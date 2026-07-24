#!/usr/bin/env python3
"""Fail-closed, exact-version NCBI assembly metadata collection freezer.

Only NCBI JSON metadata endpoints are used. This module has no sequence download code.
"""
from __future__ import annotations

import argparse, base64, csv, gzip, hashlib, json, os, platform, re, resource as usage, shutil, socket, stat, subprocess, sys, time, urllib.error, urllib.request
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

SCHEMA = "collection-release-v1"
EXPECTED_ACCESSION_SHA = "1570cbe51da509445ba1536da0a8552202d6a87566db67305903ffcced15e8c5"
EXPECTED_PROPHAGE_SHA = "6c8e3df4fce6bfa1119a594f83b758bd4ce097a22d9ee25289fdf2fbcb362996"
EXPECTED_LINES = 26078
EXPECTED_CANDIDATES = 26077
ACCESSION_RE = re.compile(rb"GCF_[0-9]{9}\.[1-9][0-9]*\Z")
CHECK_URL = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/check"
REPORT_URL = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/dataset_report"
STAGE_B = [("first",1),("10th_percentile",2609),("25th_percentile",6520),("median",13039),
           ("75th_percentile",19558),("90th_percentile",23469),("last_valid",26077),
           ("first_version_3",80),("first_version_4",4531),("median_version_2",1506)]
OCC_FIELDS = ["schema_version","record_status","record_type","record_id","parent_record_id",
 "project_release_id","source_task_id","input_source_file","input_source_sha256","input_line_number",
 "input_raw_line_b64","input_occurrence_id","normalized_accession","duplicate_group","resolution_status",
 "requested_assembly_accession_version","resolved_assembly_accession_base","assembly_version",
 "resolved_assembly_accession_version","assembly_id","artifact_state","row_sha256"]
ASM_FIELDS = ["schema_version","record_status","record_type","record_id","parent_record_id",
 "project_release_id","source_task_id","input_source_file","input_source_sha256","input_line_number",
 "input_occurrence_id","normalized_accession","resolution_status","requested_assembly_accession_version",
 "resolved_assembly_accession_base","assembly_version","resolved_assembly_accession_version",
 "paired_assembly_accession_version","assembly_status","supersedes_assembly_accession_version",
 "superseded_by_assembly_accession_version","assembly_id","biosample_accession","strain_raw_b64",
 "isolate_raw_b64","organism_name_raw_b64","annotation_available","annotation_name_raw_b64",
 "source_url","source_response_relpath","source_response_sha256","metadata_retrieved_at_utc",
 "sample_id_basis","pansn_sample","analysis_unit_id","assembly_qc_status","assembly_qc_reason",
 "artifact_state","row_sha256"]

class GateError(RuntimeError): pass

def utcnow() -> str: return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")
def sha_bytes(b: bytes) -> str: return hashlib.sha256(b).hexdigest()
def sha_file(p: Path) -> str:
    h=hashlib.sha256()
    with p.open("rb") as f:
        for block in iter(lambda:f.read(1024*1024), b""): h.update(block)
    return h.hexdigest()
def canonical(obj: Any) -> bytes: return (json.dumps(obj,sort_keys=True,separators=(",",":"),ensure_ascii=False)+"\n").encode()
def stable_id(prefix: str, *parts: str) -> str:
    return prefix+hashlib.sha256((prefix+"\0"+"\0".join(parts)).encode()).hexdigest()[:24]
def b64(value: Any) -> str:
    if value in (None, ""): return "."
    return base64.b64encode(str(value).encode()).decode()
def row_hash(row: dict[str,str], fields: list[str]) -> str:
    return sha_bytes(("\t".join(row.get(k,".") for k in fields if k!="row_sha256")+"\n").encode())
def write_tsv(path: Path, fields: list[str], rows: list[dict[str,str]]) -> None:
    with path.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter="\t",lineterminator="\n",extrasaction="ignore")
        w.writeheader()
        for row in rows:
            row={k:("." if row.get(k) in (None,"") else str(row[k])) for k in fields}
            row["row_sha256"]=row_hash(row,fields); w.writerow(row)
def append_jsonl(path: Path, obj: dict[str,Any]) -> None:
    with path.open("ab") as f: f.write(canonical(obj)); f.flush(); os.fsync(f.fileno())

def verify_external(final:Path) -> None:
    complete=(final/"COMPLETE").read_text().split()[0]
    if complete!=sha_file(final/"SHA256SUMS"): raise GateError("external COMPLETE inventory mismatch")
    for line in (final/"SHA256SUMS").read_text().splitlines():
        digest,rel=line.split("  ",1)
        if sha_file(final/rel)!=digest: raise GateError(f"external release checksum mismatch: {rel}")

def publish_tracked(final:Path, tracked:Path) -> None:
    tracked.mkdir(parents=True,exist_ok=True)
    for old in tracked.iterdir():
        if old.is_file(): old.unlink()
    for name in ("occurrences.tsv","assemblies.tsv"):
        with (final/"manifests"/name).open("rb") as src,(tracked/(name+".gz")).open("wb") as raw:
            with gzip.GzipFile(fileobj=raw,mode="wb",mtime=0) as dst: shutil.copyfileobj(src,dst)
    shutil.copyfile(final/"manifests"/"stage_b_10.tsv",tracked/"stage_b_10.tsv")
    shutil.copyfile(final/"release.json",tracked/"release.json")
    sums=[]
    for p in sorted(tracked.iterdir()):
        if p.is_file() and p.name!="SHA256SUMS": sums.append(f"{sha_file(p)}  {p.name}\n")
    (tracked/"SHA256SUMS").write_text("".join(sums))

def read_input(path: Path, expected_sha: str) -> tuple[bytes,list[bytes]]:
    data=path.read_bytes()
    if sha_bytes(data)!=expected_sha: raise GateError(f"input checksum mismatch: {path}")
    if not data.endswith(b"\n") or b"\r" in data or b"\0" in data: raise GateError("accession input is not LF-only, LF-terminated ASCII")
    try: data.decode("ascii")
    except UnicodeDecodeError as e: raise GateError("accession input is not ASCII") from e
    return data,data[:-1].split(b"\n")

@dataclass
class Resources:
    assigned_ram:int; durable_allocation:int; scratch_allocation:int; inode_allocation:int
    predicted_peak_disk:int; predicted_files:int; unfinished_write:int

def swap_free_bytes() -> int:
    for line in Path("/proc/meminfo").read_text().splitlines():
        if line.startswith("SwapFree:"): return int(line.split()[1])*1024
    raise GateError("cannot read SwapFree allocation evidence")

def preflight(path: Path, scratch: Path, r: Resources, require_floors: bool=True) -> dict[str,Any]:
    vals=[r.assigned_ram,r.durable_allocation,r.scratch_allocation,r.inode_allocation,r.predicted_peak_disk,r.predicted_files,r.unfinished_write]
    if any(v<=0 for v in vals): raise GateError("resource allocation/reservation is blank or non-positive")
    def existing_parent(p:Path) -> Path:
        while not p.exists() and p!=p.parent: p=p.parent
        return p
    target=existing_parent(path); scratch_target=existing_parent(scratch)
    st=os.statvfs(target); ss=os.statvfs(scratch_target)
    durable_free=st.f_bavail*st.f_frsize; durable_inodes=st.f_favail
    scratch_free=ss.f_bavail*ss.f_frsize; scratch_inodes=ss.f_favail
    checks={
      "durable_start_floor": durable_free>=2_400_000_000_000 and durable_inodes>=1_000_000,
      "scratch_start_floor": scratch_free>=4_000_000_000_000 and scratch_inodes>=5_000_000,
      "disk_allocation": r.predicted_peak_disk*100<=r.durable_allocation*70,
      "inode_allocation": r.predicted_files*2<=r.inode_allocation,
      "unfinished_reserve": durable_free-r.predicted_peak_disk>=2*r.unfinished_write,
    }
    if (require_floors and not all(checks.values())) or (not require_floors and not all(list(checks.values())[2:])):
        raise GateError("resource gate NO_GO: "+json.dumps(checks,sort_keys=True))
    probe=target/(f".collection-release-probe.{os.getpid()}")
    probe.write_bytes(b"probe\n"); probe.unlink()
    def mount(p:Path) -> str:
        return subprocess.run(["findmnt","-T",str(p),"-o","TARGET,SOURCE,FSTYPE,OPTIONS","-n"],check=True,text=True,capture_output=True).stdout.strip()
    return {"verdict":"PASS","path":str(path),"resolved_probe_parent":str(target),"scratch":str(scratch),
      "durable_findmnt":mount(target),"scratch_findmnt":mount(scratch_target),
      "durable_owner":{"uid":target.stat().st_uid,"gid":target.stat().st_gid,"mode":stat.S_IMODE(target.stat().st_mode)},
      "scratch_owner":{"uid":scratch_target.stat().st_uid,"gid":scratch_target.stat().st_gid,"mode":stat.S_IMODE(scratch_target.stat().st_mode)},"write_probe":"PASS",
      "durable_free_bytes":durable_free,"durable_free_inodes":durable_inodes,"scratch_free_bytes":scratch_free,
      "scratch_free_inodes":scratch_inodes,"swap_free_bytes":swap_free_bytes(),"allocations":r.__dict__,"checks":checks,"captured_at_utc":utcnow()}

def request_batch(endpoint: str, accessions:list[str], request_path:Path, response_path:Path,
                  receipt_path:Path, retries:int, sleep_seconds:float) -> dict[str,Any]:
    body=canonical({"accessions":accessions,"page_size":1000})
    if response_path.exists() and receipt_path.exists() and request_path.exists():
        receipt=json.loads(receipt_path.read_text())
        if sha_file(request_path)==receipt["request_sha256"] and sha_file(response_path)==receipt["response_sha256"]:
            return receipt
        stamp=str(time.time_ns())
        for p in (request_path,response_path,receipt_path): p.rename(p.with_name(p.name+".invalid."+stamp))
    request_path.write_bytes(body)
    last=None
    for attempt in range(1,retries+1):
        try:
            req=urllib.request.Request(endpoint,data=body,method="POST",headers={"Content-Type":"application/json","Accept":"application/json","User-Agent":"phind-collection-resolver/1.0"})
            started=utcnow()
            with urllib.request.urlopen(req,timeout=120) as resp:
                raw=resp.read(); headers=dict(resp.headers.items()); status=resp.status
            json.loads(raw); response_path.write_bytes(raw)
            receipt={"endpoint":endpoint,"method":"POST","status":status,"attempt":attempt,"started_at_utc":started,
              "completed_at_utc":utcnow(),"request_sha256":sha_file(request_path),"response_sha256":sha_file(response_path),
              "response_bytes":len(raw),"response_headers":headers,"accession_count":len(accessions)}
            receipt_path.write_bytes(canonical(receipt)); time.sleep(sleep_seconds); return receipt
        except (urllib.error.URLError,TimeoutError,json.JSONDecodeError) as e:
            last=e
            if attempt==retries: break
            time.sleep(min(60,2**(attempt-1)))
    raise GateError(f"metadata request failed after {retries} attempts: {last}")

def freeze(args: argparse.Namespace) -> dict[str,Any]:
    root=Path(args.repo_root).resolve(); accession=root/"26k_ecoli_accession.txt"; prophage=root/"26k_prophage1.csv"; template=root/"artifacts/project_manifest_template.tsv"
    _,lines=read_input(accession,EXPECTED_ACCESSION_SHA)
    if sha_file(prophage)!=EXPECTED_PROPHAGE_SHA: raise GateError("prophage input checksum mismatch")
    template_sha=sha_file(template)
    if len(lines)!=EXPECTED_LINES: raise GateError(f"physical line mismatch: {len(lines)}")
    valid=[x.decode() for x in lines if ACCESSION_RE.fullmatch(x)]; rejects=[(i+1,x) for i,x in enumerate(lines) if not ACCESSION_RE.fullmatch(x)]
    if len(valid)!=EXPECTED_CANDIDATES or rejects!=[(EXPECTED_LINES,b"genome")] or len(set(valid))!=len(valid): raise GateError("cardinality/terminal reject/uniqueness reconciliation failed")
    release_seed={"schema":SCHEMA,"accession_sha256":EXPECTED_ACCESSION_SHA,"prophage_sha256":EXPECTED_PROPHAGE_SHA,"template_sha256":template_sha,"exact_version_policy":"requested-version-is-identity","candidate_count":len(valid)}
    release_id="collection-v1-"+sha_bytes(canonical(release_seed))[:16]
    final=Path(args.durable_root)/release_id; stage=Path(args.durable_root)/(".stage."+release_id+"."+args.run_id)
    tracked=Path(args.tracked_root)
    if final.exists():
        if not (final/"COMPLETE").is_file(): raise GateError("existing external release lacks COMPLETE")
        verify_external(final); publish_tracked(final,tracked)
        return json.loads((final/"release.json").read_text())
    res=Resources(args.assigned_ram_bytes,args.durable_allocation_bytes,args.scratch_allocation_bytes,args.inode_allocation,args.predicted_peak_disk_bytes,args.predicted_files,args.unfinished_write_bytes)
    resource=preflight(Path(args.durable_root),Path(args.scratch_root),res,not args.test_relaxed_floors)
    Path(args.durable_root).mkdir(parents=True,exist_ok=True)
    resumed=False
    if stage.exists():
        if (stage/"COMPLETE").exists(): raise GateError("interrupted stage contains COMPLETE but was not promoted")
        seed_path=stage/"input_manifest.json"
        if not seed_path.exists() or seed_path.read_bytes()!=canonical(release_seed): shutil.rmtree(stage)
        else: resumed=True
    (stage/"raw"/"check").mkdir(parents=True,exist_ok=True); (stage/"raw"/"report").mkdir(parents=True,exist_ok=True); (stage/"manifests").mkdir(exist_ok=True)
    state=stage/"state.jsonl"; failures=stage/"failures.jsonl"; provenance=stage/"provenance.json"; failures.touch(exist_ok=True)
    append_jsonl(state,{"state":"RESUME_PREFLIGHT_PASS" if resumed else "PREFLIGHT_PASS","at":utcnow(),"resource":resource})
    (stage/"resource.json").write_bytes(canonical(resource))
    (stage/"input_manifest.json").write_bytes(canonical(release_seed))
    provenance.write_bytes(canonical({"schema":SCHEMA,"run_id":args.run_id,"argv":sys.argv,"python":sys.version,"platform":platform.platform(),"hostname":socket.gethostname(),"pid":os.getpid(),"environment":{k:os.environ.get(k,".") for k in ("USER","LANG","WG_TASK_ID","WG_AGENT_ID","WG_MODEL","WG_TIER")},"api_endpoints":[CHECK_URL,REPORT_URL],"sequence_payload_downloads":0,"created_at_utc":utcnow()}))
    metadata:dict[str,dict[str,Any]]={}; valid_check:set[str]=set(); response_ref:dict[str,tuple[str,str,str]]={}
    try:
      for n,start in enumerate(range(0,len(valid),args.batch_size),1):
        batch=valid[start:start+args.batch_size]; bid=f"batch-{n:04d}"
        batch_resource=preflight(Path(args.durable_root),Path(args.scratch_root),res,not args.test_relaxed_floors)
        append_jsonl(state,{"state":"BATCH_STARTED","batch_id":bid,"count":len(batch),"resource":batch_resource,"at":utcnow()})
        for kind,url in (("check",CHECK_URL),("report",REPORT_URL)):
            d=stage/"raw"/kind; req=d/(bid+".request.json"); resp=d/(bid+".response.json"); receipt=d/(bid+".receipt.json")
            rec=request_batch(url,batch,req,resp,receipt,args.retries,args.sleep_seconds)
            obj=json.loads(resp.read_bytes())
            if kind=="check": valid_check.update(obj.get("valid_assemblies",[]))
            else:
                for report in obj.get("reports",[]):
                    acc=report.get("accession")
                    if acc not in batch: raise GateError(f"API version identity violation: unrequested report {acc}")
                    if acc in metadata: raise GateError(f"duplicate API report {acc}")
                    metadata[acc]=report; response_ref[acc]=(str(resp.relative_to(stage)),rec["response_sha256"],rec["completed_at_utc"])
        append_jsonl(state,{"state":"BATCH_VALIDATED","batch_id":bid,"at":utcnow()})
      if valid_check!=set(valid): raise GateError(f"exact check mismatch missing={len(set(valid)-valid_check)} extra={len(valid_check-set(valid))}")
      append_jsonl(state,{"state":"ALL_METADATA_VALIDATED","check_count":len(valid_check),"report_count":len(metadata),"at":utcnow()})
      occurrences=[]; assemblies=[]
      for line_no,raw in enumerate(lines,1):
        token=raw.decode(); occ=stable_id("occ-v1-",EXPECTED_ACCESSION_SHA,str(line_no),base64.b64encode(raw).decode())
        common={"schema_version":"project-manifest-v1","record_status":"FROZEN","project_release_id":release_id,"source_task_id":"freeze-collection-release","input_source_file":"26k_ecoli_accession.txt","input_source_sha256":EXPECTED_ACCESSION_SHA,"input_line_number":str(line_no),"input_occurrence_id":occ,"normalized_accession":token,"artifact_state":"COMPLETE"}
        if line_no==EXPECTED_LINES:
            occurrences.append(common|{"record_type":"INPUT_OCCURRENCE","record_id":occ,"parent_record_id":".","input_raw_line_b64":b64(token),"duplicate_group":".","resolution_status":"REJECTED_MALFORMED_TERMINAL_TOKEN"}); continue
        acc=token; base,ver=acc.split("."); aid=stable_id("asm-v1-",acc); report=metadata.get(acc); status="EXACT_VERSION_RESOLVED" if report else "EXACT_VERSION_VALID_METADATA_UNAVAILABLE"
        occurrences.append(common|{"record_type":"INPUT_OCCURRENCE","record_id":occ,"parent_record_id":".","input_raw_line_b64":b64(token),"duplicate_group":".","resolution_status":status,"requested_assembly_accession_version":acc,"resolved_assembly_accession_base":base,"assembly_version":ver,"resolved_assembly_accession_version":acc,"assembly_id":aid})
        ai=(report or {}).get("assembly_info",{}); bio=ai.get("biosample",{}); infra=(report or {}).get("organism",{}).get("infraspecific_names",{}); current=(report or {}).get("current_accession",".")
        strain=bio.get("strain") or infra.get("strain"); isolate=infra.get("isolate")
        rel,digest,retrieved=response_ref.get(acc,(".",".","."))
        assemblies.append(common|{"record_type":"ASSEMBLY","record_id":aid,"parent_record_id":occ,"resolution_status":status,"requested_assembly_accession_version":acc,"resolved_assembly_accession_base":base,"assembly_version":ver,"resolved_assembly_accession_version":acc,"paired_assembly_accession_version":(report or {}).get("paired_accession","."),"assembly_status":ai.get("assembly_status","METADATA_UNAVAILABLE"),"supersedes_assembly_accession_version":".","superseded_by_assembly_accession_version":current if current not in (".",acc) else ".","assembly_id":aid,"biosample_accession":bio.get("accession","."),"strain_raw_b64":b64(strain),"isolate_raw_b64":b64(isolate),"organism_name_raw_b64":b64((report or {}).get("organism",{}).get("organism_name")),"annotation_available":"true" if (report or {}).get("annotation_info") else ("false" if report else "unknown"),"annotation_name_raw_b64":b64((report or {}).get("annotation_info",{}).get("name")),"source_url":REPORT_URL,"source_response_relpath":rel,"source_response_sha256":digest,"metadata_retrieved_at_utc":retrieved,"sample_id_basis":"exact_resolved_assembly_accession_version","pansn_sample":acc,"analysis_unit_id":aid,"assembly_qc_status":"NOT_RUN_METADATA_ONLY","assembly_qc_reason":"Phase 0 metadata resolution only"})
      write_tsv(stage/"manifests"/"occurrences.tsv",OCC_FIELDS,occurrences); write_tsv(stage/"manifests"/"assemblies.tsv",ASM_FIELDS,assemblies)
      stage_rows=[]
      byline={int(r["input_line_number"]):r for r in assemblies}
      for order,(label,line) in enumerate(STAGE_B,1):
          r=byline[line]; stage_rows.append({"stage_b_order":str(order),"selection_rule":label,"input_line_number":str(line),"input_occurrence_id":r["input_occurrence_id"],"assembly_id":r["assembly_id"],"requested_assembly_accession_version":r["requested_assembly_accession_version"],"resolved_assembly_accession_version":r["resolved_assembly_accession_version"],"resolution_status":r["resolution_status"]})
      stage_fields=list(stage_rows[0]); write_tsv(stage/"manifests"/"stage_b_10.tsv",stage_fields+["row_sha256"],stage_rows)
      def data_rows(p:Path) -> int:
          with p.open(encoding="utf-8") as f: return sum(1 for _ in f)-1
      manifest_digests={p.name:{"sha256":sha_file(p),"bytes":p.stat().st_size,"rows":data_rows(p)} for p in sorted((stage/"manifests").glob("*.tsv"))}
      if manifest_digests["occurrences.tsv"]["rows"]!=EXPECTED_LINES or manifest_digests["assemblies.tsv"]["rows"]!=EXPECTED_CANDIDATES or manifest_digests["stage_b_10.tsv"]["rows"]!=len(STAGE_B): raise GateError("publication row-count gate failed")
      if sha_file(accession)!=EXPECTED_ACCESSION_SHA or sha_file(prophage)!=EXPECTED_PROPHAGE_SHA or sha_file(template)!=template_sha: raise GateError("finish source immutability gate failed")
      append_jsonl(state,{"state":"FINISH_SOURCE_IMMUTABILITY_PASS","at":utcnow(),"input_sha256":EXPECTED_ACCESSION_SHA,"prophage_sha256":EXPECTED_PROPHAGE_SHA,"template_sha256":template_sha})
      resource_end=preflight(Path(args.durable_root),Path(args.scratch_root),res,not args.test_relaxed_floors)
      max_rss_bytes=usage.getrusage(usage.RUSAGE_SELF).ru_maxrss*1024
      resource_end["peak_rss_bytes"]=max_rss_bytes; resource_end["peak_rss_limit_bytes"]=res.assigned_ram*70//100
      resource_end["swap_growth_bytes"]=max(0,resource["swap_free_bytes"]-resource_end["swap_free_bytes"])
      if max_rss_bytes*100>res.assigned_ram*70 or resource_end["swap_growth_bytes"]!=0: raise GateError("RAM/swap end gate NO_GO")
      (stage/"resource_end.json").write_bytes(canonical(resource_end))
      release={"schema_version":SCHEMA,"release_id":release_id,"verdict":"PASS","immutable":True,"source_task_id":"freeze-collection-release","created_at_utc":utcnow(),"inputs":{"26k_ecoli_accession.txt":EXPECTED_ACCESSION_SHA,"26k_prophage1.csv":EXPECTED_PROPHAGE_SHA,"artifacts/project_manifest_template.tsv":template_sha},"counts":{"physical_occurrences":EXPECTED_LINES,"candidate_exact_versioned_assemblies":EXPECTED_CANDIDATES,"malformed_rejects":1,"exact_check_valid":len(valid_check),"detailed_metadata_reports":len(metadata),"metadata_unavailable_but_exact_check_valid":len(valid)-len(metadata),"stage_b":len(STAGE_B),"distinct_sequence_bearing_assemblies":0,"global_distinct_assembly_cap":1000},"exact_version_policy":"requested version is identity; never substitute latest","manifests":manifest_digests,"external_release_path":str(final),"raw_cache_relative_path":"raw","sequence_payload_downloads":0,"applicable_gates":{"source_immutability":"PASS","accession_version_identity":"PASS","row_accounting":"PASS","metadata_provenance":"PASS","resource":"PASS","global_cap":"PASS","bgzf_index_name_roundtrip":"NOT_APPLICABLE_METADATA_ONLY","source_coordinate_policy":"NOT_APPLICABLE_METADATA_ONLY","scale_trend":"NOT_APPLICABLE_METADATA_ONLY"}}
      (stage/"release.json").write_bytes(canonical(release)); append_jsonl(state,{"state":"READY_TO_PROMOTE","at":utcnow()})
      inventory=[]
      for p in sorted(x for x in stage.rglob("*") if x.is_file() and x.name not in ("SHA256SUMS","COMPLETE")):
          inventory.append(f"{sha_file(p)}  {p.relative_to(stage)}\n")
      (stage/"SHA256SUMS").write_text("".join(inventory))
      for p in [x for x in stage.rglob("*") if x.is_file()]:
          with p.open("rb") as f: os.fsync(f.fileno())
      (stage/"COMPLETE").write_text(sha_file(stage/"SHA256SUMS")+"  SHA256SUMS\n")
      with (stage/"COMPLETE").open("rb") as f: os.fsync(f.fileno())
      fd=os.open(stage,os.O_RDONLY); os.fsync(fd); os.close(fd)
      os.rename(stage,final); fd=os.open(final.parent,os.O_RDONLY); os.fsync(fd); os.close(fd)
      verify_external(final); publish_tracked(final,tracked)
      return release
    except Exception as e:
      append_jsonl(failures,{"state":"FAILED","at":utcnow(),"type":type(e).__name__,"message":str(e)})
      raise

def parser() -> argparse.ArgumentParser:
 p=argparse.ArgumentParser(description=__doc__); sub=p.add_subparsers(dest="command",required=True); f=sub.add_parser("freeze")
 f.add_argument("--repo-root",default="."); f.add_argument("--tracked-root",default="manifests/collection-v1"); f.add_argument("--durable-root",required=True); f.add_argument("--scratch-root",required=True); f.add_argument("--run-id",required=True)
 f.add_argument("--batch-size",type=int,default=1000); f.add_argument("--retries",type=int,default=8); f.add_argument("--sleep-seconds",type=float,default=0.34)
 f.add_argument("--assigned-ram-bytes",type=int,required=True); f.add_argument("--durable-allocation-bytes",type=int,required=True); f.add_argument("--scratch-allocation-bytes",type=int,required=True); f.add_argument("--inode-allocation",type=int,required=True); f.add_argument("--predicted-peak-disk-bytes",type=int,required=True); f.add_argument("--predicted-files",type=int,required=True); f.add_argument("--unfinished-write-bytes",type=int,required=True); f.add_argument("--test-relaxed-floors",action="store_true",help=argparse.SUPPRESS)
 return p

def main() -> int:
 try:
  args=parser().parse_args(); result=freeze(args); print(json.dumps(result,sort_keys=True)); return 0
 except GateError as e: print(f"NO_GO: {e}",file=sys.stderr); return 2
if __name__=="__main__": raise SystemExit(main())
