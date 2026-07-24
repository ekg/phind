#!/usr/bin/env python3
"""Independent semantic/checksum validator for collection-v1."""
import argparse,csv,gzip,json
from collections import Counter
from pathlib import Path
from workflow.collection_release.resolver import EXPECTED_ACCESSION_SHA,EXPECTED_PROPHAGE_SHA,sha_file,verify_external,utcnow

def rows(path:Path):
    with gzip.open(path,"rt",newline="") as f: return list(csv.DictReader(f,delimiter="\t"))
def main():
    p=argparse.ArgumentParser(); p.add_argument("--repo-root",default="."); p.add_argument("--output",default="workflow/collection_release/validation.json"); a=p.parse_args(); root=Path(a.repo_root)
    release=json.loads((root/"manifests/collection-v1/release.json").read_text()); external=Path(release["external_release_path"]); verify_external(external)
    occ=rows(root/"manifests/collection-v1/occurrences.tsv.gz"); asm=rows(root/"manifests/collection-v1/assemblies.tsv.gz")
    access=root/"26k_ecoli_accession.txt"; prop=root/"26k_prophage1.csv"; tokens=access.read_text().splitlines()
    assert sha_file(access)==EXPECTED_ACCESSION_SHA and sha_file(prop)==EXPECTED_PROPHAGE_SHA
    assert len(occ)==26078 and len(asm)==26077 and tokens[-1]=="genome"
    assert [x["normalized_accession"] for x in occ]==tokens
    assert len({x["input_occurrence_id"] for x in occ})==26078 and len({x["assembly_id"] for x in asm})==26077
    assert all(x["requested_assembly_accession_version"]==x["resolved_assembly_accession_version"] for x in asm)
    assert occ[-1]["resolution_status"]=="REJECTED_MALFORMED_TERMINAL_TOKEN"
    status=Counter(x["assembly_status"] for x in asm); resolution=Counter(x["resolution_status"] for x in asm)
    assert status==Counter({"current":25291,"suppressed":714,"METADATA_UNAVAILABLE":72})
    assert resolution==Counter({"EXACT_VERSION_RESOLVED":26005,"EXACT_VERSION_VALID_METADATA_UNAVAILABLE":72})
    with (root/"manifests/collection-v1/stage_b_10.tsv").open(newline="") as f: stage=list(csv.DictReader(f,delimiter="\t"))
    assert len(stage)==10 and [int(x["stage_b_order"]) for x in stage]==list(range(1,11))
    for line in (root/"manifests/collection-v1/SHA256SUMS").read_text().splitlines():
        digest,rel=line.split("  ",1); assert sha_file(root/"manifests/collection-v1"/rel)==digest
    result={"verdict":"PASS","validated_at_utc":utcnow(),"release_id":release["release_id"],"external_inventory_sha256":sha_file(external/"SHA256SUMS"),"counts":{"occurrences":len(occ),"input_candidates":len(asm),"current":status["current"],"suppressed_terminal":status["suppressed"],"metadata_unavailable_exact_valid":status["METADATA_UNAVAILABLE"],"stage_b":len(stage),"sequence_payloads":0},"checks":["root_input_sha256","external_COMPLETE_and_SHA256SUMS","tracked_SHA256SUMS","lossless_occurrence_order","exact_requested_equals_resolved","deterministic_unique_ids","terminal_reject","metadata_status_accounting","stage_b_order"]}
    Path(a.output).write_text(json.dumps(result,sort_keys=True,indent=2)+"\n"); print(json.dumps(result,sort_keys=True))
if __name__=="__main__": main()
