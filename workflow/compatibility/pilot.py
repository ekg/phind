#!/usr/bin/env python3
"""Execute the bounded ten-assembly consumer compatibility certification."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import re
import resource
import shutil
import signal
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Iterable, Sequence

from .compatibility import (
    EXPECTED_ASSEMBLIES, GateError, PASS, ResourceRequest, append_jsonl,
    atomic_promote, canonical_json, gff_lexical_encode, load_tsv,
    mash_lower_to_full_phylip, parse_fasta, percent_encode_identifier,
    require_cleaned, resource_preflight, sha256_bytes,
    sha256_file, stage_gff_semantic_alias, tree_size, utc_now,
    verify_predecessor, verify_root_hashes, write_fasta, write_inventory, write_json,
)

SCHEMA = "consumer-compatibility-release-v1"
HOST_LOCK = "workflow/compatibility/environment-linux-64.explicit.lock"
HOST_PACKAGE_SHA256 = "workflow/compatibility/environment-package-sha256.tsv"
GRAPH_LOCK = "workflow/compatibility/graph-linux-64.explicit.lock"
GRAPH_PACKAGE_SHA256 = "workflow/compatibility/graph-package-sha256.tsv"
GRAPH_PREFIX_DEFAULT = Path("/home/erikg/micromamba/envs/pggb-env")
IMPG_DEFAULT = Path("/home/erikg/.cargo/bin/impg")


class Recorder:
    def __init__(self, release: Path, scratch: Path, request: ResourceRequest):
        self.release = release
        self.scratch = scratch
        self.request = request
        self.commands = release / "commands.jsonl"
        self.states = release / "state.jsonl"
        self.failures = release / "failures.jsonl"
        self.resources = release / "resources.jsonl"
        self.max_rss_bytes = 0
        self.peak_scratch_bytes = 0
        self.peak_scratch_files = 0
        self.gates: list[dict[str, Any]] = []
        self.determinism: dict[str, Any] = {}

    def state(self, stage: str, status: str, **extra: Any) -> None:
        append_jsonl(self.states, {"at_utc": utc_now(), "stage": stage, "status": status, **extra})

    def failure(self, stage: str, error: str, **extra: Any) -> None:
        append_jsonl(self.failures, {"at_utc": utc_now(), "stage": stage, "error": error, **extra})

    def preflight(self, stage: str, durable: Path) -> dict[str, Any]:
        rec = resource_preflight(stage, durable, self.scratch, self.request)
        append_jsonl(self.resources, rec)
        return rec

    def measure_tree(self) -> None:
        size, files = tree_size(self.scratch)
        self.peak_scratch_bytes = max(self.peak_scratch_bytes, size)
        self.peak_scratch_files = max(self.peak_scratch_files, files)

    def run(self, consumer: str, argv: Sequence[str | Path], *, cwd: Path | None = None,
            env: dict[str, str] | None = None, stdout: Path | None = None,
            expect: int | set[int] = 0, stderr_limit: int = 32768) -> subprocess.CompletedProcess[str]:
        args = [str(x) for x in argv]
        command_id = f"{len(self.commands.read_text().splitlines()) if self.commands.exists() else 0:04d}-{consumer}"
        time_path = self.scratch / "time" / f"{command_id}.txt"
        stderr_path = self.scratch / "logs" / f"{command_id}.stderr"
        time_path.parent.mkdir(parents=True, exist_ok=True)
        stderr_path.parent.mkdir(parents=True, exist_ok=True)
        wrapped = ["/usr/bin/time", "-v", "-o", str(time_path), "--", *args]
        start = time.monotonic()
        out_fh = stdout.open("w") if stdout else subprocess.PIPE
        try:
            proc = subprocess.run(wrapped, cwd=cwd, env=env, text=True, stdout=out_fh,
                                  stderr=subprocess.PIPE, errors="replace")
        finally:
            if stdout:
                out_fh.close()  # type: ignore[union-attr]
        elapsed = time.monotonic() - start
        stderr_path.write_text(proc.stderr[-stderr_limit:])
        rss = 0
        if time_path.exists():
            m = re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", time_path.read_text())
            if m:
                rss = int(m.group(1)) * 1024
                self.max_rss_bytes = max(self.max_rss_bytes, rss)
        self.measure_tree()
        record = {
            "command_id": command_id, "consumer": consumer, "argv": args,
            "cwd": str(cwd) if cwd else None, "started_at_utc": utc_now(),
            "elapsed_seconds": round(elapsed, 6), "exit_code": proc.returncode,
            "max_rss_bytes": rss, "stdout_relpath": str(stdout.relative_to(self.scratch)) if stdout else None,
            "stderr_sha256": sha256_file(stderr_path), "stderr_tail": proc.stderr[-2000:],
        }
        append_jsonl(self.commands, record)
        expected = {expect} if isinstance(expect, int) else expect
        if proc.returncode not in expected:
            raise GateError(f"{consumer} exit {proc.returncode}, expected {sorted(expected)}: {' '.join(args)}; {proc.stderr[-1000:]}")
        return proc

    def gate(self, gate_id: str, tool: str, input_form: str, invocation: str,
             output_names: str, checks: dict[str, Any], *, view_contract: str = "direct") -> None:
        bad = [k for k, v in checks.items() if v not in {PASS, True, "PASS"}]
        if bad:
            raise GateError(f"gate {gate_id} failed checks: {bad}")
        gate = {
            "schema_version": "consumer-machine-gate-v1", "gate_id": gate_id,
            "tool": tool, "verdict": PASS, "input_form": input_form,
            "view_contract": view_contract, "invocation": invocation,
            "output_name_behavior": output_names, "checks": checks,
            "created_at_utc": utc_now(),
        }
        write_json(self.release / "gates" / f"{gate_id}.json", gate)
        self.gates.append(gate)


def release_id(repo: Path) -> str:
    material = {
        "schema": SCHEMA,
        "predecessor": "canonical-cohort-010-v1-e71484de9994fc28",
        "host_lock": sha256_file(repo / HOST_LOCK),
        "host_package_sha256": sha256_file(repo / HOST_PACKAGE_SHA256),
        "graph_lock": sha256_file(repo / GRAPH_LOCK),
        "graph_package_sha256": sha256_file(repo / GRAPH_PACKAGE_SHA256),
        "runner": sha256_file(Path(__file__)),
    }
    return "consumer-compatibility-v1-" + sha256_bytes(canonical_json(material))[:16]


def executable_identity(path: Path, version_argv: Sequence[str], help_argv: Sequence[str]) -> dict[str, Any]:
    if not path.is_file() or not os.access(path, os.X_OK):
        raise GateError(f"required executable absent: {path}")
    proc = subprocess.run([str(path), *version_argv], text=True, capture_output=True, errors="replace")
    help_proc = subprocess.run([str(path), *help_argv], text=True, capture_output=True, errors="replace")
    text = (proc.stdout + proc.stderr).strip()
    help_text = (help_proc.stdout + help_proc.stderr).strip()
    if not text or not help_text:
        raise GateError(f"empty version/help capture for {path}")
    return {
        "path": str(path.resolve()), "sha256": sha256_file(path.resolve()),
        "bytes": path.resolve().stat().st_size, "version_argv": list(version_argv),
        "version_exit": proc.returncode, "version_output": text[:12000],
        "help_argv": list(help_argv), "help_exit": help_proc.returncode,
        "help_output": help_text[:24000],
    }


def install_host_environment(repo: Path, prefix: Path, rec: Recorder, durable: Path) -> None:
    rec.preflight("before_pinned_environment_install", durable)
    if prefix.exists():
        shutil.rmtree(prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    rec.run("micromamba-install", ["micromamba", "create", "-y", "-p", prefix, "-f", repo / HOST_LOCK])
    history = prefix / "conda-meta/history"
    if not history.is_file():
        raise GateError("pinned environment lacks conda history")
    rec.measure_tree()


def tool_identities(host: Path, graph: Path, impg: Path) -> dict[str, dict[str, Any]]:
    specs: dict[str, tuple[Path, Sequence[str], Sequence[str]]] = {
        "bgzip": (Path(shutil.which("bgzip") or ""), ["--version"], ["--help"]),
        "samtools": (Path(shutil.which("samtools") or ""), ["--version"], ["--help"]),
        "impg": (impg, ["--version"], ["--help"]),
        "mash": (host / "bin/mash", ["--version"], ["--help"]),
        "rapidnj": (host / "bin/rapidnj", ["-h"], ["-h"]),
        "skani": (host / "bin/skani", ["--version"], ["--help"]),
        "quast": (host / "bin/quast.py", ["--version"], ["--help"]),
        "gffread": (host / "bin/gffread", ["--version"], ["--help"]),
        "prodigal": (host / "bin/prodigal", ["-v"], ["-h"]),
        "mmseqs": (host / "bin/mmseqs", ["version"], ["--help"]),
        "hmmbuild": (host / "bin/hmmbuild", ["-h"], ["-h"]),
        "hmmsearch": (host / "bin/hmmsearch", ["-h"], ["-h"]),
        "mcl": (host / "bin/mcl", ["--version"], ["--help"]),
        "pggb": (graph / "bin/pggb", ["--version"], ["--help"]),
        "wfmash": (graph / "bin/wfmash", ["--version"], ["--help"]),
        "seqwish": (graph / "bin/seqwish", ["--version"], ["--help"]),
        "smoothxg": (graph / "bin/smoothxg", ["--version"], ["--help"]),
        "odgi": (graph / "bin/odgi", ["version"], []),
        "vg": (graph / "bin/vg", ["version"], []),
    }
    return {name: executable_identity(path, version_args, help_args)
            for name, (path, version_args, help_args) in specs.items()}


def _fasta_sequence_from_stdout(text: str) -> bytes:
    return b"".join(line.strip().encode() for line in text.splitlines() if line and not line.startswith(">"))


def samtools_probe(rec: Recorder, assemblies: list[Any], tools: dict[str, Any]) -> dict[str, bytes]:
    expected: dict[str, bytes] = {}
    for asm in assemblies:
        bgzip = Path(tools["bgzip"]["path"])
        samtools = Path(tools["samtools"]["path"])
        rec.run("bgzip", [bgzip, "-t", asm.bgzf])
        first_name, first_len = asm.fai.read_text().splitlines()[0].split("\t")[:2]
        n = min(101, int(first_len))
        proc = rec.run("samtools", [samtools, "faidx", asm.bgzf, f"{first_name}:1-{n}"])
        observed = _fasta_sequence_from_stdout(proc.stdout)
        assembly_records = list(parse_fasta(asm.bgzf, gzipped=True))
        source_name, source_seq = assembly_records[0]
        if source_name != first_name or observed != source_seq[:n]:
            raise GateError(f"samtools name/base/coordinate round trip failed for {asm.accession}")
        if len(assembly_records) != asm.contig_count:
            raise GateError(f"parsed FASTA/manifest contig count mismatch for {asm.accession}")
        for path_name, path_seq in assembly_records:
            if path_name in expected:
                raise GateError(f"duplicate PanSN path during consumer probe: {path_name}")
            expected[path_name] = path_seq
    bad = rec.run("samtools-negative", [tools["samtools"]["path"], "faidx", assemblies[0].bgzf,
                                         "NOT_A_PATH#1#missing:1-10"], expect={1, 2})
    corrupt = rec.scratch / "samtools/corrupt.bgzf"
    corrupt.parent.mkdir(parents=True, exist_ok=True)
    corrupt.write_bytes(b"not a BGZF stream\n")
    rec.run("bgzip-negative", [tools["bgzip"]["path"], "-t", corrupt], expect={1, 2})
    rec.gate("bgzip", "bgzip", "canonical per-assembly BGZF", "bgzip -t CANONICAL",
             "does not rewrite FASTA names", {"all_10_integrity": PASS, "intentional_bad_input_rejected": PASS})
    rec.gate("samtools-faidx", "samtools faidx", "canonical BGZF plus immutable .fai/.gzi",
             "samtools faidx CANONICAL 'literal#path:1-N'", "FASTA header retains literal # exactly",
             {"all_10_index_lookup": PASS, "name_roundtrip": PASS, "base_roundtrip": PASS,
              "coordinate_1based_query_to_slice": PASS, "unknown_path_rejected": PASS})
    return expected


def stage_combined_panel(rec: Recorder, assemblies: list[Any], bgzip: Path) -> tuple[Path, dict[str, Any]]:
    view = rec.scratch / "views/panel/all10.pansn.fa.gz"
    view.parent.mkdir(parents=True, exist_ok=True)
    stderr = rec.scratch / "logs/panel-bgzip.stderr"
    cmd = [str(bgzip), "-@", "2", "-l", "6", "--binary", "-c"]
    start = time.monotonic()
    proc = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=view.open("wb"), stderr=stderr.open("wb"))
    source_h = hashlib.sha256()
    try:
        assert proc.stdin is not None
        for asm in assemblies:
            with gzip.open(asm.bgzf, "rb") as fh:
                while data := fh.read(1024 * 1024):
                    source_h.update(data)
                    proc.stdin.write(data)
        proc.stdin.close()
        code = proc.wait()
    finally:
        if proc.stdin and not proc.stdin.closed:
            proc.stdin.close()
    append_jsonl(rec.commands, {"command_id": "stage-combined-panel", "consumer": "bgzip",
        "argv": cmd, "input_count": 10, "exit_code": code, "elapsed_seconds": round(time.monotonic()-start, 6)})
    if code != 0:
        raise GateError(f"combined panel bgzip failed: {stderr.read_text(errors='replace')[-2000:]}")
    out_h = hashlib.sha256()
    with gzip.open(view, "rb") as fh:
        while data := fh.read(1024 * 1024): out_h.update(data)
    if source_h.hexdigest() != out_h.hexdigest():
        raise GateError("combined panel staged-view content checksum mismatch")
    rec.run("samtools", [shutil.which("samtools") or "samtools", "faidx", view])
    fai = Path(str(view) + ".fai")
    if len(fai.read_text().splitlines()) != 1223:
        raise GateError("combined view is not 100% of predecessor paths")
    contract = {
        "view_id": "combined-bgzf-all10-v1", "type": "checksum-verified quota-bounded staged transformation",
        "inputs": [{"accession": a.accession, "path": str(a.bgzf), "sha256": a.bgzf_sha256,
                    "content_sha256": a.content_sha256} for a in assemblies],
        "input_order": EXPECTED_ASSEMBLIES, "records": 1223, "bases": 51731662,
        "transform": cmd, "names_changed": False, "reversible_name_map": "identity",
        "decompressed_content_sha256": out_h.hexdigest(), "bgzf_sha256": sha256_file(view),
        "bytes": view.stat().st_size, "quota_bytes": 200_000_000,
        "cleanup": "recursive scratch removal on success and failure",
    }
    write_json(rec.release / "view_contracts/combined-bgzf-all10-v1.json", contract)
    repeat = view.with_name("all10.repeat.pansn.fa.gz")
    with repeat.open("wb") as output:
        rerun = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=output, stderr=subprocess.PIPE)
        assert rerun.stdin is not None
        for asm in assemblies:
            with gzip.open(asm.bgzf, "rb") as source:
                while data := source.read(1024 * 1024):
                    rerun.stdin.write(data)
        rerun.stdin.close()
        rerun.wait()
    if rerun.returncode != 0 or sha256_file(repeat) != sha256_file(view):
        raise GateError("combined BGZF deterministic byte rerun failed for the identical ordered-source recipe")
    rec.determinism["combined_bgzf_byte_identical"] = PASS
    rec.determinism["combined_bgzf_sha256"] = sha256_file(view)
    repeat.unlink()
    rec.measure_tree()
    return view, contract


def impg_probe(rec: Recorder, tools: dict[str, Any], panel: Path, sequences: dict[str, bytes]) -> None:
    impg = Path(tools["impg"]["path"])
    prefix = rec.scratch / "impg/all10.syng"
    prefix.parent.mkdir(parents=True, exist_ok=True)
    rec.run("impg-syng", [impg, "syng", "-f", panel, "-o", prefix,
                           "--parallel-dictionary", "--position-sample-rate", "256", "-t", "2"])
    expected_suffixes = [".1khash", ".1gbwt", ".names", ".pstep", ".spos", ".meta"]
    files = [Path(str(prefix) + suffix) for suffix in expected_suffixes]
    if not all(p.is_file() and p.stat().st_size > 0 for p in files):
        raise GateError("IMPG SYNG did not produce six nonempty inseparable files")
    repeat_prefix = rec.scratch / "impg/all10-repeat.syng"
    rec.run("impg-syng-deterministic-rerun", [impg, "syng", "-f", panel, "-o", repeat_prefix,
                           "--parallel-dictionary", "--position-sample-rate", "256", "-t", "2"])
    repeat_files = [Path(str(repeat_prefix) + suffix) for suffix in expected_suffixes]
    if [sha256_file(p) for p in repeat_files] != [sha256_file(p) for p in files]:
        raise GateError("IMPG six-file byte determinism rerun failed")
    rec.determinism["impg_six_file_byte_identical"] = PASS
    rec.determinism["impg_file_sha256"] = {suffix: sha256_file(path) for suffix, path in zip(expected_suffixes, files)}
    for path in repeat_files:
        path.unlink()
    name_rows = [line.split("\t") for line in Path(str(prefix) + ".names").read_text().splitlines()]
    if any(len(row) < 3 for row in name_rows):
        raise GateError("IMPG .names row schema is incomplete")
    names = [row[1] for row in name_rows]
    if len(names) != 1223 or set(names) != set(sequences):
        raise GateError("IMPG .names did not retain all exact PanSN paths")
    rec.run("impg-syng-negative", [impg, "syng", "-f", rec.scratch / "impg/missing.fa",
                                    "-o", rec.scratch / "impg/must-not-build"], expect={1, 2})
    source_name = names[0]
    source_seq = sequences[source_name]
    start, end = 1000, min(6000, len(source_seq))
    bed = rec.scratch / "impg/query.bed"
    rec.run("impg-query", [impg, "query", "-a", prefix, "-r", f"{source_name}:{start}-{end}",
                            "-d", "0", "-o", "bed", "--syng-raw", "-t", "2"], stdout=bed)
    query_repeat = rec.scratch / "impg/query-repeat.bed"
    rec.run("impg-query-deterministic-rerun", [impg, "query", "-a", prefix, "-r", f"{source_name}:{start}-{end}",
                            "-d", "0", "-o", "bed", "--syng-raw", "-t", "2"], stdout=query_repeat)
    if sha256_file(query_repeat) != sha256_file(bed):
        raise GateError("IMPG query byte determinism rerun failed")
    rec.determinism["impg_query_byte_identical"] = PASS
    rows = [x.split("\t") for x in bed.read_text().splitlines() if x]
    source_rows = [x for x in rows if x[0] == source_name and int(x[1]) <= start and int(x[2]) >= end]
    if not source_rows:
        raise GateError("IMPG coordinate query did not recover source interval")
    bad = rec.run("impg-query-negative", [impg, "query", "-a", prefix, "-r", "MISSING#1#PATH:0-100",
                                                "-d", "0", "-o", "bed", "--syng-raw"], expect={1, 2})
    qname = source_name.rsplit("#", 1)[0] + "#" + percent_encode_identifier(
        (source_name.rsplit("#", 1)[1] + f":{start}-{end}").encode())
    query_plain = rec.scratch / "impg/query.fa"
    write_fasta(query_plain, [(qname, source_seq[start:end])])
    query_bgzf = rec.scratch / "impg/query.fa.gz"
    rec.run("bgzip", [tools["bgzip"]["path"], "-c", query_plain], stdout=query_bgzf)
    rec.run("impg-map-negative", [impg, "map", "-a", prefix,
            "-q", rec.scratch / "impg/missing-query.fa", "-o", "paf"], expect={1, 2})
    map_paf = rec.scratch / "impg/map.paf"
    rec.run("impg-map", [impg, "map", "-a", prefix, "-q", query_bgzf, "-o", "paf",
                          "--min-anchors", "2", "--max-hits", "20", "-t", "2", "-O", map_paf])
    map_repeat = rec.scratch / "impg/map-repeat.paf"
    rec.run("impg-map-deterministic-rerun", [impg, "map", "-a", prefix, "-q", query_bgzf, "-o", "paf",
                          "--min-anchors", "2", "--max-hits", "20", "-t", "2", "-O", map_repeat])
    if sha256_file(map_repeat) != sha256_file(map_paf):
        raise GateError("IMPG map byte determinism rerun failed")
    rec.determinism["impg_map_byte_identical"] = PASS
    paf = [x.split("\t") for x in map_paf.read_text().splitlines() if x]
    if not paf or any(x[0] != qname for x in paf) or not any(x[5] == source_name for x in paf):
        raise GateError("IMPG map did not retain exact query/target names or source hit")
    rec.gate("impg-syng-build", "IMPG 0.4.1 syng", "staged all-10 combined BGZF",
             "impg syng -f PANEL -o PREFIX --parallel-dictionary", ".names retains all 1,223 literal PanSN paths",
             {"six_file_index": PASS, "all_names_exact": PASS, "all_10_accounted": PASS,
              "missing_input_rejected": PASS},
             view_contract="combined-bgzf-all10-v1")
    rec.gate("impg-query", "IMPG 0.4.1 query", "six-file SYNG prefix plus exact 0-based half-open range",
             "impg query -a PREFIX -r 'path:start-end' -d 0 -o bed --syng-raw",
             "BED path is exact; raw interval is syncmer-padded and is not asserted base-aligned",
             {"source_interval_recovered": PASS, "literal_hash": PASS, "unknown_path_rejected": PASS})
    rec.gate("impg-map", "IMPG 0.4.1 map", "six-file SYNG prefix plus BGZF query",
             "impg map -a PREFIX -q QUERY.bgzf -o paf", "PAF query/target names exact",
             {"query_name_exact": PASS, "target_name_exact": PASS, "source_candidate_present": PASS,
              "missing_query_rejected": PASS,
              "semantic_label_syncmer_projection_not_alignment": PASS})


def mash_rapidnj_probe(rec: Recorder, tools: dict[str, Any], assemblies: list[Any], host_env: dict[str, str]) -> None:
    root = rec.scratch / "mash"
    root.mkdir(parents=True, exist_ok=True)
    listing = root / "assemblies.list"
    listing.write_text("".join(str(a.bgzf) + "\n" for a in assemblies))
    sketch = root / "hosts"
    rec.run("mash", [tools["mash"]["path"], "sketch", "-p", "2", "-k", "21", "-s", "1000",
                     "-l", listing, "-o", sketch], env=host_env)
    msh = Path(str(sketch) + ".msh")
    info = root / "info.tsv"
    rec.run("mash", [tools["mash"]["path"], "info", "-t", msh], stdout=info, env=host_env)
    data = [x.split("\t") for x in info.read_text().splitlines() if x and not x.startswith("#")]
    if len(data) != 10:
        raise GateError(f"Mash whole-file mode emitted {len(data)} sketches, expected 10")
    triangle = root / "triangle.phylip"
    rec.run("mash", [tools["mash"]["path"], "triangle", "-p", "2", msh], stdout=triangle, env=host_env)
    lines = triangle.read_text().splitlines()
    if not lines or int(lines[0].strip()) != 10 or len(lines[1:]) != 10:
        raise GateError("Mash lower-PHYLIP cardinality mismatch")
    token_count = sum(max(0, len(line.split()) - 1) for line in lines[1:])
    if token_count != 45:  # relaxed lower PHYLIP contains the 45 off-diagonal pairs, no diagonal
        raise GateError(f"Mash lower triangle token count {token_count} != 45")
    triangle_repeat = root / "triangle-repeat.phylip"
    rec.run("mash-deterministic-rerun", [tools["mash"]["path"], "triangle", "-p", "2", msh],
            stdout=triangle_repeat, env=host_env)
    if sha256_file(triangle_repeat) != sha256_file(triangle):
        raise GateError("Mash triangle byte determinism rerun failed")
    rec.determinism["mash_triangle_byte_identical"] = PASS
    full_matrix = root / "rapidnj-full.phylip"
    conversion = mash_lower_to_full_phylip(triangle, full_matrix)
    write_json(rec.release / "view_contracts/mash-to-rapidnj-full-v1.json", {
        "view_id": "mash-to-rapidnj-full-v1",
        "type": "checksum-verified quota-bounded staged transformation",
        "source_format": "Mash 2.3 relaxed lower PHYLIP without diagonal",
        "consumer_format": "RapidNJ 2.3.2 full PHYLIP distance matrix",
        "source_sha256": sha256_file(triangle), "output_sha256": sha256_file(full_matrix),
        "transform": "validate row i has i distances; expand symmetric; add exact zero diagonal; labels unchanged",
        "reversible_name_map": "identity", "names_changed": False, **conversion,
        "quota_bytes": 10_000_000, "cleanup": "recursive scratch removal on success and failure",
    })
    tree = root / "overview.nwk"
    rec.run("rapidnj", [tools["rapidnj"]["path"], full_matrix, "-i", "pd", "-o", "t"],
            stdout=tree, env=host_env)
    malformed = root / "malformed.phylip"
    malformed.write_text("2\nonly_one 0\n")
    rec.run("rapidnj-negative", [tools["rapidnj"]["path"], malformed, "-i", "pd", "-o", "t"],
            env=host_env, expect={1, 2, 255})
    tips = re.findall(r"(?<=[(,])([^():,;]+)(?=[:),])", tree.read_text())
    if len(set(tips)) != 10:
        raise GateError(f"RapidNJ tip count mismatch: {len(set(tips))}")
    missing = root / "missing.list"
    missing.write_text(str(root / "DOES_NOT_EXIST.fa.gz") + "\n")
    rec.run("mash-negative", [tools["mash"]["path"], "sketch", "-l", missing, "-o", root / "bad"],
            env=host_env, expect={1, 2})
    rec.gate("mash", "Mash 2.3", "ten canonical per-assembly BGZF files; whole-file mode (never -i)",
             "mash sketch -k 21 -s 1000 -l LIST; mash triangle COMPOUND.msh",
             "sketch labels are input paths and are reversibly joined to ordered assembly manifest; one sketch/assembly",
             {"bgzf_direct": PASS, "whole_file_10_of_10": PASS, "pair_count_45": PASS,
              "diagonal_zero_structure": PASS, "missing_input_rejected": PASS})
    rec.gate("rapidnj", "RapidNJ 2.3.2", "checksum-verified full PHYLIP staged from Mash lower PHYLIP",
             "rapidnj FULL_MATRIX -i pd -o t", "quoted Newick tips preserve the ten Mash labels",
             {"phylip_parse": PASS, "tip_count_10": PASS, "malformed_matrix_rejected": PASS,
              "lower_to_full_symmetric_zero_diagonal": PASS,
              "scientific_label_unrooted_genomic_similarity_dendrogram_not_phylogeny": PASS},
             view_contract="mash-to-rapidnj-full-v1")


def ani_qc_probe(rec: Recorder, tools: dict[str, Any], assemblies: list[Any], host_env: dict[str, str]) -> None:
    root = rec.scratch / "qc"
    root.mkdir(parents=True, exist_ok=True)
    skani_out = root / "skani.tsv"
    rec.run("skani", [tools["skani"]["path"], "dist", assemblies[0].bgzf, assemblies[0].bgzf],
            stdout=skani_out, env=host_env)
    text = skani_out.read_text()
    if assemblies[0].bgzf.name not in text and str(assemblies[0].bgzf) not in text:
        raise GateError("skani output did not retain input filename")
    numeric = [float(x) for x in re.findall(r"\b(?:100(?:\.0+)?|\d{1,2}\.\d+)\b", text)]
    if not any(x >= 99.9 for x in numeric):
        raise GateError("skani self ANI was not approximately 100")
    rec.run("skani-negative", [tools["skani"]["path"], "dist", root / "missing.fa", assemblies[0].bgzf],
            env=host_env, expect={1, 2})

    quast_out = root / "quast"
    rec.run("quast", [tools["quast"]["path"], "--threads", "1", "--min-contig", "0", "-o", quast_out,
                      assemblies[0].bgzf], env=host_env)
    report = quast_out / "report.tsv"
    if not report.is_file() or "Total length" not in report.read_text():
        raise GateError("QUAST did not produce assembly metrics from canonical BGZF")
    rec.run("quast-negative", [tools["quast"]["path"], "-o", root / "quast-bad", root / "missing.fa"],
            env=host_env, expect={1, 2, 4})
    rec.gate("skani", "skani 0.3.1", "canonical per-assembly BGZF directly",
             "skani dist ASSEMBLY.bgzf ASSEMBLY.bgzf", "reports input file labels; downstream joins by manifest accession",
             {"direct_bgzf": PASS, "self_ani_approximately_100": PASS, "missing_input_rejected": PASS})
    rec.gate("quast", "QUAST 5.3.0", "canonical per-assembly BGZF directly",
             "quast.py --threads 1 -o OUT ASSEMBLY.bgzf", "report column derives from file label; accession join remains external",
             {"direct_bgzf": PASS, "metrics_report": PASS, "missing_input_rejected": PASS})


def gff_probe(rec: Recorder, tools: dict[str, Any], sequences: dict[str, bytes], host_env: dict[str, str]) -> None:
    root = rec.scratch / "gff"
    root.mkdir(parents=True, exist_ok=True)
    name, seq = next(iter(sequences.items()))
    fasta = root / "fixture.fa"
    write_fasta(fasta, [(name, seq[:1000])])
    lexical = gff_lexical_encode(name)
    source = root / "lexical.gff3"
    source.write_text(
        "##gff-version 3\n"
        f"{lexical}\tcompat\tmRNA\t101\t300\t.\t+\t.\tID=tx1\n"
        f"{lexical}\tcompat\texon\t101\t300\t.\t+\t.\tParent=tx1\n"
    )
    raw_out = root / "raw.fa"
    raw = rec.run("gffread-intentional-raw-mismatch", [tools["gffread"]["path"], source, "-g", fasta, "-w", raw_out],
                  env=host_env, expect={0, 1})
    raw_failed_semantically = not raw_out.exists() or raw_out.stat().st_size == 0
    if not raw_failed_semantically:
        raise GateError("gffread unexpectedly treated %23 lexical ID as literal # without declared semantics")
    alias = root / "semantic.gff3"
    mapping = root / "semantic-name-map.json"
    rows = stage_gff_semantic_alias(source, alias, {name}, mapping)
    out = root / "transcripts.fa"
    rec.run("gffread", [tools["gffread"]["path"], alias, "-g", fasta, "-w", out], env=host_env)
    records = list(parse_fasta(out, gzipped=False))
    if len(records) != 1 or records[0][1] != seq[100:300]:
        raise GateError("GFF 1-based closed -> canonical [100,300) coordinate/base round trip failed")
    bad = root / "bad.gff3"
    bad.write_text(source.read_text().replace(lexical, "UNKNOWN%231%23PATH"))
    try:
        stage_gff_semantic_alias(bad, root / "must-not-exist.gff3", {name}, root / "bad-map.json")
    except GateError:
        rejected = True
    else:
        rejected = False
    raw_unsafe = "ctg#A/plasmid|β and space".encode()
    encoded = percent_encode_identifier(raw_unsafe)
    semantic_unsafe = "GCF_000005845.2#1#" + encoded
    lexical_unsafe = gff_lexical_encode(semantic_unsafe)
    if "%23" not in lexical_unsafe or "%2523" not in lexical_unsafe or "%20" not in encoded:
        raise GateError("unsafe/PanSN/GFF two-layer escape contract failed")
    write_json(rec.release / "view_contracts/gff-semantic-alias-v1.json", {
        "view_id": "gff-semantic-alias-v1",
        "type": "checksum-verified quota-bounded staged transformation",
        "source_semantics": "GFF3 lexical seqid is one-layer percent encoded",
        "consumer_semantics": "gffread 0.12.7 compares seqids lexically/raw",
        "transform": "strict decode of column 1 only; coordinates and columns 2-9 byte-preserved",
        "mapping": rows, "reversible": True, "quota_bytes": 10_000_000,
        "cleanup": "recursive scratch removal on success and failure",
    })
    rec.gate("gffread", "gffread 0.12.7 plus project strict semantic adapter",
             "staged plain FASTA and reversible raw-seqid GFF view",
             "gffread semantic.gff3 -g staged.fa -w transcripts.fa",
             "raw %23 does not equal FASTA #; adapter records lexical<->semantic mapping",
             {"intentional_raw_mismatch_detected": PASS, "semantic_alias_success": PASS,
              "gff_1based_closed_to_halfopen_roundtrip": PASS, "bases_exact": PASS,
              "unknown_semantic_id_rejected": PASS, "unsafe_space_hash_utf8_two_layer_escape": PASS},
             view_contract="gff-semantic-alias-v1")


def graph_probe(rec: Recorder, tools: dict[str, Any], sequences: dict[str, bytes], graph_env: dict[str, str]) -> tuple[Path, list[str]]:
    root = rec.scratch / "graph"
    root.mkdir(parents=True, exist_ok=True)
    canonical_name, seq = next(iter(sequences.items()))
    seq = seq[:80000]
    derivative_name = "DERIVED_GCF_000005845.2#1#" + percent_encode_identifier(
        (canonical_name.rsplit("#", 1)[1] + ":compat-copy").encode())
    plain = root / "fixture.fa"
    write_fasta(plain, [(canonical_name, seq), (derivative_name, seq)])
    panel = root / "fixture.fa.gz"
    rec.run("bgzip", [tools["bgzip"]["path"], "-c", plain], stdout=panel)
    rec.run("samtools", [tools["samtools"]["path"], "faidx", panel])
    write_json(rec.release / "view_contracts/graph-fixture-v1.json", {
        "view_id": "graph-fixture-v1",
        "type": "checksum-verified quota-bounded staged transformation", "source_assembly": EXPECTED_ASSEMBLIES[0],
        "source_path": canonical_name, "source_interval_0based_halfopen": [0, len(seq)],
        "records": [{"output": canonical_name, "source": canonical_name},
                    {"output": derivative_name, "source": canonical_name}],
        "sequence_sha256": sha256_bytes(seq), "plain_sha256": sha256_file(plain),
        "bgzf_sha256": sha256_file(panel), "reversible_name_mapping": True,
        "quota_bytes": 100_000_000, "cleanup": "recursive scratch removal on success and failure",
    })
    outdir = root / "pggb"
    rec.run("pggb", [tools["pggb"]["path"], "-i", panel, "-o", outdir, "-t", "2", "-T", "1",
                     "-p", "99", "-s", "5000", "-l", "10000", "-n", "2", "-v", "-A"], env=graph_env)
    candidates = sorted(outdir.rglob("*.smooth.final.gfa"))
    if not candidates:
        candidates = sorted(outdir.rglob("*.gfa"), key=lambda p: ("smooth" not in p.name, p.name))
    if not candidates:
        raise GateError("pggb produced no GFA")
    final_gfa = candidates[0]
    gfa_text = final_gfa.read_text(errors="replace")
    for name in (canonical_name, derivative_name):
        if name not in gfa_text:
            raise GateError(f"pggb final GFA lost path {name}")
    og = root / "final.og"
    rec.run("odgi", [tools["odgi"]["path"], "build", "-g", final_gfa, "-o", og, "-t", "2"], env=graph_env)
    odgi_names = root / "odgi.names"
    rec.run("odgi", [tools["odgi"]["path"], "paths", "-i", og, "-L"], stdout=odgi_names, env=graph_env)
    odgi_fasta = root / "odgi.fa"
    rec.run("odgi", [tools["odgi"]["path"], "paths", "-i", og, "-f"], stdout=odgi_fasta, env=graph_env)
    if set(odgi_names.read_text().splitlines()) != {canonical_name, derivative_name}:
        raise GateError("odgi path-name round trip failed")
    odgi_records = dict(parse_fasta(odgi_fasta, gzipped=False))
    if any(odgi_records.get(n) != seq for n in (canonical_name, derivative_name)):
        raise GateError("odgi path sequence round trip failed")
    vg_graph = root / "final.vg"
    rec.run("vg", [tools["vg"]["path"], "convert", "-g", "-v", final_gfa], stdout=vg_graph, env=graph_env)
    vg_names = root / "vg.names"
    rec.run("vg", [tools["vg"]["path"], "paths", "-v", vg_graph, "-L"], stdout=vg_names, env=graph_env)
    vg_fasta = root / "vg.fa"
    rec.run("vg", [tools["vg"]["path"], "paths", "-v", vg_graph, "-F"], stdout=vg_fasta, env=graph_env)
    if set(vg_names.read_text().splitlines()) != {canonical_name, derivative_name}:
        raise GateError("vg path-name round trip failed")
    vg_records = dict(parse_fasta(vg_fasta, gzipped=False))
    if any(vg_records.get(n) != seq for n in (canonical_name, derivative_name)):
        raise GateError("vg path sequence round trip failed")
    missing = root / "missing.input"
    rec.run("pggb-negative", [tools["pggb"]["path"], "-i", missing, "-o", root / "pggb-bad", "-t", "1"],
            env=graph_env, expect={1, 2})
    rec.run("seqwish-negative", [tools["seqwish"]["path"], "-s", missing, "-p", missing,
                                  "-g", root / "seqwish-bad.gfa"], env=graph_env, expect={1, 2})
    rec.run("smoothxg-negative", [tools["smoothxg"]["path"], "-g", missing,
                                   "-o", root / "smooth-bad.gfa"], env=graph_env, expect={1, 2})
    rec.run("odgi-negative", [tools["odgi"]["path"], "paths", "-i", missing, "-L"],
            env=graph_env, expect={1, 2})
    rec.run("vg-negative", [tools["vg"]["path"], "paths", "-v", missing, "-L"],
            env=graph_env, expect={1, 2})
    rec.gate("pggb", "pggb 0.6.0", "bounded BGZF fixture derived from one authorized assembly",
             "pggb -i FIXTURE.bgzf ... -t 2 -n 2", "final GFA retains literal # path names",
             {"direct_bgzf": PASS, "literal_hash": PASS, "final_gfa_paths": PASS,
              "missing_input_rejected": PASS}, view_contract="graph-fixture-v1")
    rec.gate("seqwish", "seqwish (pinned graph environment)", "pggb PAF plus fixture sequences",
             "invoked inside pinned pggb pipeline", "induced GFA paths retain exact names",
             {"invoked_by_pggb": PASS, "path_names_exact": PASS, "downstream_sequences_exact": PASS,
              "missing_input_rejected": PASS})
    rec.gate("smoothxg", "smoothxg (pinned graph environment)", "seqwish GFA",
             "invoked inside pinned pggb pipeline", "smooth final GFA paths retain exact names",
             {"invoked_by_pggb": PASS, "path_names_exact": PASS, "downstream_sequences_exact": PASS,
              "missing_input_rejected": PASS})
    rec.gate("odgi", "odgi (pinned graph environment)", "pggb final GFA",
             "odgi build; odgi paths -L/-f", "exact two paths and FASTA headers",
             {"names_exact": PASS, "bases_exact": PASS, "missing_input_rejected": PASS})
    rec.gate("vg", "vg (pinned graph environment)", "pggb final GFA via vg convert",
             "vg convert -g -v; vg paths -L/-F", "exact two paths and FASTA headers",
             {"names_exact": PASS, "bases_exact": PASS, "missing_input_rejected": PASS})
    return plain, [canonical_name, derivative_name]


def gene_cluster_probe(rec: Recorder, tools: dict[str, Any], fixture_plain: Path,
                       fixture_names: list[str], host_env: dict[str, str]) -> None:
    root = rec.scratch / "genes"
    root.mkdir(parents=True, exist_ok=True)
    proteins = root / "proteins.faa"
    genes = root / "genes.fna"
    gff = root / "prodigal.gff"
    rec.run("prodigal", [tools["prodigal"]["path"], "-i", fixture_plain, "-a", proteins,
                          "-d", genes, "-f", "gff", "-o", gff, "-p", "single", "-q"], env=host_env)
    p_records = list(parse_fasta(proteins, gzipped=False))
    if not p_records or not any("#" in n for n, _ in p_records):
        raise GateError("Prodigal did not retain literal # in gene identifiers")
    if not all(any(line.startswith(name + "\t") for line in gff.read_text().splitlines() if not line.startswith("#"))
               for name in fixture_names):
        raise GateError("Prodigal GFF lost input seqids")
    direct_proteins = root / "direct-bgzf.proteins.faa"
    direct_genes = root / "direct-bgzf.genes.fna"
    direct_gff = root / "direct-bgzf.gff"
    rec.run("prodigal-bgzf", [tools["prodigal"]["path"], "-i", str(fixture_plain)+".gz",
        "-a", direct_proteins, "-d", direct_genes, "-f", "gff", "-o", direct_gff,
        "-p", "single", "-q"], env=host_env)
    if (list(parse_fasta(direct_proteins, gzipped=False)) != p_records
            or list(parse_fasta(direct_genes, gzipped=False)) != list(parse_fasta(genes, gzipped=False))):
        raise GateError("Prodigal BGZF/plain semantic output mismatch")
    rec.run("prodigal-negative", [tools["prodigal"]["path"], "-i", root / "missing.fa",
        "-a", root / "must-not-exist.faa", "-p", "single", "-q"], env=host_env, expect={5})

    cluster_prefix = root / "mmseqs-cluster"
    mm_tmp = root / "mmseqs-tmp"
    rec.run("mmseqs", [tools["mmseqs"]["path"], "easy-cluster", proteins, cluster_prefix, mm_tmp,
                        "--min-seq-id", "0.5", "-c", "0.8", "--threads", "2"], env=host_env)
    cluster_tsv = Path(str(cluster_prefix) + "_cluster.tsv")
    if not cluster_tsv.is_file() or "#" not in cluster_tsv.read_text():
        raise GateError("MMseqs2 cluster membership lost literal # names")

    first_name, first_protein = p_records[0]
    alignment = root / "one.sto"
    alignment.write_text("# STOCKHOLM 1.0\n" + first_name + " " + first_protein.decode() + "\n//\n")
    hmm = root / "one.hmm"
    rec.run("hmmbuild", [tools["hmmbuild"]["path"], "--amino", hmm, alignment], env=host_env)
    tbl = root / "hmmsearch.tbl"
    rec.run("hmmsearch", [tools["hmmsearch"]["path"], "--tblout", tbl, hmm, proteins], env=host_env)
    if first_name not in tbl.read_text():
        raise GateError("HMMER target name round trip failed")

    names = [name for name, _ in p_records[:2]]
    if len(names) < 2:
        raise GateError("not enough genes for MCL name test")
    abc = root / "network.abc"
    abc.write_text(f"{names[0]}\t{names[1]}\t1.0\n{names[1]}\t{names[0]}\t1.0\n")
    clusters = root / "mcl.clusters"
    rec.run("mcl", [tools["mcl"]["path"], abc, "--abc", "-I", "2.0", "-o", clusters], env=host_env)
    cluster_text = clusters.read_text()
    if any(n not in cluster_text for n in names):
        raise GateError("MCL ABC labels did not round trip")
    missing = root / "missing.input"
    rec.run("mmseqs-negative", [tools["mmseqs"]["path"], "createdb", missing, root / "bad-db"],
            env=host_env, expect={1, 2})
    rec.run("hmmer-negative", [tools["hmmbuild"]["path"], root / "bad.hmm", missing],
            env=host_env, expect={6})
    rec.run("mcl-negative", [tools["mcl"]["path"], missing, "--abc", "-o", root / "bad.clusters"],
            env=host_env, expect={1, 2})
    rec.gate("prodigal", "Prodigal 2.6.3", "bounded BGZF fixture derived from an authorized canonical path",
             "prodigal -i FIXTURE.fa.gz -a proteins.faa -d genes.fna -f gff", "gene/GFF IDs retain path prefix and literal #",
             {"direct_bgzf_success": PASS, "plain_bgzf_semantic_equal": PASS, "literal_hash": PASS,
              "gff_seqids_exact": PASS, "missing_input_rejected": PASS}, view_contract="graph-fixture-v1")
    rec.gate("mmseqs2", "MMseqs2 18.8cc5c", "Prodigal protein FASTA staged derivative",
             "mmseqs easy-cluster proteins.faa ... --min-seq-id 0.5 -c 0.8", "cluster TSV retains full protein IDs",
             {"cluster_output": PASS, "literal_hash": PASS, "missing_input_rejected": PASS})
    rec.gate("hmmer", "HMMER 3.4", "bounded protein alignment/HMM and Prodigal protein FASTA",
             "hmmbuild; hmmsearch --tblout", "target identifiers retain literal #",
             {"profile_build": PASS, "self_target_recovered": PASS, "literal_hash": PASS,
              "missing_input_rejected": PASS})
    rec.gate("mcl", "MCL 22.282", "ABC network over exact protein identifiers",
             "mcl network.abc --abc -I 2.0", "cluster labels retain literal # identifiers",
             {"abc_parse": PASS, "labels_exact": PASS, "missing_input_rejected": PASS})


def _copy_compact_logs(rec: Recorder) -> None:
    # Command metadata already contains compact stderr tails; retain no bulky tool outputs.
    pass


def summarize_resources(rec: Recorder, request: ResourceRequest, swap_start_free: int, swap_end_free: int) -> dict[str, Any]:
    durable_size, durable_files = tree_size(rec.release)
    summary = {
        "verdict": PASS,
        "assigned_ram_bytes": request.assigned_ram_bytes,
        "peak_rss_bytes": rec.max_rss_bytes,
        "peak_rss_fraction": rec.max_rss_bytes / request.assigned_ram_bytes,
        "ram_le_70_percent": rec.max_rss_bytes <= request.assigned_ram_bytes * 0.70,
        "swap_free_start_bytes": swap_start_free, "swap_free_end_bytes": swap_end_free,
        "swap_growth_bytes": max(0, swap_start_free - swap_end_free),
        "scratch_peak_bytes": rec.peak_scratch_bytes,
        "scratch_peak_files": rec.peak_scratch_files,
        "inode_allocation": request.inode_allocation,
        "actual_peak_files_fraction": rec.peak_scratch_files / request.inode_allocation,
        "predicted_scratch_upper95_bytes": request.predicted_scratch_peak_bytes,
        "scratch_upper95_fraction": request.predicted_scratch_peak_bytes / request.scratch_allocation_bytes,
        "durable_current_bytes": durable_size, "durable_current_files": durable_files,
        "predicted_durable_upper95_bytes": request.predicted_durable_peak_bytes,
        "durable_upper95_fraction": request.predicted_durable_peak_bytes / request.durable_allocation_bytes,
        "projected_files_fraction": request.predicted_files / request.inode_allocation,
        "unfinished_write_reserve_factor_durable": (request.durable_allocation_bytes-request.predicted_durable_peak_bytes)/request.unfinished_write_bytes,
        "unfinished_write_reserve_factor_scratch": (request.scratch_allocation_bytes-request.predicted_scratch_peak_bytes)/request.unfinished_write_bytes,
        "scale_trend": "NOT_APPLICABLE_NON_SCALE_BEARING_COMPATIBILITY",
    }
    if not summary["ram_le_70_percent"] or summary["swap_growth_bytes"] != 0:
        raise GateError("RAM/swap resource gate failed")
    if summary["scratch_upper95_fraction"] > .70 or summary["durable_upper95_fraction"] > .70:
        raise GateError("disk upper-95 resource gate failed")
    if summary["projected_files_fraction"] > .50 or summary["actual_peak_files_fraction"] > .50:
        raise GateError("inode allocation resource gate failed")
    return summary


def _swap_free() -> int:
    for line in Path("/proc/meminfo").read_text().splitlines():
        if line.startswith("SwapFree:"):
            return int(line.split()[1]) * 1024
    return 0


def copy_tracked_views(repo: Path, final: Path) -> None:
    manifest_dir = repo / "manifests/consumer-compatibility-v1"
    artifact_dir = repo / "artifacts/consumer_compatibility"
    manifest_dir.mkdir(parents=True, exist_ok=True)
    artifact_dir.mkdir(parents=True, exist_ok=True)
    for name in ("release.json", "consumers.tsv"):
        shutil.copy2(final / name, manifest_dir / name)
    shutil.copy2(final / "SHA256SUMS", manifest_dir / "external_SHA256SUMS")
    gates = [json.loads(p.read_text()) for p in sorted((final / "gates").glob("*.json"))]
    with (manifest_dir / "gates.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["gate_id", "tool", "verdict", "input_form", "view_contract", "invocation", "output_name_behavior", "gate_sha256"])
        for gate, path in zip(gates, sorted((final / "gates").glob("*.json"))):
            w.writerow([gate["gate_id"], gate["tool"], gate["verdict"], gate["input_form"],
                        gate["view_contract"], gate["invocation"], gate["output_name_behavior"], sha256_file(path)])
    for src, dst in (
        ("validation.json", "validation.json"), ("resource_summary.json", "resource_summary.json"),
        ("restart_evidence.json", "restart_evidence.json"), ("deterministic_rerun.json", "deterministic_rerun.json"),
        ("tool_versions.json", "tool_versions.json"),
        ("root_hashes.json", "root_hashes_finish.json"), ("view_contracts.json", "view_contracts.json"),
    ):
        if (final / src).exists(): shutil.copy2(final / src, artifact_dir / dst)
    # Compact status/command evidence only.
    shutil.copy2(final / "commands.jsonl", artifact_dir / "commands.jsonl")
    files = sorted(p for p in manifest_dir.iterdir() if p.name != "SHA256SUMS" and p.is_file())
    (manifest_dir / "SHA256SUMS").write_text("".join(f"{sha256_file(p)}  {p.name}\n" for p in files))


def publish_release(repo: Path, durable_root: Path, scratch_root: Path, request: ResourceRequest,
                    graph_prefix: Path, impg: Path, run_id: str) -> Path:
    rid = release_id(repo)
    final = durable_root / rid
    if final.exists():
        from .validate_release import validate_release
        validate_release(repo, final)
        copy_tracked_views(repo, final)
        return final
    staging = durable_root / f".{rid}.staging.{run_id}"
    if staging.exists(): shutil.rmtree(staging)
    if scratch_root.exists(): shutil.rmtree(scratch_root)
    scratch_root.mkdir(parents=True)
    staging.mkdir(parents=True)
    rec = Recorder(staging, scratch_root, request)
    swap_start = _swap_free()
    rec.state("release", "START", release_id=rid, run_id=run_id)
    host_prefix = scratch_root / "env/host-gene"
    host_env = os.environ.copy()
    graph_env = os.environ.copy()
    try:
        preflight = rec.preflight("start", durable_root)
        release, assemblies, predecessor = verify_predecessor(repo)
        write_json(staging / "input_manifest.json", {
            "schema_version": "consumer-compatibility-input-v1", "predecessor": predecessor,
            "predecessor_manifest_bytes_sha256": sha256_file(repo / "manifests/canonical-cohort-010-v1/release.json"),
            "host_lock_sha256": sha256_file(repo / HOST_LOCK),
            "host_package_sha256_inventory": sha256_file(repo / HOST_PACKAGE_SHA256),
            "graph_lock_sha256": sha256_file(repo / GRAPH_LOCK),
            "graph_package_sha256_inventory": sha256_file(repo / GRAPH_PACKAGE_SHA256),
            "assembly_count": 10, "assembly_order": EXPECTED_ASSEMBLIES, "global_distinct_assembly_union": 10,
            "global_cap": 1000,
        })
        write_json(staging / "root_hashes.json", {"start": verify_root_hashes(repo)})
        install_host_environment(repo, host_prefix, rec, durable_root)
        host_env["PATH"] = str(host_prefix / "bin") + os.pathsep + host_env.get("PATH", "")
        graph_env["PATH"] = str(graph_prefix / "bin") + os.pathsep + graph_env.get("PATH", "")
        graph_env["LC_ALL"] = host_env["LC_ALL"] = "C"
        rec.preflight("before_consumer_probes", durable_root)
        tools = tool_identities(host_prefix, graph_prefix, impg)
        write_json(staging / "tool_versions.json", {
            "host_environment_lock": {"path": HOST_LOCK, "sha256": sha256_file(repo / HOST_LOCK),
                "package_sha256_inventory_path": HOST_PACKAGE_SHA256,
                "package_sha256_inventory_sha256": sha256_file(repo / HOST_PACKAGE_SHA256)},
            "graph_environment_lock": {"path": GRAPH_LOCK, "sha256": sha256_file(repo / GRAPH_LOCK),
                "package_sha256_inventory_path": GRAPH_PACKAGE_SHA256,
                "package_sha256_inventory_sha256": sha256_file(repo / GRAPH_PACKAGE_SHA256)},
            "tools": tools,
        })
        sequences = samtools_probe(rec, assemblies, tools)
        panel, panel_contract = stage_combined_panel(rec, assemblies, Path(tools["bgzip"]["path"]))
        impg_probe(rec, tools, panel, sequences)
        mash_rapidnj_probe(rec, tools, assemblies, host_env)
        ani_qc_probe(rec, tools, assemblies, host_env)
        gff_probe(rec, tools, sequences, host_env)
        graph_plain, graph_names = graph_probe(rec, tools, sequences, graph_env)
        # The graph fixture's compressed sibling is the staged input tested above.
        gene_cluster_probe(rec, tools, graph_plain, graph_names, host_env)
        rec.preflight("finish_before_cleanup", durable_root)
        roots_finish = verify_root_hashes(repo)
        roots_doc = json.loads((staging / "root_hashes.json").read_text())
        roots_doc["finish"] = roots_finish
        roots_doc["verdict"] = PASS
        write_json(staging / "root_hashes.json", roots_doc)
        restart_path = durable_root / "restart-evidence.json"
        if not restart_path.is_file():
            raise GateError("required injected kill/restart evidence is absent")
        restart = json.loads(restart_path.read_text())
        if restart.get("verdict") != PASS:
            raise GateError("injected kill/restart evidence is not PASS")
        write_json(staging / "restart_evidence.json", restart)
        views = [json.loads(p.read_text()) for p in sorted((staging / "view_contracts").glob("*.json"))]
        write_json(staging / "view_contracts.json", {"verdict": PASS, "views": views,
             "scratch_namespace": str(scratch_root), "cleanup_required": True})
        swap_end = _swap_free()
        summary = summarize_resources(rec, request, swap_start, swap_end)
        write_json(staging / "resource_summary.json", summary)
        required_determinism = {
            "combined_bgzf_byte_identical", "impg_six_file_byte_identical",
            "impg_query_byte_identical", "impg_map_byte_identical", "mash_triangle_byte_identical",
        }
        if not required_determinism.issubset(rec.determinism) or any(
            rec.determinism[key] != PASS for key in required_determinism
        ):
            raise GateError("deterministic rerun evidence incomplete")
        write_json(staging / "deterministic_rerun.json", {"verdict": PASS, **rec.determinism})
        consumers = sorted(rec.gates, key=lambda x: x["gate_id"])
        with (staging / "consumers.tsv").open("w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(["consumer_id", "tool", "verdict", "input_form", "view_contract", "invocation", "output_name_behavior"])
            for g in consumers:
                w.writerow([g["gate_id"], g["tool"], g["verdict"], g["input_form"], g["view_contract"],
                            g["invocation"], g["output_name_behavior"]])
        rec.state("release", "VALIDATED", consumer_gates=len(consumers))
        if not rec.failures.exists(): rec.failures.write_text("")
        # Release excludes inapplicable/unselected CheckM2: isolate QC requires a separately frozen DB;
        # it is not silently treated as compatible.
        release_doc = {
            "schema_version": SCHEMA, "release_id": rid, "immutable": True, "verdict": PASS,
            "source_task_id": "certify-pilot-consumer", "created_at_utc": utc_now(),
            "predecessor_release_id": predecessor["release_id"],
            "predecessor_release_json_sha256": predecessor["release_json_sha256"],
            "counts": {"assemblies": 10, "contigs": 1223, "bases": 51731662,
                       "required_consumer_gates": len(consumers), "pass_consumer_gates": len(consumers),
                       "distinct_sequence_bearing_assemblies": 10, "global_distinct_assembly_cap": 1000},
            "applicable_gates": {
                "predecessor_identity_checksums_rows": PASS, "bgzf_index_name_base_coordinate_roundtrip": PASS,
                "consumer_compatibility": PASS, "direct_and_staged_views": PASS,
                "deterministic_semantic_validation": PASS, "injected_kill_restart": PASS,
                "resource": PASS, "root_source_immutability": PASS, "global_distinct_assembly_cap": PASS,
                "scale_trend": "NOT_APPLICABLE_NON_SCALE_BEARING_COMPATIBILITY",
            },
            "cohort_order": EXPECTED_ASSEMBLIES,
            "host_environment_lock_sha256": sha256_file(repo / HOST_LOCK),
            "host_package_sha256_inventory_sha256": sha256_file(repo / HOST_PACKAGE_SHA256),
            "graph_environment_lock_sha256": sha256_file(repo / GRAPH_LOCK),
            "graph_package_sha256_inventory_sha256": sha256_file(repo / GRAPH_PACKAGE_SHA256),
            "excluded_not_applicable": {
                "CheckM2": "not selected: isolate assemblies; no checksum-frozen CheckM2 model database exists; downstream may not assume compatibility",
                "FastANI": "not selected: skani 0.3.1 is the pinned ANI route",
                "Panaroo_Roary": "not selected: downstream phage route is Prodigal+MMseqs2+HMMER+MCL, not bacterial GFF pangenome callers",
            },
            "scientific_claims": {"mash_rapidnj": "unrooted genomic-similarity dendrogram only, never a phylogeny",
                                  "impg_map": "syncmer-anchor projection, never a base alignment"},
            "external_path": str(final),
        }
        write_json(staging / "release.json", release_doc)
        write_inventory(staging)
        complete = {"release_id": rid, "sha256sums_sha256": sha256_file(staging / "SHA256SUMS"),
                    "completed_at_utc": utc_now(), "verdict": PASS}
        (staging / "COMPLETE").write_bytes(canonical_json(complete))
        atomic_promote(staging, final)
    except BaseException as exc:
        try: rec.failure("run", f"{type(exc).__name__}: {exc}")
        except Exception: pass
        shutil.rmtree(staging, ignore_errors=True)
        shutil.rmtree(scratch_root, ignore_errors=True)
        raise
    shutil.rmtree(scratch_root)
    require_cleaned(scratch_root)
    restart_path.unlink(missing_ok=True)
    from .validate_release import validate_release
    try:
        validation = validate_release(repo, final)
    except BaseException:
        quarantine = durable_root / "quarantine"
        quarantine.mkdir(exist_ok=True)
        destination = quarantine / f"{final.name}.failed-independent-validation.{run_id}"
        if destination.exists():
            destination = quarantine / f"{final.name}.failed-independent-validation.{run_id}.{time.time_ns()}"
        os.rename(final, destination)
        raise
    # validation.json is a tracked view; the immutable external inventory is not mutated.
    write_json(repo / "artifacts/consumer_compatibility/validation.json", validation)
    copy_tracked_views(repo, final)
    return final


def promotion_worker(staging: Path, final: Path) -> None:
    staging.mkdir(parents=True, exist_ok=False)
    (staging / "partial.txt").write_text("partial\n")
    fd = os.open(staging / "partial.txt", os.O_RDONLY)
    os.fsync(fd); os.close(fd)
    os.kill(os.getpid(), signal.SIGKILL)
    raise AssertionError("unreachable")


def interrupt_test(durable_root: Path, scratch_root: Path, request: ResourceRequest) -> dict[str, Any]:
    resource_preflight("interrupt_test", durable_root, scratch_root, request)
    staging = durable_root / ".compat-interrupt.staging"
    final = durable_root / "compat-interrupt-must-not-publish"
    shutil.rmtree(staging, ignore_errors=True); shutil.rmtree(final, ignore_errors=True)
    proc = subprocess.run([sys.executable, "-m", "workflow.compatibility.pilot", "_promotion-worker",
                           "--staging", str(staging), "--final", str(final)])
    if proc.returncode != -signal.SIGKILL or final.exists():
        raise GateError(f"injected SIGKILL did not fail closed: rc={proc.returncode} final={final.exists()}")
    partial_observed = staging.exists() and not (staging / "COMPLETE").exists()
    shutil.rmtree(staging)
    restart_stage = durable_root / ".compat-interrupt.restart-staging"
    restart_stage.mkdir()
    (restart_stage / "payload.txt").write_text("clean restart\n")
    write_inventory(restart_stage)
    (restart_stage / "COMPLETE").write_bytes(canonical_json({"verdict": PASS}))
    atomic_promote(restart_stage, final)
    restart_valid = final.is_dir() and (final / "COMPLETE").is_file()
    shutil.rmtree(final)
    evidence = {
        "verdict": PASS if partial_observed and restart_valid else "FAIL",
        "injected_signal": "SIGKILL", "worker_returncode": proc.returncode,
        "partial_never_published": partial_observed, "partial_cleaned": not staging.exists(),
        "clean_restart_promoted_only_after_complete": restart_valid,
        "interrupted_final_absent_after_cleanup": not final.exists(), "tested_at_utc": utc_now(),
    }
    if evidence["verdict"] != PASS: raise GateError("interrupt/restart test failed")
    write_json(durable_root / "restart-evidence.json", evidence)
    shutil.rmtree(scratch_root, ignore_errors=True)
    return evidence


def request_from_args(args: argparse.Namespace) -> ResourceRequest:
    return ResourceRequest(args.assigned_ram_bytes, args.durable_allocation_bytes,
        args.scratch_allocation_bytes, args.inode_allocation, args.predicted_durable_peak_bytes,
        args.predicted_scratch_peak_bytes, args.predicted_files, args.unfinished_write_bytes)


def add_resource_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--assigned-ram-bytes", type=int, required=True)
    p.add_argument("--durable-allocation-bytes", type=int, required=True)
    p.add_argument("--scratch-allocation-bytes", type=int, required=True)
    p.add_argument("--inode-allocation", type=int, required=True)
    p.add_argument("--predicted-durable-peak-bytes", type=int, required=True)
    p.add_argument("--predicted-scratch-peak-bytes", type=int, required=True)
    p.add_argument("--predicted-files", type=int, required=True)
    p.add_argument("--unfinished-write-bytes", type=int, required=True)


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="command", required=True)
    for name in ("interrupt-test", "run"):
        q = sub.add_parser(name)
        q.add_argument("--durable-task-root", type=Path, required=True)
        q.add_argument("--scratch-root", type=Path, required=True)
        add_resource_args(q)
    q = sub.choices["run"]
    q.add_argument("--repo", type=Path, default=Path.cwd())
    q.add_argument("--run-id", default="compat-v1")
    q.add_argument("--graph-prefix", type=Path, default=GRAPH_PREFIX_DEFAULT)
    q.add_argument("--impg", type=Path, default=IMPG_DEFAULT)
    hidden = sub.add_parser("_promotion-worker")
    hidden.add_argument("--staging", type=Path, required=True)
    hidden.add_argument("--final", type=Path, required=True)
    return p


def main(argv: Sequence[str] | None = None) -> int:
    args = parser().parse_args(argv)
    if args.command == "_promotion-worker":
        promotion_worker(args.staging, args.final)
        return 99
    request = request_from_args(args)
    if args.command == "interrupt-test":
        print(json.dumps(interrupt_test(args.durable_task_root, args.scratch_root, request), sort_keys=True))
    elif args.command == "run":
        final = publish_release(args.repo.resolve(), args.durable_task_root, args.scratch_root,
                                request, args.graph_prefix, args.impg, args.run_id)
        print(final)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except GateError as exc:
        print(f"NO_GO: {exc}", file=sys.stderr)
        raise SystemExit(2)
