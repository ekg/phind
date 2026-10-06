#!/usr/bin/env python3
"""
build_v3_release.py — assemble the unified NTM v3 release (task ntm-v3-release).

Joins the v3 ML/ancestral phage genomes + per-clade ``release_manifest.tsv``
(task ntm-v3-ml) to the v3 host clades (task ntm-v3-host) and emits the
deliverable tables + scope-labeled FASTAs for the **single unified
Mycobacteriaceae catalog** (binding SCOPE DECISION, ``ntm/v3/PLAYBOOK.md``:
"one catalog spanning the whole family. Tag, do not filter.").

Inputs (defaults; all overridable):
  $NVME/ntm/v3/ml/all_ntm_v3_ml_phage_genomes.fa        1304 ML genomes
  $NVME/ntm/v3/ml/all_ntm_v3_ancestral_phage_genomes.fa 893 ancestral genomes
  $NVME/ntm/v3/ml/release_manifest.tsv                  per-clade ML stats
  $NVME/ntm/v3/host_clades/host_clades.tsv              accession -> host clade/species
  $NVME/ntm/v3/clades/clade_summary.tsv                 per-clade median MASH

Join rule (prophage member -> host clade):
  a member id is ``<canonical_acc>#1#<prophage_id>``; the canonical key is the
  first ``#`` token.  ``canonical_acc -> host_clades.tsv.accession`` gives the
  host species and host clade.  (The v3 ML manifest ``source_prophages`` column
  carries the canonical-keyed member list for every clade, so no separate clade
  membership file is needed.)

Scope rule (BINDING):
  The per-genome MTC classifier is a **species-prefix** test: a host genome is
  ``MTC`` when its ``host_clades.tsv.species`` **starts with**
  ``Mycobacterium tuberculosis`` / ``canetti`` / ``orygis`` / ``africanum`` /
  ``bovis`` / ``microti``, else ``NTM``.  ``startswith`` excludes substring false
  positives such as ``Mycobacterium avium subsp. paratuberculosis``.

  ``host_scope`` is then derived **per host clade**: a MASH-defined host clade is
  ``MTC`` when the majority of its genomes carry an MTC species prefix, else
  ``NTM``; every prophage member inherits its host clade's scope.  (A host clade
  is a lineage, so its scope is the defensible unit; ``Mycobacterium sp.``
  genomes that ANI-cluster inside the M. tuberculosis clade ``host_clade_0002``
  are thus correctly MTC.)  Per-genome ``host_scope`` ∈ {NTM, MTC} — binary.
  ``host_scope`` is a LABEL, never a filter — no population is dropped.

Clade-level (per prophage clade) scope label ∈ {NTM, MTC, MIXED}:
  ``MTC``   every resolved member inherits an MTC host clade
  ``NTM``   every resolved member inherits an NTM host clade
  ``MIXED`` members span both scopes (cross-boundary clade)

Outputs (default $NVME/ntm/v3/release/):
  all_ntm_v3_ml_phage_genomes.fa         one ML genome per clade, header carries
                                         host_clade_ids + host_scope + species
  all_ntm_v3_ancestral_phage_genomes.fa  one ancestral genome per alignable clade
  release_manifest.tsv                   input manifest + host_scope,
                                         host_clade_ids, n_host_clades,
                                         source_breakdown
  per_scope_summary.tsv                  NTM / MTC / MIXED roll-up
  cross_boundary_report.tsv              MIXED clades with per-clade confidence
                                         from median MASH + member counts
  release_scope_statement.md             scoping statement for the release notes
  host_clade_scope.tsv                   per host clade: derived scope, genome
                                         counts, majority fraction (audit of
                                         the clade-level derivation)

Usage:
  python3 ntm/v3/scripts/build_v3_release.py \
      --ml-dir  /mnt/nvme3n1/erikg/phind-genome-work/ntm/v3/ml \
      --host-clades /mnt/nvme3n1/erikg/phind-genome-work/ntm/v3/host_clades/host_clades.tsv \
      --out /mnt/nvme3n1/erikg/phind-genome-work/ntm/v3/release
  python3 ntm/v3/scripts/build_v3_release.py --dry-run
  python3 ntm/v3/scripts/build_v3_release.py --clade-ids 0_0000,0_0306 --out /tmp/rel-test
"""
import argparse
import csv
import datetime as dt
import os
import sys
from collections import defaultdict

V3 = "/mnt/nvme3n1/erikg/phind-genome-work/ntm/v3"
ML_DIR = f"{V3}/ml"
HOST_CLADES = f"{V3}/host_clades/host_clades.tsv"
CLADE_SUMMARY = f"{V3}/clades/clade_summary.tsv"
OUT_DIR = f"{V3}/release"

# MTC scope prefixes (BINDING rule; host_clades.tsv.species is the source).
MTC_PREFIXES = (
    "Mycobacterium tuberculosis",
    "Mycobacterium canetti",
    "Mycobacterium orygis",
    "Mycobacterium africanum",
    "Mycobacterium bovis",
    "Mycobacterium microti",
)

SCOPE_NTM = "NTM"
SCOPE_MTC = "MTC"
SCOPE_MIXED = "MIXED"

MANIFEST_NEW_COLS = ["host_scope", "host_clade_ids", "n_host_clades",
                     "source_breakdown"]


def scope_of(species):
    """MTC if the host species starts with an MTC prefix, else NTM."""
    return SCOPE_MTC if species.startswith(MTC_PREFIXES) else SCOPE_NTM


def member_accession(member_id):
    """<canonical_acc>#1#<prophage_id> -> <canonical_acc>."""
    return member_id.split("#", 1)[0]


def ml_file_for(clade_id):
    return f"ntm3_{clade_id}_ML"


def anc_file_for(clade_id):
    return f"ntm3_{clade_id}_ANCESTRAL"


def read_all_fasta(path):
    """Return {first_header_token: sequence} for a FASTA file."""
    out = {}
    cur, seq = None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if cur is not None:
                    out[cur] = "".join(seq)
                cur = line[1:].split()[0]
                seq = []
            else:
                seq.append(line.strip())
    if cur is not None:
        out[cur] = "".join(seq)
    return out


def load_host_clades(path):
    """Return ({accession: (host_clade_id, species)}, {host_clade_id: scope},
    {host_clade_id: (n_mtc_genomes, n_genomes)}).

    The host-clade scope is derived by majority vote using the per-genome
    species-prefix classifier ``scope_of``: a clade is MTC when strictly more
    than half of its genomes carry an MTC prefix, else NTM.  Every prophage
    member inherits its host clade's scope.
    """
    host = {}
    n_mtc = defaultdict(int)
    n_total = defaultdict(int)
    with open(path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            hc, sp = row["host_clade_id"], row["species"]
            host[row["accession"]] = (hc, sp)
            n_total[hc] += 1
            if scope_of(sp) == SCOPE_MTC:
                n_mtc[hc] += 1
    hc_scope = {hc: (SCOPE_MTC if n_mtc[hc] * 2 > n_total[hc] else SCOPE_NTM)
                for hc in n_total}
    hc_counts = {hc: (n_mtc[hc], n_total[hc]) for hc in n_total}
    return host, hc_scope, hc_counts


def load_clade_mash(path):
    """clade_id -> median_mash (float) or None."""
    mash = {}
    with open(path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            try:
                mash[row["clade_id"]] = float(row["median_mash"])
            except (KeyError, ValueError):
                mash[row["clade_id"]] = None
    return mash


def cross_boundary_confidence(median_mash, n_members, minority_n):
    """Ordinal confidence for a MIXED clade.

    A single minority-scope member is always ``low`` (a lone cross-scope
    prophage can be a host-assignment artefact).  Otherwise confidence rises
    with clade cohesion (low median MASH) and support (member count).
    """
    if minority_n <= 1:
        return "low"
    if median_mash is not None and median_mash <= 0.10 and n_members >= 10:
        return "high"
    if median_mash is not None and median_mash <= 0.25 and n_members >= 5:
        return "medium"
    return "low"


def write_scope_statement(path, stamp, hc_scope, hc_counts, scope_clades,
                          scope_members, scope_host_clades, n_ntm_members,
                          n_mtc_members, alt_scope_clades, clade_derived_mtc,
                          out_manifest, anc_seqs, n_clades, is_subset):
    """Write release_scope_statement.md documenting the scope derivation."""
    mtc_host_clades = sorted(h for h in hc_counts if hc_scope[h] == SCOPE_MTC)
    derived_mtc = sum(clade_derived_mtc.values())
    ref = "host_clade_0002"
    if ref in hc_counts:
        nm, nt = hc_counts[ref]
        ref_txt = f"{ref} ({nm}/{nt} = {100.0 * nm / nt:.1f}% MTC)"
    else:
        ref_txt = ref
    lines = []
    lines.append("# NTM v3 release - scoping statement")
    lines.append("")
    lines.append(f"Generated: {stamp} (unified Mycobacteriaceae catalog)")
    lines.append("")
    lines.append("`host_scope` is a **label**, not a filter: no population is "
                 "dropped. Per the binding scope decision (`ntm/v3/PLAYBOOK.md`), "
                 "the release is one catalog spanning the whole family.")
    lines.append("")
    lines.append("## Label namespaces")
    lines.append("")
    lines.append("- per genome (host) `host_scope` in {NTM, MTC} - binary only.")
    lines.append("- per prophage clade label in {NTM, MTC, MIXED} "
                 "(MIXED = cross-boundary; see `cross_boundary_report.tsv`).")
    lines.append("")
    lines.append("## Clade-level scope derivation (binding)")
    lines.append("")
    lines.append("A MASH-defined host clade is MTC when the majority of its "
                 "genomes carry an MTC species prefix (tuberculosis / canetti / "
                 "orygis / africanum / bovis / microti, via `startswith`), else "
                 "NTM. Every prophage member inherits its host clade's scope; "
                 "the per-prophage-clade label is MTC/NTM when all members "
                 "inherit that scope, and MIXED when both occur.")
    lines.append("")
    lines.append("The per-genome MTC species-prefix classifier is kept and used "
                 "only to decide each host clade's majority. Rationale: a host "
                 "clade is an ANI-defined lineage, so scope is a property of "
                 "the lineage rather than of an individual unresolved species "
                 f"label. Members labeled `Mycobacterium sp.` that sit in {ref_txt} "
                 "are unresolved-ANI MTC genomes, so the species-prefix rule "
                 "under-calls them; the clade-level rule assigns them MTC.")
    lines.append("")
    lines.append(
        "Per-genome (host) scope is auditable in `host_clade_scope.tsv`.")
    lines.append("")
    lines.append(f"MTC host clades (majority MTC): {len(mtc_host_clades)} "
                 f"({', '.join(mtc_host_clades) if mtc_host_clades else 'none'}).")
    lines.append("")
    lines.append("### Clade-derived MTC assignments (unresolved species labels)")
    lines.append("")
    lines.append(f"{derived_mtc} prophage members carry a non-MTC species label "
                 "but sit in an MTC host clade, so they are MTC under the "
                 "clade-level rule:")
    lines.append("")
    if clade_derived_mtc:
        lines.append("| prophage clade | clade-derived MTC members |")
        lines.append("|---|---:|")
        for pc, n in sorted(clade_derived_mtc.items()):
            lines.append(f"| {pc} | {n} |")
    else:
        lines.append("_(none)_")
    lines.append("")
    lines.append("### Clade-level scope vs species-prefix accounting")
    lines.append("")
    lines.append("| rule | MTC members | wholly-MTC clades | MIXED clades | fully-NTM clades |")
    lines.append("|---|---:|---:|---:|---:|")
    lines.append(f"| clade-level scope (released) | {n_mtc_members} | "
                 f"{len(scope_clades[SCOPE_MTC])} | {len(scope_clades[SCOPE_MIXED])} | "
                 f"{len(scope_clades[SCOPE_NTM])} |")
    lines.append(f"| species-prefix only | {n_mtc_members - derived_mtc} | "
                 f"{alt_scope_clades[SCOPE_MTC]} | {alt_scope_clades[SCOPE_MIXED]} | "
                 f"{alt_scope_clades[SCOPE_NTM]} |")
    lines.append("")
    lines.append(f"Members by clade-derived member-scope: NTM={n_ntm_members}, "
                 f"MTC={n_mtc_members} (of {sum(scope_members.values())}).")
    lines.append("")
    lines.append("| host_scope | clades | ML genomes | ancestral genomes | members | host clades |")
    lines.append("|---|---:|---:|---:|---:|---:|")
    for sc in (SCOPE_NTM, SCOPE_MTC, SCOPE_MIXED, "UNJOINED"):
        cids = scope_clades[sc]
        if not cids:
            continue
        n_anc_sc = sum(1 for r in out_manifest if r["host_scope"] == sc
                       and anc_file_for(r["clade_id"]) in anc_seqs)
        lines.append(f"| {sc} | {len(cids)} | {len(cids)} | {n_anc_sc} | "
                     f"{scope_members[sc]} | {len(scope_host_clades[sc])} |")
    lines.append("")
    lines.append(f"_Subset run: {n_clades} clades "
                 f"({'subset' if is_subset else 'full set'})._")
    lines.append("")
    with open(path, "w") as f:
        f.write("\n".join(lines))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ml-dir", default=ML_DIR)
    ap.add_argument("--host-clades", default=HOST_CLADES)
    ap.add_argument("--clade-summary", default=CLADE_SUMMARY)
    ap.add_argument("--out", default=OUT_DIR)
    ap.add_argument("--clade-ids", default="",
                    help="comma-separated clade ids to restrict to (debug/subset)")
    ap.add_argument("--limit", type=int, default=0,
                    help="restrict to the first N clade ids (debug/subset)")
    ap.add_argument("--no-ancestral", action="store_true",
                    help="skip writing the ancestral FASTA")
    ap.add_argument("--dry-run", action="store_true",
                    help="compute + report counts, write no files")
    args = ap.parse_args()

    manifest_path = os.path.join(args.ml_dir, "release_manifest.tsv")
    ml_fa = os.path.join(args.ml_dir, "all_ntm_v3_ml_phage_genomes.fa")
    anc_fa = os.path.join(args.ml_dir, "all_ntm_v3_ancestral_phage_genomes.fa")

    host, hc_scope, hc_counts = load_host_clades(args.host_clades)
    mash = load_clade_mash(args.clade_summary)

    with open(manifest_path) as f:
        manifest_rows = list(csv.DictReader(f, delimiter="\t"))
    manifest_cols = list(manifest_rows[0].keys())

    selected = None
    if args.clade_ids:
        selected = [c for c in args.clade_ids.split(",") if c]
    elif args.limit:
        selected = [r["clade_id"] for r in manifest_rows[:args.limit]]
    if selected is not None:
        sel = set(selected)
        manifest_rows = [r for r in manifest_rows if r["clade_id"] in sel]

    ml_seqs = read_all_fasta(ml_fa)
    anc_seqs = read_all_fasta(anc_fa) if not args.no_ancestral else {}

    out_manifest = []
    out_ml, out_anc = [], []
    scope_clades = defaultdict(list)      # scope -> [clade_id]
    scope_members = defaultdict(int)      # scope -> member count (by clade scope)
    scope_host_clades = defaultdict(set)  # scope -> host clade ids
    n_ntm_members = n_mtc_members = 0      # members by member scope (clade-derived)
    alt_scope_clades = defaultdict(int)    # species-prefix-only clade labels (report)
    clade_derived_mtc = defaultdict(int)   # prophage clade -> members reassigned MTC
    cross_rows = []
    unjoined = []
    n_ml = n_anc = 0

    for row in manifest_rows:
        cid = row["clade_id"]
        members = [m for m in row["source_prophages"].split(";") if m]
        hc_ids, hc_species = set(), set()
        n_ntm = n_mtc = 0        # clade-derived member scopes (binding)
        n_ntm_alt = n_mtc_alt = 0  # species-prefix-only member scopes (report)
        mtc_species = set()
        for m in members:
            acc = member_accession(m)
            entry = host.get(acc)
            if entry is None:
                unjoined.append((cid, m, acc))
                continue
            hc, sp = entry
            hc_ids.add(hc)
            hc_species.add(sp)
            # per-genome scope = host clade scope (clade-level rule)
            if hc_scope[hc] == SCOPE_MTC:
                n_mtc += 1
                mtc_species.add(sp)
                if scope_of(sp) != SCOPE_MTC:
                    clade_derived_mtc[cid] += 1
            else:
                n_ntm += 1
            if scope_of(sp) == SCOPE_MTC:
                n_mtc_alt += 1
            else:
                n_ntm_alt += 1

        if n_mtc and n_ntm:
            clade_scope = SCOPE_MIXED
        elif n_mtc:
            clade_scope = SCOPE_MTC
        elif n_ntm:
            clade_scope = SCOPE_NTM
        else:
            clade_scope = "UNJOINED"

        if n_mtc_alt and n_ntm_alt:
            alt_scope_clades[SCOPE_MIXED] += 1
        elif n_mtc_alt:
            alt_scope_clades[SCOPE_MTC] += 1
        elif n_ntm_alt:
            alt_scope_clades[SCOPE_NTM] += 1

        scope_clades[clade_scope].append(cid)
        scope_members[clade_scope] += len(members)
        n_ntm_members += n_ntm
        n_mtc_members += n_mtc
        scope_host_clades[clade_scope].update(hc_ids)

        hc_str = ",".join(sorted(hc_ids)) if hc_ids else "NA"
        breakdown = f"NTM={n_ntm};MTC={n_mtc}"
        sp_str = ",".join(sorted(hc_species))

        orig_ml = ml_file_for(cid)
        orig_anc = anc_file_for(cid)

        # ML genome (one per clade)
        seq = ml_seqs.get(orig_ml)
        assert seq is not None, f"ML genome missing for {cid} ({orig_ml})"
        n_ml += 1
        ml_hdr = (f"{orig_ml} status={row['status']} n_members={row['n_members']} "
                  f"length={len(seq)} host_clade_ids={hc_str} "
                  f"host_scope={clade_scope} species={sp_str}")
        out_ml.append((ml_hdr, seq))

        # Ancestral genome (alignable clades only)
        if orig_anc in anc_seqs:
            anc_seq = anc_seqs[orig_anc]
            n_anc += 1
            anc_hdr = (f"{orig_anc} status={row['status']} n_members={row['n_members']} "
                       f"length={len(anc_seq)} host_clade_ids={hc_str} "
                       f"host_scope={clade_scope} species={sp_str}")
            out_anc.append((anc_hdr, anc_seq))

        new_row = dict(row)
        new_row["host_scope"] = clade_scope
        new_row["host_clade_ids"] = hc_str
        new_row["n_host_clades"] = len(hc_ids)
        new_row["source_breakdown"] = breakdown
        out_manifest.append(new_row)

        if clade_scope == SCOPE_MIXED:
            n_members = int(row["n_members"])
            median = mash.get(cid)
            minority = min(n_ntm, n_mtc)
            cross_rows.append((
                cid, n_members, n_ntm, n_mtc,
                "MTC" if n_mtc < n_ntm else "NTM", minority,
                "NA" if median is None else f"{median:.6f}",
                cross_boundary_confidence(median, n_members, minority),
                len(hc_ids),
                ",".join(sorted(mtc_species)),
            ))

    cross_rows.sort(key=lambda r: (-r[3], r[0]))

    print(f"clades={len(manifest_rows)} ml_genomes={n_ml} ancestral_genomes={n_anc}")
    print(f"members by member-scope (clade-derived): NTM={n_ntm_members} "
          f"MTC={n_mtc_members}")
    print(f"  [species-prefix-only comparison] clades: "
          f"NTM={alt_scope_clades[SCOPE_NTM]} MTC={alt_scope_clades[SCOPE_MTC]} "
          f"MIXED={alt_scope_clades[SCOPE_MIXED]}")
    for sc in (SCOPE_NTM, SCOPE_MTC, SCOPE_MIXED, "UNJOINED"):
        print(f"  host_scope={sc}: clades={len(scope_clades[sc])} "
              f"members={scope_members[sc]} host_clades={len(scope_host_clades[sc])}")
    print(f"cross-boundary (MIXED) clades: {len(cross_rows)}")
    for r in cross_rows:
        print(f"   {r[0]} n={r[1]} NTM={r[2]} MTC={r[3]} minority={r[4]}x{r[5]} "
              f"median_mash={r[6]} confidence={r[7]}")
    if unjoined:
        print(f"WARNING unjoined members: {len(unjoined)}", file=sys.stderr)

    if args.dry_run:
        print("dry-run: no files written")
        return 0

    os.makedirs(args.out, exist_ok=True)

    with open(os.path.join(args.out, "all_ntm_v3_ml_phage_genomes.fa"), "w") as f:
        for hdr, seq in out_ml:
            f.write(f">{hdr}\n{seq}\n")
    if not args.no_ancestral:
        with open(os.path.join(args.out,
                               "all_ntm_v3_ancestral_phage_genomes.fa"), "w") as f:
            for hdr, seq in out_anc:
                f.write(f">{hdr}\n{seq}\n")

    with open(os.path.join(args.out, "release_manifest.tsv"), "w") as f:
        f.write("\t".join(manifest_cols + MANIFEST_NEW_COLS) + "\n")
        for row in out_manifest:
            f.write("\t".join(str(row[c]) for c in manifest_cols + MANIFEST_NEW_COLS)
                    + "\n")

    with open(os.path.join(args.out, "per_scope_summary.tsv"), "w") as f:
        f.write("host_scope\tn_clades\tn_ml_genomes\tn_ancestral_genomes\t"
                "n_members\tn_host_clades\n")
        for sc in (SCOPE_NTM, SCOPE_MTC, SCOPE_MIXED, "UNJOINED"):
            cids = scope_clades[sc]
            if not cids:
                continue
            n_ml_sc = len(cids)
            n_anc_sc = sum(1 for r in out_manifest
                           if r["host_scope"] == sc
                           and anc_file_for(r["clade_id"]) in anc_seqs)
            f.write(f"{sc}\t{len(cids)}\t{n_ml_sc}\t{n_anc_sc}\t"
                    f"{scope_members[sc]}\t{len(scope_host_clades[sc])}\n")

    with open(os.path.join(args.out, "cross_boundary_report.tsv"), "w") as f:
        f.write("clade_id\tn_members\tn_ntm_members\tn_mtc_members\t"
                "minority_scope\tminority_members\tmedian_mash\tconfidence\t"
                "n_host_clades\tmtc_species\n")
        for r in cross_rows:
            f.write("\t".join(str(x) for x in r) + "\n")

    with open(os.path.join(args.out, "host_clade_scope.tsv"), "w") as f:
        f.write("host_clade_id\thost_scope\tn_genomes\tn_mtc_genomes\t"
                "mtc_fraction\n")
        for hc, (nm, nt) in sorted(hc_counts.items()):
            f.write(f"{hc}\t{hc_scope[hc]}\t{nt}\t{nm}\t{nm / nt:.4f}\n")

    stamp = dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    write_scope_statement(
        os.path.join(args.out, "release_scope_statement.md"), stamp,
        hc_scope, hc_counts, scope_clades, scope_members, scope_host_clades,
        n_ntm_members, n_mtc_members, alt_scope_clades, clade_derived_mtc,
        out_manifest, anc_seqs, len(manifest_rows), selected is not None)

    print(f"release -> {args.out}/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
