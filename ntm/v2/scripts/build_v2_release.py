#!/usr/bin/env python3
"""
build_v2_release.py — assemble the NTM v2 release: join the v2 ML/ancestral
phage genomes (task ntm-v2-ml) to the v2 host clades (task ntm-v2-host) and
emit the deliverable tables + enriched FASTAs.

Inputs:
  ntm/v2/ml/all_ntm2_ml_phage_genomes.fa         (2388 ML genomes)
  ntm/v2/ml/all_ntm2_ancestral_phage_genomes.fa  (1251 ancestral genomes)
  ntm/v2/ml/release_manifest.tsv                 (per-clade ML stats)
  ntm/v2/clades/0/tight_clades.json              (clade -> prophage members)
  ntm/v2/host_clades/host_clades.tsv             (accession -> host_clade_id, species)

Join rule (prophage member -> host clade):
  member id `GCA_x_y_genomic_prophageN` -> accession = first two `_` tokens
  -> host_clades.tsv row.  If the accession is absent (GCA/GCF twin cases:
  prophage extracted from the GCF copy while the host-clade cohort kept the
  GCA copy), fall back to the same-numeric twin accession in the cohort and
  record the accession in unjoined_genomes.tsv with reason `twin_joined`.
  Clades with no resolvable member at all are listed with reason
  `no_host_clade` and get host_clades=NA.

Outputs (ntm/v2/release/):
  all_ntm_ml_phage_genomes.fa        one ML genome per clade, header carries host clades
  all_ntm_ancestral_phage_genomes.fa one ancestral genome per alignable clade
  release_manifest.tsv               per-clade: members, status, length, source genomes,
                                     host clades, species, join status
  host_range.tsv                     clade -> host_clade_ids + species spanned
  per_ntm_clade_summary.tsv          host_clade -> # ML genomes, # clades
  unjoined_genomes.tsv               every member-level join exception + reason
"""
import csv
import json
import os
import sys

ROOT = "/mnt/nvme3n1/erikg/phind-genome-work/ntm/v2"
ML = f"{ROOT}/ml"
CL = f"{ROOT}/clades/0/tight_clades.json"
HOST = f"{ROOT}/host_clades/host_clades.tsv"
REL = f"{ROOT}/release"


def member_to_accession(member):
    return "_".join(member.split("_")[:2])


def load_fasta(path):
    recs, name, seq = [], None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    recs.append((name, "".join(seq)))
                name, seq = line[1:], []
            else:
                seq.append(line)
    if name is not None:
        recs.append((name, "".join(seq)))
    return recs


def main():
    os.makedirs(REL, exist_ok=True)
    tc = json.load(open(CL))

    # host lookup + numeric twin index
    host = {}
    num2acc = {}
    with open(HOST) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            host[row["accession"]] = (row["host_clade_id"], row["species"])
            num = row["accession"].split(".")[0].split("_", 1)[1]
            num2acc.setdefault(num, []).append(row["accession"])

    ml_meta = {}
    with open(f"{ML}/release_manifest.tsv") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            ml_meta[row["clade_id"]] = row

    ml_fasta = load_fasta(f"{ML}/all_ntm2_ml_phage_genomes.fa")
    anc_fasta = load_fasta(f"{ML}/all_ntm2_ancestral_phage_genomes.fa")
    anc_by_cid = {}
    for hdr, seq in anc_fasta:
        cid = hdr.split()[0].replace("ntm2_", "").replace("_ANCESTRAL", "")
        anc_by_cid[cid] = seq

    unjoined_rows = []  # member-level exceptions
    manifest_rows, host_range_rows, combined, ancestral = [], [], [], []
    per_host = {}
    stats = dict(clades=0, ml=0, singleton=0, joined=0, twin=0, unjoined_clades=0)

    for cid, members in sorted(tc.items()):
        stats["clades"] += 1
        meta = ml_meta.get(cid, {})
        status = meta.get("status", "?")
        if status == "singleton":
            stats["singleton"] += 1
        else:
            stats["ml"] += 1

        hc_ids, hc_spp = set(), set()
        for m in members:
            acc = member_to_accession(m)
            if acc in host:
                hc, sp = host[acc]
            else:
                num = acc.split(".")[0].split("_", 1)[1]
                twins = num2acc.get(num, [])
                if twins:
                    twin = sorted(twins)[0]
                    hc, sp = host[twin]
                    unjoined_rows.append(
                        (cid, m, acc, "twin_joined", f"host clade via {twin}"))
                    stats["twin"] += 1
                else:
                    hc, sp = None, None
                    unjoined_rows.append(
                        (cid, m, acc, "no_host_clade",
                         "accession absent from host_clades.tsv and no same-numeric twin"))
                    continue
            hc_ids.add(hc)
            hc_spp.add(sp)

        src_genomes = sorted({member_to_accession(m) for m in members})
        hc_str = ",".join(sorted(hc_ids)) if hc_ids else "NA"
        sp3 = ",".join(sorted(hc_spp)[:3])
        join_status = "joined" if hc_ids else "unjoined:no_host_clade"
        if hc_ids:
            stats["joined"] += 1
        else:
            stats["unjoined_clades"] += 1

        # ML genome: preserve sequence + original header fields, append host info
        seq = None
        orig_hdr = None
        for hdr, s in ml_fasta:
            fid = hdr.split()[0]
            if fid == f"ntm2_{cid}_ML":
                orig_hdr, seq = hdr, s
                break
        assert seq is not None, f"ML genome missing for {cid}"
        parts = orig_hdr.split()
        new_hdr = " ".join(
            [parts[0]] + parts[1:] + [f"host_clades={hc_str}", f"species={sp3}"])
        combined.append((new_hdr, seq))
        host_range_rows.append((cid, len(members), hc_str,
                                ",".join(sorted(hc_spp))))
        manifest_rows.append((
            cid, len(members), status, meta.get("ancestral", ""),
            meta.get("length_bp", ""), meta.get("median_member_len_bp", ""),
            len(src_genomes), hc_str, ",".join(sorted(hc_spp)),
            join_status, meta.get("flags", "")))

        if cid in anc_by_cid:
            ancestral.append((
                f"ntm2_{cid}_ANCESTRAL host_clades={hc_str} species={sp3}",
                anc_by_cid[cid]))

        for hcid in hc_ids:
            per_host.setdefault(hcid, set()).add(cid)

    with open(f"{REL}/all_ntm_ml_phage_genomes.fa", "w") as f:
        for hdr, seq in combined:
            f.write(f">{hdr}\n{seq}\n")
    with open(f"{REL}/all_ntm_ancestral_phage_genomes.fa", "w") as f:
        for hdr, seq in ancestral:
            f.write(f">{hdr}\n{seq}\n")
    with open(f"{REL}/release_manifest.tsv", "w") as f:
        f.write("clade_id\tn_members\tstatus\tancestral\tlength_bp\t"
                "median_member_len_bp\tn_source_genomes\thost_clades\tspecies\t"
                "join_status\tflags\n")
        for r in manifest_rows:
            f.write("\t".join(str(x) for x in r) + "\n")
    with open(f"{REL}/host_range.tsv", "w") as f:
        f.write("clade_id\tn_members\thost_clade_ids\tspecies\n")
        for r in host_range_rows:
            f.write("\t".join(str(x) for x in r) + "\n")
    with open(f"{REL}/per_ntm_clade_summary.tsv", "w") as f:
        f.write("host_clade_id\tn_ml_genomes\tn_clades\n")
        for hcid, cids in sorted(per_host.items(), key=lambda kv: -len(kv[1])):
            f.write(f"{hcid}\t{len(cids)}\t{len(cids)}\n")
    with open(f"{REL}/unjoined_genomes.tsv", "w") as f:
        f.write("clade_id\tprophage_id\taccession\treason\tdetail\n")
        for r in unjoined_rows:
            f.write("\t".join(r) + "\n")

    print(f"clades={stats['clades']} ml={stats['ml']} singleton={stats['singleton']}")
    print(f"joined={stats['joined']} twin_fallback_members={stats['twin']} "
          f"unjoined_clades={stats['unjoined_clades']}")
    print(f"host clades represented: {len(per_host)}")
    print(f"release -> {REL}/")


if __name__ == "__main__":
    sys.exit(main())
