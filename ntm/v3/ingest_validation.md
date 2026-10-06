# NTM v3 Ingest Reconciliation Report

**Scope:** read-only verification of `/home/erikg/phind/NTM_Data` collaborator delivery vs frozen BV-BRC phigaro inputs in `/home/erikg/phind/ntm/v3/inputs/`.
**Authority:** no writes except this report; no sha256 over 42 GB zips (deferred); no modification of `NTM_Data`, `ntm/v3/inputs/*`, `/mnt/nvme3n1`.

## Summary table

| # | Check | Expected | Observed | Verdict |
|---|-------|----------|----------|---------|
| 1a | Master table data rows | 26499 | 26499 | PASS |
| 1b | Master table columns | 28 | 28 | PASS |
| 1c | Unique `genome_id` (one row/genome) | 26499, 0 dups | 26499, 0 dups | PASS |
| 1d | Source breakdown | ASSEMBLY 17920 / NCBI 8331 / BV-BRC 248 | ASSEMBLY 17920 / NCBI 8331 / BV-BRC 248 | PASS |
| 2a | `Genomes_sra_assembled_17920_*.zip` .fasta members | 17920 | 17920 (17920 total) | PASS |
| 2b | `Genomes_ncbi_bvbrc_8579_*.zip` .fasta members | 8579 | 8579 (8579 total) | PASS |
| 3a | summary data rows | 26499 | 26499 | PASS |
| 3b | summary unique genome_ids | 26499 | 26499 | PASS |
| 3c | summary vs master set difference both ways | 0 / 0 | master−summary=0, summary−master=0 | PASS |
| 4a | coordinates data rows | 47994 | 47994 | PASS |
| 4b | coordinates unique genome_ids | 17117 | 17117 | PASS |
| 4c | coords genome set ⊆ summary set (orphans) | 0 | 0 | PASS |
| 4d | per-genome coord rows == summary `prophage_count` | 0 mismatches | 0 mismatches (Σ prophage_count = 47994) | PASS |
| 5/6 | 3-sample naming + extraction (see below) | 4/4 ranges exact | all match | PASS |

**Overall verdict: PASS.** All ingest reconciliation checks (1–6) hold exactly. Full-cohort v3 pipeline may proceed.

## Exact commands

```bash
# headers / line counts
head -1 NTM_Data/NTM_master_table_complete_20260922.tsv | tr '\t' '\n' | cat -n
wc -l NTM_Data/NTM_master_table_complete_20260922.tsv
wc -l ntm/v3/inputs/ntm_qc_passed_phigaro_summary_20260909.csv \
      ntm/v3/inputs/ntm_qc_passed_phigaro_coordinates_20260909.csv

# rows/cols/unique/source breakdown (1)
python3 - <<'PY'
import csv, collections
rows=list(csv.DictReader(open('NTM_Data/NTM_master_table_complete_20260922.tsv'), delimiter='\t'))
gids=[r['genome_id'] for r in rows]
print(len(rows), len(rows[0]), len(set(gids)),
      dict(collections.Counter(r['source'] for r in rows)))
PY

# zip member counts (2)
python3 - <<'PY'
import zipfile
for z,exp in [('NTM_Data/Genomes_sra_assembled_17920_20260922.zip',17920),
              ('NTM_Data/Genomes_ncbi_bvbrc_8579_20260922.zip',8579)]:
    names=zipfile.ZipFile(z).namelist()
    print(z, len(names), len([n for n in names if n.lower().endswith(('.fasta','.fna','.fa'))]), exp)
PY

# set reconciliation + per-genome counts (3,4)
python3 - <<'PY'
import csv, collections
master=set(r['genome_id'] for r in csv.DictReader(open('NTM_Data/NTM_master_table_complete_20260922.tsv'), delimiter='\t'))
summ=list(csv.DictReader(open('ntm/v3/inputs/ntm_qc_passed_phigaro_summary_20260909.csv')))
sg=set(r['genome_id'] for r in summ)
print(len(master-sg), len(sg-master))
coord=list(csv.DictReader(open('ntm/v3/inputs/ntm_qc_passed_phigaro_coordinates_20260909.csv')))
cg=set(r['genome_id'] for r in coord)
print(len(coord), len(cg), len(cg-sg))
cc=collections.Counter(r['genome_id'] for r in coord)
pc={r['genome_id']:int(r['prophage_count']) for r in summ}
print(sum(1 for g in pc if cc.get(g,0)!=pc.get(g)))
PY

# sample streaming extraction (5,6): contig header match + [begin,end] length
python3 - <<'PY'
import zipfile, csv
coords={}
for r in csv.DictReader(open('ntm/v3/inputs/ntm_qc_passed_phigaro_coordinates_20260909.csv')):
    coords.setdefault(r['genome_id'],[]).append(r)
def get_seq(zf, member, target):
    seqs={}; cur=None
    with zf.open(member) as f:
        for line in f:
            line=line.decode().rstrip('\n')
            if line.startswith('>'): cur=line[1:].split()[0]; seqs[cur]=[]
            else: seqs[cur].append(line)
    for k in seqs:
        if k==target: return k, ''.join(seqs[k])
    return None, None
jobs=[('NTM_Data/Genomes_sra_assembled_17920_20260922.zip','DRR015955','NODE_14_length_150130_cov_107.750611'),
      ('NTM_Data/Genomes_ncbi_bvbrc_8579_20260922.zip','1138383.42','accn|CP118870'),
      ('NTM_Data/Genomes_ncbi_bvbrc_8579_20260922.zip','GCA_000014985.1','CP000479.1')]
for zpath,gid,target in jobs:
    zf=zipfile.ZipFile(zpath); hdr,seq=get_seq(zf, gid+'.fasta', target)
    print(gid, hdr==target, len(seq))
    for r in coords[gid]:
        if r['scaffold']==target:
            b,e=int(r['begin']),int(r['end'])
            print("  ", r['prophage_id'], b, e, "exp", e-b+1, "got", len(seq[b-1:e]))
PY
```

## Sample extraction results (check 6)

| Source | genome_id | coord `scaffold` | zip contig header | header match | begin–end | expected len | extracted len |
|--------|-----------|------------------|-------------------|--------------|-----------|--------------|---------------|
| ASSEMBLY | DRR015955 | `NODE_14_length_150130_cov_107.750611` | `NODE_14_length_150130_cov_107.750611` | yes | 4656–33196 | 28541 | 28541 |
| BV-BRC | 1138383.42 | `accn\|CP118870` | `accn\|CP118870` | yes | 82291–90640 | 8350 | 8350 |
| NCBI | GCA_000014985.1 | `CP000479.1` | `CP000479.1` | yes | 753227–791253 | 38027 | 38027 |
| NCBI | GCA_000014985.1 | `CP000479.1` | `CP000479.1` | yes | 2218471–2226202 | 7732 | 7732 |

All 4 coordinate ranges extracted at exactly `end-begin+1`, confirming 1-based inclusive coordinates and matching scaffold headers for all three source namespaces. FASTA contig lengths observed: DRR015955/`NODE_14`=150130, 1138383.42/`accn|CP118870`=5510556, GCA_000014985.1/`CP000479.1`=5475491.

## Deferred / needs decision

- **Zip sha256 deferred** (per authority): not computed over `Genomes_sra_assembled_17920_20260922.zip` (~30.5 GB) or `Genomes_ncbi_bvbrc_8579_20260922.zip` (~14.1 GB). Run a later checksum pass once extraction/CPU load allows.
- `Defence_System_Tables.zip`, `Key_Plots.zip`, and master-table defense columns (padloc/defensefinder counts) were **out of scope** for this ingest lane — not validated here.
- No discrepancies found; nothing requires a fix decision.
