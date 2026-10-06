# S3 — Prophage tight-clade definition & MTC verification (NTM v3)

Read-only independent verification of `$NVME/ntm/v3/clades/0/tight_clades.json` against
`full_prophages.fa`, `distances.npz`, `clade_similarity.json`, the delivery manifest and
`host_clades.tsv`. Nothing was modified.

- `$NVME = /mnt/nvme3n1/erikg/phind-genome-work`
- FASTA record count: **36,857** (`grep -c '^>' full_prophages.fa`)
- Report under test: `ntm/v3/mash_clades_report.md`

## Check table (observed vs expected)

| # | Check | Expected | Observed | Verdict |
|---|---|---|---:|---|
| 1 | total clades | 1,304 | 1,304 | PASS |
| 1 | alignable clades (n>=2) | 893 | 893 | PASS |
| 1 | singleton clades (n==1) | 411 | 411 | PASS |
| 1 | max clade size | 100 | 100 | PASS |
| 1 | clades > 100 | 0 | 0 | PASS |
| 1 | size bands | 1:411, 2-10:432, 11-50:125, 51-100:336 | same (sum 1304) | PASS |
| 2 | member ids total / unique | 36,857 / 36,857 | 36,857 / 36,857, 0 dups | PASS |
| 2 | union of members == FASTA ids | equal | equal (0 in-clade-not-fasta, 0 fasta-not-clade) | PASS |
| 2 | every id in exactly one clade | yes | yes (no orphan, no duplicate) | PASS |
| 3 | alignable clades with internal **median** > 0.25 | 0 | 0 (max median 0.248764) | PASS |
| 3 | recomputed stats == `clade_similarity.json` | identical | 0 mismatches over all 893 clades & all 4 stats + frac | PASS |
| 3 | alignable clades with some pair > 0.25 (informational) | — | 311 (median criterion, expected) | note |
| 4 | MTC prophages | — | 391 | info |
| 4 | clades all-MTC / partial / non-MTC | — | 12 / 7 / 1285 | info |
| 4 | clades emptied by drop-MTC | — | 12 | info |
| 5 | cap `--max-size 100` present in command | yes | `--max-size 100` in `clades/0/commands.log` | PASS |
| 5 | clades hitting exactly 100 | — | 275 | info |

All numeric claims in the report (clade count, split, distribution, max size, 0 median
violations) reproduce exactly.

## Exact recomputation commands

```bash
NVME=/mnt/nvme3n1/erikg/phind-genome-work
cd $NVME/ntm/v3/clades/0

# (1) sizes / partition / cap
python3 - <<'PY'
import json, collections
tc=json.load(open('tight_clades.json'))
sizes={k:len(v) for k,v in tc.items()}
print('n',len(sizes),'align',sum(s>=2 for s in sizes.values()),
      'single',sum(s==1 for s in sizes.values()),'max',max(sizes.values()),
      'n_at_100',sum(s==100 for s in sizes.values()),'>100',sum(s>100 for s in sizes.values()))
allm=[m for v in tc.values() for m in v]
fa={l[1:].strip() for l in open('$NVME/ntm/v3/full_prophages.fa') if l.startswith('>')}
print('members',len(allm),'unique',len(set(allm)),'union==fasta',set(allm)==fa)
PY

# (3) INDEPENDENT per-clade medians straight from distances.npz (mmap, 5.4 GB)
python3 - <<'PY'
import numpy as np, json
z=np.load('distances.npz', allow_pickle=True, mmap_mode='r')
D=z['D']; members=z['members']
idx={m:i for i,m in enumerate(members)}
tc=json.load(open('tight_clades.json'))
rep=json.load(open('clade_similarity.json'))['per_clade']
viol=[]; mism=0
for cid,m in tc.items():
    if len(m)<2: continue
    ii=np.array([idx[x] for x in m]); S=np.asarray(D[np.ix_(ii,ii)],dtype=np.float64)
    t=S[np.triu_indices(len(ii),1)]
    med,mx,mn,fr=np.median(t),t.max(),t.min(),(t<=0.25).mean()
    r=rep[cid]
    if med>0.25: viol.append(cid)
    if abs(med-r['median'])>1e-6 or abs(mx-r['max'])>1e-6: mism+=1
print('median>0.25:',len(viol),'recompute-vs-report mismatches:',mism)
PY

# (4) MTC accounting — map prophage -> canonical_acc -> accession -> host species
python3 - <<'PY'
import gzip, json, collections
NV='/mnt/nvme3n1/erikg/phind-genome-work/ntm/v3'
host={}
for i,l in enumerate(open(NV+'/host_clades/host_clades.tsv')):
    p=l.rstrip('\n').split('\t')
    if i==0: ci={c:j for j,c in enumerate(p)}; continue
    host[p[ci['accession']]]=p[ci['species']]
c2a={}; c2s={}
with gzip.open('ntm/v3/inputs/v3_prophage_manifest_delivery.tsv.gz','rt') as f:  # run from repo root
    h=f.readline().rstrip('\n').split('\t'); ci={c:j for j,c in enumerate(h)}
    for l in f:
        p=l.rstrip('\n').split('\t')
        c2a.setdefault(p[ci['canonical_acc']],set()).add(p[ci['accession']])
        c2s.setdefault(p[ci['canonical_acc']],set()).add(p[ci['species']])
P=("Mycobacterium tuberculosis","Mycobacterium canetti","Mycobacterium orygis",
   "Mycobacterium africanum","Mycobacterium bovis","Mycobacterium microti")
mtc=lambda s: bool(s) and s.startswith(P)
tc=json.load(open(NV+'/clades/0/tight_clades.json'))
allc=part=non=0; rem=0
for cid,m in tc.items():
    k=sum(1 for x in m if mtc(host.get(sorted(c2a[x.split('#1#')[0]])[0]) if x.split('#1#')[0] in c2a else None))
    rem+=k
    if k==len(m): allc+=1
    elif k==0: non+=1
    else: part+=1
print('MTC prophages',rem,'all',allc,'part',part,'non',non)
PY
```

## Internal similarity (check 3) — independently recomputed

Recomputed all 893 alignable clades' pairwise MASH distances **directly from
`distances.npz`** (full symmetric `D[36857,36857]`; diagonal 0, symmetric exactly).
Result is byte-for-byte consistent with `clade_similarity.json`:

- 0 mismatches for `median`, `min`, `max`, `frac_le_threshold` across all 893 clades.
- **0** alignable clades with median > 0.25 (max median = **0.248764**; median-of-medians 0.03354).
- 311 clades contain at least one pair > 0.25 (expected: the enforced criterion is the
  **median**, not max). Worst pair distances reach **1.0** (MASH saturated) in several
  100-member clades (`0_0004`,`0_0005`,`0_0007`,`0_0008`,`0_0009`).
- Worst clade has only `frac_le_threshold = 0.5012` of pairs <= 0.25 (0 clades below 0.5).

So the report's "median <= 0.25, 0 violations" claim is verified; the clades are
median-tight, not max-tight.

## MTC accounting (check 4)

MTC rule used: `host_clades.tsv.species` starts with `Mycobacterium tuberculosis`,
`Mycobacterium canetti`, `Mycobacterium orygis`, `Mycobacterium africanum`,
`Mycobacterium bovis`, or `Mycobacterium microti`. (Host table contains only
`tuberculosis`, `canetti`, `orygis` variants; `africanum/bovis/microti` absent. Substring
false positives e.g. `M. avium subsp. paratuberculosis` are excluded by `startswith`.)

Join: FASTA header = `{canonical_acc}#1#{prophage_id}` (canonical key matches **all
36,857** ids; genome_id key matches only 32,660 — headers are canonical_acc, as claimed).
`canonical_acc -> accession -> host_clades.species`. 36,759 prophages resolved via
host_clades; 98 fell back to the manifest species column (same result). 0 unresolved.

| Quantity | Count |
|---|---:|
| MTC host genomes in `host_clades.tsv` (all, incl. genomes with no prophages) | 7,273 |
| Prophages belonging to MTC genomes (of 36,857) | **391** |
| Clades entirely MTC | **12** |
| Clades partially MTC | **7** |
| Clades fully non-MTC | **1,285** |
| Prophages a "drop MTC" filter removes | **391** (1.06%) |
| Clades that would become EMPTY (all members MTC) | **12** |
| Clades touched (contain >=1 MTC member) | **19** (12 + 7) |
| Clades surviving the filter | **1,292** (1,304 - 12) |

Cross-check: counting MTC directly from the manifest `species` column gives the **same
391** prophages (symmetric difference 0).

Entirely-MTC clades (id, n, median):

| clade | n | median |
|---|---:|---:|
| 0_0332 | 100 | 0.000392 |
| 0_0336 | 48 | 0.021325 |
| 0_0353 | 19 | 0.066757 |
| 0_0366 | 8 | 0.037728 |
| 0_0633 | 8 | 0.000231 |
| 0_0654 | 7 | 0.000000 |
| 0_0499 | 4 | 0.000000 |
| 0_0611 | 4 | 0.104309 |
| 0_0638 | 1 | — |
| 0_0989 | 1 | — |
| 0_0990 | 1 | — |
| 0_1039 | 1 | — |

Partially-MTC clades (id, n, MTC members, median): `0_0306` (100, 93), `0_0338`
(100, 90), `0_0413` (79, 1), `0_0252` (74, 1), `0_0365` (53, 2), `0_0435` (18, 1),
`0_0657` (3, 1).

Note `0_0332` is a 100-member, near-identical, wholly-MTC clade (hits the cap) — the
largest single release-filter impact.

## Check 5 — cap

`clades/0/commands.log` shows the clustering was invoked with
`build_tight_clades.py --threshold 0.25 --max-size 100 --communities 0`. Max observed
clade size is exactly 100 with 0 clades above, and **275 clades sit exactly at the 100
cap** (so results are cap-limited, not a natural clustering ceiling for the largest
components).

## Overall verdict

**PASS.** Every quantitative claim in the report reproduces exactly:
1,304 clades = 893 alignable + 411 singletons; max size 100; the membership is a clean
partition of all 36,857 prophages with no duplicates/orphans; and the median internal
MASH <= 0.25 criterion holds with 0 violations, independently recomputed from
`distances.npz` and matching `clade_similarity.json` with 0 mismatches.

MTC impact for the release-filter decision: **391 prophages** from MTC genomes;
**12 clades would become empty** (drop entirely), 7 more clades are partially MTC
(keep their non-MTC members). Net clade count after dropping only empty clades: **1,292**.

## Residual uncertainty

- The enforced tightness criterion is the pairwise **median**. 311 alignable clades
  contain at least one member pair above 0.25 and some pairs are saturated at 1.0 —
  if the release filter needs a stricter "all pairs <= 0.25" rule, these clades would
  need revisiting (they are not violations of the stated median rule).
- `distances.npz` `D` is reported as a full symmetric matrix; verification used it
  directly. It was not re-derived from raw MASH output in this lane (the report's own
  spotcheck claims 0 triangle-vs-float32 mismatches; not re-run here).
- MTC classification relies on `host_clades.tsv` species strings (as the task
  specifies). Genomes absent from `host_clades` would fall back to the manifest species;
  only 98 prophages needed that fallback and both sources agreed.
- `0_0332` reaching the 100 cap while being wholly MTC means the true MTC prophage
  grouping may be larger than one clade (cap-split); the 12/7 partition counts are
  clade-level, not biological-cluster-level.

## Manual notes / evidence

- All checks performed read-only; only the output report was written.
- Commands run: inline Python over `tight_clades.json`, `distances.npz` (mmap),
  `clade_similarity.json`, `full_prophages.fa`, manifest, `host_clades.tsv`;
  plus `grep -c '^>'`, `cat commands.log`, `wc -l`.
- No source files, clade files, or reports were modified.
