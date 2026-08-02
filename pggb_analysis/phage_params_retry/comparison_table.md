# cluster_6 pggb parameter comparison (old vs corrected phage params)

Corrected params (this run): `p=75 s=250 l=500 k=11 ani-diff=80` via pggb 0.6.0 (`-p 75 -s 250 -l 500 -k 11 -g 80 -c 5 -n 53`).

Historical params: `wfmash -p 85 -n 5 -k 19 -l 2000 -X`, `seqwish -k 19`, `smoothxg -X 100 -r 53`, `gfaffix`.

| metric | new (p=75,s=250,l=500,k=11,ani80) | old (p=85,l=2000,k=19) | notes |
|---|---|---|---|
| Graph segments (final) | 17,313 | 17,663 | lower = less fragmentation |
| Graph edges (final) | 23,869 | 24,456 |  |
| Graph total bp (final) | 75,459 | 77,304 |  |
| Segment N50 (final, bp) | 4 | 5 | higher = longer clean nodes |
| Seqwish segments (pre-smooth) | 10,422 | 7,689 |  |
| Seqwish paths | 53 | 53 |  |
| Final graph paths (P/W lines) | 53 | 130 | old had chimeric 69kb traversal path |
| Prophage paths | 53 | 53 | should = 53 |
| Prophage path total bp | 817,520 | 1,014,303 |  |
| Prophage path N50 (bp) | 14,752 | 19,512 |  |
| Prophage path mean len (bp) | 15424.9 | 19137.8 |  |
| Consensus paths in graph | 0 | 0 | smoothxg POA consensus paths |
| Mean query coverage (PAF) | 0.9234 | 0.9138 | higher = more prophage captured |
| PAF records | 1,347 | 2,128 |  |
| Alignment columns | 93,741 | 157,426 |  |
| Ancestral consensus len (raw paths, bp) | 51,439 | 69,153 | like-for-like with old pipeline ancestral |
| Ancestral confidence (raw paths) | 0.6041 | None |  |
| Consensus len (MAFFT-majority, bp) | 93,741 | 69,153 | new-only (alignment-based) |
| Consensus method | majority-rule (MAFFT-aligned paths) | pipeline majority-rule (raw paths) |  |
