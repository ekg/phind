# NTM v3 release - scoping statement

Generated: 2026-10-06T19:23:08Z (unified Mycobacteriaceae catalog)

`host_scope` is a **label**, not a filter: no population is dropped. Per the binding scope decision (`ntm/v3/PLAYBOOK.md`), the release is one catalog spanning the whole family.

## Label namespaces

- per genome (host) `host_scope` in {NTM, MTC} - binary only.
- per prophage clade label in {NTM, MTC, MIXED} (MIXED = cross-boundary; see `cross_boundary_report.tsv`).

## Clade-level scope derivation (binding)

A MASH-defined host clade is MTC when the majority of its genomes carry an MTC species prefix (tuberculosis / canetti / orygis / africanum / bovis / microti, via `startswith`), else NTM. Every prophage member inherits its host clade's scope; the per-prophage-clade label is MTC/NTM when all members inherit that scope, and MIXED when both occur.

The per-genome MTC species-prefix classifier is kept and used only to decide each host clade's majority. Rationale: a host clade is an ANI-defined lineage, so scope is a property of the lineage rather than of an individual unresolved species label. Members labeled `Mycobacterium sp.` that sit in host_clade_0002 (7273/7467 = 97.4% MTC) are unresolved-ANI MTC genomes, so the species-prefix rule under-calls them; the clade-level rule assigns them MTC.

Per-genome (host) scope is auditable in `host_clade_scope.tsv`.

MTC host clades (majority MTC): 1 (host_clade_0002).

### Clade-derived MTC assignments (unresolved species labels)

9 prophage members carry a non-MTC species label but sit in an MTC host clade, so they are MTC under the clade-level rule:

| prophage clade | clade-derived MTC members |
|---|---:|
| 0_0306 | 1 |
| 0_0332 | 5 |
| 0_0336 | 2 |
| 0_0338 | 1 |

### Clade-level scope vs species-prefix accounting

| rule | MTC members | wholly-MTC clades | MIXED clades | fully-NTM clades |
|---|---:|---:|---:|---:|
| clade-level scope (released) | 391 | 12 | 7 | 1285 |
| species-prefix only | 382 | 10 | 9 | 1285 |

Members by clade-derived member-scope: NTM=36466, MTC=391 (of 36857).

| host_scope | clades | ML genomes | ancestral genomes | members | host clades |
|---|---:|---:|---:|---:|---:|
| NTM | 1285 | 1285 | 878 | 36228 | 289 |
| MTC | 12 | 12 | 8 | 202 | 1 |
| MIXED | 7 | 7 | 7 | 427 | 26 |

_Subset run: 1304 clades (full set)._
