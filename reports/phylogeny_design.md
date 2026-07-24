# Scalable host-genome distance/tree design for 24,000–26,000 *E. coli* assemblies

**Design date / source access date:** 2026-07-24
**Scope:** host-genome QC, distances, tree-like summaries, and the later mapping of prophage traits. This report does **not** acquire genomes, define storage policy, audit current prophage tables, test IMPG, build a pangenome, or run production-scale computation.

## Decision in one paragraph

Use a **two-resolution design**. First, after assembly/species QC, make one Mash sketch per host assembly with Mash v2.3 candidate parameters `k=21, s=10,000`, compare sensitivity to `s=1,000/50,000` and a second `k`, and—only if measured gates pass—construct an all-host **Mash-distance neighbour-joining overview** with `mash triangle` plus a memory-bounded RapidNJ implementation. Call this an **unrooted host genomic-similarity dendrogram**, not a resolved organismal phylogeny. Do not run MashTree directly at 24k–26k: v1.4.6 is useful for a bounded representative pilot and sketch-resampling support, but internally materializes a PHYLIP matrix and its own documentation explicitly says it does not infer phylogeny. Second, choose host-genetic representatives and boundary cases without looking at prophages, infer within-lineage core-genome SNP trees, mask recombination with Gubbins, and obtain model-based ML support with IQ-TREE. Use those higher-fidelity trees to validate or revise host-only clades. Retain every host in the analysis manifest; exact/near duplicates and non-representatives keep explicit representative/placement links and remain in prevalence denominators.

## 1. What “mesh distance” and “mesh triangle” mean here

### Installed evidence

`artifacts/phylogeny_probe/tool_versions.txt` records the 2026-07-24 search of `PATH`, Debian packages, and all enumerated micromamba environments. None of `mesh`, likely “mesh” variants, `mash`, `mashtree`, `rapidnj`, or the higher-fidelity tools was installed or discoverable. Therefore this report makes **no claim that any executable or version is currently available** and did not fabricate local help output. `artifacts/phylogeny_probe/commands.sh` reproduces that bounded discovery; `environment.txt` records the host. No sequence probe was run (`input_count=0`, `input_bytes=0`, `output_bytes=0`).

The proposed correction is documentation-based: the versioned Mash v2.3 source defines `mash dist` as pairwise reference/query distance and `mash triangle` as “Estimate a lower-triangular distance matrix”; it contains no `mesh` command [1–3]. Thus the user's phrases are interpreted as **Mash distance** and **Mash triangle**, contingent on installing and pinning Mash v2.3 and capturing its real `--version`/`--help` before any pilot. A literal executable called `mesh` was not found.

### Command semantics that matter

* `mash sketch` creates reduced MinHash representations. Sketching first is preferable to repeatedly sketching raw inputs [2].
* `mash dist REF QUERY` emits reference ID, query ID, Mash distance, p-value, and shared-hash count. It is rectangular (`q × r`) when sketches contain multiple entries.
* `mash triangle` computes every unordered pair and emits relaxed lower-triangular PHYLIP; `-E` emits an edge list, and `-d`/`-v` filter the **reported** edges [1]. Filtering output does not avoid computing the pairs in the v2.3 loop [1].
* Whole files are sketches by default. `-i` switches to individual FASTA records [1,3]. For a multi-contig assembly, **do not use `-i`**: one assembly file must be one host observation.

## 2. What each route can and cannot claim

| Route | Output and scaling | Appropriate claim | It cannot establish |
|---|---|---|---|
| Mash `sketch` + `dist` | Rectangular distances; useful for reference screens, blockwise comparisons, nearest representatives | Approximate whole-assembly k-mer dissimilarity; fast screening | A tree, vertical descent, SNP counts, ANI equality, direction/time, or branch support |
| Mash `sketch` + `triangle` + RapidNJ | All-pair lower matrix then NJ; distance stage and input are O(n²); external-memory RapidNJ can bound tree RAM [1,7] | Unrooted host genomic-similarity dendrogram and coarse neighbourhoods | A substitution-model phylogeny; recombination-free clonal history; trustworthy resolution of near-identical isolates; support from one matrix |
| MashTree v1.4.6 | Automates Mash `k=21,s=10,000`, all-pair distances, PHYLIP, QuickTree/BioPerl NJ; optional hash bootstrap/jackknife [4,5] | Convenient bounded representative tree and stability to MinHash resampling | Its own algorithm document says “Mashtree does not infer phylogeny” [5]; sketch support is not site bootstrap and does not cure biological/model bias |
| ANI screen (e.g. skani 0.3.2 documentation) | Query-to-reference DB can be near O(n) in the number of hosts; sampled pair validation | Species/boundary QC using ANI **and aligned fraction**; validation of Mash distance ranges | Topology or ancestry; ANI alone is not phylogeny and 95% is not a clade definition [8,9] |
| Core-genome SNP alignment → Gubbins → IQ-TREE | Costly; perform within host lineages on a preregistered representative/boundary subset | Model-based, recombination-filtered estimate of the sampled clonal frame with branch support and explicit missingness | The history of accessory loci, plasmids, horizontally transferred prophages, unsampled hosts, or a single species-wide genealogy free of uncertainty |
| Single-copy core-gene concatenation/consensus | Alternative where reference mapping is poor; gene-tree discordance can be examined | Species/lineage phylogeny conditional on orthology/alignment/model | A prophage tree or proof that every gene shares one history; no pangenome is built under this task |

Mash distance is a transform of estimated k-mer-set Jaccard,

\[
D(k,j)=-\frac{1}{k}\ln\left(\frac{2j}{1+j}\right),\quad 0<j\le1,
\]

and is set to 1 when no hashes are shared [2,10]. Mash reports a random-match p-value, not a confidence interval on topology. The Mash paper found a strong empirical relation to `1 − ANI` in a tested 90–100% ANI range, but this does not make Mash distance ANI or a nucleotide substitution model [10]. Fragmentation, missing sequence, contamination, repeats, plasmids, horizontally acquired islands, and assembly errors change the k-mer sets.

## 3. Exact scale and resource arithmetic

For `n` hosts, the number of unordered, off-diagonal pairs is

\[
P(n)={n\choose2}=\frac{n(n-1)}2.
\]

| n | Exact unordered pairs P(n) | Dense n×n float32 | Dense n×n float64 | Packed off-diagonal triangle float32 | Packed off-diagonal triangle float64 |
|---:|---:|---:|---:|---:|---:|
| 24,000 | **287,988,000** | 2,304,000,000 B = 2.304 GB = **2.145767 GiB** | 4,608,000,000 B = 4.608 GB = **4.291534 GiB** | 1,151,952,000 B = **1.072839 GiB** | 2,303,904,000 B = **2.145678 GiB** |
| 26,000 | **337,987,000** | 2,704,000,000 B = 2.704 GB = **2.518296 GiB** | 5,408,000,000 B = 5.408 GB = **5.036592 GiB** | 1,351,948,000 B = **1.259100 GiB** | 2,703,896,000 B = **2.518199 GiB** |

GB is decimal (`10^9`); GiB is binary (`2^30`). These are payload-only lower bounds. Arrays of pointers/objects, row indices, names, parsers, NJ search structures, temporary copies, and allocator overhead can dominate them.

### CPU and I/O implications

* Mash v2.3 compares sorted hashes until at most `s` hash positions have been considered [1]. `P×s` is a useful upper-scale indicator, not a runtime prediction: at `s=10,000`, it is 2.87988×10¹² and 3.37987×10¹² examined-position opportunities. Similar pairs, implementation details, threads, NUMA, and I/O change actual cost.
* If a measured pilot sustains `R` completed pairs/s, distance-only wall time is at least `P/R`. Illustrative planning values—not measurements—are: at 1,000 pairs/s, 80.00 h and 93.89 h; at 10,000 pairs/s, 8.00 h and 9.39 h; at 100,000 pairs/s, 0.800 h and 0.939 h. Production approval must substitute the lower confidence bound of pilot throughput.
* At `k=21`, Mash uses 64-bit hashes under its documented hash-width rule [2]. Full `s=10,000` sketches therefore contain at least `n×s×8` hash bytes: 1.92 GB (24k) or 2.08 GB (26k), before metadata/index overhead. Increasing `s` approximately increases sketch payload and comparison work linearly.
* A PHYLIP text token stream is larger than packed binary. If the mean serialized distance plus delimiter is `b` bytes, the distance text alone is `P×b`. For an explicitly illustrative `b=11`, that is 3.168 GB (24k) or 3.718 GB (26k), excluding labels/newlines. A five-field edge list can be much larger; at an illustrative 32 B/record it is 9.216 GB or 10.816 GB. Measure actual `wc -c`, do not reserve from these examples alone.
* Naive NJ is O(n³). RapidNJ accelerates pair selection and provides memory-efficient/external-disk modes, but still consumes a complete distance matrix and can become random-I/O bound [7]. External memory trades RAM for scratch capacity and IOPS; it does not change the evidentiary limits of Mash distances.
* MashTree v1.4.6 is especially risky at this size: it issues `mash dist -t` for each sketch, stores distances in Perl hashes/a database, serializes a full PHYLIP string, and invokes QuickTree or BioPerl [4,6]. Its convenient wrapper is not the production default.

The captured host happens to report 1.0 TiB RAM and 2.5 TiB free on its current filesystem, but that is **not** an allocation or a performance result. Scheduler limits, local scratch, concurrent users, and I/O policy must be measured.

## 4. Primary host-only workflow

### 4.1 Immutable inputs and QC (consume, do not acquire)

Consume the acquisition/storage task's manifest with one row per assembly and at least:

`assembly_id`, versioned accession, sequence path, SHA-256, retrieval/source date, assembly level, release version, length, contig count, N50, N fraction, completeness, contamination, taxonomy result, and exclusion reason.

Required gates before sketching:

1. **Identity/version:** stable `assembly_id` unique; versioned accessions only; sequence checksum agrees with manifest; annotation/prophage results later must point to the identical sequence checksum. Never silently mix `.1` and `.2` releases.
2. **Assembly QC:** record length, contigs, N50/L50, ambiguous bases, and reference-based misassembly metrics where an appropriate reference exists (QUAST documents these metrics [11]). Flag extreme length/fragmentation using robust cohort distributions rather than N50 alone.
3. **Completeness/contamination:** candidate inclusion gate ≥95% completeness and ≤5% contamination, with the estimator/version/database recorded. Those project thresholds are deliberately stricter than “medium quality”; they are not a universal definition. Investigate or quarantine failures. Contamination may create false shared k-mers; incompleteness may inflate distances. CheckM2 is MAG-oriented, so for isolate assemblies corroborate surprising values with lineage markers/read/assembly evidence rather than treating an ML score as truth [12].
4. **Species:** search every host against a curated, versioned *E. coli/Shigella/Escherichia* reference panel using ANI plus both aligned fractions. Candidate gate: best appropriate reference ANI ≥95% and aligned fraction ≥65% in each reported direction; manually adjudicate borderlines and known taxonomic ambiguity. This is a project QC gate, not a clade definition. Outgroups are exempt and marked `role=outgroup`.
5. **Sample independence:** resolve multiple assemblies for one biological isolate. Keep the best assembly as primary; preserve aliases. Technical duplicates must not masquerade as biological replication.

### 4.2 Tip IDs and the BGZF/PanSN hand-off

The tree tip is a **host assembly/sample**, never a contig. Use filesystem-safe assembly IDs (for example a versioned assembly accession) as staged basenames and maintain `tip_id ↔ assembly_id ↔ sequence checksum ↔ original path`. With multi-contig assemblies, sketch the whole file and never pass Mash `-i`.

Canonical BGZF input and literal PanSN `sample#haplotype#contig` headers are owned by `pansn-bgzip-genome-layout`. This report did not test or duplicate that compatibility matrix. **NO-GO** for production until that task confirms the pinned Mash/MashTree/other tools can read the chosen BGZF files and preserve literal `#` safely. If not, create a derived ordinary-FASTA analysis view under that task's policy while retaining the checksum/crosswalk. Shell paths must always be quoted. PanSN contig headers must not become thousands of tree tips.

### 4.3 Candidate baseline commands (templates; not executed)

Pin versions/containers and replace variables only after the ≤20-input smoke test and help capture. The list must contain one safely named assembly file per line.

```bash
# Evidence first
mash --version > versions.txt
mash sketch --help >> versions.txt 2>&1
mash dist --help >> versions.txt 2>&1
mash triangle --help >> versions.txt 2>&1
rapidnj --help >> versions.txt 2>&1

# One compound archive containing one whole-file sketch per host; no -i.
mash sketch -p "$CPUS" -k 21 -s 10000 \
  -l host_assemblies.list -o host_k21_s10000
mash info host_k21_s10000.msh > host_k21_s10000.info.txt

# Full route only after gates. Capture /usr/bin/time -v, stderr, checksums, wc -c.
/usr/bin/time -v mash triangle -p "$CPUS" \
  host_k21_s10000.msh > host_k21_s10000.lt.phylip \
  2> host_k21_s10000.triangle.time.txt

# Confirm exact PHYLIP compatibility in the tiny/pilot stages first.
/usr/bin/time -v rapidnj host_k21_s10000.lt.phylip \
  -i pd -o t -m "$RAPIDNJ_MB" -x host_mash_nj.unrooted.nwk \
  2> host_mash_nj.rapidnj.time.txt
```

`mash dist` is preferable for a representative fallback or blocks:

```bash
# Each query block and the representative archive are pre-sketched identically.
mash dist -p "$CPUS" representatives.msh query_block.msh \
  > block_to_representatives.tsv
```

Validate symmetry, diagonal zero, finite range `[0,1]`, unique IDs, exactly `P(n)` triangle tokens, and one Newick tip per accepted host. Preserve the unrooted tree and matrix checksums.

### 4.4 Sketch sensitivity

Mash v2.3 defaults are `k=21,s=1,000`; MashTree v1.4.6 defaults to `k=21,s=10,000` [3,4]. Candidate production baseline is `k=21,s=10,000`, not an unquestioned truth.

* Mash documentation gives the generic sketch error scale `sqrt(1/s)` [2]: 0.03162 at 1k, 0.01000 at 10k, and 0.004472 at 50k. This is not a topology support value or a Mash-distance confidence interval.
* Larger `k` raises specificity but lowers sensitivity; smaller `k` does the reverse [2]. Test `k=21` against `k=31` on the same pilot; heed Mash's random-match warning and p-values. Do not combine sketches with mismatched `k`; mixed `s` comparisons use the smaller sketch [2].
* Compare `s=1k,10k,50k` on preregistered pilot pairs spanning near duplicates, within-lineage, between-lineage, poor-but-accepted assemblies, and outgroups. Record rank correlation, nearest-neighbour agreement, split agreement, zero/shared-hash counts, and distance/p-value distributions.
* Close outbreak-scale isolates may collapse to identical/unstable sketches even at large `s`. That is a trigger for core-SNP analysis, not a reason to report arbitrary Mash branch order.

### 4.5 Near duplicates and keeping all hosts

1. Collapse exact sequence SHA-256 duplicates immediately into an equivalence class while preserving every sample/accession row.
2. Define a near-duplicate threshold only after the sensitivity pilot and core-SNP calibration; do not hard-code `Mash D=0` as biological identity. Choose the highest-QC medoid per class using host/assembly data only.
3. Build the high-resolution tree on medoids plus diverse and boundary representatives. Assign non-representatives to candidate clades by distances to multiple representatives, not one centroid. Mark ties/discordance `ambiguous`.
4. Retain all hosts in a membership table. For visualization, exact duplicates may be expanded as zero-length polytomies and near duplicates shown as collapsed fans; these are display/placement conventions, **not** an NJ tree recomputed with all tips. Within clusters where transmission-level order matters, run a dedicated core-SNP tree.
5. All later prophage prevalence denominators use biological hosts, with an explicit rule for repeated isolates; no host disappears merely because it was not a high-resolution representative.

## 5. Avoiding or staging explicit O(n²)

Options are ordered by fidelity to the all-host NJ target:

1. **Full triangle, staged:** simplest and exact for the chosen Mash sketches. It materializes O(n²), so use only if measured full-cost and tree gates pass.
2. **Blockwise stream/filter:** use query/reference compound sketches and immediately retain distances below a predeclared host-genetic threshold. This bounds peak RAM and permits restart/checksum by block. It still computes O(n²) unless candidate generation reduces comparisons, and a sparse threshold graph is a network/cluster result—not a complete NJ tree.
3. **Representative tree plus all-host assignment (recommended fallback):** choose `r` representatives using only host QC/genetics; compute `r(r−1)/2` for the overview tree and up to `n×r` for explicit assignment. Example only: `r=2,000` gives 49,999,000 comparisons for 24k and 53,999,000 for 26k, versus 287,988,000/337,987,000. This is not the same topology as all-host NJ, but it keeps all hosts mapped and reserves expensive inference for informative diversity.
4. **Host-genetic candidate indexing/bucketing:** a documented sketch database/LSH/cgMLST/ANI search can propose neighbours; then compute exact Mash/skani values on candidates and bridge representatives. Approximate retrieval must be recall-tested against an exact pilot. Connected components can chain through weak edges, so check bridges and do not call graph clusters a phylogeny.
5. **Hierarchical divide-and-validate:** coarse host-only groups → within-group representative/core-SNP trees → a backbone of group representatives. Report the backbone, local trees, and membership separately; do not splice branch lengths as if one homogeneous model produced them.

`mash triangle -E -d X` reduces **output** when few edges pass, not the v2.3 comparison loop [1]. MashTree is not an O(n²)-avoiding alternative.

## 6. Higher-fidelity validation/reference route

### Sampling without circularity

Freeze a representative plan before opening prophage traits: medoids, QC extremes, Mash-distance extremes, suspected boundaries, every host phylogroup/lineage represented, and verified outgroups. Stratify within provisional Mash groups, but use no prophage count, prophage gene, IMPG cluster, or phage similarity. Gubbins v3.4.3 documentation warns it is for limited-diversity samples sharing a recent ancestor, scales approximately quadratically, and should not be run across species-wide diversity [13]. Therefore run it per lineage, generally hundreds rather than all 26k, and connect interpretation through a separately validated backbone.

### Assembly core-SNP template (not executed)

Snippy v4.6.0 documents `--ctgs` and `snippy-core` for assemblies [14]. A reference-mapping route is:

```bash
# Same lineage-appropriate, versioned reference for every member of a run.
snippy --cpus "$CPUS_PER_SAMPLE" --outdir "snippy/$TIP" \
  --ref lineage_reference.gbk --ctgs "$ASSEMBLY"
snippy-core --prefix lineage_core snippy/*

# Gubbins needs a whole-genome alignment with spatial context, not a
# concatenation of core genes and not just a list of variable sites.
run_gubbins.py --prefix lineage_gubbins --threads "$CPUS" \
  --first-tree-builder rapidnj --first-model JC \
  --tree-builder iqtree --model GTR lineage_core.full.aln

mask_gubbins_aln.py --aln lineage_core.full.aln \
  --gff lineage_gubbins.recombination_predictions.gff \
  --out lineage_core.recombination_masked.full.aln

# Full alignment contains invariant sites: do not add +ASC.
iqtree3 -s lineage_core.recombination_masked.full.aln \
  -m MFP -B 1000 --alrt 1000 -T AUTO --prefix lineage_ml
```

If the final input contains **only variable SNP sites**, use a tested ascertainment-bias model such as `+ASC`; do not apply `+ASC` to a full alignment with constant sites. IQ-TREE 3.1.2 documents UFBoot and ascertainment correction [15]. Inspect composition/model warnings; bootstrap does not rescue a bad alignment or violated assumptions.

A core-gene alternative uses rigorously single-copy orthologues, per-gene codon-aware alignment/masking, partitioned ML, concordance factors or gene-tree discordance. It is a validation route, not authorized pangenome computation here.

### Recombination and reference bias

* Whole-genome Mash mixes clonal mutations, homologous recombination, mobile/accessory DNA, plasmids, and loss/gain. NJ cannot separate them.
* Gubbins identifies clustered substitutions consistent with recombination, but can confuse assembly/alignment errors, mutational hotspots, or long-branch effects; its documentation recommends subdividing diverse populations [13]. Report masked fraction, `r/m`, callable clonal-frame length, and sensitivity to lineage/reference choice.
* Mapping all assemblies to one distant reference loses lineage-specific sequence and can distort missingness. Use lineage-appropriate references and compare a core-gene/backbone result. Reject runs where reference choice moves supported boundaries.

### Rooting and outgroups

Keep primary Mash and core trees unrooted. Add one or preferably several independently verified *Escherichia* outgroups (candidate examples require taxonomic verification, such as *E. fergusonii*) selected before phage analysis. Infer the ingroup tree both with and without outgroups; root on a supported outgroup branch only if outgroups are monophyletic and do not destabilize ingroup splits. A distant outgroup may reduce core alignment and cause long-branch attraction. Midpoint rooting is an exploratory display assumption, not evidence of the ancestor. Report root sensitivity; if alternatives disagree, leave the biological result unrooted.

### Support and uncertainty

* A single Mash+NJ tree has no branch support. For representatives only, MashTree's 100 seed-bootstrap or half-hash jackknife replicates can quantify **sketch-sampling stability** [5]. They do not sample sites, assemblies, or recombination histories.
* High-fidelity branches receive UFBoot (candidate 1,000 replicates) and SH-aLRT, plus sensitivity to reference, recombination masking, alignment missingness, outgroup, representative set, `k`, and `s`. Collapse weak/unstable branches into polytomies rather than forcing bifurcation.
* Near-zero branches, negative NJ branches (if emitted), tied nearest neighbours, placement ambiguity, and topological discordance are first-class outputs.

## 7. Non-circular host clades

Host clades must be frozen from host evidence before mapping prophages:

1. Start with recombination-filtered host core topology and Mash sensitivity, not prophage traits.
2. Candidate named clade: monophyletic on the validated host tree, at least 20 independent hosts (smaller clusters remain clusters, not primary comparison clades), UFBoot ≥95 and SH-aLRT ≥80 on its defining branch, stable under the preregistered reference/mask/root/representative sensitivity analyses, and separable for ≥95% of all eligible hosts by distances to multiple clade representatives.
3. Host phylogroup, MLST/cgMLST, geography, date, and source may annotate or audit a clade but do not override topology. Prophage content is never an input, tie-breaker, or relabeling criterion.
4. Ambiguous placements remain unassigned. Lock a versioned `host_clade_membership.tsv` with tree/alignment/sketch checksums before joining any prophage table.
5. If Mash and core-SNP boundaries disagree, defer to adequately supported core inference or report alternative partitions; do not select whichever gives the strongest phage association.

These thresholds are preregistered project decision rules, not universal biological constants. A sensitivity table must show conclusions under reasonable alternatives.

## 8. Later prophage mapping—separate objects

### Map traits onto the fixed host tree

After host clades are frozen, join by stable `assembly_id`/sequence checksum:

* per-host prophage count and callable sequence denominator;
* per-host binary/count presence of each prophage protein/gene cluster, module, or element family;
* “core/accessory phage component” status with explicit denominators and uncertainty; and
* detection method/version/threshold and missing/not-callable states distinct from true absence.

Display tip rings/heatmaps and compute clade-stratified prevalence with biological-host denominators and confidence intervals. For hypothesis testing, account for host relatedness (phylogenetic logistic/mixed models or tree-aware permutation) and test sensitivity to host-clade ambiguity and prophage-calling uncertainty. Reconstruct gains/losses only on the fixed host tree under an explicit trait model and label them hypotheses: horizontal transfer, loss, detection failure, and unsampled intermediates can produce the same pattern.

### Do not confuse host and phage histories

The host core tree estimates host clonal relationships. A prophage similarity tree/network estimates relationships among prophage sequences/components. Never insert prophage sequences as host tips or use their clusters to define host clades.

Prophages are mosaics shaped by recombination, module exchange, HGT, integration/excision, and different gene histories. One bifurcating whole-prophage tree can be actively misleading. Use multiple complementary analyses: gene-content/Jaccard networks, shared-protein community graphs, synteny/module comparisons, recombination-aware alignments of conserved regions, and separate marker-gene trees where homologous. Compare phage network communities/marker trees to the host tree as an association/cophylogeny question, allowing reticulation and discordance rather than forcing one ancestry.

## 9. Explicit pilot-to-production gates

No stage advances on elapsed time alone.

| Stage | Future permitted workload | GO criteria | NO-GO / response |
|---|---|---|---|
| 0: tool/input smoke | ≤20 ordinary tiny FASTA assemblies; no BGZF claim | Pinned executable versions/help captured; one file → one tip; exact expected pair count; `dist`/`triangle` symmetry; Mash lower-PHYLIP parsed by pinned RapidNJ; outputs <100 MiB | Any parser/ID mismatch, or unresolved PanSN/BGZF task → stop and repair/pin |
| 1: QC/sensitivity | Preregistered small diversity panel, then 200-host pilot | ≥95% candidate ingroup assemblies pass completeness/contamination/species rules; `s=10k` vs `50k` Spearman ≥0.99 and nearest-neighbour agreement ≥95%; no unexplained p-value/shared-hash failures | Quarantine bad assemblies; increase `s`, revise `k`, or abandon Mash boundary at affected resolution |
| 2: scale pilot | 2k, then 5k accepted hosts; `/usr/bin/time -v`, `wc -c`, checksums | Measured pair count correct; no swapping/OOM; deterministic rerun checksums/topology conditional on pinned seed/order; lower 95% throughput bound predicts full distance within approved wall time; measured scratch projection + 30% reserve ≤70% of quota; tree peak RSS +100% safety margin fits allocation | Use representative/blockwise fallback; do not extrapolate from ≤20 probe or from theoretical rates |
| 3: high-fidelity validation | Host-only representative/boundary sets, per limited-diversity lineage | ≥90% reference callable in ≥95% samples; mean missing ≤5%; ≥100 non-recombinant informative SNPs; convergence; defining clades meet support/stability rules; ≥95% non-representatives confidently assign | Change lineage/reference/sampling; collapse unsupported clades; report ambiguity |
| 4: full overview | 24k–26k only after scheduler approval | Versioned manifest frozen; `P(n)` and tip counts exact; estimated wall/RAM/scratch/I/O approved; checkpointable blocks; independent output validation; no phage inputs | Stop on quota, I/O, pair-count, ID, symmetry, or QC failure; never silently drop hosts |
| 5: phage join | Trait tables only; no host-tree re-fitting from phage | Host memberships/checksums frozen first; callable/missing denominators explicit; host and phage results stored separately | Any clade definition influenced by phage content → invalidate and repeat blinded host definition |

Exact production caps (wall time, cores, RAM, scratch) must be filled from the actual scheduler allocation. The proposed relative margins prevent this design from pretending the captured login host is the allocation.

## 10. Audit checklist and exclusions

For every approved run retain: command line, stdout/stderr, wall/user/system time, peak RSS, filesystem and output bytes, input count/bytes, manifest and executable/container digests, random seeds, thread count, exit status, pair/tip counts, and SHA-256 for sketch/matrix/tree/membership outputs. A performance result applies only to its recorded `n,k,s`, input quality, tool build, cores, and filesystem.

This task performed **no** genome download, accession workload, BGZF conversion, production Mash comparison, tree construction, pangenome/IMPG computation, or full-table audit. It created only the four owned files. The optional `benchmark/` directory is absent because no benchmark was run; the evidence script explicitly reports zero inputs and outputs.

## Sources

All web/source material accessed 2026-07-24. Git commit links are immutable and version-applicable where a version is named.

1. Mash v2.3 `triangle` definition, output modes, and all-pair loop, commit `f228b9d9`: <https://github.com/marbl/Mash/blob/f228b9d9fc8e0f64a468d1deddb6ab9d6ac51abc/src/mash/CommandTriangle.cpp#L27-L43> and <https://github.com/marbl/Mash/blob/f228b9d9fc8e0f64a468d1deddb6ab9d6ac51abc/src/mash/CommandTriangle.cpp#L123-L213>.
2. Mash v2.3 versioned sketch documentation (specificity/sensitivity, hash width, sketch error/cost), commit `f228b9d9`: <https://github.com/marbl/Mash/blob/f228b9d9fc8e0f64a468d1deddb6ab9d6ac51abc/doc/sphinx/sketches.rst#L4-L64>.
3. Mash v2.3 option defaults and whole-file/`-i` behavior, commit `f228b9d9`: <https://github.com/marbl/Mash/blob/f228b9d9fc8e0f64a468d1deddb6ab9d6ac51abc/src/mash/Command.cpp#L168-L187> and release <https://github.com/marbl/Mash/releases/tag/v2.3>.
4. MashTree v1.4.6 README/defaults, commit `c0853a86`: <https://github.com/lskatz/mashtree/blob/c0853a86cf52ee4e47c83c6194c747cb9b9dbf5f/README.md#L12-L21> and <https://github.com/lskatz/mashtree/blob/c0853a86cf52ee4e47c83c6194c747cb9b9dbf5f/README.md#L61-L91>.
5. MashTree v1.4.6 algorithm and support semantics, including “does not infer phylogeny,” commit `c0853a86`: <https://github.com/lskatz/mashtree/blob/c0853a86cf52ee4e47c83c6194c747cb9b9dbf5f/docs/ALGORITHM.md#L21-L48>; Katz et al. 2019 JOSS: <https://doi.org/10.21105/joss.01762>.
6. MashTree v1.4.6 matrix materialization/tree implementation, commit `c0853a86`: <https://github.com/lskatz/mashtree/blob/c0853a86cf52ee4e47c83c6194c747cb9b9dbf5f/bin/mashtree#L350-L412> and <https://github.com/lskatz/mashtree/blob/c0853a86cf52ee4e47c83c6194c747cb9b9dbf5f/lib/Mashtree.pm#L307-L346>.
7. RapidNJ v2.3.2 README and documented external-memory help, commit `ed2d36e2`: <https://github.com/somme89/rapidNJ/blob/ed2d36e219d9db16778b941b5054c0fd021b528a/README#L1-L20> and <https://github.com/somme89/rapidNJ/blob/ed2d36e219d9db16778b941b5054c0fd021b528a/src/main.cpp#L535-L567>; Simonsen et al. large-NJ paper: <https://users-birc.au.dk/cstorm/software/rapidnj/papers/SimonsenOthers2011_CCIC.pdf>.
8. skani documented v0.3.2 source (ANI, aligned fractions, database search), commit `55d1dfd4`: <https://github.com/bluenote-1577/skani/blob/55d1dfd4af62ca5490f475b5596a5a7da5a8bba1/README.md#L3-L20> and <https://github.com/bluenote-1577/skani/blob/55d1dfd4af62ca5490f475b5596a5a7da5a8bba1/README.md#L69-L91>.
9. Jain et al. 2018, FastANI and ANI species-scale use/limitations: <https://pmc.ncbi.nlm.nih.gov/articles/PMC6269478/>.
10. Ondov et al. 2016, Mash formulation, accuracy, and ANI comparison: <https://pmc.ncbi.nlm.nih.gov/articles/PMC4915045/>; versioned formula source: <https://github.com/marbl/Mash/blob/f228b9d9fc8e0f64a468d1deddb6ab9d6ac51abc/doc/sphinx/distances.rst>.
11. QUAST 5.3.0 manual (assembly metrics): <https://quast.sourceforge.net/docs/manual.html>.
12. Chklovski et al. 2023, CheckM2 scope: <https://doi.org/10.1038/s41592-023-01940-w>; project documentation: <https://github.com/chklovski/CheckM2>.
13. Gubbins v3.4.3 manual, limited-diversity scope, quadratic scaling, recombination/error caveats, and tree options, commit `82ce42aa`: <https://github.com/nickjcroucher/gubbins/blob/82ce42aa657cb0d3d051c493826fb05ee9eac764/docs/gubbins_manual.md#L1-L27> and <https://github.com/nickjcroucher/gubbins/blob/82ce42aa657cb0d3d051c493826fb05ee9eac764/docs/gubbins_manual.md#L77-L146>; Croucher et al. 2015: <https://doi.org/10.1093/nar/gku1196>.
14. Snippy v4.6.0 release and version-matched assembly/core-alignment documentation, commit `86b64d09`: <https://github.com/tseemann/snippy/releases/tag/v4.6.0> and <https://github.com/tseemann/snippy/blob/86b64d09168105fe2126e65384c5cfa87448e6b6/README.md#L355-L385>.
15. IQ-TREE v3.1.2 release source, UFBoot/model features: <https://github.com/iqtree/iqtree3/tree/v3.1.2>; official ascertainment-bias documentation: <https://iqtree.github.io/doc/Substitution-Models#ascertainment-bias-correction>; official command reference: <https://iqtree.github.io/doc/Command-Reference>.
