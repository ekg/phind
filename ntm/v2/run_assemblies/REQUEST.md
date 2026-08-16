# REQUEST: NTM v2 run-assembly FASTAs (9,543 genomes)

**To:** Collaborator (Iridis cluster, `/iridisfs/scratch/ejh1g20/`)
**From:** Erik's NTM v2 pipeline (`erik@hypervolu.me`)
**Date:** 2026-08-16
**Status:** PENDING UPLOAD

## 1. What we need

The QC-passed NTM v2 cohort includes **9,543 run-assembly genomes**
(5,620 ERR + 3,915 SRR + 8 DRR) whose FASTA files exist only on the
collaborator's Iridis cluster. 5,907 of these are **prophage-bearing**
(63% of the 9,434 prophage-bearing genomes in v2); they are the top
priority — the prophage-extraction pipeline can start as soon as they
arrive.

Each genome is one file:

```
/iridisfs/scratch/ejh1g20/Mycobacteria/NTM_Db/assemblies/{ACCESSION}/contigs.fasta
```

- **Source glob:** `/iridisfs/scratch/ejh1g20/Mycobacteria/NTM_Db/assemblies/{ERR,SRR,DRR}*/contigs.fasta`
- **Expected file count:** 9,543 `contigs.fasta` files
- **Accession list (priority order):** `run_assemblies_needed.txt` in this
  directory — 9,543 accessions, one per line. **First 5,907 lines are the
  prophage-bearing accessions (highest priority); the remaining 3,636 are
  the rest.** Files may arrive in any order/batches — ingestion keys on the
  accession, not on file order — but prophage-bearing first unblocks the
  downstream prophage extraction earliest.

## 2. Destination

```
erik@hypervolu.me:~/www/phage/ntm_v2_run_assemblies/
```

SSH key auth to `erik@hypervolu.me` works. The receiving directory is
already the web root used by the pipeline; uploads land in
`~/www/phage/ntm_v2_run_assemblies/` and are automatically visible to the
ingestion watcher.

**Disk-space note:** the destination currently has **~7.9 GB free**.
Expected total for 9,543 mycobacterial assemblies is roughly 40–50 GB, so
**please upload in batches** (see §4) — do not attempt a single
full-cohort transfer. Batches are verified and fetched by ingestion as soon
as they land, so the destination never needs to hold more than one batch.

## 3. Checksum spec (required for every batch)

Each batch MUST include a checksum manifest so ingestion can verify
integrity before trusting the files. Format: standard `sha256sum -b`
output, one line per file, in the **same directory as the files**, named
`sha256sums.txt`:

```
<64 lowercase hex sha256> *<relative-path-from-batch-root>
```

Produce it with (from the assemblies dir, paths matching the destination
layout — see §4):

```bash
cd /iridisfs/scratch/ejh1g20/Mycobacteria/NTM_Db/assemblies
sha256sum -b {ERR,SRR,DRR}*/contigs.fasta > sha256sums.txt
```

(`-b` = binary mode; the `*` prefix marks binary mode so `sha256sum -c`
verifies on the receiving side without ambiguity.) Ingestion will run
`sha256sum -c sha256sums.txt` on arrival; a mismatch is treated as a
failed transfer and the file is quarantined, not ingested.

## 4. Transfer instructions

### Option A — rsync (preferred; resumable, incremental)

From Iridis, with the destination layout `ntm_v2_run_assemblies/{acc}/contigs.fasta`:

```bash
cd /iridisfs/scratch/ejh1g20/Mycobacteria/NTM_Db/assemblies

# canonical full-cohort glob (expands to all 9,543 files):
#   rsync -avz --partial {ERR,SRR,DRR}*/contigs.fasta \
#     erik@hypervolu.me:~/www/phage/ntm_v2_run_assemblies/
# ⚠️ full cohort ≈ 40–50 GB; destination has ~7.9 GB free — batch it:

# batch 1 (prophage-bearing, top priority): first 1,000 accessions
head -n 1000 /path/to/run_assemblies_needed.txt \
  | sed 's|$|/contigs.fasta|' > batch1_paths.txt
rsync -avz --partial --files-from=batch1_paths.txt \
  . erik@hypervolu.me:~/www/phage/ntm_v2_run_assemblies/

# batch 2: next 1,000 (lines 1001..2000) — same pattern with
#   sed -n '1001,2000p' run_assemblies_needed.txt | sed 's|$|/contigs.fasta|'
```

Notes:
- `--files-from` paths are relative to `assemblies/` (we `cd` there), so the
  receiving tree is exactly `ntm_v2_run_assemblies/{acc}/contigs.fasta`.
- `--partial` makes interrupted transfers resumable — re-running the same
  rsync resumes rather than restarts.
- Any batch order is fine (ingestion keys on the accession), but
  prophage-bearing accessions (first 5,907 lines) unblock extraction first.

### Checksums, every batch

Generate in the same directory on Iridis, before/after the transfer, over
paths that match the destination layout:

```bash
cd /iridisfs/scratch/ejh1g20/Mycobacteria/NTM_Db/assemblies
sha256sum -b {ERR,SRR,DRR}*/contigs.fasta > sha256sums.txt
rsync -avz sha256sums.txt erik@hypervolu.me:~/www/phage/ntm_v2_run_assemblies/
```

`sha256sums.txt` sits at the destination root next to the `{acc}/` dirs, one
line per file: `<64 hex> *{acc}/contigs.fasta`. Ingestion runs
`sha256sum -c` semantics on arrival; a mismatch quarantines the file.

### Option B — tar (if rsync is blocked by the destination's SSH config)

```bash
cd /iridisfs/scratch/ejh1g20/Mycobacteria/NTM_Db/assemblies
tar -czf ntm_v2_batch1.tar.gz $(head -n 1000 /path/to/run_assemblies_needed.txt | sed 's|$|/contigs.fasta|')
sha256sum ntm_v2_batch1.tar.gz > ntm_v2_batch1.tar.gz.sha256
# upload the two files; destination unpack:
#   tar -xzf ntm_v2_batch1.tar.gz -C ~/www/phage/ntm_v2_run_assemblies/
```

## 5. What ingestion does when files arrive

The watcher polls `~/www/phage/ntm_v2_run_assemblies/`, fetches new
batches to `/mnt/nvme3n1/erikg/phind-genome-work/ntm/v2/genomes/run_assemblies/`,
verifies:

1. `sha256sum -c sha256sums.txt` (transfer integrity),
2. per-genome **contig count** and **total length** against the v2 manifest
   columns `assembly_contigs` / `assembly_length_bp` (a strong content
   check: a mismatched file is not the assembly the manifest describes),
3. PanSN bgzip validity (`samtools faidx` on every converted file,
   header `{ACCESSION}#1#{contig}`).

Only files that pass all three are retained; anything else is itemized in
the ingestion report and quarantined.

## 6. Contact / questions

- **Questions about this request:** reply via the WG task `ntm-v2-run`
  inbox, or email erik@hypervolu.me.
- **If the full cohort cannot be uploaded:** upload at least the 5,907
  prophage-bearing accessions (first 5,907 lines of
  `run_assemblies_needed.txt`) — that covers 63% of prophage-bearing v2
  genomes. Anything less and we report partial coverage.

## 7. ENA note (probe finding, 2026-08-16)

An ENA portal probe (n=100, seed=42) found **79% of run accessions have
retrievable submitted `SEQUENCE_ASSEMBLY` analyses** (submitted_ftp
confirmed 200 OK). However, spot-checking 4 ENA assemblies against the
manifest showed **none match the collaborator assemblies** (e.g.
ERR3566243: ENA 48 contigs / 6,085,029 bp vs manifest 98 / 6,109,636 bp;
ERR330884: ENA 61 / 5,207,668 vs manifest 12 / 5,209,116). The manifest
prophage coordinates and `assembly_contigs`/`assembly_length_bp` were
computed on the collaborator's assemblies, so **ENA is not a substitute** —
it can only serve as a fallback for genomes whose collaborator files are
unavailable, and only with full re-verification/re-analysis. The
collaborator upload remains the authoritative source. Details:
`ena_probe_results.tsv` / `ena_probe_findings.md` in this directory.
