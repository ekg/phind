**Subject: NTM prophage expansion — request for your 18,055 run-assembly FASTAs**

Hi [collaborator],

Great news on our end: we've processed your BV-BRC phigaro QC-passed export (26,499 genomes, 47,994 prophage calls) and the PHIND pipeline is now running on the public portion — 16,243 genomes covering 9,446 prophages, with tight clades, alignments, and partitioning complete.

**The gap:** 18,055 genomes in the export are your SPAdes assemblies of SRA runs (ERR/SRR/DRR accessions). These aren't recoverable from any public source — BV-BRC serves them only from your private workspace, NCBI holds just the raw reads, and ENA's substitute assemblies don't match your contigs (so the phigaro coordinates can't map onto them).

**The ask:** a bulk export of those 18,055 assembly FASTAs (~90 Gbp gzipped estimate). Any of these works for us, easiest first:

1. A BV-BRC group share / read-only API token for us on your workspace, or
2. A tarball/rsync/OneDrive link to the FASTAs (gzipped, one file per assembly, named by run accession ideally), or
3. Your BV-BRC batch download of the group (they support it from the workspace UI)

Once we have them, the full pipeline re-runs on the complete 26.5k cohort — extraction through ML/ancestral phage genomes — and the per-run prophage coordinates you've already computed transfer directly, no re-calling needed.

For reference, the pending request from August (9,543 v2-era run assemblies) is subsumed by this — the 18,055 includes those.

Thanks — happy to jump on a call if the transfer needs coordinating.

[Your name]
