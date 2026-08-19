# External evidence root for validate-logan-candidates

Bulky artifacts live at
`/mnt/nvme3n1/erikg/phind-genome-work/ntm/v2/pilot_validation/`:

| file | content |
|---|---|
| `tested_genomes.fa` / `.sha256` | 17 tested ML genomes + per-genome sha256 |
| `selected_baits.json` | 20 selected NTM baits traced to `bait_manifest.tsv` |
| `baits/ntm_baits.fa`, `baits/controls.fa` | bait/control query FASTAs |
| `tested_msh.msh`, `panel_i.msh` | mash 2.3 sketches (k=21, s=20000, `-i`); panel sketch 1.0 GB |
| `tested_vs_panel_i.mash` | all genome x 6,266-record distances |
| `refs/<genome>.refs.fa`, `refs/<genome>.query.fa` | per-genome top-5 (d<=0.2) public refs |
| `paf/*.paf` | minimap2 2.31-r1302 alignments (asm20 genome-level; sr+asm20 bait-level) |
| `genome_external_evidence.tsv`, `bait_external_evidence.tsv`, `genome_vs_public_alignment.tsv` | summarized evidence (repo copies here) |
| `MANIFEST.sha256` | sha256 + size of all 63 artifacts |

Commands: `mash sketch -s 20000 -i -o X <fasta>`; `mash dist -p 8 -v 1.0 -d 1.0`;
`minimap2 -x asm20|sr -c --eqx -t 8 <panel|refs> <query>`.
