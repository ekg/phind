# Provenance — raw source artifacts

| key | retrieved_at (UTC) | bytes | sha256 | url |
|---|---|---:|---|---|
| `inphared_genomes_fa` | 2026-08-17T16:49:44Z | 714521461 | `331cb0487b5ed60863faab38c6be1496e5ec015a5c0f787dd9ce865212f110c9` | https://s3.climb.ac.uk/millardlab-inphared/2026/7Apr2026_genomes.fa.gz |
| `inphared_table` | 2026-08-17T16:48:25Z | 681798 | `312ef901680f014165597bf740cb7b3b87e1d52383f3364c7a07193ea813ba2f` | https://s3.climb.ac.uk/millardlab-inphared/2026/7Apr2026_millardlab_website_table.txt.gz |
| `ncbi_backfill` | 2026-08-17T17:22:50Z | 7221334 | `f5cbc561d63505374f45fc6dfe7594f5bf481f30e7ea579ad51132f60d159c0b` | https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary+efetch by accession id |
| `ncbi_efetch` | 2026-08-17T17:01:39Z | 146480812 | `9cbe89243f89670e91f55ad4ed119ee7b52eedde1dcfd88351724e29c712aa01` | https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&rettype=fasta |
| `ncbi_esearch` | 2026-08-17T17:01:39Z | 797 | `2a31c0374ffc156b5b777b14d8c6176d415a5d535067f92e05e19f6082c3ed0e` | https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nucleotide&term=…&usehistory=y |
| `ncbi_esummary` | 2026-08-17T17:01:39Z | 397964 | `faa7b1bef0c5326c7b7eb981bd75bd5e51e76ae1b7db7a157b9b249888300f08` | https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nuccore |
| `ncbi_gb_remainder` | 2026-08-17T17:51:21Z | 977354571 | `c6d50fd4131183caf662f93fead059dca902e788f0d87a0f9b16a3bae2c59275` | https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&rettype=gb&retmode=text |
| `phagesdb_bulk_fasta` | 2026-08-17T16:49:31Z | 375899322 | `2a612899f24120310bbc66695061875c39f83f6c48677946fdd8fbc22285eb52` | https://phagesdb.org/media/Actinobacteriophages-All.fasta |
| `phagesdb_metadata` | 2026-08-17T16:48:24Z | 2042851 | `e367ed168649b0430fe4d8ba8951459b5285b9b96a331d30df6bb3c40d0c4f5f` | https://phagesdb.org/data/?set=seq&type=full |
| `taxonomy_hosts` | 2026-08-17T17:01:45Z | 6653 | `55d8063ca5357b4b56d2b1602fe1a30c2100eb8c3410b7dbb5dcceadf1487e4a` | https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy + efetch.fcgi?db=taxonomy |

Cached raw responses live under `cache/` (external storage). Per-batch NCBI efetch/esummary files are siblings of the recorded anchor file.