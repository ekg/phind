# Full prophage FASTA artifact

`full_prophages.fa.gz` is the BGZF-compressed, Git-LFS-managed form of the full 132,393-record prophage FASTA.

## Integrity

| Form | Size | SHA-256 |
|---|---:|---|
| uncompressed `full_prophages.fa` | 3,265,823,570 bytes (3.1 GiB) | `ed85b2fb549be18bc638d8485f5b5add7c2d394f3822efe66a90ca6d979758d3` |
| BGZF `full_prophages.fa.gz` | 868,447,197 bytes (829 MiB) | `6461fa53e24a574bd7dacd2675b27b60c0b5664af589c5e2b2dc2afdc4724267` |
| BGZF index `full_prophages.fa.gz.gzi` | 1,063,848 bytes (1.0 MiB) | `3d956777b2dd57a7ec41f1fa70fbb1f894407c8d29e9ec4116dc14b6aa62351d` |

The compressed file was generated with:

```bash
bgzip -@ 16 -l 9 -c full_prophages.fa > full_prophages.fa.gz
bgzip -r full_prophages.fa.gz
```

BGZF is gzip-compatible, so tools that accept compressed FASTA can read it directly. The `.gzi` file supports indexed BGZF access.

## Materializing the historical path

Scripts that still require the uncompressed historical filename can run:

```bash
tools/materialize_full_prophages.sh
```

This creates the ignored working file `prophage_homology_survey/full_prophages.fa` and verifies its original SHA-256. The raw 3.1 GiB file is deliberately absent from rewritten Git history because GitHub Free/Pro rejects individual LFS objects larger than 2 GiB.
