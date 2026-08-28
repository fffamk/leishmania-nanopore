# Leishmania nanopore manuscript: source data and analysis code

This repository provides source data and reproducible analysis code for Figures 2–5 and Supplementary Figures S1–S2. Each figure directory contains the figure source tables, analysis inputs, and a script that regenerates the reported results.

| Folder | Analysis |
| --- | --- |
| `figure_2` | Sequencing-positive calls, positive percent agreement, and detection time |
| `figure_3` | Adaptive enrichment, target recovery, throughput, and decision outcomes |
| `figure_4` | Cohort and held-out species classification |
| `figure_5` | Marker-spanning read analysis in public datasets and biopsy samples |
| `figure_s1` | Pairwise marker-sequence divergence |
| `figure_s2` | Batch-size and platform cost calculations |

## Reproduction

```bash
python -m pip install -r requirements.txt
for script in figure_{2,3,4,5}/reproduce.py figure_s{1,2}/reproduce.py; do
    python "$script"
done
```

Each script reads `analysis_inputs` and writes the regenerated tables to `reproduced_source_data`. The values shown in the figures are stored in `source_data`.

Technical methods, decision rules, reference details, and software versions are documented in [METHODS.md](METHODS.md).
