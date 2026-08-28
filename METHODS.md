# Technical methods

## Analysis environment

| Component | Version or configuration |
| --- | --- |
| Dorado | 1.3.2+2a7f935 |
| Basecalling model | `dna_r10.4.1_e8.2_400bps_sup@v5.2.0` |
| minimap2 | 2.24-r1122 |
| samtools / htslib | 1.13 / 1.13 |
| Kraken2 | 2.17.1; k=35; minimizer length=31; confidence=0.10; minimum hit groups=2 |
| BLASTN | 2.12.0+ |
| MAFFT | 7.520 |
| Python | 3.12 |

Dorado basecalling used the SUP v5.2.0 R10.4.1 model with `--min-qscore 0`, read splitting disabled, trimming disabled, and move information emitted. Strict demultiplexing used the Native Barcoding Kit 24 V14 (`SQK-NBD114.24`; Dorado identifier `SQK-NBD114-24`) with `dorado demux --barcode-both-ends --no-trim --emit-summary`. Primary-filter reads required mean Q score ≥10 and length ≥500 bp.

Competitive mapping used minimap2 with `-x map-ont --secondary=yes -N 20`. Hits were ordered by alignment score, MAPQ, matching bases, aligned length, target name, and target start. The alignment score was taken from the PAF `AS` tag when present and from matching bases otherwise. Qualifying evidence required MAPQ ≥20, aligned length ≥500 bp, identity ≥0.85, query coverage ≥0.80, and a positive best-genus score margin.

## Figure 2: detection and positive agreement

`analysis_inputs/cohort_target_reads.tsv.gz` contains read-level competitive assignments, acquisition times, quality values, target coordinates, alignment metrics, score margins, and nuclear or kinetoplast classification. `sample_metadata.tsv` contains specimen labels and PCR status, and `run_metadata.tsv` contains the evaluable run duration.

A locus was defined by the best target and the 10-kb bin containing its 0-based alignment start. A specimen was sequencing-positive when at least three qualifying nuclear reads represented at least three loci and two chromosomes. Detection time was the acquisition time of the first chronological read prefix satisfying this rule.

Positive percent agreement was calculated among the 15 PCR-positive biopsies. The 95% confidence interval is the two-sided Clopper–Pearson exact binomial interval.

## Figure 3: adaptive sampling

Target enrichment used a combined *L. major* ASM272v2 (`GCA_000002725.2`) and *L. tropica* CDC_Ltropica216-162_v1 (`GCA_014139745.1`) reference. Host depletion used GRCh38 (`GCA_000001405.15`).

`matched_target_reads.tsv.gz` contains read-level competitive evidence for three standard libraries and their target-enrichment counterparts. Analyses use nuclear reads acquired within the first 86,400 s. Target bases are the sum of complete qualifying read lengths. Target fraction is target bases divided by all primary-filter bases in `primary_24h_totals.tsv`. Fold enrichment compares enrichment and standard target fractions; absolute recovery compares their qualifying target bases.

Detection time uses the Figure 2 three-read, three-locus, two-chromosome rule. `run_yield_checkpoints.tsv` contains cumulative primary-filter bases at 0.25, 0.5, 1, 2, 4, 8, 12, and 24 h.

`adaptive_decision_counts.tsv` links terminal MinKNOW decisions to offline competitive classes. Outcomes are grouped as follows:

- `unblock`: Ejected
- `stop_receiving` and `no_decision`: Continued
- absent from the terminal decision table: No log

Panel E percentages use all offline-classified reads within each run mode and origin.

## Figure 4: species classification

Species composition in the biopsy samples is recalculated from `cohort_target_reads.tsv.gz`. Kraken2 species composition is derived from `cohort_kraken_reports.tsv.gz`. The *Leishmania* genus clade count, taxid 5658, is the Kraken2 target total. Species counts use taxids 5664 (*L. major*), 5666 (*L. tropica*), and 5667 (*L. aethiopica*).

The held-out evaluation uses 3,000-bp genome fragments sampled every 10,000 bp. `holdout.paf.gz` and `holdout_target_manifest.tsv.gz` reproduce competitive classification using the primary alignment thresholds. `holdout.kraken.output.gz` contains per-fragment Kraken assignments, and `holdout_taxonomy.tsv` maps the nine target-species taxids. Percentages are calculated among species-level calls.

## Figure 5: marker-spanning reads

`public_marker_scores.tsv.gz` contains the best and runner-up species bit scores for each public read and marker. A read is resolved when:

```text
best_bitscore - second_bitscore >= max(15, 0.02 × best_bitscore)
```

Run-level calls use the unique plurality species among resolved reads. Reference relation and independent-accuracy eligibility are retained for every read-marker call.

`biopsy_marker_blast.tsv.gz` contains BLAST tabular hits from the study biopsy samples. Hits require identity ≥80% and the following minimum aligned spans:

| Marker | Minimum span (bp) |
| --- | ---: |
| ITS1 | 180 |
| ITS2 | 250 |
| Miniexon | 180 |
| HSP20 | 220 |
| HSP70 coding | 700 |
| HSP70 3-prime UTR | 350 |
| Cytochrome b | 650 |

Within each read and marker, the best hit for each species is ranked by bit score, aligned length, and identity. A read-level species call requires a positive bit-score margin over the runner-up. Each molecule contributes once per marker. The ITS2 input includes the separate strict-cohort ITS2 scan.

## Supplementary Figure S1: marker divergence

The analysis inputs are E-INS-i and G-INS-i marker reference alignments generated with MAFFT 7.520. E-INS-i is primary for ITS1, ITS2, miniexon, and HSP70 3-prime UTR; G-INS-i is primary for HSP20, HSP70 coding, and cytochrome b.

Pairwise comparisons use the shared homologous interval. Total differences are nucleotide substitutions plus one-sequence gap bases and are also reported per 100 comparable columns. Results are considered stable between alignment strategies when the totals differ by no more than two positions and the rates differ by no more than one difference per 100 columns. `benchmark_reference_alignments` contains one selected reference per species and marker for the public-dataset benchmark. `study_species_alignments` contains the multi-accession references for *L. major*, *L. tropica*, and *L. aethiopica*. Exact duplicate references are collapsed within marker and species before pairwise summaries are calculated.

## Supplementary Figure S2: cost model

`cost_inputs.tsv` contains the price snapshot, pack sizes, consumable allowances, platform capacities, endpoint-PCR reference costs, and batch-size scenarios used in the calculations. Nanopore cost scenarios use full-kit SQK-NBD114.96 allocations.

## Reproduction

```bash
python -m pip install -r requirements.txt
for script in figure_{2,3,4,5}/reproduce.py figure_s{1,2}/reproduce.py; do
    python "$script"
done
```

Each script writes regenerated tables to its local `reproduced_source_data` directory. The corresponding figure values are stored in `source_data`.
