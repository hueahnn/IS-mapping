# IS-mapping

Pipeline for detecting insertions of every IS element (insertion sequence) in the
ISfinder database across a large collection of *E. coli* short-read SRA accessions, and
classifying whether each insertion disrupted a gene.

At a high level: download reads for an accession, align them to the full ISfinder
IS-element database (not just one IS element — every one it contains), extract and
cluster the sequence "overhangs" flanking every IS hit, BLAST those clusters against
masked (IS-excised) *E. coli* reference genomes to classify whether the insertion landed
inside a gene, and pair up left/right overhang clusters that correspond to the same
insertion junction.

## Running the pipeline

The core per-accession pipeline is a Snakemake workflow:

```
snakemake --profile slurmprofile --rerun-incomplete --use-conda --executor slurm
```

`Snakefile` chains together: `prefetch` → `fasterq` → `bwa_isfinder` → `remove_fastqs`
→ `clip_and_cluster` (overhang extraction + position-anchored edit-distance clustering,
see `scripts/cluster_overhangs_edlib.py`) → `blast_clusters_to_ref`
(gene-disruption classification) → `add_pairing_column` (left/right junction pairing).
It runs against `atb/ecoli_atb_sra_accessions.txt` by default (override with
`--config input_path=...`).

Everything else in the repo is analysis run on top of that pipeline's output — pooling
results across accessions, plotting, summarizing, or re-running a stage at a different
scope (e.g. per IS element instead of per accession).

## Repository layout

### `Snakefile`
The per-accession pipeline described above.

### `scripts/`
Scripts invoked by the Snakefile, plus downstream analysis run standalone.

| file | purpose |
|---|---|
| `overhangs.py` | 4-stage filter (coverage → per-read quality → per-IS-element hit count → FASTA extraction) that pulls soft-clipped overhang reads out of a BAM. Used by the `clip_and_cluster` rule. |
| `cluster_overhangs_edlib.py` | Clusters one sample's overhangs per (is_element, side) using a position-anchored edit distance (edlib SHW/prefix mode) instead of CD-HIT, writing `reads.tsv` directly. Used by `clip_and_cluster`. |
| `combine_cdhit_clusters.py` | Merges one sample's cd-hit cluster + representative-sequence output into a single `reads.tsv`, dropping singleton clusters. No longer used by `clip_and_cluster` (superseded by `cluster_overhangs_edlib.py`); kept for the still-CD-HIT-based `cluster_and_align_by_is.py`. |
| `blast_clusters_to_ref.py` | BLASTs cluster representatives against the masked reference genomes and classifies gene disruption at the insertion junction. Used by the `blast_clusters_to_ref` rule. |
| `masking.py` | Builds the masked reference genomes: BLASTs each *E. coli* reference against ISfinder and excises matching regions from both FASTA and GBFF. |
| `cluster_global_reps.py` | Pools every accession's cluster representatives (by side, across all IS elements) and re-clusters at the global level. |
| `cluster_and_align_by_is.py` | Same idea scoped to specific IS elements: pools raw overhangs for those elements across all accessions, re-clusters, BLASTs, and pairs junctions. |
| `pool_reps_min_size.py` | Filters existing per-accession cluster reps by minimum cluster size and merges them into one FASTA (no re-clustering). |
| `plot_left_right_overhangs_per_is.py` | One scatter plot per IS element: genome-by-genome left vs. right deduplicated overhang counts. |
| `summarize_gene_disruption_v2.py` | Streams every `gene_disruption/*.zip` and tallies cluster counts by `hit_type` and pairing status. |
| `build_atb_species_table.py` | For a given species, builds the accession list (ATB query → BioSample IDs → SRA coverage filter) that feeds the Snakefile's `ACCESSIONS_PATH`. |
| `reads.ipynb` | Exploratory notebook from early pipeline development (minibwa hit counts, gbff parsing, insertion-site BAM inspection). |

### `pairing/`
| file | purpose |
|---|---|
| `add_pairing_column.py` | Adds a `pairing` column to gene-disruption TSVs, matching left/right overhang clusters that land at the same insertion junction within a max gap. |
| `run_pairing.sbatch` | SLURM array job running `add_pairing_column.py` over chunked zip lists. |

### `clustering/`
| file | purpose |
|---|---|
| `filter_cdhit_reps.py` | Filters a cd-hit representative FASTA down to clusters larger than a minimum size, using the matching `.clstr` file. |

### `coverage_plots/`
| file | purpose |
|---|---|
| `align_to_masked_refs.sbatch` | Array job: BWA-aligns the merged, size-filtered cluster-rep FASTA against each of the 6 masked (v1) reference genomes. |
| `align_to_masked_refs_v2.sbatch` | Same, against the v2 (ISfinder-BLAST-masked) reference genomes. |

### `is_level/`
| file | purpose |
|---|---|
| `align_is_clusters_to_masked_refs_v2.sbatch` | Array job: BWA-aligns one IS element's pooled cluster reps (from `cluster_and_align_by_is.py`) against the 6 v2 masked reference genomes. |

### `is_coverage/`
| file | purpose |
|---|---|
| `extract_is_coverage.py` | Computes per-(genome, IS element) coverage stats (mean/variance depth, full-length and edge windows) from BAMs aligned to ISfinder. |
| `coverage.ipynb` | Notebook validating `extract_is_coverage.py`'s output. |

### `insertion_event_counts/`
| file | purpose |
|---|---|
| `count_is_events.py` | Counts distinct insertion events per IS element by treating paired left/right clusters as edges in a graph and counting connected components, avoiding double-counting a junction across sides or reference genomes. |

### `ecoli_genomes_metadata/`
| file | purpose |
|---|---|
| `sra_table.ipynb` | Notebook exploring the NCBI SRA accession table and filtering it down to *E. coli* BioSamples. |

### `analysis.ipynb`
Notebook of plots/numbers assembled for presenting pipeline results.

## Not tracked in this repo

This repo holds pipeline code only. It depends on data and infrastructure that live on
the HMS O2 cluster but aren't checked in:

- **Reference data**: `atb/`, `ref_genomes/`, `ISfinder_database-master_2026/` (IS
  database + masked reference genomes and BLAST dbs).
- **Conda environments**: paths under `/home/hua575/miniconda3/envs/` (`minibwa`,
  `cd-hit`, `blast`, `bakta`, `ncbi-datasets`, `ismapper`) and `.yaml` env specs under
  `/home/hua575/baymlab/sibmi/conda-envs/`.
- **Snakemake profile**: `slurmprofile/config.yaml` (SLURM resource/group-component
  config referenced by `--profile slurmprofile`).
- **Scratch outputs**: all pipeline outputs live under
  `/n/scratch/users/h/hua575/atb_filtered/` and are not versioned.

## To Dos for the Future

- **Implement new clustering.** CD-HIT is alignment-based and clusters on percent
  identity. Want to move to Levenshtein distance / Hamming distance, incorporating
  positional start information and a count-based difference metric instead of %ID.
  Look at the Slack messages with Jim on this. Consider the `block-aligner` package
  (Rust + C).
- **Use MAPQ for repeated genes.** When mapping overhang clusters to the reference
  genome, gene repetitions cause ambiguous placement — use MAPQ score to handle this.
  See ISMapper figure 2.
- **Aggregate the output data.** Current output format is hard to parse. Look into
  concatenating across accessions without producing absurdly large files, and
  cross-referencing the cd-hit info, original reads/overhang info, and gene hits
  together rather than as separate outputs.
- **Broaden the reference genome set.** Align against a phylogenetically diverse set
  of reference genomes (e.g. ~100) and compare what's conserved vs. not, instead of
  just the current 6.
- **Biological follow-ups to look into:**
  - Megaplate data has a region mediated by transposons that causes a duplication.
  - Cluster size correlates with coverage, which gives info about duplications.
  - Genomic rearrangements: regions close together that look like inversion events.
