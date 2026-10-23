# IS-mapping

Pipeline for detecting insertions of every IS element (insertion sequence) and
transposon in user-supplied databases across a large collection of *E. coli* short-read
SRA accessions, and classifying whether each insertion disrupted a gene.

At a high level: download reads for an accession, align them to each configured
element database (by default ISfinder + ISOSDB, and TnCentral), extract and
cluster the sequence "overhangs" flanking every IS hit, BLAST those clusters against
masked (IS-excised) *E. coli* reference genomes to classify whether the insertion landed
inside a gene, and pair up left/right overhang clusters that correspond to the same
insertion junction.

## Running the pipeline

The core per-accession pipeline is a Snakemake workflow:

```
snakemake --profile slurmprofile --rerun-incomplete --use-conda --executor slurm
```

`Snakefile` chains together: `index_database` + `reference_table` (once per database)
→ `prefetch` → `fasterq` → `bwa_isfinder` → `remove_fastqs`
→ `clip_and_cluster` (overhang extraction + position-anchored edit-distance clustering,
see `scripts/cluster_overhangs_edlib.py`) → `blast_clusters_to_ref`
(gene-disruption classification) → `add_pairing_column` (left/right junction pairing)
→ `nested_is` + `summarize_nested_is` (IS elements inserted inside other elements, see
[Nested IS insertions](#nested-is-insertions)).

### Configuration

You never need to edit the Snakefile. Every setting lives in one of two files:

| file | what it sets |
|---|---|
| [`config/config.yaml`](config/config.yaml) | What to run on, where outputs go, the reference genomes, every threshold, and software paths. Each setting has a comment explaining it. |
| [`config/databases.tsv`](config/databases.tsv) | The IS / transposon databases that reads are aligned to. |

Relative paths in either file are relative to the repository root.

#### What to run on

Set **one** of these three in `config.yaml`:

| setting | use it for | what it points to |
|---|---|---|
| `accessions:` | SRA runs, downloaded automatically | a text file with one run accession per line (e.g. `SRR8275060`) |
| `samples:` | your own sequencing reads | a tab-separated table with columns `sample`, `fastq_1` and `fastq_2` (leave `fastq_2` empty for single-end). See [`config/samples.example.tsv`](config/samples.example.tsv). |
| `existing_overhangs:` | rerunning from clustering onward on a finished run | a folder of `{sample}.tar.gz` overhang archives from an earlier run. Nothing is downloaded or aligned, and the archives are only read. |

If more than one is set, `existing_overhangs` wins over `samples`, which wins over
`accessions`. Leave the ones you aren't using as `""`.

#### Databases

Each line of `databases.tsv` is one database:

```
label	fasta	description
IS	/path/to/is_elements.fa	ISfinder + ISOSDB
Tn	/path/to/transposons.fa	TnCentral
```

- **`label`**: a short name (letters, digits, `_` or `-`). It tags every read and
  overhang, so results can be split by database later.
- **`fasta`**: the element sequences.
- **`description`** is optional.

Sequence names must be unique across all databases. See
[Element databases](#element-databases) for how they are used.

#### Changing settings for one run

You don't have to edit `config.yaml` itself:

- **Several settings:** write a small YAML file with only the settings you want to
  change and pass it with `--configfile my_run.yaml`. Anything you leave out keeps
  its value from `config.yaml`. For example,
  [`config/ecoli_summer.yaml`](config/ecoli_summer.yaml) only sets
  `existing_overhangs` and `output_dir`.
- **One setting:** override it on the command line, e.g.
  `--config accessions=my_list.txt`.

#### Checks before running

The pipeline checks the configuration before it starts any job. It then lists every
problem at once: missing files, bad database labels or sample names, duplicate
names, or thresholds that aren't numbers.

### Element databases

- `index_database` checks each FASTA (no headers glued onto sequence lines, no empty or
  repeated names), copies it into `index_dir`, and builds its minibwa index there.
- `reference_table` writes `index_dir/references.tsv` (`reference`, `database`,
  `length`) and fails if any reference name appears in more than one database, since
  every downstream output is keyed by reference name alone.
- `bwa_isfinder` maps each sample to each database separately, so a read's MAPQ only
  reflects its own database, and merges the results. Each alignment carries its
  database label as its read group.

To separate evidence by database (e.g. IS vs Tn):
- **reads**: `samtools view -r Tn sample.bam`
- **overhangs**: the `database` column of `overhangs/{id}/{id}.manifest.tsv`
- **clusters / gene-disruption tables**: join `is_element` onto `references.tsv`

Everything else in the repo is analysis run on top of that pipeline's output — pooling
results across accessions, plotting, summarizing, or re-running a stage at a different
scope (e.g. per IS element instead of per accession).

## Nested IS insertions

The last pipeline step looks for IS elements inserted *into other IS elements /
transposons*. It starts from the overhang clusters that
`blast_clusters_to_ref` marks `hit_type == no_hit` (the flank matched none of the
reference genomes). The reference genomes are IS-masked, so a flank that is itself IS
sequence ends up as `no_hit`. Each such flank is BLASTed against the same IS + Tn
databases the reads were aligned to. If the flank of host element X lands inside
element Y, X may sit inside Y.

**Orientation.** Same convention as `blast_clusters_to_ref.py`: `rep_seq` is in X's
forward orientation. The junction is the last base of a `left` overhang and the first
base of a `right` one. So the insertion point inside Y is BLAST's `send` (left) or
`sstart` (right), in 1-based Y coordinates.

**Rules.**
- `nested_is_db` (once): joins the `index_database` copies of every
  `config/databases.tsv` FASTA into `index_dir/nested_is/IS_Tn` and runs
  `makeblastdb`. The nested search then always uses the databases the reads were
  aligned to.
- `nested_is` (per sample, `scripts/nested_is.py`): reads the sample's
  `gene_disruption_paired/{id}.zip` and keeps the `no_hit` rows.
  - It drops `rep_seq`s under `min_len` and BLASTs each *exactly* identical sequence
    once. Near-identical sequences stay separate, and nothing is merged at the cluster
    level.
  - The search is `-task megablast -word_size 16 -evalue 1e-5 -max_target_seqs 50
    -max_hsps 1`. On a test set, word size 16 found as many hits as full `blastn`
    (766 vs 768 queries) at about 1/5 the runtime.
  - Hits are joined back onto every cluster, which is then classified (below). The
    class depends on the cluster's own `side` and `is_element`, not only on its
    sequence.
  - About 2 s per sample. These jobs are bundled 10 per SLURM submission
    (`group-components: nested_is` in `slurmprofile/config.yaml`).
- `summarize_nested_is` (once, `scripts/summarize_nested_is.py`): combines every
  sample's table into the two summaries below.

Thresholds are in the `nested_is:` section of `config/config.yaml`. Defaults: hits
must have ≥ 90% identity (`min_pident`) and ≥ 30 bp aligned (`min_aln`). A hit
"reaches the junction" if its alignment ends within 5 bp (`junction_tol`) of the
overhang's junction end. Classes, using the best-bitscore reaching hit:

| class | meaning |
|---|---|
| `internal` | junction lands inside Y: **nested-insertion candidate** |
| `terminal` | junction within 20 bp of either end of Y: X sits next to Y's end rather than interrupting Y |
| `self` | best hit (or a tie for it) is X itself: tandem copies, IS circles/dimers, or X sequence let in by `boundary_tolerance` |
| `distal_only` | IS/Tn hits exist but none reaches the junction: host DNA lies between X and Y, so X is near Y, not in it |
| `no_is_hit` | no passing IS/Tn hit: likely host DNA missing from the 6 references |
| `too_short` | `rep_seq` under 30 bp, not BLASTed |

`paired_in_target` is the strongest evidence. It is set when the left and right
clusters of the same host in the same sample land in the same Y, on the same strand,
with junction positions ≤ 20 bp apart: both flanks of one X copy point into one site
of Y, offset by its target-site duplication. Targets hit by only one side, at a single
site, in a large share of samples are often not nesting. The host's database entry
may stop short of the real element, so its flank runs on into a related element.
Check these, and close relatives such as ISEc83/ISEc84 vs ISEc8, by hand.

**Outputs** (under `output_dir`):
- `nested_is/{id}.tsv.gz`: one row per no_hit cluster, with the class and the
  target (`target`, `target_db`, `target_len`). It also has the junction position
  inside the target (`junction_pos`, `dist_to_target_end`), `sstrand` and
  `paired_in_target`. Alignment stats are `pident`, `aln_len`, `qcov`, `bitscore`,
  `evalue`, plus `n_tied`/`tied_targets` for equal-bitscore hits.
  Header only if the sample has no no_hit clusters.
- `nested_is_summary/class_summary.tsv`: clusters and samples per class.
- `nested_is_summary/host_target_summary.tsv`: one row per (host, target, class) for
  `internal`/`terminal` clusters, with `n_samples`, `n_samples_paired` and
  `n_left`/`n_right`. It also has `n_sites` (distinct junction positions, grouped
  where ≤ `pair_tol` apart), `top_site` and `median_pident`. Sorted by
  `n_samples_paired`. Start here.

**Running it outside the pipeline.** Use this for a finished run whose intermediates
were archived, where Snakemake would want to rebuild everything upstream.
`nested_is.py --zip_dir` reads every zip in a folder and BLASTs all samples at once,
deduplicating sequences across samples:

```
python scripts/nested_is.py --zip_dir <output_dir>/gene_disruption_paired \
    --blastdb <index_dir>/nested_is/IS_Tn --references <index_dir>/references.tsv \
    --output nested_is.clusters.tsv.gz --threads 4
python scripts/summarize_nested_is.py --tables nested_is.clusters.tsv.gz --output_dir .
```

Use the `bakta` env's Python (pandas). On ecoli_summer (148,743 samples, 7.3M no_hit
clusters, 3.3M unique sequences) this took ~45 min on one core. Those results are in
`/n/scratch/users/h/hua575/ecoli_summer/nested_is/`.

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
| `nested_is.py` | BLASTs one sample's (or a folder of samples') `no_hit` overhang clusters against the IS/Tn databases and classifies nested insertions. Used by the `nested_is` rule. |
| `summarize_nested_is.py` | Combines `nested_is.py` tables into class and host→target summaries. Used by `summarize_nested_is`. |
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

- **Output every passing BLAST hit, not just the top one, per (cluster, ref_genome).**
  `blast_clusters_to_ref.py`'s `best_hit_per_genome()` currently sorts by bitscore and
  keeps only the single best hit per `(qseqid, ref_genome)`, silently dropping any
  other hit that also clears the thresholds — which is exactly what hides
  repeated/duplicated-gene ambiguity (see the MAPQ to-do below). Add a minimum-bitscore
  threshold alongside the existing `MIN_PIDENT`/`MIN_QCOV` filters, and keep every hit
  that passes all three instead of collapsing to one, so downstream analysis can see
  and reason about multiple candidate placements for the same overhang.
- **Resolve ambiguous placement in repeated reference regions.** When an overhang
  cluster's representative matches several copies of a repeat, its placement is
  ambiguous, and a left/right pair can be split across two different copies so the
  insertion never pairs. ISMapper's figure 2 handles the analogous case by anchoring on
  the flank that passed its read-depth cutoff and re-searching for a sub-threshold
  partner near it (`closestBed` of the depth-filtered flank against the *unfiltered*
  opposite side — see `create_bed_files` / `check_unpaired_hits` in the ISMapper
  source), stamping everything rescued that way with a `?` confidence suffix. Its exact
  failure mode does not apply here: ISMapper loses depth because bwa scatters
  multi-mapping reads across repeat copies, whereas this pipeline clusters overhangs by
  *sequence* before placing them, so a cluster stays intact and full-depth however
  repetitive the locus. The remedy still adapts — let the confident partner choose which
  repeat copy the ambiguous flank is assigned to, and record that the pairing rather
  than the alignment made that call. Blocked on the multi-hit output above:
  `best_hit_per_genome()` discards the alternative copies before
  `pairing/add_pairing_column.py` ever sees them. Note MAPQ is not the lever here — the
  placement is a BLAST hit, not a bwa alignment — a bitscore-ratio window over the
  retained candidates is.
- **Aggregate the output data.** Current output format is hard to parse. Look into
  concatenating across accessions without producing absurdly large files, and
  cross-referencing the cluster info, original reads/overhang info, and gene hits
  together rather than as separate outputs.
- **Broaden the reference genome set.** Align against a phylogenetically diverse set
  of reference genomes (e.g. ~100) and compare what's conserved vs. not, instead of
  just the current 6.
- **Add the ISfinder database FASTA to the repo.** The IS-element reference used for
  alignment/masking (`ISfinder_database-master_2026/IS.database.collapsed99.fa`)
  currently only lives on O2 — add it to GitHub so the exact reference the
  pipeline runs against is version-controlled and reproducible outside the cluster.
- **Decide how to handle single vs. paired-end reads beyond alignment.** Single-vs-
  paired status only matters for `fasterq`/`bwa_isfinder`/`remove_fastqs` (how many
  fastqs to hand `minibwa`, which files to delete). After that, every mate of a
  paired-end read is treated as fully independent evidence — nothing uses mate
  orientation, insert size, or requires both mates to support the same junction.
  Worth deciding whether pairing should factor in as a corroborating signal,
  especially for disambiguating low-confidence placements (ties into the MAPQ
  to-do above).
- **Biological follow-ups to look into:**
  - Megaplate data has a region mediated by transposons that causes a duplication.
  - Cluster size correlates with coverage, which gives info about duplications.
  - Genomic rearrangements: regions close together that look like inversion events.
