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
(gene-disruption classification) → `add_pairing_column` (left/right junction pairing).
### Configuration

Everything you'd change between runs is in two files; the Snakefile itself
shouldn't need editing.

- **What to run on**: either SRA accessions (`accessions:`, one per line, downloaded
  automatically) or your own reads (`samples:`, a table of `sample`, `fastq_1`,
  `fastq_2`; see [`config/samples.example.tsv`](config/samples.example.tsv)).
- **[`config/config.yaml`](config/config.yaml)**: the accession list, output
  folder, reference genomes, filtering/clustering thresholds and software paths. Each
  setting is commented.
- **[`config/databases.tsv`](config/databases.tsv)**: the IS / transposon
  databases to align reads to, one per line:

  ```
  label	fasta	description
  IS	/path/to/is_elements.fa	ISfinder + ISOSDB
  Tn	/path/to/transposons.fa	TnCentral
  ```

To run with different settings, either edit those files or copy `config.yaml`,
change what you need, and pass it with `--configfile my_config.yaml` (anything
you leave out keeps its default). One-off overrides work on the command line, e.g.
`--config accessions=my_list.txt`. Relative paths are relative to the repository root.
Before running anything, the pipeline checks the configuration and lists every
problem it finds (missing files, bad labels, non-numeric thresholds).

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

A downstream analysis, not yet a Snakemake rule, that looks for IS elements inserted
*into other IS elements / transposons*. It starts from the overhang clusters that
`blast_clusters_to_ref` marks `hit_type == no_hit` (the flank matched none of the
reference genomes). The reference genomes are IS-masked, so a flank that is itself IS
sequence ends up as `no_hit`. Each such flank is BLASTed against the same IS + Tn
databases the reads were aligned to. If the flank of host element X lands inside
element Y, X may sit inside Y.

**Orientation.** Same convention as `blast_clusters_to_ref.py`: `rep_seq` is in X's
forward orientation. The junction is the last base of a `left` overhang and the first
base of a `right` one. So the insertion point inside Y is BLAST's `send` (left) or
`sstart` (right), in 1-based Y coordinates.

**Steps** (scripts in `nested_is/scripts/`, see below):

| step | script | what it does |
|---|---|---|
| 1 | `extract_no_hit.py` | Pulls every `no_hit` row (sample, host element, side, cluster, pairing, `rep_seq`, length, `n_seqs`) out of all `gene_disruption_paired/{id}.zip` files. |
| 2 | `make_queries.py` | Drops `rep_seq`s under 30 bp and writes each *exactly* identical sequence once as a BLAST query (the same flank recurs across many samples). Near-identical sequences stay separate, and nothing is merged at the cluster level. |
| 3 | `blastn` | `-task megablast -word_size 16 -evalue 1e-5 -max_target_seqs 50 -max_hsps 1` against a BLAST db of the concatenated `config/databases.tsv` FASTAs. On a test set, word size 16 found as many hits as full `blastn` (766 vs 768 queries) at about 1/5 the runtime. |
| 4 | `classify_nested.py` | Joins hits back onto every cluster and classifies each one (below). The class depends on the cluster's own `side` and `is_element`, not only on its sequence. |
| — | `run_rest.sh` | Runs steps 2–4 in order. |

Hits must have ≥ 90% identity and ≥ 30 bp aligned. A hit "reaches the junction" if
its alignment ends within 5 bp of the overhang's junction end. Classes, using the
best-bitscore reaching hit:

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

**Outputs** (`--out_prefix nested_is`):
- `nested_is.clusters.tsv.gz`: one row per no_hit cluster, with the class and the
  target (`target`, `target_db`, `target_len`). It also has the junction position
  inside the target (`junction_pos`, `dist_to_target_end`), `sstrand` and
  `paired_in_target`. Alignment stats are `pident`, `aln_len`, `qcov`, `bitscore`,
  `evalue`, plus `n_tied`/`tied_targets` for equal-bitscore hits.
- `nested_is.class_summary.tsv`: clusters and samples per class.
- `nested_is.host_target_summary.tsv`: one row per (host, target, class) for
  `internal`/`terminal` clusters, with `n_samples`, `n_samples_paired` and
  `n_left`/`n_right`. It also has `n_sites` (distinct junction positions, grouped
  where ≤ 20 bp apart), `top_site` and `median_pident`. Start here.

**Running it.** Scripts take paths as arguments and use the `bakta` env's Python
(pandas). `blastn`/`makeblastdb` come from the `blast` env. From the output folder:

```
mkdir -p db
cat <IS fasta> <Tn fasta> > db/IS_Tn.fa                     # the FASTAs in config/databases.tsv
makeblastdb -in db/IS_Tn.fa -dbtype nucl -out db/IS_Tn -parse_seqids
# db/elements.tsv: name<TAB>database label<TAB>length, one line per element
python scripts/extract_no_hit.py <output_dir>/gene_disruption_paired no_hit_clusters.tsv.gz
./scripts/run_rest.sh                                        # THREADS=n to multithread BLAST
```

On the ecoli_summer run (148,743 samples, 7.3M no_hit clusters, 3.3M unique queries),
extraction took ~15 min and BLAST + classification ~30 min on one core. The first run
and its scripts are in `/n/scratch/users/h/hua575/ecoli_summer/nested_is/`.

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
- **Make nested-IS detection a Snakemake step.** Planned: a `build_nested_is_db`
  rule (concatenate `config/databases.tsv`, `makeblastdb` in `index_dir`), a
  per-sample `nested_is` rule after `add_pairing_column` (steps 1–4 above in one
  script, grouped with `gene_disruption`), and a `summarize_nested_is` rule writing the
  class and host→target summaries. Thresholds would go in a `nested_is:` config
  section. Per-sample BLAST loses cross-sample dedup (~2× the queries) but lets new
  samples run incrementally.
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
