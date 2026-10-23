# IS-element insertion pipeline: given a list of SRA accession IDs, downloads
# each genome's reads, aligns them to one or more IS-element / transposon
# databases, extracts and clusters the sequence "overhangs" flanking every
# element hit, BLASTs those clusters against masked E. coli reference genomes
# to classify whether the insertion disrupted a gene, and pairs up left/right
# overhang clusters that land at the same insertion junction.
# begin: 06/18/2026
# run: snakemake --profile slurmprofile --rerun-incomplete --use-conda --executor slurm
#
# All settings live in config/config.yaml (databases in config/databases.tsv);
# this file should not need editing to run on new data.

import csv
import os
import re

from snakemake.exceptions import WorkflowError

SCRIPT_DIR = workflow.basedir  # repository root; scripts and relative config paths resolve from here

configfile: os.path.join(SCRIPT_DIR, "config", "config.yaml")


# ---------------------------------------------------------------------------
# read and check the config -- every problem is collected and reported at once
# ---------------------------------------------------------------------------

_problems = []


def _resolve(path):
	"""Relative config paths are relative to the repository root."""
	path = os.path.expanduser(str(path))
	return path if os.path.isabs(path) else os.path.join(SCRIPT_DIR, path)


def _setting(*keys):
	"""config[keys[0]][keys[1]]..., recording a problem if any level is missing."""
	value = config
	for depth, key in enumerate(keys):
		if not isinstance(value, dict) or key not in value:
			_problems.append(f"missing setting '{': '.join(keys[:depth + 1])}' in config/config.yaml")
			return None
		value = value[key]
	return value


def _existing(path, what, suffix=""):
	"""Resolve a path setting and record a problem if the file is missing."""
	if path is None:
		return None
	path = _resolve(path)
	if not os.path.exists(path + suffix):
		_problems.append(f"{what} not found: {path + suffix}")
	return path


def _number(*keys, kind=float):
	value = _setting(*keys)
	if value is None:
		return None
	try:
		return kind(value)
	except (TypeError, ValueError):
		_problems.append(f"setting '{': '.join(keys)}' must be a number, got {value!r}")
		return None


def _read_databases(path):
	"""label -> fasta from databases.tsv (# comment lines and blank lines skipped)."""
	databases = {}
	if path is None or not os.path.exists(path):
		return databases
	with open(path) as fh:
		rows = csv.DictReader((l for l in fh if l.strip() and not l.startswith("#")), delimiter="\t")
		missing = {"label", "fasta"} - set(rows.fieldnames or [])
		if missing:
			_problems.append(f"{path}: header line must have columns 'label' and 'fasta' "
							 f"(tab-separated), missing {sorted(missing)}")
			return databases
		for n, row in enumerate(rows, start=1):
			label, fasta = (row["label"] or "").strip(), (row["fasta"] or "").strip()
			if not re.fullmatch(r"[A-Za-z0-9_-]+", label):
				_problems.append(f"{path}, database {n}: label {label!r} may only contain "
								 "letters, digits, '_' or '-'")
			elif label in databases:
				_problems.append(f"{path}: database label {label!r} is listed twice")
			elif not fasta:
				_problems.append(f"{path}: database {label!r} has no fasta path")
			else:
				databases[label] = _existing(fasta, f"FASTA for database {label!r}")
	if not databases and not any("databases" in p for p in _problems):
		_problems.append(f"{path}: no databases listed")
	return databases


def _read_samples(path):
	"""sample -> [fastq_1] or [fastq_1, fastq_2] from a samples table."""
	samples = {}
	if path is None or not os.path.exists(path):
		return samples
	with open(path) as fh:
		rows = csv.DictReader((l for l in fh if l.strip() and not l.startswith("#")), delimiter="\t")
		missing = {"sample", "fastq_1"} - set(rows.fieldnames or [])
		if missing:
			_problems.append(f"{path}: header line must have columns 'sample', 'fastq_1' and "
							 f"optionally 'fastq_2' (tab-separated), missing {sorted(missing)}")
			return samples
		for n, row in enumerate(rows, start=1):
			name = (row["sample"] or "").strip()
			fastqs = [(row.get(c) or "").strip() for c in ("fastq_1", "fastq_2")]
			fastqs = [f for f in fastqs if f]
			if not re.fullmatch(r"[A-Za-z0-9_-]+", name):
				_problems.append(f"{path}, sample {n}: name {name!r} may only contain "
								 "letters, digits, '_' or '-'")
			elif name in samples:
				_problems.append(f"{path}: sample {name!r} is listed twice")
			elif not fastqs:
				_problems.append(f"{path}: sample {name!r} has no fastq_1")
			else:
				samples[name] = [_existing(f, f"reads for sample {name!r}") for f in fastqs]
	if not samples and not any(path in p for p in _problems):
		_problems.append(f"{path}: no samples listed")
	return samples


# what to run on: overhang archives from an earlier run if given, else your own
# reads (samples table), else SRA accessions
EXISTING_OVERHANGS = (_existing(config["existing_overhangs"], "existing overhangs folder")
					  if config.get("existing_overhangs") else None)
SAMPLES = (_read_samples(_existing(_setting("samples"), "samples table"))
		   if config.get("samples") and not EXISTING_OVERHANGS else {})
ACCESSIONS_PATH = (None if EXISTING_OVERHANGS or SAMPLES or config.get("samples")
				   else _existing(_setting("accessions"), "accession list"))
DATABASES = _read_databases(_existing(_setting("databases"), "database table"))

# where to write. Every per-sample output lives under OUTPUT_DIR.
OUTPUT_DIR = _resolve(_setting("output_dir") or "")
# Kept out of OUTPUT_DIR on purpose: if the index were purged and rebuilt, its
# new mtime would make every existing bam look out of date (rerun-triggers is
# mtime-only, see slurmprofile/config.yaml) and re-queue every accession.
INDEX_DIR = _resolve(_setting("index_dir") or "")                        # index_database, reference_table, bwa_isfinder
REFERENCES_TSV = os.path.join(INDEX_DIR, "references.tsv")               # reference_table

# reference genomes / BLAST db for blast_clusters_to_ref (see scripts/masking.py)
REF_GENOMES_DIR = _existing(_setting("reference_genomes", "dir"), "reference genome folder")
REF_ACCESSIONS = _setting("reference_genomes", "accessions") or []
GBFF_SUFFIX = _setting("reference_genomes", "gbff_suffix")
BLASTDB = _existing(_setting("reference_genomes", "blast_db"), "BLAST database", suffix=".nsq")
REF_GBFFS = [
	_existing(os.path.join(REF_GENOMES_DIR or "", acc, f"{acc}_{GBFF_SUFFIX}.gbff"),
			  f"annotation for reference genome {acc}")
	for acc in REF_ACCESSIONS
]
REF_ACCESSIONS_STR = " ".join(REF_ACCESSIONS)  # shell {}-formatting can't eval " ".join(...) inline

# run_accession -> assembled genome length, for blast_clusters_to_ref's
# est_coverage denominator (see scripts/build_genome_size_table.py). Optional.
GENOME_SIZES = _setting("genome_sizes")
GENOME_SIZES = _existing(GENOME_SIZES, "genome size table") if GENOME_SIZES else None

# overhang extraction (scripts/overhangs.py)
MIN_ELEMENT_COVERAGE_DEPTH = _number("overhangs", "min_element_coverage_depth", kind=int)
MIN_ELEMENT_COVERAGE_FRACTION = _number("overhangs", "min_element_coverage_fraction")
MIN_MAPQ = _number("overhangs", "min_mapq", kind=int)
MIN_PIDENT = _number("overhangs", "min_pident")
MIN_CLIP_LEN = _number("overhangs", "min_clip_len", kind=int)
BOUNDARY_TOLERANCE = _number("overhangs", "boundary_tolerance", kind=int)

# clustering (scripts/cluster_overhangs_edlib.py, recruit_ambiguous_overhangs.py)
MAX_EDIT_FRAC = _number("clustering", "max_edit_frac")
MERGE_EDIT_FRAC = _number("clustering", "merge_edit_frac")
MIN_CLUSTER_SIZE = _number("clustering", "min_cluster_size", kind=int)
RECRUIT_EDIT_FRAC = _number("clustering", "recruit_edit_frac")

# Max bp between a left and right overhang's junction positions for them to count
# as the same insertion junction. Passed to BOTH blast_clusters_to_ref (which uses
# it when deciding WHERE to place each overhang) and add_pairing_column (which
# uses it to actually pair them). If the two ever disagree, the BLAST step would
# place overhangs to satisfy a rule the pairing step then rejects.
MAX_PAIR_GAP = _number("max_pair_gap", kind=int)

# nested IS insertions (scripts/nested_is.py)
NESTED_MIN_LEN = _number("nested_is", "min_len", kind=int)
NESTED_MIN_PIDENT = _number("nested_is", "min_pident")
NESTED_MIN_ALN = _number("nested_is", "min_aln", kind=int)
NESTED_JUNCTION_TOL = _number("nested_is", "junction_tol", kind=int)
NESTED_TERMINAL_TOL = _number("nested_is", "terminal_tol", kind=int)
NESTED_PAIR_TOL = _number("nested_is", "pair_tol", kind=int)

# software
SRA_TOOLS_ENV = _existing(_setting("conda_envs", "sra_tools"), "conda env file 'sra_tools'")
MINIBWA_ENV = _existing(_setting("conda_envs", "minibwa"), "conda env file 'minibwa'")
PYTHON_PYSAM = _existing(_setting("python", "pysam"), "python interpreter 'pysam'")
PYTHON_PANDAS = _existing(_setting("python", "pandas"), "python interpreter 'pandas'")
BLAST_BIN = _setting("blast_bin")
BLASTN = _existing(os.path.join(BLAST_BIN, "blastn"), "blastn") if BLAST_BIN else None
MAKEBLASTDB = _existing(os.path.join(BLAST_BIN, "makeblastdb"), "makeblastdb") if BLAST_BIN else None

if _problems:
	raise WorkflowError("Problems with the pipeline configuration:\n  - " + "\n  - ".join(_problems)
						+ "\n(settings: config/config.yaml, databases: config/databases.tsv)")

if EXISTING_OVERHANGS:
	ACCESSIONS = sorted(f[:-len(".tar.gz")] for f in os.listdir(EXISTING_OVERHANGS)
						if f.endswith(".tar.gz") and "." not in f[:-len(".tar.gz")])
	if not ACCESSIONS:
		raise WorkflowError(f"no {{sample}}.tar.gz overhang archives in {EXISTING_OVERHANGS}")
elif SAMPLES:
	ACCESSIONS = list(SAMPLES)
else:
	with open(ACCESSIONS_PATH) as f:
		ACCESSIONS = [line.strip() for line in f if line.strip()]

FASTQ_DIR = os.path.join(OUTPUT_DIR, "fastq-files")                     # fasterq, bwa_isfinder, remove_fastqs
BAM_DIR = os.path.join(OUTPUT_DIR, "bam-files")                         # bwa_isfinder, clip_and_cluster
REMOVE_DIR = os.path.join(OUTPUT_DIR, "removed")                        # remove_fastqs
READ_STATS_DIR = os.path.join(OUTPUT_DIR, "read_stats")                 # read_stats, blast_clusters_to_ref
OVERHANG_DIR = os.path.join(OUTPUT_DIR, "overhangs")                    # clip_and_cluster
# where the overhang archives are read from: this run's own, or an earlier run's
OVERHANG_ARCHIVE = os.path.join(EXISTING_OVERHANGS or OVERHANG_DIR, "{id}.tar.gz")  # recruit_ambiguous, cluster_existing_overhangs
CLUSTER_DIR = os.path.join(OUTPUT_DIR, "clusters")                      # clip_and_cluster, recruit_ambiguous
# recruit_ambiguous writes an augmented copy of clusters/{id}/{id}.reads.tsv
# here, with the ambiguous-MAPQ overhangs folded into the clusters they match.
# A separate directory because two rules cannot declare the same output file;
# blast_clusters_to_ref reads from here instead of CLUSTER_DIR.
RECRUITED_DIR = os.path.join(OUTPUT_DIR, "clusters_recruited")          # recruit_ambiguous, blast_clusters_to_ref
GENE_DISRUPTION_DIR = os.path.join(OUTPUT_DIR, "gene_disruption")       # blast_clusters_to_ref
PAIRED_DIR = os.path.join(OUTPUT_DIR, "gene_disruption_paired")         # add_pairing_column
NESTED_DIR = os.path.join(OUTPUT_DIR, "nested_is")                      # nested_is
NESTED_SUMMARY_DIR = os.path.join(OUTPUT_DIR, "nested_is_summary")      # summarize_nested_is
# BLAST db of every configured database, for nested_is. Built from the
# index_database copies, so it always matches what the reads were aligned to.
NESTED_BLASTDB = os.path.join(INDEX_DIR, "nested_is", "IS_Tn")          # nested_is_db, nested_is

# Each configured database is copied and indexed by index_database, and every
# sample is mapped to each one SEPARATELY by bwa_isfinder, then merged.
# Separate rather than one combined index because transposon databases carry
# IS entries and transposons containing IS copies: combined, reads from those
# ISs would map equally well to both and come out MAPQ 0. Each alignment is
# tagged with its database's label as its read group (RG:Z:<label>), which
# overhangs.py carries into its manifest's "database" column, and
# REFERENCES_TSV maps every reference name to its database for joining onto
# any downstream table. Reference names must be unique across databases
# (reference_table fails otherwise), since every downstream output is keyed by
# them.

# accession IDs never contain a dot, so {id} can't ambiguously match a
# "{id}.tar.gz" archive output.
wildcard_constraints:
	id = r"[^/.]+",
	db = r"[A-Za-z0-9_-]+"

rule all:
	input:
		# from-overhangs runs skip everything up to and including clipping
		[] if EXISTING_OVERHANGS else [
			expand(os.path.join(BAM_DIR, "{id}.bam"), id=ACCESSIONS), # make sure bam file is generated
			expand(os.path.join(REMOVE_DIR, "{id}.layout"), id=ACCESSIONS), # make sure fastq files are deleted
			expand(os.path.join(READ_STATS_DIR, "{id}.tsv"), id=ACCESSIONS), # read/base counts, captured before the fastqs go
			expand(os.path.join(OVERHANG_DIR, "{id}.tar.gz"), id=ACCESSIONS), # filtering, extracting, and archiving overhangs
		],
		expand(os.path.join(CLUSTER_DIR, "{id}.tar.gz"), id=ACCESSIONS), # overhang clustering + archiving
		expand(os.path.join(RECRUITED_DIR, "{id}", "{id}.cluster_coverage.tsv"), id=ACCESSIONS), # per-cluster counts after ambiguous-read recruitment
		expand(os.path.join(GENE_DISRUPTION_DIR, "{id}.zip"), id=ACCESSIONS), # BLAST clusters against masked ref genomes, classify gene disruptions
		expand(os.path.join(PAIRED_DIR, "{id}.zip"), id=ACCESSIONS), # pair up left/right overhang clusters at the same junction
		expand(os.path.join(NESTED_DIR, "{id}.tsv.gz"), id=ACCESSIONS), # unplaced overhangs BLASTed against the IS/Tn databases
		os.path.join(NESTED_SUMMARY_DIR, "host_target_summary.tsv") # nested-insertion candidates across all samples

# Reads come either from the samples table (your own fastqs, copied in
# uncompressed so every later rule sees the same {id}_1/_2.fastq layout fasterq
# produces -- remove_fastqs later deletes only these copies, never the
# originals) or from SRA (prefetch + fasterq).
if SAMPLES:
	rule import_fastqs:
		group: "align_is"
		input:
			lambda wc: SAMPLES[wc.id]
		output:
			layout=os.path.join(FASTQ_DIR, "{id}.layout")
		log:
			"logs/import_fastqs/{id}.log"
		resources:
			runtime="30m",
			mem="500M"
		shell:
			"""
			set -euo pipefail
			exec > {log} 2>&1
			mkdir -p {FASTQ_DIR}

			# zcat -f passes uncompressed input through unchanged
			set -- {input}
			if [ $# -eq 2 ]; then
				zcat -f "$1" > {FASTQ_DIR}/{wildcards.id}_1.fastq
				zcat -f "$2" > {FASTQ_DIR}/{wildcards.id}_2.fastq
				echo "paired" > {output.layout}
			else
				zcat -f "$1" > {FASTQ_DIR}/{wildcards.id}.fastq
				echo "single" > {output.layout}
			fi
			"""
else:
	rule prefetch:
		group: "align_is"
		output:
			temp(os.path.join(OUTPUT_DIR, "{id}", "{id}.sra"))
		log:
			"logs/prefetch/{id}.log"
		resources:
			runtime="5m",
			mem_mb=100
		conda: SRA_TOOLS_ENV
		shell:
			"""
			prefetch {wildcards.id} --output-directory {OUTPUT_DIR} > {log} 2>&1
			"""

	rule fasterq:
		group: "align_is"
		input:
			os.path.join(OUTPUT_DIR, "{id}", "{id}.sra")
		output:
			layout=os.path.join(FASTQ_DIR, "{id}.layout")
		log:
			"logs/fasterq/{id}.log"
		resources:
			runtime="10m",
			mem="1G"
		conda: SRA_TOOLS_ENV
		shell:
			"""
			set -euo pipefail

			fasterq-dump {input} --outdir {FASTQ_DIR} > {log} 2>&1

			if [ -f {FASTQ_DIR}/{wildcards.id}_1.fastq ]; then
				echo "paired" > {output.layout}
			else
				echo "single" > {output.layout}
			fi
			"""

# Copy each configured database into INDEX_DIR and build its minibwa index
# there, so the index always matches the exact fasta it was built from. Also
# writes one name/database/length row per reference for reference_table.
rule index_database:
	input:
		lambda wc: DATABASES[wc.db]
	output:
		fasta=os.path.join(INDEX_DIR, "{db}.fa"),
		l2b=os.path.join(INDEX_DIR, "{db}.fa.l2b"),
		mbw=os.path.join(INDEX_DIR, "{db}.fa.mbw"),
		names=os.path.join(INDEX_DIR, "{db}.references.tsv")
	log:
		"logs/index_database/{db}.log"
	resources:
		runtime="30m",
		mem="4G"
	conda: MINIBWA_ENV
	shell:
		"""
		set -euo pipefail
		exec > {log} 2>&1

		# catch the malformed-fasta failure modes that would otherwise index
		# silently: a ">" inside a sequence line (a header glued onto the end
		# of the previous record), a header with no sequence, a repeated name
		awk '
			/^>/ {{ if (name != "" && len == 0) {{ print "empty sequence: " name; bad = 1 }}
			       name = substr($1, 2); len = 0
			       if (name in seen) {{ print "duplicate name: " name; bad = 1 }}
			       seen[name] = 1; next }}
			/>/  {{ print "header glued onto a sequence line, line " NR; bad = 1 }}
			     {{ len += length($0) }}
			END  {{ if (name != "" && len == 0) {{ print "empty sequence: " name; bad = 1 }}
			       exit bad }}
		' {input}

		cp {input} {output.fasta}
		minibwa index -t 1 {output.fasta}

		awk -v db={wildcards.db} '
			/^>/ {{ if (name != "") print name "\t" db "\t" len; name = substr($1, 2); len = 0; next }}
			     {{ len += length($0) }}
			END  {{ if (name != "") print name "\t" db "\t" len }}
		' {output.fasta} > {output.names}
		"""

# reference name -> database lookup across every configured database. Fails on
# a name shared between databases: overhang files, clusters and gene-disruption
# tables are all keyed by reference name alone, so a shared name would merge
# two databases' evidence with no way to separate it afterwards.
rule reference_table:
	input:
		expand(os.path.join(INDEX_DIR, "{db}.references.tsv"), db=DATABASES)
	output:
		REFERENCES_TSV
	resources:
		runtime="5m",
		mem="500M"
	shell:
		"""
		set -euo pipefail
		DUPS=$(cut -f1 {input} | sort | uniq -d | head)
		if [ -n "$DUPS" ]; then
			echo "ERROR: reference names shared between databases:" >&2
			echo "$DUPS" >&2
			exit 1
		fi
		{{ printf "reference\tdatabase\tlength\n"; cat {input}; }} > {output}
		"""

rule bwa_isfinder:
	group: "align_is"
	input:
		layout=os.path.join(FASTQ_DIR, "{id}.layout"),
		indexes=expand(os.path.join(INDEX_DIR, "{db}.fa.{ext}"), db=DATABASES, ext=["l2b", "mbw"]),
		# not read here -- depending on it makes the cross-database name check run first
		references=REFERENCES_TSV
	output:
		bam=os.path.join(BAM_DIR, "{id}.bam")
	log:
		"logs/bwa_isfinder/{id}.log"
	resources:
		runtime="10m",
		mem="1G"
	conda: MINIBWA_ENV
	params:
		fastq_dir=FASTQ_DIR,
		index_dir=INDEX_DIR,
		dbs=" ".join(DATABASES)
	shell:
		"""
		set -euo pipefail

		LAYOUT=$(cat {input.layout})

		if [ "$LAYOUT" = "paired" ]; then
			READS="{params.fastq_dir}/{wildcards.id}_1.fastq {params.fastq_dir}/{wildcards.id}_2.fastq"
		elif [ "$LAYOUT" = "single" ]; then
			READS="{params.fastq_dir}/{wildcards.id}.fastq"
		else
			echo "ERROR: unrecognized layout '$LAYOUT'" >> {log}
			exit 1
		fi

		for f in $READS; do
			if [ ! -f "$f" ]; then
				echo "ERROR: expected read file missing: $f" >> {log}
				exit 1
			fi
		done

		# no -q MAPQ filter here on purpose: low-MAPQ reads still need to reach
		# overhangs.py, which routes them to a separate "ambiguous" FASTA
		# (--min_mapq, default 30) instead of dropping them -- filtering them out
		# here would discard MAPQ 0-19 reads before overhangs.py ever sees them,
		# rather than preserving them for later use as intended.
		#
		# one alignment per database (see DATABASES above), each MAPQ'd only
		# against its own references and tagged RG:Z:<label>, then merged into
		# the single bam the rest of the pipeline expects. A read can appear
		# once per database -- e.g. at the end of a composite transposon,
		# against both its IS and the Tn. Split by database with
		# `samtools view -r <label>`.
		WORK={resources.tmpdir}/bwa_isfinder.{wildcards.id}.$$
		mkdir -p "$WORK"
		trap 'rm -rf "$WORK"' EXIT

		: > {log}
		for DB in {params.dbs}; do
			minibwa map -R "@RG\\tID:$DB\\tSM:{wildcards.id}" "{params.index_dir}/$DB.fa" $READS 2>> {log} | \
			samtools view -b -F 0x904 - | \
			samtools sort -T "$WORK/$DB.sorttmp" -o "$WORK/$DB.bam" 2>> {log}
		done
		samtools merge -o {output.bam} "$WORK"/*.bam 2>> {log}
		"""

# Total reads/bases per sample, summed straight off the fastqs. This HAS to run
# before remove_fastqs deletes them and cannot be recovered afterwards: the bam
# holds only reads that mapped to an IS element or transposon (-F 0x904 off
# bwa_isfinder's `minibwa map`), so it carries no record of the library's total size. Feeds
# blast_clusters_to_ref's est_coverage/coverage_fraction columns.
rule read_stats:
	group: "align_is"
	input:
		layout=os.path.join(FASTQ_DIR, "{id}.layout")
	output:
		stats=os.path.join(READ_STATS_DIR, "{id}.tsv")
	log:
		"logs/read_stats/{id}.log"
	resources:
		runtime="5m",
		mem="100M"
	params:
		fastq_dir=FASTQ_DIR
	shell:
		"""
		set -euo pipefail

		LAYOUT=$(cat {input.layout})

		if [ "$LAYOUT" = "paired" ]; then
			READS="{params.fastq_dir}/{wildcards.id}_1.fastq {params.fastq_dir}/{wildcards.id}_2.fastq"
		elif [ "$LAYOUT" = "single" ]; then
			READS="{params.fastq_dir}/{wildcards.id}.fastq"
		else
			echo "ERROR: unrecognized layout '$LAYOUT'" >> {log}
			exit 1
		fi

		for f in $READS; do
			if [ ! -f "$f" ]; then
				echo "ERROR: expected read file missing: $f" >> {log}
				exit 1
			fi
		done

		mkdir -p "$(dirname {output.stats})"

		# FNR (not NR) so the every-4th-line sequence stride restarts per file,
		# which keeps this correct for the paired case regardless of the two
		# mates' line counts.
		awk 'FNR % 4 == 2 {{ n++; b += length($0) }}
		     END {{ printf "sample\\tn_reads\\ttotal_bases\\tmean_read_length\\n{wildcards.id}\\t%d\\t%d\\t%.2f\\n", n, b, (n ? b / n : 0) }}' \
		    $READS > {output.stats}.tmp 2>> {log}
		mv {output.stats}.tmp {output.stats}
		"""

rule remove_fastqs:
	group: "align_is"
	input:
		os.path.join(BAM_DIR, "{id}.bam"),
		# ordering edge: the stats must be summed while the fastqs still exist
		stats=os.path.join(READ_STATS_DIR, "{id}.tsv"),
		layout=os.path.join(FASTQ_DIR, "{id}.layout")
	output:
		removed=os.path.join(REMOVE_DIR, "{id}.layout")
	log:
		"logs/remove_fastqs/{id}.log"
	resources:
		runtime="1m",
		mem="100M"
	conda: MINIBWA_ENV
	params:
		fastq_dir=FASTQ_DIR
	shell:
		"""
		set -euo pipefail

		LAYOUT=$(cat {input.layout})

		if [ "$LAYOUT" = "paired" ]; then
			READS="{params.fastq_dir}/{wildcards.id}_1.fastq {params.fastq_dir}/{wildcards.id}_2.fastq"
		elif [ "$LAYOUT" = "single" ]; then
			READS="{params.fastq_dir}/{wildcards.id}.fastq"
		else
			echo "ERROR: unrecognized layout '$LAYOUT'" >> {log}
			exit 1
		fi

		for f in $READS; do
			if [ ! -f "$f" ]; then
				echo "ERROR: expected read file missing: $f" >> {log}
				exit 1
			fi
		done

		rm -f $READS
		echo "removed {wildcards.id}" > {output.removed}
		"""


# clip_overhangs + clustering are merged into one rule so per-IS-element
# intermediates never touch /n/scratch's inode quota: everything lives in the
# job's local $TMPDIR, and only the small manifests/reads.tsv plus one
# tar.gz per stage get persisted to scratch.
if not EXISTING_OVERHANGS:
	rule clip_and_cluster:
		group: "overhangs"
		input:
			bam=os.path.join(BAM_DIR, "{id}.bam"),
			overhang_script=os.path.join(SCRIPT_DIR, "scripts", "overhangs.py"),
			cluster_script=os.path.join(SCRIPT_DIR, "scripts", "cluster_overhangs_edlib.py")
		output:
			overhang_manifest=os.path.join(OVERHANG_DIR, "{id}", "{id}.manifest.tsv"),
			overhang_archive=os.path.join(OVERHANG_DIR, "{id}.tar.gz"),
			cluster_manifest=os.path.join(CLUSTER_DIR, "{id}", "{id}.cluster_manifest.tsv"),
			reads_tsv=os.path.join(CLUSTER_DIR, "{id}", "{id}.reads.tsv"),
			cluster_archive=os.path.join(CLUSTER_DIR, "{id}.tar.gz")
		log:
			"logs/clip_and_cluster/{id}.log"
		resources:
			runtime="15m",
			mem="500M"
		shell:
			"""
			set -euo pipefail
			exec > {log} 2>&1

			WORK={resources.tmpdir}/clip_and_cluster.{wildcards.id}.$$
			OVERHANG_WORK="$WORK/overhangs"
			CLUSTER_WORK="$WORK/clusters"
			mkdir -p "$OVERHANG_WORK" "$CLUSTER_WORK"
			trap 'rm -rf "$WORK"' EXIT

			# stage 1: clip overhangs from the bam (needs pysam/samtools -- minibwa env)
			{PYTHON_PYSAM} {input.overhang_script} \
				{input.bam} "$OVERHANG_WORK" --sample_id {wildcards.id} --skip_is_hit_filter \
				--min_is_coverage_depth {MIN_ELEMENT_COVERAGE_DEPTH} \
				--min_is_coverage_fraction {MIN_ELEMENT_COVERAGE_FRACTION} \
				--min_mapq {MIN_MAPQ} --min_pident {MIN_PIDENT} \
				--min_clip_len {MIN_CLIP_LEN} --boundary_tolerance {BOUNDARY_TOLERANCE}

			# stage 2: cluster each is_element/side pair with a position-anchored edit
			# distance (edlib SHW/prefix mode) instead of CD-HIT's percent-identity +
			# coverage-threshold approach -- see cluster_overhangs_edlib.py's module
			# docstring for the anchoring/algorithm rationale. One script call replaces
			# both the old cd-hit-est loop and the old combine_cdhit_clusters.py step.
			#
			# TODO(infra): `edlib` (pip install edlib, pure C++ ext, no Rust toolchain)
			# is not yet installed in the bakta env -- add it there before this rule
			# will run (bakta already has pandas, needed here too, so reusing it avoids
			# a new env). Set python: pandas in config/config.yaml if a different env is used.
			{PYTHON_PANDAS} {input.cluster_script} \
				"$OVERHANG_WORK/{wildcards.id}.manifest.tsv" \
				--sample_id {wildcards.id} \
				--output_dir "$CLUSTER_WORK" \
				--max_edit_frac {MAX_EDIT_FRAC} --merge_edit_frac {MERGE_EDIT_FRAC} \
				--min_cluster_size {MIN_CLUSTER_SIZE}

			# persist only the small, useful results. cluster_overhangs_edlib.py's
			# manifest (is_element/side/n_in/n_clusters) carries no tmpdir-relative
			# paths, so the sed below is a no-op pass-through today -- kept so this
			# still self-documents/rewrites correctly if a future manifest column ever
			# does reference $CLUSTER_WORK-relative paths again.
			mkdir -p "{OVERHANG_DIR}/{wildcards.id}" "{CLUSTER_DIR}/{wildcards.id}"
			cp "$OVERHANG_WORK/{wildcards.id}.manifest.tsv" {output.overhang_manifest}
			sed "s|$CLUSTER_WORK/|{CLUSTER_DIR}/{wildcards.id}/|g" \
				"$CLUSTER_WORK/{wildcards.id}.cluster_manifest.tsv" > {output.cluster_manifest}
			cp "$CLUSTER_WORK/{wildcards.id}.reads.tsv" {output.reads_tsv}

			# archive the raw intermediates -- built once, verified, then moved into place
			tar -C "$OVERHANG_WORK" -czf {output.overhang_archive}.tmp .
			tar -tzf {output.overhang_archive}.tmp > /dev/null
			mv {output.overhang_archive}.tmp {output.overhang_archive}

			tar -C "$CLUSTER_WORK" -czf {output.cluster_archive}.tmp .
			tar -tzf {output.cluster_archive}.tmp > /dev/null
			mv {output.cluster_archive}.tmp {output.cluster_archive}
			"""
else:
	# Clustering half of clip_and_cluster, run on an existing overhang archive
	# (config: existing_overhangs). Same outputs, so everything downstream is
	# unchanged. The archive is only read: it is unpacked into the job's
	# $TMPDIR and left where it is.
	rule cluster_existing_overhangs:
		group: "overhangs"
		input:
			overhang_archive=OVERHANG_ARCHIVE,
			cluster_script=os.path.join(SCRIPT_DIR, "scripts", "cluster_overhangs_edlib.py")
		output:
			cluster_manifest=os.path.join(CLUSTER_DIR, "{id}", "{id}.cluster_manifest.tsv"),
			reads_tsv=os.path.join(CLUSTER_DIR, "{id}", "{id}.reads.tsv"),
			cluster_archive=os.path.join(CLUSTER_DIR, "{id}.tar.gz")
		log:
			"logs/cluster_existing_overhangs/{id}.log"
		resources:
			runtime="15m",
			mem="500M"
		shell:
			"""
			set -euo pipefail
			exec > {log} 2>&1

			WORK={resources.tmpdir}/cluster_existing_overhangs.{wildcards.id}.$$
			OVERHANG_WORK="$WORK/overhangs"
			CLUSTER_WORK="$WORK/clusters"
			mkdir -p "$OVERHANG_WORK" "$CLUSTER_WORK"
			trap 'rm -rf "$WORK"' EXIT

			tar -xzf {input.overhang_archive} -C "$OVERHANG_WORK"

			# the manifest's fasta_path column points at the $TMPDIR of the job
			# that wrote the archive, long gone -- repoint it at the unpacked copies
			awk -F '\\t' -v OFS='\\t' -v dir="$OVERHANG_WORK" '
				NR == 1 {{ for (i = 1; i <= NF; i++) if ($i == "fasta_path") col = i
				          if (!col) {{ print "no fasta_path column in manifest" > "/dev/stderr"; exit 1 }}
				          print; next }}
				$col != "" {{ n = split($col, parts, "/"); $col = dir "/" parts[n] }}
				{{ print }}
			' "$OVERHANG_WORK/{wildcards.id}.manifest.tsv" > "$WORK/manifest.tsv"

			{PYTHON_PANDAS} {input.cluster_script} \
				"$WORK/manifest.tsv" \
				--sample_id {wildcards.id} \
				--output_dir "$CLUSTER_WORK" \
				--max_edit_frac {MAX_EDIT_FRAC} --merge_edit_frac {MERGE_EDIT_FRAC} \
				--min_cluster_size {MIN_CLUSTER_SIZE}

			mkdir -p "{CLUSTER_DIR}/{wildcards.id}"
			cp "$CLUSTER_WORK/{wildcards.id}.cluster_manifest.tsv" {output.cluster_manifest}
			cp "$CLUSTER_WORK/{wildcards.id}.reads.tsv" {output.reads_tsv}

			tar -C "$CLUSTER_WORK" -czf {output.cluster_archive}.tmp .
			tar -tzf {output.cluster_archive}.tmp > /dev/null
			mv {output.cluster_archive}.tmp {output.cluster_archive}
			"""


# Recruit the ambiguous-MAPQ overhangs into the clusters built above.
# clip_and_cluster only clusters the confident left/right buckets; overhangs.py
# stage 3 holds the low-MAPQ (ZA=1) overhangs back in a per-IS "__ambiguous"
# FASTA because their left/right call is the unreliable part. Now that
# representatives exist, each one is re-aligned against the representatives of
# its IS element (both sides, with reverse-complement handling) and joins the
# best match -- see scripts/recruit_ambiguous_overhangs.py's module docstring.
#
# The output reads.tsv keeps that filename on purpose: blast_clusters_to_ref
# globs "*.reads.tsv" and already derives n_seqs by counting rows per cluster,
# so pointing it at RECRUITED_DIR is all that is needed for the corrected
# cluster coverage to reach the final per-IS tables (and, through n_seqs, the
# coverage_fraction column).
rule recruit_ambiguous:
	group: "overhangs"
	input:
		reads_tsv=os.path.join(CLUSTER_DIR, "{id}", "{id}.reads.tsv"),
		overhang_archive=OVERHANG_ARCHIVE,
		script=os.path.join(SCRIPT_DIR, "scripts", "recruit_ambiguous_overhangs.py")
	output:
		reads_tsv=os.path.join(RECRUITED_DIR, "{id}", "{id}.reads.tsv"),
		coverage=os.path.join(RECRUITED_DIR, "{id}", "{id}.cluster_coverage.tsv")
	log:
		"logs/recruit_ambiguous/{id}.log"
	resources:
		runtime="5m",
		mem="500M"
	shell:
		"""
		set -euo pipefail
		exec > {log} 2>&1

		WORK={resources.tmpdir}/recruit_ambiguous.{wildcards.id}.$$
		mkdir -p "$WORK"
		trap 'rm -rf "$WORK"' EXIT

		# the ambiguous FASTAs only exist inside the overhang archive. NOTE the
		# script globs this directory rather than reading fasta_path out of
		# {wildcards.id}.manifest.tsv: clip_and_cluster copies that manifest
		# verbatim from its own $TMPDIR (only the CLUSTER manifest gets
		# sed-rewritten), so its paths point at a directory that is long gone.
		tar -xzf {input.overhang_archive} -C "$WORK"

		# an archive with no __ambiguous files is fine and expected for any data
		# produced before overhangs.py grew the ZA tag -- the script then just
		# passes the clustered reads through with recruited=0.
		{PYTHON_PANDAS} {input.script} \
			--reads_tsv {input.reads_tsv} \
			--overhang_dir "$WORK" \
			--sample_id {wildcards.id} \
			--output_dir {RECRUITED_DIR}/{wildcards.id} \
			--recruit_edit_frac {RECRUIT_EDIT_FRAC}
		"""


# BLASTs each accession's cluster representative sequences against the masked
# ref-genome BLAST db (BLASTDB/GBFF_SUFFIX above) and classifies whether each
# cluster's IS-insertion junction disrupted a gene (see
# scripts/blast_clusters_to_ref.py's module docstring for the full
# classification logic).
rule blast_clusters_to_ref:
	group: "gene_disruption"
	input:
		# reads from RECRUITED_DIR, not CLUSTER_DIR, so the n_seqs this script
		# computes includes the recruited ambiguous-MAPQ reads
		reads_tsv=os.path.join(RECRUITED_DIR, "{id}", "{id}.reads.tsv"),
		script=os.path.join(SCRIPT_DIR, "scripts", "blast_clusters_to_ref.py"),
		blastdb_file=BLASTDB + ".nsq",
		# no fastqs to count in a from-overhangs run (existing_overhangs), so
		# no read_stats -- est_coverage/coverage_fraction are then left empty
		read_stats=[] if EXISTING_OVERHANGS else [os.path.join(READ_STATS_DIR, "{id}.tsv")],
		genome_sizes=[GENOME_SIZES] if GENOME_SIZES and not EXISTING_OVERHANGS else [],
		ref_gbffs=REF_GBFFS
	params:
		# both optional; without them est_coverage/coverage_fraction stay empty
		read_stats_arg=lambda wc, input: f"--read_stats {input.read_stats}" if input.read_stats else "",
		genome_sizes_arg=lambda wc, input: f"--genome_sizes {input.genome_sizes}" if input.genome_sizes else ""
	output:
		os.path.join(GENE_DISRUPTION_DIR, "{id}.zip")
	log:
		"logs/blast_clusters_to_ref/{id}.log"
	resources:
		runtime="10m",
		mem="1G"
	# group-components bundles 10 of these into one SLURM submission (see
	# slurmprofile/config.yaml), summing threads across the group, so this
	# stays at 1 (10 components x 1 thread = 10 cores/submission, under the
	# short partition's 20-CPU/job cap -- blastn on this small db is fast
	# enough single-threaded that this isn't a meaningful runtime cost).
	threads: 1
	shell:
		"""
		set -euo pipefail
		exec > {log} 2>&1

		TMP_DIR={resources.tmpdir}/blast_clusters_to_ref.{wildcards.id}.$$
		mkdir -p "$TMP_DIR"
		trap 'rm -rf "$TMP_DIR"' EXIT

		{PYTHON_PANDAS} {input.script} \
			--cluster_dirs {RECRUITED_DIR}/{wildcards.id} \
			--blastdb {BLASTDB} \
			--ref_genomes_dir {REF_GENOMES_DIR} \
			--accessions {REF_ACCESSIONS_STR} \
			--gbff_suffix {GBFF_SUFFIX} \
			--output_dir {GENE_DISRUPTION_DIR} \
			--tmp_dir "$TMP_DIR" \
			--threads {threads} \
			{params.read_stats_arg} \
			{params.genome_sizes_arg} \
			--max_pair_gap {MAX_PAIR_GAP}

		# a sample with no clusters at all (every element/side under
		# min_cluster_size) gets no zip from the script -- write an empty one so
		# it reads as "no insertions" instead of failing the job (and, grouped,
		# every other sample submitted with it). add_pairing_column passes an
		# empty zip through as an empty zip.
		if [ ! -f {output} ]; then
			{PYTHON_PANDAS} -c "import zipfile, sys; zipfile.ZipFile(sys.argv[1], 'w').close()" {output}
			echo "no clusters -- wrote empty {output}"
		fi
		"""


# Adds a "pairing" column to every IS-element TSV inside a blast_clusters_to_ref
# zip, matching up left/right overhang clusters that land at the same
# insertion junction (see pairing/add_pairing_column.py's module docstring).
rule add_pairing_column:
	input:
		zip=os.path.join(GENE_DISRUPTION_DIR, "{id}.zip"),
		script=os.path.join(SCRIPT_DIR, "pairing", "add_pairing_column.py")
	output:
		os.path.join(PAIRED_DIR, "{id}.zip")
	log:
		"logs/add_pairing_column/{id}.log"
	resources:
		runtime="2m",
		mem="200M"
	shell:
		"""
		set -euo pipefail
		exec > {log} 2>&1

		{PYTHON_PYSAM} {input.script} --zip {input.zip} --out-dir {PAIRED_DIR} --max-gap {MAX_PAIR_GAP}
		"""


# BLAST db for nested_is: every configured database (the index_database
# copies, so it matches what the reads were aligned to) in one db. Depends on
# REFERENCES_TSV so reference_table's shared-name check has passed first --
# nested_is reports targets by name alone.
rule nested_is_db:
	input:
		fastas=expand(os.path.join(INDEX_DIR, "{db}.fa"), db=DATABASES),
		references=REFERENCES_TSV
	output:
		NESTED_BLASTDB + ".nsq"
	log:
		"logs/nested_is_db/nested_is_db.log"
	resources:
		runtime="15m",
		mem="2G"
	shell:
		"""
		set -euo pipefail
		exec > {log} 2>&1

		cat {input.fastas} > {NESTED_BLASTDB}.fa
		{MAKEBLASTDB} -in {NESTED_BLASTDB}.fa -dbtype nucl -out {NESTED_BLASTDB} -parse_seqids
		"""


# Nested IS insertions: the overhang clusters blast_clusters_to_ref left as
# no_hit (their flank matched none of the IS-masked reference genomes) are
# BLASTed against the IS/Tn databases, to find host elements whose flank lands
# inside another element. One per-cluster table per sample; see
# scripts/nested_is.py's header for the classes. A few seconds per sample, so
# these are bundled into one SLURM submission (slurmprofile/config.yaml,
# group-components: nested_is).
rule nested_is:
	group: "nested_is"
	input:
		zip=os.path.join(PAIRED_DIR, "{id}.zip"),
		blastdb=NESTED_BLASTDB + ".nsq",
		references=REFERENCES_TSV,
		script=os.path.join(SCRIPT_DIR, "scripts", "nested_is.py")
	output:
		os.path.join(NESTED_DIR, "{id}.tsv.gz")
	log:
		"logs/nested_is/{id}.log"
	resources:
		runtime="5m",
		mem="1G"
	threads: 1
	shell:
		"""
		set -euo pipefail
		exec > {log} 2>&1

		{PYTHON_PANDAS} {input.script} \
			--zips {input.zip} \
			--blastdb {NESTED_BLASTDB} \
			--references {input.references} \
			--output {output} \
			--blastn {BLASTN} \
			--threads {threads} \
			--tmp_dir {resources.tmpdir} \
			--min_len {NESTED_MIN_LEN} --min_pident {NESTED_MIN_PIDENT} \
			--min_aln {NESTED_MIN_ALN} --junction_tol {NESTED_JUNCTION_TOL} \
			--terminal_tol {NESTED_TERMINAL_TOL} --pair_tol {NESTED_PAIR_TOL}
		"""


# Combines every sample's nested_is table into class_summary.tsv (clusters per
# class) and host_target_summary.tsv (one row per host -> target pair of
# nested-insertion candidates; see scripts/summarize_nested_is.py). The table
# paths go through a list file because ~150k of them overflow a command line.
rule summarize_nested_is:
	input:
		tables=expand(os.path.join(NESTED_DIR, "{id}.tsv.gz"), id=ACCESSIONS),
		script=os.path.join(SCRIPT_DIR, "scripts", "summarize_nested_is.py")
	output:
		classes=os.path.join(NESTED_SUMMARY_DIR, "class_summary.tsv"),
		host_target=os.path.join(NESTED_SUMMARY_DIR, "host_target_summary.tsv")
	log:
		"logs/summarize_nested_is/summarize_nested_is.log"
	resources:
		runtime="2h",
		mem="16G"
	run:
		os.makedirs(NESTED_SUMMARY_DIR, exist_ok=True)
		table_list = os.path.join(NESTED_SUMMARY_DIR, "tables.txt")
		with open(table_list, "w") as fh:
			fh.write("\n".join(input.tables) + "\n")
		shell(
			"{PYTHON_PANDAS} {input.script} --table_list {table_list} "
			"--output_dir {NESTED_SUMMARY_DIR} --site_tol {NESTED_PAIR_TOL} > {log} 2>&1"
		)
