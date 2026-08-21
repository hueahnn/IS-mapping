# IS-element insertion pipeline: given a list of SRA accession IDs, downloads
# each genome's reads, aligns them to the ISfinder IS-element database,
# extracts and clusters the sequence "overhangs" flanking every IS hit, BLASTs
# those clusters against masked E. coli reference genomes to classify whether
# the insertion disrupted a gene, and pairs up left/right overhang clusters
# that land at the same insertion junction.
# begin: 06/18/2026
# run: snakemake --profile slurmprofile --rerun-incomplete --use-conda --executor slurm

import os

# accession list to run the pipeline over
ACCESSIONS_PATH = config.get("input_path", "/home/hua575/baymlab/mapping/atb/ecoli_atb_sra_accessions.txt")
with open(ACCESSIONS_PATH) as f:
	ACCESSIONS = [line.strip() for line in f if line.strip()]

# scratch working directory holding all per-accession intermediate/output files
PATH = "/n/scratch/users/h/hua575/atb_filtered"
FASTQ_DIR = os.path.join(PATH, "fastq-files")                     # fasterq, bwa_isfinder, remove_fastqs
BAM_DIR = os.path.join(PATH, "bam-files")                         # bwa_isfinder, clip_and_cluster
REMOVE_DIR = os.path.join(PATH, "removed")                        # remove_fastqs
OVERHANG_DIR = os.path.join(PATH, "overhangs")                    # clip_and_cluster
CLUSTER_DIR = os.path.join(PATH, "clusters")                      # clip_and_cluster, blast_clusters_to_ref
GENE_DISRUPTION_DIR = os.path.join(PATH, "gene_disruption")       # blast_clusters_to_ref
PAIRED_DIR = os.path.join(PATH, "gene_disruption_paired")         # add_pairing_column

SCRIPT_DIR = "/n/data1/hms/dbmi/baym/hue/mapping"                   # clip_and_cluster, blast_clusters_to_ref, add_pairing_column
IS_DB = "ISfinder_database-master_2026/IS.database.collapsed99.fa" # bwa_isfinder

# reference genomes / BLAST db for blast_clusters_to_ref (masking.py's
# ISfinder-BLAST-based excision -- see scripts/masking.py)
REF_GENOMES_DIR = "/home/hua575/baymlab/mapping/ref_genomes/ecoli"
REF_ACCESSIONS = [  # must match the accessions BLASTDB below was built from
	"GCA_000005845.2", "GCA_900096825.1", "GCA_000692435.1",
	"GCA_002473875.1", "GCA_000163235.1", "GCA_002966755.1",
]
REF_ACCESSIONS_STR = " ".join(REF_ACCESSIONS)  # shell {}-formatting can't eval " ".join(...) inline
BLASTDB = "/home/hua575/baymlab/mapping/ref_genomes/blastdb/ecoli_masked_combined_v2"
GBFF_SUFFIX = "no_IS_v2"
REF_GBFFS = expand(
	os.path.join(REF_GENOMES_DIR, "{acc}", "{acc}_" + GBFF_SUFFIX + ".gbff"), acc=REF_ACCESSIONS
)

# accession IDs never contain a dot, so {id} can't ambiguously match a
# "{id}.tar.gz" archive output.
wildcard_constraints:
	id = r"[^/.]+"

rule all:
	input:
		expand(os.path.join(BAM_DIR, "{id}.bam"), id=ACCESSIONS), # make sure bam file is generated
		expand(os.path.join(REMOVE_DIR, "{id}.layout"), id=ACCESSIONS), # make sure fastq files are deleted
		expand(os.path.join(OVERHANG_DIR, "{id}.tar.gz"), id=ACCESSIONS), # filtering, extracting, and archiving overhangs
		expand(os.path.join(CLUSTER_DIR, "{id}.tar.gz"), id=ACCESSIONS), # overhang clustering + archiving
		expand(os.path.join(GENE_DISRUPTION_DIR, "{id}.zip"), id=ACCESSIONS), # BLAST clusters against masked ref genomes, classify gene disruptions
		expand(os.path.join(PAIRED_DIR, "{id}.zip"), id=ACCESSIONS) # pair up left/right overhang clusters at the same junction

rule prefetch:
	group: "align_is"
	output:
		temp(os.path.join(PATH, "{id}", "{id}.sra"))
	log:
		"logs/prefetch/{id}.log"
	resources:
		runtime="5m",
		mem_mb=100
	conda: "/home/hua575/baymlab/sibmi/conda-envs/ncbi-datasets.yaml"
	shell:
		"""
		prefetch {wildcards.id} --output-directory {PATH} > {log} 2>&1
		"""

rule fasterq:
	group: "align_is"
	input:
		os.path.join(PATH, "{id}", "{id}.sra")
	output:
		layout=os.path.join(FASTQ_DIR, "{id}.layout")
	log:
		"logs/fasterq/{id}.log"
	resources:
		runtime="10m",
		mem="1G"
	conda: "/home/hua575/baymlab/sibmi/conda-envs/ncbi-datasets.yaml"
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

# for separating single vs paired ends
def get_fastqs(wildcards):
	layout_file = checkpoints.fasterq.get(id=wildcards.id).output.layout
	with open(layout_file) as f:
		layout = f.read().strip()
	if layout == "paired":
		return [
			os.path.join(FASTQ_DIR, f"{wildcards.id}_1.fastq"),
			os.path.join(FASTQ_DIR, f"{wildcards.id}_2.fastq"),
	]
	elif layout == "single":
		return [os.path.join(FASTQ_DIR, f"{wildcards.id}.fastq")]
	else:
		raise ValueError(f"Unrecognized layout '{layout}' for {wildcards.id}")

rule bwa_isfinder:
	group: "align_is"
	input:
		# fastqs=get_fastqs,
		layout=os.path.join(FASTQ_DIR, "{id}.layout")
	output:
		bam=os.path.join(BAM_DIR, "{id}.bam")
	log:
		"logs/bwa_isfinder/{id}.log"
	resources:
		runtime="5m",
		mem="500M"
	conda: "/home/hua575/baymlab/sibmi/conda-envs/minibwa.yaml"
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

		minibwa map isfinder $READS 2> {log} | \
		samtools view -b -F 0x904 -q 20 - | \
		samtools sort -o {output.bam} 2>> {log}
		"""

rule remove_fastqs:
	group: "align_is"
	input:
		os.path.join(BAM_DIR, "{id}.bam"),
		layout=os.path.join(FASTQ_DIR, "{id}.layout")
	output:
		removed=os.path.join(REMOVE_DIR, "{id}.layout")
	log:
		"logs/remove_fastqs/{id}.log"
	resources:
		runtime="1m",
		mem="100M"
	conda: "/home/hua575/baymlab/sibmi/conda-envs/minibwa.yaml"
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
		/home/hua575/miniconda3/envs/minibwa/bin/python {input.overhang_script} \
			{input.bam} "$OVERHANG_WORK" --sample_id {wildcards.id} --skip_is_hit_filter

		# stage 2: cluster each is_element/side pair with a position-anchored edit
		# distance (edlib SHW/prefix mode) instead of CD-HIT's percent-identity +
		# coverage-threshold approach -- see cluster_overhangs_edlib.py's module
		# docstring for the anchoring/algorithm rationale. One script call replaces
		# both the old cd-hit-est loop and the old combine_cdhit_clusters.py step.
		#
		# TODO(infra): `edlib` (pip install edlib, pure C++ ext, no Rust toolchain)
		# is not yet installed in the bakta env -- add it there before this rule
		# will run (bakta already has pandas, needed here too, so reusing it avoids
		# a new env). Update the interpreter path below if a different env is used.
		/home/hua575/miniconda3/envs/bakta/bin/python {input.cluster_script} \
			"$OVERHANG_WORK/{wildcards.id}.manifest.tsv" \
			--sample_id {wildcards.id} \
			--output_dir "$CLUSTER_WORK" \
			--max_edit_frac 0.10 --merge_edit_frac 0.15 --min_cluster_size 2

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


# BLASTs each accession's cluster representative sequences against the masked
# ref-genome BLAST db (BLASTDB/GBFF_SUFFIX above) and classifies whether each
# cluster's IS-insertion junction disrupted a gene (see
# scripts/blast_clusters_to_ref.py's module docstring for the full
# classification logic).
rule blast_clusters_to_ref:
	group: "gene_disruption"
	input:
		reads_tsv=os.path.join(CLUSTER_DIR, "{id}", "{id}.reads.tsv"),
		script=os.path.join(SCRIPT_DIR, "scripts", "blast_clusters_to_ref.py"),
		blastdb_file=BLASTDB + ".nsq",
		ref_gbffs=REF_GBFFS
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

		/home/hua575/miniconda3/envs/bakta/bin/python {input.script} \
			--cluster_dirs {CLUSTER_DIR}/{wildcards.id} \
			--blastdb {BLASTDB} \
			--ref_genomes_dir {REF_GENOMES_DIR} \
			--accessions {REF_ACCESSIONS_STR} \
			--gbff_suffix {GBFF_SUFFIX} \
			--output_dir {GENE_DISRUPTION_DIR} \
			--tmp_dir "$TMP_DIR" \
			--threads {threads}
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

		/home/hua575/miniconda3/envs/minibwa/bin/python {input.script} --zip {input.zip} --out-dir {PAIRED_DIR}
		"""
