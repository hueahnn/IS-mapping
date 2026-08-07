# begin: 06/18/2026
# purpose: downloading fastq files from sra accession IDs. fastq file --> bwa align to isfinder db. for every IS hit, 
# extract overhangs and cluster. align clusters to ref genome
# command to run: snakemake --profile slurmprofile --rerun-incomplete --use-conda --executor slurm

import os

### global ####
# old: pre atb db
# ACCESSIONS_PATH = "/n/scratch/users/h/hua575/ecoli_sra_accessions_v3.txt"
# PATH = "/n/scratch/users/h/hua575/v2"
# FASTQ_DIR = os.path.join(PATH, "fastq-files-v2")
# batches
# BATCH_FILE = "/n/scratch/users/h/hua575/v2/fastq-file-names-1-random5.txt"
# BATCH_FILE = "/n/scratch/users/h/hua575/v2/fastq-file-names-1-random10.txt"
# BATCH_FILE = "/n/scratch/users/h/hua575/v2/fastq-file-names-1-random20.txt"
# BATCH_FILE = "/n/scratch/users/h/hua575/v2/fastq-file-names-1-random100.txt"
# BATCH_FILE = "/n/scratch/users/h/hua575/v2/fastq-file-names-1.txt"

# new: using atb db filtered (149k IDs)
ACCESSIONS_PATH = config.get("input_path", "/home/hua575/baymlab/mapping/atb/ecoli_atb_sra_accessions.txt")
PATH = "/n/scratch/users/h/hua575/atb_filtered"
# PATH = "/n/scratch/users/h/hua575/test0708"
FASTQ_DIR = os.path.join(PATH, "fastq-files")
# batches
# BATCH_FILE = "/home/hua575/baymlab/mapping/atb/ecoli_atb_sra_accessions_random10.txt"
# BATCH_FILE = "/home/hua575/baymlab/mapping/atb/ecoli_atb_sra_accessions_random10_2.txt"
BATCH_FILE = "/home/hua575/baymlab/mapping/atb/ecoli_atb_sra_accessions_random1000.txt"
# BATCH_FILE = "/home/hua575/baymlab/mapping/atb/ecoli_atb_sra_accessions_2000.txt"
# BATCH_FILE = "/home/hua575/baymlab/mapping/atb/ecoli_atb_sra_accessions_2000_2.txt"
# BATCH_FILE = "/home/hua575/baymlab/mapping/atb/ecoli_atb_sra_accessions_tail_5000.txt"

### local rules ###
IS_DB = "ISfinder_database-master_2026/IS.database.collapsed99.fa" # bwa_isfinder
BAM_DIR = os.path.join(PATH, "bam-files") # bwa_isfinder
# BAM_DIR = os.path.join("/home/hua575/baymlab/mapping/", "bam-files") # bwa_isfinder
OVERHANG_DIR = os.path.join(PATH, "overhangs") # clip_and_cluster
SCRIPT_DIR = "/n/data1/hms/dbmi/baym/hue/mapping" # clip_and_cluster
REMOVE_DIR = os.path.join(PATH, "removed") # remove_fastqs
CDHIT_DIR = os.path.join(PATH, "cd-hit") # clip_and_cluster

GENE_DISRUPTION_DIR = os.path.join(PATH, "gene_disruption") # blast_clusters_to_ref
REF_GENOMES_DIR = "/home/hua575/baymlab/mapping/ref_genomes/ecoli" # blast_clusters_to_ref
BLASTDB = "/home/hua575/baymlab/mapping/ref_genomes/blastdb/ecoli_masked_combined" # blast_clusters_to_ref
REF_ACCESSIONS = [ # blast_clusters_to_ref -- must match the accessions the BLASTDB was built from
	"GCA_000005845.2", "GCA_900096825.1", "GCA_000692435.1",
	"GCA_002473875.1", "GCA_000163235.1", "GCA_002966755.1",
]
REF_ACCESSIONS_STR = " ".join(REF_ACCESSIONS) # blast_clusters_to_ref -- shell {}-formatting can't eval " ".join(...) inline

# v2 masked ref genomes: masking.py's ISfinder-BLAST-based excision instead of
# the original Bakta-/product-annotation regex -- see blast_clusters_to_ref_v2
GENE_DISRUPTION_DIR_V2 = os.path.join(PATH, "gene_disruption_v2") # blast_clusters_to_ref_v2
BLASTDB_V2 = "/home/hua575/baymlab/mapping/ref_genomes/blastdb/ecoli_masked_combined_v2" # blast_clusters_to_ref_v2

# change input to whatever file of accessions (if breaking into chunks, can split up original v3 accession file)
with open(ACCESSIONS_PATH) as f:
	ACCESSIONS = [line.strip() for line in f if line.strip()]

# accession IDs never contain a dot -- constrains {id} so it can't also match
# "{id}.tar.gz", which would otherwise be ambiguous with the {id}.tar.gz archive outputs
wildcard_constraints:
	id = r"[^/.]+"

rule all:
	input:
		expand(os.path.join(BAM_DIR, "{id}.bam"), id=ACCESSIONS), # make sure bam file is generated
		expand(os.path.join(REMOVE_DIR, "{id}.layout"), id=ACCESSIONS), # make sure fastq files are deleted
		expand(os.path.join(OVERHANG_DIR, "{id}.tar.gz"), id=ACCESSIONS), # filtering, extracting, and archiving overhangs
		expand(os.path.join(CDHIT_DIR, "{id}.tar.gz"), id=ACCESSIONS), # cd-hit clustering + archiving
		expand(os.path.join(GENE_DISRUPTION_DIR, "{id}.zip"), id=ACCESSIONS) # BLAST clusters against masked ref genomes, classify gene disruptions

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


# clip_overhangs and cd_hit used to be separate rules, each writing dozens of
# small per-IS-element files straight to /n/scratch, with a third pair of rules
# archiving them afterward. Since clip_overhangs runs much faster than cd_hit
# can consume its output, that left a large, unbounded "clipped but not yet
# clustered" backlog of small raw files sitting on the quota-limited
# filesystem at all times -- which is what blew the inode quota on 2026-07-22.
#
# Merged into one rule: all per-IS-element intermediates (raw overhang fastas,
# raw cd-hit cluster files) now live only in the job's local $TMPDIR, which
# never touches /n/scratch's inode quota and disappears when the job ends.
# Only the small manifest/reads.tsv results and one combined tar.gz per stage
# get written to persistent storage, per accession, atomically -- so there's
# no window where a large number of small files can accumulate on scratch.
rule clip_and_cluster:
	group: "overhangs"
	input:
		bam=os.path.join(BAM_DIR, "{id}.bam"),
		overhang_script=os.path.join(SCRIPT_DIR, "scripts", "overhangs.py"),
		combine_script=os.path.join(SCRIPT_DIR, "scripts", "combine_cdhit_clusters.py")
	output:
		overhang_manifest=os.path.join(OVERHANG_DIR, "{id}", "{id}.manifest.tsv"),
		overhang_archive=os.path.join(OVERHANG_DIR, "{id}.tar.gz"),
		cdhit_manifest=os.path.join(CDHIT_DIR, "{id}", "{id}.cdhit_manifest.tsv"),
		reads_tsv=os.path.join(CDHIT_DIR, "{id}", "{id}.reads.tsv"),
		cdhit_archive=os.path.join(CDHIT_DIR, "{id}.tar.gz")
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
		CDHIT_WORK="$WORK/cdhit"
		mkdir -p "$OVERHANG_WORK" "$CDHIT_WORK"
		trap 'rm -rf "$WORK"' EXIT

		# stage 1: clip overhangs from the bam (needs pysam/samtools -- minibwa env)
		/home/hua575/miniconda3/envs/minibwa/bin/python {input.overhang_script} \
			{input.bam} "$OVERHANG_WORK" --sample_id {wildcards.id} --skip_is_hit_filter

		# stage 2: cluster each is_element/side pair with cd-hit-est (cd-hit env)
		printf 'is_element\tside\tn_in\tcentroids_path\tclstr_path\n' > "$CDHIT_WORK/{wildcards.id}.cdhit_manifest.tsv"
		tail -n +2 "$OVERHANG_WORK/{wildcards.id}.manifest.tsv" | while IFS=$'\t' read -r is_element side n_seqs fasta_path; do
			[ -z "$fasta_path" ] && continue
			out="$CDHIT_WORK/{wildcards.id}__${{is_element}}__${{side}}.cdhit"
			/home/hua575/miniconda3/envs/cd-hit/bin/cd-hit-est -i "$fasta_path" -o "$out" \
				-c 0.9 -n 8 -G 0 -aS 0.8 -d 0 -M 0 -T {threads}
			printf '%s\t%s\t%s\t%s\t%s\n' \
				"$is_element" "$side" "$n_seqs" "$out" "$out.clstr" >> "$CDHIT_WORK/{wildcards.id}.cdhit_manifest.tsv"
		done

		/home/hua575/miniconda3/envs/cd-hit/bin/python {input.combine_script} "$CDHIT_WORK" --output_dir "$CDHIT_WORK"

		# persist only the small, useful results. centroids_path/clstr_path in the
		# saved manifest are rewritten from the (about-to-be-deleted) tmpdir to the
		# conventional persistent-style path -- valid only after extracting
		# cd-hit/{wildcards.id}.tar.gz there, but self-documenting rather than
		# pointing at a tmpdir that no longer exists once this job ends.
		mkdir -p "{OVERHANG_DIR}/{wildcards.id}" "{CDHIT_DIR}/{wildcards.id}"
		cp "$OVERHANG_WORK/{wildcards.id}.manifest.tsv" {output.overhang_manifest}
		sed "s|$CDHIT_WORK/|{CDHIT_DIR}/{wildcards.id}/|g" \
			"$CDHIT_WORK/{wildcards.id}.cdhit_manifest.tsv" > {output.cdhit_manifest}
		cp "$CDHIT_WORK/{wildcards.id}.reads.tsv" {output.reads_tsv}

		# archive the raw intermediates -- built once, verified, then moved into place
		tar -C "$OVERHANG_WORK" -czf {output.overhang_archive}.tmp .
		tar -tzf {output.overhang_archive}.tmp > /dev/null
		mv {output.overhang_archive}.tmp {output.overhang_archive}

		tar -C "$CDHIT_WORK" -czf {output.cdhit_archive}.tmp .
		tar -tzf {output.cdhit_archive}.tmp > /dev/null
		mv {output.cdhit_archive}.tmp {output.cdhit_archive}
		"""



# BLASTs each accession's cluster representative sequences (produced by
# clip_and_cluster's reads.tsv) against the combined masked-ref-genome BLAST
# db, and classifies whether each cluster's IS-insertion junction actually
# disrupted a gene. See scripts/blast_clusters_to_ref.py's module docstring
# for the full classification logic.
#
# --tmp_dir gets its own {wildcards.id}-and-$$-scoped subdir under the job's
# $TMPDIR (same idiom as clip_and_cluster's WORK dir above) so the
# intermediate BLAST query FASTA this script writes never collides with
# another {id}'s job writing into the same shared GENE_DISRUPTION_DIR at the
# same time -- only each job's own final {id}.zip lands in the shared dir.
rule blast_clusters_to_ref:
	group: "gene_disruption"
	input:
		reads_tsv=os.path.join(CDHIT_DIR, "{id}", "{id}.reads.tsv"),
		script=os.path.join(SCRIPT_DIR, "scripts", "blast_clusters_to_ref.py"),
		blastdb_file=BLASTDB + ".nsq",
		ref_gbffs=expand(
			os.path.join(REF_GENOMES_DIR, "{acc}", "{acc}_no_IS.gbff"), acc=REF_ACCESSIONS
		)
	output:
		os.path.join(GENE_DISRUPTION_DIR, "{id}.zip")
	log:
		"logs/blast_clusters_to_ref/{id}.log"
	resources:
		runtime="10m",
		mem="1G"
	# group-components bundles 10 of these into one SLURM submission (see
	# slurmprofile/config.yaml); threads sum across the group for that one
	# sbatch call, so keep this at 1 -- 10 components * 1 thread = 10 cores
	# per submission, matching the explicit "max 10 cores/job" budget (also
	# well under the short partition's 20-CPU/job cap that a higher value
	# blew past on 2026-07-23). blastn on this small combined db is fast
	# enough single-threaded that this isn't a meaningful runtime cost.
	threads: 1
	shell:
		"""
		set -euo pipefail
		exec > {log} 2>&1

		TMP_DIR={resources.tmpdir}/blast_clusters_to_ref.{wildcards.id}.$$
		mkdir -p "$TMP_DIR"
		trap 'rm -rf "$TMP_DIR"' EXIT

		/home/hua575/miniconda3/envs/bakta/bin/python {input.script} \
			--cdhit_dirs {CDHIT_DIR}/{wildcards.id} \
			--blastdb {BLASTDB} \
			--ref_genomes_dir {REF_GENOMES_DIR} \
			--accessions {REF_ACCESSIONS_STR} \
			--output_dir {GENE_DISRUPTION_DIR} \
			--tmp_dir "$TMP_DIR" \
			--threads {threads}
		"""


# same as blast_clusters_to_ref, but against the v2 (ISfinder-BLAST-masked)
# reference genomes -- kept as a separate rule/output dir rather than
# replacing blast_clusters_to_ref so v1 and v2 results can be compared
# side-by-side. Not wired into `rule all`; run by targeting
# GENE_DISRUPTION_DIR_V2/{id}.zip explicitly.
rule blast_clusters_to_ref_v2:
	group: "gene_disruption"
	input:
		reads_tsv=os.path.join(CDHIT_DIR, "{id}", "{id}.reads.tsv"),
		script=os.path.join(SCRIPT_DIR, "scripts", "blast_clusters_to_ref.py"),
		blastdb_file=BLASTDB_V2 + ".nsq",
		ref_gbffs=expand(
			os.path.join(REF_GENOMES_DIR, "{acc}", "{acc}_no_IS_v2.gbff"), acc=REF_ACCESSIONS
		)
	output:
		os.path.join(GENE_DISRUPTION_DIR_V2, "{id}.zip")
	log:
		"logs/blast_clusters_to_ref_v2/{id}.log"
	resources:
		runtime="10m",
		mem="1G"
	threads: 1
	shell:
		"""
		set -euo pipefail
		exec > {log} 2>&1

		TMP_DIR={resources.tmpdir}/blast_clusters_to_ref_v2.{wildcards.id}.$$
		mkdir -p "$TMP_DIR"
		trap 'rm -rf "$TMP_DIR"' EXIT

		/home/hua575/miniconda3/envs/bakta/bin/python {input.script} \
			--cdhit_dirs {CDHIT_DIR}/{wildcards.id} \
			--blastdb {BLASTDB_V2} \
			--ref_genomes_dir {REF_GENOMES_DIR} \
			--accessions {REF_ACCESSIONS_STR} \
			--gbff_suffix no_IS_v2 \
			--output_dir {GENE_DISRUPTION_DIR_V2} \
			--tmp_dir "$TMP_DIR" \
			--threads {threads}
		"""


# convenience aggregate target for the v2 rerun -- not wired into `rule all`.
# Target this rule name directly (with --allowed-rules blast_clusters_to_ref_v2
# gene_disruption_v2_all) rather than expanding {id}.zip paths on the command
# line: 148k+ paths blows past the shell's argv size limit.
#
# id set: reuses the v1 GENE_DISRUPTION_DIR's already-written {id}.zip names as
# a fast proxy for "has a cd-hit reads.tsv" (a v1 zip could only have been
# written if clip_and_cluster's reads.tsv existed for that id) -- a flat
# listing of ~148k files is a single readdir, versus a glob descending into
# ~148k per-id cd-hit subdirectories, which is slow enough on this filesystem
# to make DAG-building itself take minutes (confirmed empirically 2026-08-03:
# checking every id's CDHIT_DIR reads.tsv one-by-one during DAG build alone
# took ~2 min before even reaching the one id that's actually missing it).
_v1_gene_disruption_ids = {
	f[:-len(".zip")] for f in os.listdir(GENE_DISRUPTION_DIR) if f.endswith(".zip")
} if os.path.isdir(GENE_DISRUPTION_DIR) else set()
GENE_DISRUPTION_V2_IDS = [id for id in ACCESSIONS if id in _v1_gene_disruption_ids]

rule gene_disruption_v2_all:
	input:
		expand(os.path.join(GENE_DISRUPTION_DIR_V2, "{id}.zip"), id=GENE_DISRUPTION_V2_IDS)


### ismapper #########################################################################

# rule ismapper:
#     group: "mapping"
#     input:
#         layout=os.path.join(FASTQ_DIR, "{id}.layout"),
#         ref="ismapper/ncbi_dataset/data/GCF_000005845.2/genomic.gbff",
#         is_query="ISfinder_database-master/IS.database.fa"
#     output:
#         outdir=directory(os.path.join(PATH, "ismapper", "{id}"))
#     log:
#         "logs/ismapper/log/{id}.log"
#     resources:
#         runtime="2h",
#         mem="4G"
#     conda: "/home/hua575/baymlab/sibmi/conda-envs/ismapper.yaml"
#     params:
#         fastq_dir=FASTQ_DIR
#     shell:
#         """
#         LAYOUT=$(cat {input.layout})

#         if [ "$LAYOUT" = "paired" ]; then
#             R1={params.fastq_dir}/{wildcards.id}_1.fastq
#             R2={params.fastq_dir}/{wildcards.id}_2.fastq

#             if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
#                 echo "ERROR: paired layout declared but R1/R2 not found" >> {log}
#                 exit 1
#             fi

#             ismap \
#                 --reads "$R1" "$R2" \
#                 --queries {input.is_query} \
#                 --reference {input.ref} \
#                 --output_dir {output.outdir} \
# 				--log {wildcard.id} \
#                 >> {log} 2>&1

#         elif [ "$LAYOUT" = "single" ]; then
#             R1={params.fastq_dir}/{wildcards.id}.fastq

#             if [ ! -f "$R1" ]; then
#                 echo "ERROR: single layout declared but reads not found" >> {log}
#                 exit 1
#             fi

#             ismap \
#                 --reads "$R1" \
#                 --queries {input.is_query} \
#                 --reference {input.ref} \
#                 --output_dir {output.outdir} \
# 				--log {wildcard.id} \
#                 >> {log} 2>&1

#         else
#             echo "ERROR: unrecognized layout '$LAYOUT'" >> {log}
#             exit 1
#         fi
#         """



### chunking #################################################

# import os, glob

# CHUNK_DIR = "/n/scratch/users/h/hua575/chunks500"
# PATH = "/n/scratch/users/h/hua575/v2"
# FASTQ_DIR = os.path.join(PATH, "fastq-files-v2")

# CHUNKS = [os.path.splitext(os.path.basename(f))[0]
# 	for f in glob.glob(os.path.join(CHUNK_DIR, "chunk_*.txt"))][:1]

# rule all:
# 	input:
# 		expand(os.path.join(CHUNK_DIR, "{chunk}.done"), chunk=CHUNKS)

# rule process_chunk:
# 	input:
# 		os.path.join(CHUNK_DIR, "{chunk}.txt")
# 	output:
# 		os.path.join(CHUNK_DIR, "{chunk}.done")
# 	log: "logs/{chunk}.log"
# 	resources:
# 		runtime="12h",
# 		mem="500MB"
# 	conda: "/home/hua575/baymlab/sibmi/conda-envs/ncbi-datasets.yaml"
# 	shell:
# 		"""
# 		set +e
# 		while IFS= read -r id; do

# 			# prefetch: skip if .sra already exists
# 			if [ -f {PATH}/$id/$id.sra ]; then
# 				echo "[$(date)] SKIP prefetch $id (already exists)" >> {log}
# 			else
# 				echo "[$(date)] prefetch $id" >> {log}
# 				prefetch "$id" --output-directory {PATH} 2>> {log}
# 				if [ $? -ne 0 ]; then
# 					echo "[$(date)] ERROR: prefetch failed for $id" >> {log}
# 					continue
# 				fi
# 			fi

# 			# fasterq: skip if layout already exists
# 			if [ -f {FASTQ_DIR}/${{id}}.layout ]; then
# 				echo "[$(date)] SKIP fasterq $id (already exists)" >> {log}
# 			else
# 				echo "[$(date)] fasterq $id" >> {log}
# 				fasterq-dump {PATH}/$id/$id.sra --outdir {FASTQ_DIR} 2>> {log}
# 				if [ $? -ne 0 ]; then
# 					echo "[$(date)] ERROR: fasterq failed for $id" >> {log}
# 					continue
# 				fi

# 				if [ -f {FASTQ_DIR}/${{id}}_1.fastq ]; then
# 					echo "paired" > {FASTQ_DIR}/${{id}}.layout
# 				else
# 					echo "single" > {FASTQ_DIR}/${{id}}.layout
# 				fi
# 			fi

# 		done < {input}
# 		touch {output}
# 		"""