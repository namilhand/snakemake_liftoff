# Snakemake workflow for lifting over crossover sites from TAIR10 to Col-CC

# Usage
# conda env create --file environment.yaml --name GBS
# conda activate GBS
# snakemake -p --cores 30
# conda deactivate

#======= Note ============
# input cotable file name should be formatted in
# {sample}_cotable.txt
#=========================

import pandas as pd
import os

# To make a shell script invoked using the "shell" directive,
# we need to determine the base path of Snakefile since we expect the script directory to be there as well.
SRCDIR = srcdir("")

# Specify config file parameters
configfile:"config.yaml"
SAMPLE = config["SAMPLE"]
CHROM_MATCH = config["LIFTOFF"]["chroms"]
QUERY_GENOME = config["GENOME"]["query_genome"]
TARGET_GENOME = config["GENOME"]["target_genome"]
QUERY_GENOME_NAME = config["GENOME"]["query_genome_name"]
TARGET_GENOME_NAME = config["GENOME"]["target_genome_name"]
REF_FAI = config["GENOME"]["target_fai"]

binSize = config["TSV"]["binSize"]
binName = config["TSV"]["binName"]


# Specify the desired end target file(s)

rule all:
	input:
		expand("results/01_cotable_gff/{sample}.gff", sample = SAMPLE),
		expand("results/02_liftoff/{sample}_{target_genome_name}.gff", sample = SAMPLE, target_genome_name = TARGET_GENOME_NAME),
		expand("results/03_cotable/{sample}_{target_genome_name}_cotable.txt", sample=SAMPLE, target_genome_name = TARGET_GENOME_NAME),
		expand("results/04_genomeBin/{sample}_{target_genome_name}_genomeBin{binName}.tsv", sample=SAMPLE, target_genome_name=TARGET_GENOME_NAME, binName = binName)


rule cotable2gff:
	"transform cotable.txt to cotable.gff"
	output:
		"results/01_cotable_gff/{sample}.gff"
	input:
		"input/{sample}_cotable.txt"	
	shell:
		"Rscript src/cotable2bed.R {input} {output}"

rule liftoff2colcen:
	"liftoff cosite from query_genome to target-genome"
	output:
		gff="results/02_liftoff/{sample}_{TARGET_GENOME_NAME}.gff",
		unmapped="results/02_liftoff/unmapped/{sample}_{TARGET_GENOME_NAME}_unmapped.gff",
		interm=directory("results/02_liftoff/intermediate_files/{sample}_{TARGET_GENOME_NAME}_interm")
	input:
		"results/01_cotable_gff/{sample}.gff"
	shell:
		r"""
		liftoff -g {input} -o {output.gff} -u {output.unmapped} \
				-dir {output.interm} \
				-chroms {CHROM_MATCH} \
				{TARGET_GENOME} {QUERY_GENOME}
		"""

rule gff2cotable:
	"gff to cotable format"
	"Note: start and end positions set as cos +/- 1"
	output:
		"results/03_cotable/{sample}_{target_genome_name}_cotable.txt"
	input:
		"results/02_liftoff/{sample}_{target_genome_name}.gff"
	shell:
		r"""
		Rscript src/gff2cotable.R {input} {output}
		"""

rule toGenomeBin:
	output: "results/04_genomeBin/{sample}_{target_genome_name}_genomeBin{binName}.tsv"
	input: "results/03_cotable/{sample}_{target_genome_name}_cotable.txt"
	shell:
		r"""
		Rscript src/cotableToGenomeBinTSV.R {input} {binSize} {REF_FAI} {output}
		"""

