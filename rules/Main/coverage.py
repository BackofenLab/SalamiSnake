import random
import math
import os 

###############################
## COUNTING / COVERAGE FILES ##
###############################

rule extract_alignment_ends:
	input:
		DEDUPLICAITON_OUTDIR + "/{sample}_{replicate}_sorted.bam"
	output:
		COVERAGE_OUTDIR + "/{sample}_{replicate}_alignment_ends.bed"
	conda:
		config["conda_envs"] + "/bctools.yml"
	threads: 2 
	shell:
		"if [ ! -d {COVERAGE_OUTDIR} ]; then mkdir {COVERAGE_OUTDIR}; fi"
		"&& python {config[bctools]}/extract_aln_ends.py {input} > {output}"

rule extract_crosslinking_position:
	input:
		COVERAGE_OUTDIR + "/{sample}_{replicate}_alignment_ends.bed"
	output:
		COVERAGE_OUTDIR + "/{sample}_{replicate}_crosslinking_positions.bed"
	conda:
		config["conda_envs"] + "/bctools.yml"
	threads: 2
	shell:
		"if [ ! -d {COVERAGE_OUTDIR} ]; then mkdir {COVERAGE_OUTDIR}; fi"
		"&& python {config[bctools]}/coords2clnt.py {config[extract_cl]} {input} > {output}"

# Rule to check if the end coordinate is 0. It happened in a test-run.
# For a region is has to correct start and end coordinates.
rule correct_bed_format_check:
	input:
		cl=COVERAGE_OUTDIR + "/{sample}_{replicate}_crosslinking_positions.bed",
		alends=COVERAGE_OUTDIR + "/{sample}_{replicate}_alignment_ends.bed"
	output:
		cl=COVERAGE_OUTDIR + "/{sample}_{replicate}_crosslinking_positions_check.bed",
		alends=COVERAGE_OUTDIR + "/{sample}_{replicate}_alignment_ends_check.bed"
	threads: 2
	shell:
		"""awk -F "\t" 'BEGIN {{ OFS = FS }} $3 != 0 {{ print $0 }}' {input.cl} > {output.cl} """
		"""&& awk -F "\t" 'BEGIN {{ OFS = FS }} $3 != 0 {{ print $0 }}' {input.alends} > {output.alends} """

rule sort_beds:
	input:
		cl=COVERAGE_OUTDIR + "/{sample}_{replicate}_crosslinking_positions_check.bed",
		alends=COVERAGE_OUTDIR + "/{sample}_{replicate}_alignment_ends_check.bed"
	output:
		cl=COVERAGE_OUTDIR + "/{sample}_{replicate}_crosslinking_positions_sorted.bed",
		alends=COVERAGE_OUTDIR + "/{sample}_{replicate}_alignment_ends_sorted.bed"
	threads: 2
	conda:
		config["conda_envs"] + "/bedops.yml"
	shell:
		"sort-bed {input.cl} > {output.cl}"
		"&& sort-bed {input.alends} > {output.alends}"

rule calculate_coverage:
	input:
		cl=COVERAGE_OUTDIR + "/{sample}_{replicate}_crosslinking_positions_sorted.bed",
		alends=COVERAGE_OUTDIR + "/{sample}_{replicate}_alignment_ends_sorted.bed"
	output:
		cl_pos=COVERAGE_OUTDIR + "/bedgraph/{sample}_{replicate}_crosslinking_coverage_pos_strand.bedgraph",
		cl_neg=COVERAGE_OUTDIR + "/bedgraph/{sample}_{replicate}_crosslinking_coverage_neg_strand.bedgraph",
		cl_bot=COVERAGE_OUTDIR + "/bedgraph/{sample}_{replicate}_crosslinking_coverage_both_strand.bedgraph",
		alends_pos=COVERAGE_OUTDIR + "/bedgraph/{sample}_{replicate}_alignment_ends_coverage_pos_strand.bedgraph",
		alends_ned=COVERAGE_OUTDIR + "/bedgraph/{sample}_{replicate}_alignment_ends_coverage_neg_strand.bedgraph",
		alends_bot=COVERAGE_OUTDIR + "/bedgraph/{sample}_{replicate}_alignment_ends_coverage_both_strand.bedgraph"
	threads: 4
	conda:
		config["conda_envs"] + "/bedtools.yml"
	shell:
		"if [ ! -d {COVERAGE_OUTDIR}/bedgraph ]; then mkdir {COVERAGE_OUTDIR}/bedgraph; fi"
		"&& genomeCoverageBed -i {input.cl} -g {GENOME_SIZES} -bg -strand + > {output.cl_pos}" 
		"&& genomeCoverageBed -i {input.cl} -g {GENOME_SIZES} -bg -strand - > {output.cl_neg}" 
		"&& genomeCoverageBed -i {input.cl} -g {GENOME_SIZES} -bg > {output.cl_bot}" 
		"&& genomeCoverageBed -i {input.alends} -g {GENOME_SIZES} -bg -strand + > {output.alends_pos}"
		"&& genomeCoverageBed -i {input.alends} -g {GENOME_SIZES} -bg -strand - > {output.alends_ned}" 
		"&& genomeCoverageBed -i {input.alends} -g {GENOME_SIZES} -bg > {output.alends_bot}" 

rule calculate_coverage_bigwig:
	input:
		cl=COVERAGE_OUTDIR + "/bedgraph/{sample}_{replicate}_crosslinking_coverage_{type}_strand.bedgraph",
		alends=COVERAGE_OUTDIR + "/bedgraph/{sample}_{replicate}_alignment_ends_coverage_{type}_strand.bedgraph"
	output:
		cl=COVERAGE_OUTDIR + "/bigwig/{sample}_{replicate}_crosslinking_coverage_{type}_strand.bigwig",
		alends=COVERAGE_OUTDIR + "/bigwig/{sample}_{replicate}_alignment_ends_coverage_{type}_strand.bigwig"
	threads: 4
	conda:
		config["conda_envs"] + "/bedGraphToBigWig.yml"
	shell:
		"if [ ! -d {COVERAGE_OUTDIR}/bigwig ]; then mkdir {COVERAGE_OUTDIR}/bigwig; fi "
		"&& bedGraphToBigWig {input.cl} {GENOME_SIZES} {output.cl}"
		"&& bedGraphToBigWig {input.alends} {GENOME_SIZES} {output.alends}"