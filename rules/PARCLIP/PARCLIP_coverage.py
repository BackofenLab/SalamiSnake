
###############################
## COUNTING / COVERAGE FILES ##
###############################

rule bam_to_bed:
	input:
		PRE_FOR_UMI_OUTDIR + "/{sample}_{replicate}_got_umis_unlocalized_check.bam"
	output:
		COVERAGE_OUTDIR + "/{sample}_{replicate}.bed"
	threads: 2
	conda:
		config["conda_envs"] + "/bedtools.yml"
	shell:
		"if [ ! -d {COVERAGE_OUTDIR} ]; then mkdir {COVERAGE_OUTDIR}; fi"
		"&& bedtools bamtobed -i {input} > {output}"

rule correct_bed_format_check:
	input:
		COVERAGE_OUTDIR + "/{sample}_{replicate}.bed"
	output:
		COVERAGE_OUTDIR + "/{sample}_{replicate}_check.bed"
	threads: 2
	shell:
		"""awk -v OFS="\t" -v FS="\t" '$3 != 0 {{ print $0 }}' {input} > {output} """

rule sort_beds:
	input:
		COVERAGE_OUTDIR + "/{sample}_{replicate}_check.bed"
		#COVERAGE_OUTDIR + "/{sample}_{replicate}.bed"
	output:
		COVERAGE_OUTDIR + "/{sample}_{replicate}_check_sorted.bed"
	threads: 2
	conda:
		#config["conda_envs"] + "/bedtools.yml"
		config["conda_envs"] + "/bedops.yml"
	shell:
		"&& sort-bed {input} > {output}"
		#"&& bedtools sort -i {input} > {output}"

rule calculate_coverage:
	input:
		COVERAGE_OUTDIR + "/{sample}_{replicate}_check_sorted.bed"
	output:
		pos=COVERAGE_OUTDIR + "/bedgraph/{sample}_{replicate}_coverage_pos_strand.bedgraph",
		neg=COVERAGE_OUTDIR + "/bedgraph/{sample}_{replicate}_coverage_neg_strand.bedgraph",
		bot=COVERAGE_OUTDIR + "/bedgraph/{sample}_{replicate}_coverage_both_strand.bedgraph"
	threads: 4
	conda:
		config["conda_envs"] + "/bedtools.yml"
	shell:
		"if [ ! -d {COVERAGE_OUTDIR}/bedgraph ]; then mkdir {COVERAGE_OUTDIR}/bedgraph; fi"
		"&& genomeCoverageBed -i {input} -g {GENOME_SIZES} -bg -strand + > {output.pos}" 
		"&& genomeCoverageBed -i {input} -g {GENOME_SIZES} -bg -strand - > {output.neg}" 
		"&& genomeCoverageBed -i {input} -g {GENOME_SIZES} -bg > {output.bot}" 

rule calculate_coverage_bigwig:
	input:
		COVERAGE_OUTDIR + "/bedgraph/{sample}_{replicate}_coverage_{type}_strand.bedgraph"
	output:
		COVERAGE_OUTDIR + "/bigwig/{sample}_{replicate}_coverage_{type}_strand.bigwig"
	threads: 4
	conda:
		config["conda_envs"] + "/bedGraphToBigWig.yml"
	shell:
		"if [ ! -d {COVERAGE_OUTDIR}/bigwig ]; then mkdir {COVERAGE_OUTDIR}/bigwig; fi "
		"&& bedGraphToBigWig {input} {GENOME_SIZES} {output}"
