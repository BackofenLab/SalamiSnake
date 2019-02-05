import random
import math
import os 

#####################
## PEAK ANNOTATION ##
#####################

# PURECLIP
rule intersect_binding_regions_with_peaks:
	input:
		binding_regions=PEAKCALLING_OUTDIR + "/{sample_exp}_{replicate_exp}_{sample_ctl}_{replicate_ctl}_binding_regions.bed",
		annotation=GENOME_GTF
	output:
		ANNOTATION_PEAKS_OUTDIR + "/{sample_exp}_{replicate_exp}_{sample_ctl}_{replicate_ctl}_binding_regions_intersecting_peaks.gtf"
	threads: 2
	conda:
		config["conda_envs"] + "/bedtools.yml"	
	shell:
		"if [ ! -d {ANNOTATION_PEAKS_OUTDIR} ]; then mkdir {ANNOTATION_PEAKS_OUTDIR}; fi "
		"&& echo {config[annotation]} >> {file_tool_params}"
		"&& bedtools intersect -a {input.annotation} -b {input.binding_regions} {config[annotation]} > {output}"

# PEAKACHU
# rule intersect_binding_regions_with_peaks:
#     input:
#     	binding_regions=PEAKCALLING_OUTDIR + "/{sample_exp}_{replicate_exp}_{sample_ctl}_{replicate_ctl}/peakachu.gtf",
#     	annotation=GENOME_GTF
#     output:
#     	ANNOTATION_PEAKS_OUTDIR + "/{sample_exp}_{replicate_exp}_{sample_ctl}_{replicate_ctl}_binding_regions_intersecting_peaks.gtf"
#     threads: 2
#     conda:
# 		config["conda_envs"] + "/bedtools.yml"	
#     shell:
#     	"if [ ! -d {ANNOTATION_PEAKS_OUTDIR} ]; then mkdir {ANNOTATION_PEAKS_OUTDIR}; fi "
#		"&& echo {config[annotation]} >> {file_tool_params}"
#     	"&& bedtools intersect -a {input.annotation} -b {input.binding_regions} {config[annotation]} > {output}"

####################
## HTSEQ ANALYSIS ##
####################

rule peak_calling_quality:
	input:
		reads=DEDUPLICAITON_OUTDIR + "/{sample}_{replicate}_name_sorted.bam",
		annotation=GENOME_GTF
	output:
		htseq_hits=HTSEQ_COUNT_OUTDIR + "/{sample}_{replicate}_htseq_hits.txt",
		htseq_nohits=HTSEQ_COUNT_OUTDIR + "/{sample}_{replicate}_htseq_nohits.txt",
		reads_in_features=HTSEQ_COUNT_OUTDIR + "/{sample}_{replicate}_reads_summary.txt"
	threads: 2
	conda:
		config["conda_envs"] + "/htseq.yml"
	shell:
		"if [ ! -d {HTSEQ_COUNT_OUTDIR} ]; then mkdir {HTSEQ_COUNT_OUTDIR}; fi "
		"&& htseq-count --mode=union --stranded=yes --minaqual=10 --type='transcript' --idattr='gene_id' --nonunique=none --secondary-alignments=ignore "
		"--supplementary-alignments=ignore --order=name --format=bam {input.reads} {input.annotation}"
		"""| awk '{{ if ($1 ~ "no_feature|ambiguous|too_low_aQual|not_aligned|alignment_not_unique") print $0 | "cat 1>&2"; else print $0 }}' > {output.htseq_hits} 2> {output.htseq_nohits} """
		"""&& awk 'FNR==NR{{ SUM1+=$2; next }} {{ SUM2+=$2 }} END {{ print "#reads in features: "SUM1; """
		"""print "#culled reads: "SUM2; print "percentage reads in features: " SUM1/(SUM1+SUM2) }}' {output.htseq_hits} {output.htseq_nohits} > {output.reads_in_features}"""
