import random
import math
import os 

#################
## PEAKCALLING ##
#################

# BE CAREFUL when using other protocols then you might need to define a different bitflag!!
rule mate_reads_fitlering:
	input:
		DEDUPLICAITON_OUTDIR + "/{sample}_{replicate}_sorted.bam"
	output:
		bam=MATEFILTER_OUTDIR + "/{sample}_{replicate}.bam",
		sorted_bam=MATEFILTER_OUTDIR + "/{sample}_{replicate}_sorted.bam",
		sorted_bam_bai=MATEFILTER_OUTDIR + "/{sample}_{replicate}_sorted.bam.bai"
	threads: 2
	conda:
		config["conda_envs"] + "/samtools.yml"
	shell:
		"if [ ! -d {MATEFILTER_OUTDIR} ]; then mkdir {MATEFILTER_OUTDIR}; fi"
		"&& samtools view -b -f 0x0040 {input} > {output.bam}"
		"&& samtools sort {output.bam} > {output.sorted_bam}"
		"&& samtools index {output.sorted_bam}"

if ( control == "yes" ):

	rule pureclip:
		input:
			experiment=MATEFILTER_OUTDIR + "/{sample_exp}_{replicate_exp}_sorted.bam",
			experiment_bai=MATEFILTER_OUTDIR + "/{sample_exp}_{replicate_exp}_sorted.bam.bai",
			control=MATEFILTER_OUTDIR + "/{sample_ctl}_{replicate_ctl}_sorted.bam",
			control_bai=MATEFILTER_OUTDIR + "/{sample_ctl}_{replicate_ctl}_sorted.bam.bai",
			genome_fasta=GENOME_FASTA
		output:
			crosslinking_sites=PEAKCALLING_OUTDIR + "/{sample_exp}_{replicate_exp}_{sample_ctl}_{replicate_ctl}_crosslinkind_sites.bed",
			binding_regions=PEAKCALLING_OUTDIR + "/{sample_exp}_{replicate_exp}_{sample_ctl}_{replicate_ctl}_binding_regions.bed",
		threads: 10
		params:
			tmp=PEAKCALLING_OUTDIR + "/tmp/",
			parameters=PEAKCALLING_OUTDIR + "/{sample_exp}_{replicate_exp}_{sample_ctl}_{replicate_ctl}_parameters.txt"
		conda:
			config["conda_envs"] + "/pureclip.yml"
		shell:
			"if [ ! -d {PEAKCALLING_OUTDIR} ]; then mkdir {PEAKCALLING_OUTDIR}; fi "
			"&& echo {config[pureclip]} >> {file_tool_params}"
			"&& pureclip -i {input.experiment} -bai {input.experiment_bai} -g {input.genome_fasta} -o {output.crosslinking_sites} -tmp {params.tmp} "
			"-ibam {input.control} -ibai {input.control_bai} -or {output.binding_regions} -p {params.parameters} -nt {threads} {config[pureclip]} "

	rule peaks_extend_frontiers:
		input:
			bed=PEAKCALLING_OUTDIR + "/{sample_exp}_{replicate_exp}_{sample_ctl}_{replicate_ctl}_binding_regions.bed",
			genome=GENOME_SIZES
		output:
			PEAKCALLING_OUTDIR + "/{sample_exp}_{replicate_exp}_{sample_ctl}_{replicate_ctl}_peaks_extended.bed"
		threads: 2
		conda:
			config["conda_envs"] + "/bedtools.yml"
		shell:
			"echo {config[peaks_extend_frontiers]} >> {file_tool_params}"
			"&& bedtools slop {config[peaks_extend_frontiers]} -i {input.bed} -g {input.genome} > {output}"
else:

	rule piranha:
		input:
			MATEFILTER_OUTDIR + "/{sample}_{replicate}_sorted.bam"
		output:
			PEAKCALLING_OUTDIR + "/{sample}_{replicate}_binding_regions.bed"
		conda:
			config["conda_envs"] + "/piranha.yml"
		threads: 2 
		shell:
			"if [ ! -d {PEAKCALLING_OUTDIR} ]; then mkdir {PEAKCALLING_OUTDIR}; fi"
			"&& echo {config[piranha]} >> {file_tool_params}"
			"&& Piranha {config[piranha]} {input} -o {output}"	

	# rule pureclip:
	# 	input:
	# 		experiment=DEDUPLICAITON_OUTDIR + "/{sample}_{replicate}_sorted.bam",
	# 		experiment_bai=DEDUPLICAITON_OUTDIR + "/{sample}_{replicate}_sorted.bam.bai",
	# 		genome_fasta=GENOME_FASTA
	# 	output:
	# 		crosslinking_sites=PEAKCALLING_OUTDIR + "/{sample}_{replicate}_crosslinkind_sites.bed",
	# 		binding_regions=PEAKCALLING_OUTDIR + "/{sample}_{replicate}_binding_regions.bed",
	# 	threads: 10
	# 	params:
	# 		tmp=PEAKCALLING_OUTDIR + "/tmp/",
	# 		parameters=PEAKCALLING_OUTDIR + "/{sample}_{replicate}_parameters.txt"
	# 	conda:
	# 		config["conda_envs"] + "/pureclip.yml"
	# 	shell:
	# 		"if [ ! -d {PEAKCALLING_OUTDIR} ]; then mkdir {PEAKCALLING_OUTDIR}; fi "
	#		"&& echo {config[pureclip]} >> {file_tool_params}"
	# 		"&& pureclip -i {input.experiment} -bai {input.experiment_bai} -g {input.genome_fasta} -o {output.crosslinking_sites} -tmp {params.tmp} "
	# 		"-or {output.binding_regions} -p {params.parameters} -nt {threads} {config[pureclip]} "

	rule peaks_extend_frontiers:
		input:
			bed=PEAKCALLING_OUTDIR + "/{sample}_{replicate}_binding_regions.bed",
			genome=GENOME_SIZES
		output:
			PEAKCALLING_OUTDIR + "/{sample}_{replicate}_peaks_extended.bed"
		threads: 2
		conda:
			config["conda_envs"] + "/bedtools.yml"
		shell:
			"echo {config[peaks_extend_frontiers]} >> {file_tool_params}"
			"&& bedtools slop {config[peaks_extend_frontiers]} -i {input.bed} -g {input.genome} > {output}"

# rule peakachu:
#     input:
#     	clip=expand(DEDUPLICAITON_OUTDIR + "/{sample}_sorted.bam", sample=REPLICATES_CLIP),
#     	control=expand(DEDUPLICAITON_OUTDIR + "/{sample}_sorted.bam", sample=REPLICATES_CONTROL)
#     output:
#     	peaks_tsv=PEAKCALLING_OUTDIR + "/peakachu.tsv",
#     	peaks_gtf=PEAKCALLING_OUTDIR + "/peakachu.gtf",
#     	blockbuster=PEAKCALLING_OUTDIR + "/blockbuster.bed"
#     threads: 4
#     conda:
#     	config["conda_envs"] + "/peakachu.yml"
#     shell:
#     	#"--max_proc "${GALAXY_SLOTS:-1}"
#     	"if [ ! -d {PEAKCALLING_OUTDIR} ]; then mkdir {PEAKCALLING_OUTDIR}; fi"
#		"&& echo {config[peakachu]} >> {file_tool_params}"
#     	"&& source activate peakachu"
#     	"&& peakachu adaptive --exp_libs {input.clip} --ctr_libs {input.control} {config[peakachu]} --output_folder {PEAKCALLING_OUTDIR} "
#     	"&& source deactivate "
#     	"&& head -q -n 1 {PEAKCALLING_OUTDIR}/peak_tables/*.csv > tmp.tsv "
#     	"&& head -q -n 1 tmp.tsv > {output.peaks_tsv}"
#     	"&& rm tmp.tsv"
#     	"&& tail -n +2 -q {PEAKCALLING_OUTDIR}/peak_tables/*.csv >> {output.peaks_tsv} "
#     	"&& cat {PEAKCALLING_OUTDIR}/peak_annotations/*.gff | awk '/peak/ {{ print $0 }}' > {output.peaks_gtf} "
#     	"&& cat {PEAKCALLING_OUTDIR}/blockbuster_input/*bed > {output.blockbuster}"

# rule peak_calling_quality:
# 	input:
# 		peaks=PEAKCALLING_OUTDIR + "/peakachu.gtf",
# 		clip=DEDUPLICAITON_OUTDIR + "/{sample}_{replicate}_sorted.bam"
# 	output:
# 		htseq_hits=PEAKCALLING_OUTDIR + "/{sample}_{replicate}_htseq_hits.txt",
# 		htseq_nohits=PEAKCALLING_OUTDIR + "/{sample}_{replicate}_htseq_nohits.txt",
# 		reads_in_peaks=PEAKCALLING_OUTDIR + "/{sample}_{replicate}_reads_summary.txt"
# 	threads: 2
# 	conda:
# 		config["conda_envs"] + "/htseq.yml"
# 	shell:
# 		"htseq-count --mode=union --stranded=yes --minaqual=10 --type='peak_region' --idattr='ID' --order=name --format=bam {input.clip} {input.peaks} "
# 		"""| awk '{{ if ($1 ~ "no_feature|ambiguous|too_low_aQual|not_aligned|alignment_not_unique") print $0 | "cat 1>&2"; else print $0 }}' > {output.htseq_hits} 2> {output.htseq_nohits} """
# 		"""&& awk 'FNR==NR{{ SUM1+=$2; next }} {{ SUM2+=$2 }} END {{ print "#reads in peaks: "SUM1; """
# 		"""print "#culled reads: "SUM2; print "percentage reads in peaks: " SUM1/(SUM1+SUM2) }}' {output.htseq_hits} {output.htseq_nohits} > {output.reads_in_peaks}"""