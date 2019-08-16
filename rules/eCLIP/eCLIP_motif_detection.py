

#####################
## MOTIF DETECTION ##
#####################

rule extract_genomic_DNA_dreme:
	input:
		PEAKCALLING_OUTDIR + "/{sample}_{replicate}_robust_peaks.bed"
	output:
		MOTIF_DETECTION_OUTDIR + "/{sample}_{replicate}_peaks.fa"
	threads: 2
	shell:
		"if [ ! -d {MOTIF_DETECTION_OUTDIR} ]; then mkdir {MOTIF_DETECTION_OUTDIR}; fi"
		"""&& awk '$9<0.05 {{ print $1"\t"$2"\t"$3 }}' {input} > {MOTIF_DETECTION_OUTDIR}/peaks.bed3 """
		"&& python " + config["extract_genomic_dna"] + "/fetch_DNA_sequence.py -o {output} {MOTIF_DETECTION_OUTDIR}/peaks.bed3 {GENOME_FASTA}"

rule meme_chip:
	input: 
		MOTIF_DETECTION_OUTDIR + "/{sample}_{replicate}_peaks.fa"
	output:
		MOTIF_DETECTION_OUTDIR + "/{sample}_{replicate}_meme_chip/meme-chip.html"
	threads: 4
	conda:
		config["conda_envs"] + "/meme_suite.yml"
	params:
		outdir=MOTIF_DETECTION_OUTDIR + "/{sample}_{replicate}_meme_chip"
	shell:
		"echo {config[meme_chip]} >> {file_tool_params}"
		"&& meme-chip -meme-p {threads} {config[meme_chip]} -oc {params.outdir} {input}"