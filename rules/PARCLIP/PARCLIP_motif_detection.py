

#####################
## MOTIF DETECTION ##
#####################

rule extract_genomic_DNA_dreme:
	input:
		ROBUSTPEAKS_OUTDIR + "/robust_between_all.bed"
	output:
		MOTIF_DETECTION_OUTDIR + "/peaks.fa"
	threads: 2
	shell:
		"if [ ! -d {MOTIF_DETECTION_OUTDIR} ]; then mkdir {MOTIF_DETECTION_OUTDIR}; fi"
		"""&& awk '{{ print $1"\t"$2"\t"$3 }}' {input} > {MOTIF_DETECTION_OUTDIR}/peaks.bed3 """
		"&& python " + config["extract_genomic_dna"] + "/fetch_DNA_sequence.py -o {output} {MOTIF_DETECTION_OUTDIR}/peaks.bed3 {GENOME_FASTA}"

rule meme_chip:
	input: 
		MOTIF_DETECTION_OUTDIR + "/peaks.fa"
	output:
		MOTIF_DETECTION_OUTDIR + "/meme_chip/meme-chip.html"
	threads: 4
	conda:
		"envs/meme_suite.yml"
	shell:
		"meme-chip {input} -oc {MOTIF_DETECTION_OUTDIR}/meme_chip {config[meme_chip]}"