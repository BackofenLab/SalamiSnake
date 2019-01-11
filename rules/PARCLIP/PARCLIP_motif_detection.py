

#####################
## MOTIF DETECTION ##
#####################

rule extract_genomic_DNA_dreme:
	input:
		PEAKCALLING_OUTDIR + "/{sample}_{replicate}_peaks_extended.bed"
	output:
		MOTIF_DETECTION_OUTDIR + "/{sample}_{replicate}_peaks.fa"
	threads: 2
	shell:
		"if [ ! -d {MOTIF_DETECTION_OUTDIR} ]; then mkdir {MOTIF_DETECTION_OUTDIR}; fi"
		"""&& awk '{{ print $1"\t"$2"\t"$3 }}' {input} > {MOTIF_DETECTION_OUTDIR}/peaks.bed3 """
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

# ##################
# ## MOTIF SEARCH ##
# ##################

# rule fimo_for_dreme_output:
#  	input: 
#  		MOTIF_DETECTION_OUTDIR + "/dreme/dreme.xml"
#  	output:
#  		html=MOTIF_SEARCH_OUTDIR + "/fimo_dreme/fimo.html"
#  	threads: 4
#  	conda:
#  		config["conda_envs"] + "/meme_suite.yml"
#  	shell:
#  		"if [ ! -d {MOTIF_SEARCH_OUTDIR} ]; then mkdir {MOTIF_SEARCH_OUTDIR}; fi"
#  		"&& fimo --alpha 1.000000 --max-stored-scores 100000 --motif-pseudo 0.100000 --qv-thresh --thresh 0.000100 --verbosity 1 "
#  		"--oc {MOTIF_SEARCH_OUTDIR}/fimo_dreme {input} {GENOME_FASTA}"

rule fimo_for_meme_output:
 	input: 
 		MOTIF_DETECTION_OUTDIR + "/{sample}_{replicate}_meme_chip/meme_out/meme.xml"
 	output:
 		MOTIF_SEARCH_OUTDIR + "/fimo_meme/{sample}_{replicate}_meme"
 	threads: 4
 	conda:
 		config["conda_envs"] + "/meme_suite.yml"
 	shell:
 		"fimo --alpha 1.000000 --max-stored-scores 100000 --motif-pseudo 0.100000 --qv-thresh --thresh 0.000100 --verbosity 1 "
 		"--oc {output} {input} {GENOME_FASTA}"