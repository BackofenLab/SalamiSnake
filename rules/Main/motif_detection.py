
if ( control == "yes" ):

	if ( peakcaller == "PureCLIP" ):

		#####################
		## MOTIF DETECTION ##
		#####################

		rule extract_genomic_DNA_dreme:
			input:
				PEAKCALLING_OUTDIR + "/{sample_exp}_{replicate_exp}_{sample_ctl}_{replicate_ctl}_peaks_extended.bed"
			output:
				bed3=MOTIF_DETECTION_OUTDIR + "/{sample_exp}_{replicate_exp}_{sample_ctl}_{replicate_ctl}_peaks.bed3",
				fa=MOTIF_DETECTION_OUTDIR + "/{sample_exp}_{replicate_exp}_{sample_ctl}_{replicate_ctl}_peaks.fa"
			threads: 2
			shell:
				"if [ ! -d {MOTIF_DETECTION_OUTDIR} ]; then mkdir {MOTIF_DETECTION_OUTDIR}; fi"
				"""&& awk '{{ print $1"\t"$2"\t"$3 }}' {input} > {output.bed3} """
				"&& python " + config["extract_genomic_dna"] + "/fetch_DNA_sequence.py -o {output.fa} {output.bed3} {GENOME_FASTA}"

		rule meme_chip:
			input: 
				MOTIF_DETECTION_OUTDIR + "/{sample_exp}_{replicate_exp}_{sample_ctl}_{replicate_ctl}_peaks.fa"
			output:
				MOTIF_DETECTION_OUTDIR + "/{sample_exp}_{replicate_exp}_{sample_ctl}_{replicate_ctl}_meme_chip/meme-chip.html"
			threads: 4
			conda:
				config["conda_envs"] + "/meme_suite.yml"
			params:
				outdir=MOTIF_DETECTION_OUTDIR + "/{sample_exp}_{replicate_exp}_{sample_ctl}_{replicate_ctl}_meme_chip"
			shell:
				"echo {config[meme_chip]} >> {file_tool_params}"
				"&& meme-chip -meme-p {threads} {config[meme_chip]} -oc {params.outdir} {input}"

		##################
		## MOTIF SEARCH ##
		##################

		rule fimo_for_dreme_output:
		 	input: 
		 		MOTIF_DETECTION_OUTDIR + "/{sample}_{replicate}_meme_chip/dreme_out/dreme.xml"
		 	output:
		 		MOTIF_SEARCH_OUTDIR + "/fimo_dreme/{sample_exp}_{replicate_exp}_{sample_ctl}_{replicate_ctl}_dreme"
		 	threads: 4
		 	conda:
		 		config["conda_envs"] + "/meme_suite.yml"
		 	shell:
		 		"if [ ! -d {MOTIF_SEARCH_OUTDIR} ]; then mkdir {MOTIF_SEARCH_OUTDIR}; fi"
		 		"fimo {config[fimo]} --oc {output} {input} {GENOME_FASTA}"

		rule fimo_for_meme_output:
		 	input: 
		 		MOTIF_DETECTION_OUTDIR + "/{sample_exp}_{replicate_exp}_{sample_ctl}_{replicate_ctl}_meme_chip/meme_out/meme.xml"
		 	output:
		 		MOTIF_SEARCH_OUTDIR + "/fimo_meme/{sample_exp}_{replicate_exp}_{sample_ctl}_{replicate_ctl}_meme"
		 	threads: 4
		 	conda:
		 		config["conda_envs"] + "/meme_suite.yml"
		 	shell:
		 		"fimo {config[fimo]} --oc {output} {input} {GENOME_FASTA}"

else:

	if ( peakcaller == "PureCLIP" or peakcaller == "Piranha" ):
		#####################
		## MOTIF DETECTION ##
		#####################

		rule extract_genomic_DNA_dreme:
			input:
				PEAKCALLING_OUTDIR + "/{sample}_{replicate}_peaks_extended.bed"
			output:
				bed3=MOTIF_DETECTION_OUTDIR + "/{sample}_{replicate}_peaks.bed3",
				fa=MOTIF_DETECTION_OUTDIR + "/{sample}_{replicate}_peaks.fa"
			threads: 2
			shell:
				"if [ ! -d {MOTIF_DETECTION_OUTDIR} ]; then mkdir {MOTIF_DETECTION_OUTDIR}; fi"
				"""&& awk '{{ print $1"\t"$2"\t"$3 }}' {input} > {output.bed3} """
				"&& python " + config["extract_genomic_dna"] + "/fetch_DNA_sequence.py -o {output.fa} {output.bed3} {GENOME_FASTA}"

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

		##################
		## MOTIF SEARCH ##
		##################

		rule fimo_for_dreme_output:
		 	input: 
		 		MOTIF_DETECTION_OUTDIR + "/{sample}_{replicate}_meme_chip/dreme_out/dreme.xml"
		 	output:
		 		MOTIF_SEARCH_OUTDIR + "/fimo_dreme/{sample}_{replicate}_dreme"
		 	threads: 4
		 	conda:
		 		config["conda_envs"] + "/meme_suite.yml"
		 	shell:
		 		"if [ ! -d {MOTIF_SEARCH_OUTDIR} ]; then mkdir {MOTIF_SEARCH_OUTDIR}; fi"
				"fimo {config[fimo]} --oc {output} {input} {GENOME_FASTA}"

		rule fimo_for_meme_output:
		 	input: 
		 		MOTIF_DETECTION_OUTDIR + "/{sample}_{replicate}_meme_chip/meme_out/meme.xml"
		 	output:
		 		MOTIF_SEARCH_OUTDIR + "/fimo_meme/{sample}_{replicate}_meme"
		 	threads: 4
		 	conda:
		 		config["conda_envs"] + "/meme_suite.yml"
		 	shell:
		 		"fimo {config[fimo]} --oc {output} {input} {GENOME_FASTA}"