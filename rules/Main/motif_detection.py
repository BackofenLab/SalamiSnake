

if ( control == "yes" ):

	if ( peakcaller == "PanPeaker" ):


		#####################
		## MOTIF DETECTION ##
		#####################

		rule extract_genomic_DNA_dreme:
			input:
				PEAKCALLING_OUTDIR + "/robust_peaks_refined.bed"
			output:
				fa=MOTIF_DETECTION_OUTDIR + "/robust_peaks_refined.fa"
			threads: 2
			conda:
				config["conda_envs"] + "/fetch_DNA_sequence.yml"			
			shell:
				"if [ ! -d {MOTIF_DETECTION_OUTDIR} ]; then mkdir {MOTIF_DETECTION_OUTDIR}; fi"
				"&& python " + config["extract_genomic_dna"] + "/fetch_DNA_sequence.py -o {output.fa} {input} {GENOME_FASTA}"

		rule meme_chip:
			input: 
				MOTIF_DETECTION_OUTDIR + "/robust_peaks_refined.fa"
			output:
				MOTIF_DETECTION_OUTDIR + "/robust_peaks_refined_meme_chip/meme-chip.html"
			threads: 4
			conda:
				config["conda_envs"] + "/meme_suite.yml"
			params:
				outdir=MOTIF_DETECTION_OUTDIR + "/robust_peaks_refined_meme_chip"
			shell:
				"echo {config[meme_chip]} >> {file_tool_params}"
				"&& meme-chip -meme-p {threads} {config[meme_chip]} -oc {params.outdir} {input}"

	if ( peakcaller == "PureCLIP" ):

		#####################
		## MOTIF DETECTION ##
		#####################

		rule extract_genomic_DNA_dreme:
			input:
				PEAKCALLING_OUTDIR + "/{sample_exp}_{replicate_exp}_{sample_ctl}_{replicate_ctl}_peaks_extended.bed"
			output:
				fa=MOTIF_DETECTION_OUTDIR + "/{sample_exp}_{replicate_exp}_{sample_ctl}_{replicate_ctl}_peaks.fa"
			threads: 2
			conda:
				config["conda_envs"] + "/fetch_DNA_sequence.yml"			
			shell:
				"if [ ! -d {MOTIF_DETECTION_OUTDIR} ]; then mkdir {MOTIF_DETECTION_OUTDIR}; fi"
				"&& python " + config["extract_genomic_dna"] + "/fetch_DNA_sequence.py -o {output.fa} {input} {GENOME_FASTA}"

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

		#####################
		## TARGET ANALYSIS ##
		#####################

		rule sort_peaks:
			input:
				expand(PEAKCALLING_OUTDIR + "/{sample_exp}_{replicate_exp}_{sample_ctl}_{replicate_ctl}_peaks_extended.bed",
					sample_exp=SAMPLES[0], replicate_exp=REP_NAME_CLIP, sample_ctl=SAMPLES[1], replicate_ctl=REP_NAME_CONTROL)
			output:
				PEAKCALLING_OUTDIR + "/robust_peaks_sorted.bed"
			threads: 2
			conda:
				config["conda_envs"] + "/bedtools.yml"	
			params:
				tmp=PEAKCALLING_OUTDIR + "/robust_peaks_refined_tmp.bed"
			shell:
				"cat {input} > {params.tmp}"
				"&& bedtools sort -i {params.tmp} > {output}"
				"&& rm {params.tmp}"			

		rule robust_peaks_generation:
			input:
				PEAKCALLING_OUTDIR + "/robust_peaks_sorted.bed"
			output:
				PEAKCALLING_OUTDIR + "/robust_peaks_refined.bed"
			threads: 2
			conda:
				config["conda_envs"] + "/bedtools.yml"	
			params:
				tmp=PEAKCALLING_OUTDIR + "/robust_peaks_refined_tmp.bed"
			shell:
				"bedtools merge -i {input} -s -c 6 -o distinct > {params.tmp}"
				"""&& awk -v OFS="\t" '{{print $1,$2,$3,NR,NR,$4}}' {params.tmp} > {output}"""
				"&& rm {params.tmp}"
		
		bam_folder = MATEFILTER_OUTDIR
		bam_suffix = "_sorted"
		if ( PROTOCOL == "PARCLIP" ):
			bam_folder = PRE_FOR_UMI_OUTDIR
			bam_suffix = "_got_umis_unlocalized_check"

		rule stoatydive:
			input:
				peaks=PEAKCALLING_OUTDIR + "/robust_peaks_refined.bed",
				bam=bam_folder + "/{sample}_{replicate}" + bam_suffix + ".bam",
				genome=GENOME_SIZES
			output:
				MOTIF_DETECTION_OUTDIR + "/stoatydive_with_length_norm/VC_Distribution_{sample}_{replicate}" + bam_suffix + ".pdf"
			threads: 2
			conda:
				config["conda_envs"] + "/stoatydive.yml"
			params:
				output_folder=MOTIF_DETECTION_OUTDIR + "/stoatydive_with_length_norm/"
			shell:
				"echo {config[stoatydive_para]} >> {file_tool_params}"
				"&& if [ ! -d {params.output_folder} ]; then mkdir {params.output_folder}; fi"
				"&& python3 {config[stoatydive]}/bin/StoatyDive.py -a {input.peaks} -b {input.bam} -c {input.genome} "
				"{config[stoatydive_para]} -o {params.output_folder}"

		rule stoatydive2:
			input:
				peaks=PEAKCALLING_OUTDIR + "/robust_peaks_refined.bed",
				bam=bam_folder + "/{sample}_{replicate}" + bam_suffix + ".bam",
				genome=GENOME_SIZES
			output:
				MOTIF_DETECTION_OUTDIR + "/stoatydive_wo_length_norm/VC_Distribution_{sample}_{replicate}" + bam_suffix + ".pdf"
			threads: 2
			conda:
				config["conda_envs"] + "/stoatydive.yml"
			params:
				output_folder=MOTIF_DETECTION_OUTDIR + "/stoatydive_wo_length_norm/"
			shell:
				"echo {config[stoatydive_para2]} >> {file_tool_params}"
				"&& if [ ! -d {params.output_folder} ]; then mkdir {params.output_folder}; fi"
				"&& python3 {config[stoatydive]}/bin/StoatyDive.py -a {input.peaks} -b {input.bam} -c {input.genome} "
				"{config[stoatydive_para2]} -o {params.output_folder}"

else:

	if ( peakcaller == "PureCLIP" or peakcaller == "Piranha" ):
		#####################
		## MOTIF DETECTION ##
		#####################

		rule extract_genomic_DNA_dreme:
			input:
				PEAKCALLING_OUTDIR + "/{sample}_{replicate}_peaks_extended.bed"
			output:
				fa=MOTIF_DETECTION_OUTDIR + "/{sample}_{replicate}_peaks.fa"
			threads: 2
			conda:
				config["conda_envs"] + "/fetch_DNA_sequence.yml"			
			shell:
				"if [ ! -d {MOTIF_DETECTION_OUTDIR} ]; then mkdir {MOTIF_DETECTION_OUTDIR}; fi"
				"&& python " + config["extract_genomic_dna"] + "/fetch_DNA_sequence.py -o {output.fa} {input} {GENOME_FASTA}"

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