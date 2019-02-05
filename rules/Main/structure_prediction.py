
if ( control == "yes" ):

	if ( peakcaller == "PureCLIP" ):

		##########################
		## STRUCTURE PREDICTION ##
		##########################

		rule rna_structure_prediction:
			input: 
				MOTIF_DETECTION_OUTDIR + "/{sample_exp}_{replicate_exp}_{sample_ctl}_{replicate_ctl}_peaks.fa"
			output:
				STRUCTURE_PREDIC_OUTDIR + "/{sample_exp}_{replicate_exp}_{sample_ctl}_{replicate_ctl}_structures.txt"
			conda:
				config["conda_envs"] + "/viennarna.yml"
			shell:
				"if [ ! -d {STRUCTURE_PREDIC_OUTDIR} ]; then mkdir {STRUCTURE_PREDIC_OUTDIR}; fi"
				"&& echo {config[rnafold]} >> {file_tool_params}"
				"&& /home/heylf/ViennaRNA-2.4.11/src/bin/RNAfold {config[rnafold]} < {input} > {output}"

else:

	if ( peakcaller == "PureCLIP" or peakcaller == "Piranha" ):
		
		##########################
		## STRUCTURE PREDICTION ##
		##########################

		rule rna_structure_prediction:
			input: 
				MOTIF_DETECTION_OUTDIR + "/{sample}_{replicate}_peaks.fa"
			output:
				STRUCTURE_PREDIC_OUTDIR + "/{sample}_{replicate}_structures.txt"
			conda:
				config["conda_envs"] + "/viennarna.yml"
			shell:
				"if [ ! -d {STRUCTURE_PREDIC_OUTDIR} ]; then mkdir {STRUCTURE_PREDIC_OUTDIR}; fi"
				"&& echo {config[rnafold]} >> {file_tool_params}"
				"&& /home/heylf/ViennaRNA-2.4.11/src/bin/RNAfold {config[rnafold]} < {input} > {output}"