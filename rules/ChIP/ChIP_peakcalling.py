
#################
## PEAKCALLING ##
#################

if ( control == "yes" ):

	if ( peakcaller == "PEAKachu" ):

		rule peakachu:
		    input:
		    	clip=expand(PRE_FOR_UMI_OUTDIR + "/{sample}_got_umis_unlocalized_check.bam", sample=REPLICATES_CLIP),
		    	control=expand(PRE_FOR_UMI_OUTDIR + "/{sample}_got_umis_unlocalized_check.bam", sample=REPLICATES_CONTROL)
		    output:
		    	peaks_tsv=PEAKCALLING_OUTDIR + "/peakachu.tsv",
		    	peaks_gtf=PEAKCALLING_OUTDIR + "/peakachu.gtf",
		    	blockbuster=PEAKCALLING_OUTDIR + "/blockbuster.bed"
		    threads: 4
		    conda:
		    	config["conda_envs"] + "/peakachu.yml"
		    shell:
		    	#"--max_proc "${GALAXY_SLOTS:-1}"
		    	"if [ ! -d {PEAKCALLING_OUTDIR} ]; then mkdir {PEAKCALLING_OUTDIR}; fi"
		    	"&& echo {config[peakachu]} >> {file_tool_params}"
		    	#"&& source activate peakachu"
		    	"&& peakachu adaptive --exp_libs {input.clip} --ctr_libs {input.control} --paired_end --output_folder {PEAKCALLING_OUTDIR} {config[peakachu]}"
		    	#"&& source deactivate "
		    	"&& head -q -n 1 {PEAKCALLING_OUTDIR}/peak_tables/*.csv > tmp.tsv "
		    	"&& head -q -n 1 tmp.tsv > {output.peaks_tsv}"
		    	"&& rm tmp.tsv"
		    	"&& tail -n +2 -q {PEAKCALLING_OUTDIR}/peak_tables/*.csv >> {output.peaks_tsv} "
		    	"&& cat {PEAKCALLING_OUTDIR}/peak_annotations/*.gff | awk '/peak/ {{ print $0 }}' > {output.peaks_gtf} "
		    	"&& cat {PEAKCALLING_OUTDIR}/blockbuster_input/*bed > {output.blockbuster}"

		rule peaks_tsv_to_bed:
			input:
				PEAKCALLING_OUTDIR + "/peakachu.tsv"
			output:
				PEAKCALLING_OUTDIR + "/peakachu.bed"
			threads: 2
			shell:
				"""&&  awk -F "\t" 'BEGIN {{ OFS = FS }} NR>1 {{ if ($3 < $4) {{ print $1,$3,$4,"clip_peak_"NR-1,$9,$5; }} else {{ print $1,$4,$3,"clip_peak_"NR-1,$9,$5; }} }}' {input} > {output} """

		rule peaks_extend_frontiers:
			input:
				bed=PEAKCALLING_OUTDIR + "/peakachu.bed",
				genome=GENOME_SIZES
			output:
				PEAKCALLING_OUTDIR + "/peakachu_peaks_extended.bed"
			threads: 2
			conda:
				config["conda_envs"] + "/bedtools.yml"
			shell:
				"echo {config[peaks_extend_frontiers]} >> {file_tool_params}"
				"&& bedtools slop {config[peaks_extend_frontiers]} -i {input.bed} -g {input.genome} > {output}"

	if ( peakcaller == "MACS2" ):
		rule macs2:
			input:
				experiment=PRE_FOR_UMI_OUTDIR + "/{sample_exp}_{replicate_exp}_got_umis_unlocalized_check.bam",
				control=PRE_FOR_UMI_OUTDIR + "/{sample_ctl}_{replicate_ctl}_got_umis_unlocalized_check.bam",
				genome_fasta=GENOME_FASTA
			output:
				# If you want to find the motifs at the binding sites, this file is recommended (from https://github.com/taoliu/MACS)
				summits=PEAKCALLING_OUTDIR + "/{sample_exp}_{replicate_exp}_{sample_ctl}_{replicate_ctl}_summits.bed",
				# NAME_peaks.narrowPeak is BED6+4 format file which contains the peak locations together with peak summit, pvalue and qvalue. 
				peaks=PEAKCALLING_OUTDIR + "/{sample_exp}_{replicate_exp}_{sample_ctl}_{replicate_ctl}_binding_regions.bed"
			threads: 10
			params:
				output_folder=PEAKCALLING_OUTDIR,
				name="{sample_exp}_{replicate_exp}_{sample_ctl}_{replicate_ctl}",
				peaks_tmp=PEAKCALLING_OUTDIR + "/{sample_exp}_{replicate_exp}_{sample_ctl}_{replicate_ctl}_peaks.narrowPeak"
			conda:
				config["conda_envs"] + "/macs2.yml"
			shell:
				"if [ ! -d {PEAKCALLING_OUTDIR} ]; then mkdir {PEAKCALLING_OUTDIR}; fi "
				"&& echo {config[macs2]} >> {file_tool_params}"
				"&& macs2 callpeak -t {input.experiment} -c {input.control} -n {params.name} -f BAMPE --outdir {params.output_folder} {config[macs2]}"
				"&& mv {params.peaks_tmp} {output.peaks}"

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

		rule find_robust_peaks:
			input:
				expand(PEAKCALLING_OUTDIR + "/{sample_exp}_{replicate_exp}_{sample_ctl}_{replicate_ctl}_peaks_extended.bed", 
					sample_exp=SAMPLES[0], replicate_exp=REP_NAME_CLIP, sample_ctl=SAMPLES[1], replicate_ctl=REP_NAME_CONTROL)
			output:
				ROBUSTPEAKS_OUTDIR + "/robust_between_all.bed"
			threads: 2
			conda:
				config["conda_envs"] + "/bedtools.yml"
			params:
				input_folder=PEAKCALLING_OUTDIR,
				output_folder=ROBUSTPEAKS_OUTDIR
			shell:
				"if [ ! -d {ROBUSTPEAKS_OUTDIR} ]; then mkdir {ROBUSTPEAKS_OUTDIR}; fi"
				"&& {config[find_robust_intersections]}/robust_intersections.sh {params.input_folder} {params.output_folder} bed"

else:

	if ( peakcaller == "Piranha" ):
		rule piranha:
			input:
				PRE_FOR_UMI_OUTDIR + "/{sample}_{replicate}_got_umis_unlocalized_check.bam"
			output:
				PEAKCALLING_OUTDIR + "/{sample}_{replicate}_peaks.bed"
			conda:
				config["conda_envs"] + "/piranha.yml"
			threads: 2 
			shell:
				"if [ ! -d {PEAKCALLING_OUTDIR} ]; then mkdir {PEAKCALLING_OUTDIR}; fi"
				"&& echo {config[piranha]} >> {file_tool_params}"
				"&& Piranha {config[piranha]} {input} -o {output}"

		rule peaks_extend_frontiers:
			input:
				bed=PEAKCALLING_OUTDIR + "/{sample}_{replicate}_peaks.bed",
				genome=GENOME_SIZES
			output:
				PEAKCALLING_OUTDIR + "/{sample}_{replicate}_peaks_extended.bed"
			threads: 2
			conda:
				config["conda_envs"] + "/bedtools.yml"
			shell:
				"echo {config[peaks_extend_frontiers]} >> {file_tool_params}"
				"&& bedtools slop {config[peaks_extend_frontiers]} -i {input.bed} -g {input.genome} > {output}"

	elif ( peakcaller == "PureCLIP" ):
		rule pureclip:
			input:
				experiment=PRE_FOR_UMI_OUTDIR + "/{sample}_{replicate}_got_umis_unlocalized_check.bam",
				experiment_bai=PRE_FOR_UMI_OUTDIR + "/{sample}_{replicate}_got_umis_unlocalized_check.bam.bai",
				genome_fasta=GENOME_FASTA
			output:
				crosslinking_sites=PEAKCALLING_OUTDIR + "/{sample}_{replicate}_crosslinkind_sites.bed",
				binding_regions=PEAKCALLING_OUTDIR + "/{sample}_{replicate}_binding_regions.bed"
			threads: 10
			params:
				tmp=PEAKCALLING_OUTDIR + "/tmp/",
				parameters=PEAKCALLING_OUTDIR + "/{sample}_{replicate}_parameters.txt"
			conda:
				config["conda_envs"] + "/pureclip.yml"
			shell:
				"if [ ! -d {PEAKCALLING_OUTDIR} ]; then mkdir {PEAKCALLING_OUTDIR}; fi "
				"&& echo {config[pureclip]} >> {file_tool_params}"
				"&& pureclip -i {input.experiment} -bai {input.experiment_bai} -g {input.genome_fasta} -o {output.crosslinking_sites} -tmp {params.tmp} "
				"-or {output.binding_regions} -p {params.parameters} -nt {threads} {config[pureclip]} "

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

	rule find_robust_peaks:
		input:
			expand(PEAKCALLING_OUTDIR + "/{sample}_{replicate}_peaks_extended.bed", sample=SAMPLES[0], replicate=REP_NAME_CLIP)
		output:
			ROBUSTPEAKS_OUTDIR + "/robust_between_all.bed"
		threads: 2
		conda:
			config["conda_envs"] + "/bedtools.yml"
		params:
			input_folder=PEAKCALLING_OUTDIR,
			output_folder=ROBUSTPEAKS_OUTDIR
		shell:
			"if [ ! -d {ROBUSTPEAKS_OUTDIR} ]; then mkdir {ROBUSTPEAKS_OUTDIR}; fi"
			"&& {config[find_robust_intersections]}/robust_intersections.sh {params.input_folder} {params.output_folder} bed"
