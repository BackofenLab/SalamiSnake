
#####################
## MAPPING QUALITY ##
#####################

rule finger_print_plot:
	input:
		expand(PRE_FOR_UMI_OUTDIR + "/{all}_got_umis_unlocalized_check.bam", all=ALL_REPLICATES)
	output:
		MAPPING_QUALITY_OUTDIR + "/fingerprint_plot.png"
	threads: 2
	conda:
		config["conda_envs"] + "/deeptools.yml"
	shell:	
		"if [ ! -d {MAPPING_QUALITY_OUTDIR} ]; then mkdir {MAPPING_QUALITY_OUTDIR}; fi"
		"&& echo {config[fingerprint]} >> {file_tool_params}"
		"&& plotFingerprint --numberOfProcessors {threads} --bamfiles {input} {config[fingerprint]} "
		"--labels {ALL_REPLICATES} --plotFile {output}  --plotFileFormat 'png'"

rule correlating_bam_files_plot:
	input:
		expand(PRE_FOR_UMI_OUTDIR + "/{all}_got_umis_unlocalized_check.bam", all=ALL_REPLICATES)
	output:
		bamsummary=MAPPING_QUALITY_OUTDIR + "/correlating_bam_files_summary.txt",
		plot=MAPPING_QUALITY_OUTDIR + "/correlating_bam_files_plot.png"
	threads: 4
	conda:
		config["conda_envs"] + "/deeptools.yml"
	shell:	
		"if [ ! -d {MAPPING_QUALITY_OUTDIR} ]; then mkdir {MAPPING_QUALITY_OUTDIR}; fi"	
		"&& echo {config[multiBamSummary]} >> {file_tool_params}"
		"&& multiBamSummary bins --numberOfProcessors {threads} --outFileName {output.bamsummary} --bamfiles {input} "
		"--labels {ALL_REPLICATES} {config[multiBamSummary]} "
		"&& plotCorrelation {config[plotCorrelation]} --corData {output.bamsummary} --plotFile {output.plot}"

rule fastqc_after_dedup:
    input:
    	PRE_FOR_UMI_OUTDIR + "/{sample}_{replicate}_got_umis_unlocalized_check.bam"
    output:
    	FASTQC_DEDUP_OUTDIR + "/{sample}_{replicate}_got_umis_unlocalized_check_fastqc.html"
    threads: 2
    conda:
    	config["conda_envs"] + "/fastqc.yml"
    shell:
    	"if [ ! -d {FASTQC_DEDUP_OUTDIR} ]; then mkdir {FASTQC_DEDUP_OUTDIR}; fi"
    	"&& fastqc {input} --outdir {FASTQC_DEDUP_OUTDIR}"

#################
## END QUALITY ##
#################

if ( control == "yes" ):
	if ( demultiplexed == "yes" ):
		rule multiqc:
		    input:
		    	begin=expand(FASTQC_BEG_OUTDIR + "/{s}_{r}_{pair}.fastqsanger_fastqc.html", s=SAMPLES[0], r="rep1", pair=PAIR),
		    	trim=expand(FASTQC_ADAPT_OUTDIR + "/{s}_{r}_{pair}_trimmed.fastqsanger_fastqc.html", s=SAMPLES[1], r="rep1", pair=PAIR),
		    	dedup=[expand(FASTQC_DEDUP_OUTDIR + "/{sample}_{replicate}_got_umis_unlocalized_check_fastqc.html", sample=SAMPLES[0], replicate=REP_NAME_CLIP),
		    		   expand(FASTQC_DEDUP_OUTDIR + "/{sample}_{replicate}_got_umis_unlocalized_check_fastqc.html", sample=SAMPLES[1], replicate=REP_NAME_CONTROL)],
		    	mapping=[expand(MAPPING_OUTDIR + "/{sample}_{replicate}_Log.final.out", sample=SAMPLES[0], replicate=REP_NAME_CLIP),
		    			 expand(MAPPING_OUTDIR + "/{sample}_{replicate}_Log.final.out", sample=SAMPLES[1], replicate=REP_NAME_CONTROL)]
		    output:
		    	MULTIQC_OUTDIR + "/multiqc_report.html"
		    threads: 2
		    conda:
		    	config["conda_envs"] + "/multiqc.yml"
		    shell:
		    	"if [ ! -d {MULTIQC_OUTDIR} ]; then mkdir {MULTIQC_OUTDIR}; fi"
		    	"&& multiqc -s {FASTQC_BEG_OUTDIR} {FASTQC_ADAPT_OUTDIR} {FASTQC_DEDUP_OUTDIR} {input.mapping} --outdir {MULTIQC_OUTDIR}"
	else:
		rule multiqc:
		    input:
		    	begin=expand(FASTQC_BEG_OUTDIR + "/{s}_{r}_{pair}.fastqsanger_fastqc.html", s=MULTIPLEX_SAMPLE_NAME, r="rep1", pair=PAIR),
		    	trim=expand(FASTQC_ADAPT_OUTDIR + "/{s}_{r}_{pair}_trimmed.fastqsanger_fastqc.html", s=MULTIPLEX_SAMPLE_NAME, r="rep1", pair=PAIR),
		    	dedup=[expand(FASTQC_DEDUP_OUTDIR + "/{sample}_{replicate}_got_umis_unlocalized_check_fastqc.html", sample=SAMPLES[0], replicate=REP_NAME_CLIP),
		    		   expand(FASTQC_DEDUP_OUTDIR + "/{sample}_{replicate}_got_umis_unlocalized_check_fastqc.html", sample=SAMPLES[1], replicate=REP_NAME_CONTROL)],
		    	mapping=[expand(MAPPING_OUTDIR + "/{sample}_{replicate}_Log.final.out", sample=SAMPLES[0], replicate=REP_NAME_CLIP),
		    			 expand(MAPPING_OUTDIR + "/{sample}_{replicate}_Log.final.out", sample=SAMPLES[1], replicate=REP_NAME_CONTROL)]
		    output:
		    	MULTIQC_OUTDIR + "/multiqc_report.html"
		    threads: 2
		    conda:
		    	config["conda_envs"] + "/multiqc.yml"
		    shell:
		    	"if [ ! -d {MULTIQC_OUTDIR} ]; then mkdir {MULTIQC_OUTDIR}; fi"
		    	"&& multiqc -s {FASTQC_BEG_OUTDIR} {FASTQC_ADAPT_OUTDIR} {FASTQC_DEDUP_OUTDIR} {input.mapping} --outdir {MULTIQC_OUTDIR}"
else:
	if ( demultiplexed == "yes" ):
		rule multiqc:
		    input:
		    	begin=expand(FASTQC_BEG_OUTDIR + "/{sample}_{replicate}_{pair}.fastqsanger_fastqc.html", sample=SAMPLES[0], replicate=REP_NAME_CLIP, pair=PAIR),
		    	trim=expand(FASTQC_ADAPT_OUTDIR + "/{sample}_{replicate}_{pair}_trimmed.fastqsanger_fastqc.html", sample=SAMPLES[0], replicate=REP_NAME_CLIP, pair=PAIR),
		    	dedup=expand(FASTQC_DEDUP_OUTDIR + "/{sample}_{replicate}_got_umis_unlocalized_check_fastqc.html", sample=SAMPLES[0], replicate=REP_NAME_CLIP),
		    	mapping=expand(MAPPING_OUTDIR + "/{sample}_{replicate}_Log.final.out", sample=SAMPLES[0], replicate=REP_NAME_CLIP)
		    output:
		    	MULTIQC_OUTDIR + "/multiqc_report.html"
		    threads: 2
		    conda:
		    	config["conda_envs"] + "/multiqc.yml"
		    shell:
		    	"if [ ! -d {MULTIQC_OUTDIR} ]; then mkdir {MULTIQC_OUTDIR}; fi"
		    	"&& multiqc -s {FASTQC_BEG_OUTDIR} {FASTQC_ADAPT_OUTDIR} {FASTQC_DEDUP_OUTDIR} {input.mapping} --outdir {MULTIQC_OUTDIR}"
	else:
		rule multiqc:
		    input:
		    	begin=expand(FASTQC_BEG_OUTDIR + "/{s}_{r}_{pair}.fastqsanger_fastqc.html", s=MULTIPLEX_SAMPLE_NAME, r="rep1", pair=PAIR),
		    	trim=expand(FASTQC_ADAPT_OUTDIR + "/{s}_{r}_{pair}_trimmed.fastqsanger_fastqc.html", s=MULTIPLEX_SAMPLE_NAME, r="rep1", pair=PAIR),
		    	dedup=expand(FASTQC_DEDUP_OUTDIR + "/{sample}_{replicate}_got_umis_unlocalized_check_fastqc.html", sample=SAMPLES[0], replicate=REP_NAME_CLIP),
		    	mapping=expand(MAPPING_OUTDIR + "/{sample}_{replicate}_Log.final.out", sample=SAMPLES[0], replicate=REP_NAME_CLIP)
		    output:
		    	MULTIQC_OUTDIR + "/multiqc_report.html"
		    threads: 2
		    conda:
		    	config["conda_envs"] + "/multiqc.yml"
		    shell:
		    	"if [ ! -d {MULTIQC_OUTDIR} ]; then mkdir {MULTIQC_OUTDIR}; fi"
		    	"&& multiqc -s {FASTQC_BEG_OUTDIR} {FASTQC_ADAPT_OUTDIR} {FASTQC_DEDUP_OUTDIR} {input.mapping} --outdir {MULTIQC_OUTDIR}"
