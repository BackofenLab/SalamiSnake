
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
		"envs/deeptools.yml"
	shell:	
		"if [ ! -d {MAPPING_QUALITY_OUTDIR} ]; then mkdir {MAPPING_QUALITY_OUTDIR}; fi"
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
		"envs/deeptools.yml"
	shell:	
		"if [ ! -d {MAPPING_QUALITY_OUTDIR} ]; then mkdir {MAPPING_QUALITY_OUTDIR}; fi"	
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
    	"envs/fastqc.yml"
    shell:
    	"if [ ! -d {FASTQC_DEDUP_OUTDIR} ]; then mkdir {FASTQC_DEDUP_OUTDIR}; fi"
    	"&& fastqc {input} --outdir {FASTQC_DEDUP_OUTDIR}"

#################
## END QUALITY ##
#################

rule multiqc:
    input:
    	dedup=expand(FASTQC_DEDUP_OUTDIR + "/{sample}_{replicate}_got_umis_unlocalized_check_fastqc.html", sample=SAMPLES[0], replicate=REP_NAME_CLIP),
    	mapping=expand(MAPPING_OUTDIR + "/{sample}_{replicate}_Log.final.out", sample=SAMPLES[0], replicate=REP_NAME_CLIP)
    output:
    	MULTIQC_OUTDIR + "/multiqc_report.html"
    threads: 2
    conda:
    	"envs/multiqc.yml"
    shell:
    	"if [ ! -d {MULTIQC_OUTDIR} ]; then mkdir {MULTIQC_OUTDIR}; fi"
    	"&& multiqc -s {FASTQC_BEG_OUTDIR} {FASTQC_ADAPT_OUTDIR} {FASTQC_DEDUP_OUTDIR} {input.mapping} --outdir {MULTIQC_OUTDIR}"