import random
import math
import os 

SRC_PATH = os.getcwd()

configfile: "config.yml"

REF_GENOME_DIR= config['ref_genome_dir']
GENOME_FASTA = REF_GENOME_DIR + config['genome_fasta']
GENOME_2BIT = REF_GENOME_DIR + config['genome_2bit']
GENOME_GTF = REF_GENOME_DIR + config['genome_gtf']
GENOME_SIZES = REF_GENOME_DIR + config['genome_sizes']

###############
## FUNCTIONS ##
###############

def list_all_values_of_dict(dictionary):
	for values in dictionary.values():
	    if isinstance(values, dict):
	      yield from list_all_values_of_dict(values)
	    else:
	      yield values

def find_all_values_on_key(key, dictionary):
	for k, v in dictionary.items():
		if k == key:
			yield v
		if isinstance(v, dict):
			yield from find_all_values_on_key(key, v)

##############################
## WRITE TOOL PARAMS TO LOG ##
##############################

file_tool_params = open(config['sample_data_dir'] + "/tool_params.txt", "w")
for i in config['all_tool_params']:
	file_tool_params.write(i + "\n")
file_tool_params.close()

######################
## SAMPLE VARIABLES ##
######################

FIRST_READS = list(find_all_values_on_key("r1", config['clip_samples'])) + list(find_all_values_on_key("r1", config['control_samples']))
SECOND_READS = list(find_all_values_on_key("r2", config['clip_samples'])) + list(find_all_values_on_key("r2", config['control_samples']))
CLIP_SAMPLES = list(list_all_values_of_dict(config['clip_samples']))
CONTROL_SAMPLES = list(list_all_values_of_dict(config['control_samples']))
ALL_SAMPLES = CLIP_SAMPLES + CONTROL_SAMPLES

REPLICATES_CLIP = list(config['clip_samples'].keys())
REPLICATES_CONTROL = list(config['control_samples'].keys())
ALL_REPLICATES = REPLICATES_CLIP + REPLICATES_CONTROL

print("[NOTE] Loaded Samples.")
print("All samples: " + str(ALL_SAMPLES))
print("CLIP samples: " + str(CLIP_SAMPLES))
print("Conrol samples: " + str(CONTROL_SAMPLES))
print("First read samples: " + str(FIRST_READS))
print("Second read samples: " + str(SECOND_READS))

PAIR = ["r1", "r2"]
REP_NAME_CLIP = ["rep" + str(x+1) for x in range(0,len(REPLICATES_CLIP))]
REP_NAME_CONTROL = ["rep" + str(x+1) for x in range(0,len(REPLICATES_CONTROL))]
SAMPLES = [REPLICATES_CLIP[0].split("_")[0], REPLICATES_CONTROL[0].split("_")[0]]

print(SAMPLES)

######################
## PATH DEFINITIONS ##
######################

RENAMING = config['sample_data_dir'] + "/" + config['renaming_outdir']
FASTQC_BEG_OUTDIR = config['sample_data_dir'] + "/" + config['fastqc_beg_outdir']
FASTQC_ADAPT_OUTDIR = config['sample_data_dir'] + "/" + config['fastqc_adapt_outdir']
CUTADAPT_OUTDIR = config['sample_data_dir'] + "/" + config['cutadapt_outdir']
REMOVE_TAIL_OUTDIR = config['sample_data_dir'] + "/" + config['remove_tail_outdir']
TRIMGALORE_OUTDIR = config['sample_data_dir'] + "/" + config['trim_galore_outdir']
MAPPING_OUTDIR = config['sample_data_dir'] + "/" + config['mapping_outdir']
MAPPING_QUALITY_OUTDIR = config['sample_data_dir'] + "/" + config['mapping_quality_outdir']
PRE_FOR_UMI_OUTDIR = config['sample_data_dir'] + "/" + config['preprocessing_for_umi_outdir']
DEDUPLICAITON_OUTDIR = config['sample_data_dir'] + "/" + config['deduplication_outdir']
FASTQC_DEDUP_OUTDIR = config['sample_data_dir'] + "/" + config['fastqc_dedup_outdir']
COVERAGE_OUTDIR = config['sample_data_dir'] + "/" + config['coverage_outdir']
PEAKCALLING_OUTDIR = config['sample_data_dir'] + "/" + config['peakcalling_outdir']
POSTPROCESSING_OUTDIR = config['sample_data_dir'] + "/" + config['post_processing_outdir']
MOTIF_DETECTION_OUTDIR = config['sample_data_dir'] + "/" + config['motif_detection_outdir']
MOTIF_SEARCH_OUTDIR = config['sample_data_dir'] + "/" + config['motif_search_outdir']

###################
## PREPROCESSING ##
###################

rule all:
	input: 
		expand(CUTADAPT_OUTDIR + "/{sample}_{replicate}_r1.fastqsanger", sample=SAMPLES[0], replicate=REP_NAME_CLIP),
		expand(CUTADAPT_OUTDIR + "/{sample}_{replicate}_r1.fastqsanger", sample=SAMPLES[1], replicate=REP_NAME_CONTROL),
		expand(FASTQC_FIRST_OUTDIR + "/{sample}_{replicate}_{pair}.fastqsanger_fastqc.html", sample=SAMPLES[0], replicate=REP_NAME_CLIP, pair=PAIR),
		expand(FASTQC_FIRST_OUTDIR + "/{sample}_{replicate}_{pair}.fastqsanger_fastqc.html", sample=SAMPLES[1], replicate=REP_NAME_CONTROL, pair=PAIR),
		REF_GENOME_DIR + "/sjdbList.fromGTF.out.tab",
		expand(MAPPING_OUTDIR + "/{sample}_{replicate}.bam", sample=SAMPLES[0], replicate=REP_NAME_CLIP),
		expand(MAPPING_OUTDIR + "/{sample}_{replicate}.bam", sample=SAMPLES[1], replicate=REP_NAME_CONTROL),
		expand(MAPPING_OUTDIR + "/{sample}_{replicate}.bam.bai", sample=SAMPLES[0], replicate=REP_NAME_CLIP),
		expand(MAPPING_OUTDIR + "/{sample}_{replicate}.bam.bai", sample=SAMPLES[1], replicate=REP_NAME_CONTROL),
		expand(MAPPING_QUALITY_OUTDIR + "/{sample}_{replicate}_gc_bias_plot.png", sample=SAMPLES[0], replicate=REP_NAME_CLIP),
		expand(MAPPING_QUALITY_OUTDIR + "/{sample}_{replicate}_gc_bias_plot.png", sample=SAMPLES[1], replicate=REP_NAME_CONTROL),
		expand(MAPPING_QUALITY_OUTDIR + "/{sample}_{replicate}_insert_size_plot.pdf", sample=SAMPLES[0], replicate=REP_NAME_CLIP),
		expand(MAPPING_QUALITY_OUTDIR + "/{sample}_{replicate}_insert_size_plot.pdf", sample=SAMPLES[1], replicate=REP_NAME_CONTROL),
		MAPPING_QUALITY_OUTDIR + "/fingerprint_plot.png",
		#MAPPING_QUALITY_OUTDIR + "/correlating_bam_files_plot.png",
		expand(PRE_FOR_UMI_OUTDIR + "/{sample}_{replicate}_got_umis.bam", sample=SAMPLES[0], replicate=REP_NAME_CLIP),
		expand(PRE_FOR_UMI_OUTDIR + "/{sample}_{replicate}_got_umis.bam", sample=SAMPLES[1], replicate=REP_NAME_CONTROL),
		expand(FASTQC_SECOND_OUTDIR + "/{sample}_{replicate}_fastqc.html", sample=SAMPLES[0], replicate=REP_NAME_CLIP),
		expand(FASTQC_SECOND_OUTDIR + "/{sample}_{replicate}_fastqc.html", sample=SAMPLES[1], replicate=REP_NAME_CONTROL),
		PEAKCALLING_OUTDIR + "/peakachu.tsv",
		expand(PEAKCALLING_OUTDIR + "/{sample}_{replicate}_htseq_hits.txt", sample=SAMPLES[0], replicate=REP_NAME_CLIP),
		expand(PEAKCALLING_OUTDIR + "/{sample}_{replicate}_htseq_hits.txt", sample=SAMPLES[1], replicate=REP_NAME_CONTROL),
		POSTPROCESSING_OUTDIR + "/peakachu_extended.bed",
		MOTIF_DETECTION_OUTDIR + "/peaks.fa",
		expand(COVERAGE_OUTDIR + "/{sample}_{replicate}_crosslink_sites_intersecting_peaks.bed", sample=SAMPLES[0], replicate=REP_NAME_CLIP),
		expand(COVERAGE_OUTDIR + "/{sample}_crosslink_sites_quality.txt", sample=REPLICATES_CLIP),
		expand(COVERAGE_OUTDIR + "/{sample}_{replicate}_crosslinking_positions_sorted.bed", sample=SAMPLES[0], replicate=REP_NAME_CLIP),
		expand(COVERAGE_OUTDIR + "/bedgraph/{sample}_{replicate}_crosslinking_coverage_pos_strand.bedgraph", sample=SAMPLES[0], replicate=REP_NAME_CLIP),
		expand(COVERAGE_OUTDIR + "/bigwig/{sample}_{replicate}_crosslinking_coverage_{type}_strand.bigwig", sample=SAMPLES[0], replicate=REP_NAME_CLIP, type=["pos", "neg", "both"]),
		expand(COVERAGE_OUTDIR + "/bigwig/{sample}_{replicate}_crosslinking_coverage_{type}_strand.bigwig", sample=SAMPLES[1], replicate=REP_NAME_CONTROL, type=["pos", "neg", "both"]),
		MOTIF_DETECTION_OUTDIR + "/meme_chip/meme-chip.html"


ALL_NEW_FILE_NAMES = ["none"] * ((len(REPLICATES_CLIP) + len(REPLICATES_CONTROL)) * 2)
i = 0 
for j in REP_NAME_CLIP:
	for k in PAIR:
		ALL_NEW_FILE_NAMES[i] = RENAMING + "/" + SAMPLES[0] + "_" + j + "_" + k + ".fastqsanger"
		i += 1

for j in REP_NAME_CONTROL:
	for k in PAIR:
		ALL_NEW_FILE_NAMES[i] = RENAMING + "/" + SAMPLES[1] + "_" + j + "_" + k + ".fastqsanger"
		i += 1

rule renaming:
	input:
		expand(config['sample_data_dir'] + "/{all}.fastqsanger", all=ALL_SAMPLES)
	output:
		file=ALL_NEW_FILE_NAMES,
		log=RENAMING + "/renaming_log.txt"
	run:
		shell("if [ ! -d {RENAMING} ]; then mkdir {RENAMING}; fi "
		"&& echo 'original_name \t new_name' > {output.log}")
		for i in range(0, len(input)): 
			shell("cp " + input[i] + " " + output.file[i] + 
			"&& echo " + input[i] + " \t " + output.file[i] + " >> {output.log}")

rule fastqc_beginning:
    input:
    	RENAMING + "/{sample}_{replicate}_{pair}.fastqsanger"
    output:
    	FASTQC_BEG_OUTDIR + "/{sample}_{replicate}_{pair}.fastqsanger_fastqc.html",
    	FASTQC_BEG_OUTDIR + "/{sample}_{replicate}_{pair}.fastqsanger_fastqc.zip"
    threads: 2
    conda:
    	"envs/fastqc.yml"
    shell:
    	"if [ ! -d {FASTQC_BEG_OUTDIR} ]; then mkdir {FASTQC_BEG_OUTDIR}; fi"
    	"&& fastqc {input} --outdir {FASTQC_BEG_OUTDIR}"

# remove adapter from first and second read
rule cutadapt_adapter_trimming:
	input:
		first=RENAMING + "/{sample}_{replicate}_r1.fastqsanger",
		second=RENAMING + "/{sample}_{replicate}_r2.fastqsanger"
	output:
		seq_first=CUTADAPT_OUTDIR + "/{sample}_{replicate}_r1.fastqsanger",
		seq_second=CUTADAPT_OUTDIR + "/{sample}_{replicate}_r2.fastqsanger",
		log=CUTADAPT_OUTDIR + "/{sample}_{replicate}.txt"
	threads: 8 
	conda:
		"envs/cutadapt.yml"
	shell:
		"if [ ! -d {CUTADAPT_OUTDIR} ]; then mkdir {CUTADAPT_OUTDIR}; fi"
		"&& cutadapt -j {threads} {config[cutadapt]} "
		"--paired-output={output.seq_second} --output={output.seq_first} {input.first} {input.second} > {output.log}"

# remove barcode length from second read
rule remove_tail:
	input:
		first=CUTADAPT_OUTDIR + "/{samples}_{replicate}_r1.fastqsanger",
		second=CUTADAPT_OUTDIR + "/{samples}_{replicate}_r2.fastqsanger"
	output:
		first=REMOVE_TAIL_OUTDIR + "/{samples}_{replicate}_r1.fastqsanger",
		second=REMOVE_TAIL_OUTDIR + "/{samples}_{replicate}_r2.fastqsanger"
	threads: 2
	conda:
		"envs/bctools.yml"
	shell:
		"if [ ! -d {REMOVE_TAIL_OUTDIR} ]; then mkdir {REMOVE_TAIL_OUTDIR}; fi"
		"&& python {config[bctools]}/remove_tail.py {input.second} 9 > {output.second}"
		"&& cp {input.first} {output.first}"

# old version used to cut off bases of potential UMI length (here N=5 of eCLIP data)
# rule trim_galore_clip:
# 	input:
# 		first_read=REMOVE_TAIL_OUTDIR + "/{sample}_{replicate}_r1.fastqsanger",
# 		second_read=REMOVE_TAIL_OUTDIR + "/{sample}_{replicate}_r2.fastqsanger"
# 	output:
# 		first_read_out=TRIMGALORE_OUTDIR + "/{sample}_{replicate}_r1.fastqsanger",
# 		second_read_out=TRIMGALORE_OUTDIR + "/{sample}_{replicate}_r2.fastqsanger"
# 	threads: 2
# 	params:
# 		cut_n_bases=5
# 	conda:
# 		"envs/trim_galore.yml"
# 	shell:
# 		"if [ ! -d {TRIMGALORE_OUTDIR} ]; then mkdir {TRIMGALORE_OUTDIR}; fi"
# 		"&& trim_galore --phred33 --suppress_warn --three_prime_clip_R1 {params.cut_n_bases} --paired --dont_gzip "
# 		"--output_dir {TRIMGALORE_OUTDIR} {input.first_read} {input.second_read}"
# 		"&& mv {output.first_read_out}_val_1.fq {output.first_read_out}" 
# 		"&& mv {output.second_read_out}_val_2.fq {output.second_read_out}"

# Activate if you need to remove UMIs and barcodes

# rule got_umis:
# 	input:
# 		first=REMOVE_TAIL_OUTDIR + "/{samples}_{replicate}_r1.fastqsanger",
# 		second=REMOVE_TAIL_OUTDIR + "/{samples}_{replicate}_r2.fastqsanger"
# 	output:
# 		first=PRE_FOR_UMI_OUTDIR + "/{samples}_{replicate}_r1.fastqsanger",
# 		second=PRE_FOR_UMI_OUTDIR + "/{samples}_{replicate}_r2.fastqsanger",
# 		log=PRE_FOR_UMI_OUTDIR + "/{samples}_{replicate}_log.txt"
# 	params:
# 		first=PRE_FOR_UMI_OUTDIR + "/{samples}_{replicate}_buff_r1.fastqsanger",
# 		second=PRE_FOR_UMI_OUTDIR + "/{samples}_{replicate}_buff_r2.fastqsanger",
# 		log=PRE_FOR_UMI_OUTDIR + "/{samples}_{replicate}_buff__log.txt"
# 	threads: 2
# 	conda:
# 		"envs/umi.yml"
# 	shell:
# 		"if [ ! -d {PRE_FOR_UMI_OUTDIR} ]; then mkdir {PRE_FOR_UMI_OUTDIR}; fi "
# 		"&& umi_tools extract {config[got_umis_1]} -I {input.second} -S {params.second} --read2-in {input.first} --read2-out {params.first} -L {params.log}"
# 		"&& umi_tools extract {config[got_umis_2]} -I {params.second} -S {output.second} --read2-in {params.first} --read2-out {output.first} -L {output.log}"

rule fastqc_after_adapter_removal:
    input:
    	REMOVE_TAIL_OUTDIR + "/{sample}_{replicate}_{pair}.fastqsanger"
    output:
    	FASTQC_ADAPT_OUTDIR + "/{sample}_{replicate}_{pair}.fastqsanger_fastqc.html",
    	FASTQC_ADAPT_OUTDIR + "/{sample}_{replicate}_{pair}.fastqsanger_fastqc.zip"
    threads: 2
    conda:
    	"envs/fastqc.yml"
    shell:
    	"if [ ! -d {FASTQC_ADAPT_OUTDIR} ]; then mkdir {FASTQC_ADAPT_OUTDIR}; fi"
    	"&& fastqc {input} --outdir {FASTQC_ADAPT_OUTDIR}"

#############
## MAPPING ##
#############

rule star_generate_index_for_genome:
	input:
		fasta=GENOME_FASTA,
		annotation=GENOME_GTF
	output:
		REF_GENOME_DIR + "/sjdbList.fromGTF.out.tab"
	threads: 4
	conda:
		"envs/star.yml"
	shell:
		"if [ -d {config[sample_data_dir]}/STAR_tmp_Index ]; then rm -r {config[sample_data_dir]}/STAR_tmp_Index; fi "
		"&& STAR --outTmpDir {config[sample_data_dir]}/STAR_tmp_Index --runThreadN {threads} --runMode genomeGenerate --genomeDir {REF_GENOME_DIR} "
		"--genomeFastaFiles {input.fasta} --sjdbGTFfile {input.annotation}"

rule star:
	input:
		first_read=REMOVE_TAIL_OUTDIR + "/{sample}_{replicate}_r1.fastqsanger",
		second_read=REMOVE_TAIL_OUTDIR + "/{sample}_{replicate}_r2.fastqsanger",
	output:
		log=MAPPING_OUTDIR + "/{sample}_{replicate}.txt",
		bam=MAPPING_OUTDIR + "/{sample}_{replicate}.bam"
	threads: 4
	conda:
		"envs/star.yml"
	params:
		output_folder=MAPPING_OUTDIR + "/{sample}_{replicate}"	
	shell:
		"if [ ! -d {MAPPING_OUTDIR} ]; then mkdir {MAPPING_OUTDIR}; fi "
		"&& TIME=$(date +%N) "
		"&& if [ -d {config[sample_data_dir]}/STAR_tmp_$TIME ]; then rm -r {config[sample_data_dir]}/STAR_tmp_$TIME; fi "
		"&& STAR --runThreadN {threads} --genomeLoad NoSharedMemory --genomeDir {REF_GENOME_DIR} "   
		"--readFilesIn {input.first_read} {input.second_read} --outTmpDir {config[sample_data_dir]}/STAR_tmp_$TIME  --outFileNamePrefix {params.output_folder}_ "  
		"{config[star]} > {output.log} "
		"&& mv {params.output_folder}_Aligned.sortedByCoord.out.bam {output.bam} "
		"&& rm -r {config[sample_data_dir]}/STAR_tmp_$TIME"

rule indexing:
	input:
		MAPPING_OUTDIR + "/{sample}_{replicate}.bam"
	output:
		MAPPING_OUTDIR + "/{sample}_{replicate}.bam.bai"
	threads: 2
	conda:
		"envs/samtools.yml"
	shell:
		"samtools index {input}"

########################
## POST-MAP FILTERING ##
########################

# -F 0x100 do not include secondary alignments 
rule unique_reads_fitlering:
	input:
		MAPPING_OUTDIR + "/{sample}_{replicate}.bam"
	output:
		PRE_FOR_UMI_OUTDIR + "/{sample}_{replicate}_unique_reads_fitlering.bam"
	threads: 2
	conda:
		"envs/samtools.yml"
	shell:
		"if [ ! -d {PRE_FOR_UMI_OUTDIR} ]; then mkdir {PRE_FOR_UMI_OUTDIR}; fi"
		"&& samtools view -h -F 0x100 {input} "
		"| awk '$0 ~ /^@/{{ print }} ($2 == 163 || $2 == 147 || $2 == 83 || $2 == 99)&&$0!~/XS:i/{{ print }} ' "  
		"| samtools view -bSh > {output}"

rule got_umis:
	input:
		PRE_FOR_UMI_OUTDIR + "/{sample}_{replicate}_unique_reads_fitlering.bam"
	output:
		PRE_FOR_UMI_OUTDIR + "/{sample}_{replicate}_got_umis.bam"
	threads: 2
	conda:
		"envs/samtools.yml"
	shell:
		"if [ ! -d {PRE_FOR_UMI_OUTDIR} ]; then mkdir {PRE_FOR_UMI_OUTDIR}; fi"
		"&& samtools view -h {input} "
		"""| awk -F "\t" 'BEGIN {{ OFS = FS }} {{ if ($0 ~ /^@/) {{ print $0; }} else {{split($1,str,":"); $1=str[2]":"str[3]":"str[4]":"str[5]":"str[6]":"str[7]":"str[8]"_"str[1]; print $0; }} }}' """
		"| samtools view -bSh > {output}"

rule filter_out_unlocalized_regions_for_later_genome_versions:
	input:
		PRE_FOR_UMI_OUTDIR + "/{sample}_{replicate}_got_umis.bam"
	output:
		PRE_FOR_UMI_OUTDIR + "/{sample}_{replicate}_got_umis_unlocalized_check.bam"
	threads: 2
	conda:
		"envs/samtools.yml"
	shell:
		"samtools view -h {input} "
		"""| awk -F "\t" 'BEGIN {{ OFS = FS }} {{ if ($0 ~ /^@/) {{print $0;}} else {{ if ($3 ~/^chr/) {{print $0;}} }} }}' """
		"| samtools view -bSh > {output}"
		"&& samtools index {output}"
 
###################
## DEDUPLICATION ##
###################

rule deduplication:
	input:
		PRE_FOR_UMI_OUTDIR + "/{sample}_{replicate}_got_umis_unlocalized_check.bam"
	output:
		bam=DEDUPLICAITON_OUTDIR + "/{sample}_{replicate}.bam",
		log=DEDUPLICAITON_OUTDIR + "/{sample}_{replicate}_log.txt",
		sorted_bam=DEDUPLICAITON_OUTDIR + "/{sample}_{replicate}_sorted.bam"
	threads: 2
	conda:
		"envs/umi.yml"
	shell:
		"if [ ! -d {DEDUPLICAITON_OUTDIR} ]; then mkdir {DEDUPLICAITON_OUTDIR}; fi "
		"&& umi_tools dedup {config[umi_dedup]} -I {input} -S {output.bam} -L {output.log}"
		"&& samtools sort {output.bam} > {output.sorted_bam}"
		"&& samtools index {output.sorted_bam}"

#####################
## MAPPING QUALITY ##
#####################

rule compute_gc_bias_plots:
	input:
		DEDUPLICAITON_OUTDIR + "/{sample}_{replicate}.bam"
	output:
		plot=MAPPING_QUALITY_OUTDIR + "/{sample}_{replicate}_gc_bias_plot.png",
		file=MAPPING_QUALITY_OUTDIR + "/{sample}_{replicate}_gc_bias_data.txt"
	threads: 2
	conda:
		"envs/deeptools.yml"
	shell:
		"if [ ! -d {MAPPING_QUALITY_OUTDIR} ]; then mkdir {MAPPING_QUALITY_OUTDIR}; fi"
		"&& computeGCBias --numberOfProcessors {threads} --bamfile {input} --GCbiasFrequenciesFile {output.file} {config[gc_bias]} "
		"--genome {GENOME_2BIT} --biasPlot {output.plot} --plotFileFormat png"

# for picard please conda install R
rule estimate_insert_size:
	input:
		DEDUPLICAITON_OUTDIR + "/{sample}_{replicate}_got_umis.bam"
	output:
		plot=MAPPING_QUALITY_OUTDIR + "/{sample}_{replicate}_insert_size_plot.pdf",
		report=MAPPING_QUALITY_OUTDIR + "/{sample}_{replicate}_insert_size_report.txt"
	threads: 2
	conda:
		"envs/picard.yml"
	shell:
		"if [ ! -d {MAPPING_QUALITY_OUTDIR} ]; then mkdir {MAPPING_QUALITY_OUTDIR}; fi"
		"&& source activate picard"
		"&& picard CollectInsertSizeMetrics INPUT={input} OUTPUT={output.report} "
		"HISTOGRAM_FILE={output.plot} {config[estimate_insert_size]} "
		"&& source deactivate"

rule finger_print_plot:
	input:
		expand(DEDUPLICAITON_OUTDIR + "/{sample}.bam", sample=ALL_REPLICATES)
	output:
		MAPPING_QUALITY_OUTDIR + "/fingerprint_plot.png"
	threads: 2
	conda:
		"envs/deeptools.yml"
	shell:	
		"if [ ! -d {MAPPING_QUALITY_OUTDIR} ]; then mkdir {MAPPING_QUALITY_OUTDIR}; fi"
		"&& plotFingerprint --numberOfProcessors {threads} --bamfiles {input} {config[fingerprint]} "
		"--labels {ALL_REPLICATES} --plotFile {output}  --plotFileFormat 'png'"

# rule correlating_bam_files_plot:
# 	input:
# 		expand(DEDUPLICAITON_OUTDIR + "/{sample}.bam", sample=ALL_REPLICATES)
# 	output:
# 		bamsummary=MAPPING_QUALITY_OUTDIR + "/correlating_bam_files_summary.txt",
# 		plot=MAPPING_QUALITY_OUTDIR + "/correlating_bam_files_plot.png"
# 	threads: 10
#	params:
#		
# 	conda:
# 		"envs/deeptools.yml"
# 	shell:	
# 		"if [ ! -d {MAPPING_QUALITY_OUTDIR} ]; then mkdir {MAPPING_QUALITY_OUTDIR}; fi"	
# 		"&& multiBamSummary bins --numberOfProcessors {threads} --outFileName {output.bamsummary} --bamfiles {input} "
# 		"--labels {ALL_REPLICATES} --binSize '10000' --distanceBetweenBins '0'"
# 		"&& plotCorrelation --corData {output.bamsummary} --plotFile {output.plot} --corMethod 'spearman' --whatToPlot 'heatmap' "
# 		"--colorMap 'RdYlBu'  --plotTitle ''  --plotWidth 11.0 --plotHeight 9.5  --plotFileFormat 'png'"

############################
## SECOND QUALITY CONTROL ##
############################

rule fastqc_after_dedup:
    input:
    	DEDUPLICAITON_OUTDIR + "/{sample}_{replicate}.bam"
    output:
    	FASTQC_DEDUP_OUTDIR + "/{sample}_{replicate}_fastqc.html"
    threads: 2
    conda:
    	"envs/fastqc.yml"
    shell:
    	"if [ ! -d {FASTQC_DEDUP_OUTDIR} ]; then mkdir {FASTQC_DEDUP_OUTDIR}; fi"
    	"&& fastqc {input} --outdir {FASTQC_DEDUP_OUTDIR}"

#################
## PEAKCALLING ##
#################

# BE CAREFUL when using other protocols then you might need to define a different bitflag!!
# rule mate_reads_fitlering:
# 	input:
# 		DEDUPLICAITON_OUTDIR + "/{sample}_{replicate}_sorted.bam"
# 	output:
# 		bam=MATEFILTER_OUTDIR + "/{sample}_{replicate}.bam",
# 		sorted_bam=MATEFILTER_OUTDIR + "/{sample}_{replicate}_sorted.bam",
# 		sorted_bam_bai=MATEFILTER_OUTDIR + "/{sample}_{replicate}_sorted.bam.bai"
# 	threads: 2
#	conda:
#		"envs/samtools.yml"
# 	shell:
# 		"if [ ! -d {MATEFILTER_OUTDIR} ]; then mkdir {MATEFILTER_OUTDIR}; fi"
# 		"&& samtools view -b -f 0x0040 {input} > {output.bam}"
# 		"&& samtools sort {output.bam} > {output.sorted_bam}"
# 		"&& samtools index {output.sorted_bam}"

# rule piranha:
# 	input:
# 		MATEFILTER_OUTDIR + "/{sample}_{replicate}_sorted.bam"
# 	output:
# 		PEAKCALLING_OUTDIR + "/{sample}_{replicate}_peaks.bed"
# 	conda:
# 		"envs/piranha.yml"
# 	threads: 2 
# 	shell:
# 		"if [ ! -d {PEAKCALLING_OUTDIR} ]; then mkdir {PEAKCALLING_OUTDIR}; fi"
# 		"&& Piranha {config[prianha]} {input} -o {output}"

# rule pureclip:
# 	input:
# 		experiment=DEDUPLICAITON_OUTDIR + "/{sample}_{replicate}_sorted.bam",
# 		experiment_bai=DEDUPLICAITON_OUTDIR + "/{sample}_{replicate}_sorted.bam.bai",
# 		genome_fasta=GENOME_FASTA
# 	output:
# 		crosslinking_sites=PEAKCALLING_OUTDIR + "/{sample}_{replicate}_crosslinkind_sites.bed",
# 		binding_regions=PEAKCALLING_OUTDIR + "/{sample}_{replicate}_binding_regions.bed",
# 	threads: 4
# 	params:
# 		tmp=PEAKCALLING_OUTDIR + "/tmp/",
# 		parameters=PEAKCALLING_OUTDIR + "/{sample}_{replicate}_parameters.txt"
# 	conda:
# 		"envs/pureclip.yml"
# 	shell:
# 		"if [ ! -d {PEAKCALLING_OUTDIR} ]; then mkdir {PEAKCALLING_OUTDIR}; fi "
# 		"&& pureclip -i {input.experiment} -bai {input.experiment_bai} -g {input.genome_fasta} -o {output.crosslinking_sites} -tmp {params.tmp} "
# 		"-or {output.binding_regions} -p {params.parameters} -nt {threads} -nta {threads} {config[pureclip]} "

rule peakachu:
    input:
    	clip=expand(DEDUPLICAITON_OUTDIR + "/{sample}_sorted.bam", sample=REPLICATES_CLIP),
    	control=expand(DEDUPLICAITON_OUTDIR + "/{sample}_sorted.bam", sample=REPLICATES_CONTROL)
    output:
    	peaks_tsv=PEAKCALLING_OUTDIR + "/peakachu.tsv",
    	peaks_gtf=PEAKCALLING_OUTDIR + "/peakachu.gtf",
    	blockbuster=PEAKCALLING_OUTDIR + "/blockbuster.bed"
    threads: 4
    conda:
    	"envs/peakachu.yml"
    shell:
    	#"--max_proc "${GALAXY_SLOTS:-1}"
    	"if [ ! -d {PEAKCALLING_OUTDIR} ]; then mkdir {PEAKCALLING_OUTDIR}; fi"
    	"&& source activate peakachu"
    	"&& peakachu adaptive --exp_libs {input.clip} --ctr_libs {input.control} {config[peakachu]} --output_folder {PEAKCALLING_OUTDIR} "
    	"&& source deactivate "
    	"&& head -q -n 1 {PEAKCALLING_OUTDIR}/peak_tables/*.csv > tmp.tsv "
    	"&& head -q -n 1 tmp.tsv > {output.peaks_tsv}"
    	"&& rm tmp.tsv"
    	"&& tail -n +2 -q {PEAKCALLING_OUTDIR}/peak_tables/*.csv >> {output.peaks_tsv} "
    	"&& cat {PEAKCALLING_OUTDIR}/peak_annotations/*.gff | awk '/peak/ {{ print $0 }}' > {output.peaks_gtf} "
    	"&& cat {PEAKCALLING_OUTDIR}/blockbuster_input/*bed > {output.blockbuster}"

rule peak_calling_quality:
	input:
		peaks=PEAKCALLING_OUTDIR + "/peakachu.gtf",
		clip=DEDUPLICAITON_OUTDIR + "/{sample}_{replicate}_sorted.bam"
	output:
		htseq_hits=PEAKCALLING_OUTDIR + "/{sample}_{replicate}_htseq_hits.txt",
		htseq_nohits=PEAKCALLING_OUTDIR + "/{sample}_{replicate}_htseq_nohits.txt",
		reads_in_peaks=PEAKCALLING_OUTDIR + "/{sample}_{replicate}_reads_summary.txt"
	threads: 2
	conda:
		"envs/htseq.yml"
	shell:
		"htseq-count --mode=union --stranded=yes --minaqual=10 --type='peak_region' --idattr='ID' --order=name --format=bam {input.clip} {input.peaks} "
		"""| awk '{{ if ($1 ~ "no_feature|ambiguous|too_low_aQual|not_aligned|alignment_not_unique") print $0 | "cat 1>&2"; else print $0 }}' > {output.htseq_hits} 2> {output.htseq_nohits} """
		"""&& awk 'FNR==NR{{ SUM1+=$2; next }} {{ SUM2+=$2 }} END {{ print "#reads in peaks: "SUM1; """
		"""print "#culled reads: "SUM2; print "percentage reads in peaks: " SUM1/(SUM1+SUM2) }}' {output.htseq_hits} {output.htseq_nohits} > {output.reads_in_peaks}"""

###############################
## COUNTING / COVERAGE FILES ##
###############################

rule extract_alignment_ends:
	input:
		DEDUPLICAITON_OUTDIR + "/{sample}_{replicate}_sorted.bam"
	output:
		COVERAGE_OUTDIR + "/{sample}_{replicate}_alignment_ends.bed"
	conda:
		"envs/bctools.yml"
	threads: 2 
	shell:
		"if [ ! -d {COVERAGE_OUTDIR} ]; then mkdir {COVERAGE_OUTDIR}; fi"
		"&& python {config[bctools]}/extract_aln_ends.py {input} > {output}"

rule extract_crosslinking_position:
	input:
		COVERAGE_OUTDIR + "/{sample}_{replicate}_alignment_ends.bed"
	output:
		COVERAGE_OUTDIR + "/{sample}_{replicate}_crosslinking_positions.bed"
	conda:
		"envs/bctools.yml"
	threads: 2
	shell:
		"if [ ! -d {COVERAGE_OUTDIR} ]; then mkdir {COVERAGE_OUTDIR}; fi"
		"&& python {config[bctools]}/coords2clnt.py {config[extract_cl]} {input} > {output}"

rule intersect_crosslink_sites_with_peaks:
    input:
    	peaks=PEAKCALLING_OUTDIR + "/peakachu.gtf",
    	cl=COVERAGE_OUTDIR + "/{sample}_{replicate}_crosslinking_positions.bed"
    output:
    	COVERAGE_OUTDIR + "/{sample}_{replicate}_crosslink_sites_intersecting_peaks.bed"
    threads: 2
    conda:
		"envs/bedtools.yml"
    shell:
    	"bedtools intersect -a {input.cl} -b {input.peaks} -s -u -wa  > {output}"

rule crosslink_sites_quality:
	input:
		cl=expand(COVERAGE_OUTDIR + "/{sample}_crosslinking_positions.bed", sample=REPLICATES_CLIP),
		ind=expand(COVERAGE_OUTDIR + "/{sample}_crosslink_sites_intersecting_peaks.bed", sample=REPLICATES_CLIP)
	output:
		expand(COVERAGE_OUTDIR + "/{sample}_crosslink_sites_quality.txt", sample=REPLICATES_CLIP)
	run:
		for i in range(0,len(input.cl)):
			total_crosslinking_events = sum(1 for line in open(input.cl[i]))
			crosslinking_events_in_peaks = sum(1 for line in open(input.ind[i]))
			shell("echo 'total crosslink sites: ' " + str(total_crosslinking_events) + " > " + output[i] + 
			"&& echo 'crosslink sites in peaks: ' " + str(crosslinking_events_in_peaks) + " >> " + output[i] + 
			"&& echo 'percentage of crosslinking sites in peaks: ' " + str(crosslinking_events_in_peaks/total_crosslinking_events) + " >> " + output[i])

rule correct_bed_format_check:
	input:
		cl=COVERAGE_OUTDIR + "/{sample}_{replicate}_crosslinking_positions.bed",
		alends=COVERAGE_OUTDIR + "/{sample}_{replicate}_alignment_ends.bed"
	output:
		cl=COVERAGE_OUTDIR + "/{sample}_{replicate}_crosslinking_positions_check.bed",
		alends=COVERAGE_OUTDIR + "/{sample}_{replicate}_alignment_ends_check.bed"
	threads: 2
	shell:
		"""awk -F "\t" 'BEGIN {{ OFS = FS }} $3 != 0 {{ print $0 }}' {input.cl} > {output.cl} """
		"""&& awk -F "\t" 'BEGIN {{ OFS = FS }} $3 != 0 {{ print $0 }}' {input.alends} > {output.alends} """

rule sort_beds:
	input:
		cl=COVERAGE_OUTDIR + "/{sample}_{replicate}_crosslinking_positions_check.bed",
		alends=COVERAGE_OUTDIR + "/{sample}_{replicate}_alignment_ends_check.bed"
	output:
		cl=COVERAGE_OUTDIR + "/{sample}_{replicate}_crosslinking_positions_sorted.bed",
		alends=COVERAGE_OUTDIR + "/{sample}_{replicate}_alignment_ends_sorted.bed"
	threads: 2
	conda:
		"envs/bedtools.yml"
	shell:
		"bedtools sort -i {input.cl} > {output.cl}"
		"&& bedtools sort -i {input.alends} > {output.alends}"

rule calculate_coverage:
	input:
		cl=COVERAGE_OUTDIR + "/{sample}_{replicate}_crosslinking_positions_sorted.bed",
		alends=COVERAGE_OUTDIR + "/{sample}_{replicate}_alignment_ends_sorted.bed"
	output:
		cl_pos=COVERAGE_OUTDIR + "/bedgraph/{sample}_{replicate}_crosslinking_coverage_pos_strand.bedgraph",
		cl_neg=COVERAGE_OUTDIR + "/bedgraph/{sample}_{replicate}_crosslinking_coverage_neg_strand.bedgraph",
		cl_bot=COVERAGE_OUTDIR + "/bedgraph/{sample}_{replicate}_crosslinking_coverage_both_strand.bedgraph",
		alends_pos=COVERAGE_OUTDIR + "/bedgraph/{sample}_{replicate}_alignment_ends_coverage_pos_strand.bedgraph",
		alends_ned=COVERAGE_OUTDIR + "/bedgraph/{sample}_{replicate}_alignment_ends_coverage_neg_strand.bedgraph",
		alends_bot=COVERAGE_OUTDIR + "/bedgraph/{sample}_{replicate}_alignment_ends_coverage_both_strand.bedgraph"
	threads: 4
	conda:
		"envs/bedtools.yml"
	shell:
		"if [ ! -d {COVERAGE_OUTDIR}/bedgraph ]; then mkdir {COVERAGE_OUTDIR}/bedgraph; fi"
		"&& genomeCoverageBed -i {input.cl} -g {GENOME_SIZES} -bg -strand + > {output.cl_pos}" 
		"&& genomeCoverageBed -i {input.cl} -g {GENOME_SIZES} -bg -strand - > {output.cl_neg}" 
		"&& genomeCoverageBed -i {input.cl} -g {GENOME_SIZES} -bg > {output.cl_bot}" 
		"&& genomeCoverageBed -i {input.alends} -g {GENOME_SIZES} -bg -strand + > {output.alends_pos}"
		"&& genomeCoverageBed -i {input.alends} -g {GENOME_SIZES} -bg -strand - > {output.alends_ned}" 
		"&& genomeCoverageBed -i {input.alends} -g {GENOME_SIZES} -bg > {output.alends_bot}" 

rule calculate_coverage_bigwig:
	input:
		cl=COVERAGE_OUTDIR + "/bedgraph/{sample}_{replicate}_crosslinking_coverage_{type}_strand.bedgraph",
		alends=COVERAGE_OUTDIR + "/bedgraph/{sample}_{replicate}_alignment_ends_coverage_{type}_strand.bedgraph"
	output:
		cl=COVERAGE_OUTDIR + "/bigwig/{sample}_{replicate}_crosslinking_coverage_{type}_strand.bigwig",
		alends=COVERAGE_OUTDIR + "/bigwig/{sample}_{replicate}_alignment_ends_coverage_{type}_strand.bigwig"
	threads: 4
	conda:
		"envs/bedGraphToBigWig.yml"
	shell:
		"if [ ! -d {COVERAGE_OUTDIR}/bigwig ]; then mkdir {COVERAGE_OUTDIR}/bigwig; fi "
		"&& bedGraphToBigWig {input.cl} {GENOME_SIZES} {output.cl}"
		"&& bedGraphToBigWig {input.alends} {GENOME_SIZES} {output.alends}"

##############################
## POST-PROCESSING OF PEAKS ##
##############################

rule peaks_tsv_to_bed:
	input:
		PEAKCALLING_OUTDIR + "/peakachu.tsv"
	output:
		POSTPROCESSING_OUTDIR + "/peakachu.bed"
	threads: 2
	shell:
		"if [ ! -d {POSTPROCESSING_OUTDIR} ]; then mkdir {POSTPROCESSING_OUTDIR}; fi"
		"""&&  awk -F "\t" 'BEGIN {{ OFS = FS }} NR>1 {{ if ($3 < $4) {{ print $1,$3,$4,"clip_peak_"NR-1,$9,$5; }} else {{ print $1,$4,$3,"clip_peak_"NR-1,$9,$5; }} }}' {input} > {output} """

rule peaks_extend_frontiers:
	input:
		bed=POSTPROCESSING_OUTDIR + "/peakachu.bed",
		genome=GENOME_SIZES
	output:
		POSTPROCESSING_OUTDIR + "/peakachu_extended.bed"
	threads: 2
	conda:
		"envs/bedtools.yml"
	shell:
		"bedtools slop {config[peaks_extend_frontiers]} -i {input.bed} -g {input.genome} > {output}"

#####################
## MOTIF DETECTION ##
#####################

rule extract_genomic_DNA_dreme:
	input:
		POSTPROCESSING_OUTDIR + "/peakachu_extended.bed"
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
		"meme-chip {input} -noecho -dna -oc {MOTIF_DETECTION_OUTDIR}/meme_chip {config[meme_chip]}"

# for RCAS please use ensembl annotation gtfs
# rule rcas:
# 	input:
# 		POSTPROCESSING_OUTDIR + "/peakachu_extended.bed"
# 	output:
# 		MOTIF_DETECTION_OUTDIR + "/rcas/rcas_summary.html"
# 	threads: 4
# 	shell:
# 		"if [ ! -d {MOTIF_DETECTION_OUTDIR}/rcas ]; then mkdir {MOTIF_DETECTION_OUTDIR}/rcas; fi"
# 		"&& Rscript " + config["rcas"] + "/RCAS.R {input} {GENOME_GTF} ... 'hg19' {SRC_PATH}/{MOTIF_DETECTION_OUTDIR}/rcas 0 "
# 		"&& mv {MOTIF_DETECTION_OUTDIR}/rcas/*.html {MOTIF_DETECTION_OUTDIR}/rcas/rcas_summary.html"

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
#  		"envs/meme_suite.yml"
#  	shell:
#  		"if [ ! -d {MOTIF_SEARCH_OUTDIR} ]; then mkdir {MOTIF_SEARCH_OUTDIR}; fi"
#  		"&& fimo --alpha 1.000000 --max-stored-scores 100000 --motif-pseudo 0.100000 --qv-thresh --thresh 0.000100 --verbosity 1 "
#  		"--oc {MOTIF_SEARCH_OUTDIR}/fimo_dreme {input} {GENOME_FASTA}"

# rule fimo_for_meme_output:
#  	input: 
#  		MOTIF_DETECTION_OUTDIR + "/meme/{peak_file}_meme_output/meme.xml"
#  	output:
#  		MOTIF_SEARCH_OUTDIR + "/fimo_meme/{peak_file}_meme_output"
#  	threads: 4
#  	conda:
#  		"envs/meme_suite.yml"
#  	shell:
#  		"fimo --alpha 1.000000 --max-stored-scores 100000 --motif-pseudo 0.100000 --qv-thresh --thresh 0.000100 --verbosity 1 "
#  		"--oc {output} {input} {GENOME_FASTA}"
