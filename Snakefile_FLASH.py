import random
import math
import os 

SRC_PATH = os.getcwd()

configfile: "config.yml"

REF_GENOME_DIR = config['ref_genome_dir']
GENOME_FASTA = REF_GENOME_DIR + config['genome_fasta']
GENOME_2BIT = REF_GENOME_DIR + config['genome_2bit']
GENOME_GTF = REF_GENOME_DIR + config['genome_gtf']
GENOME_SIZES = REF_GENOME_DIR + config['genome_sizes']

control = "no"
demultiplexed = "no"

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

######################
## SAMPLE VARIABLES ##
######################

FIRST_READS = list()
SECOND_READS = list()
CLIP_SAMPLES = list()
CONTROL_SAMPLES = list()
REPLICATES_CONTROL = list()

if ( control == "yes" ):
	FIRST_READS = list(find_all_values_on_key("r1", config['clip_samples'])) + list(find_all_values_on_key("r1", config['control_samples']))
	SECOND_READS = list(find_all_values_on_key("r2", config['clip_samples'])) + list(find_all_values_on_key("r2", config['control_samples']))
	CONTROL_SAMPLES = list(list_all_values_of_dict(config['control_samples']))
	REPLICATES_CONTROL = list(config['control_samples'].keys())
else:
	FIRST_READS = list(find_all_values_on_key("r1", config['clip_samples']))
	SECOND_READS = list(find_all_values_on_key("r2", config['clip_samples']))

CLIP_SAMPLES = list(list_all_values_of_dict(config['clip_samples']))
ALL_SAMPLES = CLIP_SAMPLES + CONTROL_SAMPLES

REPLICATES_CLIP = list(config['clip_samples'].keys())
ALL_REPLICATES = REPLICATES_CLIP + REPLICATES_CONTROL

PAIR = ["r1", "r2"]
REP_NAME_CLIP = ["rep" + str(x+1) for x in range(0,len(REPLICATES_CLIP))]
REP_NAME_CONTROL = ["rep" + str(x+1) for x in range(0,len(REPLICATES_CONTROL))]

SAMPLES = []
if ( control == "yes" ):
	SAMPLES = [REPLICATES_CLIP[0].split("_")[0], REPLICATES_CONTROL[0].split("_")[0]]
else:
	SAMPLES = [REPLICATES_CLIP[0].split("_")[0]]

print("[NOTE] Loaded Samples.")
print("All samples: " + str(ALL_SAMPLES))
print("CLIP samples: " + str(CLIP_SAMPLES))
print("Conrol samples: " + str(CONTROL_SAMPLES))
print("First read samples: " + str(FIRST_READS))
print("Second read samples: " + str(SECOND_READS))
print(SAMPLES)

######################
## PATH DEFINITIONS ##
######################

RENAMING = config['sample_data_dir'] + "/" + config['renaming_outdir']
FASTQC_BEG_OUTDIR = config['sample_data_dir'] + "/" + config['fastqc_beg_outdir']
FASTQC_ADAPT_OUTDIR = config['sample_data_dir'] + "/" + config['fastqc_adapt_outdir']
FLEXBAR_OUTDIR = config['sample_data_dir'] + "/" + config['flexbar_outdir']
CUTADAPT_OUTDIR = config['sample_data_dir'] + "/" + config['cutadapt_outdir']
REMOVE_TAIL_OUTDIR = config['sample_data_dir'] + "/" + config['remove_tail_outdir']
MAPPING_OUTDIR = config['sample_data_dir'] + "/" + config['mapping_outdir']
MAPPING_QUALITY_OUTDIR = config['sample_data_dir'] + "/" + config['mapping_quality_outdir']
PRE_FOR_UMI_OUTDIR = config['sample_data_dir'] + "/" + config['preprocessing_for_umi_outdir']
DEDUPLICAITON_OUTDIR = config['sample_data_dir'] + "/" + config['deduplication_outdir']
FASTQC_DEDUP_OUTDIR = config['sample_data_dir'] + "/" + config['fastqc_dedup_outdir']
COVERAGE_OUTDIR = config['sample_data_dir'] + "/" + config['coverage_outdir']

PEAKCALLING_OUTDIR = config['sample_data_dir'] + "/" + config['peakcalling_outdir']
ANNOTATION_PEAKS_OUTDIR = config['sample_data_dir'] + "/" + config['annotation_peaks_outdir']
HTSEQ_COUNT_OUTDIR = config['sample_data_dir'] + "/" + config['htseq_outdir']
DEDUPLICAITON_OUTDIR = config['sample_data_dir'] + "/" + config['deduplication_outdir']
MATEFILTER_OUTDIR = config['sample_data_dir'] + "/" + config['matefilter_outdir']

MULTIQC_OUTDIR = config['sample_data_dir'] + "/" + config['multiqc_outdir']

###################
## PREPROCESSING ##
###################

ALL_NEW_FILE_NAMES = []
if ( control == "yes" ):
	rule all:
		input: 
			expand(FASTQC_BEG_OUTDIR + "/{sample}_{replicate}_{pair}.fastqsanger_fastqc.html", sample=SAMPLES[0], replicate=REP_NAME_CLIP, pair=PAIR),
			expand(FASTQC_BEG_OUTDIR + "/{sample}_{replicate}_{pair}.fastqsanger_fastqc.html", sample=SAMPLES[1], replicate=REP_NAME_CONTROL, pair=PAIR),
			expand(FASTQC_ADAPT_OUTDIR + "/{sample}_{replicate}_{pair}.fastqsanger_fastqc.html", sample=SAMPLES[0], replicate=REP_NAME_CLIP, pair=PAIR),
			expand(FASTQC_ADAPT_OUTDIR + "/{sample}_{replicate}_{pair}.fastqsanger_fastqc.html", sample=SAMPLES[1], replicate=REP_NAME_CONTROL, pair=PAIR),
			REF_GENOME_DIR + "/sjdbList.fromGTF.out.tab",
			expand(MAPPING_OUTDIR + "/{sample}_{replicate}.bam", sample=SAMPLES[0], replicate=REP_NAME_CLIP),
			expand(MAPPING_OUTDIR + "/{sample}_{replicate}.bam", sample=SAMPLES[1], replicate=REP_NAME_CONTROL),
			expand(MAPPING_OUTDIR + "/{sample}_{replicate}.bam.bai", sample=SAMPLES[0], replicate=REP_NAME_CLIP),
			expand(MAPPING_OUTDIR + "/{sample}_{replicate}.bam.bai", sample=SAMPLES[1], replicate=REP_NAME_CONTROL),
			expand(DEDUPLICAITON_OUTDIR + "/{sample}_{replicate}_name_sorted.bam", sample=SAMPLES[0], replicate=REP_NAME_CLIP),
			expand(DEDUPLICAITON_OUTDIR + "/{sample}_{replicate}_name_sorted.bam", sample=SAMPLES[1], replicate=REP_NAME_CONTROL),
			expand(FASTQC_DEDUP_OUTDIR + "/{sample}_{replicate}_fastqc.html", sample=SAMPLES[0], replicate=REP_NAME_CLIP),
			expand(FASTQC_DEDUP_OUTDIR + "/{sample}_{replicate}_fastqc.html", sample=SAMPLES[1], replicate=REP_NAME_CONTROL),
			MAPPING_QUALITY_OUTDIR + "/fingerprint_plot.png",
			expand(COVERAGE_OUTDIR + "/bigwig/{sample}_{replicate}_crosslinking_coverage_{type}_strand.bigwig", sample=SAMPLES[0], replicate=REP_NAME_CLIP, type=["pos", "neg", "both"]),
			expand(COVERAGE_OUTDIR + "/bigwig/{sample}_{replicate}_crosslinking_coverage_{type}_strand.bigwig", sample=SAMPLES[1], replicate=REP_NAME_CONTROL, type=["pos", "neg", "both"]),
			expand(ANNOTATION_PEAKS_OUTDIR + "/{sample_exp}_{replicate_exp}_{sample_ctl}_{replicate_ctl}_binding_regions_intersecting_peaks.gtf", 
					sample_exp=SAMPLES[0], replicate_exp=REP_NAME_CLIP, sample_ctl=SAMPLES[1], replicate_ctl=REP_NAME_CONTROL),
			expand(HTSEQ_COUNT_OUTDIR + "/{sample}_{replicate}_htseq_hits.txt", sample=SAMPLES[0], replicate=REP_NAME_CLIP),
			expand(HTSEQ_COUNT_OUTDIR + "/{sample}_{replicate}_htseq_hits.txt", sample=SAMPLES[1], replicate=REP_NAME_CONTROL),
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
else:
	rule all:
		input: 
			expand(FASTQC_BEG_OUTDIR + "/{sample}_{replicate}_{pair}.fastqsanger_fastqc.html", sample=SAMPLES[0], replicate=REP_NAME_CLIP, pair=PAIR),
			expand(FASTQC_ADAPT_OUTDIR + "/{sample}_{replicate}_{pair}.fastqsanger_fastqc.html", sample=SAMPLES[0], replicate=REP_NAME_CLIP, pair=PAIR),
			REF_GENOME_DIR + "/sjdbList.fromGTF.out.tab",
			expand(MAPPING_OUTDIR + "/{sample}_{replicate}.bam", sample=SAMPLES[0], replicate=REP_NAME_CLIP),
			expand(MAPPING_OUTDIR + "/{sample}_{replicate}.bam.bai", sample=SAMPLES[0], replicate=REP_NAME_CLIP),
			expand(DEDUPLICAITON_OUTDIR + "/{sample}_{replicate}_name_sorted.bam", sample=SAMPLES[0], replicate=REP_NAME_CLIP),
			expand(FASTQC_DEDUP_OUTDIR + "/{sample}_{replicate}_fastqc.html", sample=SAMPLES[0], replicate=REP_NAME_CLIP),
			MAPPING_QUALITY_OUTDIR + "/fingerprint_plot.png",
			expand(COVERAGE_OUTDIR + "/bigwig/{sample}_{replicate}_crosslinking_coverage_{type}_strand.bigwig", sample=SAMPLES[0], replicate=REP_NAME_CLIP, type=["pos", "neg", "both"]),
			expand(ANNOTATION_PEAKS_OUTDIR + "/{sample_exp}_{replicate_exp}_{sample_ctl}_{replicate_ctl}_binding_regions_intersecting_peaks.gtf", 
					sample_exp=SAMPLES[0], replicate_exp=REP_NAME_CLIP, sample_ctl=SAMPLES[1], replicate_ctl=REP_NAME_CONTROL),
			expand(HTSEQ_COUNT_OUTDIR + "/{sample}_{replicate}_htseq_hits.txt", sample=SAMPLES[0], replicate=REP_NAME_CLIP),
			MOTIF_DETECTION_OUTDIR + "/meme_chip/meme-chip.html"

	ALL_NEW_FILE_NAMES = ["none"] * (len(REPLICATES_CLIP) * 2)
	i = 0 
	for j in REP_NAME_CLIP:
		for k in PAIR:
			ALL_NEW_FILE_NAMES[i] = RENAMING + "/" + SAMPLES[0] + "_" + j + "_" + k + ".fastqsanger"
			i += 1

rule renaming:
	input:
		expand(config['sample_data_dir'] + "/{all}.fastq", all=ALL_SAMPLES)
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

# flexbar is not working on the cluster !!!
# rule flexbar:
# 	input:
# 		first=RENAMING + "/{samples}_{replicate}_r1.fastqsanger",
# 		second=RENAMING + "/{samples}_{replicate}_r2.fastqsanger",
# 		barcodes=config['sample_data_dir'] + "/" + config['barcodes']
# 	output:
# 		first=FLEXBAR_OUTDIR + "/{samples}_{replicate}_r1.fastqsanger",
# 		second=FLEXBAR_OUTDIR + "/{samples}_{replicate}_r2.fastqsanger",
# 	threads: 2
# 	conda:
# 		"envs/flexbar.yml"
# 	shell:
# 		"if [ ! -d {FLEXBAR_OUTDIR} ]; then mkdir {FLEXBAR_OUTDIR}; fi"
# 		"&& flexbar --threads {threads} --reads {input.first} --reads2 {input.second} --max-uncalled 0 --min-read-length 18 "
# 		"--adapters {input.barcodes} --adapter-trim-end RIGHT --adapter-min-overlap 1 --adapter-threshold 3"

rule cutadapt_first_read_clip:
	input:
		first=RENAMING + "/{samples}_{replicate}_r1.fastqsanger",
		second=RENAMING + "/{samples}_{replicate}_r2.fastqsanger"
	output:
		seq_first=CUTADAPT_OUTDIR + "/{samples}_{replicate}_r1.fastqsanger",
		seq_second=CUTADAPT_OUTDIR + "/{samples}_{replicate}_r2.fastqsanger",
		log=CUTADAPT_OUTDIR + "/{samples}_{replicate}.txt"
	threads: 2 
	conda:
		"envs/cutadapt.yml"
	shell:
		"if [ ! -d {CUTADAPT_OUTDIR} ]; then mkdir {CUTADAPT_OUTDIR}; fi"
		"&& cutadapt -j {threads} {config[cutadapt]} " 
		"--paired-output={output.seq_second} --output={output.seq_first} {input.first} {input.second} > {output.log}"

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
		"&& python {config[bctools]}/remove_tail.py {input.first} 13  > {output.first}"
		"&& cp {input.second} {output.second}"

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
# 		"&& umi_tools extract -p NNNNNXXXXXXNN -I {input.second} -S {params.second} --read2-in {input.first} --read2-out {params.first} -L {params.log}"
# 		"&& umi_tools extract -p NNNNNN -I {params.second} -S {output.second} --read2-in {params.first} --read2-out {output.first} -L {output.log}"

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
		"--outSAMtype BAM SortedByCoordinate --outSAMattributes All --outSAMstrandField intronMotif --outFilterIntronMotifs None --outSAMunmapped None " 
		"--outSAMprimaryFlag OneBestScore --outSAMmapqUnique '255' --outFilterType Normal --outFilterMultimapScoreRange '1' --outFilterMultimapNmax '10' "
		"--outFilterMismatchNmax '10' --outFilterMismatchNoverLmax '0.3' --outFilterMismatchNoverReadLmax '1.0' --outFilterScoreMin '0' --outFilterScoreMinOverLread '0.66' " 
		"--outFilterMatchNmin '0' --outFilterMatchNminOverLread '0.66'   --seedSearchStartLmax '50' --seedSearchStartLmaxOverLread '1.0' --seedSearchLmax '0' " 
		"--seedMultimapNmax '10000' --seedPerReadNmax '1000' --seedPerWindowNmax '50' --seedNoneLociPerWindow '10'  --alignIntronMin '21' --alignIntronMax '0' " 
		"--alignMatesGapMax '0' --alignSJoverhangMin '5' --alignSJDBoverhangMin '3' --alignSplicedMateMapLmin '0' --alignSplicedMateMapLminOverLmate '0.66' " 
		"--alignWindowsPerReadNmax '10000' --alignTranscriptsPerWindowNmax '100' --alignTranscriptsPerReadNmax '10000' --alignEndsType EndToEnd  --twopassMode 'Basic' " 
		"--twopass1readsN '-1' --limitBAMsortRAM '0' --limitOutSJoneRead '1000' --limitOutSJcollapsed '1000000' --limitSjdbInsertNsj '1000000' > {output.log} "
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

rule filter_out_unlocalized_regions_for_later_genome_versions:
	input:
		PRE_FOR_UMI_OUTDIR + "/{sample}_{replicate}_unique_reads_fitlering.bam"
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
 
# rule barcode_header_removal:
# 	input:
# 		PRE_FOR_UMI_OUTDIR + "/{sample}_{replicate}_got_umis_unlocalized_check.bam"
# 	output:
# 		PRE_FOR_UMI_OUTDIR + "/{sample}_{replicate}_got_umis_unlocalized_check_barcode_removed.bam"
# 	threads: 2
#	conda:
#		"envs/samtools.yml"
# 	shell:
# 		"samtools view -h {input} "
# 		"| sed s/\_[ATGCN]*//2 "
# 		"| samtools view -bSh > {output}"
# 		"&& samtools index {output}"
 
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
		"&& umi_tools dedup --random-seed 0 --extract-umi-method read_id --method adjacency --edit-distance-threshold 1 --paired " 
		"--soft-clip-threshold 4 --subset 1.0 -I {input} -S {output.bam} -L {output.log}"
		"&& samtools sort {output.bam} > {output.sorted_bam}"
		"&& samtools index {output.sorted_bam}"

rule sam_name_sort:
	input:
		sorted_bam=DEDUPLICAITON_OUTDIR + "/{sample}_{replicate}_sorted.bam"
	output:
		name_sorted_bam=DEDUPLICAITON_OUTDIR + "/{sample}_{replicate}_name_sorted.bam"
	threads: 2
	conda:
		"envs/samtools.yml"
	params:
		tmp=DEDUPLICAITON_OUTDIR + "/{sample}_{replicate}"
	shell:
		"samtools sort -n -T {params.tmp} -o {output.name_sorted_bam} {input.sorted_bam}" 

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
		"&& computeGCBias --numberOfProcessors {threads} --bamfile {input} --GCbiasFrequenciesFile {output.file} --fragmentLength 300 "
		"--genome {GENOME_2BIT} --effectiveGenomeSize 2451960000 --biasPlot {output.plot} --plotFileFormat png"

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
		"HISTOGRAM_FILE={output.plot} DEVIATIONS='10.0' MINIMUM_PCT='0.05' ASSUME_SORTED='false' " 
		"METRIC_ACCUMULATION_LEVEL='ALL_READS'  VALIDATION_STRINGENCY='LENIENT' QUIET=true VERBOSITY=ERROR "
		"&& source deactivate"

rule finger_print_plot:
	input:
		expand(DEDUPLICAITON_OUTDIR + "/{sample}.bam", sample=ALL_REPLICATES)
	output:
		MAPPING_QUALITY_OUTDIR + "/fingerprint_plot.png"
	threads: 2
	params:
		binsize=100
	conda:
		"envs/deeptools.yml"
	shell:	
		"if [ ! -d {MAPPING_QUALITY_OUTDIR} ]; then mkdir {MAPPING_QUALITY_OUTDIR}; fi"
		"&& plotFingerprint --numberOfProcessors {threads} --bamfiles {input} "
		"--binSize {params.binsize} --labels {ALL_REPLICATES} --plotFile {output}  --plotFileFormat 'png'"

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
		"&& python {config[bctools]}/coords2clnt.py {input} > {output}"

# Rule to check if the end coordinate is 0. It happened in a test-run.
# For a region is has to correct start and end coordinates.
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

# rule pureclip:
# 	input:
# 		experiment=MATEFILTER_OUTDIR + "/{sample_exp}_{replicate_exp}_sorted.bam",
# 		experiment_bai=MATEFILTER_OUTDIR + "/{sample_exp}_{replicate_exp}_sorted.bam.bai",
# 		control=MATEFILTER_OUTDIR + "/{sample_ctl}_{replicate_ctl}_sorted.bam",
# 		control_bai=MATEFILTER_OUTDIR + "/{sample_ctl}_{replicate_ctl}_sorted.bam.bai",
# 		genome_fasta=GENOME_FASTA
# 	output:
# 		crosslinking_sites=PEAKCALLING_OUTDIR + "/{sample_exp}_{replicate_exp}_{sample_ctl}_{replicate_ctl}_crosslinkind_sites.bed",
# 		binding_regions=PEAKCALLING_OUTDIR + "/{sample_exp}_{replicate_exp}_{sample_ctl}_{replicate_ctl}_binding_regions.bed",
# 	threads: 4
# 	params:
# 		tmp=PEAKCALLING_OUTDIR + "/tmp/",
# 		parameters=PEAKCALLING_OUTDIR + "/{sample_exp}_{replicate_exp}_{sample_ctl}_{replicate_ctl}_parameters.txt"
# 	conda:
# 		"envs/pureclip.yml"
# 	shell:
# 		"if [ ! -d {PEAKCALLING_OUTDIR} ]; then mkdir {PEAKCALLING_OUTDIR}; fi "
# 		"&& pureclip -i {input.experiment} -bai {input.experiment_bai} -g {input.genome_fasta} -o {output.crosslinking_sites} -tmp {params.tmp} "
# 		"-ibam {input.control} -ibai {input.control_bai} -or {output.binding_regions} -p {params.parameters} -nt {threads} -nta {threads} -iv 'chr1;chr2;chr3;' "

rule peakachu:
	input:
		clip=DEDUPLICAITON_OUTDIR + "/{sample_exp}_{replicate_exp}_sorted.bam",
		control=DEDUPLICAITON_OUTDIR + "/{sample_ctl}_{replicate_ctl}_sorted.bam"
	output:
		peaks_tsv=PEAKCALLING_OUTDIR + "/{sample_exp}_{replicate_exp}_{sample_ctl}_{replicate_ctl}/peakachu.tsv",
		peaks_gtf=PEAKCALLING_OUTDIR + "/{sample_exp}_{replicate_exp}_{sample_ctl}_{replicate_ctl}/peakachu.gtf",
		blockbuster=PEAKCALLING_OUTDIR + "/{sample_exp}_{replicate_exp}_{sample_ctl}_{replicate_ctl}/blockbuster.bed"
	threads: 4
	conda:
		"envs/peakachu.yml"
	params:
		insert_size=200,
		mad_multiplier=0.0,
		fold_change_cutoff=1.0,
		q_value_cutoff=0.05,
		output_folder=PEAKCALLING_OUTDIR + "/{sample_exp}_{replicate_exp}_{sample_ctl}_{replicate_ctl}"
	shell:
		#"--max_proc "${GALAXY_SLOTS:-1}"
		"if [ ! -d {PEAKCALLING_OUTDIR} ]; then mkdir {PEAKCALLING_OUTDIR}; fi"
		"&& if [ ! -d {params.output_folder} ]; then mkdir {params.output_folder} ; fi"
		"&& source activate peakachu"
		"&& peakachu adaptive --exp_libs {input.clip} --ctr_libs {input.control} "
		"--paired_end --max_insert_size {params.insert_size} --features '' --sub_features '' "
		"--output_folder {PEAKCALLING_OUTDIR} --min_cluster_expr_frac 0.01 --min_block_overlap 0.5 --min_max_block_expr 0.1 "
		"--norm_method deseq --mad_multiplier {params.mad_multiplier} --fc_cutoff {params.fold_change_cutoff} --padj_threshold {params.q_value_cutoff}"
		"&& source deactivate "
		"&& head -q -n 1 {PEAKCALLING_OUTDIR}/peak_tables/*.csv > tmp.tsv "
		"&& head -q -n 1 tmp.tsv > {output.peaks_tsv}"
		"&& rm tmp.tsv"
		"&& tail -n +2 -q {PEAKCALLING_OUTDIR}/peak_tables/*.csv >> {output.peaks_tsv} "
		"&& cat {PEAKCALLING_OUTDIR}/peak_annotations/*.gff | awk '/peak/ {{ print $0 }}' > {output.peaks_gtf} "
		"&& cat {PEAKCALLING_OUTDIR}/blockbuster_input/*bed > {output.blockbuster}"

#####################
## PEAK ANNOTATION ##
#####################

# PURECLIP
# rule intersect_binding_regions_with_peaks:
#     input:
#     	binding_regions=PEAKCALLING_OUTDIR + "/{sample_exp}_{replicate_exp}_{sample_ctl}_{replicate_ctl}_binding_regions.bed",
#     	annotation=GENOME_GTF
#     output:
#     	ANNOTATION_PEAKS_OUTDIR + "/{sample_exp}_{replicate_exp}_{sample_ctl}_{replicate_ctl}_binding_regions_intersecting_peaks.gtf"
#     threads: 2
#     conda:
# 		"envs/bedtools.yml"	
#     shell:
#     	"if [ ! -d {ANNOTATION_PEAKS_OUTDIR} ]; then mkdir {ANNOTATION_PEAKS_OUTDIR}; fi "
#     	"&& bedtools intersect -a {input.annotation} -b {input.binding_regions} -s -u -wa > {output}"

# PEAKACHU
rule intersect_binding_regions_with_peaks:
    input:
    	binding_regions=PEAKCALLING_OUTDIR + "/{sample_exp}_{replicate_exp}_{sample_ctl}_{replicate_ctl}/peakachu.gtf",
    	annotation=GENOME_GTF
    output:
    	ANNOTATION_PEAKS_OUTDIR + "/{sample_exp}_{replicate_exp}_{sample_ctl}_{replicate_ctl}_binding_regions_intersecting_peaks.gtf"
    threads: 2
    conda:
		"envs/bedtools.yml"	
    shell:
    	"if [ ! -d {ANNOTATION_PEAKS_OUTDIR} ]; then mkdir {ANNOTATION_PEAKS_OUTDIR}; fi "
    	"&& bedtools intersect -a {input.annotation} -b {input.binding_regions} -s -u -wa -f 0.1 > {output}"

####################
## HTSEQ ANALYSIS ##
####################

rule peak_calling_quality:
	input:
		reads=DEDUPLICAITON_OUTDIR + "/{sample}_{replicate}_name_sorted.bam",
		annotation=GENOME_GTF
	output:
		htseq_hits=HTSEQ_COUNT_OUTDIR + "/{sample}_{replicate}_htseq_hits.txt",
		htseq_nohits=HTSEQ_COUNT_OUTDIR + "/{sample}_{replicate}_htseq_nohits.txt",
		reads_in_features=HTSEQ_COUNT_OUTDIR + "/{sample}_{replicate}_reads_summary.txt"
	threads: 2
	conda:
		"envs/htseq.yml"
	shell:
		"if [ ! -d {HTSEQ_COUNT_OUTDIR} ]; then mkdir {HTSEQ_COUNT_OUTDIR}; fi "
		"&& htseq-count --mode=union --stranded=yes --minaqual=10 --type='transcript' --idattr='gene_id' --nonunique=none --secondary-alignments=ignore "
		"--supplementary-alignments=ignore --order=name --format=bam {input.reads} {input.annotation}"
		"""| awk '{{ if ($1 ~ "no_feature|ambiguous|too_low_aQual|not_aligned|alignment_not_unique") print $0 | "cat 1>&2"; else print $0 }}' > {output.htseq_hits} 2> {output.htseq_nohits} """
		"""&& awk 'FNR==NR{{ SUM1+=$2; next }} {{ SUM2+=$2 }} END {{ print "#reads in features: "SUM1; """
		"""print "#culled reads: "SUM2; print "percentage reads in features: " SUM1/(SUM1+SUM2) }}' {output.htseq_hits} {output.htseq_nohits} > {output.reads_in_features}"""

#################
## END QUALITY ##
#################

rule multiqc:
    input:
    	beginning=expand(RENAMING + "/{sample}_{replicate}_{pair}.fastqsanger", sample=SAMPLES[0], replicate=REP_NAME_CLIP, pair=PAIR),
    	trimming=expand(FASTQC_ADAPT_OUTDIR + "/{sample}_{replicate}_{pair}.fastqsanger_fastqc.html", sample=SAMPLES[0], replicate=REP_NAME_CLIP, pair=PAIR),
    	dedup=expand(FASTQC_DEDUP_OUTDIR + "/{sample}_{replicate}_got_umis_unlocalized_check_fastqc.html", sample=SAMPLES[0], replicate=REP_NAME_CLIP),
    	mapping=expand(MAPPING_OUTDIR + "/{sample}_{replicate}_Log.final.out", sample=SAMPLES[0], replicate=REP_NAME_CLIP)
    output:
    	MULTIQC_OUTDIR + "multiqc_report.html"
    threads: 2
    conda:
    	"envs/multiqc.yml"
    shell:
    	"if [ ! -d {MULTIQC_OUTDIR} ]; then mkdir {MULTIQC_OUTDIR}; fi"
    	"&& multiqc {input.beginning} {input.trimming} {input.dedup} {input.mapping} --outdir {MULTIQC_OUTDIR}"