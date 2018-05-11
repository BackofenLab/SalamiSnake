import random
import math

configfile: "config.yml"

REF_MAPPING="genomes/hg19"

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

num_peak_files = 5
PEAK_FILES_FOR_MEME = [""] * num_peak_files
for i in range(0,5):
	peak_file = "peak_random_sample_" + str(i+1)
	PEAK_FILES_FOR_MEME[i] = peak_file + "/" + peak_file

######################
## PATH DEFINITIONS ##
######################

FASTQC_FIRST_OUTDIR = config['sample_data_dir'] + "/" + config['fastqc_first_outdir']
CUTADAPT_OUTDIR = config['sample_data_dir'] + "/" + config['cutadapt_outdir']
TRIMGALORE_OUTDIR = config['sample_data_dir'] + "/" + config['trim_galore_outdir']
MAPPING_OUTDIR = config['sample_data_dir'] + "/" + config['mapping_outdir']
PRE_FOR_UMI_OUTDIR = config['sample_data_dir'] + "/" + config['preprocessing_for_umi_outdir']
DEDUPLICAITON_OUTDIR = config['sample_data_dir'] + "/" + config['deduplication_outdir']
FASTQC_SECOND_OUTDIR = config['sample_data_dir'] + "/" + config['fastqc_second_outdir']
PEAKCALLING_OUTDIR = config['sample_data_dir'] + "/" + config['peakcalling_outdir']
POSTPROCESSING_OUTDIR = config['sample_data_dir'] + "/" + config['post_processing_outdir']
MOTIF_DETECTION_OUTDIR = config['sample_data_dir'] + "/" + config['motif_detection_outdir']
MOTIF_SEARCH_OUTDIR = config['sample_data_dir'] + "/" + config['motif_search_outdir']

###################
## PREPROCESSINg ##
###################

# rule all:
#   input: MAPPING_OUTDIR + "/{sample}.sam"

rule fastqc_first:
    input:
    	expand(config['sample_data_dir'] + "/{sample}.fastqsanger", sample=ALL_SAMPLES)
    output:
    	expand(FASTQC_FIRST_OUTDIR + "/{sample}.fastqsanger_fastqc.html", sample=ALL_SAMPLES),
    	expand(FASTQC_FIRST_OUTDIR + "/{sample}.fastqsanger_fastqc.zip", sample=ALL_SAMPLES)
    shell:
    	"if [ ! -d {FASTQC_FIRST_OUTDIR} ]; then mkdir {FASTQC_FIRST_OUTDIR}; else rm -r {FASTQC_FIRST_OUTDIR} && mkdir {FASTQC_FIRST_OUTDIR}; fi"
    	"&& fastqc {input} --outdir {FASTQC_FIRST_OUTDIR}"

rule cutadapt_first_read_clip:
	input:
		expand(config['sample_data_dir'] + "/{sample}.fastqsanger", sample=FIRST_READS)
	output:
		seq=expand(CUTADAPT_OUTDIR + "/{sample}.fastqsanger", sample=FIRST_READS),
		log=expand(CUTADAPT_OUTDIR + "/{sample}.txt", sample=FIRST_READS)
	run:
		shell("if [ ! -d {CUTADAPT_OUTDIR} ]; then mkdir {CUTADAPT_OUTDIR}; else rm -r {CUTADAPT_OUTDIR} && mkdir {CUTADAPT_OUTDIR}; fi")
		for i in range(0,len(input)):
			shell("cutadapt --format=fastq  --error-rate=0.1 --times=1 --overlap=5 " 
			"--front='CTTCCGATCTACAAGTT' --front='CTTCCGATCTTGGTCCT' " 
			"--output=" + output.seq[i] + " " + input[i] + " > " + output.log[i])

rule cutadapt_second_read_clip:
	input:
		expand(config['sample_data_dir'] + "/{sample}.fastqsanger", sample=SECOND_READS)
	output:
		seq=expand(CUTADAPT_OUTDIR + "/{sample}.fastqsanger", sample=SECOND_READS),
		log=expand(CUTADAPT_OUTDIR + "/{sample}.txt", sample=SECOND_READS)
	run:
		for i in range(0,len(input)):
			shell("cutadapt --format=fastq --error-rate=0.1 --times=1 --overlap=5 " 
		"--adapter='AACTTGTAGATCGGA' --adapter='AGGACCAAGATCGGA' --adapter='ACTTGTAGATCGGAA' --adapter='GGAAGAGCGTCGTGT' "
		"--adapter='GGACCAAGATCGGAA' --adapter='CTTGTAGATCGGAAG' --adapter='GACCAAGATCGGAAG' --adapter='TTGTAGATCGGAAGA' " 
		"--adapter='ACCAAGATCGGAAGA' --adapter='TGTAGATCGGAAGAG' --adapter='CCAAGATCGGAAGAG' --adapter='GTAGATCGGAAGAGC' " 
		"--adapter='CAAGATCGGAAGAGC' --adapter='TAGATCGGAAGAGCG' --adapter='AAGATCGGAAGAGCG' --adapter='AGATCGGAAGAGCGT' " 
		"--adapter='GATCGGAAGAGCGTC' --adapter='ATCGGAAGAGCGTCG' --adapter='TCGGAAGAGCGTCGT' --adapter='CGGAAGAGCGTCGTG' "      
		"--output=" + output.seq[i] + " " + input[i] + " > " + output.log[i])

rule trim_galore_clip:
	input:
		first_read=expand(CUTADAPT_OUTDIR + "/{sample}.fastqsanger", sample=FIRST_READS),
		second_read=expand(CUTADAPT_OUTDIR + "/{sample}.fastqsanger", sample=SECOND_READS)
	params:
		cut_n_bases=5
	output:
		first_read_out=expand(TRIMGALORE_OUTDIR + "/{sample}.fastqsanger", sample=FIRST_READS),
		second_read_out=expand(TRIMGALORE_OUTDIR + "/{sample}.fastqsanger", sample=SECOND_READS)
	run:
		shell("if [ ! -d {TRIMGALORE_OUTDIR} ]; then mkdir {TRIMGALORE_OUTDIR}; else rm -r {TRIMGALORE_OUTDIR} && mkdir {TRIMGALORE_OUTDIR}; fi")
		for i in range(0,len(input.first_read)):
			shell("trim_galore --phred33 --suppress_warn --three_prime_clip_R1 {params.cut_n_bases} --paired --dont_gzip "
			"--output_dir {TRIMGALORE_OUTDIR} " + input.first_read[i] + " " + input.second_read[i] + " " + 
			"&& mv " + output.first_read_out[i] + "_val_1.fq " + output.first_read_out[i] + " " + 
			"&& mv " + output.second_read_out[i] + "_val_2.fq " + output.second_read_out[i])

#############
## MAPPING ##
#############

# rule star:
# 	input:
# 		first_read=expand(TRIMGALORE_OUTDIR + "/{sample}.fastqsanger", sample=FIRST_READS),
# 		second_read=expand(TRIMGALORE_OUTDIR + "/{sample}.fastqsanger", sample=SECOND_READS)
# 	params:
# 		cut_n_bases=5
# 	resources:
# 		mem=100
# 	threads: 12
# 	output:
# 		expand(MAPPING_OUTDIR + "/.bam", sample=FIRST_READS)
# 	run:
# 		shell("if [ ! -d {MAPPING_OUTDIR} ]; then mkdir {MAPPING_OUTDIR}; else rm -r {MAPPING_OUTDIR} && mkdir {MAPPING_OUTDIR}; fi")
# 		for i in range(0,len(input.first_read)):
# 			shell("STAR --runThreadN {threads} --genomeLoad NoSharedMemory --genomeDir {REF_MAPPING} "   
# 			"--readFilesIn " + input.first_read[i] + " " + input.second_read[i] + " "  
# 			"--outSAMtype BAM SortedByCoordinate --outSAMattributes All --outSAMstrandField intronMotif --outFilterIntronMotifs None --outSAMunmapped None " 
# 			"--outSAMprimaryFlag OneBestScore --outSAMmapqUnique '255' --outFilterType Normal --outFilterMultimapScoreRange '1' --outFilterMultimapNmax '10' "
# 			"--outFilterMismatchNmax '10' --outFilterMismatchNoverLmax '0.3' --outFilterMismatchNoverReadLmax '1.0' --outFilterScoreMin '0' --outFilterScoreMinOverLread '0.66' " 
# 			"--outFilterMatchNmin '0' --outFilterMatchNminOverLread '0.66'   --seedSearchStartLmax '50' --seedSearchStartLmaxOverLread '1.0' --seedSearchLmax '0' " 
# 			"--seedMultimapNmax '10000' --seedPerReadNmax '1000' --seedPerWindowNmax '50' --seedNoneLociPerWindow '10'  --alignIntronMin '21' --alignIntronMax '0' " 
# 			"--alignMatesGapMax '0' --alignSJoverhangMin '5' --alignSJDBoverhangMin '3' --alignSplicedMateMapLmin '0' --alignSplicedMateMapLminOverLmate '0.66' " 
# 			"--alignWindowsPerReadNmax '10000' --alignTranscriptsPerWindowNmax '100' --alignTranscriptsPerReadNmax '10000' --alignEndsType Local  --twopassMode 'Basic' " 
# 			"--twopass1readsN '-1'   --limitBAMsortRAM '0' --limitOutSJoneRead '1000' --limitOutSJcollapsed '1000000' --limitSjdbInsertNsj '1000000' ")

########################
## POST-MAP FILTERING ##
########################

rule unique_reads_fitlering:
	input:
		expand(MAPPING_OUTDIR + "/{sample}.bam", sample=ALL_REPLICATES)
	output:
		expand(PRE_FOR_UMI_OUTDIR + "/{sample}_unique_reads_fitlering.bam", sample=ALL_REPLICATES)
	run:
		shell("if [ ! -d {PRE_FOR_UMI_OUTDIR} ]; then mkdir {PRE_FOR_UMI_OUTDIR}; else rm -r {PRE_FOR_UMI_OUTDIR} && mkdir {PRE_FOR_UMI_OUTDIR}; fi")
		for i in range(0,len(input)):
			shell("samtools view -h " +  input[i] + " | "
			"awk '$0 ~ /^@/{{ print }} ($2 == 163 || $2 == 147 || $2 == 83 || $2 == 99)&&$0!~/XS:i/{{ print }} ' | " +  
			"samtools view -bSh > " + output[i])

rule got_umis:
	input:
		expand(PRE_FOR_UMI_OUTDIR + "/{sample}_unique_reads_fitlering.bam", sample=ALL_REPLICATES)
	output:
		expand(PRE_FOR_UMI_OUTDIR + "/{sample}_got_umis.bam", sample=ALL_REPLICATES)
	run:
		for i in range(0,len(input)):
			shell("samtools view -h " +  input[i] + " | "
			""" awk -F "\t" 'BEGIN {{ OFS = FS }} {{ if ($0 ~ /^@/) {{ print $0; }} else {{split($1,str,":"); $1=str[2]":"str[3]":"str[4]":"str[5]":"str[6]":"str[7]":"str[8]"_"str[1]; print $0; }} }}' | """
			"samtools view -bSh > " + output[i])
 
###################
## DEDUPLICATION ##
###################

rule deduplication:
	input:
		expand(PRE_FOR_UMI_OUTDIR + "/{sample}_got_umis.bam", sample=ALL_REPLICATES)
	output:
		bam=expand(DEDUPLICAITON_OUTDIR + "/{sample}.bam", sample=ALL_REPLICATES),
		log=expand(DEDUPLICAITON_OUTDIR + "/{sample}_log.txt", sample=ALL_REPLICATES),
		sorted_bam=expand(DEDUPLICAITON_OUTDIR + "/{sample}_sorted.bam", sample=ALL_REPLICATES)
	run:
		shell("if [ ! -d {DEDUPLICAITON_OUTDIR} ]; then mkdir {DEDUPLICAITON_OUTDIR}; else rm -r {DEDUPLICAITON_OUTDIR} && mkdir {DEDUPLICAITON_OUTDIR}; fi")
		for i in range(0,len(input)):
			shell("samtools index " + input[i] + " "
				"&& umi_tools dedup --random-seed 0 --extract-umi-method read_id --method adjacency --edit-distance-threshold 1 --paired " 
				"--soft-clip-threshold 4 --subset 1.0 -I " + input[i] + " -S " + output.bam[i] + " -L " + output.log[i] + " "
				"&& samtools sort " + output.bam[i] + " > " +  output.sorted_bam[i] + " "
				"&& samtools index " + output.sorted_bam[i])

############################
## SECOND QUALITY CONTROL ##
############################

rule fastqc_second:
    input:
    	expand(DEDUPLICAITON_OUTDIR + "/{sample}.bam", sample=ALL_REPLICATES)
    output:
    	expand(FASTQC_SECOND_OUTDIR + "/{sample}.fastqsanger_fastqc.html", sample=ALL_REPLICATES),
    	expand(FASTQC_SECOND_OUTDIR + "/{sample}.fastqsanger_fastqc.zip", sample=ALL_REPLICATES)
    shell:
    	"if [ ! -d {FASTQC_SECOND_OUTDIR} ]; then mkdir {FASTQC_SECOND_OUTDIR}; else rm -r {FASTQC_SECOND_OUTDIR} && mkdir {FASTQC_SECOND_OUTDIR}; fi"
    	"&& fastqc {input} --outdir {FASTQC_SECOND_OUTDIR}"

#################
## PEAKCALLING ##
#################

rule peakachu:
    input:
    	clip=expand(DEDUPLICAITON_OUTDIR + "/{sample}_sorted.bam", sample=REPLICATES_CLIP),
    	control=expand(DEDUPLICAITON_OUTDIR + "/{sample}_sorted.bam", sample=REPLICATES_CONTROL)
    output:
    	peaks_tsv=PEAKCALLING_OUTDIR + "/peakachu.tsv",
    	peaks_gtf=PEAKCALLING_OUTDIR + "/peakachu.gtf",
    	blockbuster=PEAKCALLING_OUTDIR + "/blockbuster.bed"
    shell:
    	#"--max_proc "${GALAXY_SLOTS:-1}"
    	"if [ ! -d {PEAKCALLING_OUTDIR} ]; then mkdir {PEAKCALLING_OUTDIR}; else rm -r {PEAKCALLING_OUTDIR} && mkdir {PEAKCALLING_OUTDIR}; fi "
    	"&& peakachu adaptive --exp_libs {input.clip} --ctr_libs {input.control} "
    	"--paired_end --max_insert_size 50 --features '' --sub_features '' "
    	"--output_folder {PEAKCALLING_OUTDIR} --min_cluster_expr_frac 0.01 --min_block_overlap 0.5 --min_max_block_expr 0.1 "
    	"--norm_method deseq --mad_multiplier 0.0 --fc_cutoff 0.0 --padj_threshold 0.05 "
    	"&& head -q -n 1 {PEAKCALLING_OUTDIR}/peak_tables/*.csv > tmp.tsv "
    	"&& head -q -n 1 tmp.tsv > {output.peaks_tsv}"
    	"&& rm tmp.tsv"
    	"&& tail -n +2 -q {PEAKCALLING_OUTDIR}/peak_tables/*.csv >> {output.peaks_tsv} "
    	"&& cat {PEAKCALLING_OUTDIR}/peak_annotations/*.gff | awk '/peak/ {{ print $0 }}' > {output.peaks_gtf} "
    	"&& cat {PEAKCALLING_OUTDIR}/blockbuster_input/*bed > {output.blockbuster}"

##############################
## POST-PROCESSING OF PEAKS ##
##############################

rule peak_tsv_to_bed:
	input:
		PEAKCALLING_OUTDIR + "/peakachu.tsv"
	output:
		POSTPROCESSING_OUTDIR + "/peakachu.bed"
	shell:
		"if [ ! -d {POSTPROCESSING_OUTDIR} ]; then mkdir {POSTPROCESSING_OUTDIR}; else rm -r {POSTPROCESSING_OUTDIR} && mkdir {POSTPROCESSING_OUTDIR}; fi"
		"""&&  awk -F "\t" 'BEGIN {{ OFS = FS; print "Chrom", "Start","End","Name","Score","Strand" }} NR>1 {{ if ($3 < $4) {{ print $1,$3,$4,"clip_peak_"NR-1,$9,$5; }} else {{ print $1,$4,$3,"clip_peak_"NR-1,$9,$5; }} }}' {input} > {output} """

#####################
## MOTIF DETECTION ##
#####################

rule extract_genomic_DNA_dreme:
	input:
		POSTPROCESSING_OUTDIR + "/peakachu.bed"
	output:
		MOTIF_DETECTION_OUTDIR + "/peaks.fa"
	run:
		shell("if [ ! -d {MOTIF_DETECTION_OUTDIR} ]; then mkdir {MOTIF_DETECTION_OUTDIR}; else rm -r {MOTIF_DETECTION_OUTDIR} && mkdir {MOTIF_DETECTION_OUTDIR}; fi"
		"""&& awk 'NR>1 {{ print $1"\t"$2"\t"$3 }}' {input} > {MOTIF_DETECTION_OUTDIR}/peaks.bed3 """
		"&& python /home/flow/Documents/PycharmProjects/ClIP_QC/fetch_DNA_sequence.py -o {output} {MOTIF_DETECTION_OUTDIR}/peaks.bed3 {REF_MAPPING}/hg19_2.fa")

rule dreme:
	input: 
		MOTIF_DETECTION_OUTDIR + "/peaks.fa"
	output:
		MOTIF_DETECTION_OUTDIR + "/dreme/dreme.html"
	run:
		shell("python2 /home/flow/miniconda3/bin/dreme -p {input} -dna -s '1' -e 0.05 -m 100 -g 100 -mink 5 -maxk 20 "
			"-oc {MOTIF_DETECTION_OUTDIR}/dreme")

rule cross_random_subsampling_for_meme:
	input:
		POSTPROCESSING_OUTDIR + "/peakachu.bed"
	output:
		expand(MOTIF_DETECTION_OUTDIR + "/meme/{peak_file}.bed", peak_file=PEAK_FILES_FOR_MEME)
	params:
		n_peaks=1000,
		seed_value=123
	run:
		lines = open(str(input)).read().splitlines()
		header = lines[0]
		lines = lines[1:]
		if params.seed_value != -1:
			random.seed(params.seed_value)
		shell("if [ ! -d {MOTIF_DETECTION_OUTDIR}/meme ]; then mkdir {MOTIF_DETECTION_OUTDIR}/meme; else rm -r {MOTIF_DETECTION_OUTDIR}/meme && mkdir {MOTIF_DETECTION_OUTDIR}/meme; fi")
		for i in range(0, len(PEAK_FILES_FOR_MEME)):
			peak_file = "peak_random_sample_" + str(i+1) 
			complet_path_to_peak_file = MOTIF_DETECTION_OUTDIR + "/meme/" + peak_file
			shell("if [ ! -d {complet_path_to_peak_file} ]; then mkdir {complet_path_to_peak_file}; else rm -r {complet_path_to_peak_file} && mkdir {complet_path_to_peak_file}; fi")
			outfile = open(str(complet_path_to_peak_file + "/" + peak_file + ".bed"),"w") 
			outfile.write(header + "\n")
			choices = random.sample(lines, params.n_peaks)
			for l in choices[:-1]:			
				outfile.write(l + "\n")
			outfile.write(choices[-1])
			outfile.close()

rule extract_genomic_DNA_meme:
	input:
		expand(MOTIF_DETECTION_OUTDIR + "/meme/{peak_file}.bed", peak_file=PEAK_FILES_FOR_MEME)
	output:
		bed3=expand(MOTIF_DETECTION_OUTDIR + "/meme/{peak_file}.bed3", peak_file=PEAK_FILES_FOR_MEME),
		fa=expand(MOTIF_DETECTION_OUTDIR + "/meme/{peak_file}.fa", peak_file=PEAK_FILES_FOR_MEME)
	run:
		for i in range(0,len(input)):
			shell("""awk 'NR>1 {{ print $1"\t"$2"\t"$3 }}' """ + input[i] + " > " + output.bed3[i] + " "
			"&& python /home/flow/Documents/PycharmProjects/ClIP_QC/fetch_DNA_sequence.py -o " + output.fa[i] + " " + 
			output.bed3[i] + " {REF_MAPPING}/hg19_2.fa")

rule meme:
	input: 
		expand(MOTIF_DETECTION_OUTDIR + "/meme/{peak_file}.fa", peak_file=PEAK_FILES_FOR_MEME)
	output:
		output_dir=expand(MOTIF_DETECTION_OUTDIR + "/meme/{peak_file}_meme_output", peak_file=PEAK_FILES_FOR_MEME)
	run:
		for i in range(0, len(input)):
			shell("meme " + input[i] + " -dna -revcomp -pal -prior dirichlet -b 0.01 -spmap uni -spfuzz 0.5 -nmotifs 5 -mod zoops -wnsites 0.8 -minw 5 -maxw 20 "
			"-wg 11 -ws 1 -maxiter 50 -distance 0.001 -oc " + output.output_dir[i])

##################
## MOTIF SEARCH ##
##################

rule fimo_for_dreme_output:
	input: 
		MOTIF_DETECTION_OUTDIR + "/dreme/dreme.gff"
	output:
		html=MOTIF_SEARCH_OUTDIR + "/fimo_dreme/{peak_file}.html"
		interval=MOTIF_SEARCH_OUTDIR + "/fimo_dreme/{peak_file}.interval"
		txt=MOTIF_SEARCH_OUTDIR + "/fimo_dreme/{peak_file}.txt"
		xml=MOTIF_SEARCH_OUTDIR + "/fimo_dreme/{peak_file}.xml"
	run:
		shell("if [ ! -d {MOTIF_SEARCH_OUTDIR} ]; then mkdir {MOTIF_SEARCH_OUTDIR}; else rm -r {MOTIF_SEARCH_OUTDIR} && mkdir {MOTIF_SEARCH_OUTDIR}; fi"
		"&& if [ ! -d {MOTIF_SEARCH_OUTDIR}/fimo_dreme ]; then mkdir {MOTIF_SEARCH_OUTDIR}/fimo_dreme; else rm -r {MOTIF_SEARCH_OUTDIR}/fimo_dreme && mkdir {MOTIF_SEARCH_OUTDIR}/fimo_dreme; fi" 
		"&& fimo --input_motifs {input} --input_fasta {REF_MAPPING}/hg19_2.fa "
		"--options_type advanced --alpha '1.0'  --max_stored_scores '100000' --output_separate_motifs yes --motif_pseudo '0.1' "  
		"--qv_thresh --thresh 0.0001 --output_path {MOTIF_SEARCH_OUTDIR}/fimo_dreme "
		"--html_output {output.html} --interval_output {output.interval} --txt_output {output.txt} --xml_output {output.xml} --gff_output 'None' ")


# rule motif_search_for_meme_output:
# 	input: 
# 		expand(MOTIF_DETECTION_OUTDIR + "/meme/{peak_file}.fa", peak_file=PEAK_FILES_FOR_MEME)
# 	output:
# 		expand(MOTIF_DETECTION_OUTDIR + "/meme/{peak_file}.html", peak_file=PEAK_FILES_FOR_MEME)
# 	run:
# 		shell("if [ ! -d {MOTIF_SEARCH_OUTDIR} ]; then mkdir {MOTIF_SEARCH_OUTDIR}; else rm -r {MOTIF_SEARCH_OUTDIR} && mkdir {MOTIF_SEARCH_OUTDIR}; fi")
# 		for i in range(0, len(input)):
# 			shell(mkdir -p output && python '/usr/local/galaxy/shed_tools/toolshed.g2.bx.psu.edu/repos/iuc/meme_fimo/c470b36b592d/meme_fimo/fimo_wrapper.py'
# 			 --input_motifs '/data/dnb01/galaxy_db/files/006/541/dataset_6541255.dat' --input_fasta '/data/db/reference_genomes/hg19/seq/hg19.fa' 
# 			 --options_type advanced --alpha '1.0'  --max_stored_scores '100000' --output_separate_motifs yes --motif_pseudo '0.1'   
# 			 --qv_thresh --thresh 0.0001 --output_path '/data/dnb01/galaxy_db/job_working_directory/003/635/3635323/dataset_6532511_files' 
# 			 --html_output '/data/dnb01/galaxy_db/files/006/532/dataset_6532511.dat' --interval_output '/data/dnb01/galaxy_db/files/006/532/dataset_6532514.dat' 
# 			 --txt_output '/data/dnb01/galaxy_db/files/006/532/dataset_6532512.dat' --xml_output '/data/dnb01/galaxy_db/files/006/532/dataset_6532513.dat' 
# 			 --gff_output 'None')