import random
import math
import itertools as iter
import os 
from snakemake.utils import validate, min_version
min_version("5.3.0")

SRC_PATH = os.getcwd()

configfile: "config.yml"

REF_GENOME_DIR = config['ref_genome_dir']
GENOME_FASTA = REF_GENOME_DIR + config['genome_fasta']
GENOME_2BIT = REF_GENOME_DIR + config['genome_2bit']
GENOME_GTF = REF_GENOME_DIR + config['genome_gtf']
GENOME_SIZES = REF_GENOME_DIR + config['genome_sizes']

PROTOCOL = config["protocol"]
paired = config["paired"]
control = config["control"]
demultiplexed = config["demultiplexed"]
finished_demultiplexed = config["finished_demultiplexed"]

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

def name_generation_samples(path, sample_name, replicate_names, pair_names, suffix):
	list_new_names = ["none"] * len(replicate_names) * len(pair_names)
	i = 0 
	for j in replicate_names:
		for k in pair_names:
			list_new_names[i] = path + "/" + sample_name + "_" + j + "_" + k + suffix
			i += 1
	return(list_new_names)

##############################
## WRITE TOOL PARAMS TO LOG ##
##############################

file_tool_params = open(config['sample_data_dir'] + "/tool_params.txt", "w")
for i in config['all_tool_params']:
	file_tool_params.write(i + "\n")
file_tool_params.close()

######################
## PATH DEFINITIONS ##
######################

RENAMING = config['sample_data_dir'] + "/" + config['renaming_outdir']
FASTQC_BEG_OUTDIR = config['sample_data_dir'] + "/" + config['fastqc_beg_outdir']
FASTQC_ADAPT_OUTDIR = config['sample_data_dir'] + "/" + config['fastqc_adapt_outdir']
CUTADAPT_OUTDIR = config['sample_data_dir'] + "/" + config['cutadapt_outdir']
REMOVE_TAIL_OUTDIR = config['sample_data_dir'] + "/" + config['remove_tail_outdir']
DEMULTI_OUTDIR = config['sample_data_dir'] + "/" + config['demulti_outdir']

MAPPING_OUTDIR = config['sample_data_dir'] + "/" + config['mapping_outdir']
MAPPING_QUALITY_OUTDIR = config['sample_data_dir'] + "/" + config['mapping_quality_outdir']
PRE_FOR_UMI_OUTDIR = config['sample_data_dir'] + "/" + config['preprocessing_for_umi_outdir']
DEDUPLICAITON_OUTDIR = config['sample_data_dir'] + "/" + config['deduplication_outdir']
FASTQC_BEFORE_DEDUP_OUTDIR = config['sample_data_dir'] + "/" + config['fastqc_before_dedup_outdir']
FASTQC_DEDUP_OUTDIR = config['sample_data_dir'] + "/" + config['fastqc_dedup_outdir']
COVERAGE_OUTDIR = config['sample_data_dir'] + "/" + config['coverage_outdir']

PEAKCALLING_OUTDIR = config['sample_data_dir'] + "/" + config['peakcalling_outdir']
ANNOTATION_PEAKS_OUTDIR = config['sample_data_dir'] + "/" + config['annotation_peaks_outdir']
HTSEQ_COUNT_OUTDIR = config['sample_data_dir'] + "/" + config['htseq_outdir']
DEDUPLICAITON_OUTDIR = config['sample_data_dir'] + "/" + config['deduplication_outdir']
MATEFILTER_OUTDIR = config['sample_data_dir'] + "/" + config['matefilter_outdir']

MULTIQC_OUTDIR = config['sample_data_dir'] + "/" + config['multiqc_outdir']

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

PAIR = []
if ( paired == "yes" ):
	PAIR = ["r1", "r2"]
else:
	PAIR = ["r1"]

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

#### For Demultiplexing ####
DEMULTIPLEX_SAMPLES = [config["multiplex_setup"]["r1"],  config["multiplex_setup"]["r2"]]
MULTIPLEX_SAMPLE_NAME = "multiplexed"
BARCODE_FILES = []
BARCODE_NEWFILES = []

barcode_sample_dict = {}

# generate buffer names 
if ( demultiplexed == "no" ):
	print("[NOTE] DEMULTIPLEX: " + demultiplexed)
	print(DEMULTIPLEX_SAMPLES)

	barcode_file = open(config["barcodes"],"r")
	new_barcode_file = open(config["sample_data_dir"] + "/new_barcode_file.txt","w")

	for line in barcode_file:
		sample = line.split("\t")[0]
		barcode =  line.split("\t")[1]
		if sample in barcode_sample_dict:
			sys.exit("[ERROR] Barcode Error, Sample with two barcodes.")
		else:
			barcode_sample_dict[sample] = str.strip(barcode)

		new_barcode_file.write(MULTIPLEX_SAMPLE_NAME + "_rep1_" + sample + "\t" + str.strip(barcode) + "\n")
	new_barcode_file.close()

	for i in FIRST_READS:
		if ( i in barcode_sample_dict ):
			# BARCODE_FILES.append(DEMULTI_OUTDIR + "/" + i + "_" + barcode_sample_dict[i] + "_1.txt")
			# BARCODE_FILES.append(DEMULTI_OUTDIR + "/" + i + "_" + barcode_sample_dict[i] + "_2.txt")
			BARCODE_FILES.append("/" + i + "_" + barcode_sample_dict[i] + "_1.txt")
			BARCODE_FILES.append("/" + i + "_" + barcode_sample_dict[i] + "_2.txt")
		else:
			print(i)
			sys.exit("[ERROR] Sample name in config is not coherent with sample name in the barcode file.")

	print(BARCODE_FILES)
	BARCODE_NEWFILES = name_generation_samples(DEMULTI_OUTDIR, SAMPLES[0], REP_NAME_CLIP, PAIR, "_trimmed.fastqsanger") 
	if ( control == "yes" ):
		BARCODE_NEWFILES = BARCODE_NEWFILES + name_generation_samples(DEMULTI_OUTDIR, SAMPLES[1], REP_NAME_CONTROL, PAIR, "_trimmed.fastqsanger")
	print(BARCODE_NEWFILES)

###########
## RULES ##
###########

if ( PROTOCOL == "FLASH" ):
	if( demultiplexed == "yes" ):
		# Preprocessing (no demultiplexing)
		include: config["rules"] + "/FLASH/FLASH_preprocess.py"
	else:
		# Preprocessing (with demultiplexing)
		include: config["rules"] + "/FLASH/FLASH_preprocess_demult.py"
	# Mapping
	include: config["rules"] + "/Main/mapping.py"
	# Postmap Filtering
	include: config["rules"] + "/FLASH/FLASH_postmap_filtering.py"
	# Deduplication
	include: config["rules"] + "/FLASH/FLASH_deduplication.py"
	# Mapping Quality
	include: config["rules"] + "/FLASH/FLASH_mapping_quality.py"
	# Coverage
	include: config["rules"] + "/Main/coverage.py"
	# Peakcalling
	include: config["rules"] + "/FLASH/FLASH_peakcalling.py"
	# # Peak Annotation
	# #include: config["rules"] + "/FLASH/FLASH.misc.py"
elif ( PROTOCOL == "PARCLIP" ):
	# Preprocessing (without demultiplexing)
	include: config["rules"] + "/PARCLIP/PARCLIP_preprocessing_singleend.py"
	# Mapping
	include: config["rules"] + "/Main/mapping_singleend.py"
	# Postmap Filtering
	include: config["rules"] + "/PARCLIP/PARCLIP_postmap_filtering.py"
	# Mapping Quality
	include: config["rules"] + "/PARCLIP/PARCLIP_mapping_quality.py"
	# Coverage
	include: config["rules"] + "/PARCLIP/PARCLIP_coverage.py"
	# Peakcalling
	include: config["rules"] + "/PARCLIP/PARCLIP_peakcalling.py"
	# Motif Detection
	include: config["rules"] + "/PARCLIP/PARCLIP_motif_detection.py"
else:
	sys.exit("[ERROR] Protocol not provided yet.")

############################
## INVOCATION (RULE ALL) ##
############################

ALL_NEW_FILE_NAMES = []

if ( control == "yes" ):
	if ( demultiplexed == "yes" ):
		rule all:
			input: 
				expand(FASTQC_BEG_OUTDIR + "/{sample}_{replicate}_{pair}.fastqsanger_fastqc.html", sample=SAMPLES[0], replicate=REP_NAME_CLIP, pair=PAIR),
				
		ALL_NEW_FILE_NAMES = name_generation_samples(RENAMING, SAMPLES[0], REP_NAME_CLIP, PAIR, ".fastqsanger") + name_generation_samples(RENAMING, SAMPLES[1], REP_NAME_CONTROL, PAIR, ".fastqsanger")
	else:
		# if ( finished_demultiplexed == "no" ):
		# 	rule all:
		# 		input: 
		# 			expand(DEMULTI_OUTDIR + "/{sample}_{replicate}_diag.log", sample=MULTIPLEX_SAMPLE_NAME, replicate="rep1")
		# else:
		# 	rule all:
		# 		input: 
		# 			BARCODE_NEWFILES,
		# 			expand(MAPPING_OUTDIR + "/{sample}_{replicate}.bam", sample=SAMPLES[0], replicate=REP_NAME_CLIP),
		# 			expand(MAPPING_OUTDIR + "/{sample}_{replicate}.bam", sample=SAMPLES[1], replicate=REP_NAME_CONTROL),
		# 			MULTIQC_OUTDIR + "/multiqc_report.html",
		# 			expand(PEAKCALLING_OUTDIR + "/{sample_exp}_{replicate_exp}_{sample_ctl}_{replicate_ctl}_peaks_extended.bed",
		# 			sample_exp=SAMPLES[0], replicate_exp=REP_NAME_CLIP, sample_ctl=SAMPLES[1], replicate_ctl=REP_NAME_CONTROL),
		# 			expand(COVERAGE_OUTDIR + "/bigwig/{sample}_{replicate}_alignment_ends_coverage_{type}_strand.bigwig", sample=SAMPLES[0], replicate=REP_NAME_CLIP, type=["pos", "neg", "both"]),
		# 			expand(COVERAGE_OUTDIR + "/bigwig/{sample}_{replicate}_alignment_ends_coverage_{type}_strand.bigwig", sample=SAMPLES[1], replicate=REP_NAME_CONTROL, type=["pos", "neg", "both"]),
		
		rule all:
			input: 
				expand(DEMULTI_OUTDIR + "/{sample}_{replicate}_diag.log", sample=MULTIPLEX_SAMPLE_NAME, replicate="rep1"),
				BARCODE_NEWFILES,
				expand(MAPPING_OUTDIR + "/{sample}_{replicate}.bam", sample=SAMPLES[0], replicate=REP_NAME_CLIP),
				expand(MAPPING_OUTDIR + "/{sample}_{replicate}.bam", sample=SAMPLES[1], replicate=REP_NAME_CONTROL),
				expand(MAPPING_OUTDIR + "/{sample}_{replicate}.bam.bai", sample=SAMPLES[0], replicate=REP_NAME_CLIP),
				expand(MAPPING_OUTDIR + "/{sample}_{replicate}.bam.bai", sample=SAMPLES[1], replicate=REP_NAME_CONTROL),
				# MAPPING_QUALITY_OUTDIR + "/fingerprint_plot.png",
				#MAPPING_QUALITY_OUTDIR + "/correlating_bam_files_plot.png",
				MULTIQC_OUTDIR + "/multiqc_report.html",
				expand(PEAKCALLING_OUTDIR + "/{sample_exp}_{replicate_exp}_{sample_ctl}_{replicate_ctl}_peaks_extended.bed",
				sample_exp=SAMPLES[0], replicate_exp=REP_NAME_CLIP, sample_ctl=SAMPLES[1], replicate_ctl=REP_NAME_CONTROL),
				expand(COVERAGE_OUTDIR + "/bigwig/{sample}_{replicate}_alignment_ends_coverage_{type}_strand.bigwig", sample=SAMPLES[0], replicate=REP_NAME_CLIP, type=["pos", "neg", "both"]),
				expand(COVERAGE_OUTDIR + "/bigwig/{sample}_{replicate}_alignment_ends_coverage_{type}_strand.bigwig", sample=SAMPLES[1], replicate=REP_NAME_CONTROL, type=["pos", "neg", "both"]),
	
		ALL_SAMPLES = DEMULTIPLEX_SAMPLES
		ALL_NEW_FILE_NAMES = name_generation_samples(RENAMING, MULTIPLEX_SAMPLE_NAME, ["rep1"], PAIR, ".fastqsanger")

else:
	rule all:
		input: 
			expand(MAPPING_OUTDIR + "/{sample}_{replicate}.bam", sample=SAMPLES[0], replicate=REP_NAME_CLIP),
			expand(MAPPING_OUTDIR + "/{sample}_{replicate}.bam.bai", sample=SAMPLES[0], replicate=REP_NAME_CLIP),
			MULTIQC_OUTDIR + "/multiqc_report.html",
			MAPPING_QUALITY_OUTDIR + "/fingerprint_plot.png",
			MAPPING_QUALITY_OUTDIR + "/correlating_bam_files_plot.png",
			expand(PEAKCALLING_OUTDIR + "/{sample}_{replicate}_peaks_extended.bed", sample=SAMPLES[0], replicate=REP_NAME_CLIP),
			expand(COVERAGE_OUTDIR + "/bigwig/{sample}_{replicate}_alignment_ends_coverage_{type}_strand.bigwig", sample=SAMPLES[0], replicate=REP_NAME_CLIP, type=["pos", "neg", "both"])#,
			# MOTIF_DETECTION_OUTDIR + "/meme_chip/meme-chip.html"

	ALL_NEW_FILE_NAMES = name_generation_samples(RENAMING, SAMPLES[0], REP_NAME_CLIP, PAIR, ".fastqsanger")

print("[NOTE] New names:")
print(ALL_NEW_FILE_NAMES)

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