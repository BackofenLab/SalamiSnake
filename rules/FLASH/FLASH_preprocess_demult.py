import random
import math
import os 

rule fastqc_beginning:
    input:
    	RENAMING + "/{sample}_{replicate}_{pair}.fastqsanger"
    output:
    	FASTQC_BEG_OUTDIR + "/{sample}_{replicate}_{pair}.fastqsanger_fastqc.html",
    	FASTQC_BEG_OUTDIR + "/{sample}_{replicate}_{pair}.fastqsanger_fastqc.zip"
    threads: 2
    conda:
    	config["conda_envs"] + "/fastqc.yml"
    shell:
    	"if [ ! -d {FASTQC_BEG_OUTDIR} ]; then mkdir {FASTQC_BEG_OUTDIR}; fi"
    	"&& fastqc {input} --outdir {FASTQC_BEG_OUTDIR}"

rule cutadapt_first_read_clip:
	input:
		first=RENAMING + "/{sample}_{replicate}_r1.fastqsanger",
		second=RENAMING + "/{sample}_{replicate}_r2.fastqsanger"
	output:
		seq_first=CUTADAPT_OUTDIR + "/{sample}_{replicate}_r1.fastqsanger",
		seq_second=CUTADAPT_OUTDIR + "/{sample}_{replicate}_r2.fastqsanger",
		log=CUTADAPT_OUTDIR + "/{sample}_{replicate}.txt"
	threads: 2 
	conda:
		config["conda_envs"] + "/cutadapt.yml"
	shell:
		"if [ ! -d {CUTADAPT_OUTDIR} ]; then mkdir {CUTADAPT_OUTDIR}; fi"
		"&& echo {config[cutadapt]} >> {file_tool_params}"
		"&& cutadapt -j {threads} {config[cutadapt]} " 
		"--paired-output={output.seq_second} --output={output.seq_first} {input.first} {input.second} > {output.log}"

rule remove_tail:
	input:
		first=CUTADAPT_OUTDIR + "/{sample}_{replicate}_r1.fastqsanger",
		second=CUTADAPT_OUTDIR + "/{sample}_{replicate}_r2.fastqsanger"
	output:
		first=REMOVE_TAIL_OUTDIR + "/{sample}_{replicate}_r1_trimmed.fastqsanger",
		second=REMOVE_TAIL_OUTDIR + "/{sample}_{replicate}_r2_trimmed.fastqsanger"
	threads: 2
	conda:
		config["conda_envs"] + "/bctools.yml"
	shell:
		"if [ ! -d {REMOVE_TAIL_OUTDIR} ]; then mkdir {REMOVE_TAIL_OUTDIR}; fi"
		"&& echo {config[remove_tail]} >> {file_tool_params}"
		"&& python {config[bctools]}/remove_tail.py {input.first} {config[remove_tail]} > {output.first}"
		"&& cp {input.second} {output.second}"

rule fastqc_after_adapter_removal:
    input:
    	REMOVE_TAIL_OUTDIR + "/{sample}_{replicate}_{pair}_trimmed.fastqsanger"
    output:
    	FASTQC_ADAPT_OUTDIR + "/{sample}_{replicate}_{pair}_trimmed.fastqsanger_fastqc.html",
    	FASTQC_ADAPT_OUTDIR + "/{sample}_{replicate}_{pair}_trimmed.fastqsanger_fastqc.zip"
    threads: 2
    conda:
    	config["conda_envs"] + "/fastqc.yml"
    shell:
    	"if [ ! -d {FASTQC_ADAPT_OUTDIR} ]; then mkdir {FASTQC_ADAPT_OUTDIR}; fi"
    	"&& fastqc {input} --outdir {FASTQC_ADAPT_OUTDIR}"

# Activate if you need to remove UMIs and barcodes
rule got_umis:
	input:
		first=REMOVE_TAIL_OUTDIR + "/{sample}_{replicate}_r1_trimmed.fastqsanger",
		second=REMOVE_TAIL_OUTDIR + "/{sample}_{replicate}_r2_trimmed.fastqsanger"
	output:
		first=PRE_FOR_UMI_OUTDIR + "/{sample}_{replicate}_r1_trimmed.fastqsanger",
		second=PRE_FOR_UMI_OUTDIR + "/{sample}_{replicate}_r2_trimmed.fastqsanger",
		log=PRE_FOR_UMI_OUTDIR + "/{sample}_{replicate}_log.txt"
	threads: 2
	conda:
		config["conda_envs"] + "/umi.yml"
	shell:
		"if [ ! -d {PRE_FOR_UMI_OUTDIR} ]; then mkdir {PRE_FOR_UMI_OUTDIR}; fi "
		"&& echo {config[got_umis_1]} >> {file_tool_params}"
		"&& umi_tools extract {config[got_umis_1]} -I {input.second} -S {output.second} --read2-in {input.first} --read2-out {output.first} -L {output.log}"

# Convert Binary barcodes nulceotides
rule create_binary_barcode:
	input:
		PRE_FOR_UMI_OUTDIR + "/{sample}_{replicate}_r2_trimmed.fastqsanger"
	output:
		PRE_FOR_UMI_OUTDIR + "/{sample}_{replicate}_r2_trimmed_bbconverted.fastqsanger"
	threads: 2
	shell: 
		"""cat {input} | paste - - - - | awk -v FS=" " -v OFS=" " '{{bc=substr($3,0,2); rest=substr($3,3,length($3)); gsub("[AG]","A",bc); gsub("[TCU]","C",bc); $3=bc rest; print $1,$2"\t"$3"\t"$4"\t"$5}}' | tr -s "\t" "\n" > {output}"""

# Please use FORCE=true to allow overwrite.
rule demultiplex:
	input:
		first=PRE_FOR_UMI_OUTDIR + "/{sample}_{replicate}_r1_trimmed.fastqsanger",
		second=PRE_FOR_UMI_OUTDIR + "/{sample}_{replicate}_r2_trimmed_bbconverted.fastqsanger",
		barcodes=NEW_BARCODE_FILE
	output:
		diag=DEMULTI_OUTDIR + "/{sample}_{replicate}_diag.log",
		files=[DEMULTI_OUTDIR + "/{sample}_{replicate}_" + x for x in BARCODE_FILES]
	threads: 4
	conda:
		config["conda_envs"] + "/je_demultiplex.yml"
	params:
		outdir=DEMULTI_OUTDIR  
	shell:
		"if [ ! -d {DEMULTI_OUTDIR} ]; then mkdir {DEMULTI_OUTDIR}; fi "
		"&& echo {config[demultiplex]} >> {file_tool_params}"
		"&& je demultiplex FORCE=true F1={input.first} F2={input.second} BARCODE_FILE={input.barcodes} {config[demultiplex]} "
		"OUTPUT_DIR={params.outdir} BARCODE_DIAG_FILE={output.diag} METRICS_FILE_NAME={params.outdir}/metric.txt"

rule renaming_demultiplex:
	input:
		[DEMULTI_OUTDIR + "/" + MULTIPLEX_SAMPLE_NAME + "_rep1_" + x for x in BARCODE_FILES]
	output:
		BARCODE_NEWFILES
	run:
		for i in range(0, len(BARCODE_FILES)): 
			shell("mv " + input[i] + " " + output[i])