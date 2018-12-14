
###################
## PREPROCESSING ##
###################

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
		first=RENAMING + "/{sample}_{replicate}_r1.fastqsanger"
	output:
		seq_first=CUTADAPT_OUTDIR + "/{sample}_{replicate}_r1.fastqsanger",
		log=CUTADAPT_OUTDIR + "/{sample}_{replicate}.txt"
	threads: 8 
	conda:
		"envs/cutadapt.yml"
	shell:
		"if [ ! -d {CUTADAPT_OUTDIR} ]; then mkdir {CUTADAPT_OUTDIR}; fi"
		"&& cutadapt -j {threads} {config[cutadapt]} " 
		"--output={output.seq_first} {input.first} > {output.log}"

# trim of low quality bases at the ned of the parclip data 
rule remove_tail:
	input:
		first=CUTADAPT_OUTDIR + "/{samples}_{replicate}_r1.fastqsanger"
	output:
		first=REMOVE_TAIL_OUTDIR + "/{samples}_{replicate}_r1_trimmed.fastqsanger"
	threads: 2
	conda:
		"envs/bctools.yml"
	shell:
		"if [ ! -d {REMOVE_TAIL_OUTDIR} ]; then mkdir {REMOVE_TAIL_OUTDIR}; fi"
		"&& python {config[bctools]}/remove_tail.py {input.first} {config[remove_tail]} > {output.first}"

rule fastqc_after_adapter_removal:
    input:
    	REMOVE_TAIL_OUTDIR + "/{sample}_{replicate}_{pair}_trimmed.fastqsanger"
    output:
    	FASTQC_ADAPT_OUTDIR + "/{sample}_{replicate}_{pair}_trimmed.fastqsanger_fastqc.html",
    	FASTQC_ADAPT_OUTDIR + "/{sample}_{replicate}_{pair}_trimmed.fastqsanger_fastqc.zip"
    threads: 2
    conda:
    	"envs/fastqc.yml"
    shell:
    	"if [ ! -d {FASTQC_ADAPT_OUTDIR} ]; then mkdir {FASTQC_ADAPT_OUTDIR}; fi"
    	"&& fastqc {input} --outdir {FASTQC_ADAPT_OUTDIR}"