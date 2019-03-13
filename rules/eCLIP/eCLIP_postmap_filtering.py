import random
import math
import os 

########################
## POST-MAP FILTERING ##
########################

# -F 0x100 do not include secondary alignments 
# 0!~/XS:i/{{ print }} for bowtie flag of nonunqiue mapped reads
rule unique_reads_fitlering:
	input:
		MAPPING_OUTDIR + "/{sample}_{replicate}.bam"
	output:
		PRE_FOR_UMI_OUTDIR + "/{sample}_{replicate}_unique_reads_fitlering.bam"
	threads: 2
	conda:
		config["conda_envs"] + "/samtools.yml"
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
		config["conda_envs"] + "/samtools.yml"
	shell:
		"samtools view -h {input} "
		"""| awk -F "\t" 'BEGIN {{ OFS = FS }} {{ if ($0 ~ /^@/) {{print $0;}} else {{ if ($3 ~/^chr/) {{print $0;}} }} }}' """
		"| samtools view -bSh > {output}"
		"&& samtools index {output}"
 