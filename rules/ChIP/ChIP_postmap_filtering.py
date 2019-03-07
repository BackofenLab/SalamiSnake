import random
import math
import os 

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
		config["conda_envs"] + "/samtools.yml"
	shell:
		"if [ ! -d {PRE_FOR_UMI_OUTDIR} ]; then mkdir {PRE_FOR_UMI_OUTDIR}; fi"
		"&& samtools view -h -F 0x100 {input} "
		"| awk '$0 ~ /^@/{{ print }} ($2 == 163 || $2 == 147 || $2 == 83 || $2 == 99)&&$0!~/XS:i/{{ print }} ' "  
		"| samtools view -bSh > {output}"

rule filter_out_unlocalized_regions_for_later_genome_versions:
	input:
		PRE_FOR_UMI_OUTDIR + "/{sample}_{replicate}_unique_reads_fitlering.bam"
	output:
		bam=PRE_FOR_UMI_OUTDIR + "/{sample}_{replicate}_got_umis_unlocalized_check.bam",
		bai=PRE_FOR_UMI_OUTDIR + "/{sample}_{replicate}_got_umis_unlocalized_check.bam.bai",
		bam_pos=PRE_FOR_UMI_OUTDIR + "/{sample}_{replicate}_got_umis_unlocalized_check_pos.bam",
		bai_pos=PRE_FOR_UMI_OUTDIR + "/{sample}_{replicate}_got_umis_unlocalized_check_pos.bam.bai",
		bam_neg=PRE_FOR_UMI_OUTDIR + "/{sample}_{replicate}_got_umis_unlocalized_check_neg.bam",
		bai_neg=PRE_FOR_UMI_OUTDIR + "/{sample}_{replicate}_got_umis_unlocalized_check_neg.bam.bai"
	threads: 2
	conda:
		config["conda_envs"] + "/samtools.yml"
	shell:
		"samtools view -h {input} "
		"""| awk -F "\t" 'BEGIN {{ OFS = FS }} {{ if ($0 ~ /^@/) {{print $0;}} else {{ if ($3 ~/^chr/) {{print $0;}} }} }}' """
		"| samtools view -bSh > {output.bam}"
		"&& samtools index {output.bam}"
		"&& samtools view -h -F 0x10 {output.bam} | samtools view -bSh > {output.bam_pos}"
		"&& samtools view -h -f 0x10 {output.bam} | samtools view -bSh > {output.bam_neg}"
		"&& samtools index {output.bam_pos}"
		"&& samtools index {output.bam_neg}"
 