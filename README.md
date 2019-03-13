# SalamiSnake_CLIP

## Background
Snakemake pipeline for FLASH, eCLIP, ChIP-Seq and PARCLIP data. Changes needs to be done in order to adapt to your data.
The file **Snakefile** includes evey step up to the point of the quality check of the mapping. The **cluster.yml** can be amended to give some tools more threads and memory. The **config.yml** need to be adapted for your data. The **.qsub** scripts can then be called with ``qsub *.qsub`` to queue it into SGE (Son of Grid Engine).
Conda environments are placed into **envs/**. Other additional scirpts can be placed into **scripts/**. Rules that are imported into the **Snakefile** are located in **rules/**.

## Input Data
Data needs to be provided in **.fastqsanger** and it must be **de-multiplexed**.

## Config
To apply the pipleine to your data you need to specifiy the following things inseide the **config.yml**. <br/>
- Define the names of your flash and control experiments:
  - first define a key-name for each replicate (this is your choice)
  - then insert the name of your fastqsanger files **without** the **.fastqsanger** extension
  - files need to be correctly sorreted into **r1 = forward** and **r2 = reverse** order
  - **NOTE!: The key-name for the files needs to be consistent, i.e., the pipeline do not allow for something like: EBF1_rep1 and HR_EBF1_rep2. It need to be EBF1_rep1 and EBF1_rep2.**
- Define the path where the data is located in **sample_data_dir:**
- Change output directory names if you want.

## Qsub Files
You need to define the following things in the two **qsub** files:
- Define the output and error log path (**-o** and **-e**)
- Miniconda path **MINICONDA=**
- Email to message when job fails **-M** and !! also in the snakemake command !!

## Cluster.yml
You can specifiy the resources of various tools. You also need to define the output and error log path in **out** and **err** for each tool.
- **cores** = number of threads
- **mem** = memory

## Conda environments
You can place conda environemnts inside envs. Take a look into the [Snakemake guide for conda environments](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html).
To automatically install the tool itself in the conda environment, place the name of the tool inside the **dependencies**. The dependencies for the tools can be found in the [conda github for recipies](https://github.com/bioconda/bioconda-recipes), or in the diretory of each tool in your miniconda. <br/>

Tools can be linked to a conda envs inside the Snakemake pipline with **conda:**.

## Snakefile
The Snakefile is mainly responsible for the pre-processing of your data. It can be executed with the **snakemake_qsub.sh** script. <br/>

You need to define the following things:
- Path to your reference genome **REF_GENOME_DIR**. Following data is needed:
  - **genome** in fasta format
  - **transcriptome** annotation in .gtf format
  - **chromosome sizes** in .txt format
- Path to your config.yml **configfile**
- Adapter sequences inside the **Cutadapt** calls <br/>

Furthermore you need to keep in mind:
- Activate **got_umis**  and **barcode_header_removal** to extract UMIs and barcodes from your reads if not already done.
- Activate **estimate_insert_size** to estimate the insert size if you are using PEAKachu for peak calling. It will give you information about this parameter.
- Deactivate **filter_out_unlocalized_regions_for_later_genome_versions** if you are using versions of genomes which do not include these regions.

## SnakefilePeakcalling.py
The SnakefilePeakcalling.py file is mainly responsible for the peak calling. It can be executed with the **snakemake_peakcalling_qsub.sh** script. **The data needs to be pre-processed with the Snakefile beforehand**.

You need to define the following things:
- Path to your reference genome **REF_GENOME_DIR**. Following data is needed:
  - **genome** in fasta format
  - **transcriptome** annotation in .gtf format
  - **chromosome sizes** in .txt format
- Path to your config.yml **configfile**
- When using PureCLIP only one side of the read pair is used. Based on the protocol change the samtools flag inside **mate_reads_fitlering**. Normally the crosslinking site is located on the forward (first) read in a standard FLASH protocol, thus the flag is **0x0040**. Change it to **0x0080** for the reverse (second) read.

You can additionally change the following things:
- Either use PureCLIP or PEAKachu for peak calling. When using PureCLIP activate the following things:
  - **mate_reads_fitlering**
  - **pureclip**
  - **intersect_binding_regions_with_peaks** for PureCLIP

## Further Remarks
Keep the following thins in mind:
- Make yourself familiar with the tools.
- **Cutadapt** is very important to change
- **FASTQC** protocols are very important to look at.
- **PEAKachu** `Pairwise Replicates` set it on `true` if all replicates from the signal and control come from the same sample. It models batch effects. `true` = (~ samples + conditions); `false` = (~ conditions)
- **PureCLIP** works best with only one side of the paired reads.
- **PureCLIP** works best if you only estimate the HMM parameters for three to four chromosomes. Else you get problewms with runtime and memory.
- **PureCLIP** binding regions are often only one nucleotide long apart form the crosslinking site file. Further analysis require **SlopBED** to extend these regions.
