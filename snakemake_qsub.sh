#!/bin/bash
#$ -N SalamiSnake
#$ -cwd
#$ -pe smp 1
#$ -l h_vmem=4G
#$ -o /scratch/bi03/heylf/Snakemake_Cluster_Logs/
#$ -e /scratch/bi03/heylf/Snakemake_Cluster_Logs/
#$ -M heylf@informatik.uni-freiburg.de
#$ -m a
#$ -R y
MINICONDA="/scratch/bi01/heylf/miniconda3/bin/:$PATH"
export PATH=$MINICONDA
source activate snakemake5.3
snakemake -s Snakefile.py --use-conda -j 20 --directory ${PWD} --cluster-config cluster.yml --cluster "qsub -v PATH=$MINICONDA -R y -N {cluster.jobname} -cwd -pe {cluster.cores} -l {cluster.mem} -o {cluster.out} -e {cluster.err} -M heylf@informatik.uni-freiburg.de -m a" --latency-wait 60
