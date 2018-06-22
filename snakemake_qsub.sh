#!/bin/bash
#$ -N SalamiSnake
#$ -cwd
#$ -pe smp 1
#$ -l h_vmem=4G
#$ -o /scratch/bi01/heylf/SalamiSnake/
#$ -e /scratch/bi01/heylf/SalamiSnake/
#$ -M heylf@informatik.uni-freiburg.de
#$ -m a
#$ -R y
export PATH="/home/heylf/miniconda2/bin/:$PATH"
source activate snakemake
snakemake --use-conda -j 20 --cluster-config cluster.yml --cluster "qsub -R y -N {cluster.jobname} -cwd -pe {cluster.cores} -l {cluster.mem} -o {cluster.out} -e {cluster.err} -M heylf@informatik.uni-freiburg.de -m a" --latency-wait 60
