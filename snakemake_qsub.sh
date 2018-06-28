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
MINICONDA="/scratch/bi01/heylf/miniconda3/bin/:$PATH"
export PATH=$MINICONDA
source activate snakemake
snakemake --use-conda -j 20 --directory ${PWD} --cluster-config cluster.yml --cluster "qsub -v PATH=$MINICONDA -R y -N {cluster.jobname} -cwd -pe {cluster.cores} -l {cluster.mem} -o {cluster.out} -e {cluster.err} -M heylf@informatik.uni-freiburg.de -m a" --latency-wait 60
