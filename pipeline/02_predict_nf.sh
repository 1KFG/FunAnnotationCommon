#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 2 --mem 4G --out logs/funannotate_nf.log -J funannotate_nf -t 7-00:00:00

module load nextflow

mkdir -p logs

nextflow run pipeline/nextflow/funannotate.nf \
    -profile slurm \
    -resume \
    "$@"
