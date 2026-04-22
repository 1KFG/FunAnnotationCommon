#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 2 --mem 4G --out logs/antismash_nf.log -J antismash_nf -t 7-00:00:00

module load nextflow

mkdir -p logs

nextflow run pipeline/nextflow/antismash.nf \
    -profile slurm \
    -resume \
    "$@"
