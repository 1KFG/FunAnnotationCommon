#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 2 -p batch --mem 4G --out logs/funannotate_nf.log -J funannotate_nf -t 14-00:00:00

# use local one
#module load nextflow

mkdir -p logs

nextflow run pipeline/nextflow/funannotate.nf \
    -profile slurm,funannotate \
    --suppress suppress.txt
    -resume \
    "$@"
