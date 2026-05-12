#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 2 -p batch --mem 4G --out logs/interproscan6_nf.log -J interproscan6_nf -t 14-00:00:00

# Run interproscan6 (ebi-pf-team/interproscan6) on all predicted proteomes.
#
# Prerequisites (run once on head node before submitting this job):
#   nextflow pull ebi-pf-team/interproscan6 -r 6.0.0
#   mkdir -p /bigdata/stajichlab/shared/lib/interproscan6_data
#   mkdir -p /bigdata/stajichlab/shared/containers/interproscan6
#
# The first sample run will download InterPro data into --ips6_datadir.
# Subsequent samples reuse it.  Pin --ips6_interpro to a release number
# (e.g. 107.0) for reproducibility once the data is downloaded.

mkdir -p logs work.interproscan6

nextflow run pipeline/nextflow/interproscan6.nf \
    -profile slurm,interproscan6 \
    -w work.interproscan6 \
    --suppress suppress.txt \
    --ips6_interpro 108.0 \
    -resume \
    "$@"
