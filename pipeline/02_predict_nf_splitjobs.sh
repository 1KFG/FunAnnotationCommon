#!/usr/bin/bash -l
#SBATCH --mem 16gb -c 1 -N 1 -n 1 -a 1 --out logs/nf_predict_split.%a.log

hostname
CPU=$SLURM_CPUS_ON_NODE
if [ -z "$CPU" ]; then
	CPU=24
fi
N=${SLURM_ARRAY_TASK_ID}
if [ -z "$N" ]; then
	N=$1
	if [ -z "${N}" ]; then
		echo "need to provide and array-id or cmdline arg for number to run (1-3)."
		N=1
		echo "defaulting to 1"
	fi
fi

cd run${N}
mkdir -p logs

nextflow run pipeline/nextflow/funannotate.nf \
    -profile slurm,funannotate \
    --suppress ../suppress.txt \
    --samples ../samples_${N}_7500.csv  \
    -resume \
    "$@"


