#!/usr/bin/bash -l
#SBATCH -p epyc --mem 500gb -c 16 -N 1 -n 1 -a 1 --out logs/nf_clean.%a.log

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
	-profile local,funannotate \
	-resume  \
	--only_clean \
	--samples ../samples_${N}_7500.csv \
	-w ../work.funannotate \
	--max_cpus $CPU \
	--suppress ../suppress.txt
