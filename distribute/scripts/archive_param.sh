#!/usr/bin/bash -l
#SBATCH -N 1 -n 32 --mem 24gb  --out archive_params.%A.log

mkdir -p $SCRATCH/gene_prediction_params
pushd parameters
parallel -j 32 tar cfhz $SCRATCH/gene_prediction_params/{}.tar.gz {} ::: $(ls -U)
pushd $SCRATCH
tar cf /bigdata/stajichlab/shared/projects/1KFG/common_annotate/distribute_tar/gene_prediction_params.tar gene_prediction_params
popd
popd
