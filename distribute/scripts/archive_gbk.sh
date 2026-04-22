#!/usr/bin/bash -l
#SBATCH -N 1 -c 48 --mem 24gb 

mkdir -p $SCRATCH/gbk
cd gbk
parallel -j 24 pigz -p 2 -c {} \> $SCRATCH/gbk/{}.gz ::: $(ls -U)
pushd $SCRATCH
tar cf /bigdata/stajichlab/shared/projects/1KFG/common_annotate/distribute_tar/gbk.tar gbk
