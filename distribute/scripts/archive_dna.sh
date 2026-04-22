#!/usr/bin/bash -l
#SBATCH -p hpcc_default -c 64 -N 1 -n 1 --out archive_dna.log

module load workspace/scratch
module load parallel
FOLDER=DNA
rsync -aL $FOLDER $SCRATCH/
pushd $SCRATCH/$FOLDER
parallel -j 16 pigz -p 4 {} ::: $(ls -U)
pushd $SCRATCH
tar cf /bigdata/stajichlab/shared/projects/1KFG/common_annotate/distribute_tar/$FOLDER.tar $FOLDER
