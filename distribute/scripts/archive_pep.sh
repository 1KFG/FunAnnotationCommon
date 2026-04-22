#!/usr/bin/bash -l
#SBATCH -p hpcc_default

FOLDER=pep
module load workspace/scratch
rsync -aL $FOLDER $SCRATCH/
pushd $SCRATCH
tar cf $FOLDER.tar $FOLDER
pigz $FOLDER.tar
rsync -a $FOLDER.tar.gz /bigdata/stajichlab/shared/projects/1KFG/common_annotate/distribute_tar
