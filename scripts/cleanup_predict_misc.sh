#!/usr/bin/bash -l
#SBATCH -c 16 -N 1 -n 1 --mem 2gb
CPU=2
if [ $SLURM_CPUS_ON_NODE ]; then
  CPU=$SLURM_CPUS_ON_NODE
fi

deldirs() {
	in=$1
	taxaup=$(dirname $in)
	taxa=$(dirname $taxaup)
	echo $taxa
	rm -rf $taxa/predict_misc/busco
	rm -rf $taxa/predict_misc/EVM
	rm -rf $taxa/predict_misc/glimmerhmm
	rm -rf $taxa/predict_misc/busco_proteins
	rm -rf $taxa/predict_misc/busco
	rm -rf $taxa/predict_misc/proteins.combined.fa
	echo "done with $taxa."
}
module load parallel
export -f deldirs
ls annotate/*/predict_results/*.gbk > todel.txt

parallel -j $CPU deldirs ::: $(cat todel.txt)
