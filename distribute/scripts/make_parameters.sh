#!/usr/bin/bash -l

mkdir -p parameters
pushd parameters

for a in $(ls -d ../../annotate/*); do 
	n=$(ls $a/predict_misc/ab_initio_parameters/*.genemark.mod)
	if [[ ! -z "$n" ]]; then 
		echo "($a)"
		target=$(basename $n .genemark.mod)
		if [ ! -d $target ]; then 
			ln -s $(dirname $n) $target
		fi 
	fi
done

popd
