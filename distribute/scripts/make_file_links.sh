#!/usr/bin/bash -l

mkdir -p cds  DNA  gbk  gff pep

pushd cds
ln -s ../../annotate/*/predict_results/*.cds-transcripts.fa .
popd

pushd pep
ln -s ../../annotate/*/predict_results/*.proteins.fa .
popd

pushd gbk
ln -s ../../annotate/*/predict_results/*.gbk .
popd

pushd DNA
ln -s ../../annotate/*/predict_results/*.scaffolds.fa .
popd

pushd gff
ln -s ../../annotate/*/predict_results/*.gff3 .
popd

