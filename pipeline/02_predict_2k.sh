#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 16 --mem 48gb --out logs/funannotate_predict_2k.%a.log -a 1-2000 --time 24:00:00

module load funannotate
SAMPLES=samples.csv
SOURCE=/bigdata/stajichlab/shared/projects/1KFG/2021/NCBI_fungi/source/NCBI_ASM
TARGET=annotate
SEQCENTER=NCBI
export AUGUSTUS_CONFIG_PATH=$(realpath lib/augustus/3.5/config)
mkdir -p $TARGET
CPU=2
if [ ! -z $SLURM_CPUS_ON_NODE ]; then
    CPU=$SLURM_CPUS_ON_NODE
fi
N=$SLURM_ARRAY_TASK_ID
if [ -z $N ]; then
    N=$1
    if [ -z $N ]; then
        echo "need to provide a number by --array or cmdline"
        exit
    fi
fi
N=$(expr $N + 2000)
export FUNANNOTATE_DB=/bigdata/stajichlab/shared/lib/funannotate_db
IFS=,
tail -n +2 $SAMPLES | sed -n ${N}p | while read ASMID SPECIES STRAIN BIOPROJECT NCBI_TAXONID BUSCO_LINEAGE PHYLUM SUBPHYLUM CLASS SUBCLASS ORDER FAMILY GENUS SPECIES LOCUSTAG
do
    LOCUSTAG=$(echo -n "$LOCUSTAG" | perl -p -e 's/[\r\n]//g')
    echo "Running $ASM for $SPECIES $STRAIN ( $BUSCO_LINEAGE, $LOCUSTAG )"
    GENOMEGZ=$SOURCE/$ASMID/$ASMID\_genomic.fna.gz
    GENOME=$SCRATCH/$ASMID.fa
    #pigz -dc $GENOMEGZ | perl -p -e 's/>(\S+)\s+.+/>$1/' > $GENOME
    pigz -dc $GENOMEGZ | ./scripts/clean_genome_fa.py --len 2000 > $GENOME
    OUT=$(echo -n $SPECIES | perl -p -e 'chomp; s/\s+/_/g')
    time funannotate predict --name $LOCUSTAG -i $GENOME --strain "$STRAIN" -o $TARGET/$OUT -s "$SPECIES" --cpu $CPU --busco_db $BUSCO_LINEAGE \
        --AUGUSTUS_CONFIG_PATH $AUGUSTUS_CONFIG_PATH -w codingquarry:0 --min_training_models 50 --tmpdir $SCRATCH --SeqCenter $SEQCENTER --keep_no_stops --header_length 24

    F=$(ls $TARGET/$OUT/predict_results/*.gbk | head -n 1)
    if [ ! -z $F ]; then
        rm -rf $TARGET/$OUT/predict_misc/EVM $TARGET/$OUT/predict_misc/proteins.combined.fa
        rm -rf $TARGET/$OUT/predict_misc/glimmerhmm
        rm -rf $TARGET/$OUT/predict_misc/busco
        rm -rf $TARGET/$OUT/predict_misc/busco_proteins
    fi
done
