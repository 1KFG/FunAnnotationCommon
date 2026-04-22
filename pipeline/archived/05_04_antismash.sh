#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 8 --mem 16G --out logs/antismash_4k.%a.log -J antismash -a 1-1889

module load antismash
which antismash
hostname
CPU=1
if [ ! -z $SLURM_CPUS_ON_NODE ]; then
  CPU=$SLURM_CPUS_ON_NODE
fi
TARGET=annotate
SAMPLES=samples.csv
N=${SLURM_ARRAY_TASK_ID}
if [ ! $N ]; then
  N=$1
  if [ ! $N ]; then
    echo "need to provide a number by --array or cmdline"
    exit
  fi
fi
N=$(expr $N + 4000)
MAX=`wc -l $SAMPLES | awk '{print $1}'`

if [ $N -gt $MAX ]; then
  echo "$N is too big, only $MAX lines in $SAMPLES"
  exit
fi

INPUTFOLDER=predict_results
tail -n +2 $SAMPLES | sed -n ${N}p | while IFS=, read ASMID SPECIESIN STRAIN BIOPROJECT NCBI_TAXONID BUSCO_LINEAGE PHYLUM SUBPHYLUM CLASS SUBCLASS ORDER FAMILY GENUS SPECIES LOCUSTAG
do
    OUT=$(echo -n $SPECIES | perl -p -e 'chomp; s/\s+/_/g')
    F=$(ls $TARGET/$OUT/$INPUTFOLDER/*.gbk | head -n 1)
    if [ -z $F ]; then
	    echo "No annotate dir for ${OUT} in $TARGET/$OUT"
	    exit
     fi
     echo "processing $OUT"
     JSON=($TARGET/$OUT/antismash_local/*.json)
     JSONGZ=($TARGET/$OUT/antismash_local/*.json.gz)
     if [[ -f ${JSON[0]} || -f ${JSONGZ[0]} ]]; then
	     echo "already produced ${JSON[0]} ${JSONGZ[0]}, skipping"
	     continue
     fi
     if [[ ! -d $TARGET/$OUT/antismash_local || ! -s $TARGET/$OUT/antismash_local/index.html ]]; then
	     rm -rf $TARGET/$OUT/antismash_local
	    time antismash --taxon fungi --output-dir $TARGET/$OUT/antismash_local \
		 --genefinding-tool none --fullhmmer --clusterhmmer --cb-general \
		 --pfam2go -c $CPU $F
     else
	     echo "antiSMASH already run for $name, skipping"
     fi
done
