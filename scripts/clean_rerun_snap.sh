#!/usr/bin/ksh
SAMPLES=samples.csv
TARGET=annotate
IFS=,
TORUN=()
TORUN2=()
TORUN4=()
N=1

tail -n +2 $SAMPLES | while read ASMID SPECIES STRAIN BIOPROJECT NCBI_TAXONID BUSCO_LINEAGE PHYLUM SUBPHYLUM CLASS SUBCLASS ORDER FAMILY GENUS SPECIES LOCUSTAG
do
    SPECIEASNOSPACE=$(echo -n "$SPECIES" | perl -p -e 's/\s+/_/g')
    ADD=0
    if [[ -s $TARGET/$SPECIESNOSPACE/predict_misc/snap-predictions.gff3 ]]; then
        SNAPCOUNT=$(wc -l $TARGET/$SPECIESNOSPACE/predict_misc/snap-predictions.gff3 | awk '{print $1}')
        if [[ $SNAPCOUNT -lt 5 ]]; then  # we had successful snap predictions
            rm -rf $TARGET/$SPECIESNOSPACE/predict_misc/snap*
	    ADD=1
        fi
    elif [[ ! -f $TARGET/$SPECIESNOSPACE/predict_misc/snap-predictions.gff3 ]]; then
	ADD=1
    fi
    
    if [[ $ADD -eq 1 ]]; then
	if [[ $N -gt 4000 ]]; then
	    ADJ=$(expr $N - 4000)
	    TORUN4+=($ADJ)
	elif [[ $N -gt 2000 ]]; then
	    ADJ=$(expr $N - 2000)
	    TORUN2+=($ADJ)
	else
            TORUN+=($N)
	fi
    fi

    N=$(expr $N + 1)
done
RUNSET=$(echo "${TORUN[@]}" | perl -p -e 's/ /,/g')

echo "sbatch -a $RUNSET pipeline/01_predict.sh" 


RUNSET=$(echo "${TORUN2[@]}" | perl -p -e 's/ /,/g')

echo "sbatch -a $RUNSET pipeline/02_predict_2k.sh"

RUNSET=$(echo "${TORUN4[@]}" | perl -p -e 's/ /,/g')

echo "sbatch -a $RUNSET pipeline/03_predict_4k.sh" 
