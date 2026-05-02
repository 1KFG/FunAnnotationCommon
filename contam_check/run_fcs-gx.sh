#!/usr/bin/bash -l
#SBATCH --mem 512gb -N 1 -n 1 -c 64 --out contam_purge.log

module load AAFTF
hostname
rsync -a --progress /srv/projects/db/ncbi-fcs/0.5.4/gxdb /dev/shm/
#Fungi
TAXID=4751
INDIR=../distribute/DNA
OUTDIR=contam_reports
mkdir -p $OUTDIR
run_fcs() {
    local f="$1"
    local base="${f%.*}"
    local INPUT="$INDIR/$f"
    local REPORT="$OUTDIR/$base.report"
    if [ ! -f "$REPORT" ] || [ "$INPUT" -nt "$REPORT" ]; then
        AAFTF fcs_gx_purge --db /dev/shm/gxdb/all -i "$INPUT" --cpus 8 -o "$OUTDIR/$base.purge.fasta" -t "${TAXID}" -w "$REPORT"
    fi
}
export -f run_fcs
export INDIR OUTDIR TAXID

parallel -j 8 run_fcs ::: $(ls -U $INDIR)

rm -rf /dev/shm/gxdb
