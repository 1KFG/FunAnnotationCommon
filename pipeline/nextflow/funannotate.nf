#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.samples         = "${launchDir}/samples.csv"
params.target          = "${launchDir}/annotate"
params.taxondb	       = "/bigdata/stajichlab/shared/projects/1KFG/2026/NCBI_fungi/tmp/taxa/"
params.proteins	       = "${launchDir}/lib/swissprot_fungi.faa"
params.source          = "/bigdata/stajichlab/shared/projects/1KFG/2026/NCBI_fungi/source/NCBI_ASM"
params.seqcenter       = "NCBI"
params.augustus_config = "${launchDir}/lib/augustus/3.5/config"
params.funannotate_db  = "/bigdata/stajichlab/shared/lib/funannotate_db"
params.min_contig_len  = 2000
params.clean_script    = "${launchDir}/scripts/clean_genome_fa.py"
params.sbt_template    = "${launchDir}/lib/template.sbt"
params.debug           = false   // --debug: verbose logging in script + channel views
params.n_test          = 0       // --n_test N: limit to first N samples (0 = all)
params.max_cpus        = 64 // --max_cpus N: total CPUs for local executor
params.suppress        = ""      // --suppress path/to/file: ASMID per line, optional comma+comment; # = comment line
params.max_rnaseq_runs = 4       // --max_rnaseq_runs N: max paired-end SRA sets to download per species
params.run_antismash   = false   // --run_antismash: run antiSMASH BGC detection after predict
params.run_interpro    = false   // --run_interpro: run InterProScan before funannotate annotate
params.run_signalp     = false   // --run_signalp: run SignalP 6 (requires GPU node) before funannotate annotate
params.antismash_taxon = "fungi" // --antismash_taxon: antiSMASH --taxon value
params.only_clean      = false   // --only_clean: stop after GENOME_CLEAN (skip prediction and all post-predict steps)
params.run_sra_fetch   = true    // --run_sra_fetch false: skip SRA download + funannotate train
params.skip_repeatmasker = false // --skip_repeatmasker: skip RepeatModeler+RepeatMasker steps
params.pasa_mysql    = false   // --pasa_mysql: start a per-task MariaDB instance for PASA
params.mariadb_sif   = "/bigdata/stajichlab/shared/lib/mariadb/mariadb.sif"
params.mysql_datadir = ""      // path to template MySQL data dir (required with --pasa_mysql)
params.pasa_conf_dir = ""      // path to dir with my.cnf + conf.txt (required with --pasa_mysql)


// Metadata tuple order used throughout:
//   val(out), val(asmid), val(species), val(strain), val(locustag),
//   val(busco_lineage), val(header_length), val(transl_table)
// GENOME_CLEAN receives: ..., path(genome_gz), val(taxonid), val(taxondb)
//   → emits: ..., path(genome_fa), val(taxonid)   [storeDir moves .fa; workflow maps to abs string]
//   → writes <asmid>.fa to input_clean_genomes/ (storeDir; skip check targets this file)
//   → purge/FCS intermediates written as side effects to input_clean_genomes/clean/
// REPEATMODELER_RUN receives: val(species_tag), val(asmid), val(genome_fa)  [one per species]
//   → emits: val(species_tag), path(rmlib)   [storeDir caches repeat_library/<species_tag>.RMlib.fasta]
// REPEATMASKER_RUN receives: val(species_tag), ..., val(genome_fa), val(taxonid), path(rmlib)
//   → emits: ..., path(masked_fa), val(taxonid)   [storeDir caches input_clean_genomes/<asmid>.masked.fasta]
//   [skipped when --skip_repeatmasker; masked_fa falls back to unmasked .fa if .masked.fasta absent]
// SRA_FETCH receives: ..., val(genome_fa), val(taxonid)   [only when --run_sra_fetch]
//   → emits: ..., val(genome_fa), path(reads_dir)   [reads_dir may be empty]
// FUNANNOTATE_TRAIN receives: ..., val(genome_fa), path(reads_dir)   [only when --run_sra_fetch]
//   → emits: ..., val(genome_fa)   [reads deleted after training]
// FUNANNOTATE_PREDICT receives: ..., val(genome_fa)   [from TRAIN or directly after masking/clean]

// Download and extract NCBI taxdump once; storeDir caches it at params.taxondb so
// subsequent runs skip this entirely.
process SETUP_TAXONDB {
    storeDir params.taxondb

    cpus   1
    memory '4 GB'
    time   '1h'

    output:
    path "names.dmp",    emit: ready
    path "nodes.dmp"
    path "merged.dmp"
    path "delnodes.dmp"
    path "division.dmp"
    path "gencode.dmp"
    path "citations.dmp"

    script:
    """
    set -euo pipefail
    wget --no-verbose https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
    tar zxf taxdump.tar.gz
    rm taxdump.tar.gz
    """

    stub:
    """
    for f in names.dmp nodes.dmp merged.dmp delnodes.dmp division.dmp gencode.dmp citations.dmp; do
        touch \$f
    done
    """
}

process GENOME_CLEAN {
    tag "$asmid"

    container '/rhome/jstajich/projects/AAFTF/AAFTF_v0.6.1-signed.sif'

    // Nextflow skips this task when input_clean_genomes/<asmid>.fa already exists.
    storeDir "${launchDir}/input_clean_genomes"

    cpus   16
    memory '450 GB'
    time   '6h'

    input:
    tuple val(out), val(asmid), val(species), val(strain), val(locustag),
          val(busco_lineage), val(header_length), val(transl_table),
          path(genome_gz), val(taxonid), val(taxondb)

    output:
    tuple val(out), val(asmid), val(species), val(strain), val(locustag),
          val(busco_lineage), val(header_length), val(transl_table),
          path("${asmid}.fa"), val(taxonid), emit: genome

    script:
    """
    if [ ! -f "${genome_gz}" ]; then
        echo "ERROR: genome_gz not found at path: ${genome_gz}" >&2
        exit 1
    fi
    module load AAFTF

    # Ensure /dev/shm/gxdb is present on this node; register for cleanup when done.
    source ${launchDir}/scripts/setup_fcs_shm.sh
    SCRATCH=\$(printf '%s' "\${SCRATCH}" | tr -d '\\n\\r')
    TAXONKIT_DB=${taxondb}
    module load taxonkit
    phylum=\$(echo ${taxonid} | taxonkit --data-dir \$TAXONKIT_DB lineage | taxonkit --data-dir \$TAXONKIT_DB reformat -f "{p}" | cut -f3 | taxonkit --data-dir \$TAXONKIT_DB name2taxid | cut -f2 | uniq | head -n 1)
    if [ -z "\$phylum" ]; then
    	phylum=\$(echo ${taxonid} | taxonkit --data-dir \$TAXONKIT_DB lineage | taxonkit --data-dir \$TAXONKIT_DB reformat -f "{K}" | cut -f3 | taxonkit --data-dir \$TAXONKIT_DB name2taxid | uniq | cut -f2 | head -n 1)
	# weird we are getting 2 lines from name2taxid when input is Fungi add the uniq/head -n 1 to ensure only one line
    fi
    module unload taxonkit
    echo "[INFO] Phylum for ${asmid} (taxonid=${taxonid}): \$phylum"
    echo "[INFO] Decompressing and cleaning genome for ${asmid}..."
    pigz -dc ${genome_gz} > \$SCRATCH/${asmid}.raw.fa
    AAFTF fcs_gx_purge --db /dev/shm/gxdb/all \
        -i \$SCRATCH/${asmid}.raw.fa --cpus ${task.cpus} \
        -o \$SCRATCH/${asmid}.purge.fasta \
        -t "\$phylum" -w \$SCRATCH/${asmid}.fcs_report
    mkdir -p ${launchDir}/input_clean_genomes/clean
    cat \$SCRATCH/${asmid}.purge.fasta | \
        ${params.clean_script} --len ${params.min_contig_len} > ${asmid}.fa
    echo "[INFO] Clean genome written: ${asmid}.fa (\$(du -sh ${asmid}.fa | cut -f1))"
    pigz \$SCRATCH/${asmid}.purge.fasta 
    pigz \$SCRATCH/${asmid}.purge.fcs_gx-taxonomy.tsv 
    mv \$SCRATCH/${asmid}.purge.fasta.gz \$SCRATCH/${asmid}.purge.fcs_gx-taxonomy.tsv.gz ${launchDir}/input_clean_genomes/clean/
    """

    stub:
    """
    echo ">stub_${asmid}" > ${asmid}.fa
    mkdir -p ${launchDir}/input_clean_genomes/clean
    touch ${launchDir}/input_clean_genomes/clean/${asmid}.purge.fasta
    touch ${launchDir}/input_clean_genomes/clean/${asmid}.purge.fcs_gx-taxonomy.tsv
    """
}

// Run RepeatModeler on one representative assembly per species to build a de-novo
// repeat library.  storeDir caches the library so the process is skipped on re-runs
// or when another assembly of the same species already produced it.
//
// Output naming convention:
//   repeat_library/{species_tag}.{asmid}.RM_lib.fasta  — the canonical file (storeDir output)
//   repeat_library/{species_tag}.RMlib.fasta           — symlink to the above (for downstream use)
//   repeat_library/library_manifest.tsv                — records which asmid built each library
process REPEATMODELER_RUN {
    tag "$species_tag"

    storeDir "${launchDir}/repeat_library"

    cpus   16
    memory '32 GB'
    time   '72h'

    input:
    tuple val(species_tag), val(asmid), val(genome_fa)

    output:
    tuple val(species_tag), path("${species_tag}.${asmid}.RM_lib.fasta"), emit: rmlib

    script:
    """
    module load RepeatModeler
    DBNAME=${asmid}_rmdb
    BuildDatabase -name \$DBNAME -engine ncbi ${genome_fa}
    RepeatModeler -database \$DBNAME -threads ${task.cpus} -LTRStruct
    if [ -f "\${DBNAME}-families.fa" ]; then
        cp \${DBNAME}-families.fa ${species_tag}.${asmid}.RM_lib.fasta
    else
        echo "[WARN] RepeatModeler found no families for ${asmid}; creating empty library"
        touch ${species_tag}.${asmid}.RM_lib.fasta
    fi
    pigz -dc ${launchDir}/lib/fungi_repeat.20170127.lib.gz >> ${species_tag}.${asmid}.RM_lib.fasta

    # Copy to storeDir so the symlink target exists before Nextflow moves the declared output
    mkdir -p ${launchDir}/repeat_library
    cp ${species_tag}.${asmid}.RM_lib.fasta \
        ${launchDir}/repeat_library/${species_tag}.${asmid}.RM_lib.fasta
    # (Re-)create the convenience symlink; relative so it resolves within repeat_library/
    ln -sf ${species_tag}.${asmid}.RM_lib.fasta \
        ${launchDir}/repeat_library/${species_tag}.RMlib.fasta

    # Append provenance record; initialise header on first entry
    MANIFEST="${launchDir}/repeat_library/library_manifest.tsv"
    if [ ! -f "\$MANIFEST" ]; then
        printf "species_tag\tasmid\trmlib_file\ttimestamp\n" > "\$MANIFEST"
    fi
    printf "%s\t%s\t%s\t%s\n" \
        "${species_tag}" "${asmid}" \
        "${species_tag}.${asmid}.RM_lib.fasta" \
        "\$(date -Iseconds)" >> "\$MANIFEST"
    """

    stub:
    """
    echo ">stub_repeat_${species_tag}" > ${species_tag}.${asmid}.RM_lib.fasta
    pigz -dc ${launchDir}/lib/fungi_repeat.20170127.lib.gz >> ${species_tag}.${asmid}.RM_lib.fasta
    mkdir -p ${launchDir}/repeat_library
    cp ${species_tag}.${asmid}.RM_lib.fasta \
        ${launchDir}/repeat_library/${species_tag}.${asmid}.RM_lib.fasta
    ln -sf ${species_tag}.${asmid}.RM_lib.fasta \
        ${launchDir}/repeat_library/${species_tag}.RMlib.fasta
    MANIFEST="${launchDir}/repeat_library/library_manifest.tsv"
    if [ ! -f "\$MANIFEST" ]; then
        printf "species_tag\tasmid\trmlib_file\ttimestamp\n" > "\$MANIFEST"
    fi
    printf "%s\t%s\t%s\t%s\n" \
        "${species_tag}" "${asmid}" \
        "${species_tag}.${asmid}.RM_lib.fasta" \
        "\$(date -Iseconds)" >> "\$MANIFEST"
    """
}

// Soft-mask each assembly genome using the per-species repeat library produced by
// REPEATMODELER_RUN.  storeDir caches the masked FASTA alongside the clean genome.
// The full RepeatMasker output folder is written as a side effect to repeat_masker/<asmid>/.
process REPEATMASKER_RUN {
    tag "$asmid"

    storeDir "${launchDir}/input_clean_genomes"

    cpus   32
    memory '24 GB'
    time   '24h'

    input:
    tuple val(species_tag), val(out), val(asmid), val(species), val(strain), val(locustag),
          val(busco_lineage), val(header_length), val(transl_table),
          val(genome_fa), val(taxonid), path(rmlib)

    output:
    tuple val(out), val(asmid), val(species), val(strain), val(locustag),
          val(busco_lineage), val(header_length), val(transl_table),
          path("${asmid}.masked.fasta"), val(taxonid), emit: masked

    script:
    """
    module load RepeatMasker
    mkdir -p rm_out
    RepeatMasker -lib ${rmlib} -pa ${task.cpus} -dir rm_out -noisy -xsmall ${genome_fa}
    MASKED=rm_out/\$(basename ${genome_fa}).masked
    if [ ! -f "\$MASKED" ]; then
        echo "[WARN] RepeatMasker produced no .masked file; copying unmasked genome"
        cp ${genome_fa} ${asmid}.masked.fasta
    else
        cp \$MASKED ${asmid}.masked.fasta
    fi
    mkdir -p ${launchDir}/repeat_masker/${asmid}
    mv rm_out/* ${launchDir}/repeat_masker/${asmid}/ 2>/dev/null || true
    """

    stub:
    """
    echo ">stub_${asmid}_masked" > ${asmid}.masked.fasta
    mkdir -p ${launchDir}/repeat_masker/${asmid}
    """
}

// Search NCBI SRA for paired-end RNA-seq runs for this taxon and download up to
// params.max_rnaseq_runs sets.  Emits a reads/ directory (possibly empty) so
// FUNANNOTATE_TRAIN always has a consistent input regardless of SRA availability.
process SRA_FETCH {
    tag "$asmid"

    cpus   8
    memory '16 GB'
    time   '6h'

    input:
    tuple val(out), val(asmid), val(species), val(strain), val(locustag),
          val(busco_lineage), val(header_length), val(transl_table),
          val(genome_fa), val(taxonid)

    output:
    tuple val(out), val(asmid), val(species), val(strain), val(locustag),
          val(busco_lineage), val(header_length), val(transl_table),
          val(genome_fa), path("reads"), emit: reads

    script:
    """
    mkdir -p reads
    module load ncbi_edirect
    module load sratoolkit
    module load parallel-fastq-dump

    ACCESSIONS=\$(esearch -db sra \
        -query "txid${taxonid}[Organism:noexp] AND RNA-Seq[Strategy] AND PAIRED[Layout] AND Illumina[Platform]" | \
        efetch -format runinfo | \
        awk -F',' 'NR>1 && \$13=="RNA-Seq" && \$16=="PAIRED" && \$1~/^[SDE]RR/ {print \$1}' | \
        sort -u | head -n ${params.max_rnaseq_runs} || true)

    if [ -z "\$ACCESSIONS" ]; then
        echo "[INFO] No paired-end RNA-seq runs found for ${species} (taxonid=${taxonid})"
        exit 0
    fi

    echo "[INFO] SRA accessions for ${species}: \$ACCESSIONS"
    TMPDIR=\${SCRATCH:-/tmp}

    for ACC in \$ACCESSIONS; do
        echo "[INFO] Downloading \$ACC ..."
        parallel-fastq-dump --sra-id \$ACC --threads ${task.cpus} -X 500 \\
            --outdir reads/ --split-files --gzip --tmpdir \$TMPDIR || {
            echo "[WARN] Download failed for \$ACC, skipping"
            #rm -f reads/\${ACC}*.fastq.gz
        }
    done

    NPAIRS=\$(ls reads/*_1.fastq.gz 2>/dev/null | wc -l)
    echo "[INFO] Downloaded \$NPAIRS paired read set(s) for ${species}"
    """

    stub:
    """
    mkdir -p reads
    echo "[STUB] SRA_FETCH noop for ${out} (taxonid=${taxonid})"
    """
}

process FUNANNOTATE_TRAIN {
    tag "$out"

    cpus   16
    memory '96 GB'
    time   '120h'

    input:
    tuple val(out), val(asmid), val(species), val(strain), val(locustag),
          val(busco_lineage), val(header_length), val(transl_table),
          val(genome_fa), path(reads_dir)

    output:
    tuple val(out), val(asmid), val(species), val(strain), val(locustag),
          val(busco_lineage), val(header_length), val(transl_table),
          val(genome_fa)

    script:
    def pasa_db_arg = params.pasa_mysql ? "--pasa_db mysql" : ""
    """
    # ── Skip if training output already present ───────────────────────────────
    TRAIN_GFF3="${params.target}/${out}/training/funannotate_train.pasa.gff3"
    if [ -f "\$TRAIN_GFF3" ]; then
        echo "[INFO] Training already complete for ${out}; skipping"
        exit 0
    fi

    # ── Detect reads ──────────────────────────────────────────────────────────
    R1=(\$(ls ${reads_dir}/*_1.fastq.gz 2>/dev/null || true))
    R2=(\$(ls ${reads_dir}/*_2.fastq.gz 2>/dev/null || true))
    SE=(\$(ls ${reads_dir}/*.fastq.gz 2>/dev/null | grep -v '_[12]\\.fastq\\.gz\$' || true))

    if [ \${#R1[@]} -eq 0 ] && [ \${#SE[@]} -eq 0 ]; then
        echo "[INFO] No RNAseq reads for ${out}, skipping funannotate train"
        exit 0
    fi
    if [ \${#R1[@]} -gt 0 ] && [ \${#R1[@]} -ne \${#R2[@]} ]; then
        echo "[WARN] Unequal R1/R2 counts for ${out} (R1=\${#R1[@]}, R2=\${#R2[@]}), skipping train"
        exit 0
    fi

    source /etc/profile.d/modules.sh 2>/dev/null || true
    module load miniconda3
    eval "\$(conda shell.bash hook)"
    module load funannotate

    export AUGUSTUS_CONFIG_PATH=${params.augustus_config}
    export FUNANNOTATE_DB=${params.funannotate_db}
    TMPDIR=\${SCRATCH:-/tmp}

    # ── Optional per-task MariaDB for PASA ────────────────────────────────────
    if [ "${params.pasa_mysql}" = "true" ]; then
        RUNID=\$\$
        MYSQL_SCRATCH=\$TMPDIR/mysql_\${RUNID}
        mkdir -p \$MYSQL_SCRATCH/db \$MYSQL_SCRATCH/conf
        rsync -a ${params.mysql_datadir}/mysql \$MYSQL_SCRATCH/db/ || \
            { echo "ERROR: Failed to copy mysql data from ${params.mysql_datadir}" >&2; exit 1; }
        cp ${params.pasa_conf_dir}/my.cnf \$MYSQL_SCRATCH/conf/my.cnf || \
            { echo "ERROR: Failed to copy my.cnf" >&2; exit 1; }
        MYHOSTNAME=\$(hostname -s)
        PORT=\$(shuf -i3000-4999 -n1)
        PASACONF=\$MYSQL_SCRATCH/conf/pasa-local-\${MYHOSTNAME}.config.txt
        cp ${params.pasa_conf_dir}/conf.txt \$PASACONF
        sed -i "s/^MYSQLSERVER.*\$/MYSQLSERVER=\${MYHOSTNAME}:\${PORT}/" \$PASACONF
        perl -i -p -e "s/port = \\d+/port = \${PORT}/" \$MYSQL_SCRATCH/conf/my.cnf
        export SINGULARITY_BINDPATH=\$TMPDIR
        export PASACONF
        stop_mysqldb() { singularity instance stop mysqldb\${RUNID} 2>/dev/null || true; }
        trap "stop_mysqldb; exit 130" SIGHUP SIGINT SIGTERM
        trap "stop_mysqldb" EXIT
        module load singularity
        singularity instance start --writable-tmpfs \\
            -B \$MYSQL_SCRATCH/conf/my.cnf:/etc/mysql/my.cnf,\$MYSQL_SCRATCH/db/:/var/lib/mysql,\$MYSQL_SCRATCH/conf:/usr/conf \\
            ${params.mariadb_sif} mysqldb\${RUNID} /usr/bin/mysqld_safe
        sleep 5
    fi

    # ── Run funannotate train ─────────────────────────────────────────────────
    if [ \${#R1[@]} -gt 0 ]; then
        echo "[INFO] Running funannotate train (paired-end) for ${out} with \${#R1[@]} read pair(s)"
        funannotate train -i ${genome_fa} -o ${params.target}/${out} \\
            --left \${R1[@]} --right \${R2[@]} \\
            --species "${species}" --strain "${strain}" \\
            --cpus ${task.cpus} --memory "${task.memory}" \\
            --jaccard_clip --no-progress --min_coverage 4 \\
            --header_length ${header_length} \\
            --tmpdir \$TMPDIR ${pasa_db_arg}
    else
        echo "[INFO] Running funannotate train (single-end) for ${out} with \${#SE[@]} read file(s)"
        funannotate train -i ${genome_fa} -o ${params.target}/${out} \\
            --single \${SE[@]} \\
            --species "${species}" --strain "${strain}" \\
            --cpus ${task.cpus} --memory "${task.memory}" \\
            --jaccard_clip --no-progress --min_coverage 4 \\
            --header_length ${header_length} \\
            --tmpdir \$TMPDIR ${pasa_db_arg}
    fi

    # ── Reclaim SRA reads immediately ─────────────────────────────────────────
    REAL_READS=\$(readlink -f ${reads_dir})
    if [ "\$REAL_READS" != "${reads_dir}" ]; then
        rm -rf "\$REAL_READS"
        echo "[INFO] Removed SRA reads at \$REAL_READS"
    fi
    """

    stub:
    """
    echo "[STUB] FUNANNOTATE_TRAIN stub for ${out} (reads_dir=${reads_dir})"
    mkdir -p ${params.target}/${out}/training
    touch ${params.target}/${out}/training/funannotate_train.pasa.gff3
    """
}

process FUNANNOTATE_PREDICT {
    tag "$out"

    cpus   16
    memory '32 GB'
    time   '32h'

    publishDir "${params.target}", mode: 'copy', overwrite: true

    input:
    tuple val(out), val(asmid), val(species), val(strain), val(locustag),
          val(busco_lineage), val(header_length), val(transl_table),
          val(genome_fa)

    output:
    tuple val(out), path("${out}/**"), emit: files
    tuple val(out), val(asmid), val(species), val(strain), val(locustag),
          val(busco_lineage), val(header_length), val(transl_table), emit: metadata

    script:
    """
    source /etc/profile.d/modules.sh 2>/dev/null || true
    module load miniconda3
    eval "\$(conda shell.bash hook)"
    module load funannotate

    export AUGUSTUS_CONFIG_PATH=${params.augustus_config}
    export FUNANNOTATE_DB=${params.funannotate_db}
    TMPDIR=\${SCRATCH:-/tmp}

    if [ "${params.debug}" = "true" ]; then
        echo "[DEBUG] out          = ${out}"
        echo "[DEBUG] asmid        = ${asmid}"
        echo "[DEBUG] species      = ${species}"
        echo "[DEBUG] strain       = ${strain}"
        echo "[DEBUG] locustag     = ${locustag}"
        echo "[DEBUG] busco        = ${busco_lineage}"
        echo "[DEBUG] transl_table = ${transl_table}"
        echo "[DEBUG] proteins     = ${params.proteins}"
        echo "[DEBUG] genome_fa    = ${genome_fa}"
        echo "[DEBUG] TMPDIR       = \$TMPDIR"
        echo "[DEBUG] pwd          = \$(pwd)"
    fi

    TBL2ASN_PARAMS="-l paired-ends"

    funannotate predict --name ${locustag} -i ${genome_fa} --strain "${strain}" \\
        -o ${out} -s "${species}" --cpu ${task.cpus} --busco_db ${busco_lineage} \\
        --AUGUSTUS_CONFIG_PATH \$AUGUSTUS_CONFIG_PATH -w codingquarry:0 \\
        --min_training_models 30 --tmpdir \$TMPDIR --SeqCenter ${params.seqcenter} \\
        --keep_no_stops --header_length ${header_length} --protein_evidence ${params.proteins} \\
        --tbl2asn "\$TBL2ASN_PARAMS" --table ${transl_table}

    EXPECTED_GBK="${out}/predict_results/${out}.gbk"
    if [ ! -f "\$EXPECTED_GBK" ]; then
        echo "ERROR: funannotate predict did not produce expected GBK: \$EXPECTED_GBK" >&2
        exit 1
    fi
    mv ${out}/predict_misc/ab_initio_parameters ${out}
    rm -rf ${out}/predict_misc
    mkdir -p ${out}/predict_misc
    mv ${out}/ab_initio_parameters ${out}/predict_misc
    pigz ${out}/predict_results/*.txt ${out}/predict_results/*.mrna-transcripts.fa
    """

    stub:
    """
    echo "[STUB] Would run funannotate predict for ${out} using ${genome_fa}"
    [ -f "${genome_fa}" ] || { echo "ERROR: genome not found at ${genome_fa}" >&2; exit 1; }
    mkdir -p ${out}/predict_results ${out}/predict_misc
    touch ${out}/predict_results/${out}.gbk ${out}/predict_results/${out}.proteins.fa
    """
}

process ANTISMASH_RUN {
    tag "$out"

    cpus   8
    memory '16 GB'
    time   '60h'

    publishDir "${params.target}", mode: 'copy', overwrite: true

    input:
    tuple val(out), val(asmid), val(species), val(strain), val(locustag),
          val(busco_lineage), val(header_length), val(transl_table)

    output:
    tuple val(out), path("${out}/antismash_local/**")

    script:
    def gbk = "${params.target}/${out}/predict_results/${out}.gbk"
    """
    if [ ! -f "${gbk}" ]; then
        echo "ERROR: predict GBK not found: ${gbk}" >&2
        exit 1
    fi
    source /etc/profile.d/modules.sh 2>/dev/null || true
    module load miniconda3
    eval "\$(conda shell.bash hook)"
    module load antismash
    antismash --taxon ${params.antismash_taxon} \\
        --output-dir ${out}/antismash_local \\
        --genefinding-tool none \\
        --fullhmmer --clusterhmmer --cb-general --pfam2go \\
        -c ${task.cpus} \\
        ${gbk}
    pigz ${out}/antismash_local/*.json
    """

    stub:
    """
    mkdir -p ${out}/antismash_local
    touch ${out}/antismash_local/${out}.json.gz
    touch ${out}/antismash_local/index.html
    """
}

// IPRSCAN5
process INTERPROSCAN_RUN {
    tag "$out"

    cpus   8
    memory '32 GB'
    time   '60h'

    publishDir "${params.target}", mode: 'copy', overwrite: true

    input:
    tuple val(out), val(asmid), val(species), val(strain), val(locustag),
          val(busco_lineage), val(header_length), val(transl_table)

    output:
    tuple val(out), path("${out}/annotate_misc/iprscan.xml")

    script:
    def proteins = "${params.target}/${out}/predict_results/${out}.proteins.fa"
    """
    if [ ! -f "${proteins}" ]; then
        echo "ERROR: protein FASTA not found: ${proteins}" >&2
        exit 1
    fi
    mkdir -p ${out}/annotate_misc
    module load interproscan
    interproscan.sh -i ${proteins} -f XML -o ${out}/annotate_misc/iprscan.xml \\
        -dp -goterms -pa -t p -cpu ${task.cpus}
    """

    stub:
    """
    mkdir -p ${out}/annotate_misc
    touch ${out}/annotate_misc/iprscan.xml
    """
}

process SIGNALP_RUN {
    tag "$out"

    publishDir "${params.target}", mode: 'copy', overwrite: true

    input:
    tuple val(out), val(asmid), val(species), val(strain), val(locustag),
          val(busco_lineage), val(header_length), val(transl_table)

    output:
    tuple val(out), path("${out}/annotate_misc/signalp.results.txt")

    script:
    def proteins = "${params.target}/${out}/predict_results/${out}.proteins.fa"
    """
    if [ ! -f "${proteins}" ]; then
        echo "ERROR: protein FASTA not found: ${proteins}" >&2
        exit 1
    fi
    module load signalp/6-gpu
    TMPDIR=\${SCRATCH:-/tmp}
    signalp6 -od \$TMPDIR/${out}_signalp \\
        -org euk --mode fast -format txt \\
        -fasta ${proteins} \\
        --write_procs ${task.cpus} -bs 16
    mkdir -p ${out}/annotate_misc
    cp \$TMPDIR/${out}_signalp/prediction_results.txt ${out}/annotate_misc/signalp.results.txt
    rm -rf \$TMPDIR/${out}_signalp
    """

    stub:
    """
    mkdir -p ${out}/annotate_misc
    touch ${out}/annotate_misc/signalp.results.txt
    """
}

process FUNANNOTATE_ANNOTATE {
    tag "$out"

    cpus   16
    memory '32 GB'
    time   '48h'

    publishDir "${params.target}", mode: 'copy', overwrite: true

    input:
    tuple val(out), val(asmid), val(species), val(strain), val(locustag),
          val(busco_lineage), val(header_length), val(transl_table)

    output:
    tuple val(out), path("${out}/**")

    script:
    def ipr     = file("${params.target}/${out}/annotate_misc/iprscan.xml")
    def iprArg  = ipr.exists()     ? "--iprscan ${ipr}"    : ""
    def sp      = file("${params.target}/${out}/annotate_misc/signalp.results.txt")
    def spArg   = sp.exists()      ? "--signalp ${sp}"     : ""
    """
    source /etc/profile.d/modules.sh 2>/dev/null || true
    module load miniconda3
    eval "\$(conda shell.bash hook)"
    module load funannotate

    export AUGUSTUS_CONFIG_PATH=${params.augustus_config}
    export FUNANNOTATE_DB=${params.funannotate_db}
    TMPDIR=\${SCRATCH:-/tmp}

    funannotate annotate -i ${params.target}/${out} \\
        --species "${species}" --strain "${strain}" \\
        --busco_db ${busco_lineage} --rename ${locustag} \\
        --sbt ${params.sbt_template} \\
        ${iprArg} ${spArg} \\
        --cpu ${task.cpus} --tmpdir \$TMPDIR

    EXPECTED_GBK="${out}/annotate_results/${out}.gbk"
    if [ ! -f "\$EXPECTED_GBK" ]; then
        echo "ERROR: funannotate annotate did not produce expected GBK: \$EXPECTED_GBK" >&2
        exit 1
    fi
    """

    stub:
    """
    echo "[STUB] Would run funannotate annotate for ${out}"
    mkdir -p ${out}/annotate_results
    touch ${out}/annotate_results/${out}.gbk
    """
}

workflow {
    def suppressSet = (params.suppress && file(params.suppress).exists())
        ? file(params.suppress).readLines()
              .collect { it.trim().split(',')[0].trim() }
              .findAll { it && !it.startsWith('#') }
              .toSet()
        : ([] as Set)
    if (suppressSet) {
        log.info "Suppress list loaded: ${suppressSet.size()} ASMIDs will be skipped"
    }

    // ── Prediction pipeline ───────────────────────────────────────────────────
    def jobs = channel.fromPath(params.samples)
        .splitCsv(header: true)
        .map { row ->
            def species       = row.SPECIES?.trim()?.replaceAll(/['"]/, '')
            def strain        = row.STRAIN?.trim()?.replaceAll(/['"]/, '')
            strain = strain.replaceAll(/;.*$/, '').trim()
            def out           = [species, strain].findAll { it }.join('_').replaceAll(/\s+/, '_')
            def asmid         = row.ASMID?.trim()
            def locustag      = row.LOCUSTAG?.replaceAll(/[\r\n]/, '')?.trim()
            def busco         = row.BUSCO_LINEAGE?.trim()
            def header_length = 24
            def transl_table  = row.TRANSL_TABLE?.trim() ?: '1'
            def taxonid       = row.NCBI_TAXONID?.trim()
            tuple(out, asmid, species, strain, locustag, busco, header_length, transl_table, taxonid)
        }
        .filter { out, asmid, _sp, _st, _lt, _bl, _hl, _tt, _tid -> out && asmid }
        .take((params.n_test as int) > 0 ? params.n_test as int : -1)
        .filter { out, asmid, _sp, _st, _lt, _bl, _hl, _tt, _tid ->
            if (suppressSet.contains(asmid)) {
                log.info "Suppressing ${out} (asmid=${asmid})"
                return false
            }
            return true
        }
        .map { out, asmid, species, strain, locustag, busco, header_length, transl_table, taxonid ->
            def gz = file("${params.source}/${asmid}/${asmid}_genomic.fna.gz")
            tuple(out, asmid, species, strain, locustag, busco, header_length, transl_table, gz, taxonid)
        }
        .filter { out, asmid, _sp, _st, _lt, _bl, _hl, _tt, gz, _tid ->
            if (!gz.exists()) {
                log.warn "Missing genome for ${out} (asmid=${asmid}): ${gz}"
                return false
            }
            if (params.debug) {
                log.info "Queuing ${out}: genome=${gz} (${gz.size()} bytes)"
            }
            return true
        }

    if (params.debug) {
        jobs.view { t -> "[CHANNEL] Submitting: out=${t[0]}, asmid=${t[1]}, transl_table=${t[7]}, gz=${t[8]}" }
    }

    // Ensure taxondb is populated before any GENOME_CLEAN task starts.
    // SETUP_TAXONDB uses storeDir so it runs at most once across all pipeline runs.
    SETUP_TAXONDB()
    def taxondb_ch = SETUP_TAXONDB.out.ready.map { params.taxondb }
    GENOME_CLEAN(jobs.combine(taxondb_ch))

    if (!params.only_clean) {
        // Convert path output to absolute-path string so downstream val(genome_fa) processes
        // can reference the file directly without Nextflow re-staging it per-process.
        def clean_genome_ch = GENOME_CLEAN.out.genome
            .map { out, asmid, species, strain, locustag, busco, hlen, ttable, genome_fa, taxonid ->
                tuple(out, asmid, species, strain, locustag, busco, hlen, ttable,
                      genome_fa.toAbsolutePath().toString(), taxonid)
            }

        // ── Repeat masking ────────────────────────────────────────────────────────
        // predict_genome_ch carries the genome path to use for prediction — either
        // the soft-masked genome (default) or the clean unmasked genome (--skip_repeatmasker).
        def predict_genome_ch
        if (!params.skip_repeatmasker) {
            // Tag each assembly with species_tag for grouping and RM output naming.
            def tagged_ch = clean_genome_ch
                .map { out, asmid, species, strain, locustag, busco, hlen, ttable, genome_fa, taxonid ->
                    def species_tag = species.replaceAll(/\s+/, '_')
                    tuple(species_tag, out, asmid, species, strain, locustag, busco, hlen, ttable, genome_fa, taxonid)
                }

            // Run RepeatModeler once per species: group assemblies, pick the first.
            def rm_model_input = tagged_ch
                .map { species_tag, out, asmid, species, strain, locustag, busco, hlen, ttable, genome_fa, taxonid ->
                    tuple(species_tag, asmid, genome_fa)
                }
                .groupTuple(by: 0)
                .map { species_tag, asmids, genome_fas ->
                    tuple(species_tag, asmids[0], genome_fas[0])
                }
            REPEATMODELER_RUN(rm_model_input)

            // Join per-species RM library back to every assembly, then run RepeatMasker.
            REPEATMASKER_RUN(tagged_ch.join(REPEATMODELER_RUN.out.rmlib, by: 0))

            predict_genome_ch = REPEATMASKER_RUN.out.masked
                .map { out, asmid, species, strain, locustag, busco, hlen, ttable, masked_fa, taxonid ->
                    tuple(out, asmid, species, strain, locustag, busco, hlen, ttable,
                          masked_fa.toAbsolutePath().toString(), taxonid)
                }
        } else {
            // --skip_repeatmasker: use masked genome if a prior run produced it, else unmasked.
            predict_genome_ch = clean_genome_ch
                .map { out, asmid, species, strain, locustag, busco, hlen, ttable, genome_fa, taxonid ->
                    def masked = file("${launchDir}/input_clean_genomes/${asmid}.masked.fasta")
                    def use_fa = masked.exists() ? masked.toString() : genome_fa
                    if (params.debug) {
                        log.info "[DEBUG] ${asmid}: genome_fa=${use_fa} (masked=${masked.exists()})"
                    }
                    tuple(out, asmid, species, strain, locustag, busco, hlen, ttable, use_fa, taxonid)
                }
        }

        // FUNANNOTATE_PREDICT input tuple drops taxonid (not needed after masking/clean).
        // When SRA is enabled, TRAIN feeds PREDICT; otherwise PREDICT draws directly from
        // predict_genome_ch so it always runs regardless of RNA-seq availability.
        def predict_input_ch
        if (params.run_sra_fetch) {
            SRA_FETCH(predict_genome_ch)
            FUNANNOTATE_TRAIN(SRA_FETCH.out.reads)
            predict_input_ch = FUNANNOTATE_TRAIN.out
        } else {
            predict_input_ch = predict_genome_ch
                .map { out, asmid, species, strain, locustag, busco, hlen, ttable, genome_fa, _taxonid ->
                    tuple(out, asmid, species, strain, locustag, busco, hlen, ttable, genome_fa)
                }
        }

        def predict_ch = predict_input_ch
            .filter { out, _asmid, _sp, _st, _lt, _bl, _hl, _tt, _gfa ->
                !file("${params.target}/${out}/predict_results/${out}.gbk").exists()
            }
        FUNANNOTATE_PREDICT(predict_ch)

        // ── Post-predict steps ────────────────────────────────────────────────────
        // Independent channel parse so already-predicted species still reach these steps.
        def postpredict = channel.fromPath(params.samples)
            .splitCsv(header: true)
            .map { row ->
                def species       = row.SPECIES?.trim()?.replaceAll(/['"]/, '')
                def strain        = row.STRAIN?.trim()?.replaceAll(/['"]/, '')
                strain = strain.replaceAll(/;.*$/, '').trim()
                def out           = [species, strain].findAll { it }.join('_').replaceAll(/\s+/, '_')
                def asmid         = row.ASMID?.trim()
                def locustag      = row.LOCUSTAG?.replaceAll(/[\r\n]/, '')?.trim()
                def busco         = row.BUSCO_LINEAGE?.trim()
                def header_length = 24
                def transl_table  = row.TRANSL_TABLE?.trim() ?: '1'
                tuple(out, asmid, species, strain, locustag, busco, header_length, transl_table)
            }
            .filter { out, asmid, _sp, _st, _lt, _bl, _hl, _tt -> out && asmid }
            .take((params.n_test as int) > 0 ? params.n_test as int : -1)
            .filter { out, asmid, _sp, _st, _lt, _bl, _hl, _tt -> !suppressSet.contains(asmid) }
            .filter { out, _asmid, _sp, _st, _lt, _bl, _hl, _tt -> file("${params.target}/${out}/predict_results/${out}.gbk").exists() }

        if (params.run_antismash) {
            ANTISMASH_RUN(postpredict.filter { out, _asmid, _sp, _st, _lt, _bl, _hl, _tt ->
                def asDir = file("${params.target}/${out}/antismash_local")
                !(asDir.isDirectory() && asDir.list()?.any { it.endsWith('.json') || it.endsWith('.json.gz') })
            })
        }
        if (params.run_interpro) {
            INTERPROSCAN_RUN(postpredict.filter { out, _asmid, _sp, _st, _lt, _bl, _hl, _tt ->
                !file("${params.target}/${out}/annotate_misc/iprscan.xml").exists()
            })
        }
        if (params.run_signalp) {
            SIGNALP_RUN(postpredict.filter { out, _asmid, _sp, _st, _lt, _bl, _hl, _tt ->
                !file("${params.target}/${out}/annotate_misc/signalp.results.txt").exists()
            })
        }
    }
}
