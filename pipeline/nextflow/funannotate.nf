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


// Metadata tuple order used throughout:
//   val(out), val(asmid), val(species), val(strain), val(locustag),
//   val(busco_lineage), val(header_length), val(transl_table)
// GENOME_CLEAN receives: ..., path(genome_gz), val(taxonid)
//   → emits: ..., path(genome_fa), val(taxonid)
//   → writes <asmid>.fa to input_clean_genomes/ (storeDir; skip check targets this file)
//   → purge/FCS intermediates written as side effects to input_clean_genomes/clean/
// SRA_FETCH receives: ..., path(genome_fa), val(taxonid)
//   → emits: ..., val(genome_fa_abs), path(reads_dir)   [reads_dir may be empty]
// FUNANNOTATE_TRAIN receives: ..., val(genome_fa), path(reads_dir)
//   → emits: ..., val(genome_fa)   [reads deleted after training]
// FUNANNOTATE_PREDICT receives: ..., val(genome_fa)

process GENOME_CLEAN {
    tag "$asmid"

    // Nextflow skips this task when input_clean_genomes/<asmid>.fa already exists.
    storeDir "${launchDir}/input_clean_genomes"

    cpus   16
    memory '450 GB'
    time   '6h'

    input:
    tuple val(out), val(asmid), val(species), val(strain), val(locustag),
          val(busco_lineage), val(header_length), val(transl_table),
          path(genome_gz), val(taxonid)

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
    TAXONKIT_DB=${params.taxondb}
    module load taxonkit
    phylum=\$(echo ${taxonid} | taxonkit --data-dir \$TAXONKIT_DB lineage | taxonkit --data-dir \$TAXONKIT_DB reformat -f "{p}" | cut -f3 | taxonkit --data-dir \$TAXONKIT_DB name2taxid | cut -f2)
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
    pigz \$SCRATCH/${asmid}.purge.fasta  \$SCRATCH/${asmid}.purge.fcs_gx-taxonomy.tsv 
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
          path(genome_fa), val(taxonid)

    output:
    tuple val(out), val(asmid), val(species), val(strain), val(locustag),
          val(busco_lineage), val(header_length), val(transl_table),
          val("${launchDir}/input_clean_genomes/${asmid}.fa"), path("reads"), emit: reads

    script:
    """
    mkdir -p reads
    module load entrez-direct
    module load sratoolkit

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
        parallel-fastq-dump --sra-id \$ACC --threads ${task.cpus} \\
            --outdir reads/ --split-files --gzip --tmpdir \$TMPDIR || {
            echo "[WARN] Download failed for \$ACC, skipping"
            rm -f reads/\${ACC}*.fastq.gz
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
    memory '64 GB'
    time   '24h'

    input:
    tuple val(out), val(asmid), val(species), val(strain), val(locustag),
          val(busco_lineage), val(header_length), val(transl_table),
          val(genome_fa), path(reads_dir)

    output:
    tuple val(out), val(asmid), val(species), val(strain), val(locustag),
          val(busco_lineage), val(header_length), val(transl_table),
          val(genome_fa)

    script:
    """
    R1=(\$(ls ${reads_dir}/*_1.fastq.gz 2>/dev/null || true))

    if [ \${#R1[@]} -eq 0 ]; then
        echo "[INFO] No RNAseq reads for ${out}, skipping funannotate train"
        exit 0
    fi

    R2=(\$(ls ${reads_dir}/*_2.fastq.gz 2>/dev/null || true))
    if [ \${#R1[@]} -ne \${#R2[@]} ]; then
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

    echo "[INFO] Running funannotate train for ${out} with \${#R1[@]} read pair(s)"

    funannotate train -i ${genome_fa} -o ${params.target}/${out} \\
        --left \${R1[@]} --right \${R2[@]} \\
        --species "${species}" --strain "${strain}" \\
        --cpus ${task.cpus} --memory "${task.memory}" \\
        --tmpdir \$TMPDIR

    # Resolve and delete the original fastq files in the SRA_FETCH work directory
    # to reclaim disk immediately; scratch dir itself is auto-cleaned by Nextflow.
    REAL_READS=\$(readlink -f ${reads_dir})
    if [ "\$REAL_READS" != "${reads_dir}" ]; then
        rm -rf "\$REAL_READS"
        echo "[INFO] Removed SRA reads at \$REAL_READS"
    fi
    """

    stub:
    """
    echo "[STUB] FUNANNOTATE_TRAIN stub for ${out} (reads_dir=${reads_dir})"
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
        --keep_no_stops --header_length 24 --protein_evidence ${params.proteins} \\
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


    GENOME_CLEAN(jobs)

    if (!params.only_clean) {
        SRA_FETCH(GENOME_CLEAN.out.genome)
        FUNANNOTATE_TRAIN(SRA_FETCH.out.reads)
        def predict_ch = FUNANNOTATE_TRAIN.out
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
