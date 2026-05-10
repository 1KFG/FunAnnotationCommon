#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.samples         = "${launchDir}/samples.csv"
params.target          = "${launchDir}/annotate"
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
params.max_cpus        = 32      // --max_cpus N: total CPUs for local executor

// Metadata tuple order used throughout:
//   val(out), val(asmid), val(species), val(strain), val(locustag),
//   val(busco_lineage), val(header_length), val(transl_table)
// GENOME_CLEAN additionally receives: path(genome_gz), val(taxonid)
// All steps from GENOME_CLEAN onward pass:  ..., path(genome_fa)

process GENOME_CLEAN {
    tag "$asmid"

    // Skip the task when the output file already exists in this directory.
    storeDir "${launchDir}/input_clean_genomes/clean"

    cpus   8
    memory '500 GB'
    time   '6h'

    input:
    tuple val(out), val(asmid), val(species), val(strain), val(locustag),
          val(busco_lineage), val(header_length), val(transl_table),
          path(genome_gz), val(taxonid)

    output:
    tuple val(out), val(asmid), val(species), val(strain), val(locustag),
          val(busco_lineage), val(header_length), val(transl_table),
          path("${asmid}.fa"), emit: genome
    path("${asmid}.purge.fasta"), emit: purge_fasta
    path("${asmid}.purge.fcs_gx-taxonomy.tsv"), emit: purge_tsv, optional: true

    script:
    """
    if [ ! -f "${genome_gz}" ]; then
        echo "ERROR: genome_gz not found at path: ${genome_gz}" >&2
        exit 1
    fi
    module load AAFTF

    # Ensure /dev/shm/gxdb is present on this node; register for cleanup when done.
    source ${launchDir}/scripts/setup_fcs_shm.sh

    echo "[INFO] Decompressing and cleaning genome for ${asmid}..."
    pigz -dc ${genome_gz} > \$SCRATCH/${asmid}.raw.fa
    AAFTF fcs_gx_purge --db /dev/shm/gxdb/all \
        -i \$SCRATCH/${asmid}.raw.fa --cpus ${task.cpus} \
        -o \$SCRATCH/${asmid}.purge.fasta \
        -t "${taxonid}" -w \$SCRATCH/${asmid}.fcs_report
    cp \$SCRATCH/${asmid}.purge.fasta .
    cp \$SCRATCH/${asmid}.purge.fcs_gx-taxonomy.tsv . 2>/dev/null || true
    cat \$SCRATCH/${asmid}.purge.fasta | \
        ${params.clean_script} --len ${params.min_contig_len} > ${asmid}.fa
    echo "[INFO] Clean genome written: ${asmid}.fa (\$(du -sh ${asmid}.fa | cut -f1))"
    """

    stub:
    """
    echo ">stub_${asmid}" > ${asmid}.fa
    touch ${asmid}.purge.fasta ${asmid}.purge.fcs_gx-taxonomy.tsv
    """
}

// Placeholder — implement RNAseq-based training evidence here when ready.
process FUNANNOTATE_TRAIN {
    tag "$out"

    input:
    tuple val(out), val(asmid), val(species), val(strain), val(locustag),
          val(busco_lineage), val(header_length), val(transl_table),
          path(genome_fa)

    output:
    tuple val(out), val(asmid), val(species), val(strain), val(locustag),
          val(busco_lineage), val(header_length), val(transl_table),
          path(genome_fa)

    script:
    """
    echo "[STUB] FUNANNOTATE_TRAIN noop for ${out}"
    """

    stub:
    """
    echo "[STUB] FUNANNOTATE_TRAIN stub for ${out}"
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
          path(genome_fa)

    output:
    tuple val(out), path("${out}/**")

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

    F=\$(ls ${out}/predict_results/*.gbk 2>/dev/null | head -n 1)
    if [ -z "\$F" ]; then
        echo "ERROR: funannotate predict did not produce a .gbk in ${out}/predict_results" >&2
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
    mkdir -p ${out}/predict_results
    touch ${out}/predict_results/${out}.gbk
    """
}

// TODO: interpro scan with nextflow sub?
// TODO: signalp on gpu?

process FUNANNOTATE_ANNOTATE {
    tag "$out"

    cpus   16
    memory '32 GB'
    time   '48h'

    publishDir "${params.target}", mode: 'copy', overwrite: true

    input:
    tuple val(out), val(asmid), val(species), val(strain), val(locustag), val(busco_lineage), val(header_length)

    output:
    tuple val(out), path("${out}/**")

    script:
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
        --cpu ${task.cpus} --tmpdir \$TMPDIR
    """

    stub:
    """
    echo "[STUB] Would run funannotate annotate for ${out}"
    mkdir -p ${out}/annotate_results
    touch ${out}/annotate_results/${out}.gbk
    """
}

// Returns true if predict_results already has a .gbk for this sample.
def hasExistingGbk(targetDir, out) {
    def dir = new File("${targetDir}/${out}/predict_results")
    if (!dir.exists()) return false
    return dir.list()?.any { f -> f.endsWith('.gbk') } ?: false
}

workflow {
    def target = file(params.target)

    def jobs = channel
        .fromPath(params.samples)
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
            [out, asmid, species, strain, locustag, busco, header_length, transl_table, taxonid]
        }
        .filter { out, asmid, _sp, _st, _lt, _bl, _hl, _tt, _tid ->
            out && asmid
        }
        // n_test > 0 limits to first N samples; -1 means take all
        .take((params.n_test as int) > 0 ? params.n_test as int : -1)
        .map { out, asmid, species, strain, locustag, busco, header_length, transl_table, taxonid ->
            def gz = file("${params.source}/${asmid}/${asmid}_genomic.fna.gz")
            // Reorder so genome_gz precedes taxonid, matching GENOME_CLEAN input declaration
            tuple(out, asmid, species, strain, locustag, busco, header_length, transl_table, gz, taxonid)
        }
        .filter { out, asmid, _sp, _st, _lt, _bl, _hl, _tt, gz, _tid ->
            if (hasExistingGbk(target, out)) {
                log.info "Skipping ${out}: predict_results gbk already present"
                return false
            }
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
    // FUNANNOTATE_TRAIN(GENOME_CLEAN.out.genome)
    // FUNANNOTATE_PREDICT(FUNANNOTATE_TRAIN.out)
}
