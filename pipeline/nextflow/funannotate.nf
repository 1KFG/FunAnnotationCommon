#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.samples         = "${launchDir}/samples.csv"
params.target          = "${launchDir}/annotate"
params.source          = "/bigdata/stajichlab/shared/projects/1KFG/2021/NCBI_fungi/source/NCBI_ASM"
params.seqcenter       = "NCBI"
params.augustus_config = "${launchDir}/lib/augustus/3.5/config"
params.funannotate_db  = "/bigdata/stajichlab/shared/lib/funannotate_db"
params.min_contig_len  = 2000
params.clean_script    = "${launchDir}/scripts/clean_genome_fa.py"

process FUNANNOTATE_PREDICT {
    tag "$out"

    cpus   16
    memory '32 GB'
    time   '32h'

    publishDir "${params.target}", mode: 'copy', overwrite: true

    input:
    tuple val(out), val(asmid), val(species), val(strain), val(locustag), val(busco_lineage), path(genome_gz)

    output:
    tuple val(out), path("${out}/**")

    script:
    def tmpdir = System.getenv('SCRATCH') ?: '/tmp'
    """
    source /etc/profile.d/modules.sh 2>/dev/null || true
    [ -f \$HOME/.bashrc ] && source \$HOME/.bashrc
    module load funannotate

    export AUGUSTUS_CONFIG_PATH=${params.augustus_config}
    export FUNANNOTATE_DB=${params.funannotate_db}

    GENOME=${tmpdir}/${asmid}.fa
    pigz -dc ${genome_gz} | ${params.clean_script} --len ${params.min_contig_len} > \$GENOME

    funannotate predict --name ${locustag} -i \$GENOME --strain "${strain}" \\
        -o ${out} -s "${species}" --cpu ${task.cpus} --busco_db ${busco_lineage} \\
        --AUGUSTUS_CONFIG_PATH \$AUGUSTUS_CONFIG_PATH -w codingquarry:0 \\
        --min_training_models 30 --tmpdir ${tmpdir} --SeqCenter ${params.seqcenter} \\
        --keep_no_stops --header_length 24

    F=\$(ls ${out}/predict_results/*.gbk 2>/dev/null | head -n 1)
    if [ ! -z "\$F" ]; then
        rm -rf ${out}/predict_misc/EVM ${out}/predict_misc/proteins.combined.fa
        rm -rf ${out}/predict_misc/glimmerhmm
        rm -rf ${out}/predict_misc/busco
        rm -rf ${out}/predict_misc/busco_proteins
    fi

    rm -f \$GENOME
    """

    stub:
    """
    mkdir -p ${out}/predict_results
    touch ${out}/predict_results/${out}.gbk
    """
}

workflow {
    def target = file(params.target)

    Channel
        .fromPath(params.samples)
        .splitCsv(header: true)
        .map { row ->
            def species  = row.SPECIES?.trim()
            def out      = species?.replaceAll(/\s+/, '_')
            def asmid    = row.ASMID?.trim()
            def strain   = row.STRAIN?.trim()
            def locustag = row.LOCUSTAG?.replaceAll(/[\r\n]/, '')?.trim()
            def busco    = row.BUSCO_LINEAGE?.trim()
            [out, asmid, species, strain, locustag, busco]
        }
        .filter { out, asmid, species, strain, locustag, busco ->
            out && asmid
        }
        .map { out, asmid, species, strain, locustag, busco ->
            def gz       = file("${params.source}/${asmid}/${asmid}_genomic.fna.gz")
            def existing = file("${target}/${out}/predict_results/*.gbk")
            [out, asmid, species, strain, locustag, busco, gz, existing]
        }
        .filter { out, asmid, species, strain, locustag, busco, gz, existing ->
            if (existing) {
                log.info "Skipping ${out}: predict_results gbk already present"
                return false
            }
            if (!gz.exists()) {
                log.warn "Missing genome for ${out}: ${gz}"
                return false
            }
            return true
        }
        .map { out, asmid, species, strain, locustag, busco, gz, existing ->
            tuple(out, asmid, species, strain, locustag, busco, gz)
        }
        .set { jobs }

    FUNANNOTATE_PREDICT(jobs)
}
