process SPECIMMUNE_TYPING {
    tag "${meta.id}:${gene_class}"
    label 'process_medium'

    container 'ghcr.io/dylanmhouse/specimmune:0.0.3'

    input:
    tuple val(meta), val(gene_class), path(fastq)
    path(db)
    path(hg38, stageAs: 'hg38/*') // no-alt hg38 fasta with its .fai and bwa index files; only used for CYP and IG_TR, can be []

    output:
    tuple val(meta), val(gene_class), path("${gene_class}/${prefix}/${result}"), emit: result
    tuple val(meta), val(gene_class), path("${gene_class}/")                    , emit: outdir
    tuple val("${task.process}"), val('specimmune'), val(VERSION)              , emit: versions_specimmune, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // main.py --version reports 1.0.0 for every release, so this follows the image tag
    VERSION    = '0.0.3'
    def args   = task.ext.args ?: ''
    prefix     = task.ext.prefix ?: "${meta.id}"
    result     = resultName(prefix, gene_class)
    def hg38_fasta = (hg38 instanceof List ? hg38 : [hg38]).find { f -> f && f.name ==~ /.*\.(fa|fasta|fna)$/ }
    def hg38_arg   = hg38_fasta ? "--hg38 ${hg38_fasta}" : ''
    """
    specimmune \\
        -r ${fastq} \\
        -n ${prefix} \\
        -i ${gene_class} \\
        -o ${gene_class} \\
        -j ${task.cpus} \\
        --db ${db} \\
        ${hg38_arg} \\
        ${args}
    """

    stub:
    VERSION  = '0.0.3'
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    result   = resultName(prefix, gene_class)
    """
    echo $args
    mkdir -p ${gene_class}/${prefix}
    touch ${gene_class}/${prefix}/${result}
    """
}

def resultName(prefix, gene_class) {
    gene_class == 'CYP'   ? "${prefix}.CYP.merge.type.result.txt"
    : gene_class == 'IG_TR' ? "${prefix}.IG_TR_typing_result.txt"
    : "${prefix}.${gene_class}.final.type.result.formatted.txt"
}
