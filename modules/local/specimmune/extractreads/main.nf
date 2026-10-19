process SPECIMMUNE_EXTRACTREADS {
    tag "${meta.id}:${gene_class}"
    label 'process_low'

    // SpecImmune is not on bioconda, so there is no Seqera Containers image; this one is built from ../Dockerfile
    container 'ghcr.io/dylanmhouse/specimmune:0.0.3'

    input:
    tuple val(meta), path(bam), path(bai), val(gene_class)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), val(gene_class), path("${prefix}.${gene_class}.fastq.gz"), emit: fastq, optional: true
    tuple val(meta), val(gene_class), path("${prefix}.${gene_class}.bam")     , emit: bam
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d; s/samtools //'"), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "SPECIMMUNE_EXTRACTREADS does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args      = task.ext.args ?: ''
    prefix        = task.ext.prefix ?: "${meta.id}"
    def reference = bam.name.endsWith('.cram') ? "-r ${fasta}" : ''
    def index     = bai ? '' : "samtools index -@ ${task.cpus} ${bam}"
    """
    ${index}

    ExtractReads.sh \\
        -s ${prefix} \\
        -i ${bam} \\
        -g ${gene_class} \\
        -o . \\
        ${reference} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.${gene_class}.fastq.gz
    touch ${prefix}.${gene_class}.bam
    """
}
