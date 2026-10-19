// Type HLA, KIR, CYP and IG/TR genes with SpecImmune from reads pulled out of the haplotagged BAM

include { SPECIMMUNE_EXTRACTREADS } from '../../../modules/local/specimmune/extractreads/main'
include { SPECIMMUNE_TYPING       } from '../../../modules/local/specimmune/typing/main'

workflow IMMUNETYPING {

    take:
    ch_bam          // channel: [ val(meta), path(bam), path(bai) ] bai can be []; reads must be mapped to hg38 with alt (and HLA) contigs
    ch_fasta        // channel: [ val(meta), path(fasta) ] only read for CRAM input
    ch_gene_classes // channel: val(gene_class) any of HLA, KIR, CYP, IG_TR
    ch_db           // channel: path(db) folder built by SpecImmune's make_db.py
    ch_hg38         // channel: [ path(fasta), path(index files) ] no-alt hg38 for CYP and IG_TR, can be []

    main:

    // One extraction and one typing run per sample and gene class
    SPECIMMUNE_EXTRACTREADS(
        ch_bam.combine(ch_gene_classes),
        ch_fasta
    )

    SPECIMMUNE_TYPING(
        SPECIMMUNE_EXTRACTREADS.out.fastq,
        ch_db,
        ch_hg38
    )

    emit:
    result = SPECIMMUNE_TYPING.out.result  // channel: [ val(meta), val(gene_class), path(txt) ]
    outdir = SPECIMMUNE_TYPING.out.outdir  // channel: [ val(meta), val(gene_class), path(dir) ] <gene_class>/<meta.id>/
    fastq  = SPECIMMUNE_EXTRACTREADS.out.fastq // channel: [ val(meta), val(gene_class), path(fastq.gz) ]
}
