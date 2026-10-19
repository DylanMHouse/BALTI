// Co-phase SNVs, SVs and 5mC modifications with LongPhase, then tag reads by haplotype (HP/PS tags)

include { LONGPHASE_MODCALL  } from '../../../modules/local/longphase/modcall/main'
include { LONGPHASE_PHASE    } from '../../../modules/nf-core/longphase/phase/main'
include { LONGPHASE_HAPLOTAG } from '../../../modules/nf-core/longphase/haplotag/main'

workflow PHASEHAPLOTAG {

    take:
    ch_input // channel: [ val(meta), path(bam), path(bai), path(snvs), path(svs) ] svs can be []
    ch_fasta // channel: [ val(meta), path(fasta) ]
    ch_fai   // channel: [ val(meta), path(fai) ]

    main:

    ch_bam = ch_input.map { meta, bam, bai, _snvs, _svs -> [ meta, bam, bai ] }
    ch_fasta_fai = ch_fasta.join(ch_fai, by: 0) // [ meta, fasta, fai ]

    // Modifications come from the MM/ML tags carried through from the Dorado BAM
    LONGPHASE_MODCALL(
        ch_bam,
        ch_fasta_fai
    )

    LONGPHASE_PHASE(
        ch_input.join(LONGPHASE_MODCALL.out.vcf, by: 0),
        ch_fasta,
        ch_fai
    )

    // sv_vcf and mod_vcf are optional outputs, so fill any missing one with [] for haplotag
    ch_haplotag_input = ch_bam
        .join(LONGPHASE_PHASE.out.snv_vcf, by: 0)
        .join(LONGPHASE_PHASE.out.sv_vcf, by: 0, remainder: true)
        .join(LONGPHASE_PHASE.out.mod_vcf, by: 0, remainder: true)
        .map { meta, bam, bai, snvs, svs, mods -> [ meta, bam, bai, snvs, svs ?: [], mods ?: [] ] }

    LONGPHASE_HAPLOTAG(
        ch_haplotag_input,
        ch_fasta,
        ch_fai
    )

    emit:
    bam         = LONGPHASE_HAPLOTAG.out.bam     // channel: [ val(meta), path(bam) ]
    snv_vcf     = LONGPHASE_PHASE.out.snv_vcf    // channel: [ val(meta), path(vcf.gz) ]
    sv_vcf      = LONGPHASE_PHASE.out.sv_vcf     // channel: [ val(meta), path(vcf.gz) ]
    mod_vcf     = LONGPHASE_PHASE.out.mod_vcf    // channel: [ val(meta), path(vcf.gz) ]
    modcall_vcf = LONGPHASE_MODCALL.out.vcf      // channel: [ val(meta), path(vcf.gz) ]
}
