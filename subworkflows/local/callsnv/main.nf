include { CLAIR3 } from '../../../modules/nf-core/clair3/main'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_CLAIR3 } from '../../../modules/nf-core/bcftools/view/main'


workflow CALLSNV {

    take:
    ch_input // channel: [ val(meta), [ bam ], [bai], [ model ], [ user model ] [ platform ]]
    ch_ref
    ch_ref_index

    main:

    ch_clair3_input = ch_input.map { meta, bam, bai ->
        tuple(meta, bam, bai, params.clair3_model, [], 'ont')
        }

    CLAIR3(
        ch_clair3_input,
        ch_ref,
        ch_ref_index
    )
    ch_vcf = CLAIR3.out.vcf
    ch_tbi = CLAIR3.out.tbi

    def regions = params.regions ?: []
    def targets = params.targets ?: []
    def sample_names = params.sample_names ?: []

    ch_clair3_output = ch_vcf.join(ch_tbi, by: 0)
                            .map { meta, vcf, index ->
                                tuple(meta, vcf, index)
                            }

    BCFTOOLS_VIEW_CLAIR3(
        ch_clair3_output,
        regions,
        targets,
        sample_names
    )


    emit:

    vcf     = BCFTOOLS_VIEW_CLAIR3.out.vcf
    index   = BCFTOOLS_VIEW_CLAIR3.out.index
}
