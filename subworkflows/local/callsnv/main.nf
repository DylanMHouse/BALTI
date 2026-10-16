include { CLAIR3 } from '../../../modules/nf-core/clair3/main'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_CLAIR3 } from '../../../modules/nf-core/clair3/main'


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
    ch_idx = CLAIR3.out.idx

    ch_clair3_output = ch_vcf.join(ch_idx, by: 0)
                            .map { meta, vcf, index ->
                                tuple(meta, vcf, index)
                            }

    BCFTOOLS_VIEW_CLAIR3(
        ch_clair3_output,
        params.regions ?: [],
        params.targets ?: [],
        params.sample_names ?: []
    )


    emit:

    vcf     = BCFTOOLS_VIEW_CLAIR3.out.vcf
    index   = BCFTOOLS_VIEW_CLAIR3.out.index
}
