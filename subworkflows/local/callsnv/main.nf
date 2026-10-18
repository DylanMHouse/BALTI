include { CLAIR3 } from '../../../modules/nf-core/clair3/main'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_CLAIR3 } from '../../../modules/nf-core/bcftools/view/main'


workflow CALLSNV {

    take:
    ch_input // channel
    ch_fasta
    ch_fai

    main:

    CLAIR3(
        ch_input,
        ch_fasta,
        ch_fai
    )
    ch_vcf = CLAIR3.out.vcf
    ch_tbi = CLAIR3.out.tbi



    def regions = params.regions
        ? channel.fromPath(params.regions, checkIfExists: true).first()
        : channel.value([])


    def targets = params.targets
        ? channel.fromPath(params.targets, checkIfExists: true).first()
        : channel.value([])

    def sample_names = params.sample_names
        ? channel.fromPath(params.sample_names, checkIfExists: true).first()
        : channel.value([])


    ch_clair3_output = ch_vcf.combine(ch_tbi)

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
