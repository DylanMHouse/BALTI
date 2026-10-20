include { BCFTOOLS_SORT as BCFTOOLS_SORT_CPSR } from '../../../modules/nf-core/bcftools/sort/main'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_PASS } from '../../../modules/nf-core/bcftools/view/main'
include { CPSR_RUN                            } from '../../../modules/local/cpsr/main'


workflow CPSR {

    take:
    ch_vcf       // channel: [ val(meta), path(vcf) ] germline SNV calls
    ch_pcgr_dir  // channel: path(pcgr_refdata)
    ch_vep_cache // channel: path(vep_cache)

    main:

    // Sort into a bgzipped VCF with a .tbi (--temp-dir . is set by the module)
    BCFTOOLS_SORT_CPSR(ch_vcf)

    // Keep PASS variants; output is bgzipped and tabixed. bcftools view keeps the sorted order
    BCFTOOLS_VIEW_PASS(
        BCFTOOLS_SORT_CPSR.out.vcf.join(BCFTOOLS_SORT_CPSR.out.index, by: 0),
        [],
        [],
        []
    )
    ch_pass = BCFTOOLS_VIEW_PASS.out.vcf.join(BCFTOOLS_VIEW_PASS.out.index, by: 0)

    CPSR_RUN(
        ch_pass,
        ch_pcgr_dir,
        ch_vep_cache
    )

    emit:
    vcf     = BCFTOOLS_VIEW_PASS.out.vcf
    index   = BCFTOOLS_VIEW_PASS.out.index
    reports = CPSR_RUN.out.cpsr_reports
}
