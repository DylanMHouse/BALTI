// Structural variant calling: Sniffles for SVs from split/spanning reads, Spectre for large CNVs (>100kb) from binned coverage

include { SNIFFLES                      } from '../../../modules/nf-core/sniffles/main'
include { MOSDEPTH as MOSDEPTH_SPECTRE  } from '../../../modules/nf-core/mosdepth/main'
include { SPECTRE_CNVCALLER             } from '../../../modules/local/spectre/cnvcaller/main'

workflow CALLSV {

    take:
    ch_input // channel: [ val(meta), path(bam), path(bai) ]
    ch_snv   // channel: [ val(meta), path(vcf) ] small variant calls for the same samples, used by Spectre
    ch_fasta // channel: [ val(meta), path(fasta) ]
    ch_spectre_metadata  // channel: path(mdr) or [] (built-in metadata is set by name in conf/modules.config)
    ch_spectre_blacklist // channel: path(bed) or [] (built-in blacklist is set by name in conf/modules.config)

    main:

    SNIFFLES(
        ch_input.map { meta, bam, bai -> [ meta, bam, bai, [] ] },
        ch_fasta,
        [ [:], [] ],
        true,
        false
    )

    // Spectre needs binned coverage; the window size is set with ext.args (--by) in conf/modules.config
    MOSDEPTH_SPECTRE(
        ch_input.map { meta, bam, bai -> [ meta, bam, bai, [] ] },
        ch_fasta,
        []
    )

    ch_spectre_input = MOSDEPTH_SPECTRE.out.summary_txt
        .join(MOSDEPTH_SPECTRE.out.regions_bed, by: 0)
        .join(MOSDEPTH_SPECTRE.out.regions_csi, by: 0)
        .join(ch_snv, by: 0)

    SPECTRE_CNVCALLER(
        ch_spectre_input,
        ch_fasta,
        ch_spectre_metadata,
        ch_spectre_blacklist
    )

    emit:
    sv_vcf    = SNIFFLES.out.vcf                  // channel: [ val(meta), path(vcf.gz) ]
    sv_tbi    = SNIFFLES.out.tbi                  // channel: [ val(meta), path(vcf.gz.tbi) ]
    cnv_vcf   = SPECTRE_CNVCALLER.out.vcf         // channel: [ val(meta), path(vcf) ]
    cnv_bed   = SPECTRE_CNVCALLER.out.bed         // channel: [ val(meta), path(bed) ]
    karyotype = SPECTRE_CNVCALLER.out.karyotype   // channel: [ val(meta), path(txt) ]
}
