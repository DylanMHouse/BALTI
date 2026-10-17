
// Bam alignment for unaligned bam input files from Dorado basecaller. Workflow will collate and convert bam files to fastq. Subsequent fastqc and nanoplots will be made to assess initial quality

include { SAMTOOLS_CAT          } from '../../../modules/nf-core/samtools/cat/main'
include { SAMTOOLS_FASTQ        } from '../../../modules/nf-core/samtools/fastq/main'
include { FASTQC                } from '../../../modules/nf-core/fastqc/main'
include { NANOPLOT              } from '../../../modules/nf-core/nanoplot/main'


workflow PREPAREBAM {

    take:
    ch_bam // channel: [ val(meta), [ bam ] ]

    main:

    ch_versions = channel.empty()
    ch_png      = channel.empty()
    ch_txt      = channel.empty()
    ch_zip      = channel.empty()

    SAMTOOLS_CAT(ch_bam)


    SAMTOOLS_FASTQ(SAMTOOLS_CAT.out.bam,
                    false)

    ch_fastq = SAMTOOLS_FASTQ.out.other



    // qc checks for fastq file
    if (params.qc) {
        FASTQC(ch_fastq)
        NANOPLOT(ch_fastq)
        ch_png = NANOPLOT.out.png
        ch_txt = NANOPLOT.out.txt
        ch_zip = FASTQC.out.zip
    }



    emit:
    fastq    = ch_fastq                  // channel: [ val(meta), [ fastq ] ]
    png      = ch_png                    // channel: [ val(meta), [ png ] ]
    txt      = ch_txt                    // channel: [ val(meta), [ txt ] ]
    zip      = ch_zip                    // channel: [ val(meta), [ zip ] ]
    versions = ch_versions
}
