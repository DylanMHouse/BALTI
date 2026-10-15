
// Bam alignment for unaligned bam input files from Dorado basecaller. Workflow will collate and convert bam files to fastq. Subsequent fastq and nanoplots will be made to assess initial quality

include { SAMTOOLS_CAT          } from '../../../modules/nf-core/samtools/cat/main'
include { SAMTOOLS_FASTQ        } from '../../../modules/nf-core/samtools/fastq/main'
include { FASTQC                } from '../../../modules/nf-core/fastqc/main'
include { NANOPLOT              } from '../../../modules/nf-core/nanoplot/main'


workflow PREPAREBAM {

    take:
    ch_bam // channel: [ val(meta), [ bam ] ]
    ch_fasta // channel: [ val(meta), [ ref ] ]

    main:

    ch_versions = channel.empty()

    SAMTOOLS_CAT(ch_bam)


    SAMTOOLS_FASTQ(SAMTOOLS_CAT.out.bam,
                    false)

    ch_fastq = SAMTOOLS_FASTQ.out.other



    // qc checks for fastq file
    FASTQC(ch_fastq)
    NANOPLOT(ch_fastq)



    emit:
    fastq    = ch_fastq                  // channel: [ val(meta), [ fastq ] ]
    png      = NANOPLOT.out.png          // channel: [ val(meta), [ bai ] ]
    zip      = FASTQC.out.zip            // channel: [ val(meta), [ csi ] ]
    versions = ch_versions
}
