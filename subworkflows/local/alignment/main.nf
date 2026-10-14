include { MINIMAP2_ALIGN          } from '../../../modules/nf-core/minimap2/align/main'
include { PARABRICKS_MINIMAP2     } from '../../../modules/nf-core/parabricks/minimap2/main'

workflow ALIGNMENT {

    take:
    ch_reads // channel: [ val(meta), [ fast(a).gz ] ]
    ch_fasta // channel: [ val(meta), [ reference ] ]

    main:

        if (params.alignment_software == 'minimap_nogpu') {
	    MINIMAP2_ALIGN(
	        ch_reads,
		ch_fasta,
		true,
		'bai',
		[],
		[]
	    )
	ch_bam      = MINIMAP2_ALIGN.out.bam
	ch_bai      = MINIMAP2_ALIGN.out.index
	ch_versions = MINIMAP2_ALIGN.out.versions_minimap2
	} else if (params.alignment_software == 'minimap_gpu') {
	    PARABRICKS_MINIMAP2(
	        ch_reads,
		ch_fasta,
		[[:], []],
		[[:], []],
		'bam'
	    )
	    ch_bam = PARABRICKS_MINIMAP2.out.bam
	    ch_bai = PARABRICKS_MINIMAP2.out.bai
	    ch_versions = PARABRICKS_MINIMAP2.out.versions_parabricks
	} else {
	    error("Unknown --alignment_software '${params.alignment_software}'. Please choose either minimap_nogpu, or minimap_gpu")
	}


    emit:
    bam           = ch_bam               // channel: [ val(meta), [ bam ] ]
    bai           = ch_bai          	 // channel: [ val(meta), [ bai ] ]
    versions      = ch_versions		 // channel: [ path(versions.yml) ]
}
