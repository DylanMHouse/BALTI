#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/balti
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/balti
    Website: https://nf-co.re/balti
    Slack  : https://nfcore.slack.com/channels/balti
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BALTI  } from './workflows/balti'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_balti_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_balti_pipeline'
include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_balti_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// TODO nf-core: Remove this line if you don't need a FASTA file
//   This is an example of how to use getGenomeAttribute() to fetch parameters
//   from igenomes.config using `--genome`
params.fasta = getGenomeAttribute('fasta')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_BALTI {

    take:
    samplesheet // channel: samplesheet read in from --input

    main:

    //
    // Some defaults among all pipelines

    ch_samples = samplesheet.flatMap { meta, tumour, normal ->
    [ [ meta + [patient: meta.id, id: "${meta.id}_T", type: 'tumour'], tumour ],
      [ meta + [patient: meta.id, id: "${meta.id}_N", type: 'normal'], normal ] ]
}



    def valid_input_types = ['fastq', 'ubam', 'bam' ]
    if (!(params.input_type in valid_input_types)) {
        error("Invalid --input_type '${params.input_type}'. Valid options: ${valid_input_types.join(', ')}")
    }

    def valid_aligners = [
        'fastq': ['minimap_nogpu', 'minimap_gpu'],
	'ubam': ['minimap_nogpu', 'minimap_gpu'],
    ]

    if (valid_aligners.containsKey(params.input_type) && !(params.alignment_software in valid_aligners[params.input_type])) {
        error("--input_type ${params.input_type} requires alignment_software to be one of ${valid_aligners[params.input_type].join(', ')} (got '${params.alignment_software}')")
	}

    //
    // WORKFLOW: Run pipeline
    //
        BALTI (
            ch_samples,
            params.multiqc_config,
            params.multiqc_logo,
            params.multiqc_methods_description,
            params.outdir,
        )

    emit:
    multiqc_report = BALTI.out.multiqc_report // channel: /path/to/multiqc_report.html
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_BALTI (
        PIPELINE_INITIALISATION.out.samplesheet
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        NFCORE_BALTI.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
