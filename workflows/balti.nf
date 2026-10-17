/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap            } from 'plugin/nf-schema'

include { MOSDEPTH                    } from '../modules/nf-core/mosdepth/main'
include { BAM_STATS_SAMTOOLS          } from '../subworkflows/nf-core/bam_stats_samtools'

include { SAMTOOLS_FAIDX              } from '../modules/nf-core/samtools/faidx/main'
include { ALIGNMENT		      } from '../subworkflows/local/alignment'
include { PREPAREBAM                  } from '../subworkflows/local/preparebam'


include { paramsSummaryMultiqc        } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML      } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText      } from '../subworkflows/local/utils_nfcore_balti_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow BALTI {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:



    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()
    def ch_multiqc_report = channel.empty()
    def ch_aligned = channel.empty()
    def ch_aligned_index = channel.empty()

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name: 'nf_core_'  +  'balti_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )


    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    			REFERENCES
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // --fasta takes priority over --genome if provided of course.
    def fasta_path = params.fasta
        ?: params.genome
            ? params.genomes.containsKey(params.genome)
                ? params.genomes[params.genome].fasta
                : error("Genome '${params.genome}' not found in iGenomes config. Check --genome or set --igenomes_ignore false")
            : error("Provide a reference with --fasta or --genome")

    def ch_fasta = channel.fromPath(fasta_path, checkIfExists: true)
        .map { fasta -> tuple([id: "ref"], fasta) }
        .first()

    // Indexing the reference to .fai is necessary here as it is mandatory for some steps.
    SAMTOOLS_FAIDX ( ch_fasta.map{ meta, fasta -> [meta, fasta, [] ] }, false)
    def ch_fasta_fai = ch_fasta
        .join(SAMTOOLS_FAIDX.out.fai, by: 0)
        .first()

    // Optional inputs: an empty list stands in for a missing file.
    def ch_intervals = params.intervals
        ? channel.fromPath(params.intervals, checkIfExists: true).first()
        : channel.value([])

    // minimap2 accepts a prebuilt .mmi in place of the fasta; parabricks always needs the fasta
    def ch_minimap2_ref = params.minimap2_index && params.alignment_software == 'minimap_nogpu'
        ? channel.fromPath(params.minimap2_index, checkIfExists: true).map { mmi -> tuple([id: "ref"], mmi) }.first()
        : ch_fasta

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    			READ QC
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    if (params.qc && params.input_type == 'fastq') {
        FASTQC(ch_samplesheet)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.map{ _meta, file -> file })
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    			INDEXING,  ALIGNMENT & QC CHECKS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    if (params.minimap2) {

        if (params.input_type == 'fastq') {

            ALIGNMENT(ch_samplesheet, ch_minimap2_ref)

        } else if (params.input_type == 'bam') {
            PREPAREBAM(ch_samplesheet)
            ch_multiqc_files = ch_multiqc_files.mix(PREPAREBAM.out.zip.map{ _meta, file -> file })
            ch_multiqc_files = ch_multiqc_files.mix(PREPAREBAM.out.txt.map{ _meta, file -> file })

            ALIGNMENT(PREPAREBAM.out.fastq, ch_minimap2_ref)
        }

        // Keep the full meta (id, patient, type) so tumour and germline can be told apart downstream. I've had this issue in previous commits :(
        ch_aligned       = ALIGNMENT.out.bam
        ch_aligned_index = ALIGNMENT.out.bai

        def ch_input_bam = ch_aligned.join(ch_aligned_index, by: 0)

        // Generate some statistics of aligned reads
        BAM_STATS_SAMTOOLS(
            ch_input_bam,
            ch_fasta_fai
        )
        ch_multiqc_files = ch_multiqc_files
            .mix(BAM_STATS_SAMTOOLS.out.stats.map{ _meta, file -> file })
            .mix(BAM_STATS_SAMTOOLS.out.flagstat.map{ _meta, file -> file })
            .mix(BAM_STATS_SAMTOOLS.out.idxstats.map{ _meta, file -> file })

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                COVERAGE
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */

        if (params.coverage) {
            MOSDEPTH(
                // wrap the bed so an empty [] still fills the bed slot: [ meta, bam, bai, bed ]
                ch_input_bam.combine(ch_intervals.map { bed -> [bed] }),
                ch_fasta,
                []
            )
            ch_multiqc_files = ch_multiqc_files
                .mix(MOSDEPTH.out.global_txt.map{ _meta, file -> file })
                .mix(MOSDEPTH.out.summary_txt.map{ _meta, file -> file })
        }
    }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))

    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'balti'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )
    ch_multiqc_report = MULTIQC.out.report.map { _meta, report -> report }

    emit:
    multiqc_report = ch_multiqc_report
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
