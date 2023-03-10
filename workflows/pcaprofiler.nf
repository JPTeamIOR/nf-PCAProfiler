/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowPcaprofiler.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

ch_dummy_file = file("$projectDir/assets/dummy_file.txt", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'
include { PREPARE_GENOME } from '../subworkflows/local/prepare_genome'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { TRIMGALORE } from '../modules/nf-core/trimgalore/main'
include { STAR_ALIGN } from '../modules/nf-core/star/align/main'
include { SALMON_QUANT } from '../modules/nf-core/salmon/quant/main'
include { SORTMERNA       } from '../modules/nf-core/sortmerna/main'

include { IRFINDER_PROCESS_BAM } from '../modules/local/irfinder/process_bam'
include { IMOKA_EXTRACT_KMERS } from '../modules/local/imoka/extract_kmers'
include { DECONTAMINER } from '../modules/local/decontaminer/main'
include { WHIPPET_PROCESS_FASTQ } from '../modules/local/whippet/process_fastq/main.nf'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow PCAPROFILER {
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    //
    // SUBWORKFLOW: prepare genome
    //
    PREPARE_GENOME()
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)

    // Check rRNA databases for sortmerna

    ch_ribo_db = file(params.ribo_database_manifest, checkIfExists: true)
    if (ch_ribo_db.isEmpty()) { exit 1, "File provided with --ribo_database_manifest is empty: ${ch_ribo_db.getName() }!" }

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK(
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    //
    // MODULE: Run FastQC
    //
    FASTQC(
        INPUT_CHECK.out.reads
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // MODULE: Run Trimgalore
    //
    ch_filtered_reads = TRIMGALORE(INPUT_CHECK.out.reads).reads

    //
    // MODULE: Remove ribosomal RNA reads
    //

    ch_sortmerna_fastas = Channel.from(ch_ribo_db.readLines()).map { row -> file(row, checkIfExists: true) }.collect()
    SORTMERNA(
            ch_filtered_reads,
            ch_sortmerna_fastas
        )
        .reads
        .set { ch_filtered_reads }

    ch_versions = ch_versions.mix(SORTMERNA.out.versions.first())

    //
    // MODULE : iMOKA extract k-mers
    //
    IMOKA_EXTRACT_KMERS(ch_filtered_reads, params.k_len)
    ch_versions   = ch_versions.mix(IMOKA_EXTRACT_KMERS.out.versions.first())

    //
    // MODULE : Whippet quant
    //

    WHIPPET_PROCESS_FASTQ(ch_filtered_reads, PREPARE_GENOME.out.whippet_index)
    ch_versions   = ch_versions.mix(WHIPPET_PROCESS_FASTQ.out.versions.first())

    //
    // MODULE: STAR alignment
    //

    STAR_ALIGN(ch_filtered_reads,
        PREPARE_GENOME.out.star_index,
        PREPARE_GENOME.out.gtf,
        true,
        '',
        '')
    ch_orig_bam       = STAR_ALIGN.out.bam
    ch_log_final      = STAR_ALIGN.out.log_final
    ch_log_out        = STAR_ALIGN.out.log_out
    ch_log_progress   = STAR_ALIGN.out.log_progress
    ch_bam_sorted     = STAR_ALIGN.out.bam_sorted
    ch_bam_transcript = STAR_ALIGN.out.bam_transcript
    ch_fastq          = STAR_ALIGN.out.fastq
    ch_tab            = STAR_ALIGN.out.tab
    ch_versions       = ch_versions.mix(STAR_ALIGN.out.versions.first())

    ///
    /// MODULE: Decontaminer
    ///
    ch_human_ribo_db = params.human_ribo  ? file(params.human_ribo, checkIfExists: true) : file('NO_RIBO')
    ch_bacterial_db = params.bacterial_db ?  file(params.bacterial_db, checkIfExists: true) : file('NO_BACTERIAL')
    ch_fungi_db = params.fungi_db ? file(params.fungi_db, checkIfExists: true) : file('NO_FUNGI')
    ch_virus_db =  params.virus_db ? file(params.virus_db, checkIfExists: true) : file('NO_VIRUS')
    DECONTAMINER(ch_fastq, ch_human_ribo_db, ch_bacterial_db, ch_fungi_db, ch_virus_db, false)

    ///
    /// MODULE: Quantify using Salmon
    ///
    //
    SALMON_QUANT(
        ch_bam_transcript,
        ch_dummy_file,
        PREPARE_GENOME.out.gtf,
        PREPARE_GENOME.out.transcript_fasta,
        true,
        'A'
    )
    ch_versions = ch_versions.mix(SALMON_QUANT.out.versions)

    ///
    /// MODULE: IRFinder BAM
    ///

    IRFINDER_PROCESS_BAM(ch_orig_bam, PREPARE_GENOME.out.irfinder_ref)

    CUSTOM_DUMPSOFTWAREVERSIONS(
        ch_versions.unique { it.text }.collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowPcaprofiler.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowPcaprofiler.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect { it[1] }.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(TRIMGALORE.out.zip.collect { it[1] }.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(SORTMERNA.out.log.collect { it[1] }.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(STAR_ALIGN.out.log_final.collect { it[1] }.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(SALMON_QUANT.out.results.collect { it[1] }.ifEmpty([]))

    MULTIQC(
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
