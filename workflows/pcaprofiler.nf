/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowPcaprofiler.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.ref_dir ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

if (params.conatminant_method ==~ /all|kraken2|decontaminer/) {
    error 'Contaminant method must be one of all, kraken2, decontaminer'
}

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
include { CHECK_REFERENCE } from '../subworkflows/local/check_reference'

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
include { KRAKEN2_KRAKEN2 } from '../modules/nf-core/kraken2/kraken2/main'
include { BRACKEN_BRACKEN } from '../modules/nf-core/bracken/bracken/main'

include { IRFINDER_PROCESS_BAM } from '../modules/local/irfinder/process_bam'
include { IMOKA_EXTRACT_KMERS } from '../modules/local/imoka/extract_kmers'
include { DECONTAMINER } from '../modules/local/decontaminer/main'
include { WHIPPET_PROCESS_FASTQ } from '../modules/local/whippet/process_fastq/main.nf'
include { STARFUSION_PROCESS } from '../modules/local/star-fusion/main.nf'
include { CTAT_MUTATION_PROCESS_BAM } from '../modules/local/ctat-mutation/main.nf'

include { SAMTOOLS_INDEX } from '../modules/nf-core/samtools/index/main'
include { SAMTOOLS_SORT } from '../modules/nf-core/samtools/sort/main'

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
    CHECK_REFERENCE(params.ref_dir)
    ch_versions = ch_versions.mix(CHECK_REFERENCE.out.versions)
    ch_sortmerna_fastas = CHECK_REFERENCE.out.sortmerna_fastas
    ch_whippet_index = CHECK_REFERENCE.out.whippet_index
    ch_star_index = CHECK_REFERENCE.out.star_index
    ch_gtf = CHECK_REFERENCE.out.gtf
    ch_decontaminer_ref = CHECK_REFERENCE.out.decontaminer
    ch_ctat_ref = CHECK_REFERENCE.out.ctat
    ch_transcript_fasta = CHECK_REFERENCE.out.transcript_fasta
    ch_irfinder_ref = CHECK_REFERENCE.out.irfinder_ref
    ch_kraken2_ref = CHECK_REFERENCE.out.kraken2

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
    if (params.extract_kmers) {
        IMOKA_EXTRACT_KMERS(ch_filtered_reads, params.k_len)
        ch_versions   = ch_versions.mix(IMOKA_EXTRACT_KMERS.out.versions.first())
    }

    //
    // MODULE : Whippet quant
    //

    WHIPPET_PROCESS_FASTQ(ch_filtered_reads, ch_whippet_index)
    ch_versions   = ch_versions.mix(WHIPPET_PROCESS_FASTQ.out.versions.first())

    //
    // MODULE: STAR alignment
    //

    STAR_ALIGN(ch_filtered_reads,
        ch_star_index,
        ch_gtf,
        true,
        '',
        '')
    ch_orig_bam       = STAR_ALIGN.out.bam
    ch_bam_transcript = STAR_ALIGN.out.bam_transcript
    ch_unmapped_fastq = STAR_ALIGN.out.fastq
    ch_junction       = STAR_ALIGN.out.junction

    ch_versions       = ch_versions.mix(STAR_ALIGN.out.versions.first())

    ///
    /// MODULE: Decontaminer
    ///
    if (params.conatminant_method == 'all' || params.conatminant_method == 'decontaminer') {
        ch_human_ribo_db = file("${ch_decontaminer_ref}/HUMAN_RNA").exists  ? file("${ch_decontaminer_ref}/HUMAN_RNA") : file('NO_RIBO')
        ch_bacterial_db = file("${ch_decontaminer_ref}/HUMAN_RNA").exists  ? file("${ch_decontaminer_ref}/BACTERIA") : file('NO_BACTERIAL')
        ch_fungi_db = file("${ch_decontaminer_ref}/HUMAN_RNA").exists  ? file("${ch_decontaminer_ref}/FUNGI") : file('NO_FUNGI')
        ch_virus_db =  file("${ch_decontaminer_ref}/HUMAN_RNA").exists  ? file("${ch_decontaminer_ref}/VIRUSES") : file('NO_VIRUS')
        DECONTAMINER(ch_unmapped_fastq, ch_human_ribo_db, ch_bacterial_db, ch_fungi_db, ch_virus_db, false)
        ch_versions = ch_versions.mix(DECONTAMINER.out.versions)
    }

    if (params.conatminant_method == 'all' || params.conatminant_method == 'kraken2') {
        KRAKEN2_KRAKEN2(ch_unmapped_fastq, ch_kraken2_ref, false, false)
        BRACKEN_BRACKEN(KRAKEN2_KRAKEN2.out.report, ch_kraken2_ref)
        ch_versions = ch_versions.mix(BRACKEN_BRACKEN.out.versions.first())
        ch_versions = ch_versions.mix(KRAKEN2_KRAKEN2.out.versions.first())
    }

    ///
    /// MODULE: STAR-Fusion
    ///

    STARFUSION_PROCESS(
        ch_junction,
        ch_ctat_ref
    )
    ch_version = ch_version.mix(STARFUSION_PROCESS.out.versions)

    ///
    /// MODULE: Quantify using Salmon
    ///

    SALMON_QUANT(
        ch_bam_transcript,
        ch_dummy_file,
        ch_gtf,
        ch_transcript_fasta,
        true,
        'A'
    )
    ch_versions = ch_versions.mix(SALMON_QUANT.out.versions)

    ///
    /// MODULE: IRFinder BAM
    ///

    IRFINDER_PROCESS_BAM(ch_orig_bam, ch_irfinder_ref)
    ch_versions = ch_versions.mix(IRFINDER_PROCESS_BAM.out.versions)

    ///
    /// MODULE: CTAT-Mutation
    ///
    SAMTOOLS_SORT(ch_orig_bam)
    ch_sorted_bam = SAMTOOLS_SORT.out.bam
    SAMTOOLS_INDEX(ch_sorted_bam)
    ch_sorted_bam_bai = SAMTOOLS_INDEX.out.bai

    CTAT_MUTATION_PROCESS_BAM(ch_sorted_bam.join(ch_sorted_bam_bai) , ch_ctat_ref)

    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)
    ch_versions = ch_versions.mix(CTAT_MUTATION_PROCESS_BAM.out.versions)

    ///

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
