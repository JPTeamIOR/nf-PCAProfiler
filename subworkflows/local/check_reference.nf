
//
// Check input reference and get channels
//

include { GTF_GENE_FILTER } from '../../modules/local/gtf_gene_filter'
include { GFFREAD_GTF_TO_FASTA } from '../../modules/local/gffread_gtf_to_fasta'

workflow CHECK_REFERENCE {
    take:
        ref_dir
    main:

    // Check few files to ensure that they are all present
    file("$projectDir/assets/ref_check.txt").readLines().each { v ->
        if (! file("${ref_dir}/${v}").exists()) {
            error "The file ${ref_dir}/${v} in the reference does not exists."
        }
    }

    ch_sortmerna_fastas = files("${ref_dir}/rRNA/*.fasta")

    ch_irfinder_ref = file("${ref_dir}/IRFinder/", checkIfExists : true)
    ch_ctat_ref = file("${ref_dir}/ctat_genome_lib_build_dir/", checkIfExists : true)

    ch_decontaminer_human_ribo_db = Channel.empty()
    ch_decontaminer_bacterial_db = Channel.empty()
    ch_decontaminer_fungi_db = Channel.empty()
    ch_decontaminer_virus_db = Channel.empty()
    if ( params.contaminant_method == 'all' ||  params.contaminant_method == 'decontaminer' ) {
        base_decontaminer_path = "${ref_dir}/Decontaminer"
        if ( ! file(base_decontaminer_path).exists() ){
            error "No Decontaminer found in the reference directory."
        }
        ch_decontaminer_human_ribo_db = file("${base_decontaminer_path}/HUMAN_RNA").exists()  ? file("${base_decontaminer_path}/HUMAN_RNA") : file('NO_RIBO')
        ch_decontaminer_bacterial_db = file("${base_decontaminer_path}/BACTERIA").exists()  ? file("${base_decontaminer_path}/BACTERIA") : file('NO_BACTERIAL')
        ch_decontaminer_fungi_db = file("${base_decontaminer_path}/FUNGI").exists()  ? file("${base_decontaminer_path}/FUNGI") : file('NO_FUNGI')
        ch_decontaminer_virus_db = file("${base_decontaminer_path}/VIRUSES").exists()  ? file("${base_decontaminer_path}/VIRUSES") : file('NO_VIRUS')
    }

    ch_star_ref = file("${ref_dir}/ctat_genome_lib_build_dir/ref_genome.fa.star.idx", checkIfExists : true)
    ch_kraken2_ref = file("${ref_dir}/Kraken2/", checkIfExists : true)
    ch_whippet_ref = file("${ref_dir}/whippet_index.jls", checkIfExists : true)

    ch_fasta = file("${ref_dir}/IRFinder/genome.fa", checkIfExists : true)
    ch_gtf = file("${ref_dir}/IRFinder/transcripts.gtf", checkIfExists : true)

    ch_filter_gtf = file("${ref_dir}/filtered_genes.gtf", checkIfExists : true)
    ch_transcript_fasta = file("${ref_dir}/transcripts.fa", checkIfExists : true)

    ch_versions = Channel.empty()

    emit:
    fasta            = ch_fasta            //    path: genome.fasta
    gtf              = ch_gtf              //    path: genome.gtf
    transcript_fasta = ch_transcript_fasta //    path: transcript.fasta
    irfinder_ref     = ch_irfinder_ref     //    path: IRFinder
    star_index       = ch_star_ref         //    path: ctat_genome_lib_build_dir/ref_genome.fa.star.idx
    kraken2          = ch_kraken2_ref      //    path: Kraken2
    whippet_index    = ch_whippet_ref      //    path: whippet_index.jls
    ctat             = ch_ctat_ref         //    path: ctat_genome_lib_build_dir
    decontaminer_human_ribo = ch_decontaminer_human_ribo_db
    decontaminer_bacterial = ch_decontaminer_bacterial_db
    decontaminer_fungi = ch_decontaminer_fungi_db
    decontaminer_virus = ch_decontaminer_virus_db
    sortmerna_fastas = ch_sortmerna_fastas //    path: rRNA/*.fasta

    versions = ch_versions                     // channel: [ versions.yml ]
}
