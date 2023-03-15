
//
// Check input reference and get channels
//

workflow CHECK_REFERENCE {
    take:
        ref_dir
    main:

    // Check few files to ensure that they are all present
    file("$projectDir/assets/ref_check.txt").each((v)=>{
        if (! file("${ref_dir}/${v}").exists ) {
            error "The file ${} in the reference does not exists."
        }
    })

    ch_sortmerna_fastas = files("${ref_dir}/rRNA/*.fasta")

    ch_irfinder_ref = file("${ref_dir}/IRFinder/", checkIfExists : true)
    ch_ctat_ref = file("${ref_dir}/ctat_genome_lib_build_dir/", checkIfExists : true)
    ch_decontaminer_ref= Channel.empty()
    if (params.conatminant_method == 'all' ||  params.conatminant_method == 'decontaminer' ){
        ch_decontaminer_ref = file("${ref_dir}/Decontaminer/", checkIfExists : true)
    }

    ch_star_ref = file("${ref_dir}/ctat_genome_lib_build_dir/ref_genome.fa.star.idx", checkIfExists : true)
    ch_kraken2_ref = file("${ref_dir}/Kraken2/", checkIfExists : true)
    ch_whippet_ref=file("${ref_dir}/whippet_index.jls", checkIfExists : true)

    ch_fasta=file("${ref_dir}/IRFinder/genome.fa", checkIfExists : true)
    ch_gtf=file("${ref_dir}/IRFinder/transcripts.gtf", checkIfExists : true)

    ch_filter_gtf = GTF_GENE_FILTER(ch_fasta, ch_gtf).gtf
    ch_transcript_fasta = GFFREAD_GTF_TO_FASTA(ch_fasta, ch_filter_gtf).fasta

    ch_versions         = ch_versions.mix(GTF_GENE_FILTER.out.versions)
    ch_versions         = ch_versions.mix(GFFREAD_GTF_TO_FASTA.out.versions)

    emit:
    fasta            = ch_fasta            //    path: genome.fasta
    gtf              = ch_gtf              //    path: genome.gtf
    transcript_fasta = ch_transcript_fasta //    path: transcript.fasta
    irfinder_ref     = ch_irfinder_ref     //    path: IRFinder
    star_index       = ch_star_ref         //    path: ctat_genome_lib_build_dir/ref_genome.fa.star.idx
    kraken2          = ch_kraken2_ref      //    path: Kraken2
    whippet_index    = ch_whippet_ref      //    path: whippet_index.jls
    ctat             = ch_ctat_ref         //    path: ctat_genome_lib_build_dir
    decontaminer     = ch_decontaminer_ref //    path: decontaminer
    sortmerna_fastas = ch_sortmerna_fastas //    path: rRNA/*.fasta

    versions = ch_versions                     // channel: [ versions.yml ]
}
