include { GUNZIP as GUNZIP_FASTA  } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GTF  } from '../../modules/nf-core/gunzip/main'
include { STAR_GENOMEGENERATE } from '../../modules/nf-core/star/genomegenerate/main'
include { GFFREAD_GTF_TO_FASTA   } from '../../modules/local/gffread_gtf_to_fasta'
include { GTF_GENE_FILTER  } from '../../modules/local/gtf_gene_filter'
include { IRFINDER_GENERATE_REF } from '../../modules/local/irfinder/generate_reference'

workflow PREPARE_GENOME {
    take:

    main:

    ch_versions = Channel.empty()
    if (params.fasta.endsWith('.gz')) {
        ch_fasta    = GUNZIP_FASTA([ [:], params.fasta ]).gunzip.map { it[1] }
        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
    } else {
        ch_fasta = file(params.fasta)
    }
    if (params.gtf.endsWith('.gz')) {
        ch_gtf      = GUNZIP_GTF([ [:], params.gtf ]).gunzip.map { it[1] }
        ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)
    } else {
        ch_gtf = file(params.gtf)
    }

    //
    // Link to STAR index or generate from scratch if required
    //
    ch_star_index = Channel.empty()
    if (params.star_index) {
        ch_star_index = file(params.star_index)
    } else {
        ch_star_index = STAR_GENOMEGENERATE(ch_fasta, ch_gtf).index
        ch_versions   = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)
    }

    //
    // Link to IRFinder index or generate from scratch if required
    //
    ch_irfinder_ref = Channel.empty()
    if (params.irfinder_ref) {
        ch_irfinder_ref = file(params.irfinder_ref)
    } else {
        ch_irfinder_ref = IRFINDER_GENERATE_REF(
                ch_star_index,
                ch_fasta,
                ch_gtf,
                params.mapability ? Channel.fromPath(params.mapability, checkIfExists:true) : file('NoMapability'),
                params.mapability_len ? params.mapability_len : 75  ,
                params.extraref ?  Channel.fromPath(params.extraref, checkIfExists:true) : file('NoExtraRef'),
                params.blacklist? Channel.fromPath(params.blacklist, checkIfExists:true) : file('NoBlackList'),
                params.roi? Channel.fromPath(params.roi, checkIfExists:true) : file('NoROI')
                ).reference
        ch_versions   = ch_versions.mix(IRFINDER_GENERATE_REF.out.versions)
    }
    ch_filter_gtf = GTF_GENE_FILTER(ch_fasta, ch_gtf).gtf
    ch_transcript_fasta = GFFREAD_GTF_TO_FASTA(ch_fasta, ch_filter_gtf).fasta
    ch_versions         = ch_versions.mix(GTF_GENE_FILTER.out.versions)
    ch_versions         = ch_versions.mix(GFFREAD_GTF_TO_FASTA.out.versions)



    emit:
    fasta            = ch_fasta            //    path: genome.fasta
    gtf              = ch_gtf              //    path: genome.gtf
    star_index       = ch_star_index       //    path: star/index/
    transcript_fasta = ch_transcript_fasta //    path: transcript.fasta
    irfinder_ref     = ch_irfinder_ref     //    path: IRFinder

    versions = ch_versions                     // channel: [ versions.yml ]
}
