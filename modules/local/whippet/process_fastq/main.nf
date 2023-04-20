process WHIPPET_PROCESS_FASTQ {
    tag "$meta.id"
    label 'process_medium'

     container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://github.com/JPTeamIOR/nf-PCAProfiler/releases/download/v0.1/Whippet.sif' :
        'cloxd/whippet:1.6.1' }"

    input:
    tuple val(meta), path(reads)
    path whippet_ref

    output:
    path "*.psi.gz"    , emit: psi
    path "*.gene.tpm.gz", emit : gene_tpm
    path "*.isoform.tpm.gz", emit : isoform_tpm
    path "*.jnc.gz", emit : jnc
    path "*.map.gz", emit : map
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ? task.ext.args : ''
    def reads_arg = ""
    def prefix = task.ext.prefix ?: "${meta.id}"

    /// bug when the third line is not a simple + symbol.
    if (meta.single_end) {
        if ( "$reads".endsWith('.gz') ){
            reads_arg = "<( zcat ${reads} | grep . | awk '{if ( NR %4 == 3 ) { print \"+\" } else { print } }' )"
        } else {
            reads_arg = "<( cat ${reads} | grep . | awk '{if ( NR %4 == 3 ) { print \"+\" } else { print } }' )"
        }
    } else {
        if ( "${reads[0]}".endsWith('.gz') ){
            reads_arg = "<( zcat ${reads[0]} | grep . |  awk '{if ( NR %4 == 3 ) { print \"+\" } else { print } }' ) <( zcat ${reads[1]} | grep . | awk '{if ( NR %4 == 3 ) { print \"+\" } else { print } }' )"
        } else {
            reads_arg = "<( cat ${reads[0]} | grep . | awk '{if ( NR %4 == 3 ) { print \"+\" } else { print } }' ) <( cat ${reads[1]} | grep . | awk '{if ( NR %4 == 3 ) { print \"+\" } else { print } }' )"
        }
    }

    """


    whippet-quant.jl -x $whippet_ref -o ./$prefix $args $reads_arg

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Whippet: \$(whippet-index.jl -h 2>&1 | head -1 | awk '{print \$2}' )
        Julia : \$(julia --version | awk '{print \$3}' )
    END_VERSIONS
    """
}
