process WHIPPET_GENERATE_REF {
    tag "$fasta"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://github.com/JPTeamIOR/nf-PCAProfiler/releases/download/v0.1/Whippet.sif' :
        'cloxd/whippet:1.6.1' }"

    input:
    path fasta
    path gtf
    val TSL // Transcript Support Level ( http://www.ensembl.org/Help/Glossary?id=492 )

    output:
    path 'whippet_index.jls' , emit: index
    path 'whippet_index.exons.tab.gz' , emit: exons
    path 'versions.yml' , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    def args = task.ext.args ? task.ext.args : ''


    """
    gawk -v lev="$TSL" ' \$3 == "exon" { if ( match(\$0, /transcript_support_level "([0-9]+)"/ , ary) ) { if ( ary[1] <= lev ) { print }  }   } ' $gtf > ${gtf}.filtered.gtf
    whippet-index.jl --fasta $fasta --gtf ${gtf}.filtered.gtf -x ./whippet_index $args
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Whippet: \$(whippet-index.jl -h 2>&1 | head -1 | awk '{print \$2}' )
        Julia : \$(julia --version | awk '{print \$3}' )
    END_VERSIONS
    """
}
