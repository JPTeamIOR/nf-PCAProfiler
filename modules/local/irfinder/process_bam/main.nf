process IRFINDER_PROCESS_BAM {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://github.com/RitchieLabIGH/IRFinder/releases/download/v2.0.1/IRFinder' :
        'cloxd/irfinder:2.0.1' }"

    input:
    tuple val(meta), path(unsorted_bam)
    path irfinder_ref

    output:
    path "./results"    , emit: results
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ? task.ext.args : ''

    """
    IRFinder BAM  \\
        -r $irfinder_ref \\
        -d ./results/ \\
        $args $unsorted_bam

    gzip ./results/*.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        IRFinder: \$(IRFinder --version | sed -e "s/IRFinder version: //g" )
    END_VERSIONS
    """
}
