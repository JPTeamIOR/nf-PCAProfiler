process STARFUSION_PROCESS {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://data.broadinstitute.org/Trinity/CTAT_SINGULARITY/STAR-Fusion/star-fusion.v1.12.0.simg' :
        'trinityctat/starfusion:1.12.0' }"

    input:
    tuple val(meta), path(chimeric_junction)
    path ctat_ref

    output:
    tuple val(meta), path("${prefix}/*.tsv") , emit: results
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ? task.ext.args : ''
    prefix   = task.ext.prefix ?: "${meta.id}"

    """

    STAR-Fusion --genome_lib_dir $ctat_ref \
             -J $chimeric_junction \
             --output_dir $prefix \
             --CPU $task.cpus $args


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        STAR-Fusion: \$(STAR-Fusion --version | awk '!/^\$/{ print \$NF }' )
    END_VERSIONS
    """
}
