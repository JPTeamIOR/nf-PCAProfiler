process STARFUSION_PROCESS {
    tag "$star_bam"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://data.broadinstitute.org/Trinity/CTAT_SINGULARITY/STAR-Fusion/star-fusion.v1.12.0.simg' :
        'trinityctat/starfusion:1.12.0' }"

    input:
    path chimeric_junction
    path ctat_ref

    output:
    path "./results"    , emit: results
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ? task.ext.args : ''

    """
     STAR-Fusion --genome_lib_dir $ctat_ref \
             -J $chimeric_junction \
             --output_dir results \
             --CPU $task.cpus $args


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        STAR-Fusion: \$(STAR-Fusion --version | awk '{ print \$NF }' )
    END_VERSIONS
    """
}
