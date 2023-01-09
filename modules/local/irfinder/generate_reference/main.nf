process IRFINDER_GENERATE_REF {
    tag "$fasta"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://github.com/RitchieLabIGH/IRFinder/releases/download/v2.0.1/IRFinder' :
        'cloxd/irfinder:2.0.1' }"

    input:
    path star_ref
    path fasta
    path gtf
    path mapability
    val mapability_len
    path extraref
    path blacklist
    path roi

    output:
    path 'irfinder'    , emit: reference
    path 'irfinder/STAR', emit : star_ref
    path 'irfinder/Mapability', emit : mapability
    path 'irfinder/IRFinder', emit : irfinder_inner
    path 'irfinder/genome.fa', emit : genome
    path 'irfinder/transcripts.gtf', emit : transcripts
    path 'versions.yml' , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    def args = task.ext.args ? task.ext.args : ''
    def irfinder_process_name = 'BuildRef'
    def star_ref_arg = ''
    if (star_ref.exists() && star_ref.isDirectory()) {
        irfinder_process_name = 'BuildRefFromSTARRef'
        star_ref_arg = "-x $star_ref"
    }
    def mapability_arg = mapability.exists() ? "-M $mapability" : mapability_len ? "-m $mapability_len" : ''
    def extraref_arg = extraref.exists()  ? "-e $extraref" : ''
    def blacklist_arg = blacklist.exists()   ?  "-b $blacklist " : ''
    def roi_arg = roi.exists()  ?  "-R $roi" : ''

    """
    IRFinder $irfinder_process_name -l \\
        -f $fasta \\
        -g $gtf \\
        -t $task.cpus \\
        $mapability_arg $roi_arg $extraref_arg $blacklist_arg \\
        -r ./irfinder \\
        $args \\
        $star_ref_arg

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        IRFinder: \$(IRFinder --version )
    END_VERSIONS
    """
}
