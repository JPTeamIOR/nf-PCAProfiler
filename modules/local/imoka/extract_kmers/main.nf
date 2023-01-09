process IMOKA_EXTRACT_KMERS {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://github.com/RitchieLabIGH/iMOKA/releases/download/v1.1/iMOKA' :
        'cloxd/imoka:1.1' }"

    input:
    tuple val(meta), path(reads)
    val k_len

    output:
    tuple val(meta), path("preprocess/*/*.tsv.sorted.bin"), emit: sorted_bin
    path  'versions.yml'    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def read_col = meta.single_end ? "\$(realpath $reads)" : "\$(realpath ${reads[0]});\$(realpath ${reads[1]})"
    def max_mem = task.memory ? " ${task.memory.toBytes().intdiv(1000000000)}" : '12'
    """
    echo -e "${meta.id}\tNA\t${read_col}" > ./${meta.id}.tsv
    preprocess.sh -i ./${meta.id}.tsv -k $k_len -t $task.cpus -r $max_mem -o ./preprocess/ $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iMOKA: v1.1
        KMC3: \$(echo \$(kmc --version | head -n 1 | awk '{print \$4 }' ) )
    END_VERSIONS
    """
}
