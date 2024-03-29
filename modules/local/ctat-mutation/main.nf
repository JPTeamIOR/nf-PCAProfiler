process CTAT_MUTATION_PROCESS_BAM {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://pcaprofilertest.tk/static/ctat_mutations.v3.3.1.simg' :
        'cloxd/ctat-mutations:3.3.1' }"
    containerOptions "${workflow.containerEngine == 'singularity' ? '-B ' : '--volume '}${ctat_ref}/ctat_mutation_lib/cravat/:/mnt/modules/"

    input:
    tuple val(meta), path(bam) , path(bai)
    path ctat_ref

    output:
    tuple val(meta), path("${prefix}.vcf.gz")    , emit: vcf
    tuple val(meta), path("*classifier.vcf.gz")    , emit: vcf_filtered
    tuple val(meta), path("*cancer.vcf")    , emit: vcf_cancer
    tuple val(meta), path("*cancer.tsv")    , emit: tsv_cancer
    path 'versions.yml' , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ? task.ext.args : ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    /usr/local/src/ctat-mutations/ctat_mutations \
            --genome_lib_dir $ctat_ref \
            --bam $bam \
            --sample_id $meta.id \
            --cpu $task.cpus \
            -O ./ $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CTAT-Mutations: \$( /usr/local/src/ctat-mutations/ctat_mutations --version 2>&1 )
    END_VERSIONS
    """
}
