process CTAT_MUTATION_PROCESS_BAM {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://data.broadinstitute.org/Trinity/CTAT_SINGULARITY/CTAT_MUTATIONS/ctat_mutations.v3.3.1.simg' :
        'trinityctat/ctat_mutations:3.3.1' }"

    input:
    tuple val(meta), path(bam) , path(bai)
    path ctat_ref

    output:
    tuple val(meta), path("${prefix}/${prefix}.vcf.gz")    , emit: vcf
    tuple val(meta), path("${prefix}/*classifier.vcf.gz")    , emit: vcf_filtered
    tuple val(meta), path("${prefix}/*cancer.vcf")    , emit: vcf_cancer
    tuple val(meta), path("${prefix}/*cancer.tsv")    , emit: tsv_cancer
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
            -O $prefix $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CTAT-Mutations: \$( /usr/local/src/ctat-mutations/ctat_mutations --version 2>&1 )
    END_VERSIONS
    """
}
