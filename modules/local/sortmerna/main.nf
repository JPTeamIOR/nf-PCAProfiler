process SORTMERNA {
    tag "$meta.id"
    label "process_high"

    conda "bioconda::sortmerna=2.1b"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sortmerna:2.1b--he860b03_4' :
        'quay.io/biocontainers/sortmerna:2.1b--he860b03_4' }"

    input:
    tuple val(meta), path(reads)
    path  fastas
    path  db

    output:
    tuple val(meta), path("*non_rRNA.fastq.gz"), emit: reads
    tuple val(meta), path("*.log")     , emit: log
    path  "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def references = ""
    for (ref in fastas) {
        references = references == "" ? "${ref},${ref}.db" : "${references}:${ref},${ref}.db"
    }
    if (meta.single_end) {
        """
        zcat $reads > ${prefix}.fastq
        sortmerna \\
            --ref $references \\
            --reads ${meta.id}.fastq \\
            -a $task.cpus \\
            --aligned rRNA_reads \\
            --fastx \\
            --other non_rRNA_reads \\
            --log \\
            $args

        mv non_rRNA_reads.fastq ${prefix}.non_rRNA.fastq
        gzip ${prefix}.non_rRNA.fastq
        mv rRNA_reads.log ${prefix}.sortmerna.log

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sortmerna: \$(echo \$(sortmerna --version 2>&1) | sed 's/^.*SortMeRNA version //; s/ Build Date.*\$//')
        END_VERSIONS
        """
    } else {
        """
        merge_fq.sh ${reads[0]} ${reads[1]} ${prefix}.fastq
        sortmerna \\
            --ref $references \\
            --reads ${prefix}.fastq \\
            -a $task.cpus \\
            --aligned rRNA_reads \\
            --fastx \\
            --other non_rRNA_reads \\
            --paired_in \\
            --log \\
            $args
        split_fq.sh non_rRNA_reads.fastq ${prefix}_1.non_rRNA.fastq ${prefix}_2.non_rRNA.fastq
        gzip ${prefix}_1.non_rRNA.fastq
        gzip ${prefix}_2.non_rRNA.fastq
        mv rRNA_reads.log ${prefix}.sortmerna.log

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sortmerna: \$(echo \$(sortmerna --version 2>&1) | sed 's/^.*SortMeRNA version //; s/ Build Date.*\$//')
        END_VERSIONS
        """
    }
}
