process DECONTAMINER {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://github.com/JPTeamIOR/nf-PCAProfiler/releases/download/v0.1/decontaminer.img' :
        'cloxd/decontaminer:1.4' }"

    input:
    tuple val(meta), path(reads)
    path human_ribo_db // 'NO_RIBO' to disable
    path bacterial_db // 'NO_BACTERIAL' to disable
    path fungi_db // 'NO_FUNGI' to disable
    path virus_db // 'NO_VIRUS' to disable
    val qualityFilter // false to disable

    output:

    path  'versions.yml'           , emit: versions
    path 'results/RESULTS/', emit : results

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def human_ribo = "$human_ribo_db" != 'NO_RIBO'
    def bacterial = "$bacterial_db" != 'NO_BACTERIAL'
    def fungi = "$fungi_db" != 'NO_FUNGI'
    def virus = "$virus_db" != 'NO_VIRUS'
    def ribo_arg = human_ribo ? 'y' : 'n'
    def qualityFilter_arg = qualityFilter ? 'y' : 'n'
    def org_arg = "-${bacterial? 'b' : ''}${fungi ? 'f' : ''}${virus ? 'v' : ''}"
    if (!(bacterial || fungi || virus)) {
        error 'One of the organism must be selected.'
    }
    def pairing_arg = 'S'
    def init_input_1 = ''
    def init_input_2 = ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        init_input_1 = "zcat $reads > ./input/${prefix}.fastq "
    } else {
        init_input_1 = "zcat ${reads[0]} | awk 'BEGIN {tot=1} {if ( NR % 4 == 1 ) { print \"@READ_\" tot \"/1\"; tot=tot+1 } else { print } }' > ./input/${prefix}_1.fastq"
        init_input_2 = "zcat ${reads[1]} | awk 'BEGIN {tot=1} {if ( NR % 4 == 1 ) { print \"@READ_\" tot \"/2\"; tot=tot+1 } else { print } }' > ./input/${prefix}_2.fastq"
        pairing_arg = 'P'
    }
    def bacterial_process = ''
    def fungi_process = ''
    def virus_process = ''
    if (bacterial) {
        bacterial_process = "filterBlastInfo.sh -i ./results/RESULTS/BACTERIA/ -s $pairing_arg "
        if (!bacterial_db.exists()) {
            error "$bacterial_db doesn't exists"
        }
    }
    if (fungi) {
        fungi_process = "filterBlastInfo.sh -i ./results/RESULTS/FUNGI/ -s $pairing_arg "
        if (!fungi_db.exists()) {
            error "$fungi_db doesn't exists"
        }
    }
    if (virus) {
        virus_process = "filterBlastInfo.sh -i ./results/RESULTS/VIRUSES/ -s $pairing_arg -V V "
        if (!virus_db.exists()) {
            error "$virus_db doesn't exists"
        }
    }
    """
    echo -e '[ext_soft]
    SAMTOOLS_EXEC=\$(which samtools)
    FASTX_EXEC=\$(which fastq_quality_filter)
    BLASTN_EXEC=\$(which blastn)
    SORTMERNA_EXEC=\$(which sortmerna)
    [cont_db]
    RIBO_DB=\$(realpath ${human_ribo ? human_ribo_db : './' } )
    RIBO_NAME=rRNA
    BACTERIA_DB=\$(realpath ${bacterial ? bacterial_db : './' } )
    BACTERIA_NAME=Bacteria_newblast
    FUNGI_DB=\$(realpath ${fungi ? fungi_db : './' } )
    FUNGI_NAME=Fungi_new
    VIRUSES_DB=\$(realpath ${virus ? virus_db : './' } )
    VIRUSES_NAME=Virus_new
    ' > ./config.txt
    mkdir -p ./input/
    $init_input_1
    $init_input_2
    decontaMiner.sh -i ./input/ -o ./results/ -c config.txt -F fastq -s $pairing_arg \
         -Q $qualityFilter_arg -R $ribo_arg $org_arg $args
    rm -fr ./input/
    $virus_process
    $bacterial_process
    $fungi_process
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            decontaminer: 1.4
            \$(echo \$(blastn -version | head -1 ))
    END_VERSIONS
    """
}
