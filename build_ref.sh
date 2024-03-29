#!/bin/bash
### Use this script to build a reference for the pipeline
set -e

BR_VERSION=1.0
BF_NAME="build_ref.sh"
SCRIPT_DIR=$(realpath $(dirname -- "$0" ))
THREADS=4
WHIPPET_TSL=1
BASE_DIR=$(realpath ./)

CTAT_LIB_SUPP="https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/MUTATION_LIB_SUPPLEMENT/GRCh38.mutation_lib_supplement.Jul272020.tar.gz"
DFAM_DB_URL="https://www.dfam.org/releases/Dfam_3.7/infrastructure/dfamscan/homo_sapiens_dfam.hmm"


STAR_FUSION_IMG="https://data.broadinstitute.org/Trinity/CTAT_SINGULARITY/STAR-Fusion/star-fusion.v1.12.0.simg"
CTAT_MUT_IMG="https://pcaprofilertest.tk/static/ctat_mutations.v3.3.1.simg"
K2_PLUSPF="https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_20221209.tar.gz"
IRFINDER_IMG="https://github.com/RitchieLabIGH/IRFinder/releases/download/v2.0.1/IRFinder"
WHIPPET_IMG="https://github.com/JPTeamIOR/nf-PCAProfiler/releases/download/v0.1/Whippet.sif"
PYTHON_IMG='https://depot.galaxyproject.org/singularity/python:3.9--1'
GFFREAD_IMG='https://depot.galaxyproject.org/singularity/gffread:0.12.1--h8b12597_0'
SORTMERNA_IMG="https://depot.galaxyproject.org/singularity/sortmerna:2.1b--he860b03_4"

help(){
    echo "Usage: ${BF_NAME} [-g|--gencode-version][-o|--outdir][-h|--help][-v|--version]
    The current script aims to build the references for the pipeline nf-PCaProfiler for any version of GENCODE annotation and hg38 reference genome."
    exit 2
}

SHORT="g:,o:,t:,h,v"
LONG="gencode-version:,outdir:,threads:,help,version"

OPTS=$(getopt -a -n ${BF_NAME} --options $SHORT --longoptions $LONG -- "$@" )

if [ $# -eq 0 ]; then
    help
fi

eval set -- "$OPTS"


GENCODE_VERSION="42"
OUTDIR=$(realpath ./REF/ )

while :
do
case "$1" in
    -g | --gencode-version )
        GENCODE_VERSION="$2"
        if ! [[ "${GENCODE_VERSION}" =~ ^[0-9]+$ ]]; then
            echo "Error! --gencode-version must be a positive integer number"
            exit 1
        fi
        shift 2
        ;;
    -t | --threads )
        THREADS="$2"
        if ! [[ "${THREADS}" =~ ^[0-9]+$ ]]; then
            echo "Error! --threads must be a positive integer number"
            exit 1
        fi
        shift 2
        ;;
    -o | --outdir )
        OUTDIR=$(realpath $2)
        shift 2
        ;;
    -- )
        shift;
        break
        ;;
    *)
        echo "Unexpected option: $1"
        help
        ;;
    esac
done

filter_rules=$(realpath "${SCRIPT_DIR}/assets/AnnotFilterRule.pm")
annot_lib=$(realpath "${SCRIPT_DIR}/assets/fusion_lib.dat.gz")
irf_mapability=$(realpath "${SCRIPT_DIR}/assets/MapabilityExclusion.100bp.bed.gz")
irf_spike=$(realpath "${SCRIPT_DIR}/assets/RNA.SpikeIn.ERCC.fasta.gz")

for fname in ${filter_rules} ${annot_lib} ${irf_mapability} ${irf_spike} ; do
    if ! [ -f ${fname} ]; then
        echo "Error! missing " $(basename ${fname}) " in the asset of the pipeline"
        help
    fi
done


reference_gtf="${OUTDIR}/gencode_v${GENCODE_VERSION}/gencode.v${GENCODE_VERSION}.annotation.gtf"
reference_genome="${OUTDIR}/gencode_v${GENCODE_VERSION}/GRCh38.primary_assembly.genome.fa"


mkdir -p ${OUTDIR}/gencode_v${GENCODE_VERSION}/

if [ ! -f ${reference_gtf} ]; then
    wget -O ${reference_gtf}.gz https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE_VERSION}/gencode.v${GENCODE_VERSION}.annotation.gtf.gz
    gzip -d ${reference_gtf}.gz
fi

if [ ! -f ${reference_genome} ]; then
    wget -O ${reference_genome}.gz https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE_VERSION}/GRCh38.primary_assembly.genome.fa.gz
    gzip -d ${reference_genome}.gz
fi


mkdir -p ${OUTDIR}/singularity/

SF_IMG=${OUTDIR}/singularity/star-fusion.v1.12.0.simg
CM_IMG=${OUTDIR}/singularity/ctat_mut.img
IR_IMG=${OUTDIR}/singularity/irfinder.img
W_IMG=${OUTDIR}/singularity/whippet.img
PTS_IMG=${OUTDIR}/singularity/pathoscope2.img
PY_IMG=${OUTDIR}/singularity/python.img
GFR_IMG=${OUTDIR}/singularity/gffread.img
SMR_IMG=${OUTDIR}/singularity/sortmerna.img




cd $OUTDIR
LOG_DIR="${OUTDIR}/logs/"
mkdir -p ${LOG_DIR}

DFAM_DB=${OUTDIR}/homo_sapiens_dfam.hmm
if ! [ -f ${LOG_DIR}/.hmmpress.done ]; then
    ### Need to download the dfamdb due to issues with url ( in package is without www and it doens't work anymore)
    [ -f $DFAM_DB ] || wget -O $DFAM_DB $DFAM_DB_URL
    [ -f ${SF_IMG} ] || wget -O ${SF_IMG} $STAR_FUSION_IMG
    echo "STEP 1 - hmmpress"
    singularity exec -e ${SF_IMG} hmmpress $DFAM_DB > ${LOG_DIR}/hmmpress.stdout 2> ${LOG_DIR}/hmmpress.stderr
    touch ${LOG_DIR}/.hmmpress.done
else
    echo "STEP 1 - hmmpress skipped"
fi


if ! [ -f ${LOG_DIR}/.prep_genome_lib.done ]; then
    echo "STEP 2 - CTAT genome lib preparation"
    [ -f ${SF_IMG} ] || wget -O ${SF_IMG} $STAR_FUSION_IMG
    singularity exec -e ${SF_IMG} \
        /usr/local/src/STAR-Fusion/ctat-genome-lib-builder/prep_genome_lib.pl \
            --genome_fa ${reference_genome} \
            --gtf ${reference_gtf} \
            --fusion_annot_lib ${annot_lib} \
            --annot_filter_rule ${filter_rules} \
            --pfam_db current \
            --dfam_db $DFAM_DB \
            --CPU ${THREADS} \
            --human_gencode_filter > /dev/null 2> ${LOG_DIR}/STAR-Fusion.stderr

    CTAT_REF=$(realpath ./ctat_genome_lib_build_dir)
    STAR_REF=$(realpath ./STAR/ )
    ln -s ${CTAT_REF}/ref_genome.fa.star.idx ${STAR_REF}
    touch ${LOG_DIR}/.prep_genome_lib.done
else
    echo "STEP 2 - CTAT genome lib preparation skipped"
fi

reference_gtf=$(realpath ${OUTDIR}/gencode_v${GENCODE_VERSION}/gencode.v${GENCODE_VERSION}.annotation.gtf.revised.custom.gtf )
reference_genome=$(realpath ${OUTDIR}/gencode_v${GENCODE_VERSION}/*.genome.fa.pseudo_masked.fa )

CTAT_REF=$(realpath ./ctat_genome_lib_build_dir)
STAR_REF=$(realpath ./STAR/ )

### Add CTAT Mutation Lib Supplement
if ! [ -f ${LOG_DIR}/.mutation_lib_supplement.done ]; then
    echo "STEP 3 - CTAT genome lib supplement preparation"
    cd ${CTAT_REF}
    wget -O ./mutation_lib_supplement.tar.gz $CTAT_LIB_SUPP
    tar xvf ./mutation_lib_supplement.tar.gz
    rm -fr ./mutation_lib_supplement.tar.gz
    cd ${OUTDIR}
    touch ${LOG_DIR}/.mutation_lib_supplement.done
else
    echo "STEP 3 - CTAT genome lib supplement preparation skipped"
fi

if ! [ -f ${LOG_DIR}/.ctat_lib_integration.done ]; then
    echo "STEP 4 - CTAT genome lib supplement integration"
    [ -f ${CM_IMG} ] || wget -O ${CM_IMG} $CTAT_MUT_IMG
    singularity exec -e ${CM_IMG} /usr/local/src/ctat-mutations/mutation_lib_prep/ctat-mutation-lib-integration.py \
        --genome_lib_dir ${CTAT_REF} > ${LOG_DIR}/ctat-mutation-lib-integration.stdout 2> ${LOG_DIR}/ctat-mutation-lib-integration.stderr
    touch ${LOG_DIR}/.ctat_lib_integration.done
else
    echo "STEP 4 - CTAT genome lib supplement integration skipped"
fi

### Add Open-CRAVAT - created directly on the image
OPEN_CRAVAT="${CTAT_REF}/ctat_mutation_lib/cravat"
if ! [ -f ${LOG_DIR}/.open_cravat.done ]; then
    echo "STEP 5 - CTAT cravat preparation"
   mkdir -p ${OPEN_CRAVAT}
   [ -f ${CM_IMG} ] || wget -O ${CM_IMG} $CTAT_MUT_IMG
   singularity exec -B ${OPEN_CRAVAT}:/mnt/modules/ -e ${CM_IMG} oc module install-base > ${LOG_DIR}/cravat-base.stdout 2> ${LOG_DIR}/cravat-base.stderr
   singularity exec -B ${OPEN_CRAVAT}:/mnt/modules/ -e ${CM_IMG} oc module install --yes vest chasmplus vcfreporter mupit clinvar > ${LOG_DIR}/cravat-module.stdout 2> ${LOG_DIR}/cravat-module.stderr
   touch ${LOG_DIR}/.open_cravat.done
else
    echo "STEP 5 - CTAT cravat preparation skipped"
fi

### Cleanup of CTAT-ref

if ! [ -f ${LOG_DIR}/.ctat_cleanup.done ]; then
    echo "STEP 6 - CTAT cleanup"
   rm -fr ${OUTDIR}/ref_annot* ${OUTDIR}/pipeliner*  ${OUTDIR}/homo_sapiens_dfam* ${OUTDIR}/Pfam-A.*
   rm -fr ${OUTDIR}/_* ${OUTDIR}/PFAM* ${OUTDIR}/para.mask.list
   touch ${LOG_DIR}/.ctat_cleanup.done
else
    echo "STEP 6 - CTAT cleanup skipped"
fi

if ! [ -f ${LOG_DIR}/.irfinder.done ]; then
    echo "STEP 7 - IRFinder reference creation"
    [ -f ${IR_IMG} ] || wget -O ${IR_IMG} $IRFINDER_IMG
    singularity exec -e ${IR_IMG} IRFinder \
        BuildRefFromSTARRef -l \
          -f ${reference_genome} \
          -g ${reference_gtf} \
          -t ${THREADS} \
          -r ./IRFinder/ \
          -e ${irf_spike} \
          -M ${irf_mapability} \
          -x ${STAR_REF} > ${LOG_DIR}/IRFinder.stdout 2> ${LOG_DIR}/IRFinder.stderr
    touch ${LOG_DIR}/.irfinder.done
else
    echo "STEP 7 - IRFinder reference creation skipped"
fi



### Whippet ref generation

if ! [ -f ${LOG_DIR}/.whippet.done ]; then
    echo "STEP 8 - Whippet reference creation"
    [ -f ${W_IMG} ] || wget -O ${W_IMG} $WHIPPET_IMG
    singularity exec -e ${W_IMG} gawk -v lev="$WHIPPET_TSL" ' $3 == "exon" { if ( match($0, /transcript_support_level "([0-9]+)"/ , ary) ) { if ( ary[1] <= lev ) { print }  }   } ' $reference_gtf > ${reference_gtf}.filtered.gtf
    singularity exec -e ${W_IMG} whippet-index.jl --fasta $reference_genome --gtf ${reference_gtf}.filtered.gtf -x ./whippet_index > ${LOG_DIR}/Whippet.stdout 2> ${LOG_DIR}/Whippet.stderr
    rm -fr ${reference_gtf}.filtered.gtf
    touch ${LOG_DIR}/.whippet.done
else
    echo "STEP 8 - Whippet reference creation skipped"
fi
### Decontaminer ref is to be downloaded manually from gdrive ( probably we will exclude it)


### Kraken2 db
if ! [ -f ${LOG_DIR}/.kraken2.done ]; then
    echo "STEP 9 - Kraken2 reference download"
    mkdir ${OUTDIR}/Kraken2/
    wget -O ${OUTDIR}/Kraken2/k2_pluspf.tar.gz $K2_PLUSPF
    cd ${OUTDIR}/Kraken2/
    tar xvf k2_pluspf.tar.gz
    rm k2_pluspf.tar.gz
    touch ${LOG_DIR}/.kraken2.done
else
    echo "STEP 9 - Kraken2 reference download skipped"
fi

cd ${OUTDIR}

if ! [ -f ${LOG_DIR}/.rRNA.done ]; then
    echo "STEP 10 - rRNA reference download"
    [ -f ${SMR_IMG} ] || wget -O ${SMR_IMG} $SORTMERNA_IMG
    mkdir -p ${OUTDIR}/rRNA/
    cd ${OUTDIR}/rRNA/
    while read line ; do
        wget $line
        fname=$(echo $line | awk -F '/' '{print $NF}')
        singularity exec -e ${SMR_IMG} indexdb_rna --ref ${fname},${fname}.db >> ${LOG_DIR}/sortmerna.stdout 2>> ${LOG_DIR}/sortmerna.stderr
    done < ${SCRIPT_DIR}/assets/rrna-db-defaults.txt
    touch ${LOG_DIR}/.rRNA.done
else
    echo "STEP 10 - rRNA reference download skipped"
fi
## produce the filtered gtf and fasta genes

if ! [ -f ${LOG_DIR}/.pipeline_files.done ]; then
    echo "STEP 11 - Transcriptome preparation"
    [ -f ${PY_IMG} ] || wget -O ${PY_IMG} $PYTHON_IMG
    [ -f ${GFR_IMG} ] || wget -O ${GFR_IMG} $GFFREAD_IMG
    cd ${OUTDIR}
    cp ${SCRIPT_DIR}/bin/filter_gtf_for_genes_in_genome.py ./
    singularity exec -e ${PY_IMG} ./filter_gtf_for_genes_in_genome.py --gtf $reference_gtf --fasta $reference_genome -o ./filtered_genes.gtf > ${LOG_DIR}/filter_gtf_for_genes_in_genome.stdout 2> ${LOG_DIR}/filter_gtf_for_genes_in_genome.stderr
    rm filter_gtf_for_genes_in_genome.py
    singularity exec -e ${GFR_IMG} gffread -w ./transcripts.fa -g $reference_genome ./filtered_genes.gtf > ${LOG_DIR}/gffread.stdout 2> ${LOG_DIR}/gffread.stderr
    touch ${LOG_DIR}/.pipeline_files.done
else
    echo "STEP 11 - Transcriptome preparation skipped"
fi


## cleanup
echo "All Done."
rm -fr ${OUTDIR}/singularity

cd $BASE_DIR


