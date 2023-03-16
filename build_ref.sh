#!/bin/bash
### Use this script to build a reference for the pipeline
set -e

BR_VERSION=1.0
BF_NAME="build_ref.sh"
SCRIPT_DIR=$(dirname -- "$0" )
THREADS=4
WHIPPET_TSL=1
BASE_DIR=$(realpath ./)

CTAT_LIB_SUPP="https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/MUTATION_LIB_SUPPLEMENT/GRCh38.mutation_lib_supplement.Jul272020.tar.gz"

STAR_FUSION_IMG="https://data.broadinstitute.org/Trinity/CTAT_SINGULARITY/STAR-Fusion/star-fusion.v1.12.0.simg"
CTAT_MUT_IMG="https://data.broadinstitute.org/Trinity/CTAT_SINGULARITY/CTAT_MUTATIONS/ctat_mutations.v3.3.1.simg"
K2_PLUSPF="https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_20221209.tar.gz"
IRFINDER_IMG="https://github.com/RitchieLabIGH/IRFinder/releases/download/v2.0.1/IRFinder"
WHIPPET_IMG="https://github.com/JPTeamIOR/nf-PCAProfiler/releases/download/v0.1/Whippet.sif"

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

[ -f ${IR_IMG} ] || wget -O ${IR_IMG} $IRFINDER_IMG
[ -f ${W_IMG} ] || wget -O ${W_IMG} $WHIPPET_IMG
[ -f ${CM_IMG} ] || wget -O ${CM_IMG} $CTAT_MUT_IMG
[ -f ${SF_IMG} ] || wget -O ${SF_IMG} $STAR_FUSION_IMG

cd $OUTDIR

singularity exec -e ${SF_IMG} \
   /usr/local/src/STAR-Fusion/ctat-genome-lib-builder/prep_genome_lib.pl \
      --genome_fa ${reference_genome} \
      --gtf ${reference_gtf} \
      --fusion_annot_lib ${annot_lib} \
      --annot_filter_rule ${filter_rules} \
      --pfam_db current \
      --dfam_db human \
      --CPU ${THREADS} \
      --human_gencode_filter

CTAT_REF=$(realpath ./ctat_genome_lib_build_dir)
STAR_REF=$(realpath ./STAR/ )
ln -s ${CTAT_REF}/ref_genome.fa.star.idx ${STAR_REF}


### Add CTAT Mutation Lib Supplement

cd ${CTAT_REF}
wget -O ./mutation_lib_supplement.tar.gz $CTAT_LIB_SUPP
tar xvf ./mutation_lib_supplement.tar.gz
rm -fr ./mutation_lib_supplement.tar.gz
cd ${OUTDIR}


singularity exec -e ${CM_IMG} /usr/local/src/ctat-mutations/mutation_lib_prep/ctat-mutation-lib-integration.py \
      --genome_lib_dir ${CTAT_REF}

### Add Open-CRAVAT
OPEN_CRAVAT="${CTAT_REF}/ctat_mutation_lib/cravat"
mkdir -p ${OPEN_CRAVAT}
singularity exec -e ${CM_IMG} oc config md ${OPEN_CRAVAT}
singularity exec -e ${CM_IMG} oc module install-base
singularity exec -e ${CM_IMG} oc module install --yes vest chasmplus vcfreporter mupit clinvar



singularity exec -e ${IR_IMG} \
   BuildRefFromSTARRef -l \
      -f ${reference_genome} \
      -g ${reference_gtf} \
      -t ${THREADS} \
      -r ./IRFinder/ \
      -e ${irf_spike} \
      -M ${irf_mapability} \
      -x ${STAR_REF}



### Whippet ref generation


singularity exec -e ${W_IMG} gawk -v lev="$WHIPPET_TSL" ' $3 == "exon" { if ( match($0, /transcript_support_level "([0-9]+)"/ , ary) ) { if ( ary[1] <= lev ) { print }  }   } ' $reference_gtf > ${reference_gtf}.filtered.gtf
singularity exec -e ${W_IMG} whippet-index.jl --fasta $reference_genome --gtf ${reference_gtf}.filtered.gtf -x ./whippet_index

rm -fr ${reference_gtf}.filtered.gtf

### Decontaminer ref is to be downloaded manually from gdrive ( probably we will exclude it)


### Kraken2 db
mkdir ${OUTDIR}/Kraken2/
wget -O ${OUTDIR}/Kraken2/k2_pluspf.tar.gz $K2_PLUSPF
cd ${OUTDIR}/Kraken2/
tar k2_pluspf.tar.gz
rm k2_pluspf.tar.gz

cd ${OUTDIR}
mkdir -p ${OUTDIR}/rRNA/

cd ${OUTDIR}/rRNA/
while read line ; do
    wget $line
done < ${SCRIPT_DIR}/assets/rrna-db-defaults.txt


cd $BASE_DIR

