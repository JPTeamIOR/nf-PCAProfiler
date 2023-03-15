#!/bin/bash

export NXF_SINGULARITY_CACHEDIR=$(realpath ./test_data/singularity/)

if [[ $1 == "dc" ]]; then
    nf-test test tests/modules/local/decontaminer/main.nf.test --profile singularity
    exit
fi

if [[ $1 == "wg" ]]; then
    nf-test test tests/modules/local/whippet/generate_reference/main.nf.test --profile singularity
    exit
fi

if [[ $1 == "wp" ]]; then
    nf-test test tests/modules/local/whippet/process_fastq/main.nf.test --profile singularity
    exit
fi

if [[ $1 == "imoka" ]]; then
    nf-test test tests/modules/local/imoka/extract_kmer/main.nf.test --profile singularity
    exit
fi


if [[ $1 == "align" ]]; then
    nf-test test tests/modules/nf-core/star/align/main.nf.test --profile singularity
    exit
fi


if [[ $1 == "sf" ]]; then
    nf-test test tests/modules/local/star-fusion/main.nf.test --profile singularity
    exit
fi

if [[ $1 == "cm" ]]; then
    nf-test test tests/modules/local/ctat-mutation/main.nf.test --profile singularity
    exit
fi

echo "test.sh [dc|wg|wp|imoka|align|sf|cm]"

