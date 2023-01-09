#!/bin/bash

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

echo "test.sh [dc|wg|wp]"

