#!/usr/bin/env bash
input=$1
genome_fasta=$2
output=$3
threads=$4

picard -Xmx110g -XX:ParallelGCThreads=${threads} CollectWgsMetrics \
       INPUT=$input \
       OUTPUT=$output \
       READ_LENGTH=150 \
       REFERENCE_SEQUENCE=$genome_fasta
