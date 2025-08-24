#!/bin/bash

# Loop over all BAM files
for bam in *.bam; do
    start=$(date +%s)  # start time

    echo "Processing $bam ..."
    featureCounts -s 0 -a gencode.v48.primary_assembly.annotation.gtf \
        -o quants/${bam%.bam}_featurecounts.txt \
        "$bam"

    end=$(date +%s)  # end time
    runtime=$(( (end - start) / 60 ))  # in minutes

    echo "✅ Completed $bam in $runtime minutes."
    echo "------------------------------------"
done
