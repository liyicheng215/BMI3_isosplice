#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <input_file.fq>"
    exit 1
fi

input_fq=$1
input_basename=$(basename "$input_fq" .fq)
genome="ce11.fa"

minimap2 -ax splice --MD --secondary=no $genome $input_fq > ${input_basename}.sam


samtools view -bS ${input_basename}.sam > ${input_basename}.bam
samtools sort -o ${input_basename}_sorted.bam ${input_basename}.bam
samtools index ${input_basename}_sorted.bam

rm ${input_basename}.sam ${input_basename}.bam

echo "BAM and BAI files have been created: ${input_basename}_sorted.bam and ${input_basename}_sorted.bam.bai"
