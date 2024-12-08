#!/bin/bash

# Define input files and parameters
FASTA_FILE="ce11.chr1.fa"
GTF_FILE="ce11.ncbiRefSeq.chr1.gtf"
AS_GTF_FILE="ce11.ncbiRefSeq.chr1.as.gtf"
GENE_DEPTH_FILE_PREFIX="ce11_gene_depth"
ISOFORM_DEPTH_FILE_PREFIX="ce11_isoform_depth"
TRANS_FA_FILE="ce11_trans_as.fa"
COVERAGES=(10 20 30)
THREADS=40
OUTPUT_DIR="simulated_data"

# Create the directory to store results
mkdir -p $OUTPUT_DIR

# 1. Create AS events
echo "Creating AS events..."
python -m yasim generate_as_events \
    -f $FASTA_FILE \
    -g $GTF_FILE \
    -o $AS_GTF_FILE \
    -c 5

# Iterate over each coverage
for COVERAGE in "${COVERAGES[@]}"; do
    echo "Processing coverage=$COVERAGE..."

    # 2. Create gene coverage file
    GENE_DEPTH_FILE="${GENE_DEPTH_FILE_PREFIX}_${COVERAGE}.tsv"
    echo "Generating gene depth file: $GENE_DEPTH_FILE"
    python -m yasim generate_gene_depth \
        -g $AS_GTF_FILE \
        -o $GENE_DEPTH_FILE \
        -d $COVERAGE

    # 3. Create isoform coverage file
    ISOFORM_DEPTH_FILE="${ISOFORM_DEPTH_FILE_PREFIX}_${COVERAGE}.tsv"
    echo "Generating isoform depth file: $ISOFORM_DEPTH_FILE"
    python -m yasim generate_isoform_depth \
        -g $AS_GTF_FILE \
        -d $GENE_DEPTH_FILE \
        -o $ISOFORM_DEPTH_FILE

    # 4. Transcribe GTF to FASTA
    echo "Transcribing GTF to FASTA..."
    python -m labw_utils.bioutils transcribe \
        -f $FASTA_FILE \
        -g $AS_GTF_FILE \
        -o $TRANS_FA_FILE

    # 5. Call LLRG to generate simulated data - ONT
    ONT_OUTPUT="ce11_ONT_${COVERAGE}"
    echo "Simulating ONT reads: $ONT_OUTPUT"
    python -m yasim pbsim3 \
        -e pbsim \
        -F ${TRANS_FA_FILE}.d \
        --hmm_method errhmm \
        --hmm_model ONT \
        -j $THREADS \
        -d $ISOFORM_DEPTH_FILE \
        -o $ONT_OUTPUT
    # Move the generated ONT fq files to the target directory
    mv ${ONT_OUTPUT}/*.fq $OUTPUT_DIR

    # 6. Call LLRG to generate simulated data - PacBio
    PACBIO_OUTPUT="ce11_Pacbio_${COVERAGE}"
    echo "Simulating PacBio reads: $PACBIO_OUTPUT"
    python -m yasim pbsim3 \
        -e pbsim \
        -F ${TRANS_FA_FILE}.d \
        --hmm_method errhmm \
        --hmm_model RSII \
        -j $THREADS \
        -d $ISOFORM_DEPTH_FILE \
        -o $PACBIO_OUTPUT

    # Move the generated PacBio fq files to the target directory
    mv ${PACBIO_OUTPUT}/*.fq $OUTPUT_DIR
done

echo "All tasks completed! All .fq files are stored in the '$OUTPUT_DIR' directory."
