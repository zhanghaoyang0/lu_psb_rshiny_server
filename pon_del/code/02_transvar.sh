#!/bin/bash

set -x

task_id="$1"
input_type="$2"

# task_id="protein_del";input_type="protein"
# task_id="genomic_del";input_type="genomic"
# task_id="transcript_del";input_type="transcript"

BASE_DIR="/srv/shiny-server/pon_del"
TEMP_DIR="${BASE_DIR}/temp"

echo "Running transvar for task ID: $task_id"
echo "Input type: $input_type"

REF_FASTA="/ref/hg38.fa"

if [ "$input_type" = "protein" ]; then
  transvar panno -l "${TEMP_DIR}/${task_id}_transvarInput.txt" --ensembl --reference "${REF_FASTA}" > "${TEMP_DIR}/${task_id}_transvarOutput.txt"
fi

if [ "$input_type" = "genomic" ]; then
  transvar ganno -l "${TEMP_DIR}/${task_id}_transvarInput.txt" --ensembl --reference "${REF_FASTA}" > "${TEMP_DIR}/${task_id}_transvarOutput.txt"
fi

if [ "$input_type" = "transcript" ]; then
  transvar canno -l "${TEMP_DIR}/${task_id}_transvarInput.txt" --ensembl --reference "${REF_FASTA}" > "${TEMP_DIR}/${task_id}_transvarOutput.txt"
fi


# Check if transvar output was created
if [ -f "${TEMP_DIR}/${task_id}_transvarOutput.txt" ]; then
    echo "Transvar output created successfully"
    echo "Output file size: $(wc -l < ${TEMP_DIR}/${task_id}_transvarOutput.txt) lines"
else
    echo "Error: Transvar output file was not created"
    exit 1
fi

echo "Done!"