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

# Use full path for docker command
DOCKER_CMD="/usr/bin/docker"

if [ "$input_type" = "protein" ]; then
  ${DOCKER_CMD} exec transvar_forPonDel \
    bash -c "transvar panno -l /temp/${task_id}_transvarInput.txt --ensembl --reference /ref/hg38.fa > /temp/${task_id}_transvarOutput.txt"
fi

if [ "$input_type" = "genomic" ]; then
  ${DOCKER_CMD} exec transvar_forPonDel \
    bash -c "transvar ganno -l /temp/${task_id}_transvarInput.txt --ensembl --reference /ref/hg38.fa > /temp/${task_id}_transvarOutput.txt"
fi

if [ "$input_type" = "transcript" ]; then
  ${DOCKER_CMD} exec transvar_forPonDel \
    bash -c "transvar canno -l /temp/${task_id}_transvarInput.txt --ensembl --reference /ref/hg38.fa > /temp/${task_id}_transvarOutput.txt"
fi


# run transvar and redirect output inside the container


# Check if transvar output was created
if [ -f "${TEMP_DIR}/${task_id}_transvarOutput.txt" ]; then
    echo "Transvar output created successfully"
    echo "Output file size: $(wc -l < ${TEMP_DIR}/${task_id}_transvarOutput.txt) lines"
else
    echo "Error: Transvar output file was not created"
    exit 1
fi

echo "Done!"