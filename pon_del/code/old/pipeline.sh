#!/bin/bash

# Check if task_id and input_type are provided
if [ $# -lt 2 ]; then
    echo "Error: Please provide a task ID and input type"
    echo "Usage: $0 <task_id> <input_type>"
    exit 1
fi

task_id="$1"
input_type="$2"
BASE_DIR="/srv/shiny-server/pon_del"
TEMP_DIR="${BASE_DIR}/temp"

echo "Running task id: $task_id"
echo "Input type: $input_type"

# make transvar input
echo "Adding pchange ..."
Rscript "${BASE_DIR}/code/01_transvarInput.r" "$task_id" "$input_type"
# head temp/${task_id}_transvarInput.txt

# run transvar
echo "Running transvar ..."
sh "${BASE_DIR}/code/02_transvar.sh" "$task_id" "$input_type"
# head temp/${task_id}_transvarOutput.txt

# add features
echo "Adding features ..."
Rscript "${BASE_DIR}/code/03_addfeat.r" "$task_id" "$input_type"
# head temp/${task_id}_withfeat.csv

# predict
echo "Running prediction ..."
python3 "${BASE_DIR}/code/04_predict.py" "$task_id"
