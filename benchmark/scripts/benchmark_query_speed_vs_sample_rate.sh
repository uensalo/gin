#!/bin/bash

# Executables
FMD_DIR=../bin
QUERY_SCRIPT=generate_queries.py

# Set the hard-coded directories
LOG_DIR=../log/query_speed_vs_sample_rate
QUERY_DIR=../res/query/query_speed_vs_sample_rate
INDEX_OUTPUT_DIR=../res/index/query_speed_vs_sample_rate

mkdir $LOG_DIR
mkdir $INDEX_OUTPUT_DIR
mkdir $QUERY_DIR

QUERY_LEN=8
NO_QUERIES=65536
SEED=420
QUERY_NUM_THREADS=4
BATCH_SIZE=64
INDEX_CONSTRUCTION_PARALLELISM=4

# Get the input file name and the sample rates
INPUT_FILE=$1
shift
SAMPLE_RATES=("$@")

# Get the base name of the input file
BASENAME=$(basename $INPUT_FILE .fmdg)

# Generate the queries
QUERY_FILE="$QUERY_DIR/${BASENAME}_query_speed_vs_sample_rate.fmdq"
python3 $QUERY_SCRIPT $INPUT_FILE $QUERY_LEN $NO_QUERIES "$QUERY_DIR/_" $SEED > $QUERY_FILE

# Keep track of all background PIDs
PIDS=()

# Loop over each sample rate and perform the benchmark
for SAMPLE_RATE in "${SAMPLE_RATES[@]}"
do
    echo "[fmd:benchmark] Running benchmark for sample rate $SAMPLE_RATE:"

    # Set the output file names based on the sample rate
    INDEX_FILE="$INDEX_OUTPUT_DIR/${BASENAME}_index_${SAMPLE_RATE}.fmdi"
    LOG_FILE="$LOG_DIR/${BASENAME}_log_sample_rate_${SAMPLE_RATE}.txt"
    touch $LOG_FILE

    # Run the index operation in background, redirecting stderr to the log file
    $FMD_DIR/fmd index -i $INPUT_FILE -s $SAMPLE_RATE -r $SAMPLE_RATE -o $INDEX_FILE -v 2>> $LOG_FILE &

    # Store the background job PID
    PIDS+=("$!")

    # If we have reached our limit for concurrent jobs, wait for them to finish
    if (( ${#PIDS[@]} % INDEX_CONSTRUCTION_PARALLELISM == 0 )); then
        wait "${PIDS[@]}"
        PIDS=()
    fi
done

# Wait for any remaining background jobs to finish
wait "${PIDS[@]}"

# Now perform the queries one by one
for SAMPLE_RATE in "${SAMPLE_RATES[@]}"
do
    echo "[fmd:benchmark] Running queries for sample rate $SAMPLE_RATE:"

    # Set the output file names based on the sample rate
    INDEX_FILE="$INDEX_OUTPUT_DIR/${BASENAME}_index_${SAMPLE_RATE}.fmdi"
    LOG_FILE="$LOG_DIR/${BASENAME}_log_sample_rate_${SAMPLE_RATE}.txt"

    # Run the query operation, redirecting stderr to the log file
    $FMD_DIR/fmd query find -r $INDEX_FILE -i $QUERY_FILE -j $QUERY_NUM_THREADS -b $BATCH_SIZE -v 2>> $LOG_FILE
done
