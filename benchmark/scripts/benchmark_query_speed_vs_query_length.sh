#!/bin/bash

# Executables
FMD_DIR=../bin
QUERY_SCRIPT=generate_queries.py

# Set the hard-coded directories
LOG_DIR=../log/query_speed_vs_query_length
QUERY_DIR=../res/query/query_speed_vs_query_length
INDEX_OUTPUT_DIR=../res/index/query_speed_vs_query_length
PERMUTATION_DIR=../res/permutation/query_speed_vs_query_length
COMMON_INDEX_DIR=../res/index/common
COMMON_PERMUTATION_DIR=../res/permutation/common

mkdir -p $LOG_DIR
mkdir -p $INDEX_OUTPUT_DIR
mkdir -p $QUERY_DIR
mkdir -p $PERMUTATION_DIR

NO_QUERIES=65536
SEED=420
QUERY_NUM_THREADS=16
TEMPERATURE=1e2
COOLING=0.99
TIME=3600
DEPTH=6
BATCH_SIZE=256
PERMUTATION_NUM_THREADS=128

# Get the input file name and the query lengths
INPUT_FILE=$1
shift
QUERY_LENGTHS=("$@")

# Get the base name of the input file
BASENAME=$(basename $INPUT_FILE .fmdg)

# Check if a common permutation file already exists, otherwise create one
COMMON_PERMUTATION_FILE="$COMMON_PERMUTATION_DIR/${BASENAME}_permutation_time_${TIME}_depth_${DEPTH}_threads_${PERMUTATION_NUM_THREADS}.fmdp"
if [[ ! -f $COMMON_PERMUTATION_FILE ]]; then
    mkdir -p $COMMON_PERMUTATION_DIR
    # No common permutation file exists, need to create one
    # Run the permutation operation
    $FMD_DIR/fmd permutation -i $INPUT_FILE -t $TIME -u $TIME -e $TEMPERATURE -c $COOLING -d $DEPTH -o $COMMON_PERMUTATION_FILE -j $PERMUTATION_NUM_THREADS
fi

# Use the common permutation file for all the indexing benchmarks
PERMUTATION_FILE=$COMMON_PERMUTATION_FILE

# Check if a common index file already exists, otherwise create one
COMMON_INDEX_FILE="$COMMON_INDEX_DIR/${BASENAME}.${TIME}.fmdi"
if [[ ! -f $COMMON_INDEX_FILE ]]; then
    mkdir -p $COMMON_INDEX_DIR
    # No common index file exists, need to create one
    # Run the index operation and save the index file to the common directory
    $FMD_DIR/fmd index -i $INPUT_FILE -p $PERMUTATION_FILE -o $COMMON_INDEX_FILE
fi

# Use the common index file for all the query benchmarks
INDEX_FILE=$COMMON_INDEX_FILE

# Generate the queries with specified lengths, and perform the benchmark for each
for QUERY_LENGTH in "${QUERY_LENGTHS[@]}"
do
    echo "[fmd:benchmark] Running benchmark for query length $QUERY_LENGTH:"

    # Generate the queries
    QUERY_FILE="$QUERY_DIR/${BASENAME}_query_length_${QUERY_LENGTH}.fmdq"
    python3 $QUERY_SCRIPT $INPUT_FILE $QUERY_LENGTH $NO_QUERIES "$QUERY_DIR/_" $SEED > $QUERY_FILE

    # Set the log file name based on the query length
    LOG_FILE="$LOG_DIR/${BASENAME}_log_query_length_${QUERY_LENGTH}.txt"
    touch $LOG_FILE

    # Benchmark the index with the query set, redirecting stderr to the log file
    $FMD_DIR/fmd query breadth -r $INDEX_FILE -i $QUERY_FILE -j $QUERY_NUM_THREADS -b $BATCH_SIZE -v 2>> $LOG_FILE

    # Log permutation parameters
    echo "[fmd:benchmark] Permutation Parameters:" >> $LOG_FILE
    echo "[fmd:benchmark] Permutation Time: $TIME" >> $LOG_FILE
    echo "[fmd:benchmark] Permutation Depth: $DEPTH" >> $LOG_FILE
    echo "[fmd:benchmark] Permutation Num Threads: $PERMUTATION_NUM_THREADS" >> $LOG_FILE
done