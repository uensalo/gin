#!/bin/bash

# Executables
FMD_DIR=../bin
QUERY_SCRIPT=generate_queries.py

# Set the hard-coded directories
LOG_DIR=../log/pipeline
QUERY_DIR=../res/query/pipeline
INDEX_OUTPUT_DIR=../res/index/pipeline
PERMUTATION_DIR=../res/permutation/pipeline
COMMON_INDEX_DIR=../res/index/common
COMMON_PERMUTATION_DIR=../res/permutation/common

mkdir -p $LOG_DIR
mkdir -p $INDEX_OUTPUT_DIR
mkdir -p $QUERY_DIR
mkdir -p $PERMUTATION_DIR

# Query generation parameters
QUERY_LEN=10
NO_QUERIES=1048576
SEED=420
TEMPERATURE=1e2
COOLING=0.99
BATCH_SIZE=4096

# Get the input file name, time, depth and number of threads
INPUT_FILE=$1
TIME=$2
DEPTH=$3
QUERY_NUM_THREADS=$4
PERMUTATION_NUM_THREADS=$4

# Get the base name of the input file
BASENAME=$(basename $INPUT_FILE .fmdg)

LOG_FILE="$LOG_DIR/${BASENAME}_log.txt"
QUERY_FILE="$QUERY_DIR/${BASENAME}_query.fmdq"

# Generate the queries
python3 $QUERY_SCRIPT $INPUT_FILE $QUERY_LEN $NO_QUERIES "$QUERY_DIR/_" $SEED > $QUERY_FILE

# Check if a common permutation file already exists, otherwise create one
COMMON_PERMUTATION_FILE="$COMMON_PERMUTATION_DIR/${BASENAME}_permutation.fmdp"
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

# Benchmark the index with the query set
echo "Running query benchmark" | tee -a $LOG_FILE
$FMD_DIR/fmd query enumerate -r $INDEX_FILE -i $QUERY_FILE -j $QUERY_NUM_THREADS -b $BATCH_SIZE -v 2>> $LOG_FILE
