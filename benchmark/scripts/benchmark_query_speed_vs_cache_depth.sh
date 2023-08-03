#!/bin/bash

# Executables
FMD_DIR=../bin
QUERY_SCRIPT=generate_queries.py

# Set the hard-coded directories
LOG_DIR=../log/query_length_vs_cache_depth
QUERY_DIR=../res/query/query_length_vs_cache_depth
COMMON_INDEX_DIR=../res/index/common
COMMON_PERMUTATION_DIR=../res/permutation/common

mkdir -p $LOG_DIR
mkdir -p $QUERY_DIR

NO_QUERIES=65536
SEED=420
QUERY_LENGTH=64
QUERY_NUM_THREADS=4
BATCH_SIZE=64
TEMPERATURE=1e2
COOLING=0.99
TIME=3600
DEPTH=6
PERMUTATION_NUM_THREADS=128

# Get the input file name and the cache depths
INPUT_FILE=$1
shift
CACHE_DEPTHS=("$@")

# Get the base name of the input file
BASENAME=$(basename $INPUT_FILE .fmdg)

# Check if a common index file already exists, otherwise create one
COMMON_INDEX_FILE="$COMMON_INDEX_DIR/${BASENAME}.${TIME}.fmdi"
if [[ ! -f $COMMON_INDEX_FILE ]]; then
    mkdir -p $COMMON_INDEX_DIR
    # No common index file exists, need to create one
    # Run the index operation and save the index file to the common directory
    $FMD_DIR/fmd index -i $INPUT_FILE -o $COMMON_INDEX_FILE
fi

# Use the common index file for all the query benchmarks
INDEX_FILE=$COMMON_INDEX_FILE

# Generate the queries
QUERY_FILE="$QUERY_DIR/${BASENAME}_query_length_${QUERY_LENGTH}.fmdq"
python3 $QUERY_SCRIPT $INPUT_FILE $QUERY_LENGTH $NO_QUERIES "$QUERY_DIR/_" $SEED > $QUERY_FILE

# Perform the benchmark for each cache depth
for CACHE_DEPTH in "${CACHE_DEPTHS[@]}"
do
    echo "[fmd:benchmark] Running benchmark for cache depth $CACHE_DEPTH:"

    # Set the log file name based on the cache depth
    LOG_FILE="$LOG_DIR/${BASENAME}_log_cache_depth_${CACHE_DEPTH}.txt"
    touch $LOG_FILE

    # Benchmark the index with the query set, redirecting stderr to the log file
    $FMD_DIR/fmd query find -r $INDEX_FILE -i $QUERY_FILE -j $QUERY_NUM_THREADS -b $BATCH_SIZE -c $CACHE_DEPTH -v 2>> $LOG_FILE

    # Log cache depth
    echo "[fmd:benchmark] Cache Depth: $CACHE_DEPTH" >> $LOG_FILE
done
