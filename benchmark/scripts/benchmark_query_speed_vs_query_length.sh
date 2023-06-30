#!/bin/bash

# Executables
FMD_DIR=../bin
QUERY_SCRIPT=generate_queries.py

# Set the hard-coded directories
LOG_DIR=../log/query_speed_vs_query_length
QUERY_DIR=../res/query/query_speed_vs_query_length
INDEX_OUTPUT_DIR=../res/index/query_speed_vs_query_length
PERMUTATION_DIR=../res/permutation/query_speed_vs_query_length

mkdir -p $LOG_DIR
mkdir -p $INDEX_OUTPUT_DIR
mkdir -p $QUERY_DIR
mkdir -p $PERMUTATION_DIR

NO_QUERIES=1000
SEED=420
NUM_THREADS=4
TEMPERATURE=1e2
COOLING=0.99
TIME=30
DEPTH=2

# Get the input file name, query length start, end, and increment
INPUT_FILE=$1
QUERY_LENGTH_START=$2
QUERY_LENGTH_END=$3
QUERY_LENGTH_INC=$4

# Get the base name of the input file
BASENAME=$(basename $INPUT_FILE .fmdg)

# Set the output file names
PERMUTATION_FILE="$PERMUTATION_DIR/${BASENAME}_permutation.fmdp"
INDEX_FILE="$INDEX_OUTPUT_DIR/${BASENAME}_index.fmdi"

# Run the permutation and index operations
$FMD_DIR/fmd permutation -i $INPUT_FILE -t $TIME -u $TIME -e $TEMPERATURE -c $COOLING -d $DEPTH -o $PERMUTATION_FILE -j $NUM_THREADS
$FMD_DIR/fmd index -i $INPUT_FILE -p $PERMUTATION_FILE -o $INDEX_FILE

# Generate the queries with varying length, and perform the benchmark for each
for QUERY_LENGTH in $(seq $QUERY_LENGTH_START $QUERY_LENGTH_INC $QUERY_LENGTH_END)
do
    QUERY_LENGTH=$(printf "%.0f" $QUERY_LENGTH) # Round up to the nearest integer
    echo "Running benchmark for query length $QUERY_LENGTH"

    # Generate the queries
    QUERY_FILE="$QUERY_DIR/${BASENAME}_query_length_${QUERY_LENGTH}.fmdq"
    python3 $QUERY_SCRIPT $INPUT_FILE $QUERY_LENGTH $NO_QUERIES "$QUERY_DIR/_" $SEED > $QUERY_FILE

    # Set the log file name based on the query length
    LOG_FILE="$LOG_DIR/${BASENAME}_log_query_length_${QUERY_LENGTH}.txt"
    touch $LOG_FILE

    # Benchmark the index with the query set, redirecting stderr to the log file
    $FMD_DIR/fmd query enumerate -r $INDEX_FILE -i $QUERY_FILE -j $NUM_THREADS -v 2>> $LOG_FILE
done
