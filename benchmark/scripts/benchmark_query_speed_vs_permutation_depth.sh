#!/bin/bash

# Executables
FMD_DIR=../bin
QUERY_SCRIPT=generate_queries.py

# Set the hard-coded directories
LOG_DIR=../log/query_speed_vs_permutation_depth
QUERY_DIR=../res/query/query_speed_vs_permutation_depth
INDEX_OUTPUT_DIR=../res/index/query_speed_vs_permutation_depth
PERMUTATION_DIR=../res/permutation/query_speed_vs_permutation_depth

mkdir -p $LOG_DIR
mkdir -p $INDEX_OUTPUT_DIR
mkdir -p $QUERY_DIR
mkdir -p $PERMUTATION_DIR

QUERY_LEN=25
NO_QUERIES=1000
SEED=420
NUM_THREADS=4
TEMPERATURE=1e2
COOLING=0.99
TIME=30

# Get the input file name and the depth values
INPUT_FILE=$1
shift
DEPTHS=("$@")

# Get the base name of the input file
BASENAME=$(basename $INPUT_FILE .fmdg)

# Generate the queries
QUERY_FILE="$QUERY_DIR/${BASENAME}_query_speed_vs_permutation_depth.fmdq"
python3 $QUERY_SCRIPT $INPUT_FILE $QUERY_LEN $NO_QUERIES "$QUERY_DIR/_" $SEED > $QUERY_FILE

# Loop over each depth value and perform the benchmark
for DEPTH in "${DEPTHS[@]}"
do
    echo "Running benchmark for permutation depth $DEPTH"

    # Set the output file names based on the depth
    PERMUTATION_FILE="$PERMUTATION_DIR/${BASENAME}_permutation_depth_${DEPTH}.fmdp"
    INDEX_FILE="$INDEX_OUTPUT_DIR/${BASENAME}_index_depth_${DEPTH}.fmdi"
    LOG_FILE="$LOG_DIR/${BASENAME}_log_depth_${DEPTH}.txt"
    touch $LOG_FILE

    # Run the permutation and index operations, redirecting stderr to the log file
    $FMD_DIR/fmd permutation -i $INPUT_FILE -t $TIME -u $TIME -e $TEMPERATURE -c $COOLING -d $DEPTH -o $PERMUTATION_FILE -v -j $NUM_THREADS 2>> $LOG_FILE
    $FMD_DIR/fmd index -i $INPUT_FILE -p $PERMUTATION_FILE -o $INDEX_FILE -v 2>> $LOG_FILE

    # Benchmark the index with the query set
    $FMD_DIR/fmd query enumerate -r $INDEX_FILE -i $QUERY_FILE -j $NUM_THREADS -v 2>> $LOG_FILE
done
