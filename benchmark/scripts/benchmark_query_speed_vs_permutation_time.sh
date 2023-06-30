#!/bin/bash

# Executables
FMD_DIR=../bin
QUERY_SCRIPT=generate_queries.py

# Set the hard-coded directories
LOG_DIR=../log/query_speed_vs_permutation_time
QUERY_DIR=../res/query/query_speed_vs_permutation_time
INDEX_OUTPUT_DIR=../res/index/query_speed_vs_permutation_time
PERMUTATION_DIR=../res/permutation/query_speed_vs_permutation_time

mkdir -p $LOG_DIR
mkdir -p $INDEX_OUTPUT_DIR
mkdir -p $QUERY_DIR
mkdir -p $PERMUTATION_DIR

QUERY_LEN=25
NO_QUERIES=1000
SEED=420
NUM_THREADS=4
DEPTH=6
TEMPERATURE=1e2
COOLING=0.99

# Get the input file name and the time parameters
INPUT_FILE=$1
shift
TIME_PARAMETERS=("$@")

# Get the base name of the input file
BASENAME=$(basename $INPUT_FILE .fmdg)

# Generate the queries
QUERY_FILE="$QUERY_DIR/${BASENAME}_query_speed_vs_permutation_time.fmdq"
python3 $QUERY_SCRIPT $INPUT_FILE $QUERY_LEN $NO_QUERIES "$QUERY_DIR/_" $SEED > $QUERY_FILE

# Loop over each time parameter and perform the benchmark
for TIME in "${TIME_PARAMETERS[@]}"
do
    echo "Running benchmark for time parameter $TIME"

    # Set the output file names based on the time parameter
    PERMUTATION_FILE="$PERMUTATION_DIR/${BASENAME}_permutation_time_${TIME}.fmdp"
    INDEX_FILE="$INDEX_OUTPUT_DIR/${BASENAME}_index_time_${TIME}.fmdi"
    LOG_FILE="$LOG_DIR/${BASENAME}_log_time_${TIME}.txt"
    touch $LOG_FILE

    # Run the permutation operation
    $FMD_DIR/fmd permutation -i $INPUT_FILE -t $TIME -d $DEPTH -u $TIME -e $TEMPERATURE -c $COOLING -o $PERMUTATION_FILE -j $NUM_THREADS -v 2>> $LOG_FILE

    # Run the index operation
    $FMD_DIR/fmd index -i $INPUT_FILE -p $PERMUTATION_FILE -o $INDEX_FILE -v 2>> $LOG_FILE

    # Benchmark the index with the query set
    $FMD_DIR/fmd query enumerate -r $INDEX_FILE -i $QUERY_FILE -j $NUM_THREADS -v 2>> $LOG_FILE
done
