#!/bin/bash

# Executables
FMD_DIR=../bin
QUERY_SCRIPT=generate_queries.py

# Set the hard-coded directories
LOG_DIR=../log/pipeline
QUERY_DIR=../res/query/pipeline
INDEX_OUTPUT_DIR=../res/index/pipeline
PERMUTATION_DIR=../res/permutation/pipeline

mkdir -p $LOG_DIR
mkdir -p $INDEX_OUTPUT_DIR
mkdir -p $QUERY_DIR
mkdir -p $PERMUTATION_DIR

# Query generation parameters
QUERY_LEN=10
NO_QUERIES=10000
SEED=420
TEMPERATURE=1e2
COOLING=0.99

# Get the input file name, time, depth and number of threads
INPUT_FILE=$1
TIME=$2
DEPTH=$3
NUM_THREADS=$4

# Get the base name of the input file
BASENAME=$(basename $INPUT_FILE .fmdg)

# Set the output file names
PERMUTATION_FILE="$PERMUTATION_DIR/${BASENAME}_permutation.fmdp"
INDEX_FILE="$INDEX_OUTPUT_DIR/${BASENAME}_index.fmdi"
LOG_FILE="$LOG_DIR/${BASENAME}_log.txt"
QUERY_FILE="$QUERY_DIR/${BASENAME}_query.fmdq"

# Generate the queries
python3 $QUERY_SCRIPT $INPUT_FILE $QUERY_LEN $NO_QUERIES "$QUERY_DIR/_" $SEED > $QUERY_FILE

# Run the permutation operation
echo "Running permutation" | tee -a $LOG_FILE
$FMD_DIR/fmd permutation -i $INPUT_FILE -t $TIME -u $TIME -e $TEMPERATURE -c $COOLING -d $DEPTH -o $PERMUTATION_FILE -v -j $NUM_THREADS 2>> $LOG_FILE

# Run the index operation
echo "Running index" | tee -a $LOG_FILE
$FMD_DIR/fmd index -i $INPUT_FILE -p $PERMUTATION_FILE -o $INDEX_FILE -v 2>> $LOG_FILE

# Benchmark the index with the query set
echo "Running query benchmark" | tee -a $LOG_FILE
$FMD_DIR/fmd query enumerate -r $INDEX_FILE -i $QUERY_FILE -j $NUM_THREADS -v 2>> $LOG_FILE
