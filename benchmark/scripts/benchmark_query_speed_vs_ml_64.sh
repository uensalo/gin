#!/bin/bash

# Executables
FMD_DIR=../bin
QUERY_SCRIPT=generate_queries.py

# Set the hard-coded directories
LOG_DIR=../log/query_speed_vs_ml_64
QUERY_DIR=../res/query/query_speed_vs_ml
INDEX_OUTPUT_DIR=../res/index/query_speed_vs_ml
PERMUTATION_DIR=../res/permutation/query_speed_vs_ml
COMMON_INDEX_DIR=../res/index/common
COMMON_PERMUTATION_DIR=../res/permutation/common

mkdir -p $LOG_DIR
mkdir -p $INDEX_OUTPUT_DIR
mkdir -p $QUERY_DIR
mkdir -p $PERMUTATION_DIR
mkdir -p COMMON_PERMUTATION_DIR

NO_QUERIES=10000
SEED=420
QUERY_NUM_THREADS=1
TEMPERATURE=1e2
COOLING=0.99
TIME=3600
DEPTH=6
BATCH_SIZE=8
PERMUTATION_NUM_THREADS=64
INDEX_SAMPLING_RATE=64

# Get the input file name, max-forks and lengths
INPUT_FILE=$1
shift
while getopts ":m:l:" opt; do
  case $opt in
    m) IFS=' ' read -r -a MAX_FORKS <<< "$OPTARG";;
    l) IFS=' ' read -r -a LENGTHS <<< "$OPTARG";;
    \?) echo "Invalid option -$OPTARG" >&2;;
  esac
done

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
COMMON_INDEX_FILE="$COMMON_INDEX_DIR/${BASENAME}_index_time_${TIME}_depth_${DEPTH}_sampling_rate_${INDEX_SAMPLING_RATE}.fmdi"
if [[ ! -f $COMMON_INDEX_FILE ]]; then
    mkdir -p $COMMON_INDEX_DIR
    # No common index file exists, need to create one
    # Run the index operation and save the index file to the common directory
    $FMD_DIR/fmd index -i $INPUT_FILE -p $PERMUTATION_FILE -o $COMMON_INDEX_FILE -s $INDEX_SAMPLING_RATE -r $INDEX_SAMPLING_RATE
fi

# Use the common index file for all the query benchmarks
INDEX_FILE=$COMMON_INDEX_FILE

# Pre-generate the query files in parallel (4 at a time)
for LENGTH in "${LENGTHS[@]}"
do
    ((i=i%4)); ((i++==0)) && wait

    QUERY_FILE="$QUERY_DIR/${BASENAME}_query_length_${LENGTH}.fmdq"
    if [[ ! -f $QUERY_FILE ]]; then
        python3 $QUERY_SCRIPT $INPUT_FILE $LENGTH $NO_QUERIES "$QUERY_DIR/_" $SEED > $QUERY_FILE &
    fi
done
wait

# Iterate over all combinations of max-forks and lengths
for FORK in "${MAX_FORKS[@]}"
do
  for LENGTH in "${LENGTHS[@]}"
  do
    echo "[fmd:benchmark] Running benchmark for max forks $FORK and query length $LENGTH:"

    # Set the log file name based on the max forks and length
    LOG_FILE="$LOG_DIR/${BASENAME}_log_max_forks_${FORK}_query_length_${LENGTH}.txt"
    touch $LOG_FILE

    # Benchmark the index with the query set, redirecting stderr to the log file
    QUERY_FILE="$QUERY_DIR/${BASENAME}_query_length_${LENGTH}.fmdq"
    $FMD_DIR/fmd query find -r $INDEX_FILE -i $QUERY_FILE -j $QUERY_NUM_THREADS -b $BATCH_SIZE -m $FORK -v 2>> $LOG_FILE
  done
done
