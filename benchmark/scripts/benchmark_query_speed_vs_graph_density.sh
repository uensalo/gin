#!/bin/bash

# Executables
FMD_DIR=../bin
QUERY_SCRIPT=generate_queries.py
GRAPH_SCRIPT=generate_random_fmdg.py

# Set the hard-coded directories
LOG_DIR=../log/query_speed_vs_graph_density
QUERY_DIR=../res/query/query_speed_vs_graph_density
INDEX_OUTPUT_DIR=../res/index/query_speed_vs_graph_density
PERMUTATION_DIR=../res/permutation/query_speed_vs_graph_density
GRAPH_DIR=../res/graph/query_speed_vs_graph_density

mkdir -p $LOG_DIR
mkdir -p $INDEX_OUTPUT_DIR
mkdir -p $QUERY_DIR
mkdir -p $PERMUTATION_DIR
mkdir -p $GRAPH_DIR

NO_VERTICES=$1
VERTEX_LABEL_LEN=$2
QUERY_LEN=10
NO_QUERIES=65536
SEED=420
QUERY_NUM_THREADS=16
DEPTH=6
TEMPERATURE=1e2
COOLING=0.99
TIME=300
BATCH_SIZE=256
PERMUTATION_NUM_THREADS=128

# Get the graph density values
shift 2
GRAPH_DENSITIES=("$@")

# Loop over each graph density and perform the benchmark
for DENSITY in "${GRAPH_DENSITIES[@]}"
do
    echo "[fmd:benchmark] Running benchmark for graph density $DENSITY:"

    # Calculate the number of edges
    NO_EDGES=$(printf "%.0f" $(echo "$NO_VERTICES * $DENSITY" | bc))

    # Generate the graph
    GRAPH_FILE="$GRAPH_DIR/graph_density_${DENSITY}.fmdg"
    python3 $GRAPH_SCRIPT $NO_VERTICES $NO_EDGES $VERTEX_LABEL_LEN  $SEED > $GRAPH_FILE

    # Generate the queries
    QUERY_FILE="$QUERY_DIR/query_speed_vs_graph_density_density_${DENSITY}.fmdq"
    python3 $QUERY_SCRIPT $GRAPH_FILE $QUERY_LEN $NO_QUERIES "$QUERY_DIR/_" $SEED > $QUERY_FILE

    # Set the output file names based on the density
    PERMUTATION_FILE="$PERMUTATION_DIR/permutation_density_${DENSITY}.fmdp"
    INDEX_FILE="$INDEX_OUTPUT_DIR/index_density_${DENSITY}.fmdi"
    LOG_FILE="$LOG_DIR/log_density_${DENSITY}.txt"
    touch $LOG_FILE

    # Run the permutation operation
    $FMD_DIR/fmd permutation -i $GRAPH_FILE -d $DEPTH -t $TIME -u $TIME -e $TEMPERATURE -c $COOLING -o $PERMUTATION_FILE -v -j $PERMUTATION_NUM_THREADS 2>> $LOG_FILE

    # Run the index operation
    $FMD_DIR/fmd index -i $GRAPH_FILE -p $PERMUTATION_FILE -o $INDEX_FILE -v 2>> $LOG_FILE

    # Benchmark the index with the query set
    $FMD_DIR/fmd query enumerate -r $INDEX_FILE -i $QUERY_FILE -j $QUERY_NUM_THREADS -b $BATCH_SIZE -v 2>> $LOG_FILE

done