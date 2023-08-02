#!/bin/bash

# Executables
FMD_DIR=../bin
QUERY_SCRIPT=generate_queries.py
GRAPH_SCRIPT=generate_graph.py

# Set the hard-coded directories
LOG_DIR=../log/query_speed_vs_no_connected_vertices
QUERY_DIR=../res/query/query_speed_vs_no_connected_vertices
INDEX_OUTPUT_DIR=../res/index/query_speed_vs_no_connected_vertices
GRAPH_DIR=../res/graph/query_speed_vs_no_connected_vertices

mkdir -p $LOG_DIR
mkdir -p $INDEX_OUTPUT_DIR
mkdir -p $QUERY_DIR
mkdir -p $GRAPH_DIR

NO_QUERIES=1000
QUERY_NUM_THREADS=4
BATCH_SIZE=64

# Get the number of vertices and the query lengths
NUM_VERTICES=$1
shift
QUERY_LENGTHS=("$@")

# Generate the graph if it doesn't exist
GRAPH_FILE="$GRAPH_DIR/connected_graph_${NUM_VERTICES}.fmdg"
if [[ ! -f $GRAPH_FILE ]]; then
    python3 $GRAPH_SCRIPT "connected" $NUM_VERTICES $GRAPH_FILE
fi

# Check if an index file already exists in the INDEX_OUTPUT_DIR, otherwise create one
INDEX_FILE="$INDEX_OUTPUT_DIR/connected_graph_${NUM_VERTICES}.fmdi"
if [[ ! -f $INDEX_FILE ]]; then
    # No index file exists, need to create one
    # Run the index operation and save the index file to the INDEX_OUTPUT_DIR
    $FMD_DIR/fmd index -i $GRAPH_FILE -o $INDEX_FILE
fi

# Generate the queries with specified lengths, and perform the benchmark for each
for QUERY_LENGTH in "${QUERY_LENGTHS[@]}"
do
    # Set the log file name based on the query length and number of vertices
    LOG_FILE="$LOG_DIR/log_no_vertices_${NUM_VERTICES}_query_length_${QUERY_LENGTH}.txt"
    touch $LOG_FILE
    echo "[fmd:benchmark] Fully connected graph experiment parameters:" >> $LOG_FILE
    echo "[fmd:benchmark] Number of vertices: ${NUM_VERTICES}" >> $LOG_FILE
    echo "[fmd:benchmark] Query length: ${QUERY_LENGTH}" >> $LOG_FILE

    # Generate the queries
    QUERY_FILE="$QUERY_DIR/query_length_${QUERY_LENGTH}.fmdq"
    python3 $QUERY_SCRIPT $GRAPH_FILE $QUERY_LENGTH $NO_QUERIES "$QUERY_DIR/_" > $QUERY_FILE

    # Benchmark the index with the query set, redirecting stderr to the log file
    $FMD_DIR/fmd query find -r $INDEX_FILE -i $QUERY_FILE -j $QUERY_NUM_THREADS -b $BATCH_SIZE -v 2>> $LOG_FILE

done
