#!/bin/bash
LOG_DIR=../log/query_speed_vs_no_bipartite_vertices
QUERY_DIR=../res/query/query_speed_vs_no_bipartite_vertices
INDEX_OUTPUT_DIR=../res/index/query_speed_vs_no_bipartite_vertices
GRAPH_DIR=../res/graph/query_speed_vs_no_bipartite_vertices

rm -rf $LOG_DIR
rm -rf $QUERY_DIR
rm -rf $INDEX_OUTPUT_DIR
rm -rf $PERMUTATION_DIR
rm -rf $GRAPH_DIR