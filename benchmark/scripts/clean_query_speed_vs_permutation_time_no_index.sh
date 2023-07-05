#!/bin/bash

LOG_DIR=../log/query_speed_vs_permutation_time_no_index
QUERY_DIR=../res/query/query_speed_vs_permutation_time_no_index
INDEX_OUTPUT_DIR=../res/index/query_speed_vs_permutation_time_no_index
PERMUTATION_DIR=../res/permutation/query_speed_vs_permutation_time_no_index

rm -rf $LOG_DIR
rm -rf $QUERY_DIR
rm -rf $INDEX_OUTPUT_DIR
rm -rf $PERMUTATION_DIR