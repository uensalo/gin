#!/bin/bash

LOG_DIR=../log/query_speed_vs_cached_query_length
QUERY_DIR=../res/query/query_speed_vs_cached_query_length
INDEX_OUTPUT_DIR=../res/index/query_speed_vs_cached_query_length
PERMUTATION_DIR=../res/permutation/query_speed_vs_cached_query_length

rm -rf $LOG_DIR
rm -rf $QUERY_DIR
rm -rf $INDEX_OUTPUT_DIR
rm -rf $PERMUTATION_DIR