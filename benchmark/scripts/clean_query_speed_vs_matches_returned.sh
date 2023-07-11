#!/bin/bash

LOG_DIR=../log/query_speed_vs_matches_returned
QUERY_DIR=../res/query/query_speed_vs_matches_returned
INDEX_OUTPUT_DIR=../res/index/query_speed_vs_matches_returned
PERMUTATION_DIR=../res/permutation/query_speed_vs_matches_returned

rm -rf $LOG_DIR
rm -rf $QUERY_DIR
rm -rf $INDEX_OUTPUT_DIR
rm -rf $PERMUTATION_DIR