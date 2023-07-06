#!/bin/bash
rm -rf ../res/index/common
rm -rf ../res/permutation/common
./clean_query_speed_vs_graph_density.sh
./clean_query_speed_vs_num_threads.sh
./clean_query_speed_vs_permutation_depth.sh
./clean_query_speed_vs_permutation_depth_variable_time.sh
./clean_query_speed_vs_permutation_time_no_index.sh
./clean_query_speed_vs_query_length.sh
./clean_query_speed_vs_sample_rate.sh