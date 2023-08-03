#!/bin/bash
./benchmark_query_speed_vs_cache_depth.sh ../res/graph/gencode.v40.fmdg 2 4 6 8
./benchmark_query_speed_vs_query_length.sh ../res/graph/gencode.v40.fmdg 2 4 6 8 10 75 100 200 300 400 500 600 700 800 900 1000
./benchmark_query_speed_vs_matches_returned.sh ../res/graph/gencode.v40.fmdg 1 2 4 8 16 32 64 128 256 512 1024 -1
./benchmark_query_speed_vs_num_threads.sh ../res/graph/gencode.v40.fmdg 1 2 4 8 16 32 64 128
./benchmark_query_speed_vs_permutation_time.sh ../res/graph/gencode.v40.fmdg 1 60 120 180 360 720 1440 3600 7200
./benchmark_query_speed_vs_permutation_depth.sh ../res/graph/gencode.v40.fmdg 2 4 6 8 10
./benchmark_query_speed_vs_permutation_depth_variable_time.sh ../res/graph/gencode.v40.fmdg 2 4 6 8 10
./benchmark_query_speed_vs_sample_rate.sh ../res/graph/gencode.v40.fmdg 16 32 64 128
./pipeline.sh ../res/graph/gencode.v40.fmdg 3600 4 32 128 20