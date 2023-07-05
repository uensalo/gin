#!/bin/bash
./benchmark_query_speed_vs_query_length.sh ../res/graph/GRCh38-20-0.10b.fmdg 8 10 75 100 200 300 400 500 600 700 800 900 1000
./benchmark_query_speed_vs_num_threads.sh ../res/graph/GRCh38-20-0.10b.fmdg 1 2 4 8 16 32 64 128
./benchmark_query_speed_vs_sample_rate.sh ../res/graph/GRCh38-20-0.10b.fmdg 1 4 16 32 64 1024 16384 65536
./benchmark_query_speed_vs_permutation_time.sh ../res/graph/GRCh38-20-0.10b.fmdg 1 60 120 180 360 720 1440 3600 7200
./benchmark_query_speed_vs_permutation_depth.sh ../res/graph/GRCh38-20-0.10b.fmdg 2 4 6 8
./benchmark_query_speed_vs_permutation_depth_variable_time.sh ../res/graph/GRCh38-20-0.10b.fmdg 2 4 6 8
./pipeline.sh ../res/graph/GRCh38-20-0.10b.fmdg 3600 4 128
