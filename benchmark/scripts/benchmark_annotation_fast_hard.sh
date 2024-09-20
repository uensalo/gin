#!/bin/bash
#########################################################
# Transcriptomics
#########################################################
# permutation benchmark
./gin_run_fast.sh ../input/gencode.v40.ging \
-L "gencode.v40.ging_permutation_fast" \
-c "12" \
-r "32" \
-m "-1" \
-M "-1" \
-j "1" \
-l "16 32 64 128 256 512 1024 2048 4096" \
-t "225 450 900 1800 3600" \
-p "12" \
-h

# principal benchmark
./gin_run_fast.sh ../input/gencode.v40.ging \
-L "gencode.v40.ging_principal_fast" \
-c "0 1 2 3 4 5 6 7 8 9 10 11 12" \
-r "32" \
-m "-1" \
-M "-1" \
-j "1" \
-l "16 32 64 128 256 512 1024 2048 4096" \
-t "3600" \
-p "12" \
-h

# decode benchmark
./gin_run_fast.sh ../input/gencode.v40.ging \
-L "gencode.v40.ging_rate_decode_fast" \
-c "12" \
-r "32" \
-m "-1" \
-M "-1" \
-j "1" \
-l "16 32 64 128 256 512 1024 2048 4096" \
-t "3600" \
-p "12" \
-d \
-h

# plain benchmark
./gin_run_fast.sh ../input/gencode.v40.ging \
-L "gencode.v40.ging_baseline_plain_fast" \
-c "0" \
-r "32" \
-m "-1" \
-M "-1" \
-j "1" \
-l "16 32 64 128 256 512 1024 2048 4096" \
-t "0" \
-p "0" \
-h

# permutation baseline benchmark
./gin_run_fast.sh ../input/gencode.v40.ging  \
-L "gencode.v40.ging_baseline_permutation_fast" \
-c "0" \
-r "32" \
-m "-1" \
-M "-1" \
-j "1" \
-l "16 32 64 128 256 512 1024 2048 4096" \
-t "3600" \
-p "12" \
-h

# cache baseline benchmark
./gin_run_fast.sh ../input/gencode.v40.ging \
-L "gencode.v40.ging_baseline_cache_fast" \
-c "12" \
-r "32" \
-m "-1" \
-M "-1" \
-j "1" \
-l "16 32 64 128 256 512 1024 2048 4096" \
-t "0" \
-p "0" \
-h

# cache permutation baseline benchmark
./gin_run_fast.sh ../input/gencode.v40.ging \
-L "gencode.v40.ging_baseline_cache_permutation_fast" \
-c "12" \
-r "32" \
-m "-1" \
-M "-1" \
-j "1" \
-l "16 32 64 128 256 512 1024 2048 4096" \
-t "3600" \
-p "12" \
-h
