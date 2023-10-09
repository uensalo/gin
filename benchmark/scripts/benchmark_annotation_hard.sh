#!/bin/bash
#########################################################
# Transcriptomics
#########################################################
# permutation benchmark
./fmd_run.sh ../input/gencode.v40.fmdg \
-L "gencode.v40.fmdg_permutation" \
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
./fmd_run.sh ../input/gencode.v40.fmdg \
-L "gencode.v40.fmdg_principal" \
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
./fmd_run.sh ../input/gencode.v40.fmdg \
-L "gencode.v40.fmdg_rate_decode" \
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
./fmd_run.sh ../input/gencode.v40.fmdg \
-L "gencode.v40.fmdg_baseline_plain" \
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
./fmd_run.sh ../input/gencode.v40.fmdg  \
-L "gencode.v40.fmdg_baseline_permutation" \
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
./fmd_run.sh ../input/gencode.v40.fmdg \
-L "gencode.v40.fmdg_baseline_cache" \
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
./fmd_run.sh ../input/gencode.v40.fmdg \
-L "gencode.v40.fmdg_baseline_cache_permutation" \
-c "12" \
-r "32" \
-m "-1" \
-M "-1" \
-j "1" \
-l "16 32 64 128 256 512 1024 2048 4096" \
-t "3600" \
-p "12" \
-h
