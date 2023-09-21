#!/bin/bash
#########################################################
# Transcriptomics
#########################################################
# query length grouped by cache depth
# ./fmd_run.sh ../input/gencode.v40.fmdg \
# -L "gencode.v40.fmdg_c_l_hard" \
# -c "0 2 4 6 8 10 12" \
# -r "16" \
# -m "-1" \
# -M "-1" \
# -j "1" \
# -l "10 30 100 300 1000" \
# -t "3600" \
# -p "4" \
# -h

# query length grouped by matches returned
# ./fmd_run.sh ../input/gencode.v40.fmdg \
# -L "gencode.v40.fmdg_m_l_hard" \
# -c "10" \
# -r "16" \
# -m "-1 1 4 16 64 512 4096" \
# -M "-1" \
# -j "1" \
# -l "10 30 100 300 1000" \
# -t "3600" \
# -p "4" \
# -h

# thread benchmark
# ./fmd_run.sh ../input/gencode.v40.fmdg \
# -L "gencode.v40.fmdg_j_l_hard" \
# -c "10" \
# -r "16" \
# -m "-1" \
# -M "-1" \
# -j "1 2 4 8 16 32" \
# -l "10 30 100 300 1000" \
# -t "3600" \
# -p "4" \
# -h

# decoding speed benchmark
# ./fmd_run.sh ../input/gencode.v40.fmdg \
# -L "gencode.v40.fmdg_decode_hard" \
# -c "10" \
# -r "16 32 64" \
# -m "-1" \
# -M "65536" \
# -j "1" \
# -l "10 30 100 300 500" \
# -t "3600" \
# -p "4" \
# -d \
# -h

# extensive cache benchmark
#./fmd_run.sh ../input/gencode.v40.fmdg \
# -L "gencode.v40.fmdg_c_l2_hard" \
# -c "0 1 2 3 4 5 6 7 8 9 10 11 12" \
# -r "16" \
# -m "-1" \
# -M "-1" \
# -j "1" \
# -l "2 4 8 16 32 64 128 256 512" \
# -t "3600" \
# -p "4" \
# -h

# long cache benchmark
#./fmd_run.sh ../input/gencode.v40.fmdg \
#-L "gencode.v40.fmdg_c_l2_hard" \
#-c "0 1 2 3 4 5 6 7 8 9 10 11 12" \
#-r "16" \
#-m "-1" \
#-M "-1" \
#-j "1" \
#-l "1024 2048 4096" \
#-t "3600" \
#-p "4" \
#-h

# decoding speed benchmark
#./fmd_run.sh ../input/gencode.v40.fmdg \
#-L "gencode.v40.fmdg_decode2_hard" \
#-c "12" \
#-r "16 32 64" \
#-m "-1" \
#-M "65536" \
#-j "1" \
#-l "16 32 64 128 256 512 1024 2048 4096" \
#-t "3600" \
#-p "4" \
#-d \
#-h

# threading benchmark
#./fmd_run.sh ../input/gencode.v40.fmdg \
#-L "gencode.v40.fmdg_j_l2_hard" \
#-c "12" \
#-r "16" \
#-m "-1" \
#-M "-1" \
#-j "1 2 4 8 16 32" \
#-l "16 32 64 128 256 512 1024 2048 4096" \
#-t "3600" \
#-p "4" \
#-h

# permutation benchmark
#./fmd_run.sh ../input/gencode.v40.fmdg \
#-L "gencode.v40.fmdg_j_p_hard" \
#-c "12" \
#-r "64" \
#-m "-1" \
#-M "-1" \
#-j "1" \
#-l "16 32 64 128 256 512 1024 2048 4096" \
#-t "225 450 900 1800 3600" \
#-p "4" \
#-h

# long cache benchmark 2
#./fmd_run.sh ../input/gencode.v40.fmdg \
#-L "gencode.v40.fmdg_principal" \
#-c "0 1 2 3 4 5 6 7 8 9 10 11 12" \
#-r "16" \
#-m "-1" \
#-M "-1" \
#-j "1" \
#-l "16 32 64 128 256 512 1024 2048 4096" \
#-t "3600" \
#-p "4" \
#-h

#./fmd_run.sh ../input/gencode.v40.fmdg \
#-L "gencode.v40.fmdg_rate_decode" \
#-c "12" \
#-r "16 32 64" \
#-m "-1" \
#-M "-1" \
#-j "1" \
#-l "16 32 64 128 256 512 1024 2048 4096" \
#-t "3600" \
#-p "4" \
#-d \
#-h

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
-p "4" \
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
-p "4" \
-h
