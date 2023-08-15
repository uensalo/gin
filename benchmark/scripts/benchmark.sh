#!/bin/bash
#########################################################
# Transcriptomics
#########################################################
# query length grouped by cache depth
./fmd_run.sh ../input/gencode.v40.fmdg \
-L "gencode.v40.fmdg_c_l" \
-c "0 2 4 6 8 10 12" \
-r "16" \
-m "-1" \
-M "-1" \
-j "1" \
-l "10 30 100 300 1000" \
-t "3600" \
-p "4"

# query length grouped by matches returned
./fmd_run.sh ../input/gencode.v40.fmdg \
-L "gencode.v40.fmdg_m_l" \
-c "10" \
-r "16" \
-m "-1 1 4 16 64 512 4096" \
-M "-1" \
-j "1" \
-l "10 30 100 300 1000" \
-t "3600" \
-p "4"

# thread benchmark
./fmd_run.sh ../input/gencode.v40.fmdg \
-L "gencode.v40.fmdg_j_l" \
-c "10" \
-r "16" \
-m "-1" \
-M "-1" \
-j "1 2 4 8 16 32" \
-l "10 30 100 300 1000" \
-t "3600" \
-p "4"

# decoding speed benchmark
./fmd_run.sh ../input/gencode.v40.fmdg \
-L "gencode.v40.fmdg_decode" \
-c "10" \
-r "16 32 64" \
-m "-1" \
-M "65536" \
-j "1" \
-l "10 30 100 300 1000" \
-t "3600" \
-p "4" \
-d

#########################################################
# Pangenomics
#########################################################
# query length grouped by cache depth
./fmd_run.sh ../input/GRCh38-20-0.10b.fmdg \
-L "GRCh38-20-0.10b.fmdg_c_l" \
-c "0 2 4 6 8 10 12" \
-r "16" \
-m "-1" \
-M "-1" \
-j "1" \
-l "10 30 100 300 1000" \
-t "3600" \
-p "4"

# query length grouped by matches returned
./fmd_run.sh ../input/GRCh38-20-0.10b.fmdg \
-L "GRCh38-20-0.10b.fmdg_m_l" \
-c "10" \
-r "16" \
-m "-1 1 4 16 64 512 4096" \
-M "-1" \
-j "1" \
-l "10 30 100 300 1000" \
-t "3600" \
-p "4"

# thread benchmark
./fmd_run.sh ../input/GRCh38-20-0.10b.fmdg \
-L "GRCh38-20-0.10b.fmdg_j_l" \
-c "10" \
-r "16" \
-m "-1" \
-M "-1" \
-j "1 2 4 8 16 32" \
-l "10 30 100 300 1000" \
-t "3600" \
-p "4"

# decoding speed benchmark
./fmd_run.sh ../input/GRCh38-20-0.10b.fmdg \
-L "GRCh38-20-0.10b.fmdg_decode" \
-c "10" \
-r "16 32 64" \
-m "-1" \
-M "65536" \
-j "1" \
-l "10 30 100 300 1000" \
-t "3600" \
-p "4" \
-d

#########################################################
# Duplicate Benchmarks
#########################################################

# partial match benchmark
./fmd_run.sh ../input/gencode.v40.fmdg \
-L "gencode.v40.fmdg_c_l2" \
-c "0 1 2 3 4 5 6 7 8 9 10 11 12" \
-r "16" \
-m "-1" \
-M "-1" \
-j "1" \
-l "2 4 8 16 32 64 128 256 512 1024" \
-t "3600" \
-p "4"

# partial match benchmark
./fmd_run.sh ../input/GRCh38-20-0.10b.fmdg \
-L "GRCh38-20-0.10b.fmdg_c_l2" \
-c "0 1 2 3 4 5 6 7 8 9 10 11 12" \
-r "16" \
-m "-1" \
-M "-1" \
-j "1" \
-l "2 4 8 16 32 64 128 256 512 1024" \
-t "3600" \
-p "4"
