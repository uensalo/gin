#!/bin/bash
#########################################################
# Transcriptomics
#########################################################
# query length grouped by cache depth
./fmd_run.sh ../input/gencode.v40.fmdg \
-L "gencode.v40.fmdg_c_l" \
-c "0 2 4 6 8 10" \
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
-m "-1 1 4 64 512 4096" \
-M "-1" \
-j "1" \
-l "10 30 100 300 1000" \
-t "3600" \
-p "4"

# query length grouped by permutation depth
./fmd_run.sh ../input/gencode.v40.fmdg \
-L "gencode.v40.fmdg_p_l"\
-c "10" \
-r "16" \
-m "-1" \
-M "-1" \
-j "1" \
-l "10 30 100 300 1000" \
-t "3600" \
-p "2 4 6 8"

# query length grouped by sampling rate & threads
./fmd_run.sh ../input/gencode.v40.fmdg \
-L "gencode.v40.fmdg_r_l"\
-c "10" \
-r "16 32 64" \
-m "-1" \
-M "-1" \
-j "1 4 8 16 32" \
-l "10 30 100 300 1000" \
-t "3600" \
-p "4"

#########################################################
# Pangenomics
#########################################################
# query length grouped by cache depth
./fmd_run.sh ../input/GRCh38-20-0.10b.fmdg \
-L "GRCh38-20-0.10b.fmdg_c_l" \
-c "0 2 4 6 8 10" \
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
-m "-1 1 4 64 512 4096" \
-M "-1" \
-j "1" \
-l "10 30 100 300 1000" \
-t "3600" \
-p "4"

# query length grouped by permutation depth
./fmd_run.sh ../input/GRCh38-20-0.10b.fmdg \
-L "GRCh38-20-0.10b.fmdg_p_l"\
-c "10" \
-r "16" \
-m "-1" \
-M "-1" \
-j "1" \
-l "10 30 100 300 1000" \
-t "3600" \
-p "2 4 6 8"

# query length grouped by sampling rate & threads
./fmd_run.sh ../input/GRCh38-20-0.10b.fmdg \
-L "GRCh38-20-0.10b.fmdg_r_l"\
-c "10" \
-r "16 32 64" \
-m "-1" \
-M "-1" \
-j "1 4 8 16 32" \
-l "10 30 100 300 1000" \
-t "3600" \
-p "4"
