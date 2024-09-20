#!/bin/bash
# Parallelism on the pangenome
./gin_run_fast.sh ../input/GRCh38-20-0.10b.ging \
-L "GRCh38-20-0.10b.ging_parallelism" \
-c "12" \
-r "32" \
-m "-1" \
-M "-1" \
-j "1 2 4 8 10 16 32 48 64" \
-l "16 32 64 128 256 512 1024 2048 4096" \
-t "3600" \
-p "12" \
-d \
-h

# Parallelism on the transcriptome
./gin_run_fast.sh ../input/gencode.v40.ging \
-L "gencode.v40.ging_parallelism" \
-c "12" \
-r "32" \
-m "-1" \
-M "-1" \
-j "1 2 4 8 10 16 32 48 64" \
-l "16 32 64 128 256 512 1024 2048 4096" \
-t "3600" \
-p "12" \
-d \
-h