#!/bin/bash
cleanup() {
    jobs -p | xargs kill -s SIGTERM
}
trap "cleanup" SIGINT SIGTERM

if [[ "$#" -lt 8 ]]; then
    echo "Usage: $0 -d <depth> -a <alphabet> -i <input_file> -c <cache_depth> -D (optional decode flag) -j <num_threads> -L <experiment_name>"
    exit 1
fi

DECODE_FLAG=""

while getopts ":d:a:i:c:Dj:L:" opt; do
    case $opt in
        d) DEPTH="$OPTARG";;
        a) ALPHABET="$OPTARG";;
        i) INPUT_FILE="$OPTARG";;
        c) CACHE_DEPTH="$OPTARG";;
        D) DECODE_FLAG="--decode";;
        j) THREAD="$OPTARG";;
        L) EXPERIMENT_NAME="$OPTARG";;
        \?) echo "Invalid option -$OPTARG" >&2; exit 1;;
        :) echo "Option -$OPTARG requires an argument." >&2; exit 1;;
    esac
done

# Directories, Constants and Executables
BASENAME=$(basename "$INPUT_FILE" .fmdg)
LOG_DIR=../log/$EXPERIMENT_NAME
mkdir -p $LOG_DIR

FMD_DIR=../bin
QUERY_SCRIPT=generate_all_strings.pl
INDEX_DIR=../res/index
PERMUTATION_DIR=../res/permutation
CACHE_DIR=../res/cache
PERMUTATION_TIME=3600
PERMUTATION_DEPTH=4
PERMUTATION_NUM_THREADS=16
SAMPLE_RATE=32

# Log file paths
PERM_LOG_FILE="$LOG_DIR/perm_log_ptime_${PERMUTATION_TIME}_pdepth_${PERMUTATION_DEPTH}_threads_${PERMUTATION_NUM_THREADS}.txt"
INDEX_LOG_FILE="$LOG_DIR/index_log_ptime_${PERMUTATION_TIME}_pdepth_${PERMUTATION_DEPTH}_sampling_rate_${SAMPLE_RATE}.txt"
CACHE_LOG_FILE="$LOG_DIR/cache_log_ptime_${PERMUTATION_TIME}_pdepth_${PERMUTATION_DEPTH}_depth_${CACHE_DEPTH}.txt"

# Checks
INDEX_FILE="$INDEX_DIR/${BASENAME}_index_ptime_${PERMUTATION_TIME}_pdepth_${PERMUTATION_DEPTH}_sampling_rate_${SAMPLE_RATE}.fmdi"
if [[ ! -f $INDEX_FILE ]]; then
    PERMUTATION_FILE="$PERMUTATION_DIR/${BASENAME}_ptime_${PERMUTATION_TIME}_pdepth_${PERMUTATION_DEPTH}_threads_${PERMUTATION_NUM_THREADS}.fmdp"
    if [[ ! -f $PERMUTATION_FILE ]]; then
        $FMD_DIR/fmd permutation -i "$INPUT_FILE" -t "$PERMUTATION_TIME" -u "$PERMUTATION_TIME" -d "$PERMUTATION_DEPTH" -o "$PERMUTATION_FILE" -j $PERMUTATION_NUM_THREADS 2>> "$PERM_LOG_FILE"
    fi
    $FMD_DIR/fmd index -i "$INPUT_FILE" -p "$PERMUTATION_FILE" -o "$INDEX_FILE" -s "$SAMPLE_RATE" -r "$SAMPLE_RATE" 2>> "$INDEX_LOG_FILE"
fi

CACHE_FILE="$CACHE_DIR/${BASENAME}_cache_ptime_${PERMUTATION_TIME}_pdepth_${PERMUTATION_DEPTH}_depth_${CACHE_DEPTH}.fmdc"
if [[ ! -f $CACHE_FILE ]]; then
    $FMD_DIR/fmd query cache -r "$INDEX_FILE" -o "$CACHE_FILE" -j $THREAD -c "$CACHE_DEPTH" 2>> "$CACHE_LOG_FILE"
fi

# Actual benchmarking script
MAX_CONCURRENT_BENCHMARKS=4
running=0

for (( l=1; l<=$DEPTH; l++ )); do
    (perl $QUERY_SCRIPT -l $l -a "$ALPHABET" | $FMD_DIR/fmd query find -r "$INDEX_FILE" -C "$CACHE_FILE" -j "$THREAD" $DECODE_FLAG -v 2>> "$LOG_DIR/stats_kmer_length_${l}.txt" | awk '/^[^\t]/ && !/c:0/ {print}' > "$LOG_DIR/kmer_length_$l.txt") &
    running=$(($running + 1))

    # If the number of concurrent benchmarks hits the maximum limit, wait for any one of them to finish
    if [[ $running -eq $MAX_CONCURRENT_BENCHMARKS ]]; then
        wait -n
        running=$(($running - 1))
    fi
done

wait
