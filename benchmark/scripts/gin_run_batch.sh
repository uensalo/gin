#!/bin/bash

# Executables
GIN_DIR=../bin
QUERY_SCRIPT=generate_queries.pl

# If hard queries are requested
HARD_SUFFIX=""

# Set the hard-coded directories
LOG_DIR=../log
QUERY_DIR=../res/query
INDEX_DIR=../res/index
PERMUTATION_DIR=../res/permutation
CACHE_DIR=../res/cache

# Hard-coded parameters
NO_QUERIES=65536
SEED=420
TEMPERATURE=1e2
COOLING=0.99
BATCH_SIZE_MIN=8
PERMUTATION_NUM_THREADS=16

# Variables to manage the decode option
DECODE_FLAG=""
DECODE_SUFFIX=""
EXPERIMENT_NAME=""

# Read input parameters
INPUT_FILE=$1
shift
while getopts ":c:r:m:M:j:l:t:p:dL:h" opt; do
  case $opt in
    c) IFS=' ' read -r -a CACHE_DEPTHS <<< "$OPTARG";;
    r) IFS=' ' read -r -a SAMPLE_RATES <<< "$OPTARG";;
    m) IFS=' ' read -r -a MAX_FORKS <<< "$OPTARG";;
    M) IFS=' ' read -r -a MAX_MATCHES <<< "$OPTARG";;
    j) IFS=' ' read -r -a QUERY_THREADS <<< "$OPTARG";;
    l) IFS=' ' read -r -a LENGTHS <<< "$OPTARG";;
    t) IFS=' ' read -r -a PERMUTATION_TIMES <<< "$OPTARG";;
    p) IFS=' ' read -r -a PERMUTATION_DEPTHS <<< "$OPTARG";;
    d)
       DECODE_FLAG="--decode";
       DECODE_SUFFIX="_decode";;
    L) EXPERIMENT_NAME="$OPTARG";;
    h)
       QUERY_SCRIPT=generate_queries_hard.pl; # Switch scripts.
       HARD_SUFFIX="_hard";;
    \?) echo "Invalid option -$OPTARG" >&2; exit 1;;
  esac
done

# Get the base name of the input file
BASENAME=$(basename "$INPUT_FILE" .ging)

# Create log subdirectory for the experiment
LOG_DIR="$LOG_DIR/$EXPERIMENT_NAME"
mkdir -p $LOG_DIR

# Pre-generate the query files (if they don't exist)
for LENGTH in "${LENGTHS[@]}"
do
  if [[ -z "$HARD_SUFFIX" ]]; then
    CURRENT_QUERY_SCRIPT="generate_queries.pl"
  elif [[ "$LENGTH" -le 512 ]]; then
    CURRENT_QUERY_SCRIPT="generate_queries_hard.pl"
  else
    CURRENT_QUERY_SCRIPT="generate_queries_hard_long.pl"
  fi
  QUERY_FILE="$QUERY_DIR/${BASENAME}${HARD_SUFFIX}_query_length_${LENGTH}.ginq"
  if [[ ! -f $QUERY_FILE ]]; then
    ./$CURRENT_QUERY_SCRIPT -i "$INPUT_FILE" -q "$LENGTH" -N $NO_QUERIES -s $SEED > "$QUERY_FILE" &
  fi
done
wait

# Compute all the permutations in parallel (2 at a time)
for PERMUTATION_TIME in "${PERMUTATION_TIMES[@]}"
do
  for PERMUTATION_DEPTH in "${PERMUTATION_DEPTHS[@]}"
  do
    if [[ "$PERMUTATION_TIME" -ne 0 && "$PERMUTATION_DEPTH" -ne 0 ]]; then
      PERMUTATION_FILE="$PERMUTATION_DIR/${BASENAME}_ptime_${PERMUTATION_TIME}_pdepth_${PERMUTATION_DEPTH}_threads_${PERMUTATION_NUM_THREADS}.ginp"
      PERM_LOG_FILE="$LOG_DIR/perm_log_ptime_${PERMUTATION_TIME}_pdepth_${PERMUTATION_DEPTH}_threads_${PERMUTATION_NUM_THREADS}.txt"
      if [[ ! -f $PERMUTATION_FILE ]]; then
        $GIN_DIR/gin permutation -i "$INPUT_FILE" -t "$PERMUTATION_TIME" -u "$PERMUTATION_TIME" -e $TEMPERATURE -c $COOLING -d "$PERMUTATION_DEPTH" -o "$PERMUTATION_FILE" -j $PERMUTATION_NUM_THREADS -v 2>> "$PERM_LOG_FILE" &
      fi
      # Make sure only two tasks are running at the same time
      while [ $(jobs | wc -l) -ge 2 ]; do
        sleep 1
      done
    fi
  done
done

# Compute all the indexes in parallel over sample rates (4 at a time)
for SAMPLE_RATE_IDX in $(seq 0 3 $((${#SAMPLE_RATES[@]} - 1)))
do
  for PERMUTATION_TIME in "${PERMUTATION_TIMES[@]}"
  do
    for PERMUTATION_DEPTH in "${PERMUTATION_DEPTHS[@]}"
    do
      PERMUTATION_FLAG=""
      PERMUTATION_FILE=""
      if [[ "$PERMUTATION_TIME" -ne 0 && "$PERMUTATION_DEPTH" -ne 0 ]]; then
        PERMUTATION_FLAG="-p"
        PERMUTATION_FILE="$PERMUTATION_DIR/${BASENAME}_ptime_${PERMUTATION_TIME}_pdepth_${PERMUTATION_DEPTH}_threads_${PERMUTATION_NUM_THREADS}.ginp"
      fi
      for IDX in $(seq "$SAMPLE_RATE_IDX" $(($SAMPLE_RATE_IDX + 3)))
      do
        SAMPLE_RATE=${SAMPLE_RATES[$IDX]}
        INDEX_FILE="$INDEX_DIR/${BASENAME}_index_ptime_${PERMUTATION_TIME}_pdepth_${PERMUTATION_DEPTH}_sampling_rate_${SAMPLE_RATE}.gini"
        INDEX_LOG_FILE="$LOG_DIR/index_log_ptime_${PERMUTATION_TIME}_pdepth_${PERMUTATION_DEPTH}_sampling_rate_${SAMPLE_RATE}.txt"
        if [[ ! -f $INDEX_FILE && $IDX -lt ${#SAMPLE_RATES[@]} ]]; then
          $GIN_DIR/gin index -i "$INPUT_FILE" $PERMUTATION_FLAG "$PERMUTATION_FILE" -o "$INDEX_FILE" -s "$SAMPLE_RATE" -r "$SAMPLE_RATE" -v 2>> "$INDEX_LOG_FILE" &
        fi
      done
      wait
    done
  done
done

# Compute the caches in parallel (4 at a time)
for CACHE_DEPTH_IDX in $(seq 0 3 $((${#CACHE_DEPTHS[@]} - 1)))
do
  for IDX in $(seq "$CACHE_DEPTH_IDX" $(($CACHE_DEPTH_IDX + 3)))
  do
    CACHE_DEPTH=${CACHE_DEPTHS[$IDX]}
    if [[ "$CACHE_DEPTH" -ne "0" && $IDX -lt ${#CACHE_DEPTHS[@]} ]]; then
      for PERMUTATION_TIME in "${PERMUTATION_TIMES[@]}"
      do
        for PERMUTATION_DEPTH in "${PERMUTATION_DEPTHS[@]}"
        do
          FIRST_SAMPLE_RATE=${SAMPLE_RATES[0]}
          CACHE_FILE="$CACHE_DIR/${BASENAME}_cache_ptime_${PERMUTATION_TIME}_pdepth_${PERMUTATION_DEPTH}_depth_${CACHE_DEPTH}.ginc"
          CACHE_LOG_FILE="$LOG_DIR/cache_log_ptime_${PERMUTATION_TIME}_pdepth_${PERMUTATION_DEPTH}_depth_${CACHE_DEPTH}.txt"
          INDEX_FILE="$INDEX_DIR/${BASENAME}_index_ptime_${PERMUTATION_TIME}_pdepth_${PERMUTATION_DEPTH}_sampling_rate_${FIRST_SAMPLE_RATE}.gini"
          if [[ ! -f $CACHE_FILE ]]; then
            $GIN_DIR/gin query cache -r "$INDEX_FILE" -o "$CACHE_FILE" -j 8 -c "$CACHE_DEPTH" -v 2>> "$CACHE_LOG_FILE" &
          fi
        done
      done
    fi
  done
  wait
done

# Querying
for CACHE_DEPTH in "${CACHE_DEPTHS[@]}"
do
  for SAMPLE_RATE in "${SAMPLE_RATES[@]}"
  do
    for PERMUTATION_TIME in "${PERMUTATION_TIMES[@]}"
    do
      for PERMUTATION_DEPTH in "${PERMUTATION_DEPTHS[@]}"
      do
        CACHE_FLAG=""
        CACHE_FILE=""
        if [[ "$CACHE_DEPTH" -ne "0" ]]; then
          CACHE_FLAG="-C"
          CACHE_FILE="$CACHE_DIR/${BASENAME}_cache_ptime_${PERMUTATION_TIME}_pdepth_${PERMUTATION_DEPTH}_depth_${CACHE_DEPTH}.ginc"
        fi
        for FORK in "${MAX_FORKS[@]}"
        do
          for MATCH in "${MAX_MATCHES[@]}"
          do
            for THREAD in "${QUERY_THREADS[@]}"
            do
              for LENGTH in "${LENGTHS[@]}"
              do
                LOG_FILE="$LOG_DIR/find_log_ptime_${PERMUTATION_TIME}_pdepth_${PERMUTATION_DEPTH}_sampling_rate_${SAMPLE_RATE}_cache_${CACHE_DEPTH}_fork_${FORK}_match_${MATCH}_threads_${THREAD}_length_${LENGTH}${DECODE_SUFFIX}.txt"
                # Run the find command with the queries, redirecting stderr to the log file
                QUERY_FILE="$QUERY_DIR/${BASENAME}${HARD_SUFFIX}_query_length_${LENGTH}.ginq"
                INDEX_FILE="$INDEX_DIR/${BASENAME}_index_ptime_${PERMUTATION_TIME}_pdepth_${PERMUTATION_DEPTH}_sampling_rate_${SAMPLE_RATE}.gini"
                BATCH_SIZE=$(( THREAD > BATCH_SIZE_MIN ? THREAD : BATCH_SIZE_MIN ))
                $GIN_DIR/gin query find -r "$INDEX_FILE" $CACHE_FLAG "$CACHE_FILE" -i "$QUERY_FILE" -j "$THREAD" -b $BATCH_SIZE -m "$FORK" -M "$MATCH" $DECODE_FLAG -o "/dev/null" -v 2>> "$LOG_FILE"
              done
            done
          done
        done
      done
    done
  done
done