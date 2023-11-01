#!/bin/bash

while getopts "i:L:d:ve" arg; do
    case $arg in
        i) INPUT_GRAPH="$OPTARG";;
        L) EXPERIMENT_NAME="$OPTARG";;
        d) DEPTH="$OPTARG";;
        v) VERTEX_COMPLEXITY=true;;
        e) EDGE_COMPLEXITY=true;;
        *) exit 1;;
    esac
done

LOG_DIR="../log/$EXPERIMENT_NAME"
mkdir -p "$LOG_DIR"
BASENAME=$(basename "$INPUT_GRAPH" .ging)

# Define the inner loop as a function
compute() {
    local k="$1"
    local INPUT_GRAPH="$2"
    local LOG_DIR="$3"
    local BASENAME="$4"
    local VERTEX_COMPLEXITY="$5"
    local EDGE_COMPLEXITY="$6"

    PATH_LOG="$LOG_DIR/paths_${BASENAME}_k$k.txt"
    perl ./generate_paths.pl -i "$INPUT_GRAPH" -q "$k" > "$PATH_LOG"

    if [ "$VERTEX_COMPLEXITY" = true ]; then
        VERTEX_LOG="$LOG_DIR/vertex_complexity_${BASENAME}_k$k.txt"
        cat "$PATH_LOG" | perl ./vertex_complexity.pl > "$VERTEX_LOG"
    fi

    if [ "$EDGE_COMPLEXITY" = true ]; then
        EDGE_LOG="$LOG_DIR/edge_complexity_${BASENAME}_k$k.txt"
        cat "$PATH_LOG" | perl ./edge_complexity.pl > "$EDGE_LOG"
    fi

    HISTOGRAM_LOG="$LOG_DIR/histogram_${BASENAME}_k$k.txt"
    perl ./generate_path_histogram.pl -i "$PATH_LOG" > "$HISTOGRAM_LOG"
}
export -f compute


count=0
for ((k=1; k<=$DEPTH; k++)); do
    compute "$k" "$INPUT_GRAPH" "$LOG_DIR" "$BASENAME" "$VERTEX_COMPLEXITY" "$EDGE_COMPLEXITY" &
    ((count++))
    if ((count % 4 == 0)); then
        wait
    fi
done
wait

PLOT_COMMAND="set terminal pngcairo size 800,600; \
set output '$LOG_DIR/plot_${BASENAME}_k${DEPTH}.png'; \
set xlabel 'Path Length'; \
set ylabel 'Frequency'; \
plot "

for ((k=2; k<=$DEPTH; k++)); do
    HISTOGRAM_LOG="$LOG_DIR/histogram_${BASENAME}_k$k.txt"
    PLOT_COMMAND="$PLOT_COMMAND '$HISTOGRAM_LOG' using 1:2 with linespoints title 'k=$k', "
done

# Remove the last comma and space
PLOT_COMMAND=${PLOT_COMMAND%, }
echo "$PLOT_COMMAND" | gnuplot
