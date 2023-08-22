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
BASENAME=$(basename "$INPUT_GRAPH" .fmdg)

PLOT_COMMAND="set terminal pngcairo size 800,600; \
set output '$LOG_DIR/plot_${BASENAME}_k${DEPTH}.png'; \
set xlabel 'Path Length'; \
set ylabel 'Frequency'; \
plot "

for ((k=1; k<=$DEPTH; k++)); do
    PATH_LOG="$LOG_DIR/paths_${BASENAME}_k$k.txt"
    ./path_generation.pl -i "$INPUT_GRAPH" -q "$k" > "$PATH_LOG"

    if [ "$VERTEX_COMPLEXITY" = true ]; then
        VERTEX_LOG="$LOG_DIR/vertex_complexity_${BASENAME}_k$k.txt"
        cat "$PATH_LOG" | ./vertex_complexity.pl > "$VERTEX_LOG"
    fi

    if [ "$EDGE_COMPLEXITY" = true ]; then
        EDGE_LOG="$LOG_DIR/edge_complexity_${BASENAME}_k$k.txt"
        cat "$PATH_LOG" | ./edge_complexity.pl > "$EDGE_LOG"
    fi

    HISTOGRAM_LOG="$LOG_DIR/histogram_${BASENAME}_k$k.txt"
    ./path_histogram.pl -i "$PATH_LOG" > "$HISTOGRAM_LOG"

    # Append to the gnuplot commands
    PLOT_COMMAND="$PLOT_COMMAND '$HISTOGRAM_LOG' using 1:2 with linespoints title 'k=$k', "
done

# Remove the last comma and space
PLOT_COMMAND=${PLOT_COMMAND%, }
echo "$PLOT_COMMAND" | gnuplot
