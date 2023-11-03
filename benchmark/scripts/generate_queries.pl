#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

my @graph;  # List of lists representation
my @labels;
my %vertex_to_index;  # Map from vertex ID to index
my ($input_file, $query_length, $no_queries, $seed);

# Parse command-line arguments
GetOptions(
    'i=s' => \$input_file,
    'q=i' => \$query_length,
    'N=i' => \$no_queries,
    's=i' => \$seed
) or die("Error in command line arguments\n");

srand($seed) if defined $seed;

# Read the graph from the input file
open my $fh, '<', $input_file or die "Cannot open $input_file: $!";
while (<$fh>) {
    chomp;
    my @fields = split /\t/;
    if ($fields[0] eq 'V') {
        my $index = scalar(@labels);  # Current index
        $vertex_to_index{$fields[1]} = $index;
        push @labels, $fields[2];
        push @graph, [];  # Initialize an empty adjacency list for this vertex
    } elsif ($fields[0] eq 'E') {
        my $src_index = $vertex_to_index{$fields[1]};
        my $dst_index = $vertex_to_index{$fields[2]};
        push @{$graph[$src_index]}, $dst_index;
    }
}
close $fh;

# Function to generate a substring starting from a given vertex and offset
sub generate_substring {
    my ($start_vertex, $start_offset, $length) = @_;
    my @stack;
    push @stack, { vertex => $start_vertex, offset => $start_offset, substring => '' };

    while (@stack) {
        my $current = pop @stack;
        my $vertex = $current->{vertex};
        my $offset = $current->{offset};
        my $substring = $current->{substring};

        # If the substring length is reached, return it
        if (length($substring) == $length) {
            return $substring;
        }

        # If the current label contains 'N', backtrack
        if (index($labels[$vertex], 'N') != -1) {
            next;
        }

        # Extract the portion of the label from the offset
        my $remaining_label = substr($labels[$vertex], $offset);

        # If the remaining label can be added entirely
        if (length($substring) + length($remaining_label) <= $length) {
            $substring .= $remaining_label;
            if (defined $graph[$vertex]) {
                # Push all outgoing vertices to the stack
                for my $neighbor (@{$graph[$vertex]}) {
                    push @stack, { vertex => $neighbor, offset => 0, substring => $substring };
                }
            }
        } else {
            # If only a part of the remaining label is needed
            $substring .= substr($remaining_label, 0, $length - length($substring));
            return $substring;
        }
    }
    return '';
}

# Binary search to find the interval for the random number
sub binary_search {
    my ($arr, $rand) = @_;
    my ($low, $high) = (0, scalar(@$arr) - 1);
    while ($low <= $high) {
        my $mid = int(($low + $high) / 2);
        if ($rand < $arr->[$mid]) {
            $high = $mid - 1;
        } elsif ($mid == scalar(@$arr) - 1 || ($rand >= $arr->[$mid] && $rand < $arr->[$mid + 1])) {
            return $mid;
        } else {
            $low = $mid + 1;
        }
    }
    return -1;  # Should not reach here if the random number is valid
}

# Uniformly sample positions on the graph to generate substrings
sub sample_substrings {
    my ($n, $q) = @_;
    my @vertices = 0..$#labels;
    my $total_length = 0;
    $total_length += length($labels[$_]) for @vertices;

    my @cumulative_probabilities;
    my $cumulative = 0;
    for my $vertex (@vertices) {
        $cumulative += length($labels[$vertex]) / $total_length;
        push @cumulative_probabilities, $cumulative;
    }

    for (1..$n) {
        my $rand = rand();
        my $index = binary_search(\@cumulative_probabilities, $rand);
        my $selected_offset = int(rand(length($labels[$index])));

        my $substring = generate_substring($index, $selected_offset, $q);
        if ($substring) {
            print "$substring\n";
        } else {
            redo;
        }
    }
}

sample_substrings($no_queries, $query_length);

print "exit();\n"