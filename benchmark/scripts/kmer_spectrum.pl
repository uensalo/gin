#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

my @graph;  # List of lists representation
my @labels;
my %vertex_to_index;  # Map from vertex ID to index
my ($input_file, $query_length);

# Parse command-line arguments
GetOptions(
    'i=s' => \$input_file,
    'q=i' => \$query_length
) or die("Error in command line arguments\n");

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

# Generate all possible substrings of length q
sub generate_all_substrings {
    my ($q) = @_;
    my %unique_strings;

    for my $vertex_index (0..$#labels) {
        for my $offset (0..length($labels[$vertex_index]) - 1) {
            my $substring = generate_substring($vertex_index, $offset, $q);
            $unique_strings{$substring} = 1 if length($substring) == $q;
        }
    }

    return keys %unique_strings;
}

my @all_unique_strings = generate_all_substrings($query_length);
foreach my $string (sort @all_unique_strings) {
    print "$string\n";
}
