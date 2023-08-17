#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

my @graph;  # List of lists representation
my @labels;
my %vertex_to_index;  # Map from vertex ID to index
my ($input_file, $query_length);
my %memo;  # For memoization

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
        my $index = scalar(@labels);
        $vertex_to_index{$fields[1]} = $index;
        push @labels, $fields[2];
        push @graph, [];
    } elsif ($fields[0] eq 'E') {
        my $src_index = $vertex_to_index{$fields[1]};
        my $dst_index = $vertex_to_index{$fields[2]};
        push @{$graph[$src_index]}, $dst_index;
    }
}
close $fh;

# Function to generate a substring starting from a given vertex and offset
sub generate_substring {
    my ($vertex, $offset, $length, $current_string) = @_;

    return [$current_string] if length($current_string) == $length;
    return [] if index($labels[$vertex], 'N') != -1;

    my $key = "$vertex:$offset:" . length($current_string);
    return $memo{$key} if exists $memo{$key};

    my $remaining_label = substr($labels[$vertex], $offset);
    my @results;

    if (length($current_string) + length($remaining_label) <= $length) {
        $current_string .= $remaining_label;
        if (defined $graph[$vertex]) {
            for my $neighbor (@{$graph[$vertex]}) {
                push @results, @{generate_substring($neighbor, 0, $length, $current_string)};
            }
        }
    } else {
        $current_string .= substr($remaining_label, 0, $length - length($current_string));
        push @results, $current_string;
    }

    $memo{$key} = \@results;
    return \@results;
}

# Generate all possible substrings of length q
sub generate_all_substrings {
    my ($q) = @_;
    my %unique_strings;

    for my $vertex_index (0..$#labels) {
        for my $offset (0..length($labels[$vertex_index]) - 1) {
            my $substrings = generate_substring($vertex_index, $offset, $q, "");
            for my $substring (@$substrings) {
                $unique_strings{$substring} = 1 if length($substring) == $q;
            }
        }
    }

    return keys %unique_strings;
}

my @all_unique_strings = generate_all_substrings($query_length);
foreach my $string (sort @all_unique_strings) {
    print "$string\n";
}
