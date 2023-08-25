#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use File::Temp qw/ tempfile /;

my @graph;
my @labels;
my %vertex_to_index;
my ($input_file, $query_length, $num_queries, $seed);

GetOptions(
    'i=s' => \$input_file,
    'q=i' => \$query_length,
    'N=i' => \$num_queries,
    's=i' => \$seed   # For backward compatibility
) or die("Error in command line arguments\n");

if ($query_length <= 512) {
    die("Please use generate_queries_hard.pl for query lengths <= 512.\n");
}

open my $fh, '<', $input_file or die "Cannot open $input_file: $!";
while (<$fh>) {
    chomp;
    next if /^\s*$/;

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

# Temporary file for writing paths
my ($tmp_fh, $tmp_filename) = tempfile();

sub generate_paths {
    my ($vertex, $offset_start, $offset_end, $consumable_start, $consumable_end, $current_path, $current_weight) = @_;
    my $label_length = length($labels[$vertex]);
    my $path_length = scalar(split(/:/, $current_path));

    return if index($labels[$vertex], 'N') != -1;  # Skip paths containing 'N'

    $current_path = "$current_path:$vertex";
    $current_weight += $label_length / ($path_length);

    if ($label_length >= $consumable_end) {
        # Entire vertex label is consumed
        print $tmp_fh "$offset_start,$offset_end;$current_path;$current_weight\n";
        return;
    } elsif ($consumable_end > $label_length && $label_length >= $consumable_start) {
        # Vertex label partially satisfies the path length
        my $push_end = $offset_start + $label_length - $consumable_start;
        print $tmp_fh "$offset_start,$push_end;$current_path;$current_weight\n";
        $offset_start = $push_end + 1;
        $consumable_start = 1;
        $consumable_end -= $label_length;
    } else {
        # Entire vertex label is needed but doesn't fulfill required length
        $consumable_start -= $label_length;
        $consumable_end -= $label_length;
    }

    # Recurse
    foreach my $neighbor (@{$graph[$vertex]}) {
        generate_paths($neighbor, $offset_start, $offset_end, $consumable_start, $consumable_end, $current_path, $current_weight);
    }
}

# Generate paths with a hardcoded query_length of 512
for my $vertex (0..$#labels) {
    my $label_length = length($labels[$vertex]);
    my ($offset_start, $offset_end, $consumable_start, $consumable_end);

    if ($label_length >= 512) {
        $offset_start = $label_length - 512 + 1;
        $offset_end = $label_length - 1;
        $consumable_start = 1;
        $consumable_end = 512 - 1;
    } else {
        $offset_start = 0;
        $offset_end = $label_length - 1;
        $consumable_start = 512 - $label_length;
        $consumable_end = 512 - 1;
    }

    # Skip if the first label already contains an 'N'
    next if index(substr($labels[$vertex], $offset_start, $offset_end - $offset_start + 1), 'N') != -1;

    # Call generate_paths on each neighbor of the vertex
    foreach my $neighbor (@{$graph[$vertex]}) {
        generate_paths($neighbor, $offset_start, $offset_end, $consumable_start, $consumable_end, "$vertex", $label_length);
    }
}

# Reset the temporary file's position to the beginning for reading
seek $tmp_fh, 0, 0;

# Read back from the temporary file, sort by weight in ascending order, and extract extended paths
open my $sorted_fh, "-|", "sort -t';' -k3,3n $tmp_filename" or die "Failed to sort: $!";
my %unique_strings;

while (<$sorted_fh>) {
    my ($offsets, $vertices, $weight) = split /;/;
    my ($start, $end) = split /,/, $offsets;
    my @vertex_ids = split /:/, $vertices;

    # Construct the query string from the path
    my $query_string = substr($labels[$vertex_ids[0]], $start, $end-$start+1);
    for my $i (1..$#vertex_ids) {
        $query_string .= $labels[$vertex_ids[$i]];
    }

    # Extend the query string until it reaches the required length
    my $current_vertex = $vertex_ids[-1];
    while (length($query_string) < $query_length) {
        my @sorted_neighbors = sort { length($labels[$a]) <=> length($labels[$b]) } @{$graph[$current_vertex]};
        my $selected_neighbor = shift @sorted_neighbors;
        $query_string .= $labels[$selected_neighbor];
        $current_vertex = $selected_neighbor;
    }
    $query_string = substr($query_string, 0, $query_length);
    $unique_strings{$query_string} = 1;

    last if scalar(keys %unique_strings) >= $num_queries;
}

foreach my $string (keys %unique_strings) {
    print "$string\n";
}

# Cleanup: Remove the temporary file
unlink $tmp_filename;

print "exit();\n"
