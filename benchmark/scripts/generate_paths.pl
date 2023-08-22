#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

my @graph;
my @labels;
my %vertex_to_index;
my ($input_file, $query_length);
my %memo;

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

# Function to generate paths for a given vertex and offset
sub generate_paths {
    my ($vertex, $offset, $length, $current_path) = @_;

    return [$current_path] if $length <= 0;

    my $key = "$vertex:$offset:$length";
    return $memo{$key} if exists $memo{$key};

    my $remaining_label_length = length($labels[$vertex]) - $offset;
    my @paths;

    if ($length <= $remaining_label_length) {
        push @paths, $current_path;
    } else {
        if (defined $graph[$vertex]) {
            for my $neighbor (@{$graph[$vertex]}) {
                my $new_length = $length - $remaining_label_length;
                push @paths, @{generate_paths($neighbor, 0, $new_length, "$current_path:$neighbor")};
            }
        }
    }

    $memo{$key} = \@paths;
    return \@paths;
}

for my $vertex_index (0..$#labels) {
    my $start_offset = length($labels[$vertex_index]) - $query_length + 1;
    $start_offset = 0 if $start_offset < 0; # Make sure the offset isn't negative

    for my $offset ($start_offset .. length($labels[$vertex_index]) - 1) {
        my $paths = generate_paths($vertex_index, $offset, $query_length, "$offset,$vertex_index");
        for my $path (@$paths) {
            print "$path\n";
        }
    }
}
