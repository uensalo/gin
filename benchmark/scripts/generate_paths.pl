#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

my @graph;
my @labels;
my %vertex_to_index;
my ($input_file, $query_length);

GetOptions(
    'i=s' => \$input_file,
    'q=i' => \$query_length
) or die("Error in command line arguments\n");

open my $fh, '<', $input_file or die "Cannot open $input_file: $!";
while (<$fh>) {
    chomp;
    next if /^\s*$/; # Skip empty lines

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

sub generate_paths {
    my ($vertex, $offset_start, $offset_end, $consumable_start, $consumable_end, $current_path) = @_;
    my $label_length = length($labels[$vertex]);

    $current_path = "$current_path:$vertex";

    if ($label_length >= $consumable_end) {
        # Everything is consumed
        print "$offset_start,$offset_end;$current_path\n";
        return;
    } elsif ($consumable_end > $label_length && $label_length >= $consumable_start) {
        # Compute the end of the consumed interval
        my $push_end = $offset_start + $label_length - $consumable_start;
        print "$offset_start,$push_end;$current_path\n";
        $offset_start = $push_end + 1;
        $consumable_start = 1;
        $consumable_end -= $label_length;
    } else {
        # Consume symbols
        $consumable_start -= $label_length;
        $consumable_end -= $label_length;
    }

    # Recurse
    foreach my $neighbor (@{$graph[$vertex]}) {
        generate_paths($neighbor, $offset_start, $offset_end, $consumable_start, $consumable_end, "$current_path");
    }
}

for my $vertex (0..$#labels) {
    my $label_length = length($labels[$vertex]);
    my ($offset_start, $offset_end, $consumable_start, $consumable_end);

    if ($label_length >= $query_length) {
        $offset_start = $label_length - $query_length + 1;
        $offset_end = $label_length - 1;
        $consumable_start = 1;
        $consumable_end = $query_length - 1;
    } else {
        $offset_start = 0;
        $offset_end = $label_length - 1;
        $consumable_start = $query_length - $label_length;
        $consumable_end = $query_length - 1;
    }

    # Call generate_paths on each neighbor of the vertex
    foreach my $neighbor (@{$graph[$vertex]}) {
        generate_paths($neighbor, $offset_start, $offset_end, $consumable_start, $consumable_end, "$vertex");
    }
}
