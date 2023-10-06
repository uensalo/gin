#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

# A hash to store the complexity of each vertex
my %complexity;
my $input_file;

GetOptions(
    'i=s' => \$input_file
) or die("Error in command line arguments\n");

# Function to calculate the width of an interval
sub get_interval_width {
    my ($interval) = @_;
    my ($start, $end) = split /,/, $interval;
    $end = $start unless $end; # If only one number is given, start = end
    return $end - $start + 1;
}

my $fh;
if ($input_file) {
    open($fh, '<', $input_file) or die "Cannot open $input_file: $!";
} else {
    $fh = *STDIN;
}

while (<$fh>) {
    chomp;
    my ($offsets, $path) = split /;/;
    my @offset_intervals = split /:/, $offsets;
    my @vertices = split /:/, $path;

    # Calculate the total width of the intervals
    my $total_width = 0;
    for my $interval (@offset_intervals) {
        $total_width += get_interval_width($interval);
    }

    # Update the complexity of each vertex
    for my $vertex (@vertices) {
        $complexity{$vertex} += $total_width;
    }
}

close($fh) if $input_file;

# Print the complexity of each vertex
for my $vertex (sort {$a <=> $b} keys %complexity) {
    print "$vertex:$complexity{$vertex}\n";
}
