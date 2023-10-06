#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

my $input_file;

GetOptions('i=s' => \$input_file) or die("Error in command line arguments\n");

my $fh;
if (defined $input_file) {
    open $fh, '<', $input_file or die "Cannot open $input_file: $!";
} else {
    $fh = \*STDIN;
}

my %histogram;

while (<$fh>) {
    chomp;
    my ($offsets, $path) = split(/;/, $_, 2);
    my $length = (() = $path =~ /:/g) + 1;
    $histogram{$length}++;
}

close $fh if defined $input_file;

for my $length (sort { $a <=> $b } keys %histogram) {
    print "$length $histogram{$length}\n";
}
