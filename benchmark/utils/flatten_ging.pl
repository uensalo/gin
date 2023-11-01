#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

my ($input_file, $output_file);

# Process command line options
GetOptions(
    'i=s' => \$input_file,
    'o=s' => \$output_file
) or die("Error in command line arguments\n");

unless ($input_file && $output_file) {
    die "Usage: $0 -i <input_file> -o <output_file>\n";
}

open my $in_fh, '<', $input_file or die "Cannot open $input_file: $!";
open my $out_fh, '>', $output_file or die "Cannot open $output_file: $!";

while (<$in_fh>) {
    chomp;
    my @fields = split /\t/;
    if ($fields[0] eq 'V') {
        print $out_fh $fields[2];
    }
}
print $out_fh "\n";

close $in_fh;
close $out_fh;

