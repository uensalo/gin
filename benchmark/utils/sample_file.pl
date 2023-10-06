#!/usr/bin/perl

use strict;
use warnings;

if (@ARGV != 1) {
    die "Usage: $0 <input_file>\n";
}

my $input_file = $ARGV[0];

if (!-e $input_file) {
    die "Error: File '$input_file' not found!\n";
}

open(my $fh, '<', $input_file) or die "Cannot open $input_file: $!";
my $s = <$fh>;
chomp $s;
close $fh;

for my $k (16, 32, 64, 128, 256, 512, 1024, 2048, 4096) {
    my $output_file = "${input_file}.${k}.txt";
    open(my $out, '>', $output_file) or die "Cannot open $output_file for writing: $!";
    
    for (1..65536) {
        my $sample;
        do {
            my $pos = int(rand(length($s) - $k + 1));
            $sample = substr($s, $pos, $k);
        } while ($sample =~ /^N+$/); # Resample if all N's
        
        print $out $sample . "\n";
    }
    
    close $out;
}
