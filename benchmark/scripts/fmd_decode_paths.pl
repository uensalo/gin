#!/usr/bin/perl
use strict;
use warnings;

my %graph;
my %labels;

sub parse_graph {
    my ($filename) = @_;
    open my $fh, '<', $filename or die "Cannot open $filename: $!";
    while (<$fh>) {
        chomp;
        my @fields = split /\t/;
        if ($fields[0] eq 'V') {
            $labels{$fields[1]} = $fields[2];
        } elsif ($fields[0] eq 'E') {
            push @{$graph{$fields[1]}}, $fields[2];
        }
    }
    close $fh;
}

sub find_paths {
    my ($string, $v, $o, $path, $start_o, $string_pos) = @_;
    $string_pos ||= 0;

    my $label = substr($labels{$v}, $o);
    my $label_len = length($label);
    my $remaining_len = length($string) - $string_pos;

    if (substr($string, $string_pos, $label_len) eq $label) {
        if ($string_pos + $label_len == length($string)) {
            print OUT "\t($start_o,$o);" . join(':', @$path, $v) . "\n";
        } else {
            foreach my $neighbor (@{$graph{$v}}) {
                find_paths($string, $neighbor, 0, [@$path, $v], $start_o, $string_pos + $label_len);
            }
        }
    } elsif (substr($label, 0, $remaining_len) eq substr($string, $string_pos)) {
        my $end_o = $o + $remaining_len;
        print OUT "\t($start_o,$end_o);" . join(':', @$path, $v) . "\n";
    }
}

sub usage {
    print << "END_USAGE";
Description:
    This program serves as an example path decoder. It processes the output of the 'fmd query find --decode' command. For each path root output for a given query, the script traverses the input graph file and identifies all walks that produce the corresponding string, starting from a specific offset.

    Input Format (fmd query find --decode output format):
    The input should adhere to the following semantics and syntax:
    <string>:
    \t(v:<integer>,o:<integer>)
    ... (repeated for each vertex-offset pair)
    ... (repeated for each string)

    Output Format:
    The script produces output in the following format:
    <string>:
    \t(o1,o2);v1:v2:...:vN
    ... (repeated for each path)
    ... (repeated for each string)

    Semantics:
    o1: Represents the offset into the label of the first vertex (v1) where the string starts.
    o2: Represents the offset into the label of the last vertex (vN) where the string ends.

    Example decoding pipeline call:
    ./fmd query find -r <fmd_index> -i <query_file> -C <cache_file> --decode 2>/dev/null | fmd_decode_paths.pl -r <fmdg_file> -o <decoded_paths_file>

Usage:
    $0 -r <graph_file_name> [-i <input_file_name>] [-o <output_file_name>]

Options:
    -r <graph_file_name>    : Specify the graph file name. (Required)
    -i <input_file_name>    : Specify the input file name. Default is stdin.
    -o <output_file_name>   : Specify the output file name. Default is stdout.
    --help                  : Display this help message.

END_USAGE
    exit(1);
}

my ($graph_file, $input_file, $output_file) = ('', '-', '-');
while (@ARGV) {
    my $arg = shift @ARGV;
    if ($arg eq '-r') {
        $graph_file = shift @ARGV;
    } elsif ($arg eq '-i') {
        $input_file = shift @ARGV || '-';
    } elsif ($arg eq '-o') {
        $output_file = shift @ARGV || '-';
    } elsif ($arg eq '--help') {
        usage();
    } else {
        print "Unknown option: $arg\n";
        usage();
    }
}

unless ($graph_file) {
    print "Error: -r <graph_file_name> is a required parameter.\n";
    usage();
}

parse_graph($graph_file);

open my $fh, ($input_file eq '-' ? '<-' : '<' . $input_file) or die "Cannot open $input_file: $!";
open OUT, ($output_file eq '-' ? '>-' : '>' . $output_file) or die "Cannot open $output_file: $!";

my $current_string;
while (<$fh>) {
    chomp;
    if (/^(\w+):$/) {
        $current_string = $1;
        print OUT "$current_string:\n";
    } elsif (my ($v, $o) = /\(v:(\d+),o:(\d+)\)/) {
        find_paths($current_string, $v, $o, [], $o);
    }
}

close $fh;
close OUT;
