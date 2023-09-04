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
    }
    elsif (substr($label, 0, $remaining_len) eq substr($string, $string_pos)) {
        my $end_o = $o + $remaining_len;
        print OUT "\t($start_o,$end_o);" . join(':', @$path, $v) . "\n";
    }
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
    }
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
