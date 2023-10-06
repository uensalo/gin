#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use List::Util qw(sum);
use File::stat;
use File::Basename;

my ($input_file, $latex_table);
GetOptions('i=s' => \$input_file, 'T' => \$latex_table);

my (@lengths, %indegree, %outdegree, $size, $vertices, $edges);
my ($base, $dir, $ext) = fileparse($input_file, qr/\.[^.]*/);

my $fh;
if ($input_file) {
    open $fh, '<', $input_file or die "Cannot open $input_file: $!";
    $size = sprintf("%.3f", stat($fh)->size / 1024 / 1024);  # Size in MB
} else {
    $fh = *STDIN;
}

while (<$fh>) {
    chomp;
    next if /^\s*$/;

    my @fields = split /\t/;
    if ($fields[0] eq 'V') {
        push @lengths, length($fields[2]);
        $vertices++;
    } elsif ($fields[0] eq 'E') {
        $outdegree{$fields[1]}++;
        $indegree{$fields[2]}++;
        $edges++;
    }
}
close $fh;

# Sort the arrays for median calculation
@lengths = sort {$a <=> $b} @lengths;
my @indegrees = sort {$a <=> $b} values %indegree;
my @outdegrees = sort {$a <=> $b} values %outdegree;

my $min_length = $lengths[0];
my $max_length = $lengths[-1];
my $average_length = sprintf("%.3f", sum(@lengths) / @lengths);
my $median_length = $lengths[@lengths/2];
my $std_dev_length = sprintf("%.3f", sqrt(sum(map {($_ - $average_length) ** 2} @lengths) / @lengths));

my $min_indegree = $indegrees[0];
my $max_indegree = $indegrees[-1];
my $average_indegree = sprintf("%.3f", sum(@indegrees) / @indegrees);
my $median_indegree = $indegrees[@indegrees/2];
my $std_dev_indegree = sprintf("%.3f", sqrt(sum(map {($_ - $average_indegree) ** 2} @indegrees) / @indegrees));

my $min_outdegree = $outdegrees[0];
my $max_outdegree = $outdegrees[-1];
my $average_outdegree = sprintf("%.3f", sum(@outdegrees) / @outdegrees);
my $median_outdegree = $outdegrees[@outdegrees/2];
my $std_dev_outdegree = sprintf("%.3f", sqrt(sum(map {($_ - $average_outdegree) ** 2} @outdegrees) / @outdegrees));

if ($latex_table) {
    print "\\begin{table*}[h]\n";
    print "\\centering\n";
    print "\\small\n";
    print "\\scalebox{0.8}{\n";
    print "\\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|}\n";
    print "\\hline\n";
    print "\\multirow{2}{*}{Input Name} & \\multirow{2}{*}{Size (MB)} & \\multirow{2}{*}{V} & \\multirow{2}{*}{E} & \\multicolumn{5}{c|}{Label Length (bp)} & \\multicolumn{5}{c|}{Indegree} & \\multicolumn{5}{c|}{Outdegree} \\\\ \n";
    print "\\cline{5-19}\n";
    print " & & & & Avg & Med & Std & Min & Max & Avg & Med & Std & Min & Max & Avg & Med & Std & Min & Max \\\\ \n";
    print "\\hline\n";
    print "$base & $size & $vertices & $edges & $average_length & $median_length & $std_dev_length & $min_length & $max_length & $average_indegree & $median_indegree & $std_dev_indegree & $min_indegree & $max_indegree & $average_outdegree & $median_outdegree & $std_dev_outdegree & $min_outdegree & $max_outdegree \\\\ \n";
    print "\\hline\n";
    print "\\end{tabular}}\n";
    print "\\caption{Input Graph Statistics}\n";
    print "\\label{table:input_graph_statistics}\n";
    print "\\end{table*}\n";
} else {
    print "File size (MB): $size\n";
    print "Number of vertices: $vertices\n";
    print "Number of edges: $edges\n";
    print "Min label length: $min_length\n";
    print "Max label length: $max_length\n";
    print "Average label length: $average_length\n";
    print "Median label length: $median_length\n";
    print "Label length standard deviation: $std_dev_length\n";
    print "Min indegree: $min_indegree\n";
    print "Max indegree: $max_indegree\n";
    print "Average indegree: $average_indegree\n";
    print "Median indegree: $median_indegree\n";
    print "Indegree standard deviation: $std_dev_indegree\n";
    print "Min outdegree: $min_outdegree\n";
    print "Max outdegree: $max_outdegree\n";
    print "Average outdegree: $average_outdegree\n";
    print "Median outdegree: $median_outdegree\n";
    print "Outdegree standard deviation: $std_dev_outdegree\n";
}
