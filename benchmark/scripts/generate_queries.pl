#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

# Command-line arguments
my $input_file;
my $length;
my $num_samples;
my $output_file;
my $seed;

GetOptions(
    'input_file=s'  => \$input_file,
    'length=i'      => \$length,
    'num_samples=i' => \$num_samples,
    'output_file=s' => \$output_file,
    'seed=i'        => \$seed,
);

# Set seed if provided
srand($seed) if defined $seed;

my %vertices;
my %adjacency_list;

open(my $fh, '<', $input_file) or die "Can't open $input_file: $!";
while (<$fh>) {
    chomp;
    my ($type, $id, $data) = split(/\t/);
    if ($type eq 'V') {
        $vertices{$id} = $data;
    } elsif ($type eq 'E') {
        my (undef, $from_vertex, $to_vertex) = split(/\t/);
        push @{$adjacency_list{$from_vertex}}, $to_vertex;
    }
}
close $fh;

open(my $ofh, '>', $output_file) or die "Can't open $output_file: $!";
for (1 .. $num_samples) {
    while (1) {
        my $vertex = (keys %vertices)[rand keys %vertices];
        my $position = int(rand(length($vertices{$vertex})));
        my $first_offset = $position;
        my $string = substr($vertices{$vertex}, $position, $length);
        my $needed = $length - length($string);

        while (length($string) < $length) {
            last unless exists $adjacency_list{$vertex};
            $vertex = $adjacency_list{$vertex}[rand @{$adjacency_list{$vertex}}];
            $needed = $length - length($string);
            $string .= substr($vertices{$vertex}, 0, $needed);
        }

        if (length($string) == $length && $string !~ /N/) {
            print $ofh join("\t", $vertex), ",$first_offset\t$needed\n";
            print "$string\n";
            last;
        }
    }
}
print "exit();\n";
close $ofh;