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

# Pre-calculate vertex keys
my @vertex_keys = keys %vertices;

open(my $ofh, '>', $output_file) or die "Can't open $output_file: $!";
for (1 .. $num_samples) {
    my @stack;
    my $string = '';
    my $vertex = $vertex_keys[int(rand(@vertex_keys))];
    my $position = int(rand(length($vertices{$vertex})));

    while (length($string) < $length) {
        my $remaining = $length - length($string);
        my $segment = substr($vertices{$vertex}, $position, $remaining);
        $string .= $segment;
        push @stack, $vertex;

        if (length($string) == $length) {
            last if $string =~ /N/;
            # Pop the last vertex and backtrack if the string contains 'N'
            $string = substr($string, 0, length($string) - length($segment));
            pop @stack;
        }

        # Choose the next vertex from the adjacency list
        if (exists $adjacency_list{$vertex}) {
            $vertex = $adjacency_list{$vertex}[rand @{$adjacency_list{$vertex}}];
            $position = 0;
        } else {
            # If no next vertex, backtrack
            $string = substr($string, 0, length($string) - length($segment));
            pop @stack;
            last unless @stack;
            $vertex = $stack[-1];
        }
    }

    if (length($string) == $length) {
        my $start_vertex = $stack[0];
        print $ofh join("\t", $start_vertex, $position, $length), "\n";
        print $ofh "$string\n";
    }
}

print "exit();\n";
close $ofh;
