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
    my @history_stack;
    my $string = '';

    # Start with a random vertex
    my $initial_vertex = $vertex_keys[int(rand(@vertex_keys))];
    push @history_stack, { vertex => $initial_vertex, idx => -1 };

    while (1) {
        my $current_data = $history_stack[-1]; # Look at the top of the stack
        $current_data->{idx}++; # Move to the next adjacent vertex

        # If we've tried all adjacent vertices, backtrack
        if ($current_data->{idx} >= (exists $adjacency_list{$current_data->{vertex}} ? @{$adjacency_list{$current_data->{vertex}}} : 0)) {
            pop @history_stack;
            $string = substr($string, 0, -$length); # Remove the last segment from the string
            next if @history_stack; # Continue with the next vertex in the stack
            last; # Exit if the stack is empty
        }

        # Extract a segment from the current vertex data
        my $next_vertex = $adjacency_list{$current_data->{vertex}}[$current_data->{idx}];
        my $segment = substr($vertices{$next_vertex}, 0, $length - length($string));
        $string .= $segment;

        # Check the length and content of the generated string
        if (length($string) == $length && $string !~ /N/) {
            print $ofh join("\t", $initial_vertex), "\t$string\n";
            last;
        } elsif (length($string) < $length) {
            # If the string is shorter than required, continue building it
            push @history_stack, { vertex => $next_vertex, idx => -1 };
        }
    }
}

print "exit();\n";
close $ofh;
