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
    my $initial_vertex;
    while (1) {
        $initial_vertex = $vertex_keys[int(rand(@vertex_keys))];
        last if exists $adjacency_list{$initial_vertex};
    }
    push @history_stack, { vertex => $initial_vertex, idx => -1 };

    while (1) {
        my $current_data = $history_stack[-1]; # Look at the last element in the stack
        $current_data->{idx}++; # Move to the next adjacent vertex

        # If the current vertex has no adjacency or we've tried all adjacent vertices, backtrack
        if (!exists $adjacency_list{$current_data->{vertex}} || $current_data->{idx} >= @{$adjacency_list{$current_data->{vertex}}}) {
            pop @history_stack;
            $string = substr($string, 0, -$length); # Remove the last segment from the string
            next;
        }

        # Ensure that the current vertex exists in the adjacency list and that the index is within bounds
        if (exists $adjacency_list{$current_data->{vertex}} && defined $adjacency_list{$current_data->{vertex}}[$current_data->{idx}]) {
            my $next_vertex = $adjacency_list{$current_data->{vertex}}[$current_data->{idx}];
            my $segment = substr($vertices{$next_vertex}, 0, $length - length($string));
            $string .= $segment;

            # If the generated string is of the required length and doesn't contain 'N', finish this sample
            if (length($string) == $length && $string !~ /N/) {
                my $starting_vertex = $history_stack[0]{vertex};
                my $starting_offset = int(rand(length($vertices{$starting_vertex})));
                print $ofh join("\t", $starting_vertex), ",$starting_offset\t$length\n";
                print "$string\n";
                last;
            } elsif (length($string) < $length) {
                push @history_stack, { vertex => $next_vertex, idx => -1 }; # Push the next vertex onto the stack and continue
            }
        }
    }
}
print "exit();\n";
close $ofh;
