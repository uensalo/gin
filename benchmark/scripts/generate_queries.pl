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

    # Start with a random vertex and an initial substring
    my $initial_vertex = $vertex_keys[int(rand(@vertex_keys))];
    my $starting_offset = int(rand(length($vertices{$initial_vertex}) - $length + 1));
    $string = substr($vertices{$initial_vertex}, $starting_offset, $length);

    # If the starting string is already of the desired length and has no 'N', write and continue
    if (length($string) == $length && $string !~ /N/) {
        print $ofh join("\t", $initial_vertex), ",$starting_offset\t$length\n";
        print $ofh "$string\n";
        next;
    }

    push @history_stack, { vertex => $initial_vertex, idx => -1 };

    while (1) {
        my $current_data = $history_stack[-1];
        $current_data->{idx}++;

        if ($current_data->{idx} >= (exists $adjacency_list{$current_data->{vertex}} ? @{$adjacency_list{$current_data->{vertex}}} : 0)) {
            pop @history_stack;
            $string = substr($string, 0, -$length);
            next if @history_stack;
            last;
        }

        my $next_vertex = $adjacency_list{$current_data->{vertex}}[$current_data->{idx}];
        my $segment = substr($vertices{$next_vertex}, 0, $length - length($string));
        $string .= $segment;

        if (length($string) == $length && $string !~ /N/) {
            print $ofh join("\t", $initial_vertex), ",$starting_offset\t$length\n";
            print $ofh "$string\n";
            last;
        } elsif (length($string) < $length) {
            push @history_stack, { vertex => $next_vertex, idx => -1 };
        }
    }
}

print "exit();\n";
close $ofh;
