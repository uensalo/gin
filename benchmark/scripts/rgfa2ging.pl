#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

# Compute reverse complement of a sequence
sub reverse_complement {
    my $seq = shift;
    $seq = reverse($seq);
    $seq =~ tr/ACGTacgt/TGCAtgca/;
    return $seq;
}

my $gfa_file;
my $ging_file;

# Parsing command-line options
GetOptions(
    'i=s' => \$gfa_file,
    'o=s' => \$ging_file
) or die "Usage: $0 -i <input.rGFA> -o <output.GING>\n";

# Check if input and output files are provided
if (not defined $gfa_file or not defined $ging_file) {
    die "Usage: $0 -i <input.rGFA> -o <output.GING>\n";
}

# Data structures to hold sequences and links
my %sequences;
my @links;

# Open the rGFA file and read it line by line
open my $gf, '<', $gfa_file or die "Could not open rGFA file: $!\n";
while (my $line = <$gf>) {
    chomp $line;
    my @parts = split "\t", $line;
    if ($parts[0] eq 'S') {
        # Store both the sequence and its reverse complement
        $sequences{$parts[1]} = $parts[2];
        $sequences{"${parts[1]}'"} = reverse_complement($parts[2]);
    } elsif ($parts[0] eq 'L' || $parts[0] eq 'E') {
        my ($src, $src_dir, $dest, $dest_dir);
        if ($parts[0] eq 'L') {
            ($src, $src_dir, $dest, $dest_dir) = ($parts[1], $parts[2], $parts[3], $parts[4]);
        } else {
            ($src, $src_dir, $dest, $dest_dir) = ($parts[2], $parts[4], $parts[3], $parts[5]);
        }

        $src .= "'" if $src_dir eq '-';
        $dest .= "'" if $dest_dir eq '-';
        push @links, [$src, $dest];
    }
}
close $gf;

# Create a mapping for sequences to unique integers
my %sequence_map;
my $counter = 0;
for my $seq_id (sort keys %sequences) {  # Sorting to maintain consistency
    $sequence_map{$seq_id} = $counter++;
}

# Write the GING file using the remapping
open my $fd, '>', $ging_file or die "Could not open GING file: $!\n";
for my $seq_id (sort keys %sequences) {  # Sorting to maintain order
    my $new_id = $sequence_map{$seq_id};
    print $fd "V\t$new_id\t$sequences{$seq_id}\n";
}
for my $link (@links) {
    my $src_new_id = $sequence_map{$link->[0]};
    my $dest_new_id = $sequence_map{$link->[1]};
    print $fd "E\t$src_new_id\t$dest_new_id\n";
}
close $fd;
