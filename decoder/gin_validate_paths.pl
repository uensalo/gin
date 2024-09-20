##
# gin: FM-Index-like graph indexing algorithm toolkit.
# Copyright (C) 2024, Unsal Ozturk
#
# gin_validate_paths.pl is part of gin
#
# gin is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# gin is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

sub parse_graph {
    my ($graph_file) = @_;
    my %labels;
    my %graph;

    open my $fh, '<', $graph_file or die "Cannot open $graph_file: $!";
    while (<$fh>) {
        chomp;
        next unless $_;
        my @fields = split /\t/;
        if ($fields[0] eq 'V') {
            my $vertex_id = $fields[1];
            my $label = $fields[2];
            $labels{$vertex_id} = $label;
        } elsif ($fields[0] eq 'E') {
            my $from_vertex = $fields[1];
            my $to_vertex = $fields[2];
            push @{$graph{$from_vertex}}, $to_vertex;
        }
    }
    close $fh;
    return (\%labels, \%graph);
}

sub parse_output {
    my ($output_file) = @_;
    my %patterns;
    my $current_pattern;

    open my $fh, '<', $output_file or die "Cannot open $output_file: $!";
    while (<$fh>) {
        chomp;
        if (/^(.*):$/) {
            # This is a pattern
            $current_pattern = $1;
            $patterns{$current_pattern} = [];
        } elsif ($_ eq '-' or $_ eq "\t-") {
            # No paths found for this pattern
            $patterns{$current_pattern} = [];
        } elsif (/^\t\((\d+),(\d+)\);(.*)$/) {
            # This is a path
            my $o1 = $1;
            my $oN = $2;
            my $path_str = $3;
            my @vertices = split /:/, $path_str;
            push @{$patterns{$current_pattern}}, {
                'o1'      => $o1,
                'oN'      => $oN,
                'vertices' => \@vertices,
            };
        }
    }
    close $fh;
    return \%patterns;
}

sub validate_paths {
    my ($labels_ref, $patterns_ref, $verbose) = @_;
    my %labels = %$labels_ref;
    my %patterns = %$patterns_ref;

    my $total_paths = 0;
    my $max_path_length = 0;
    my $max_paths_query = '';
    my $max_paths_count = 0;

    my %invalid_paths;

    foreach my $pattern (keys %patterns) {
        my $paths = $patterns{$pattern};
        my $num_paths = scalar @$paths;
        $total_paths += $num_paths;

        # Update the query with the most paths
        if ($num_paths > $max_paths_count) {
            $max_paths_count = $num_paths;
            $max_paths_query = $pattern;
        }

        foreach my $path_info (@$paths) {
            my $o1 = $path_info->{'o1'};
            my $oN = $path_info->{'oN'};
            my @vertices = @{$path_info->{'vertices'}};
            my $path_length = scalar @vertices;

            # Update maximum path length
            if ($path_length > $max_path_length) {
                $max_path_length = $path_length;
            }

            my $reconstructed = '';
            if ($path_length == 1) {
                # Single-vertex path
                my $vertex = $vertices[0];
                my $label = $labels{$vertex} // '';
                my $length = $oN - $o1;
                my $label_part = substr($label, $o1, $length);
                $reconstructed = $label_part;
            } else {
                for (my $idx = 0; $idx < $path_length; $idx++) {
                    my $vertex = $vertices[$idx];
                    my $label = $labels{$vertex} // '';
                    my $label_part = '';
                    if ($idx == 0) {
                        # First vertex
                        $label_part = substr($label, $o1);
                    } elsif ($idx == $path_length - 1) {
                        # Last vertex
                        $label_part = substr($label, 0, $oN);
                    } else {
                        # Intermediate vertices
                        $label_part = $label;
                    }
                    $reconstructed .= $label_part;
                }
            }
            if ($reconstructed ne $pattern) {
                push @{$invalid_paths{$pattern}}, {
                    'vertices' => \@vertices,
                    'expected' => $pattern,
                    'got'      => $reconstructed,
                };
            }
        }
    }

    # Print invalid paths
    foreach my $pattern (sort keys %invalid_paths) {
        print "Pattern '$pattern':\n";
        foreach my $path_info (@{$invalid_paths{$pattern}}) {
            my @vertices = @{$path_info->{'vertices'}};
            my $path_str = join '->', @vertices;
            print "\tInvalid path: $path_str\n";
            print "\tExpected: $path_info->{'expected'}\n";
            print "\tGot:      $path_info->{'got'}\n";
        }
    }

    # Print statistics if verbose is enabled
    if ($verbose) {
        print "\nStatistics:\n";
        print "Total number of paths: $total_paths\n";
        print "Maximum path length (in vertices): $max_path_length\n";
        if ($max_paths_query ne '') {
            print "Query with the most paths: '$max_paths_query' ($max_paths_count paths)\n";
        } else {
            print "No paths found for any query.\n";
        }
    }
}

sub main {
    my ($graph_file, $input_file, $verbose);
    GetOptions(
        'g|graph=s'  => \$graph_file,
        'i|input=s'  => \$input_file,
        'v|verbose'  => \$verbose,
    ) or die "Error in command line arguments.\n";

    unless ($graph_file && $input_file) {
        die "Usage: $0 -g <graph_file> -i <input_file> [-v]\n";
    }

    my ($labels_ref, $graph_ref) = parse_graph($graph_file);
    my $patterns_ref = parse_output($input_file);
    validate_paths($labels_ref, $patterns_ref, $verbose);
}

main();
