#!/usr/bin/perl
use strict;
use warnings;

sub parse_file {
    my ($filename) = @_;
    open my $fh, '<', $filename or die "Cannot open $filename: $!";

    my %data;
    my $current_string;

    while (<$fh>) {
        chomp;
        if (/^([^:]+):$/) {
            $current_string = $1;
        } elsif (my ($offsets, $path) = /\t\((\d+,\d+)\);(.+)$/) {
            push @{$data{$current_string}}, { offsets => $offsets, path => [split /:/, $path] };
        }
    }
    close $fh;

    return \%data;
}

my ($file1, $file2) = @ARGV;
my $data1 = parse_file($file1);
my $data2 = parse_file($file2);

foreach my $string (keys %$data1) {
    unless (exists $data2->{$string}) {
        print "String '$string' exists in $file1 but not in $file2\n";
        next;
    }

    my @paths1 = sort { $a->{offsets} cmp $b->{offsets} } @{$data1->{$string}};
    my @paths2 = sort { $a->{offsets} cmp $b->{offsets} } @{$data2->{$string}};

    if (@paths1 != @paths2) {
        print "Different number of paths for string '$string' between files\n";
        next;
    }

    for my $i (0 .. $#paths1) {
        unless ($paths1[$i]->{offsets} eq $paths2[$i]->{offsets} && join(":", @{$paths1[$i]->{path}}) eq join(":", @{$paths2[$i]->{path}})) {
            print "Paths for string '$string' differ between files\n";
            last;
        }
    }
}

foreach my $string (keys %$data2) {
    unless (exists $data1->{$string}) {
        print "String '$string' exists in $file2 but not in $file1\n";
    }
}
