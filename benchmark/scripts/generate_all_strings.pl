#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

my $alphabet;
my $length;

GetOptions(
    'a=s' => \$alphabet,
    'l=i' => \$length
) or die("Error in command line arguments\n");

if (!defined $alphabet || !defined $length) {
    die("Usage: $0 -a <alphabet> -l <length>\n");
}

generate_strings("", $length);
print "exit();\n";

sub generate_strings {
    my ($prefix, $len) = @_;

    if ($len == 0) {
        print "$prefix\n";
        return;
    }

    for my $char (split //, $alphabet) {
        generate_strings($prefix . $char, $len - 1);
    }
}
