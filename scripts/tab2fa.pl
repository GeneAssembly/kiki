#! /usr/bin/env perl

use strict;

my $usage = "Usage: $0 < contigs.fa > contigs.tab\n\n";

while (<STDIN>) {
    my ($id, $len, $cov, $stdev, $gc, $seed, $seq) = split /\t/;
    print join("_", ">$id len", $len, 'cov', $cov, 'stdev', $stdev, 'GC', $gc, 'seed', $seed) . "\n";
    print $seq;
}
