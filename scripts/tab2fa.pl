#! /usr/bin/env perl

use strict;
use Getopt::Long;

my $usage = "Usage: $0 [ -m min_contig_length ] < contigs.fa > contigs.tab\n\n";

my ($help, $min_contig_length);

GetOptions("h|help"  => \$help,
           "m|min=i" => \$min_contig_length);

$help and die $usage;  

while (<STDIN>) {
    my ($id, $len, $cov, $stdev, $gc, $seed, $seq) = split /\t/;
    next if $len < $min_contig_length;
    print join("_", ">$id len", $len, 'cov', $cov, 'stdev', $stdev, 'GC', $gc, 'seed', $seed) . "\n";
    print $seq;
}
