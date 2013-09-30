#! /usr/bin/env perl

use strict;
use Getopt::Long;

my $usage = "Usage: $0 [ -m min_contig_length -c min_coverage ] < contigs.fa > contigs.tab\n\n";

my ($help, $min_contig_length, $min_coverage);

GetOptions("h|help"  => \$help,
           "m|min=i" => \$min_contig_length,
           "c|cov=f" => \$min_coverage );

$help and die $usage;  

while (<STDIN>) {
    my ($id, $len, $cov, $stdev, $gc, $seed, $seq) = split /\t/;
    next if $len < $min_contig_length;
    next if $cov < $min_coverage;
    print join("_", ">$id len", $len, 'cov', $cov, 'stdev', $stdev, 'GC', $gc, 'seed', $seed) . "\n";
    print $seq;
}
