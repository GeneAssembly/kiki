#! /usr/bin/env perl

use strict;
use Getopt::Long;

my $usage = "Usage: $0 < contigs.tab > contigs.fa\n\n";

my ($help, $min_contig_length);

GetOptions("h|help"  => \$help,
           "m|min=i" => \$min_contig_length);

$help and die $usage;  

while (<STDIN>) {
    my $hdr = $_;
    my $seq = <STDIN>;
    
    if ($hdr = /^>(\S+)\s+len_(\d+)_cov_([0-9.]+)_stdev_([0-9.]+)_GC_([0-9.]+)_seed_(\S+)/ && $seq) {
        print join("\t", $1, $2, $3, $4, $5, $6, $seq) if $2 >= $min_contig_length;
    }
}
