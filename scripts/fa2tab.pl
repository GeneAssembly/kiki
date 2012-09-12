#! /usr/bin/env perl

use strict;

my $usage = "Usage: $0 < contigs.tab > contigs.fa\n\n";

while (<STDIN>) {
    my $hdr = $_;
    my $seq = <STDIN>;
    
    if ($hdr = /^>(\S+)\s+len_(\d+)_cov_([0-9.]+)_stdev_([0-9.]+)_GC_([0-9.]+)_seed_(\S+)/ && $seq) {
        print join("\t", $1, $2, $3, $4, $5, $6, $seq);
    }
}
