Kiki
====

De novo metagenomic assembler


Installation
------------

1. Clone Kiki repo

    git clone https://github.com/GeneAssembly/kiki.git 

    (or 'git clone git://github.com/GeneAssembly/kiki.git' for read-only access)

2. Make sure CMake is installed

3. Build and test

    cd kiki<br/>
    mkdir bin<br/>
    cd bin<br/>
    cmake ..<br/>
    make ki<br/>
    ./ki -i ../test/short.fa -o test<br/>
    cat test.contig* | ../scripts/tab2fa.pl >test.fa <br/>

4. Scripts for converting 

   Example 1: sort contigs by length and create a fasta file containing contigs longer than 1kb. <br/>
     cat test.contig* | sort -k2nr | ../scripts/tab2fa.pl -m 1000 >test.fa <br/><br/>

   Example 2: convert the fasta contig file to the tab-delimited format [ id, len, cov, stdev, GC, seed, seq ] <br/>
     cat test.contig* | sort -k2nr | ../scripts/tab2fa.pl -m 1000 >test.fa <br/>


Usage
-----

    Usage: ki (-i single.file | -I file.list) [options] 
         -k int        kmer length, also used as minimum overlap (D = 25) 
         -o outname    output prefix 
         -u float|int  fraction of nodes used as users or the number of user nodes 
         -v [level]    verbosity level (D = 3)  
         -persist [t]  checkpointing interval (eg, 30m, 6h, D = 1h) 
         -dot file     generate graphviz output 
         -h            show this usage information 



License
-------

Copyright (c) 2011-2012, University of Chicago All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.