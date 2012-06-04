Kiki
====

De novo metagenomic assembler


Installation
====

1. Clone Kiki repo
git https://github.com/GeneAssembly/kiki.git
(or git git://github.com/GeneAssembly/kiki.git for read-only access)

2. Make sure CMake is installed

3. Build and test
cd kiki
mkdir bin
cd bin
cmake ..
make ki
./ki -i ../test/short.fa -o test

Usage
====

Usage: ki (-i single.file | -I file.list) [options] 
         -k int        kmer length, also used as minimum overlap (D = 25) 
         -o outname    output prefix 
         -u float|int  fraction of nodes used as users or the number of user nodes 
         -v [level]    verbosity level (D = 3)  
         -persist [t]  checkpointing interval (eg, 30m, 6h, D = 1h) 
         -dot file     generate graphviz output 
         -h            show this usage information 


License
====

Copyright (c) 2011-2012, University of Chicago All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.