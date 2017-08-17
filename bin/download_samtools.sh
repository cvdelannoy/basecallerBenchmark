#!/bin/bash

# Download and unpack samtools
wget https://sourceforge.net/projects/samtools/files/samtools/1.3.1/samtools-1.3.1.tar.bz2
tar xf samtools-1.3.1.tar.bz2
cd samtools-1.3.1
./configure --enable-plugins --enable-libcurl --with-plugin-path=$PWD/htslib-1.3.1
make all all-htslib
make install install-htslib
