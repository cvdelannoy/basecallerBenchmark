#!/bin/bash

# basecallerBenchmark Setup script

echo "Setting up basecaller benchmark for the first time"



# 1. Check availability of binaries, attempt download if necessary
echo "Checking availability required binaries..."
installAll=false
mkdir -p bin

binChecks=( "bwa" "samtools" )
# binChecks=( "bwa" "samtools" "minimap" "miniasm" "quast" "jellyfish" )
for cb in ${binChecks[@]}; do
if hash $cb; then
        echo "$cb detected...";
        else
                if [ -v  $installAll ]; then
                        echo "$cb not detected. Would you like to download $cb binaries? [y/n/a]"
                        read installChoice
                fi
                if [ $installAll ] || [ "$installChoice"="y" ] || [ "$installChoice"="a" ]
                        then echo "Downloading $cb..."
                        bash bin/download_$cb.sh
                fi
                if [ $installChoice -eq "a" ]
                        then installAll=true
                fi
        fi
done
