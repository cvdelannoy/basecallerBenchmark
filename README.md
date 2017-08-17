# basecallerBenchmark
Code for a MinION sequencer benchmark

Starting from a directoyr of fast5-files, assess basecalling accuracy, assemble the reads using minimap/miniasm and assess the quality of the assembly. 

Prerequisites:
- Samtools (v1.3.1!)
- Jellyfish
- Minimap
- Miniasm
- BWA
- Quast

Download the repository, cd to the basecallerBenchmark folder and run using bash basecallerBenchmark.sh. A help option (-h) is available. Note that running the entire routine on many reads may take a lot of memory! I plan to find a solution for this in the future.
