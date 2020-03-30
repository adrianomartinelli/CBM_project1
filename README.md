# CBM_project1
Project 1: Efficient Search and Read Alignment

To run the mapping for the data_small file use:
main_time.py --genome ./data_small/genome.chr22.5K.fa --read1 ./data_small/output_tiny_30xCov1.fq --read2 ./data_small/output_tiny_30xCov2.fq --outpath result.sam --k 16

Please adjust the paths accordingly to run other genomes or pair-end read files.
The mapping parameters are optimised for the data_small.

NOTE:
The code works perfectly for the small data set. For the large data set the read mapping is rather slow. We build the genome index in about 200 seconds for kmers of length 16 and map 100 reads in 5 sec.
This may be imporved by parameter tuning. For example increase k (from terminal), max_distance (in script) or shit (in script)