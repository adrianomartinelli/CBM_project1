# Computational Biomedicine Project 1: Efficient Search and Read Alignment

## Installation

- First clone and go to the repository:

```bash
git clone https://github.com/adrianomartinelli/CBM_project1.git
cd CBM_project1
```

- Then download the small data file.

```bash
wget -r -np -nH --cut-dirs=3 -R "index.html*" https://public.bmi.inf.ethz.ch/teaching/cbm1_2019/project1/data_small/ && cd data_small && gunzip *.gz && cd ..
```

- Optionally, download the big data file.

```bash
wget -r -np -nH --cut-dirs=3 -R "index.html*" https://public.bmi.inf.ethz.ch/teaching/cbm1_2019/project1/data/ && cd data && gunzip *.gz && cd ..
```

## Running the code

- To run the mapping for the small data file, use the following:
main_time.py --genome ./data_small/genome.chr22.5K.fa --read1 ./data_small/output_tiny_30xCov1.fq --read2 ./data_small/output_tiny_30xCov2.fq --outpath result.sam --k 16

- Please adjust the paths accordingly to run other genomes or pair-end reads for big data file. The mapping parameters are optimized for the small data file.

- NOTE:
The code works perfectly for the small dataset. For the large dataset, the read mapping is rather slow.
We build the genome index in about 200 seconds for k-mers of length 16 and map 100 reads in 5 sec.
This may be improved by parameter tuning. For example increase k (from terminal) or max_distance (in script).
