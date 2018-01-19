#!/bin/bash

~/software/kallisto_linux-v0.43.1/kallisto quant -i ./ref/transcripts.idx -o CRI1 --threads=2 -b 100 ./raw_reads/raw_reads/2528-01-2-1_S24_R1_001.fastq.gz ./raw_reads/raw_reads/2528-01-2-1_S24_R2_001.fastq.gz
~/software/kallisto_linux-v0.43.1/kallisto quant -i ./ref/transcripts.idx -o CRI2 --threads=2 -b 100 ./raw_reads/raw_reads/2528-02-2-1_S25_R1_001.fastq.gz ./raw_reads/raw_reads/2528-02-2-1_S25_R2_001.fastq.gz
~/software/kallisto_linux-v0.43.1/kallisto quant -i ./ref/transcripts.idx -o CRI3 --threads=2 -b 100 ./raw_reads/raw_reads/2528-03-2-1_S26_R1_001.fastq.gz ./raw_reads/raw_reads/2528-03-2-1_S26_R2_001.fastq.gz
~/software/kallisto_linux-v0.43.1/kallisto quant -i ./ref/transcripts.idx -o CRI4 --threads=2 -b 100 ./raw_reads/raw_reads/2528-04-2-1_S27_R1_001.fastq.gz ./raw_reads/raw_reads/2528-04-2-1_S27_R2_001.fastq.gz
~/software/kallisto_linux-v0.43.1/kallisto quant -i ./ref/transcripts.idx -o KO1 --threads=2 -b 100 ./raw_reads/raw_reads/2528-05-2-1_S28_R1_001.fastq.gz ./raw_reads/raw_reads/2528-05-2-1_S28_R2_001.fastq.gz
~/software/kallisto_linux-v0.43.1/kallisto quant -i ./ref/transcripts.idx -o KO2 --threads=2 -b 100 ./raw_reads/raw_reads/2528-06-2-1_S29_R1_001.fastq.gz ./raw_reads/raw_reads/2528-06-2-1_S29_R2_001.fastq.gz
~/software/kallisto_linux-v0.43.1/kallisto quant -i ./ref/transcripts.idx -o KO3 --threads=2 -b 100 ./raw_reads/raw_reads/2528-07-2-1_S30_R1_001.fastq.gz ./raw_reads/raw_reads/2528-07-2-1_S30_R2_001.fastq.gz
~/software/kallisto_linux-v0.43.1/kallisto quant -i ./ref/transcripts.idx -o KO4 --threads=2 -b 100 ./raw_reads/raw_reads/2528-08-2-1_S31_R1_001.fastq.gz ./raw_reads/raw_reads/2528-08-2-1_S31_R2_001.fastq.gz
