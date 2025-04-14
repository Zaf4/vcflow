# Variant Calling flow with Python

## Requirements

bwa
samtools
gatk
    - Java 17

python=>3.10
    - requests
    - tqdm

## Flow

1. unzip the zip file
2. get the reference genome hg38 from UCSC
3. Index the reference genome (bwa index)
3. align the reads to the reference genome (bwa mem)
4. post-process the alignment
5. call variants (GATK)