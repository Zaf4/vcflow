# Variant Calling flow with Python

## Requirements

bwa   
samtools   
gatk==4.6.1.0
- Java 17   

python=>3.10   
- requests  
- tqdm   

### With Conda

```shell
conda env create -f env.yaml
```

## Flow

1. unzip the zip file
2. get the reference genome hg38 from UCSC
3. Index the reference genome (`bwa index`, `samtools faidx`)
3. align the reads to the reference genome (`bwa mem`)
4. post-process the alignment
    - convert sam to bam (`samtools view -bhS`)
    - add read groups (`samtools addreplacerg`)
    - sort the bam file (`samtools sort`)
    - generate stats (`samtools stats`)
    - index the bam file (`samtools index`)
5. call variants (GATK)
    - create dict file (`gatk CreateSequenceDictionary`)
    - call variants (`gatk HaplotypeCaller`)
6. hard filter variants
    - filter vairants with (`VariantFiltration`)


Note that:
`CNNScoreVariants` could not be used because it is no longer included in GATK as of version 4.6.1.0.
Suggested alternative `NVScoreVariants` was not compatible with my computer.

Therefore, downstream `FilterVariantTranches` could not be used.

Instead, I used `VariantFiltration` to hard filter variants.

## Results

Resulting Files and Directories
```shell
ref
├── hg38.dict
├── hg38.fa
├── hg38.fa.amb
├── hg38.fa.ann
├── hg38.fa.bwt
├── hg38.fa.fai
├── hg38.fa.gz
├── hg38.fa.pac
└── hg38.fa.sa
output
├── POOL-37.bam
├── POOL-37.filtered.vcf
├── POOL-37.filtered.vcf.idx
├── POOL-37.rg.bam
├── POOL-37.rg.sorted.bam
├── POOL-37.rg.sorted.bam.bai
├── POOL-37.rg.sorted.bam.stats
├── POOL-37.sam
├── POOL-37.vcf
└── POOL-37.vcf.idx
Task_Files
├── POOL-37.bed
├── POOL-37_S182_L004_R1_001.fastq.gz
└── POOL-37_S182_L004_R2_001.fastq.gz

3 directories, 22 files
```