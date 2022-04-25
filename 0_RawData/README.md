# Raw Sequences

# Note:
Raw FASTQ files are too large and are stored on NCBI/EBI. This page is mostly the metadata.
SRA:SRP021107
BioProject:PRJNA196801

## Files
sample_metadata.txt : Metadata for both the bulk and viral raw reads.

## Directories

## Scripts

## Description

One healthy individual was sampled at 16 time points. For 3 time points bulk metagenomics was sequenced (bacteria metagenomics) and for all 16 time points viral metagenomics was completed.
For 8 time points two seperate viral metagenomic read runs were done. So there are 24 viral sequence runs.

For more information on the metadata and how samples were collected see orginal manuscript.


### Bacteria (Bulk Metagenomics)
See Fig S2 from Minot 2013 paper for days.
| SRA	  | SampleID	      | Day	| Week |		
| ------- | ----------------- | ------- | ---- |
|SRR828645| BLS020308_1014-01 | 182 	|  1   |
|SRR828661| BLS020309_1014-02 | 852	|  2   |  
|SRR828660| BLS020310_1014-03 | 882	|  3   |

### Viral Metagenomics
See Fig 1A & Table S1 (note there is a typo for day 851-2)
Unfortunately Table S1 has the '# of reads' which appears to be after QC, as they do no match raw reads.
I followed the labels from Sample_ID where 'd00-1' where the first digit is the week. Day zero is 00, and week 1 is 01-5.
The second digit is the day of the week.

| SRA       | SampleID        | Week | Day | Replicate |
|-----------|-----------------|------|-----|-----------|
| SRR829033 | BLS020813_d00-1 | 0    | 0   | 1         |
| SRR829867 | BLS020814_d00-2 | 0    | 0   | 2         |
| SRR935337 | BLS020815_d01-1 | 1    | 180 | 1         |
| SRR935339 | BLS020816_d02-1 | 1    | 181 | 1         |
| SRR935340 | BLS020817_d03-1 | 1    | 182 | 1         |
| SRR935341 | BLS020818_d03-2 | 1    | 182 | 2         |
| SRR935342 | BLS020819_d04-1 | 1    | 183 | 1         |
| SRR829034 | BLS020820_d05-1 | 1    | 184 | 1         |
| SRR935364 | BLS020821_d11-1 | 2    | 851 | 1         |
| SRR935366 | BLS020822_d11-2 | 2    | 851 | 2         |
| SRR935367 | BLS020823_d12-1 | 2    | 852 | 1         |
| SRR935361 | BLS026407_d12-2 | 2    | 852 | 2         |
| SRR935348 | BLS020824_d13-1 | 2    | 853 | 1         |
| SRR935349 | BLS020825_d13-2 | 2    | 853 | 2         |
| SRR935350 | BLS020826_d14-1 | 2    | 854 | 1         |
| SRR935351 | BLS020827_d15-1 | 2    | 855 | 1         |
| SRR935352 | BLS020828_d21-1 | 3    | 879 | 1         |
| SRR935354 | BLS020829_d21-2 | 3    | 879 | 2         |
| SRR935356 | BLS020830_d22-1 | 3    | 880 | 1         |
| SRR935363 | BLS020831_d22-2 | 3    | 880 | 2         |
| SRR935357 | BLS020832_d23-1 | 3    | 881 | 1         |
| SRR935358 | BLS020833_d23-2 | 3    | 881 | 2         |
| SRR935359 | BLS020834_d24-1 | 3    | 882 | 1         |
| SRR935360 | BLS020835_d25-1 | 3    | 883 | 1         |
