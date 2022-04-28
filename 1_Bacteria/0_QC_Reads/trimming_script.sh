#!/bin/bash

#SBATCH --time=10:00:00
#SBATCH --mem=40G
#SBATCH --cpus-per-task=16

module load trimmomatic/0.36

raw_dic=../../0_RawData/0_Bacteria_Reads
#Note make sure reads are unzipped

java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE $raw_dic/SRR828645_1.fastq \
 $raw_dic/SRR828645_2.fastq output_SRR828645_1_paired.fastq \
output_SRR828645_1_unpaired.fastq output_SRR828645_2_paired.fastq \
output_SRR828645_2_unpaired.fastq \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:8:keepBothReads LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE  $raw_dic/SRR828660_1.fastq \
 $raw_dic/SRR828660_2.fastq output_SRR828660_1_paired.fastq \
output_SRR828660_1_unpaired.fastq output_SRR828660_2_paired.fastq \
output_SRR828660_2_unpaired.fastq \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:8:keepBothReads LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE  $raw_dic/SRR828661_1.fastq \
 $raw_dic/SRR828661_2.fastq output_SRR828661_1_paired.fastq \
output_SRR828661_1_unpaired.fastq output_SRR828661_2_paired.fastq \
output_SRR828661_2_unpaired.fastq \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:8:keepBothReads LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

