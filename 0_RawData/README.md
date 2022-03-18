# Raw Sequences

# Note:
Due to file storage issues (files being too large) I've implented protocols to move them.
SRA:SRP021107
BioProject:PRJNA196801

## Files

## Directories

* 0_Bacteria : Bulk metagenomics for which we assembled bacterial genomes
* 1_Viral : VLPs were seperated and sequenced.

## Scripts

## Description

One healthy individual was sampled at 16 time points. For 3 time points bulk metagenomics was sequenced (bacteria metagenomics) and for all 16 time points viral metagenomics was completed.
For 8 time points two seperate viral metagenomic read runs were done. So there are 24 viral sequence runs.

For more information on the metadata and how samples were collected see orginal manuscript.

### Bacteria (Bulk Metagenomics)

SRR828645 - BLS020308_1014-01 - Day 182 - Week 1
SRR828661 - BLS020309_1014-02 - Day 852 - Week 2
SRR828660 - BLS020310_1014-03 - Day 882 - Week 3

### Viral Metagenomics
SRA		Sample_ID	Day	Replicate
SRR829867	BLS020814_d00-2	0	2
SRR935337	BLS020815_d01-1	180	1
SRR935339	BLS020816_d02-1	181	1
SRR935340	BLS020817_d03-1	182	1
SRR935341	BLS020818_d03-2	182	2
SRR935342	BLS020819_d04-1	183	1
SRR829034	BLS020820_d05-1	184	1
SRR935364	BLS020821_d11-1	851	1
SRR935366	BLS020822_d11-2	851	2
SRR935367	BLS020823_d12-1	852	1
SRR935361	BLS026407_d12-2	852	2
SRR935348	BLS020824_d13-1	853	1
SRR935349	BLS020825_d13-2	853	2
SRR935350	BLS020826_d14-1	854	1
SRR935351	BLS020827_d15-1	855	1
SRR935352	BLS020828_d21-1	879	1
SRR935354	BLS020829_d21-2	879	2
SRR935356	BLS020830_d22-1	880	1
SRR935363	BLS020831_d22-2	880	2
SRR935357	BLS020832_d23-1	881	1
SRR935358	BLS020833_d23-2	881	1
SRR935359	BLS020834_d24-1	882	1
SRR935360	BLS020835_d25-1	883	1

Note: Table S1 does not match up with the data they uploaded or from Fig.1 of the Minot et al dataset.
For example day 851 should be in duplicate based on Fig 1. and sequences uploaded on SRA but is only one in Table S1. So I do not trust Table S1
The pattern for d##_# works out to be d<week><day of week>-<replicate>

