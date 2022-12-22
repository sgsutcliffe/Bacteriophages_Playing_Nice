#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 15:54:15 2020

@author: steven
"""

from Bio import SeqIO

import csv

#loading in the file with three columns 1 bacterial bin contig name 2 base site 3 coverage output from mpileup
path_mpileup="Data_For_Scripts/three_columns_mpileup_DAScheckbins"
tsv_file = open(path_mpileup)

read_tsv = csv.reader(tsv_file, delimiter="\t")

#turns tsv into a list of lists, each list is a row
mpileup = []
for row in read_tsv:
    
    mpileup.append(row)

contig_list = []


#makes a list of all contigs in tsv
for n in mpileup:
    if n[0] not in contig_list:
        contig_list.append(n[0])

print('made contig list')

#iterate through all the contigs and see if there is a prophage region and add them to the dictionary

prophage_regions = []

#contig by contig
for contig in contig_list:
    
    
    #this takes all lines that belong to this specifc contig and makes a contig specific pileup list
    mpileup_contig = []
    
    #collect all base-pairs or rows for a contig
    for i in mpileup:
        
        if i[0] == contig:
            
            mpileup_contig.append(i)
    
    #For regions where there were no hits mpileup skips these lines, so I added a zero coverage value to fill in the 'gaps'
    for k in range(1,len(mpileup_contig)):
       
        if (int(mpileup_contig[k-1][1]) + 1) != (int(mpileup_contig[k][1])):
           
            mpileup_contig.insert(int(k), [contig, (int(mpileup_contig[k-1][1]) + 1), 0])
    
    #Now we calculate coverage for a sliding window I can define
    sliding_window = 1000
    
    #Set a min coverage to alert me 
    min_coverage = 10
    contig_coverage = 0
    
    #Go over contig measure number of reads that map to each base
    for d in range(0,(len(mpileup_contig) - sliding_window)):
        contig_reads = 0
        for q in range(d,(d + sliding_window)):
            contig_reads += int(mpileup_contig[q][2])
        
        #calculate coverage
        contig_coverage = contig_reads/sliding_window
     
        if contig_coverage >= min_coverage:
             # print(contig)
             # print(str(d) + ":" + str(d + sliding_window))
             # print(contig_coverage)
             # print(d)
             prophage_regions.append([contig, d, (d+sliding_window)])

for b in range(0,len(prophage_regions)):
    outF = open("prophage_regions_step1.txt", "a")
    outF.writelines(str(prophage_regions[b][0]) + ' ' + str(prophage_regions[b][1]) + ' ' + str(prophage_regions[b][2]) + '\n')
    outF.close()


print('found prophage regions')

for n in range(0, len(prophage_regions)):
    
    
    while (len(prophage_regions) > (n + 1)) and (prophage_regions[n][0] == prophage_regions[n+1][0]) and (int(prophage_regions[n][2]) > int(prophage_regions[n+1][1])):
        prophage_regions[n][2] = prophage_regions[n+1][2]
        prophage_regions.pop(n+1) 


for op in range(0,len(prophage_regions)):
    outF = open("prophage_regions_step1.txt", "a")
    outF.writelines(str(prophage_regions[op][0]) + ' ' + str(prophage_regions[op][1]) + ' ' + str(prophage_regions[op][2]) + '\n')
    outF.close()

print('Regions')

#Here I load in the concatenated fasta file of all the bacterial bins involved in the study
bacterial_bin = {}
for seq_record1 in SeqIO.parse("Data_For_Scripts/Final_Collection_DAS_check_MAGs.fa"", "fasta"):
        bacterial_bin[seq_record1.id] = seq_record1.seq

for b1 in range(0, len(prophage_regions)):
    if prophage_regions[b1][0] in bacterial_bin:
                   
        outF = open("bowtie_prophages_DASins.fa", "a")
        outF.writelines(">" + prophage_regions[b1][0] + '\n')
        outF.writelines(bacterial_bin[prophage_regions[b1][0]][(int(prophage_regions[b1][1])):(int(prophage_regions[b1][2]))] + '\n')
        outF.close()

print('End')

