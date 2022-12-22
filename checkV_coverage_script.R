#Open packages
library(dplyr) #Version ‘1.0.9’
library(ggplot2) #Version ‘3.3.6’
library(tidyr) #Version ‘1.2.0’
library(stringr) #Version ‘1.4.0’
library(reshape2) #Version ‘1.4.4’


setwd("/Users/Sutcliffe/Library/CloudStorage/OneDrive-McGillUniversity/PhD_OneDrive/Publication_3_Viral_Metagenomics_Collab/Writing/12_draft_revisions/Scripts")

#Load in the data! 

SRR <- read.csv("Data_For_Scripts/sequence_runs.csv", header=TRUE) #File has the sequence run meta-data
runinfo <- read.csv("Data_For_Scripts/sequence_runs.csv", header = TRUE)
bacterial_bin <- read.csv(file = 'Data_For_Scripts/GTDB-Tk_DAS_check_bin_modified.csv', header = FALSE)
bacphlip <- read.delim('Data_For_Scripts/CheckV_Coverage/checkv_contigs.fasta.bacphlip', header = TRUE, sep = '\t')
setwd("Data_For_Scripts/CheckV_Coverage/")
temp <-  list.files(pattern="*_RPKM_coverage") #Make a temp variable for each of the coverage file per sequence run. The all viral reads was not included unless I need it.
myfiles <-  lapply(temp, function(x) read.delim(x, header = TRUE)) #Opens up each file
coverage_data <- do.call(rbind, myfiles) #Binds all the files together in one dataframe

##### DESEQ normalized reads ####

#I will run the analysis on DESEQ normalized reads
#Load in the normalized reads
DESEQ_normalized_abundance <- read.csv("../../Data_For_Scripts/norDESEQ_counts_norbylength.tsv", sep = "\t", header=TRUE) #DESEQ normalized read counts

#Then select the 462 contigs we have lifestyle prediction for by Bacphlip
bacphlip_contigs <- match(DESEQ_normalized_abundance$contig, bacphlip$X)
bacphlip_contigs <- bacphlip_contigs[!is.na(bacphlip_contigs)]
bacphlip_norm_abun <- DESEQ_normalized_abundance[bacphlip_contigs,]

#Now I will add the 48 prophages not in the list as we know they are temperate too
bacphlip_norm_abun <- rbind(DESEQ_normalized_abundance[1:52,], bacphlip_norm_abun)
#Now we have some duplicates because I added the 52 prophages but there were already four there
bacphlip_norm_abun <- bacphlip_norm_abun[!duplicated(bacphlip_norm_abun$contig),]

#Now add a column for lifestyle

#Use the bacphlip data to determine the lifestyle 
temperate_phages <- bacphlip$X[(which(bacphlip$Temperate > 0.5))]
bacphlip_norm_abun <- bacphlip_norm_abun %>% mutate(lifestyle = case_when((!(contig %in% temperate_phages) & !(str_detect(contig, "DAS"))) ~ 'lytic', ((contig %in% temperate_phages) | ( str_detect(contig, "DAS"))) ~ 'temperate' ) )
row.names(bacphlip_norm_abun) <- bacphlip_norm_abun[,1]
bacphlip_norm_abun <- bacphlip_norm_abun[,-1]

temperate_contigs3 <- bacphlip_norm_abun %>% filter(lifestyle == "temperate")
lytic_contigs3 <- bacphlip_norm_abun %>% filter(lifestyle == "lytic")

relative_by_lifstyl <- c()
for (column in 1:24){
  temperate <- sum(temperate_contigs3[,column])/sum(bacphlip_norm_abun[,column])
  lytic <- sum(lytic_contigs3[,column])/sum(bacphlip_norm_abun[,column])
  relative_by_lifstyl <- rbind(relative_by_lifstyl, c(temperate, lytic))
}

relative_by_lifstyl <- as.data.frame(relative_by_lifstyl)
relative_by_lifstyl <- cbind(colnames(bacphlip_norm_abun)[1:24], relative_by_lifstyl)
colnames(relative_by_lifstyl) <- c('SRR','temperate', 'lytic') 

relative_by_lifstyl$ID <- runinfo$ID[match(relative_by_lifstyl$SRR,runinfo$SRR)]
relative_by_lifstyl_melt <- reshape2::melt(relative_by_lifstyl, id.vars=c("SRR", 'ID'))

DESEQ_lifestyle_abund <- ggplot(relative_by_lifstyl_melt, aes(fill=variable, y=value, x=ID)) + geom_bar(position="stack", stat="identity") + theme(axis.text.x = element_text(angle = 90)) + labs( y = 'Normalized Relative Abundance (%)')
DESEQ_lifestyle_abund <- DESEQ_lifestyle_abund + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                          panel.grid.minor = element_blank())

#Chnage titles of days on x-axis + reformatting to make more clear to read
day_sampled <- c(runinfo$Time_days)
DESEQ_lifestyle_abund3 <- DESEQ_lifestyle_abund + scale_x_discrete(labels= day_sampled) + xlab(label = 'Day') + theme(axis.text.x = element_text(size = 25), axis.text.y = element_text(size = 25)) + theme(strip.text = element_text(size =35)) + theme(axis.title = element_text(size = 35)) + scale_fill_discrete(name = "")
DESEQ_lifestyle_abund3 


##### STATS #####
#Note couldn't get the stats to work in R so I finished it in PRISM
relative_by_lifstyl_melt <- as_tibble(relative_by_lifstyl_melt)
relative_by_lifstyl_melt$day <- runinfo$ID[match(relative_by_lifstyl_melt$SRR, runinfo$SRR)]
relative_by_lifstyl_melt$week <- runinfo$week[match(relative_by_lifstyl_melt$SRR, runinfo$SRR)]
relative_by_lifstyl_melt$day_group <- runinfo$Time_days[match(relative_by_lifstyl_melt$SRR, runinfo$SRR)]
stats <- relative_by_lifstyl_melt %>% select(value, day_group, variable) %>% filter(day_group %in% c(0, 182, 851, 852, 853,879, 880, 881))
stats <- filter(variable == 'temperate')
