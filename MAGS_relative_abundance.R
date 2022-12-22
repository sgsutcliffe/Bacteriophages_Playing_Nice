#Open packages
library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(stringr)
library("ggsci")



#Relative abundances from CheckM
SRR828645_profile <- read.delim("SRR828645_profile", header = TRUE, col.names = c('BinId', 'BinSize', 'MappedReads', 'Perc_Mapped_Reads', 'Percen_binned_population', 'Perc_community'))
SRR828660_profile <- read.delim("SRR828660_profile", header = TRUE, col.names = c('BinId', 'BinSize', 'MappedReads', 'Perc_Mapped_Reads', 'Percen_binned_population', 'Perc_community'))
SRR828661_profile <- read.delim("SRR828661_profile", header = TRUE, col.names = c('BinId', 'BinSize', 'MappedReads', 'Perc_Mapped_Reads', 'Percen_binned_population', 'Perc_community'))

#Elsewhere I turned the GTDB_tk into a phylogentically sorted list
phylo_GTDB_tk <- read.csv("phylogeny_of_bins.csv", header = TRUE)
phylo_GTDB_tk$GTDB_taxa <- as.factor(phylo_GTDB_tk$GTDB_taxa)

#Create a normalized for genome size percent of mapped reads

SRR828645_profile$norm_reads <- (SRR828645_profile$MappedReads/(SRR828645_profile$BinSize*100000000))
SRR828645_norm_reads <- sum(SRR828645_profile$norm_reads)
SRR828645_profile$norm_perc_reads <- (SRR828645_profile$norm_reads/SRR828645_norm_reads)*100
SRR828645_mini <- cbind(SRR828645_profile$BinId, SRR828645_profile$norm_reads, (replicate(25, "week1")))
colnames(SRR828645_mini) <- c('BinID', 'NormPercReads', 'Week')
SRR828645_mini <- as.data.frame(SRR828645_mini)

SRR828661_profile$norm_reads <- (SRR828661_profile$MappedReads/(SRR828661_profile$BinSize*100000000))
SRR828661_norm_reads <- sum(SRR828661_profile$norm_reads)
SRR828661_profile$norm_perc_reads <- (SRR828661_profile$norm_reads/SRR828661_norm_reads)*100
SRR828661_mini <- cbind(SRR828661_profile$BinId, SRR828661_profile$norm_reads, (replicate(25, "week2")))
colnames(SRR828661_mini) <- c('BinID', 'NormPercReads', 'Week')
SRR828661_mini <- as.data.frame(SRR828661_mini)

SRR828660_profile$norm_reads <- (SRR828660_profile$MappedReads/(SRR828660_profile$BinSize*100000000))
SRR828660_norm_reads <- sum(SRR828660_profile$norm_reads)
SRR828660_profile$norm_perc_reads <- (SRR828660_profile$norm_reads/SRR828660_norm_reads)*100
SRR828660_mini <- cbind(SRR828660_profile$BinId, SRR828660_profile$norm_reads, (replicate(25, "week3")))
colnames(SRR828660_mini) <- c('BinID', 'NormPercReads', 'Week')
SRR828660_mini <- as.data.frame(SRR828660_mini)

all_weeks <- rbind(SRR828645_mini,SRR828660_mini,SRR828661_mini)
all_weeks <-  transform(all_weeks, NormPercReads = as.numeric(NormPercReads))
all_weeks$BinID <- phylo_GTDB_tk[(match(all_weeks$BinID,phylo_GTDB_tk$Bin)),9]
all_weeks$BinID <- factor(all_weeks$BinID, levels = phylo_GTDB_tk$GTDB_taxa)

stacked_norm_perc <- ggplot(all_weeks, aes(fill=BinID, y=NormPercReads, x=Week,)) +
  geom_bar(position="fill", stat="identity") + xlab(label = element_blank()) + ylab(label = "Relative Abundance") + scale_fill_igv() + theme_bw() +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +theme(panel.border = element_blank())

stacked_norm_perc  <- stacked_norm_perc + labs(fill = 'Bacteria') 
stacked_norm_perc
#save as .svg
ggsave(file="stacked_norm_perc.svg", plot=stacked_norm_perc, width=10, height=8)

#Alejandro wanted me to try an alternative order, specifically in order of the most abundant bacteria.

#I interpreted this as the most abundant bacteria over the three days, the noramlized by each day or percentage of reads.

abundance_order <- rbind(SRR828645_mini,SRR828660_mini,SRR828661_mini)
abundance_order$NormPercReads <- as.numeric(abundance_order$NormPercReads)
binorder <- abundance_order %>% group_by(BinID) %>% summarise(total_abund = sum(NormPercReads)) %>% arrange(desc(total_abund)) %>% select(BinID)
abundance_order <-  transform(abundance_order, NormPercReads = as.numeric(NormPercReads))
abundance_order$BinID <- factor(abundance_order$BinID, levels = binorder$BinID)
abundance_order$Bacteria <- phylo_GTDB_tk[(match(abundance_order$BinID,phylo_GTDB_tk$Bin)),9]
abundance_order$Bacteria <- as.character(abundance_order$Bacteria)
abundance_order$Bacteria <- factor(abundance_order$Bacteria, levels = phylo_GTDB_tk[(match(binorder$BinID,phylo_GTDB_tk$Bin)),9])

stacked_norm_perc2 <- ggplot(abundance_order, aes(fill=Bacteria, y=NormPercReads, x=Week)) +
  geom_bar(position="fill", stat="identity") + xlab(label = element_blank()) + ylab(label = "Relative Abundance") + scale_fill_igv() + theme_bw() +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.border = element_blank()) +
      + theme(text = element_text(size = 10))
stacked_norm_perc2 + theme(text = element_text(size = 10)) + theme(axis.title = element_text(size = 20)) + theme(legend.text = element_text(size = 20))     
stacked_norm_perc2 <- stacked_norm_perc2 + theme(text = element_text(size = 10)) + theme(axis.title = element_text(size = 20)) + theme(legend.text = element_text(size = 10)) + theme(axis.title.x = element_text(size = 15, angle=90, hjust=1, vjust=1))    
stacked_norm_perc2
ggsave(file="stacked_norm_perc2.svg", plot=stacked_norm_perc2 )#, width=10, height=8)

#I will try it based on abundance of week 1

abundance_order2 <- rbind(SRR828645_mini,SRR828660_mini,SRR828661_mini)
abundance_order2$NormPercReads <- as.numeric(abundance_order2$NormPercReads)
binorder2 <- abundance_order2 %>% group_by(BinID) %>% filter(Week == 'week1') %>% arrange(desc(NormPercReads))
abundance_order2 <-  transform(abundance_order2, NormPercReads = as.numeric(NormPercReads))
abundance_order2$BinID <- factor(abundance_order2$BinID, levels = binorder2$BinID)
abundance_order2$Bacteria <- phylo_GTDB_tk[(match(abundance_order2$BinID,phylo_GTDB_tk$Bin)),9]
abundance_order2$Bacteria <- as.character(abundance_order2$Bacteria)
abundance_order2$Bacteria <- factor(abundance_order2$Bacteria, levels = phylo_GTDB_tk[(match(binorder2$BinID,phylo_GTDB_tk$Bin)),9])
stacked_norm_perc3 <- ggplot(abundance_order2, aes(fill=Bacteria, y=NormPercReads, x=Week,)) +
  geom_bar(position="fill", stat="identity") + xlab(label = element_blank()) + ylab(label = "Relative Abundance") + scale_fill_igv() + theme_bw() +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +theme(panel.border = element_blank())

stacked_norm_perc3 

#The order for organizing bacteria for figures is 
bacteria_order <- as.data.frame(binorder$BinID)
bacteria_order$`binorder$BinID` <- phylo_GTDB_tk[(match(bacteria_order$`binorder$BinID`,phylo_GTDB_tk$Bin)),9]
write_csv(bacteria_order, file = 'bacteria_order.csv')

#Alternative colour schemes I could use
library(pals)
stacked_norm_perc + scale_fill_manual(values = rev(cols25(n=25)))
stacked_norm_perc + scale_fill_manual(values = cols25(n=26))
stacked_norm_perc + scale_fill_manual(values = (glasbey(n=25)))
stacked_norm_perc + scale_fill_manual(values = rev(glasbey(n=25)))
stacked_norm_perc + scale_fill_manual(values = rev(unname(polychrome(n=25))))
stacked_norm_perc + scale_fill_manual(values = rev(unname(alphabet(n=25))))
stacked_norm_perc + scale_fill_manual(values = (unname(alphabet(n=25))))

