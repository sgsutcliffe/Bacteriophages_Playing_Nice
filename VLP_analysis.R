#Open packages
library(dplyr) #Version ‘1.0.9’
library(ggplot2) #Version ‘3.3.6’
library(tidyr) #Version ‘1.2.0’
library(stringr) #Version ‘1.4.0’
library(permute) #Version ‘0.9.7’
library(vegan) #Version ‘2.6.2’
library(data.table) #Version ‘1.14.2’
library(lattice) #Version ‘0.20.45’
library(rstatix) #Version ‘0.7.0’
library(grid) #Version ‘4.2.1’
library(DESeq2) #Version 1.36.0’
library(car) #Version ‘3.1.0’
library(ape) #Version ‘5.6.2’
library(reshape2) #Version ‘1.4.4’
library(IRanges) #Version ‘2.30.0’
library(GenomicRanges) #Version ‘1.48.0’

#Working from data
setwd("VLP_analysis/")
#Open all files in a data directory
temp <-  list.files(pattern="*_readcount") #Make a temp variable for each of the readcount files
SRR_names <- substring(temp,1,9) #Make a vector of the SRR run by name in order files of temp
myfiles <-  lapply(temp, function(x) read.delim(x, header = FALSE)) #Opens up each file
myfiles2 <- do.call(cbind, myfiles) #Binds all the files together in one dataframe
colnames(myfiles2) <- 1:48 #Give each column a number
counts_df <- myfiles2[,seq(2,48,2)] #Only odd numbered columns have read data for some reason
colnames(counts_df) <- SRR_names #The runs are in ordered based on temp file so name them
Contig <- myfiles2[,1] #use first column to get contigs name
bacterial_bin <- read.csv(file = '../GTDB-Tk_DAS_check_bin_modified.csv', header = FALSE)
demovir <- read.table(file = "demovir_taxonomy.tabular", sep = '\t', header = TRUE)
runinfo <- read.csv("../sequence_runs.csv", header = TRUE)

count_matrix <- cbind(Contig, counts_df)

totalNumReads <- read.delim("totalNumReads", header = FALSE, col.names = c("SRR", "Num_Reads"))
contigLength <- read.delim("ContigLength", header = FALSE, col.names = c("Contig", "Length"))

### RPKM ####

#reads per million
totalNumReads <- totalNumReads %>%
  mutate(readsPM = Num_Reads/1000000)

#contig per KB
contigLength <-  contigLength %>%
  mutate(contigLengthperKB = Length/1000)

#Merge based on contig data
merged_count_matrix <- merge(contigLength, count_matrix, by="Contig")

#Divide counts by contigLength per Kb
RPK <- merged_count_matrix[,4:27]/merged_count_matrix$contigLengthperKB
contignames <- merged_count_matrix$Contig
temp <-  RPK
for (row in 1:nrow(totalNumReads)){
  index_name <- (totalNumReads$SRR[row])
  index_number <- match(index_name, colnames(RPK))
  new_column <- RPK[,index_number]/totalNumReads$readsPM[row]
  temp[,index_number] <- new_column}

RPKM <- cbind(contignames,temp)

setwd('../results')
write.table(RPKM, 'RPKM.tsv', sep = "\t", row.names=FALSE)


#Heatmap

#1 Get it in the right format (x-axis: day/seqrun y-axis: contig z: RPKM)
RPKM <- read.table('RPKM.tsv', sep = "\t", header = TRUE)

RPKM.long <- pivot_longer(data = RPKM, cols = -c(1), names_to = 'SRR', values_to = 'RPKM')
colnames(RPKM.long)[which(names(RPKM.long) == 'contignames')] <- 'contig'
RPKM.long <- RPKM.long %>% 
  mutate(lifestyle = case_when(str_detect(contig, "NODE") ~ 'lytic', str_detect(contig, "DAS") ~ 'prophage' ) )


write.table(RPKM.long, 'RPKM_long.tsv', sep = "\t")
RPKM.long <- read.table('RPKM_long.tsv', sep = "\t")

##### log-transform RPKM abundances #####

#removing all zeros for log transformation 
RPKM.long[RPKM.long == 0] <- NA

RPKM.long$log.abundance <- log(RPKM.long$RPKM)

write.table(RPKM.long, 'logRPKM_long.tsv', sep = "\t")
RPKM.long <- read.table('logRPKM_long.tsv', sep = "\t")

logRPKM.heatmap <- ggplot(data = RPKM.long, aes(x = SRR, y = contig, fill = log.abundance)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue",
                       mid = "white",
                       high ="red",
                       midpoint = 0,
                       space = "Lab",
                       na.value = "blue") +
  theme(axis.text.x = element_text(angle = 90), axis.text.y = element_text(size = 0.4)) +
  xlab(label = "Sequence Run")


logRPKM.heatmap

##### RELATIVE ABUNDANCE ####

relative_abundance <- read.table('RPKM.tsv', sep = "\t", header = TRUE)


for (column in 2:ncol(relative_abundance)){
  total_normalized_counts_per_sequence_run <- sum(relative_abundance[,column])
  relative_abundance[,column] <- relative_abundance[,column]/total_normalized_counts_per_sequence_run
}

#### RELATIVE ABUNDANCE ###
write.table(relative_abundance, 'relative_abundance.tsv', sep = "\t")
relative_abundance <- read.table('relative_abundance.tsv', sep = "\t")

#### 1% RELATIVE ABUNDANCE ####

#I am going to run a check to see how many contigs actually contribute more than 1% of reads to a sample
number_of_one_percent_contigs <- 0
list_contigs_greater_than_1percent <- c()
for (row in 1:nrow(relative_abundance)){
  
  if (max(relative_abundance[row,-1]) > 0.01) {
    number_of_one_percent_contigs <- number_of_one_percent_contigs + 1
    list_contigs_greater_than_1percent <- c(list_contigs_greater_than_1percent, relative_abundance[row,1] )
  }
  
}

#### 0.1% RELATIVE ABUNDANCE ####

#I am going to run a check to see how many contigs actually contribute more than 1% of reads to a sample
number_of_0.1_percent_contigs <- 0
list_contigs_greater_than_0.1percent <- c()

for (row in 1:nrow(relative_abundance)){
  
  if (max(relative_abundance[row,-1]) > 0.001) {
    number_of_0.1_percent_contigs <- number_of_0.1_percent_contigs + 1
    list_contigs_greater_than_0.1percent <- c(list_contigs_greater_than_0.1percent, relative_abundance[row,1] )
  }
  
}

number_of_0.1_percent_contigs 
list_contigs_greater_than_0.1percent

#0.1%
logRPKM.long <- read.table('logRPKM_long.tsv', sep = "\t")
logRPKM.long <- as_tibble(logRPKM.long)
zero_one_percent_RPKM <- logRPKM.long[logRPKM.long$contig %in% list_contigs_greater_than_0.1percent,]
write.table(zero_one_percent_RPKM, 'zero_one_percent_sam.tsv', sep = "\t")

##Make the sequence runs in order using factors
zero_one_percent_RPKM$SRR <- factor(zero_one_percent_RPKM$SRR, levels = c(runinfo$SRR))


logRPKM.long.heatmap <- ggplot(data = logRPKM.long, aes(x = SRR, y = contig, fill = log.abundance)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue",
                       mid = "white",
                       high ="red",
                       midpoint = 0,
                       space = "Lab",
                       na.value = "blue") +
  theme(axis.text.x = element_text(angle = 90), axis.text.y = element_text(size = 0.4)) +
  xlab(label = "Sequence Run")


#logRPKM.long.heatmap

##Look at only prophages and 0.1%
zero_one_percent_RPKM <- logRPKM.long[logRPKM.long$contig %in% list_contigs_greater_than_0.1percent,]
zero_one_percent_RPKM$SRR <- factor(zero_one_percent_RPKM$SRR, levels = c(runinfo$SRR))
prophage_zero_one_percent_RPKM <- zero_one_percent_RPKM %>% filter(zero_one_percent_RPKM$lifestyle == "prophage")

prophage_zero_one_percent_RPKM.heatmap <- ggplot(data = prophage_zero_one_percent_RPKM, aes(x = SRR, y = contig, fill = log.abundance)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue",
                       mid = "white",
                       high ="red",
                       midpoint = 0,
                       space = "Lab",
                       na.value = "blue") +
  theme(axis.text.x = element_text(angle = 90), axis.text.y = element_text(size = 0.4)) +
  xlab(label = "Sequence Run")


#prophage_zero_one_percent_RPKM.heatmap

#Diversity

#change names of SRR to days
x <- runinfo$ID[match(colnames(RPKM), runinfo$SRR)]
colnames(RPKM) <- replace_na(x,'contig')
shannon_diversity <-  c()
for (column in 2:(ncol(RPKM)-1)) {
  temp <- diversity(RPKM[,column], MARGIN = 2, index = 'shannon')
  shannon_diversity <- c(shannon_diversity, temp)
}

#Bray-curtis difference
RPKM <- read.table('RPKM.tsv', sep = "\t", header = TRUE)
row.names(RPKM) <- RPKM[,1]
RPKM <- RPKM[,-1]
colnames(RPKM) <- runinfo$ID[match(colnames(RPKM), runinfo$SRR)]
RPKM_transposed <- as.data.frame(t(RPKM))
beta_dist <- vegdist(RPKM_transposed, index = "bray")
beta_dist_clust <- hclust(beta_dist, method = "average")
plot(beta_dist_clust)

colnames(RPKM) <- runinfo$Time_days[match(colnames(RPKM), runinfo$SRR)]
RPKM_transposed <- as.data.frame(t(RPKM))

#Euclidean distance
RPKM <- read.table('RPKM.tsv', sep = "\t", header = TRUE)
row.names(RPKM) <- RPKM[,1]
RPKM <- RPKM[,-1]
colnames(RPKM) <- runinfo$Time_days[match(colnames(RPKM), runinfo$SRR)]
RPKM_transposed <- as.data.frame(t(RPKM))
euclid_dist <- as.matrix(dist(RPKM_transposed, method = "euclidean"))
colnames(euclid_dist) <- colnames(RPKM)
row.names(euclid_dist) <- colnames(RPKM)
euclid_matrix <- as.data.frame(euclid_dist)
euclid_matrix <- cbind(colnames(euclid_matrix), euclid_matrix)
euclid_long <- pivot_longer(euclid_matrix, cols = euclid_matrix$`colnames(euclid_matrix)`)
colnames(euclid_long) <- c('time1', 'time2', 'euc_dis')
euclid_long$time1 <- as.numeric(euclid_long$time1)
euclid_long$time2 <- as.numeric(euclid_long$time2)
euclid_long$srt_time_lag <- sqrt(abs((euclid_long$time1 - euclid_long$time2)))
q <- plot(euclid_long$srt_time_lag, euclid_long$euc_dis)
q <- lm(euclid_long$euc_dis ~ euclid_long$srt_time_lag)
q <- abline(9921, 1557)
q
##the results were that there was 'directional change' in the virome unlike what I would expect
#So I will look at the top most abundant 

RPKM <- read.table('RPKM.tsv', sep = "\t", header = TRUE)
top100 <- c()
for (column in 2:25){
  x <- head(arrange(RPKM, desc(RPKM[,column])), n =175)[,1]
  top100 <- c(top100, x)
}

most_abundant <- unique(top100)
rpkm_most_abundant <- RPKM_transposed[, most_abundant]
most_abdnt_euclid_dist <- as.matrix(dist(rpkm_most_abundant, method = "euclidean"))

RPKM <- read.table('RPKM.tsv', sep = "\t", header = TRUE)
row.names(RPKM) <- RPKM[,1]
RPKM <- RPKM[,-1]
colnames(RPKM) <- runinfo$Time_days[match(colnames(RPKM), runinfo$SRR)]

colnames(most_abdnt_euclid_dist) <- colnames(RPKM)
row.names(most_abdnt_euclid_dist) <- colnames(RPKM)
most_abd_euclid_matrix <- as.data.frame(most_abdnt_euclid_dist)
most_abd_euclid_matrix <- cbind(colnames(most_abd_euclid_matrix), most_abd_euclid_matrix)
most_abd_euclid_long <- pivot_longer(most_abd_euclid_matrix, cols = most_abd_euclid_matrix$`colnames(most_abd_euclid_matrix)`)
colnames(most_abd_euclid_long) <- c('time1', 'time2', 'euc_dis')
most_abd_euclid_long$time1 <- as.numeric(most_abd_euclid_long$time1)
most_abd_euclid_long$time2 <- as.numeric(most_abd_euclid_long$time2)
most_abd_euclid_long$srt_time_lag <- sqrt(abs((most_abd_euclid_long$time1 - most_abd_euclid_long$time2)))
plot(most_abd_euclid_long$srt_time_lag, most_abd_euclid_long$euc_dis)
lm(most_abd_euclid_long$euc_dis ~ most_abd_euclid_long$srt_time_lag)

abline(9919, 1557)


##Next I will look at how the active prophages compare

RPKM <- read.table('RPKM.tsv', sep = "\t", header = TRUE)
RPKM <- as_tibble(RPKM)
RPKM <- RPKM %>% mutate(lifestyle = case_when(str_detect(contignames, "NODE") ~ 'lytic', str_detect(contignames, "DAS") ~ 'temperate' ) )

#I know it's a bad move but first 52 coloumns of RPKM_transposed are active prophages and the rest are the regular ones
active_prophage_RPKM <- RPKM_transposed[,1:52]
non_prophage_RPKM <- RPKM_transposed[,53:5890]

active_prophage_euclid_dist <- as.matrix(dist(active_prophage_RPKM, method = "euclidean"))

RPKM <- read.table('RPKM.tsv', sep = "\t", header = TRUE)
row.names(RPKM) <- RPKM[,1]
RPKM <- RPKM[,-1]
colnames(RPKM) <- runinfo$Time_days[match(colnames(RPKM), runinfo$SRR)]

colnames(active_prophage_euclid_dist) <- colnames(RPKM)
row.names(active_prophage_euclid_dist) <- colnames(RPKM)
active_prophage_euclid_matrix <- as.data.frame(active_prophage_euclid_dist)
active_prophage_euclid_matrix <- cbind(colnames(active_prophage_euclid_matrix), active_prophage_euclid_matrix)
active_prophage_euclid_long <- pivot_longer(active_prophage_euclid_matrix, cols = active_prophage_euclid_matrix$`colnames(active_prophage_euclid_matrix)`)
colnames(active_prophage_euclid_long) <- c('time1', 'time2', 'euc_dis')
active_prophage_euclid_long$time1 <- as.numeric(active_prophage_euclid_long$time1)
active_prophage_euclid_long$time2 <- as.numeric(active_prophage_euclid_long$time2)
active_prophage_euclid_long$srt_time_lag <- sqrt(abs((active_prophage_euclid_long$time1 - active_prophage_euclid_long$time2)))
plot(active_prophage_euclid_long$srt_time_lag, active_prophage_euclid_long$euc_dis)
lm(active_prophage_euclid_long$euc_dis ~ active_prophage_euclid_long$srt_time_lag)
abline(139.8442,-0.7071)

#Non-active prophages
non_prophage_euclid_dist <- as.matrix(dist(non_prophage_RPKM, method = "euclidean"))

RPKM <- read.table('RPKM.tsv', sep = "\t", header = TRUE)
row.names(RPKM) <- RPKM[,1]
RPKM <- RPKM[,-1]
colnames(RPKM) <- runinfo$Time_days[match(colnames(RPKM), runinfo$SRR)]

colnames(non_prophage_euclid_dist) <- colnames(RPKM)
row.names(non_prophage_euclid_dist) <- colnames(RPKM)
non_prophage_euclid_matrix <- as.data.frame(non_prophage_euclid_dist)
non_prophage_euclid_matrix <- cbind(colnames(non_prophage_euclid_matrix), non_prophage_euclid_matrix)
non_prophage_euclid_long <- pivot_longer(non_prophage_euclid_matrix, cols = non_prophage_euclid_matrix$`colnames(non_prophage_euclid_matrix)`)
colnames(non_prophage_euclid_long) <- c('time1', 'time2', 'euc_dis')
non_prophage_euclid_long$time1 <- as.numeric(non_prophage_euclid_long$time1)
non_prophage_euclid_long$time2 <- as.numeric(non_prophage_euclid_long$time2)
non_prophage_euclid_long$srt_time_lag <- sqrt(abs((non_prophage_euclid_long$time1 - non_prophage_euclid_long$time2)))
plot(non_prophage_euclid_long$srt_time_lag, non_prophage_euclid_long$euc_dis)
lm(non_prophage_euclid_long$euc_dis ~ non_prophage_euclid_long$srt_time_lag)
abline(9919,1557)


##Now try and use bray-curtis distance

#Bray-curtis distance
RPKM <- read.table('RPKM.tsv', sep = "\t", header = TRUE)
row.names(RPKM) <- RPKM[,1]
RPKM <- RPKM[,-1]
colnames(RPKM) <- runinfo$Time_days[match(colnames(RPKM), runinfo$SRR)]
RPKM_transposed <- as.data.frame(t(RPKM))
beta_dist <- as.matrix(vegdist(RPKM_transposed, index = "bray"))
colnames(beta_dist) <- colnames(RPKM)
row.names(beta_dist) <- colnames(RPKM)
beta_matrix <- as.data.frame(beta_dist)
beta_matrix <- cbind(colnames(beta_matrix), beta_matrix)
beta_long <- pivot_longer(beta_matrix, cols = beta_matrix$`colnames(beta_matrix)`)
colnames(beta_long) <- c('time1', 'time2', 'beta_dis')
beta_long$time1 <- as.numeric(beta_long$time1)
beta_long$time2 <- as.numeric(beta_long$time2)
beta_long$srt_time_lag <- sqrt(abs((beta_long$time1 - beta_long$time2)))
plot(beta_long$srt_time_lag, beta_long$beta_dis)
lm(beta_long$beta_dis ~ beta_long$srt_time_lag)
abline(0.18549, 0.02483)



#Look at the most abundant contigs

RPKM <- read.table('RPKM.tsv', sep = "\t", header = TRUE)
top100 <- c()
for (column in 2:25){
  x <- head(arrange(RPKM, desc(RPKM[,column])), n =100)[,1]
  top100 <- cbind(top100, x)
  }

colnames(top100) <- colnames(RPKM)[-1]

#Look at abundance of active prophges compared to non-active.

#Make a column for prophage or not
RPKM <- as_tibble(RPKM)
RPKM <- as.tibble(RPKM)
RPKM <- RPKM %>% mutate(lifestyle = case_when(str_detect(contignames, "NODE") ~ 'lytic', str_detect(contignames, "DAS") ~ 'temperate' ) )
test <- aggregate(RPKM[,2:25], by=list(Category=RPKM$lifestyle), FUN=sum) #This sums the RPKM for all the contigs based on being temperate or lytic
test <- t(test)
colnames(test) <- test[1,]
test <- test[-1,]
test <- as.data.frame(test)
test$lytic <- as.numeric(test$lytic)
test$temperate <- as.numeric(test$temperate)
test$total_RPKM <- test$lytic + test$temperate
test$perc_lytic <- test$lytic / test$total_RPKM
test$perc_temp <- test$temperate / test$total_RPKM

#The number of contigs is a lot more in lytic so I decided to explore dividing it by the number of contigs
amnt_temp <- sum(RPKM$lifestyle == "temperate")
amnt_lytic <- sum(RPKM$lifestyle == "lytic")

test$lytic_averag <- test$lytic / amnt_lytic
test$temp_averag <- test$temperate / amnt_temp

test$day <- runinfo$Time_days[match(row.names(test), runinfo$SRR)]
test$week <- runinfo$week[match(row.names(test), runinfo$SRR)]
test_graph <- ggplot(test, aes(x=day, y=perc_temp)) + geom_point(aes(size = 30)) + theme(axis.text.x = element_text(size = 20))  
test_graph <- test_graph +  facet_wrap(~ week, scales="free_x") + ylab(label = 'Relative Abundance of Active Prophages') + xlab(label = 'Day') + theme(strip.text = element_text(size =20)) + theme(axis.title = element_text(size = 30))
  
test_graph
ggsave(file="/Users/steven/GutMicrobiome/Publication_3_Viral_Metagenomics_Collab/Minot_Data/Third_Round_bins/Viral_Analysis/Reads_Per_Contig/figures/active_prophage_relative_abundance.svg", plot=test_graph, width=10, height=10)


##### DESEQ2 #####

#format count data for DESEQ2
#DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ condition)
#countData
cts <- count_matrix[,-1]
rownames(cts) <- count_matrix[,1]
head(cts,2)

#colData
orderSamples <- match(colnames(cts), runinfo$SRR)
days <- c((runinfo$Time_days[orderSamples]))
weeks <- c(runinfo$week[orderSamples])
coldata <- cbind(days, weeks)
colnames(coldata) <- c('days', 'week')
rownames(coldata) <- runinfo$SRR[orderSamples]

coldata <- as.data.frame(coldata)

coldata$days <- factor(coldata$days, levels = unique(sort(coldata$days)))
coldata$week <- factor(coldata$week, levels = unique(sort(coldata$week)))
#Test that order is correct which is important for DESEQ
all(rownames(coldata) %in% colnames(cts))
all(rownames(coldata) == colnames(cts))

#Now we can run it as stated above
dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ days)

#Now we can run it
test <- DESeq(dds)
res <- results(test)
res
summary(res)
subset(res, padj<0.05) %>% order(row.names(res))
head(res[order(res$padj),])
by_contig <- subset((subset(res[order(row.names(res)),], padj<0.05)), log2FoldChange>1)

##this last step highlights the active prophages 
by_prophage <- subset(by_contig, str_detect(row.names(by_contig), "DAS"))
prophages_induced <- data.frame(
  row.names(by_prophage), 
  by_prophage$baseMean,
  by_prophage$log2FoldChange,
  by_prophage$lfcSE,
  by_prophage$stat,
  by_prophage$pvalue,
  by_prophage$padj)
colnames(prophages_induced) <- c("Phage", "baseMean", "log2FoldChange", "lfcSE", 
                   "stat", "P-value","Adjusted p-value")
#write.table(prophages_induced, 'DESEQ_induced_prophages.tsv', sep = "\t")

vsdata <- varianceStabilizingTransformation(dds)
day_PCA <- plotPCA(vsdata, intgroup="days")
week_PCA <- plotPCA(vsdata, intgroup="week") + geom_point(size=0.5) +
  theme(axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 7), axis.title = element_text(size = 7)) + 
  theme(legend.key.size = unit(0.25, 'cm'), #change legend key size
        legend.key.height = unit(0.25, 'cm'), #change legend key height
        legend.key.width = unit(0.25, 'cm'), #change legend key width
        legend.title = element_text(size=7), #change legend title font size
        legend.text = element_text(size=7)) #change legend text font size
week_PCA
ggsave(file="day_PCA2.svg", plot=day_PCA, width=10, height=10)
ggsave(file="week_PCA2_2.svg", plot=week_PCA, width=8, height=8, units = c("cm"), dpi = 600)

#Now look at only the active prophages

dds_prophages <- dds[1:52]
vsdata_prophages <- varianceStabilizingTransformation(dds_prophages)

plot_PCA_prophages <- plotPCA(vsdata_prophages, intgroup="week") +
  geom_point(size=0.5) +
  theme(axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 7), axis.title = element_text(size = 7)) + 
  theme(legend.key.size = unit(0.25, 'cm'), #change legend key size
        legend.key.height = unit(0.25, 'cm'), #change legend key height
        legend.key.width = unit(0.25, 'cm'), #change legend key width
        legend.title = element_text(size=7), #change legend title font size
        legend.text = element_text(size=7)) #change legend text font size
plot_PCA_prophages 
  

ggsave(file="plot_PCA_all_phages.svg", plot=plot_PCA_all_phages, width=10, height=10)
ggsave(file="plot_PCA_prophages2.svg", plot=plot_PCA_prophages, width=8, height=10, units = c("cm"), dpi = 600)


#Normalizing counts by DESEQ count
dds <- estimateSizeFactors(dds)
sizeFactors(dds)

norDESEQ_counts <- counts(dds, normalized=TRUE)
norDESEQ_counts <- as.data.frame(norDESEQ_counts)
norDESEQ_counts <-  setDT(norDESEQ_counts, keep.rownames = 'contig')[]
norDESEQ_counts <- norDESEQ_counts[-5891,]

#Before also normalizing for contig length I will make sure both contig list and norDESEQ_counts are in same order of contigs
all(contigLength$Contig %in% norDESEQ_counts$contig)
all(contigLength$Contig == norDESEQ_counts$contig)

#Normalize counts per contig length
temp_dataframe <- norDESEQ_counts[,2:25]/contigLength$Length

norDESEQ_counts <- cbind(norDESEQ_counts$contig, temp_dataframe)
x <- colnames(norDESEQ_counts)
x[1] <- 'contig'
colnames(norDESEQ_counts) <- x

#Get relative abundance
norDESEQ_counts <- as.data.frame(norDESEQ_counts)
for (x in 2:ncol(norDESEQ_counts)){
  total_normalized_counts_per_sequence_run <- sum(norDESEQ_counts[,x])
  norDESEQ_counts[,x] <- norDESEQ_counts[,x]/total_normalized_counts_per_sequence_run
}


norDESEQ_counts <- as.data.table(norDESEQ_counts)
temperate_DESEQ <- norDESEQ_counts %>% mutate(lifestyle = case_when(str_detect(contig, "NODE") ~ 'lytic', str_detect(contig, "DAS") ~ 'temperate' ) )
temperate_DESEQ <- temperate_DESEQ %>% filter(lifestyle == 'temperate') %>% select(!(lifestyle))
temperate_DESEQ_melt <- melt(temperate_DESEQ, id.vars = c("contig"))


temperate_DESEQ_melt$variable <- factor(temperate_DESEQ_melt$variable, levels = c(runinfo$SRR))
temperate_DESEQ_melt$day <- as.factor(runinfo$ID[match(temperate_DESEQ_melt$variable, runinfo$SRR)])

prophage_induction <- temperate_DESEQ_melt[,c('contig','value','variable')]
colnames(prophage_induction) <- c("contig", "read_count", "SRR_columns")

prophage_induction <- spread(prophage_induction, contig, read_count)
prophage_induction <- as.data.frame(prophage_induction)
prophage_induction <- prophage_induction[-25,]
normalized_coverage <- cbind(runinfo$Time_days, runinfo$week, runinfo$ID, prophage_induction)

x <- colnames(normalized_coverage)
x[1] <- 'day'
x[2] <-  'week'
x[3] <- 'day-replicate'
colnames(normalized_coverage) <- x
normalized_coverage$`day-replicate` <- factor(normalized_coverage$`day-replicate`)
normalized_coverage$SRR_columns <- factor(normalized_coverage$SRR_columns)
normalized_coverage <- as.data.table(normalized_coverage)
normalized_coverage_melted <- melt(normalized_coverage, id.vars= c('day', 'week', 'day-replicate', 'SRR_columns'))

normalized_coverage_melted$variable <- str_match(normalized_coverage_melted$variable, "DAS_Check_Bin_..")
bins_position <- match(normalized_coverage_melted$variable, bacterial_bin$V1)
bin_names <- bacterial_bin$V9[bins_position]
normalized_coverage_melted$variable <- bin_names
p <- ggplot(normalized_coverage_melted, aes(x=day, y=value, col=variable)) + geom_point(aes(size = 15)) + theme(axis.text.x = element_text(size = 20))  
p <- p + facet_wrap(~ week, scales="free_x") + ylab(label = 'Normalized Coverage') + theme(strip.text = element_text(size =20)) + theme(axis.title = element_text(size = 30))
p

#Now I will revisit the DESEQ normalized counts and look at what percentage is active prophages.
norDESEQ_counts <- counts(dds, normalized=TRUE)
norDESEQ_counts <- as.data.frame(norDESEQ_counts)
norDESEQ_counts <-  setDT(norDESEQ_counts, keep.rownames = 'contig')[]
norDESEQ_counts <- norDESEQ_counts[-5891,]
write.table(norDESEQ_counts, 'norDESEQ_counts.tsv', sep = "\t", row.names=FALSE)
norDESEQ_counts <- read.table('norDESEQ_counts.tsv', sep = "\t", header =TRUE)

#Before also normalizing for contig length I will make sure both contig list and norDESEQ_counts are in same order of contigs
all(contigLength$Contig %in% norDESEQ_counts$contig)
all(contigLength$Contig == norDESEQ_counts$contig)

#Normalize counts per contig length
temp_dataframe <- norDESEQ_counts[,2:25]/contigLength$Length

norDESEQ_counts <- cbind(norDESEQ_counts$contig, temp_dataframe)
x <- colnames(norDESEQ_counts)
x[1] <- 'contig'
colnames(norDESEQ_counts) <- x
write.table(norDESEQ_counts, 'norDESEQ_counts_norbylength.tsv', sep = "\t", row.names=FALSE)
norDESEQ_counts <- read_delim('norDESEQ_counts_norbylength.tsv', delim = "\t")
norDESEQ_counts <- norDESEQ_counts %>% mutate(lifestyle = case_when(str_detect(contig, "NODE") ~ 'lytic', str_detect(contig, "DAS") ~ 'temperate' ) )

perc_active_proph <- aggregate(norDESEQ_counts [,2:25], by=list(Category=norDESEQ_counts$lifestyle), FUN=sum) #This sums the norDESEQ_counts for all the contigs based on being temperate or lytic
perc_active_proph <- t(perc_active_proph)
colnames(perc_active_proph) <- perc_active_proph[1,]
perc_active_proph <- perc_active_proph[-1,]
perc_active_proph <- as.data.frame(perc_active_proph, stringsAsFactors = FALSE)
perc_active_proph$lytic <- as.numeric(perc_active_proph$lytic)
perc_active_proph$temperate <- as.numeric(perc_active_proph$temperate)
perc_active_proph$total_norDESEQ_counts <- perc_active_proph$lytic + perc_active_proph$temperate
perc_active_proph$perc_lytic <- perc_active_proph$lytic / perc_active_proph$total_norDESEQ_counts
perc_active_proph$perc_temp <- perc_active_proph$temperate / perc_active_proph$total_norDESEQ_counts

#The number of contigs is a lot more in lytic so I decided to explore dividing it by the number of contigs
amnt_temp <- sum(norDESEQ_counts$lifestyle == "temperate")
amnt_lytic <- sum(norDESEQ_counts$lifestyle == "lytic")

perc_active_proph$lytic_averag <- perc_active_proph$lytic / amnt_lytic
perc_active_proph$temp_averag <- perc_active_proph$temperate / amnt_temp

perc_active_proph$day <- runinfo$ID[match(row.names(perc_active_proph), runinfo$SRR)]
perc_active_proph$week <- runinfo$week[match(row.names(perc_active_proph), runinfo$SRR)]
perc_active_proph$day_group <- runinfo$Time_days[match(row.names(perc_active_proph), runinfo$SRR)]
perc_active_proph$day <- as.factor(runinfo$ID)

perc_active_proph_graph <- ggplot(perc_active_proph, aes(x=day_group, y=perc_temp, color = perc_temp < 0.010)) + geom_point(size = 1, show.legend=FALSE) + theme(axis.text.x = element_text(size = 25), axis.text.y = element_text(size = 25))  
perc_active_proph_graph <- perc_active_proph_graph + ylab(label = 'Relative Abundance of Active Prophages') + xlab(label = 'Day') + theme(strip.text = element_text(size =35)) + theme(axis.title = element_text(size = 35))+
  geom_hline(yintercept=0.010, linetype='dashed', color='#F8766D', size = 0.5) + facet_grid(~ week, scales="free_x", space="free")

#perc_active_proph_graph 

#### STATS ON PERCENT ACTIVE ####
perc_active_stats <-  as_tibble(perc_active_proph)
perc_active_stats <- perc_active_stats %>% select(perc_temp, day_group, temperate) %>% filter(day_group %in% c(0, 182, 851, 852, 853,879, 880, 881))
perc_active_stats$day_group <- as.factor(perc_active_stats$day_group)
perc_active_stats %>% group_by(day_group) %>% get_summary_stats(temperate, type = "common")
perc_active_stats$id <- as.factor(rep(c(1,2), 8))
friedman_test(perc_active_stats, temperate ~ day_group | id)

#Add days names on x-axis
#Changle-axis labels to the day they were sampled.
#day_sampled <- c(runinfo$Time_days)
#perc_active_proph_graphp_2 <- perc_active_proph_graph + scale_x_discrete(labels= day_sampled)
perc_active_proph_graphp_2 <- perc_active_proph_graph  + 
  facet_grid(~ week, scales="free_x", space="free_x", shrink = TRUE) +
  theme(axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 7), strip.text = element_text(size=7)) + 
  geom_point(size=0.5, position = position_dodge(0.05), show.legend=FALSE) + 
  theme(legend.text = element_text(size = 7)) + 
  theme(legend.title = element_text(size = 7)) + 
  theme(axis.text.x = element_text(size = 7, angle = 270), axis.title = element_text(size = 7)) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#perc_active_proph_graphp_2  
qt <- ggplot_gtable(ggplot_build(perc_active_proph_graphp_2 ))
qt$widths[5] = 12*qt$widths[5]
grid.draw(qt)
ggsave(file="perc_active_proph_graph7.svg", plot=qt, width=7, height=11, units = c("cm"), dpi = 600)

#### DESEQ normalized read abundance Bray-Curtis ####
dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ days)
dds <- estimateSizeFactors(dds)

norDESEQ_counts <- counts(dds, normalized=TRUE)
norDESEQ_counts <- as.data.frame(norDESEQ_counts)
norDESEQ_counts <-  setDT(norDESEQ_counts, keep.rownames = 'contig')[]
norDESEQ_counts <- norDESEQ_counts[-5891,]

#Before also normalizing for contig length I will make sure both contig list and norDESEQ_counts are in same order of contigs
all(contigLength$Contig %in% norDESEQ_counts$contig)
all(contigLength$Contig == norDESEQ_counts$contig)

#Normalize counts per contig length
temp_dataframe <- norDESEQ_counts[,2:25]/contigLength$Length

norDESEQ_counts <- cbind(norDESEQ_counts$contig, temp_dataframe)
x <- colnames(norDESEQ_counts)
x[1] <- 'contig'
colnames(norDESEQ_counts) <- x

#Get relative abundance
norDESEQ_counts <- as.data.frame(norDESEQ_counts)
for (x in 2:ncol(norDESEQ_counts)){
  total_normalized_counts_per_sequence_run <- sum(norDESEQ_counts[,x])
  norDESEQ_counts[,x] <- norDESEQ_counts[,x]/total_normalized_counts_per_sequence_run
}


norDESEQ_counts <- as.data.table(norDESEQ_counts)

temperate_DESEQ <- norDESEQ_counts %>% mutate(lifestyle = case_when(str_detect(contig, "NODE") ~ 'lytic', str_detect(contig, "DAS") ~ 'temperate' ) )
temperate_DESEQ <- temperate_DESEQ %>% filter(lifestyle == 'temperate') %>% select(!(lifestyle))


#Bray-curtis difference

norDESEQ_counts_transformed <- as.data.frame(t(norDESEQ_counts))
colnames(norDESEQ_counts_transformed) <- norDESEQ_counts_transformed[1,]
norDESEQ_counts_transformed <- norDESEQ_counts_transformed[-1,]
temp_row_names <- row.names(norDESEQ_counts_transformed)
norDESEQ_counts_transformed <- sapply(norDESEQ_counts_transformed, as.numeric) #This turns lifestyle row into NA
row.names(norDESEQ_counts_transformed) <- runinfo$ID[match(temp_row_names, runinfo$SRR)]
#norDESEQ_counts_transformed <- norDESEQ_counts_transformed[1:24,1:5] #drop the NA row
beta_dist <- vegdist(norDESEQ_counts_transformed, index = "bray")
beta_dist_clust <- hclust(beta_dist, method = "average")
plot(beta_dist_clust)

adonis(beta_dist ~ week + Time_days, data= runinfo)

PCOA <- pcoa(beta_dist)
biplot.pcoa(PCOA)
beta_dist

beta_dist_matrix <- as.matrix(beta_dist)
beta_dist_matrix[lower.tri(beta_dist_matrix, diag = FALSE)] <- NA
beta_dist_matrix[beta_dist_matrix==0.0000000] <- NA

beta_dist_matrix <- as.data.frame(beta_dist_matrix)
beta_dist_matrix <- cbind(row.names(beta_dist_matrix), beta_dist_matrix)
beta_dist_matrix_long <- melt(as.data.table(beta_dist_matrix), id.vars = colnames(beta_dist_matrix)[1])
colnames(beta_dist_matrix_long) <- c('time1', 'time2', 'beta_dis')
beta_dist_matrix_long$time1 <- runinfo$Time_days[match(beta_dist_matrix_long$time1, runinfo$ID)]
beta_dist_matrix_long$time2 <- runinfo$Time_days[match(beta_dist_matrix_long$time2, runinfo$ID)]
beta_dist_matrix_long$srt_time_lag <- sqrt(abs((beta_dist_matrix_long$time1 - beta_dist_matrix_long$time2)))
#### BETADIST FIGURE ####
plot(beta_dist_matrix_long$srt_time_lag, beta_dist_matrix_long$beta_dis)
linearmodel_beta_all <- lm(beta_dist_matrix_long$beta_dis ~ beta_dist_matrix_long$srt_time_lag)
abline(a=0.1849, b=0.0245)
summary(linearmodel_beta_all)
resid(linearmodel_beta_all)

beta_dist_plot <- ggplot(beta_dist_matrix_long, aes(x=srt_time_lag, y=beta_dis)) +
  geom_point(size = 0.75) + 
  geom_abline(intercept=0.1849, slope=0.0245, linetype='dashed', color='#F8766D') + 
  ylab(label = 'Bray Curtis Dissimilarity') + 
  xlab(label = bquote(Time (Days^-2))) + 
  theme(axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 7), axis.title = element_text(size = 7 )) 

beta_dist_plot

### REVIEWERS COMMENT ####
#Transform BC into log(1-BC)
beta_dist_matrix_long$z <- log(1-beta_dist_matrix_long$beta_dis)
plot(beta_dist_matrix_long$srt_time_lag, beta_dist_matrix_long$z)
reviewer_lm_beta_all <- lm(beta_dist_matrix_long$z ~beta_dist_matrix_long$srt_time_lag)
reviewer_lm_beta_all <- lm(beta_dist_matrix_long$z * beta_dist_matrix_long$srt_time_lag)
summary(reviewer_lm_beta_all)
resid(reviewer_lm_beta_all)
x <- beta_dist_matrix_long$srt_time_lag
newfunc2 <- function(x,a){x*a}
newFunc <- function(x,a,b){b-exp(a*x)}
curve(newfunc2(x,-0.075), add=TRUE)
linearmodel_beta_all <- lm(beta_dist_matrix_long$z ~ beta_dist_matrix_long$srt_time_lag)
abline(a=-0.0649, b=-0.075)


beta_dist_plot <- ggplot(beta_dist_matrix_long, aes(x=srt_time_lag, y=z)) +
  geom_point(size = 0.75) + 
  geom_abline(intercept=-0.0649, slope=-0.075, linetype='dashed', color='#F8766D') + 
  ylab(label = 'Bray Curtis Dissimilarity') + 
  xlab(label = bquote(Time (Days^-2))) + 
  theme(axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 7), axis.title = element_text(size = 7 )) 


#Missing points are simply because I am made a triangle matrix, and did so by making one side NAs, there is probably a better way
ggsave(file="beta_dist_plot3_reviewers.png", plot=beta_dist_plot, width=5, height=6, units = c("cm"), dpi = 600, device = 'png')

#Look at prophages

prophage_norDESEQ_counts_transformed <- norDESEQ_counts_transformed[,1:52]

prophage_beta_dist <- vegdist(prophage_norDESEQ_counts_transformed, index = "bray")
prophage_beta_dist_clust <- hclust(prophage_beta_dist, method = "average")
plot(prophage_beta_dist_clust)
adonis(prophage_beta_dist ~ week + Time_days, data= runinfo)
prophage_PCOA <- pcoa(prophage_beta_dist)
biplot.pcoa(prophage_PCOA)

prophage_beta_dist_matrix <- as.matrix(prophage_beta_dist)
prophage_beta_dist_matrix[lower.tri(prophage_beta_dist_matrix, diag = FALSE)] <- NA
prophage_beta_dist_matrix[prophage_beta_dist_matrix==0.0000000] <- NA

prophage_beta_dist_matrix <- as.data.frame(prophage_beta_dist_matrix)
prophage_beta_dist_matrix <- cbind(row.names(prophage_beta_dist_matrix), prophage_beta_dist_matrix)
prophage_beta_dist_matrix_long <- melt(as.data.table(prophage_beta_dist_matrix), id.vars = colnames(prophage_beta_dist_matrix)[1])
colnames(prophage_beta_dist_matrix_long) <- c('time1', 'time2', 'beta_dis')
prophage_beta_dist_matrix_long$time1 <- runinfo$Time_days[match(prophage_beta_dist_matrix_long$time1, runinfo$ID)]
prophage_beta_dist_matrix_long$time2 <- runinfo$Time_days[match(prophage_beta_dist_matrix_long$time2, runinfo$ID)]
prophage_beta_dist_matrix_long$srt_time_lag <- sqrt(abs((prophage_beta_dist_matrix_long$time1 - prophage_beta_dist_matrix_long$time2)))

plot(prophage_beta_dist_matrix_long$srt_time_lag, prophage_beta_dist_matrix_long$beta_dis)
linearmodel_beta_prophages <- lm(prophage_beta_dist_matrix_long$beta_dis ~ prophage_beta_dist_matrix_long$srt_time_lag)
abline(0.584793, 0.003089)
summary(linearmodel_beta_prophages)
prophage_beta_dist_plot <- ggplot(prophage_beta_dist_matrix_long, aes(x=srt_time_lag, y=beta_dis)) + geom_point(size = 0.75) + 
  geom_abline(intercept=0.58479, slope=0.003089, linetype='dashed', color='#F8766D') + 
  ylab(label = 'Bray Curtis Dissimilarity') + 
  xlab(label = bquote(Time (Days^-2))) + 
  theme(axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 7), axis.title = element_text(size = 7 )) 
prophage_beta_dist_plot
ggsave(file="prophage_beta_dist_plot3.svg", plot=prophage_beta_dist_plot, width=5, height=6, units = c("cm"), dpi = 600)

#### REVIEWERS COMMENT ####

#Transform BC into log(1-BC)
prophage_beta_dist_matrix_long$z <- log(1-prophage_beta_dist_matrix_long$beta_dis)
plot(prophage_beta_dist_matrix_long$srt_time_lag, prophage_beta_dist_matrix_long$z)
reviewer_lm__prophage_beta_all <- lm(prophage_beta_dist_matrix_long$z ~prophage_beta_dist_matrix_long$srt_time_lag)

summary(reviewer_lm__prophage_beta_all)
resid(reviewer_lm__prophage_beta_all)
x <- prophage_beta_dist_matrix_long$srt_time_lag
newfunc2 <- function(x,a){x*a}
newFunc <- function(x,a,b){b-exp(a*x)}
curve(newfunc2(x,-0.0075), add=TRUE)
abline(a=-1.037166, b=-0.008292)

beta_dist_plot <- ggplot(prophage_beta_dist_matrix_long, aes(x=srt_time_lag, y=z)) +
  geom_point(size = 0.75) + 
  geom_abline(intercept=-1.037166, slope=-0.008292, linetype='dashed', color='#F8766D') + 
  ylab(label = 'Bray Curtis Dissimilarity') + 
  xlab(label = bquote(Time (Days^-2))) + 
  theme(axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 7), axis.title = element_text(size = 7 )) 


#Missing points are simply because I am made a triangle matrix, and did so by making one side NAs, there is probably a better way
ggsave(file="prophage_beta_dist_plot3_reviewers.png", plot=beta_dist_plot, width=5, height=6, units = c("cm"), dpi = 600, device = 'png')


#Now I want to compare beta_dist_matrix_long (all phages) to prophage_beta_dist_matrix_long (prophages)
comaparing_beta_dist_matrix_long <- beta_dist_matrix_long
names(comaparing_beta_dist_matrix_long)[names(comaparing_beta_dist_matrix_long) == 'beta_dis'] <- 'beta_dis_all'
comaparing_beta_dist_matrix_long$beta_dis_prophage <- prophage_beta_dist_matrix_long$beta_dis
compared_model <- lm(beta_dis_all ~ srt_time_lag + beta_dis_prophage, data = comaparing_beta_dist_matrix_long)
summary(compared_model)
anova(compared_model)

linearHypothesis(compared_model, "beta_dis_all = beta_dis_prophage")

#Next I want to see if the high beta-diversity over time is due to the larger number of contigs
random_contigs <- floor(runif(59, min=0, max=5890))
random_norDESEQ_counts_transformed <-  norDESEQ_counts_transformed[,random_contigs]
random_beta_dist <- vegdist(random_norDESEQ_counts_transformed, index = "bray")

random_beta_dist_matrix <- as.matrix(random_beta_dist)
random_beta_dist_matrix[lower.tri(random_beta_dist_matrix, diag = FALSE)] <- NA
random_beta_dist_matrix[random_beta_dist_matrix==0.0000000] <- NA

random_beta_dist_matrix <- as.data.frame(random_beta_dist_matrix)
random_beta_dist_matrix <- cbind(row.names(random_beta_dist_matrix), random_beta_dist_matrix)
random_beta_dist_matrix_long <- melt(random_beta_dist_matrix, id.vars = colnames(random_beta_dist_matrix)[1])

colnames(random_beta_dist_matrix_long) <- c('time1', 'time2', 'beta_dis')
random_beta_dist_matrix_long$time1 <- runinfo$Time_days[match(random_beta_dist_matrix_long$time1, runinfo$ID)]
random_beta_dist_matrix_long$time2 <- runinfo$Time_days[match(random_beta_dist_matrix_long$time2, runinfo$ID)]
random_beta_dist_matrix_long$srt_time_lag <- sqrt(abs((random_beta_dist_matrix_long$time1 - random_beta_dist_matrix_long$time2)))

plot(random_beta_dist_matrix_long$srt_time_lag, random_beta_dist_matrix_long$beta_dis)
linearmodel_beta_subsample <- lm(random_beta_dist_matrix_long$beta_dis ~ random_beta_dist_matrix_long$srt_time_lag)
#abline(a=0.1849, b=0.0245)
summary(linearmodel_beta_subsample)


##### DESEQ taxnomic abundance #####
norDESEQ_counts <- read.table('../norDESEQ_counts_norbylength.tsv', sep = "\t", header =TRUE)

#I want relative abundance so I will divide each value in column

for (column in 2:ncol(norDESEQ_counts)){
  total_normalized_counts_per_sequence_run <- sum(norDESEQ_counts[,column])
  norDESEQ_counts[,column] <- norDESEQ_counts[,column]/total_normalized_counts_per_sequence_run
}


norDESEQ_melt <- melt(as.data.table(norDESEQ_counts), id.vars = c("contig"))

#Add in the Demovir taxonomy for phages that have it

norDESEQ_melt$Order <- demovir$Order[match(norDESEQ_melt$contig, demovir$Sequence_ID)]
norDESEQ_melt$Family <- demovir$Family[match(norDESEQ_melt$contig, demovir$Sequence_ID)]
#Make NAs "Unknown"
norDESEQ_melt$Order <- as.character(norDESEQ_melt$Order)
norDESEQ_melt$Order[is.na(norDESEQ_melt$Order)] <- "Unknown"

norDESEQ_melt$Family <- as.character(norDESEQ_melt$Family)
norDESEQ_melt$Family[is.na(norDESEQ_melt$Family)] <- "Unknown"

#I still prefer to have samples order by sequence run ID

norDESEQ_melt$day <- runinfo$ID[match(norDESEQ_melt$variable, runinfo$SRR)]
norDESEQ_melt$week <- runinfo$week[match(norDESEQ_melt$variable, runinfo$SRR)]

norDESEQ_melt$Family <- as.factor(norDESEQ_melt$Family)
#Making a side-step so that I can later on look at the number of unclassified phages
contig_with_taxonomy <- norDESEQ_melt

norDESEQ_family_melt <- norDESEQ_melt %>% group_by(day, Family, week) %>% summarise(value = sum(value))

family_relaytive_abund <- ggplot(norDESEQ_family_melt, aes(fill=Family, y=value, x=day)) + 
  geom_bar(position = "stack", stat="identity", colour=NA) + 
  theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust=1))  
family_relaytive_abund <- family_relaytive_abund +  
  facet_wrap(~ week, scales="free_x") + ylab(label = 'Relative Abundance by Order') + 
  xlab(label = 'Day') + 
  theme(strip.text = element_text(size =20)) + 
  theme(axis.title = element_text(size = 30))

family_relaytive_abund
ggsave(file="family_relaytive_abund.svg", plot=family_relaytive_abund, width=10, height=10)


norDESEQ_melt %>% filter(Family == "Unknown") -> unknown_phages
norDESEQ_melt %>% filter(Family == "Unknown") %>% filter(value >= 1e-1) -> abundant_unknown_phages

#Make a graph for %Known vs %Unknown
w <- norDESEQ_melt
w$Family <-  as.character(w$Family)
w$Family[w$Family != 'Unknown'] <- "Known"
#Another side table to see how many phages are known vs unknown
contig_with_taxonomy2 <- as_tibble(w)

w <-  as_tibble(w %>% group_by(variable, Family) %>% summarise((sum(value))))
w$day <- runinfo$ID[match(w$variable, runinfo$SRR)]
w$week <- runinfo$week[match(w$variable, runinfo$SRR)]
names(w)[names(w) == colnames(w[3])] <- "value"

perc_known <- ggplot(w, aes(fill=Family, y=value, x=day)) + geom_bar(position = "stack", stat="identity") #+ theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust=1))  
perc_known <- perc_known + 
  #facet_wrap(~ week, scales="free_x") + 
  ylab(label = 'Relative Abundance') + xlab(label = 'Day') + 
  theme(strip.text = element_text(size =20)) + theme(axis.title = element_text(size = 30))

#Remove background ggplot defeault stuff
perc_known <- perc_known + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                                              panel.grid.minor = element_blank())


#Relabel days, as days, and reformat the look
day_sampled <- c(runinfo$Time_days)
perc_known2 <- perc_known + scale_x_discrete(labels= day_sampled) + xlab(label = 'Day') + theme(axis.text.x = element_text(size = 25), axis.text.y = element_text(size = 25)) + theme(strip.text = element_text(size =35)) + theme(axis.title = element_text(size = 35))
perc_known2
ggsave(file="perc_known2.svg", plot=perc_known2, width=15, height=10)

#I will make a version that shows a line per contig 
contig_with_taxonomy2$day <- runinfo$ID[match(contig_with_taxonomy2$variable, runinfo$SRR)]
contig_with_taxonomy2$week <- runinfo$week[match(contig_with_taxonomy2$variable, runinfo$SRR)]

perc_known3 <- ggplot(contig_with_taxonomy2, aes(fill=Family, y=value, x=day)) + geom_bar(position = "stack", stat="identity") #+ theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust=1))  
perc_known3 <- perc_known3 +  
  #facet_wrap(~ week, scales="free_x") + 
  ylab(label = 'Relative Abundance') + xlab(label = 'Day') + 
  theme(strip.text = element_text(size =20)) + 
  theme(axis.title = element_text(size = 30))

#Remove background ggplot defeault stuff
perc_known3 <- perc_known3 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                              panel.grid.minor = element_blank())

perc_known3
#Relabel days, as days, and reformat the look
day_sampled <- c(runinfo$Time_days)
perc_known4 <- perc_known3 + 
  scale_x_discrete(labels= day_sampled) + 
  xlab(label = 'Day') + 
  theme(axis.text.x = element_text(size = 7, angle = 270), axis.text.y = element_text(size = 7)) + 
  theme(strip.text = element_text(size =7)) + theme(axis.title = element_text(size = 7)) +
  theme(legend.key.size = unit(0.25, 'cm'), #change legend key size
        legend.key.height = unit(0.25, 'cm'), #change legend key height
        legend.key.width = unit(0.25, 'cm'), #change legend key width
        legend.title = element_text(size=7), #change legend title font size
        legend.text = element_text(size=7)) #change legend text font size
perc_known4

ggsave(file="perc_known_w_lines2.svg", plot=perc_known4, width=11, height=11, units = c("cm"), dpi = 600)



#I will now try run the analysis without the Unknown Family/Order category to remove non-demovir hits bias in week 2 / week 3

norDESEQ_melt_known_phages <- norDESEQ_melt %>% filter(Family != "Unknown")
known_phages <- unique(norDESEQ_melt_known_phages$contig)
known_norDESEQ <- norDESEQ_counts[match(known_phages, norDESEQ_counts$contig), ]


for (column in 2:ncol(known_norDESEQ)){
  total_normalized_counts_per_sequence_run <- sum(known_norDESEQ[,column])
  known_norDESEQ[,column] <- known_norDESEQ[,column]/total_normalized_counts_per_sequence_run
}

known_norDESEQ_melt <- melt(as.data.table(known_norDESEQ), id.vars = c("contig"))

#Add in the Demovir taxonomy for phages that have it

known_norDESEQ_melt$Order <- demovir$Order[match(known_norDESEQ_melt$contig, demovir$Sequence_ID)]
known_norDESEQ_melt$Family <- demovir$Family[match(known_norDESEQ_melt$contig, demovir$Sequence_ID)]
#Make NAs "Unknown"
known_norDESEQ_melt$Order <- as.character(known_norDESEQ_melt$Order)
known_norDESEQ_melt$Order[is.na(known_norDESEQ_melt$Order)] <- "Unknown"

known_norDESEQ_melt$Family <- as.character(known_norDESEQ_melt$Family)
known_norDESEQ_melt$Family[is.na(known_norDESEQ_melt$Family)] <- "Unknown"

#### Reviewers Comment ####

known_norDESEQ_melt <- known_norDESEQ_melt %>% filter(Family == "Myoviridae" | Family == "Siphoviridae" | Family ==  "Unassigned" | Family == "Podoviridae" | Family =="Microviridae")
Family_Colours <- c("Myoviridae" = "#EEA236", "Siphoviridae" = "#46B8DA", "Unassigned" = "#B8B8B8" , "Podoviridae" = "#357EBD", Microviridae = "#5CB85C")
#I still prefer to have samples order by sequence run ID

known_norDESEQ_melt$day <- runinfo$ID[match(known_norDESEQ_melt$variable, runinfo$SRR)]
known_norDESEQ_melt$week <- runinfo$week[match(known_norDESEQ_melt$variable, runinfo$SRR)]


#Relative abundance by family is fragmented into bars for each contig. So I will sum them by family
known_norDESEQ_melt2 <-  as_tibble(known_norDESEQ_melt %>% group_by(variable, Family) %>% summarise((sum(value))))
known_norDESEQ_melt2$day <- runinfo$ID[match(known_norDESEQ_melt2$variable, runinfo$SRR)]
known_norDESEQ_melt2$week <- runinfo$week[match(known_norDESEQ_melt2$variable, runinfo$SRR)]
known_norDESEQ_melt2$day_group <- runinfo$Time_days[match(known_norDESEQ_melt2$variable, runinfo$SRR)]
names(known_norDESEQ_melt2)[names(known_norDESEQ_melt2) == colnames(known_norDESEQ_melt2[3])] <- "value"
known_norDESEQ_melt2$day <- factor(known_norDESEQ_melt2$day, levels= runinfo$ID)


day_sampled <- c(runinfo$Time_days)
ID_sampled <- c(runinfo$ID)
knownfamily_relative_abund <- ggplot(known_norDESEQ_melt2, aes(fill=Family, y=value, x=day)) + 
  geom_bar(position = "stack", stat="identity", show.legend=TRUE) + scale_x_discrete(breaks = ID_sampled, labels = day_sampled) + scale_fill_manual(values = Family_Colours)  #Add in when generating the no label figure
knownfamily_relative_abund <- knownfamily_relative_abund +  
  #facet_wrap(~ week, scales="free_x") + 
  ylab(label = 'Relative Abundance by Order') + 
  xlab(label = 'Day') + theme(strip.text = element_text(size =20)) + 
  theme(axis.title = element_text(size = 30))

#Remove background ggplot defeault stuff
knownfamily_relative_abund <- knownfamily_relative_abund + 
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
#Relabel days, as days, and reformat the look

knownfamily_relative_abund2 <- knownfamily_relative_abund + 
  xlab(label = 'Day') + 
  theme(axis.text.x = element_text(size = 25), axis.text.y = element_text(size = 25)) + 
  theme(strip.text = element_text(size =35)) + 
  theme(axis.title = element_text(size = 35))
knownfamily_relative_abund2 <- knownfamily_relative_abund2 + 
  facet_grid(~ week, scales="free_x", space="free_x", shrink = TRUE) +
  theme(axis.text.x = element_text(size = 7, angle = 270), axis.text.y = element_text(size = 7), strip.text = element_text(size=7)) + 
  theme(legend.title = element_text(size = 7)) + 
  theme(axis.text.x = element_text(size = 7), axis.title = element_text(size = 7)) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.text = element_text(size = 7)) +
  theme(legend.key.size = unit(0.25, 'cm'), #change legend key size
        legend.key.height = unit(0.25, 'cm'), #change legend key height
        legend.key.width = unit(0.25, 'cm'), #change legend key width
        legend.title = element_text(size=7), #change legend title font size
        legend.text = element_text(size=7)) #change legend text font size
knownfamily_relative_abund2
ggsave(file="knownfamily_relative_abund5.svg", plot=knownfamily_relative_abund2, width=17, height=13, units = c("cm"), dpi = 600)
ggsave(file="knownfamily_relative_abund_nolegend5.svg", plot=knownfamily_relative_abund2, width=8, height=10, units = c("cm"), dpi = 600)

#Make a figure for known vs unkown
norDESEQ_family_melt

# # Alejandro's request to see a figure of how many contigs belong to known vs unknown, to see if it is one or two highly abundant
# 
# contig_with_taxonomy2 <- as_tibble(w)
# counts_per_category <- contig_with_taxonomy2 %>% group_by(variable, Family) %>% summarise(contig_number = n())
# 
# w <-  as_tibble(w %>% group_by(variable, Family) %>% summarise((sum(value))))
# w$day <- runinfo$ID[match(w$variable, runinfo$SRR)]
# w$week <- runinfo$week[match(w$variable, runinfo$SRR)]
# names(w)[names(w) == colnames(w[3])] <- "value"
# 
# perc_known <- ggplot(w, aes(fill=Family, y=value, x=day)) + geom_bar(position = "stack", stat="identity") #+ theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust=1))  
# perc_known <- perc_known +  facet_wrap(~ week, scales="free_x") + ylab(label = 'Relative Abundance') + xlab(label = 'Day') + theme(strip.text = element_text(size =20)) + theme(axis.title = element_text(size = 30))
# 
# #Remove background ggplot defeault stuff
# perc_known <- perc_known + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                                               panel.grid.minor = element_blank())
# #Relabel days, as days, and reformat the look
# day_sampled <- c(runinfo$Time_days)
# perc_known2 <- perc_known + scale_x_discrete(labels= day_sampled) + xlab(label = 'Day') + theme(axis.text.x = element_text(size = 25), axis.text.y = element_text(size = 25)) + theme(strip.text = element_text(size =35)) + theme(axis.title = element_text(size = 35))
