#Open packages
library(dplyr) #Version ‘1.0.9’
library(ggplot2) #Version ‘3.3.6’
library(tidyr) #Version ‘1.2.0’
library(stringr) #Version ‘1.4.0’
library(reshape2) #Version ‘1.4.4’
library(cartography) #Version ‘3.0.1’
library(grid) #Version ‘4.2.1’
library(ggrepel) #Version ‘0.9.1’

#Each sequence run information
setwd('Data_For_Scripts/Prophage_Coverage/')
SRR <- read.csv("../sequence_runs.csv", header=TRUE)
contig_bed <- read.delim("merged_phage_contigs.bed", header=FALSE, blank.lines.skip = TRUE)
contig_bed <- contig_bed[,1:3]
colnames(contig_bed) <- c('contig', 'start', 'length')

setwd('samtools_coverage_coverage')
temp <-  list.files(pattern="*_samtool_coverage") #Make a temp variable for each of the readcount files
myfiles <-  lapply(temp, function(x) read.delim(x, header = TRUE)) #Opens up each file
coverage_data <- do.call(rbind, myfiles) #Binds all the files together in one dataframe
setwd('..')
runinfo <- read.csv("../sequence_runs.csv", header = TRUE)
bacterial_bin <- read.csv(file = '../GTDB-Tk_DAS_check_bin_modified.csv', header = FALSE)
DESEQ_size_factors <- read.csv("../DESEQ_sizefactors.csv", header = TRUE)
#An important file this matches each prophage contig name to a prophage id based on bacterial bin
prophage_host_match <- read.csv(file = "../prophage_host_naming.csv", header = TRUE)

#I need to make a column with the SRR it is from. 
SRR_names <- substring(temp,1,9) #Make a vector of the SRR run by name in order files of temp.
#E.g this replicate(5,SRR_names[1]) would replicate the first sequence run 5 times but I need it 5890 times
SRR_columns <- c()
for (i in 1:24)
{
  SRR_columns <- c(SRR_columns,(replicate(5890,SRR_names[i])))
}

runinfo$permillionreads <- runinfo$totalNumReads/1000000

#Now I can add it into my coverage data the SRR run
coverage_data <- cbind(coverage_data, SRR_columns)

#Change name of contig column
x <- colnames(coverage_data)
x[1] <- "contig"
colnames(coverage_data) <- x

#Make a weak cut-off for present (i.e any reads mapping to contig)
coverage_data$minot_present <- (coverage_data$coverage > 0)

#Make a column for prophage or not
coverage_data <- coverage_data %>% mutate(lifestyle = case_when(str_detect(contig, "NODE") ~ 'lytic', str_detect(contig, "DAS") ~ 'temperate' ) )

coverage_data$present <- (coverage_data$coverage >= 75) #First results I did this with 25 because I thought it was coverage of 0s not other way around

present_contigs <- coverage_data %>% filter(present == TRUE)
present_contigs <- unique(present_contigs$contig)
length(present_contigs)

#Check prophages 

present_temperate <- coverage_data %>% filter(present == TRUE)
present_temperate <- present_temperate %>% filter(lifestyle == 'temperate')
length(unique(present_temperate$contig))

#Make presence absence heatmap
present_prophages <- coverage_data %>% filter(lifestyle == 'temperate')
coverage_data$days <- runinfo$ID[(match(coverage_data$SRR_columns, runinfo$SRR))]
present_prophages$SRR_columns <- factor(present_prophages$SRR_columns, levels = c(runinfo$SRR))
#This lumps the days together but some days are sequenced twice so I will try it (below) with the day-sequence run
#present_prophages$day <- as.factor(runinfo$Time_days[match(present_prophages$SRR_columns, runinfo$SRR)])
present_prophages$day <- as.factor(runinfo$ID[match(present_prophages$SRR_columns, runinfo$SRR)])
present_prophages$contig <- prophage_host_match$prophage_names[match(present_prophages$contig, prophage_host_match$prophage_ID)]
prophage_presence.heatmap <- ggplot(data = present_prophages, aes(x = day, y = contig, fill = present)) +
  geom_tile() + scale_fill_manual(values = c('white', "#CE3D32FF")) +
  theme(axis.text.x = element_text(size = 5, angle = 90, ), axis.text.y = element_text(size = 5, face = "italic")) +
  xlab(label = "Day") + ylab(label = 'Active Prophage') + theme(axis.text.x = element_text(size = 10))
prophage_presence.heatmap

# I want to get some metrics, like how many prophages were found over multiple weeks, all time points, etc.

present_temperate_sorted <- present_temperate[ , c("contig","SRR_columns")]
present_temperate_sorted$week <- runinfo$week[match(present_temperate_sorted$SRR_columns, runinfo$SRR)]
present_temperate_sorted$day <- runinfo$Time_days[match(present_temperate_sorted$SRR_columns, runinfo$SRR)]

active_prophage_list <- unique(present_temperate_sorted$contig)
metrics_present_temperate <- c()

for (i in active_prophage_list){
  temp_matrix <- present_temperate_sorted %>% filter(contig == i)
  weeks <- length(unique(temp_matrix$week))
  days <- length(unique(temp_matrix$day))
  temp_variable <- c(i, weeks, days)
  metrics_present_temperate <- rbind(metrics_present_temperate, temp_variable)
}
colnames(metrics_present_temperate) <- c('contig', 'Number_of_Weeks' ,"Number_of_Days")
metrics_present_temperate <- as.data.frame(metrics_present_temperate)
metrics_present_temperate$Number_of_Weeks <- as.numeric(metrics_present_temperate$Number_of_Weeks)
metrics_present_temperate$Number_of_Days <- as.numeric(metrics_present_temperate$Number_of_Days)
summary(metrics_present_temperate)

#Lets see how non-prophages

#Stop phase, output the 75% coverage table
write.table(coverage_data, 'prophage_coverage_data.tsv', sep = "\t")
coverage_data <- read.table('prophage_coverage_data.tsv', sep = "\t")

#replace the SRR with days
coverage_data$days <- SRR$ID[(match(coverage_data$SRR_columns, SRR$SRR))]





#To compare the coverages of each prophage between dates to see if prophage induction has occured

prophage_induction <- present_prophages[,c('contig','meandepth','SRR_columns')]
prophage_induction <- spread(prophage_induction, contig, meandepth)

##### DESEQ SIZE FACTOR ####
#Note run at least the first 201 lines to get here

runinfo$size_factors <-  DESEQ_size_factors$sizeFactors[(match(runinfo$SRR,DESEQ_size_factors$SRR ))]

normalized_coverage2 <- c()
for (row in 1:nrow(prophage_induction))
{
  x <- (prophage_induction[row,-1]/runinfo$size_factors[row])
  normalized_coverage2 <- rbind(normalized_coverage2, x)
}

normalized_coverage2 <- cbind(runinfo$Time_days, runinfo$week, runinfo$ID, normalized_coverage2)
x <- colnames(normalized_coverage2)
x[1] <- 'day'
x[2] <-  'week'
x[3] <- 'day-replicate'
colnames(normalized_coverage2) <- x
normalized_coverage2$`day-replicate` <- factor(normalized_coverage2$`day-replicate`)
normalized_coverage2$sample_num <- c(1:24)
normalized_coverage2_melted <- melt(normalized_coverage2, id.vars= c('day', 'week', 'day-replicate', 'sample_num'))

normalized_coverage3_melted <- normalized_coverage2_melted
#normalized_coverage3_melted$variable <- prophage_host_match$prophage_names[match(normalized_coverage3_melted$variable, prophage_host_match$prophage_ID)]
normalized_coverage3_melted$`day-replicate` <- as.factor(runinfo$ID)

#Add bacterial names of prophages
normalized_coverage3_melted$phage_name <- prophage_host_match$prophage_names[match(normalized_coverage3_melted$variable, prophage_host_match$prophage_ID)]

#I will try and label all those above 100x coverage
p <- ggplot(normalized_coverage3_melted, aes(x=day, y=value, col=variable)) + geom_point(size=1, position = position_dodge(0.7)) + theme(axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20)) + geom_text(aes(label=ifelse(value>100,as.character(variable),'')),hjust=0,vjust=0, size=4, show.legend=FALSE)
p <- p + facet_grid(~ week, scales="free_x", space="free") + ylab(label = 'Normalized Coverage') + theme(strip.text = element_text(size =20)) + theme(axis.title = element_text(size = 30)) + theme(legend.text =element_text(face = "italic") )
p 

#Now I will try with z-score that is significant
nozero_coverage <- subset(normalized_coverage2_melted, value>0)
nozero_coverage$logCov <- log(nozero_coverage$value)
nozero_coverage$zscore <- scale(nozero_coverage$logCov)
#nozero_coverage$variable <- prophage_host_match$prophage_names[match(nozero_coverage$variable, prophage_host_match$prophage_ID)]
one_hundred_coverage <- subset(nozero_coverage, value>100)
significant_coverage <- subset(nozero_coverage, zscore>1.96)

normalized_coverage_w_zscore <- left_join(normalized_coverage3_melted, nozero_coverage[,c(1,2,3,4,5,7,8)], by =c('day', "week", "day-replicate", "sample_num", "variable"))
q <- ggplot(normalized_coverage_w_zscore, aes(x=day, y=value, col=variable)) + 
  geom_point(size=1, position = position_dodge(0.1), show.legend=FALSE) + 
  theme(axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 7)) + 
  geom_text(aes(label=ifelse(zscore>1.95,as.character(variable),'')),hjust=0,vjust=0, size=2, show.legend=FALSE)
q <- q + facet_grid(~ week, scales="free_x", space="free_x", shrink = TRUE) + ylab(label = 'Normalized Coverage') + theme(strip.text = element_text(size =7)) + theme(axis.title = element_text(size = 7)) + theme(legend.text =element_text(face = "italic") )
q <- q + theme(legend.text = element_text(size = 7)) + theme(legend.title = element_text(size = 7)) + theme(axis.text.x = element_text(size = 7), axis.title = element_text(size = 7)) + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
q
qt <- ggplot_gtable(ggplot_build(q))
qt$widths[5] = 12*qt$widths[5]
grid.draw(qt)

ggsave(file="prophage_coverage_DESEQ_norm6.svg", plot=q,  width=30, height=13, units = c("cm"), dpi = 600)

#Note, so I can't figure out how make a day 0 which has just one time point, look like the other weeks so I am going to make a fake dataset where
#During week has other days 1,2,3,4 so that I can later remove them in Inkskape in a way that the data looks the same. AFter trying different things
#With facet_grid or ggplot_gtable

day1 <- filter(normalized_coverage_w_zscore, week == 'day0')
day1$day <- 1

day2 <- filter(normalized_coverage_w_zscore, week == 'day0')
day2$day <- 2

day3 <- filter(normalized_coverage_w_zscore, week == 'day0')
day3$day <- 3

day4 <- filter(normalized_coverage_w_zscore, week == 'day0')
day4$day <- 4

normalized_coverage_w_zscore2 <- rbind(normalized_coverage_w_zscore, day1, day2, day3, day4)

qt <- ggplot(normalized_coverage_w_zscore2, aes(x=day, y=value, col=variable)) + 
  geom_point(size=1, position = position_dodge(0.1), show.legend=FALSE) + 
  theme(axis.text.x = element_text(size = 7, angle = 270), axis.text.y = element_text(size = 7)) + 
  geom_text_repel(aes(label=ifelse(zscore>1.95,as.character(variable),'')),hjust=0,vjust=0, size=2, show.legend=FALSE)
qt <- qt + facet_grid(~ week, scales="free_x", space="free_x", shrink = TRUE) + ylab(label = 'Normalized Coverage') + theme(strip.text = element_text(size =7)) + theme(axis.title = element_text(size = 7)) + theme(legend.text =element_text(face = "italic") )
qt <- qt + theme(legend.text = element_text(size = 7)) + theme(legend.title = element_text(size = 7)) + theme(axis.text.x = element_text(size = 7), axis.title = element_text(size = 7)) + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
qt
ggsave(file="prophage_coverage_DESEQ_norm7.svg", plot=qt,  width=8, height=11, units = c("cm"), dpi = 600)


#Make a heatmap with normalized coverage
normalized_coverage4 <- c()
for (row in 1:nrow(present_prophages)){
  if (present_prophages$coverage[row] < 75){
    coverage = 0
  } else {
    x <- match(present_prophages$SRR_columns[row], runinfo$SRR)
    coverage = present_prophages$meandepth[row]/runinfo$size_factors[x]
  }
  y <- c(present_prophages$contig[row], as.character(present_prophages$day[row]), as.numeric(coverage))
  normalized_coverage4 <- rbind(normalized_coverage4, y)
}

normalized_coverage4 <- as.data.frame(normalized_coverage4)
colnames(normalized_coverage4) <- c("contig", "day", "coverage") 
normalized_coverage4$day <- as.factor(present_prophages$day)
normalized_coverage4$coverage <- as.numeric(normalized_coverage4$coverage)
normalized_coverage4[normalized_coverage4 == 0] <- NA
normalized_coverage4$week <- runinfo$week[(match(normalized_coverage4$day, runinfo$ID))]
normalized_coverage4$number <- prophage_host_match$prophage_numbers[(match(normalized_coverage4$contig, prophage_host_match$prophage_names))]


prophage_coverage.heatmap <- ggplot(data = normalized_coverage4, aes(x = day, y = contig, fill = coverage)) +
  geom_tile() + scale_fill_gradientn(colours = carto.pal(pal1 = "orange.pal" ,n1 = 5), na.value = 'white') + 
  theme(axis.text.x = element_text(size = 7 ), axis.text.y = element_text(size = 7, face = "italic", )) +
  xlab(label = "Day") + ylab(label = 'Active Prophages') + theme(axis.text.x = element_text(size = 7), axis.title = element_text(size = 7))
prophage_coverage.heatmap_grided <- prophage_coverage.heatmap +facet_grid(~ week) + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
prophage_coverage.heatmap_grided
prophage_coverage.heatmap 

#Changle-axis labels to the day they were sampled.
day_sampled <- c(runinfo$Time_days)
prophage_coverage.heatmap_2 <- prophage_coverage.heatmap + scale_x_discrete(labels= day_sampled) + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
prophage_coverage.heatmap_2 <- prophage_coverage.heatmap_2 + theme(legend.text = element_text(size = 7)) + theme(legend.title = element_text(size = 7))  
prophage_coverage.heatmap_2 <- prophage_coverage.heatmap_2 + theme(legend.key.size = unit(2, 'mm'), #change legend key size
                                  legend.key.height = unit(2, 'mm'), #change legend key height
                                  legend.key.width = unit(2, 'mm'), #change legend key width
                                  legend.title = element_text(size=7), #change legend title font size
                                  legend.text = element_text(size=7)) #change legend text font size
prophage_coverage.heatmap_2 
ggsave(file='prophage_coverage.heatmap_3.svg', plot=prophage_coverage.heatmap_2, width=21, height=13, units = c("cm"), dpi = 600)
ggsave(file='prophage_coverage.heatmap_grided.svg', plot=prophage_coverage.heatmap_grided, width=13, height=8)

##### Check Induced Prophages ####
#Important the prophages that went though a signficant fold increase by DESEQ analysis
#Method1 use the significantly increased DESEQ prophages
DESEQ_induced_prophages <- read.csv( file = '../DESEQ_induced_prophages.tsv', sep = "\t")
DESEQ_induced_prophages$Phage <- prophage_host_match$prophage_names[match(DESEQ_induced_prophages$Phage, prophage_host_match$prophage_ID)]

#Method2 use the normalized coverage with a z-score over 1.96
nozero_coverage <- subset(normalized_coverage2_melted, value>0)
nozero_coverage$logCov <- log(nozero_coverage$value)
nozero_coverage$zscore <- scale(nozero_coverage$logCov)
#nozero_coverage$variable <- prophage_host_match$prophage_names[match(nozero_coverage$variable, prophage_host_match$prophage_ID)]
one_hundred_coverage <- subset(nozero_coverage, value>100)
significant_coverage <- subset(nozero_coverage, zscore>1.96)

# bins_with_induced <- unique((significant_coverage$variable))
# bins_position <- match(bins_with_induced, bacterial_bin$V1)
# bin_names <- bacterial_bin$V9[bins_position]
# bin_names

#Compare the different methods
method1 <- DESEQ_induced_prophages$Phage
method2 <- significant_coverage$variable
intersect(method1,method2)

##### Missclassified Prophage-Host #####

coverage_data %>% filter(str_detect(contig, "DAS")) %>% filter(coverage > 95) %>% select(contig) -> high_covered_prophages
unique(sort(high_covered_prophages[,]))

#We see that 28 of 52 prophages could fall into this 'miss-classified to host' area. 

