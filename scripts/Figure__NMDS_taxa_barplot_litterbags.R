# (A) NMDS, Grass samples for BACTERIA only & 
# (B) taxa barplot over time underneath (for ____ treatments)
# Kendra E. Walters
# April 9th 2020

# Room of ()ments
require(data.table)
require(ggplot2)
require(vegan)
require(anchors)


#####__________ For Panel A = NMDS of grass samples__________######
# Setting up and reading in
metadata <- as.data.frame(fread("metadata_q2.txt"))
table <- as.data.frame(fread("14_table_filtered.tsv"))

###__ CLEANING THE OTU TABLE________
# Filter out mock, negative, and replicates by name
mock_and_neg <- c("Negative3", "PCR Mock1", "PCR Mock3", "PCR Mock4", "PCR mock5", 
                  "PCR Neg3", "PCR Neg4", "PCR neg5", "RCT4R1 (Rep1)", "RCT4R1 (Rep 2)", 
                  "LitterT4R2 (Rep 1)", "LitterT4R2 (Rep 2)", "LitterT1R1 (Rep 2)", 
                  "LitterT1R1 (Rep 1)", "ROT5R1 (Rep1)", "ROT5R1 (Rep 2)", "Negative1","PCR Neg1")
table <- table[,!(names(table) %in% mock_and_neg)]

# Make the OTU ID the row name
row.names(table) <- table$`#OTU ID`
table <- table[ , !(names(table) %in% c("#OTU ID"))]

# and transform the dataframes because vegan expects rows = samples and columns = species
table <- as.data.frame(t(table))
table <- table[ , !(names(table) %in% c("SampleID"))]

###______ Making NMDS:
# First, calculate the BC matrix so you can do your rarefactions and sqrt transformations for EACH BC matrix
bray.dist <- avgdist(table, sample = 1000, meanfun = median, transf = sqrt, iterations = 999) #default is Bray-Curtis calculated by vegdist()

#Next, we will run an NMDS, which is a form of ordination and can be used to visualize beta diversity.
NMDS1 <- metaMDS(bray.dist, autotransform = FALSE, k = 2)

#This will make the first two columns as the x and y coordinates, so that they may be plotted.
coordinates <- data.frame(NMDS1$points[,1:2])

# Update metadata
nmds_plus_metadata <- merge(coordinates, metadata, by.x = "row.names", by.y = "SampleID")
nmds_plus_metadata$Dispersal_Route <- ifelse(nmds_plus_metadata$Dispersal_Route == "Air", "Environmental Air", 
                                             ifelse(nmds_plus_metadata$Dispersal_Route == "Soil", "Environmental Soil", 
                                                    ifelse(nmds_plus_metadata$Dispersal_Route == "Elevated", "Elevated (Air)", 
                                                           ifelse(nmds_plus_metadata$Dispersal_Route == "Up", "Overhead (Air + Vegetation)", 
                                                                  ifelse(nmds_plus_metadata$Dispersal_Route == "Open", "Open (Air + Vegetation + Soil)", nmds_plus_metadata$Dispersal_Route)))))
nmds_plus_metadata$Environmental <- ifelse(grepl("Environmental", nmds_plus_metadata$Dispersal_Route), "Yes", "No")


#Time to make a new plot colored with metadata.
cols <- c("Closed" = "#663A87",
          "Elevated (Air)" = "#dbbd29", 
          "Overhead (Air + Vegetation)" = "#d94188", 
          "Open (Air + Vegetation + Soil)" = "#a41d2c",
          "Environmental Litter" = "#d94188")

nmds_plus_metadata$Dispersal_Route <- factor(nmds_plus_metadata$Dispersal_Route, levels = c(
  "Closed", "Elevated (Air)", "Overhead (Air + Vegetation)", "Open (Air + Vegetation + Soil)", "Environmental Litter"))

nmds_plus_metadata$Environmental <- factor(nmds_plus_metadata$Environmental, levels = c("Yes", "No"))

pdf("/Users/walters_kendra/Google Drive/Dispersal/Figures/For_Manuscript/Figure_4A.pdf", width = 9, height = 4)
ggplot(data = nmds_plus_metadata) +
  geom_point(aes(x = MDS1, y = MDS2, color = as.factor(nmds_plus_metadata$Dispersal_Route), 
                 shape=as.factor(Environmental)), size = 4) + #Creates and colors legend to match, modify after $ here.
  labs(col =  "Dispersal Route") +
  theme_classic() + 
  scale_color_manual(values = cols) +
  scale_shape_manual(values = c(21,19))
dev.off()









#####__________ For Panel B = taxa barplot over time underneath for [Open?] treatment__________######
# August 26th 2020
# Kendra Walters

# A = Closed
# B (across from A) = Elevated
# C = Overhead
# D = Open
# E = Environmental Litter

# room of ()ments
require(data.table)
require(ggplot2)
require(vegan)
require(anchors)
library(stringr)
library(tidyr)
library(tidyverse)

# Setting up and reading in data
q2.table.rarefied <- as.data.frame(fread("OTU_filtered_plus_rarefied__ROUNDED.tsv"))

# Doing the part of the taxa bar plot script that doesn't need to be run over and over and over
row.names(q2.table.rarefied) <- q2.table.rarefied$V1 # set sample names as row names
q2.table.rarefied$V1 <- NULL # delete column of sample names
rowSums(q2.table.rarefied) # check to see if it's the right dataset


q2.table.rarefied$Dispersal_Route <- ifelse(grepl("^RO", row.names(q2.table.rarefied)), "Open",
                                            ifelse(grepl("^RC", row.names(q2.table.rarefied)), "Closed",
                                                   ifelse(grepl("^RU", row.names(q2.table.rarefied)), "Overhead",
                                                          ifelse(grepl("^RE", row.names(q2.table.rarefied)), "Elevated",
                                                                 ifelse(grepl("Litter", row.names(q2.table.rarefied)), "Environmental Litter", "ooops.....")))))

q2.table.rarefied$TimePoint <- str_extract(string = row.names(q2.table.rarefied), pattern = "T[1-5]")
q2.table.rarefied$Route_AND_Time <- paste0(q2.table.rarefied$Dispersal_Route, "_" ,q2.table.rarefied$TimePoint)


q2.family <- q2.table.rarefied[, -which(names(q2.table.rarefied) %in% c("Dispersal_Route", "TimePoint"))] # Remove the Dispersal Route column, since you've already chosen your variable

q2.family.agg <- aggregate(. ~ Route_AND_Time, data = q2.family, FUN = mean) # Take the average of the communities for each timepoint
rowSums(q2.family.agg[2:length(names(q2.family.agg))]) # check to make sure each agg sample still adds up to rarefaction depth
q2.family.agg <- as.data.frame(t(q2.family.agg)) # transpose table so factor names are on tp
names(q2.family.agg) <- as.character(unlist(q2.family.agg["Route_AND_Time", ])) # Set factor names as column names
q2.family.agg <- q2.family.agg[-1,] # get rid of row that gives factor names


# set the data type to be numberic 
for(i in c(1:ncol(q2.family.agg))) {
  q2.family.agg[,i] <- as.numeric(as.character(q2.family.agg[,i]))
}


# merge taxonomy with this new OTU table
q2.taxonomy <- as.data.frame(fread("15_taxonomy_filtered.tsv")) # import taxonomy file (3 columns --> OTU ID, taxonomy, and confidence)

colSums(q2.family.agg) # check to make sure it's still rarefied and nothing weird happened
q2.fam.agg.tax <- merge(q2.family.agg, q2.taxonomy[1:2], by.x = "row.names", by.y = "#OTUID") # merge taxonomy and otu table

q2.fam.agg.tax$family <- str_extract(q2.fam.agg.tax$taxonomy, "g__[:alpha:]*") # Extract genus names
q2.fam.agg.tax$family <- gsub("g__", "", q2.fam.agg.tax$family) # Make them pretty

q2.fam.agg.tax$family <- factor(q2.fam.agg.tax$family)

q2.fam.agg.tax$OTUids <- row.names(q2.fam.agg.tax)
q2.fam.agg.tax$Row.names <- NULL
q2.fam.agg.tax$taxonomy <- NULL


q2.fam.agg.tax.agg <- aggregate(. ~ family, data = q2.fam.agg.tax[c(1:(length(names(q2.fam.agg.tax))-2), match("family", names(q2.fam.agg.tax)))], FUN = sum) # aggregate by fammily
q2.fam.agg.tax.agg <- q2.fam.agg.tax.agg[!(q2.fam.agg.tax.agg$family == ""), ] # get rid of all sequences that weren't identified to family


library(vegan)
rownames(q2.fam.agg.tax.agg) <- q2.fam.agg.tax.agg$family
q2.fam.agg.tax.agg$family <- NULL
q2.fam.agg.tax.agg.tot <- decostand(q2.fam.agg.tax.agg, method="total", MARGIN=2) 
q2.fam.agg.tax.agg.tot$family <- rownames(q2.fam.agg.tax.agg.tot) 

q2.fam.agg.tax.agg.tot.long <- gather(q2.fam.agg.tax.agg.tot, Route_AND_Time, Abundance, Closed_T1:Overhead_T5) # you have to change the column range for a different dataset


#call everything below 0.09 (9%) "Other 'classification level'"
keep <- list()
PERCENTAGE <- 0.09 # fill this out to to the cutoff you want

samp <- unique(q2.fam.agg.tax.agg.tot.long$Route_AND_Time)

for (i in 1:length(samp)) {
  sample_type <- q2.fam.agg.tax.agg.tot.long[q2.fam.agg.tax.agg.tot.long$Route_AND_Time == samp[i], ]
  high_abund <- sample_type[sample_type$Abundance > PERCENTAGE, ]
  vec <- c(as.character(high_abund$family))
  keep[[i]] <- vec
}
keep.unique <- unique(unlist(keep))
# keep.unique

# Rename all families NOT in "keep.unique" to "Other Family"
q2.fam.agg.tax.agg.tot.long$family <- ifelse(q2.fam.agg.tax.agg.tot.long$family %in% keep.unique, as.character(q2.fam.agg.tax.agg.tot.long$family), "Other Genera")

q2.FINAL <- aggregate(x = q2.fam.agg.tax.agg.tot.long$Abundance, 
                      by = list(q2.fam.agg.tax.agg.tot.long$family, q2.fam.agg.tax.agg.tot.long$Route_AND_Time), 
                      FUN = sum)
names(q2.FINAL) <- c("Genus", "Route_AND_Time", "Relative Abundance")

# Checking to make sure it adds to 1 precisely
for (i in 1:length(samp)) {
  print(paste0("For Route_AND_Time ", samp[i], " the relative abundance adds to: ", sum(q2.fam.agg.tax.agg.tot.long[q2.fam.agg.tax.agg.tot.long$Route_AND_Time == samp[i], ]$Abundance), "."))
}


q2.FINAL$Genus <- factor(q2.FINAL$Genus, levels = c("Other Genera", sort(keep.unique)))
q2.FINAL$Dispersal_Route <- str_extract(string = q2.FINAL$Route_AND_Time, pattern = "^[a-zA-Z\ ]*")
q2.FINAL$TimePoint <- str_extract(string = q2.FINAL$Route_AND_Time, pattern = "T[1-5]")
q2.FINAL.double <- q2.FINAL
q2.FINAL.double$TimePoint <- c(rep("Avg", nrow(q2.FINAL.double)))
q2.FINAL.double$`Relative Abundance` <- q2.FINAL.double$`Relative Abundance` / length(unique(q2.FINAL$TimePoint))
q2.FINAL.together <- rbind(q2.FINAL, q2.FINAL.double)
q2.FINAL.together$width <- ifelse(q2.FINAL.together$TimePoint == "Avg", 1, 0.5)

#q2.FINAL.together$Genus <- factor(q2.FINAL.together$Genus, levels = c("Other Genera", sort(keep.unique))) # First, "Other Genera", then in order of abundance
q2.FINAL.together$Genus <- factor(q2.FINAL.together$Genus, levels = c("Other Genera", "Erwinia", "Janthinobacterium", "Hymenobacter",
                                                                      "Pseudomonas", "Paenibacillus", "Curtobacterium", "Frigoribacterium", 
                                                                      "Sphingomonas", "Methylobacterium", "Sanguibacter"))

q2.FINAL.together$Dispersal_Route <- factor(q2.FINAL.together$Dispersal_Route, levels = 
                                              c("Closed", "Elevated", "Overhead", "Open", "Environmental Litter"))

color = c("#a9a9a9", '#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', 
          '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')

ggplot(q2.FINAL.together) +
  geom_bar(aes(x = TimePoint, y = `Relative Abundance`, fill = Genus, width = width), stat = "identity") +
  labs(x = "Time Point", y = "Proportion of Total Community (identified to Genus)") +
  facet_wrap( ~ Dispersal_Route, ncol = 2) +
  theme_test() +
  #scale_color_viridis(discrete=TRUE) 
  scale_fill_manual(values = color, name = "Genera") 