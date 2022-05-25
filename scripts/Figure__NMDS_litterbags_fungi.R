# Fungi NMDS + taxa barplot
# Kendra E Walters
# December 3rd 2020

# Room of ()ments
require(data.table)
require(ggplot2)
require(vegan)
require(anchors)
library(stringr)
library(tidyr)
library(tidyverse)
library(gridExtra)
require(ggvenn)

#####__________ Panel A = NMDS _________________________ 
# Read in data
fungi <- as.data.frame(fread("metadata_for_PRIMER_because_stupid_sample_order_q2_single.txt"))

# Make better names
fungi$Environmental <- ifelse(fungi$Dispersal_Route == "Environmental Litter", "Yes", "No")

treatment.key <- c("Environmental Litter" = "Environmental Litter", 
                   "Phyllosphere" = "Phyllosphere", 
                   "Closed" = "Closed", 
                   "Elevated" = "Elevated (Air)", 
                   "Up" = "Overhead (Air + Vegetation)", 
                   "Open" = "Open (Air + Vegetation + Soil)")
fungi$Treatments <- treatment.key[fungi$Dispersal_Route]

# Make the graph
cols <- c("Closed" = "#474748",
          "Elevated (Air)" = "#dbbd28", 
          "Overhead (Air + Vegetation)" = "#662f90", 
          "Open (Air + Vegetation + Soil)" = "#9b2a31",
          "Environmental Litter" = "#662f90", 
          "Phyllosphere" = "black")

fungi$Treatments <- factor(fungi$Treatments, levels = c(
  "Closed", "Elevated (Air)", "Overhead (Air + Vegetation)", "Open (Air + Vegetation + Soil)", "Environmental Litter", "Phyllosphere"))

fungi$Environmental <- factor(fungi$Environmental, levels = c("Yes", "No"))


ggplot(data = fungi) +
  geom_point(aes(x = MDS1, y = MDS2, color = Treatments, 
                 shape=as.factor(Environmental)), size = 4) + #Creates and colors legend to match, modify after $ here.
  labs(col =  "Dispersal Route") +
  theme_classic() + 
  scale_color_manual(values = cols) +
  scale_shape_manual(values = c(21,19))


#####__________ Panel B = taxa barplot _________________ #####

### ______________________________________ PHYLLOSPHERE
# Setting up and reading in data
q2.table.rarefied <- as.data.frame(fread("OTU_filtered_plus_rarefied_3500__ROUND_IS_TRUE.tsv"))

# Doing the part of the taxa bar plot script that doesn't need to be run over and over and over
row.names(q2.table.rarefied) <- q2.table.rarefied$V1 # set sample names as row names
q2.table.rarefied$V1 <- NULL # delete column of sample names
rowSums(q2.table.rarefied) # check to see if it's the right dataset


q2.table.rarefied$Dispersal_Route <- ifelse(grepl("^RO", row.names(q2.table.rarefied)), "Open",
                                            ifelse(grepl("^RC", row.names(q2.table.rarefied)), "Closed",
                                                   ifelse(grepl("^RU", row.names(q2.table.rarefied)), "Overhead",
                                                          ifelse(grepl("^RE", row.names(q2.table.rarefied)), "Elevated",
                                                                 ifelse(grepl("Litter", row.names(q2.table.rarefied)), "Environmental Litter", 
                                                                        ifelse(grepl("Phyllo", row.names(q2.table.rarefied)), "Phyllosphere","ooops....."))))))
phyllo <- q2.table.rarefied[q2.table.rarefied$Dispersal_Route == "Phyllosphere", ]
phyllo$Dispersal_Route <- NULL

q2.family.agg <- as.data.frame(t(phyllo)) # transpose table so factor names are on top

# set the data type to be numeric 
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
q2.fam.agg.tax.agg <- q2.fam.agg.tax.agg[!(q2.fam.agg.tax.agg$family == ""), ] # get rid of all sequences that weren't identified to genus

rownames(q2.fam.agg.tax.agg) <- q2.fam.agg.tax.agg$family
q2.fam.agg.tax.agg$family <- NULL
q2.fam.agg.tax.agg.tot <- decostand(q2.fam.agg.tax.agg, method="total", MARGIN=2)
q2.fam.agg.tax.agg.tot$family <- rownames(q2.fam.agg.tax.agg.tot) 

q2.fam.agg.tax.agg.tot.long <- gather(q2.fam.agg.tax.agg.tot, Sample, Abundance, `PhylloT0R2 (Rep 1)`:PhylloT0R7) # you have to change the column range for a different dataset


#call everything below 0.04 (4%) "Other 'classification level'"
keep <- list()
PERCENTAGE <- 0.04 # fill this out to to the cutoff you want

samp <- unique(q2.fam.agg.tax.agg.tot.long$Sample)

for (i in 1:length(samp)) {
  sample_type <- q2.fam.agg.tax.agg.tot.long[q2.fam.agg.tax.agg.tot.long$Sample == samp[i], ]
  high_abund <- sample_type[sample_type$Abundance > PERCENTAGE, ]
  vec <- c(as.character(high_abund$family))
  keep[[i]] <- vec
}
keep.unique <- unique(unlist(keep))
#keep.unique <- keep.unique[keep.unique != "unidentified"]
keep.unique

# Rename all families NOT in "keep.unique" to "Other Family"
q2.fam.agg.tax.agg.tot.long$family <- ifelse(q2.fam.agg.tax.agg.tot.long$family %in% keep.unique, as.character(q2.fam.agg.tax.agg.tot.long$family), "Other Genera")

q2.FINAL <- aggregate(x = q2.fam.agg.tax.agg.tot.long$Abundance, 
                      by = list(q2.fam.agg.tax.agg.tot.long$family, q2.fam.agg.tax.agg.tot.long$Sample), 
                      FUN = sum)
names(q2.FINAL) <- c("Genus", "Sample", "Relative Abundance")

# Checking to make sure it adds to 1 precisely
for (i in 1:length(samp)) {
  print(paste0("For Sample ", samp[i], " the relative abundance adds to: ", sum(q2.fam.agg.tax.agg.tot.long[q2.fam.agg.tax.agg.tot.long$Sample == samp[i], ]$Abundance), "."))
}


q2.FINAL$Genus <- factor(q2.FINAL$Genus, levels = c("Other Genera", sort(keep.unique)))

#color <- c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628','#f781bf','#999999')

color = c("#a9a9a9", '#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', 
          '#fabebe', '#008080',  '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', 
          '#808080', '#ffffff', '#000000')
A <- ggplot(q2.FINAL) +
  geom_bar(aes(x = Sample, y = `Relative Abundance`, fill = Genus), stat = "identity") +
  labs(x = "Sample", y = "Proportion of Total Community (identified to Genus)") +
  theme_test() +
  #scale_color_viridis(discrete=TRUE) 
  scale_fill_manual(values = color, name = "Genera") 


########_________________________ EVERYTHING ELSE
# Setting up and reading in data
q2.table.rarefied <- as.data.frame(fread("OTU_filtered_plus_rarefied_3500__ROUND_IS_TRUE.tsv"))

# Doing the part of the taxa bar plot script that doesn't need to be run over and over and over
row.names(q2.table.rarefied) <- q2.table.rarefied$V1 # set sample names as row names
q2.table.rarefied$V1 <- NULL # delete column of sample names
rowSums(q2.table.rarefied) # check to see if it's the right dataset


q2.table.rarefied$Dispersal_Route <- ifelse(grepl("^RO", row.names(q2.table.rarefied)), "Open",
                                            ifelse(grepl("^RC", row.names(q2.table.rarefied)), "Closed",
                                                   ifelse(grepl("^RU", row.names(q2.table.rarefied)), "Overhead",
                                                          ifelse(grepl("^RE", row.names(q2.table.rarefied)), "Elevated",
                                                                 ifelse(grepl("Litter", row.names(q2.table.rarefied)), "Environmental Litter", 
                                                                        ifelse(grepl("Phyllo", row.names(q2.table.rarefied)), "Phyllosphere","ooops....."))))))
q2.table.rarefied <- q2.table.rarefied[q2.table.rarefied$Dispersal_Route != "Phyllosphere", ]

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

q2.fam.agg.tax.agg <- aggregate(. ~ family, data = q2.fam.agg.tax[c(1:(length(names(q2.fam.agg.tax))-2), match("family", names(q2.fam.agg.tax)))], FUN = sum) # aggregate by genus
q2.fam.agg.tax.agg <- q2.fam.agg.tax.agg[!(q2.fam.agg.tax.agg$family == ""), ] # get rid of all sequences that weren't identified to genus

rownames(q2.fam.agg.tax.agg) <- q2.fam.agg.tax.agg$family
q2.fam.agg.tax.agg$family <- NULL
q2.fam.agg.tax.agg.tot <- decostand(q2.fam.agg.tax.agg, method="total", MARGIN=2)
q2.fam.agg.tax.agg.tot$family <- rownames(q2.fam.agg.tax.agg.tot) 

q2.fam.agg.tax.agg.tot.long <- gather(q2.fam.agg.tax.agg.tot, Route_AND_Time, Abundance, Closed_T1:Overhead_T5) # you have to change the column range for a different dataset


#call everything below 0.07 (7%) "Other 'classification level'"
keep <- list()
PERCENTAGE <- 0.07 # fill this out to to the cutoff you want

samp <- unique(q2.fam.agg.tax.agg.tot.long$Route_AND_Time)

for (i in 1:length(samp)) {
  sample_type <- q2.fam.agg.tax.agg.tot.long[q2.fam.agg.tax.agg.tot.long$Route_AND_Time == samp[i], ]
  high_abund <- sample_type[sample_type$Abundance > PERCENTAGE, ]
  vec <- c(as.character(high_abund$family))
  keep[[i]] <- vec
}
keep.unique <- unique(unlist(keep))
keep.unique <- keep.unique[keep.unique != "unidentified"]
keep.unique

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
q2.FINAL.together$Genus <- factor(q2.FINAL.together$Genus, levels = c("Other Genera", "Alternaria", "Cladosporium", "Paraconiothyrium", "Pyrenophora", 
                                                                      "Leptosphaeria", "Crocicreas", "Aureobasidium", "Filobasidium"))

q2.FINAL.together$Dispersal_Route <- factor(q2.FINAL.together$Dispersal_Route, levels = 
                                              c("Closed", "Elevated", "Overhead", "Open", "Environmental Litter"))

# color = c("#a9a9a9", '#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', 
#           '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')

#color <- c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628','#f781bf','#999999')
color = c("#a9a9a9", '#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', 
          '#fabebe', '#008080',  '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', 
          '#808080', '#ffffff', '#000000')

B <- ggplot(q2.FINAL.together) +
  geom_bar(aes(x = TimePoint, y = `Relative Abundance`, fill = Genus, width = width), stat = "identity") +
  labs(x = "Time Point", y = "Proportion of Total Community (identified to Genus)") +
  facet_wrap( ~ Dispersal_Route, ncol = 2) +
  theme_test() +
  #scale_color_viridis(discrete=TRUE) 
  scale_fill_manual(values = color, name = "Genera") 


####### ________ put them together
grid.arrange(A, B)