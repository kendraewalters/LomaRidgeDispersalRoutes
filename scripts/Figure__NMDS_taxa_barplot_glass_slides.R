# Panel A = NMDS, Glass slides, include: elevated, open, up, soil samples, leaf litter, and air
# Panel B = taxa barplot for average, and each time point for each treatment + environmental sample

# Kendra E. Walters
# April 6th 2020

# Room of ()ments
require(data.table)
require(ggplot2)
require(vegan)
require(ggvenn)
require(anchors)
library(stringr)
library(tidyr)
library(tidyverse)


######## _____________ PANEL A = NMDS  _____________________ #################
# Loading data
metadata <- as.data.frame(fread("metadata_q1.txt"))

q1.table <- as.data.frame(fread("13_table_filtered.tsv")) # now you DON"T want to use rarefied data here because you can rarefy within your BC calculation
#   --p-trim-left 5 \
#   --p-trunc-len 279 \

#First, filter out mock, negative community standards, and replicates by name.
mock_and_neg <- c("Negative2", "PCR Mock1", "PCR Mock2", "PCR Neg1", "PCR Neg2", "Negative1", 
                  "LitterT1R1 (Rep 1)", "LitterT1R1 (Rep 2)", "LitterT4R2 (Rep 1)", "LitterT4R2 (Rep 2)")
q1.table <- q1.table[,!(names(q1.table) %in% mock_and_neg)]

# Then, filter out any treatment you DON'T want in the graph -- (i.e. death)
q1.table <- q1.table[ , !(grepl(pattern = "LD", names(q1.table)))]

# Make the OTU ID the row name
row.names(q1.table) <- q1.table$`#OTU ID`
q1.table <- q1.table[ , !(names(q1.table) %in% c("#OTU ID"))]

# and transform the dataframes because vegan expects rows = samples and columns = species
q1.table <- as.data.frame(t(q1.table))


###______ Making NMDS:
# First, calculate the BC matrix so you can do your rarefactions and sqrt transformations for EACH BC matrix
bray.dist <- avgdist(q1.table, sample = 1000, meanfun = median, transf = sqrt, iterations = 999) #default is Bray-Curtis calculated by vegdist()

#Next, we will run an NMDS, which is a form of ordination and can be used to visualize beta diversity.
NMDS1 <- metaMDS(bray.dist, autotransform = FALSE, k = 2)

#This will make the first two columns as the x and y coordinates, so that they may be plotted.
coordinates <- data.frame(NMDS1$points[,1:2])

#To make a more sophisticated plot, we will merge the stress scores with the metadata.
nmds_plus_metadata <- merge(coordinates, metadata, by.x = "row.names", by.y = "SampleID")
nmds_plus_metadata$Dispersal_Route <- ifelse(nmds_plus_metadata$Dispersal_Route == "Air", "Environmental Air", 
                                             ifelse(nmds_plus_metadata$Dispersal_Route == "Soil", "Environmental Soil", 
                                                    ifelse(nmds_plus_metadata$Dispersal_Route == "Elevated", "Elevated (Air)", 
                                                           ifelse(nmds_plus_metadata$Dispersal_Route == "Up", "Overhead (Air + Vegetation)", 
                                                                  ifelse(nmds_plus_metadata$Dispersal_Route == "Open", "Open (Air + Vegetation + Soil)", nmds_plus_metadata$Dispersal_Route)))))
nmds_plus_metadata$Slides <- ifelse(grepl("Environmental", nmds_plus_metadata$Dispersal_Route), "No", "Yes")


#Time to make a new plot colored with metadata.
cols <- c("Elevated (Air)" = "#dbbd29", 
          "Overhead (Air + Vegetation)" = "#d94188", 
          "Open (Air + Vegetation + Soil)" = "#a41d2c", 
          "Environmental Air" = "#dbbd29", 
          "Environmental Litter" = "#d94188", 
          "Environmental Soil" = "#a41d2c")

nmds_plus_metadata$Dispersal_Route <- factor(nmds_plus_metadata$Dispersal_Route, levels = c(
  "Elevated (Air)", "Overhead (Air + Vegetation)", "Open (Air + Vegetation + Soil)", 
  "Environmental Air", "Environmental Litter", "Environmental Soil"))

ggplot(data = nmds_plus_metadata) +
  geom_point(aes(x = MDS1, y = MDS2, color = as.factor(Dispersal_Route), 
                                                       shape=as.factor(Slides)), 
             size = 4, stroke = 1) + #Creates and colors legend to match, modify after $ here.
  labs(col =  "Dispersal Route") +
  theme_classic() + 
  scale_shape_manual(values = c(21,19)) + 
  scale_color_manual(values = cols)





######## _____________ PANEL B = taxa bar plot _____________ ############
# A = elevated
# B = up
# C = open
# D = air
# E = environmental litter
# F = soil 
# Kendra E. Walters
# May 21st 2020


# Setting up and reading in data
q1.table.rarefied <- as.data.frame(fread("OTU_filtered_plus_rarefied__ROUNDED.tsv"))

# Doing the part of the taxa bar plot script that doesn't need to be run over and over and over
row.names(q1.table.rarefied) <- q1.table.rarefied$V1 # set sample names as row names
q1.table.rarefied$V1 <- NULL # delete column of sample names
rowSums(q1.table.rarefied) # check to see if it's the right dataset


q1.table.rarefied$Dispersal_Route <- ifelse(grepl("^LO", row.names(q1.table.rarefied)), "Open",
                                            ifelse(grepl("^LC", row.names(q1.table.rarefied)), "Closed",
                                                   ifelse(grepl("^LU", row.names(q1.table.rarefied)), "Overhead",
                                                          ifelse(grepl("^LE", row.names(q1.table.rarefied)), "Elevated",
                                                                 ifelse(grepl("Litter", row.names(q1.table.rarefied)), "Environmental Litter",
                                                                        ifelse(grepl("LD", row.names(q1.table.rarefied)), "Death",
                                                                               ifelse(grepl("Soil", row.names(q1.table.rarefied)), "Soil",
                                                                                      ifelse(grepl("Air", row.names(q1.table.rarefied)), "Air", "ooops....."))))))))

q1.table.rarefied$TimePoint <- str_extract(string = row.names(q1.table.rarefied), pattern = "T[1-5]")
q1.table.rarefied$Route_AND_Time <- paste0(q1.table.rarefied$Dispersal_Route, "_" ,q1.table.rarefied$TimePoint)

q1.family <- q1.table.rarefied[q1.table.rarefied$Dispersal_Route != "Death", ]
q1.family <- q1.family[, -which(names(q1.family) %in% c("Dispersal_Route", "TimePoint"))] # Remove the Dispersal Route column, since you've already chosen your variable

q1.family.agg <- aggregate(. ~ Route_AND_Time, data = q1.family, FUN = mean) # Take the average of the communities for each timepoint
rowSums(q1.family.agg[2:length(names(q1.family.agg))]) # check to make sure each agg sample still adds up to rarefaction depth
q1.family.agg <- as.data.frame(t(q1.family.agg)) # transpose table so factor names are on tp
# good here
names(q1.family.agg) <- as.character(unlist(q1.family.agg["Route_AND_Time", ])) # Set factor names as column names
q1.family.agg <- q1.family.agg[-1,] # get rid of row that gives factor names


# set the data type to be numberic 
for(i in c(1:ncol(q1.family.agg))) {
  q1.family.agg[,i] <- as.numeric(as.character(q1.family.agg[,i]))
}


# merge taxonomy with this new OTU table
q1.taxonomy <- as.data.frame(fread("14_taxonomy_filtered.tsv")) # import taxonomy file (3 columns --> OTU ID, taxonomy, and confidence)

colSums(q1.family.agg) # check to make sure it's still rarefied and nothing weird happened
q1.fam.agg.tax <- merge(q1.family.agg, q1.taxonomy[1:2], by.x = "row.names", by.y = "#OTUID") # merge taxonomy and otu table

q1.fam.agg.tax$family <- str_extract(q1.fam.agg.tax$taxonomy, "g__[:alpha:]*") # Extract family names
q1.fam.agg.tax$family <- gsub("g__", "", q1.fam.agg.tax$family) # Make them pretty

q1.fam.agg.tax$family <- factor(q1.fam.agg.tax$family)

q1.fam.agg.tax$OTUids <- row.names(q1.fam.agg.tax)
q1.fam.agg.tax$Row.names <- NULL
q1.fam.agg.tax$taxonomy <- NULL


q1.fam.agg.tax.agg <- aggregate(. ~ family, data = q1.fam.agg.tax[c(1:(length(names(q1.fam.agg.tax))-2), match("family", names(q1.fam.agg.tax)))], FUN = sum) # aggregate by fammily
# good here
q1.fam.agg.tax.agg <- q1.fam.agg.tax.agg[!(q1.fam.agg.tax.agg$family == ""), ] # get rid of all sequences that weren't identified to family


library(vegan)
rownames(q1.fam.agg.tax.agg) <- q1.fam.agg.tax.agg$family
q1.fam.agg.tax.agg$family <- NULL
q1.fam.agg.tax.agg.tot <- decostand(q1.fam.agg.tax.agg, method="total", MARGIN=2) 
q1.fam.agg.tax.agg.tot$family <- rownames(q1.fam.agg.tax.agg.tot) 

q1.fam.agg.tax.agg.tot.long <- gather(q1.fam.agg.tax.agg.tot, Route_AND_Time, Abundance, Air_T1:Soil_T5)


#call everything below 0.09 (9%) "Other 'classification level'"
keep <- list()
PERCENTAGE <- 0.09 # fill this out to to the cutoff you want

samp <- unique(q1.fam.agg.tax.agg.tot.long$Route_AND_Time)

for (i in 1:length(samp)) {
  sample_type <- q1.fam.agg.tax.agg.tot.long[q1.fam.agg.tax.agg.tot.long$Route_AND_Time == samp[i], ]
  high_abund <- sample_type[sample_type$Abundance > PERCENTAGE, ]
  vec <- c(as.character(high_abund$family))
  keep[[i]] <- vec
}
keep.unique <- unique(unlist(keep))
# keep.unique

# Rename all families NOT in "keep.unique" to "Other Genera"
q1.fam.agg.tax.agg.tot.long$family <- ifelse(q1.fam.agg.tax.agg.tot.long$family %in% keep.unique, as.character(q1.fam.agg.tax.agg.tot.long$family), "Other Genera")

q1.FINAL <- aggregate(x = q1.fam.agg.tax.agg.tot.long$Abundance, 
                      by = list(q1.fam.agg.tax.agg.tot.long$family, q1.fam.agg.tax.agg.tot.long$Route_AND_Time), 
                      FUN = sum)
names(q1.FINAL) <- c("Genus", "Route_AND_Time", "Relative Abundance")

# Checking to make sure it adds to 1 precisely
for (i in 1:length(samp)) {
  print(paste0("For Route_AND_Time ", samp[i], " the relative abundance adds to: ", sum(q1.fam.agg.tax.agg.tot.long[q1.fam.agg.tax.agg.tot.long$Route_AND_Time == samp[i], ]$Abundance), "."))
}


q1.FINAL$Genus <- factor(q1.FINAL$Genus, levels = c("Other Genera", sort(keep.unique)))
q1.FINAL$Dispersal_Route <- str_extract(string = q1.FINAL$Route_AND_Time, pattern = "^[a-zA-Z\ ]*")
q1.FINAL$TimePoint <- str_extract(string = q1.FINAL$Route_AND_Time, pattern = "T[1-5]")
q1.FINAL.double <- q1.FINAL
q1.FINAL.double$TimePoint <- c(rep("Avg", nrow(q1.FINAL.double)))
q1.FINAL.double$`Relative Abundance` <- q1.FINAL.double$`Relative Abundance` / length(unique(q1.FINAL$TimePoint))
q1.FINAL.together <- rbind(q1.FINAL, q1.FINAL.double)
q1.FINAL.together$width <- ifelse(q1.FINAL.together$TimePoint == "Avg", 1, 0.5)

#q1.FINAL.together$Genus <- factor(q1.FINAL.together$Genus, levels = c("Other Genera", sort(keep.unique))) # First, "Other Genera", then in order of abundance
q1.FINAL.together$Genus <- factor(q1.FINAL.together$Genus, levels = c("Other Genera", "Hymenobacter", "Sphingomonas",
                                                                      "Janthinobacterium", "Curtobacterium", "Methylobacterium", 
                                                                      "Bacillus", "Balneimonas", "Skermanella", 
                                                                      "Bradyrhizobium", "Fervidobacterium", "Paracoccus"
)) # First, "Other Genera", then in order of abundance

q1.FINAL.together$Dispersal_Route <- factor(q1.FINAL.together$Dispersal_Route, levels = 
                                              c("Air", "Elevated", "Environmental Litter", "Overhead", "Soil", "Open"))

color = c("#a9a9a9", '#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', 
          '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
color.99.accessible <- c("#a9a9a9", "#f58231", "#4363d8", "#3cb44b", "#ffe119", "#000075", 
                         "#e6194B", "#f032e6", "#e6beff", "#42d4f4", "#aaffc3", "#800000")

ggplot(q1.FINAL.together) +
  geom_bar(aes(x = TimePoint, y = `Relative Abundance`, fill = Genus, width = width), stat = "identity") +
  labs(x = "Time Point", y = "Proportion of Total Community (identified to Genus)") +
  facet_wrap( ~ Dispersal_Route, ncol = 2) +
  theme_test() +
  #scale_color_viridis(discrete=TRUE) 
  scale_fill_manual(values = color, name = "Genera")