# A = glass slide Shannon diversity by route
# B = litterbag bacterial Shannon diversity by route
# C = litterbag fungal Shannon diversity by route
# D = glass slide beta-diversity by route
# E = litterbag bacterial beta-diversity by route
# F = litterbag fungal beta-diversity by route

# Kendra E. Walters
# September 23rd, 2020


# Room of ()ments
require(ggplot2)
require(data.table)
require(vegan)
require(gridExtra)
require(stringr)
require(ggpattern)


## Set up treatment labels for x-axis (want environmental litter to be on two lines)
treatments <- c("Air", "Environmental Litter", "Soil", "Elevated", "Overhead", "Open")
treatment_labels <- setNames(str_replace(treatments, " ", "\n "), treatments) # put in line break at the space


##### __________________ Glass Shannon = A  _________________ #########
q1.table.rarefied <- as.data.frame(fread("OTU_filtered_plus_rarefied__ROUNDED.tsv"))
row.names(q1.table.rarefied) <- q1.table.rarefied$V1
q1.table.rarefied$V1 <- NULL

# Read in metadata!
metadata <- as.data.frame(fread("metadata_q1.txt"))

# Calculate alpha diversity!
shannon <- as.data.frame(diversity(q1.table.rarefied, index = "shannon"))
names(shannon) <- c("Shannon")
richness <- as.data.frame(apply(q1.table.rarefied[,]>0, 1, sum))
names(richness) <- c("Richness")
simpson <- as.data.frame(diversity(q1.table.rarefied, index = "simpson"))
names(simpson) <- c("Simpson")

# And merge them all together!
tmp1 <- merge(richness, shannon, by.x = "row.names", by.y = "row.names")
all.alpha <- merge(tmp1, simpson, by.x = "Row.names", by.y = "row.names")
names(all.alpha) <- c("SampleID", names(all.alpha)[2:length(names(all.alpha))])

#Merge with metadata to create a plot.
merged_alpha <- merge(all.alpha, metadata, by.x = "SampleID", by.y = "SampleID")

# Make-the-plot
merged_alpha <- merged_alpha[merged_alpha$Dispersal_Route != "Death", ]
merged_alpha$Dispersal_Route_new_names <- ifelse(merged_alpha$Dispersal_Route == "Up", "Overhead", merged_alpha$Dispersal_Route)
merged_alpha$Dispersal_Route_new_names <- factor(merged_alpha$Dispersal_Route_new_names, 
                                       levels = c("Air", "Environmental Litter", "Soil", "Elevated", "Overhead", "Open"))

cols <- c("Air" = "#dbbd28", 
          "Environmental Litter" = "#662f90",
          "Soil" = "#9b2a31", 
          "Elevated" = "#dbbd28", 
          "Overhead" = "#662f90", 
          "Open" = "#9b2a31")

A <- ggplot(data = merged_alpha) + 
  geom_boxplot_pattern(aes(x = Dispersal_Route_new_names,
                           y = Shannon, 
                           fill = Dispersal_Route_new_names,
                           pattern_density = Dispersal_Route_new_names), 
                       color = "black", 
                       pattern_fill = "white", 
                       pattern_spacing = 0.05) +
  theme_classic() + 
  labs(x = "Dispersal Route", y = "Shannon") +
  scale_fill_manual(values = cols) + 
  scale_pattern_density_manual(values = c(0.5, 0.5, 0.5, 0, 0, 0)) + 
  theme(legend.position = "none") +
  guides(fill = FALSE) + 
  xlab("Treatment") + ylab("Shannon Diversity Index") +
  scale_x_discrete(labels = treatment_labels) + 
  scale_y_continuous(limits = c(0, 6))
  #coord_fixed(ratio = 2/3)

##### __________________ Glass BC = B _______________________ #####
## Load data
bray.dist.df <- as.data.frame(fread("median_bray_curtis_sqrt_q1_single_rarefied_to_1000.txt"))
bray.dist.df <- bray.dist.df[!(grepl("LD", bray.dist.df$V1)), !(grepl("LD", names(bray.dist.df)))]
row.names(bray.dist.df) <- bray.dist.df$V1
bray.dist.df$V1 <- NULL

bray.dist <- as.dist(as.matrix(bray.dist.df)) # turning into a distance object

## make our groups
groups <- ifelse(grepl("Air", names(bray.dist.df)), "Air", 
                 ifelse(grepl("Litter", names(bray.dist.df)), "Environmental Litter", 
                        ifelse(grepl("Soil", names(bray.dist.df)), "Soil", 
                               ifelse(grepl("LE", names(bray.dist.df)), "Elevated", 
                                      ifelse(grepl("LU", names(bray.dist.df)), "Overhead", 
                                             ifelse(grepl("LO", names(bray.dist.df)), "Open",
                                                    ifelse(grepl("LC", names(bray.dist.df)), "Closed", "there was a problem..")))))))

## do the test! 
disper <- betadisper(bray.dist, group = groups, type = c("centroid"))
permutest(disper, pairwise = TRUE)

## make a dataframe of distances of each sample within our factors!
distcen.df <- data.frame("Treatments" = groups, "Distance_to_Centroid" = disper$distances)

## make the graph
distcen.df$Treatments <- factor(distcen.df$Treatments,levels = c("Air", "Environmental Litter", "Soil",
                                                                 "Elevated", "Overhead", "Open") )

cols <- c("Air" = "#dbbd28", 
          "Environmental Litter" = "#662f90",
          "Soil" = "#9b2a31", 
          "Elevated" = "#dbbd28", 
          "Overhead" = "#662f90", 
          "Open" = "#9b2a31")

B <- ggplot(data = distcen.df) + 
  geom_boxplot_pattern(aes(x = Treatments, 
                           y = Distance_to_Centroid, 
                           fill = Treatments, 
                           pattern_density = Treatments), 
                       color = "black",
                       pattern_fill = "white", 
                       pattern_spacing = 0.05) + 
  theme_classic() + 
  labs(y = "Distance to Centroid", x = "Treatments") + 
  scale_fill_manual(values = cols) + 
  scale_pattern_density_manual(values = c(0.5, 0.5, 0.5, 0, 0, 0)) + 
  guides(fill = FALSE) +
  theme(legend.position = "none") +
  scale_x_discrete(labels = treatment_labels) + 
  scale_y_continuous(limits = c(0.05, 0.75))



##### __________________ Grass-B Shannon = C ________________ #########
table.rarefied <- data.frame(fread("OTU_filtered_plus_rarefied__ROUNDED.tsv"))
row.names(table.rarefied) <- table.rarefied$V1
table.rarefied$V1 <- NULL

# Calculate alpha diversity!
shannon2 <- as.data.frame(diversity(table.rarefied, index = "shannon"))
names(shannon2) <- c("Shannon")
richness2 <- as.data.frame(apply(table.rarefied[,]>0, 1, sum))
names(richness2) <- c("Richness")
simpson2 <- as.data.frame(diversity(table.rarefied, index = "simpson"))
names(simpson2) <- c("Simpson")

# And merge them all together!
tmp12 <- merge(richness2, shannon2, by.x = "row.names", by.y = "row.names")
all.alpha2 <- merge(tmp12, simpson2, by.x = "Row.names", by.y = "row.names")
names(all.alpha2) <- c("SampleID", names(all.alpha2)[2:length(names(all.alpha2))])

#Merge with metadata to create a plot.
metadata2 <- as.data.frame(fread("metadata_q2_cleaned.txt"))
merged_alpha2 <- merge(all.alpha2, metadata2, by.x = "SampleID", by.y = "SampleID")


# LOOKING AT DISPERSAL ROUTE:
merged_alpha2$Dispersal_Route <- str_replace(merged_alpha2$Dispersal_Route, "Up", "Overhead")
merged_alpha2$Dispersal_Route <- factor(merged_alpha2$Dispersal_Route, levels = c("Environmental Litter", "Closed", "Elevated", "Overhead", "Open"))

cols <- c("Air" = "#dbbd28", 
          "Environmental Litter" = "#662f90",
          "Soil" = "#9b2a31", 
          "Closed" = "#474748",
          "Elevated" = "#dbbd28", 
          "Overhead" = "#662f90", 
          "Open" = "#9b2a31")



C <- ggplot(data = merged_alpha2) +
  geom_boxplot_pattern(aes(x = Dispersal_Route, 
                           y = Shannon,
                           fill = Dispersal_Route, 
                           pattern_density = Dispersal_Route), 
                       color = "black", 
                       outlier.shape = NA, 
                       pattern_fill = "white", 
                       pattern_spacing = 0.05) +
  labs(x = 'Treatment', y = 'Shannon Diversity Index') + 
  theme_classic() + 
  scale_fill_manual(values = cols) +
  scale_pattern_density_manual(values = c(0.5, 0, 0, 0, 0)) + 
  theme(legend.position = "none") +
  guides(fill = FALSE) +
  scale_x_discrete(labels = treatment_labels) + 
  scale_y_continuous(limits = c(0, 6))



##### __________________ Grass-B BC = D _____________________ #####
## Load data
bray.dist.df <- as.data.frame(fread("median_bray_curtis_sqrt_q2_single_rarefied_to_1000.txt"))
row.names(bray.dist.df) <- bray.dist.df$V1
bray.dist.df$V1 <- NULL

bray.dist <- as.dist(as.matrix(bray.dist.df))

## make our groups
groups <- ifelse(grepl("Air", names(bray.dist.df)), "Air", 
                 ifelse(grepl("Litter", names(bray.dist.df)), "Environmental Litter", 
                        ifelse(grepl("Soil", names(bray.dist.df)), "Soil", 
                               ifelse(grepl("RE", names(bray.dist.df)), "Elevated", 
                                      ifelse(grepl("RU", names(bray.dist.df)), "Overhead", 
                                             ifelse(grepl("RO", names(bray.dist.df)), "Open",
                                                    ifelse(grepl("RC", names(bray.dist.df)), "Closed", "there was a problem..")))))))

## do the test! 
disper <- betadisper(bray.dist, group = groups, type = c("centroid"))
permutest(disper, pairwise = TRUE)

## make a dataframe of distances of each sample within our factors!
distcen.df <- data.frame("Treatments" = groups, "Distance_to_Centroid" = disper$distances)

## make the graph
distcen.df$Treatments <- factor(distcen.df$Treatments,levels = c("Environmental Litter", "Closed", 
                                                                 "Elevated", "Overhead", "Open") )

cols <- c("Air" = "#dbbd28", 
          "Environmental Litter" = "#662f90",
          "Soil" = "#9b2a31", 
          "Closed" = "#474748",
          "Elevated" = "#dbbd28", 
          "Overhead" = "#662f90", 
          "Open" = "#9b2a31")

D <- ggplot(data = distcen.df) + 
  geom_boxplot_pattern(aes(x = Treatments, 
                           y = Distance_to_Centroid,
                           fill = Treatments, 
                           pattern_density = Treatments), 
                       color = "black", 
                       pattern_fill = "white", 
                       pattern_spacing = 0.05) + 
  theme_classic() + 
  labs(y = "Distance to Centroid", x = "Treatments") + 
  scale_fill_manual(values = cols) + 
  scale_pattern_density_manual(values = c(0.5, 0, 0, 0, 0)) + 
  guides(fill = FALSE) +
  theme(legend.position = "none") +
  scale_x_discrete(labels = treatment_labels) + 
  scale_y_continuous(limits = c(0.05, 0.75))


##### __________________ Grass-F Shannon = E ________________ #########
table.rarefied <- data.frame(fread("OTU_filtered_plus_rarefied_3500__ROUND_IS_TRUE.tsv"))
row.names(table.rarefied) <- table.rarefied$V1
table.rarefied$V1 <- NULL

# Calculate alpha diversity!
shannon2 <- as.data.frame(diversity(table.rarefied, index = "shannon"))
names(shannon2) <- c("Shannon")
richness2 <- as.data.frame(apply(table.rarefied[,]>0, 1, sum))
names(richness2) <- c("Richness")
simpson2 <- as.data.frame(diversity(table.rarefied, index = "simpson"))
names(simpson2) <- c("Simpson")

# And merge them all together!
tmp12 <- merge(richness2, shannon2, by.x = "row.names", by.y = "row.names")
all.alpha2 <- merge(tmp12, simpson2, by.x = "Row.names", by.y = "row.names")
names(all.alpha2) <- c("SampleID", names(all.alpha2)[2:length(names(all.alpha2))])

#Merge with metadata to create a plot.
metadata2 <- as.data.frame(fread("metadata_q2_ITS_cleaned.txt"))
merged_alpha2 <- merge(all.alpha2, metadata2, by.x = "SampleID", by.y = "SampleID")


###______Plotting alpha diversity:
# LOOKING AT DISPERSAL ROUTE:
routes_keep <- c("Environmental Litter", "Closed", "Elevated", "Up", "Open")
merged_alpha2 <- merged_alpha2[merged_alpha2$Dispersal_Route %in% routes_keep, ]
merged_alpha2$Dispersal_Route <- str_replace(merged_alpha2$Dispersal_Route, "Up", "Overhead")
merged_alpha2$Dispersal_Route <- factor(merged_alpha2$Dispersal_Route, levels = c("Phyllosphere", "Environmental Litter", "Closed", "Elevated", "Overhead", "Open"))

cols <- c("Air" = "#dbbd28", 
          "Environmental Litter" = "#662f90",
          "Soil" = "#9b2a31", 
          "Closed" = "#474748",
          "Elevated" = "#dbbd28", 
          "Overhead" = "#662f90", 
          "Open" = "#9b2a31")


E <- ggplot(data = merged_alpha2) +
  geom_boxplot_pattern(aes(x = Dispersal_Route, 
                           y = Shannon, 
                           fill = Dispersal_Route, 
                           pattern_density = Dispersal_Route), 
                       color = "black", 
                       outlier.shape = NA, 
                       pattern_fill = "white", 
                       pattern_spacing = 0.05) +
  labs(x = "Treatment", y = 'Shannon Diversity Index') + 
  theme_classic() + 
  scale_fill_manual(values = cols) + 
  scale_pattern_density_manual(values = c(0.5, 0, 0, 0, 0)) + 
  theme(legend.position = "none") +
  guides(fill = FALSE) +
  scale_x_discrete(labels = treatment_labels) + 
  scale_y_continuous(limits = c(0, 6))



##### __________________ Grass-F BC = F _____________________ #####
## Load data
bray.dist.df <- as.data.frame(fread("median_bray_curtis_sqrt_q2_ITS_single_rarefied_to_3500.txt"))
row.names(bray.dist.df) <- bray.dist.df$V1 
bray.dist.df$V1 <- NULL

bray.dist <- as.dist(as.matrix(bray.dist.df))

## make our groups
groups <- ifelse(grepl("Litter", names(bray.dist.df)), "Environmental Litter", 
                               ifelse(grepl("RE", names(bray.dist.df)), "Elevated", 
                                      ifelse(grepl("RU", names(bray.dist.df)), "Overhead", 
                                             ifelse(grepl("RO", names(bray.dist.df)), "Open",
                                                    ifelse(grepl("RC", names(bray.dist.df)), "Closed", "there was a problem..")))))

## do the test! 
disper.fun <- betadisper(bray.dist, group = groups, type = c("centroid"))
permutest(disper.fun, permutations = 50,  pairwise = TRUE)

## make a dataframe of distances of each sample within our factors!
distcen.df <- data.frame("Treatments" = groups, "Distance_to_Centroid" = disper.fun$distances)

## make the graph
distcen.df$Treatments <- factor(distcen.df$Treatments,levels = c("Environmental Litter", "Closed", 
                                                                 "Elevated", "Overhead", "Open") )


cols <- c("Air" = "#dbbd28", 
          "Environmental Litter" = "#662f90",
          "Soil" = "#9b2a31", 
          "Closed" = "#474748",
          "Elevated" = "#dbbd28", 
          "Overhead" = "#662f90", 
          "Open" = "#9b2a31")

F <- ggplot(data = distcen.df) + 
  geom_boxplot_pattern(aes(x = Treatments, 
                           y = Distance_to_Centroid, 
                           fill = Treatments, 
                           pattern_density = Treatments), 
                       color = "black", 
                       pattern_fill = "white", 
                       pattern_spacing = 0.05) + 
  theme_classic() + 
  labs(y = "Distance to Centroid", x = "Treatments") + 
  scale_fill_manual(values = cols) + 
  scale_pattern_density_manual(values = c(0.5, 0, 0, 0, 0)) + 
  theme(legend.position = "none") +
  guides(fill = FALSE) +
  scale_x_discrete(labels = treatment_labels) + 
  scale_y_continuous(limits = c(0.05, 0.75))



#####

##### __________________ Put them altogether ________________ ######
grid.arrange(A, C, E, B, D, F, nrow = 2)