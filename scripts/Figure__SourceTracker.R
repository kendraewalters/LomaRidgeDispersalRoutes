# SourceTracker graph
# November 2019
# Kendra E. Walters

# load
require(data.table)
require(ggplot2)
require(tidyr)
require(ggpattern)

# Setting up
proportions <- as.data.frame(fread("mixing_proportions.txt"))
names(proportions) <- c("SampleID", "Air", "Environmental Litter", "Soil", "Unknown")


# Creating metadata
proportions$Dispersal_Route <- ifelse(grepl("^LO", proportions$SampleID), "Open", 
                                                              ifelse(grepl("^LC", proportions$SampleID), "Closed", 
                                                                     ifelse(grepl("^LU", proportions$SampleID), "Overhead", 
                                                                            ifelse(grepl("^LE", proportions$SampleID), "Elevated",
                                                                                   ifelse(grepl("Litter", proportions$SampleID), "Environmental Litter", 
                                                                                          ifelse(grepl("LD", proportions$SampleID), "Death",
                                                                                                 ifelse(grepl("Soil", proportions$SampleID), "Soil",
                                                                                                        ifelse(grepl("Air", proportions$SampleID), "Air", "ooops....."))))))))

# Formatting dataframe for ggplot
proportions.long <- gather(proportions, Source, Proportion, Air:Unknown, factor_key=TRUE)

# Creating the amazing boxplot
cols <- c("Air" = "#dbbd28", 
          "Environmental Litter" = "#662f90",
          "Soil" = "#9b2a31", 
          "Unknown" = "#474748")

proportions.long$Dispersal_Route <- factor(proportions.long$Dispersal_Route, 
                                              levels = c("Elevated", "Overhead", "Open"))

ggplot(data = proportions.long, aes(x = Dispersal_Route, y = Proportion)) + 
  geom_boxplot_pattern(aes(pattern_fill = Source, pattern_density = Source), 
                       pattern = 'stripe',
                       colour  = 'black', 
                       pattern_spacing = 0.04) + 
  theme_classic() + 
  xlab("Treatment") +
  theme(legend.key.size = unit(0.8, 'cm')) + 
  scale_pattern_fill_manual(values = cols) + 
  scale_pattern_density_manual(values = c(Air = 0.4, `Environmental Litter` = 0.4, Soil = 0.4, Unknown = 0.0))