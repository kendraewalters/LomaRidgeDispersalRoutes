# Script specifically for Figure increase in community similarity of litterbags to environmental litter
# Kendra E. Walters
# September 23rd, 2020

# Room of ()ments
require(data.table)
require(ggplot2)

# Load data + remove unwanted samples
distances <- as.data.frame(fread("distances_each_sample_to_env_by_timepoint.txt"))

# Set up metadata for visualization
months.key <- c("T1" = "May", "T2" = "June", "T3" = "July", "T4" = "September", "T5" = "October")
distances$Month <- months.key[str_extract(distances$Sample, "T[1-5]")]

treatment.key <- c("RC" = "Closed", "RE" = "Elevated", "RU" = "Overhead", "RO" = "Open")
distances$Treatment <- treatment.key[str_extract(distances$Sample, "R[EOUC]")]

# Create graph
cols <- c("Closed" = "#474748",
          "Elevated" = "#dbbd28", 
          "Overhead" = "#662f90", 
          "Open" = "#9b2a31")
distances$Treatment <- factor(distances$Treatment, 
                                levels = c("Closed", "Elevated", "Overhead", "Open"))
distances$Month <- factor(distances$Month, 
                          levels = c("May", "June", "July", "September", "October"))

ggplot(data = distances, aes(x = Month, y = 1-Distance_to_Env_Centroid, fill = Treatment)) + 
  geom_boxplot(color = "black") + 
  theme_classic() + 
  scale_fill_manual(values = cols) +
  ylab("Similarity to Environmental Litter") +
  theme(axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black")) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.55))