# April 5 2021
# Kendra E Walters

# Room of ()ments
require(data.table)
require(ggplot2)
require(ggpattern)
require(tidyverse)

# Load the data
remove <- c("LUT4R1", "LOT4R1", "LOT5R4") # Remove LUT4R1 (252), LOT4R1 (257), and LOT5R4 (326), two of which were literally COVERED in soil

glass <- as.data.frame(fread("glass_slides_abundance_all_timepoints.csv")) %>% 
  mutate(Treatment = ifelse(grepl("LC", `Bag Label`), "Closed", 
                             ifelse(grepl("LO", `Bag Label`), "Open", 
                                    ifelse(grepl("LU", `Bag Label`), "Overhead", 
                                           ifelse(grepl("LE", `Bag Label`), "Elevated", "unknown")))), 
         Month = ifelse(`Time Point` == 1, "May", 
                        ifelse(`Time Point` == 2, "June", 
                               ifelse(`Time Point` == 3, "July",
                                      ifelse(`Time Point` == 4, "September", 
                                             ifelse(`Time Point` == 5, "October", "error... :("))))), 
         Abundance_per_cm2 = Abundance/(2.5*7.5), 
         Month = factor(Month, levels = c("May", "June", "July", "September", "October")), 
         Treatment = factor(Treatment, levels = c("Closed", "Elevated", "Overhead", "Open"))) %>% 
  filter(!(`Bag Label` %in% remove))


# Create plot
cols <- c("Closed" = "#474748", 
          "Elevated" = "#dbbd28",
          "Overhead" = "#662f90",
          "Open" = "#9b2a31")

ggplot(glass, aes(x = Month, y = Abundance_per_cm2)) + 
  geom_boxplot(aes(fill = Treatment), color = "black") + 
  theme_classic() +
  scale_fill_manual(values = cols) + 
  theme(axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black")) + 
  ylab("Bacterial Abundance / cm2")