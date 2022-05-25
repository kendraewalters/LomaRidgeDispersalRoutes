# April 5 2021
# Kendra E Walters

# Room of ()ments
require(data.table)
require(ggplot2)
require(ggpattern)

# Load the data
grass <- as.data.frame(fread("grass_litterbag_abundance_all_timepoints.csv"))

# Treatment names
grass$Treatment <- ifelse(grepl(pattern = "Phyllosphere", x = grass$`Bag Label`), "Phyllosphere", 
                          ifelse(grepl("RC", grass$`Bag Label`), "Closed", 
                                 ifelse(grepl("RO", grass$`Bag Label`), "Open", 
                                        ifelse(grepl("RU", grass$`Bag Label`), "Overhead", 
                                               ifelse(grepl("RE", grass$`Bag Label`), "Elevated", 
                                                      ifelse(grepl("Litter", grass$`Bag Label`), "Environmental Litter", "unknown"))))))
grass$Month <- ifelse(grass$`Time Point` == 1, "May", 
                      ifelse(grass$`Time Point` == 2, "June", 
                             ifelse(grass$`Time Point` == 3, "July",
                                    ifelse(grass$`Time Point` == 4, "September", 
                                           ifelse(grass$`Time Point` == 5, "October", 
                                                  ifelse(grepl(pattern = "Phyllosphere", x = grass$`Bag Label`), "Before Dispersal", "error... :("))))))

# Exclude treatments not wanted
t_not_0 <- grass[grass$`Time Point` != 0, ]

# Set up factors
t_not_0$Treatment <- factor(t_not_0$Treatment, levels = c("Environmental Litter", "Closed", 
                                                          "Elevated", "Overhead", 
                                                          "Open"))
t_not_0$Month <- factor(t_not_0$Month, levels = c("May", "June", "July", "September", "October"))
# Create plot
cols <- c("Environmental Litter" = "#662f90",
          "Closed" = "#474748", 
          "Elevated" = "#dbbd28",
          "Overhead" = "#662f90",
          "Open" = "#9b2a31")

ggplot(t_not_0, aes(x = Month, y = `Abundance Per G`)) + 
  geom_boxplot_pattern(aes(fill = Treatment, pattern_density = Treatment), 
               color = "black", 
               pattern = "stripe", 
               pattern_fill = "white", 
               pattern_spacing = 0.022) + 
  theme_classic() +
  theme(legend.key.size = unit(0.8, 'cm')) +
  scale_fill_manual(values = cols) + 
  scale_pattern_density_manual(values = c(0.5, 0, 0, 0, 0)) + 
  theme(axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black")) + 
  ylab("Bacterial Abundance / g(dw)") 