# Box plot of  MASS LOSS, taking out the elevated treatment
# Kendra E. Walters
# Nov 17th 2020

# Room of ()ments
require(data.table)
require(ggplot2)

# Input data
mass <- as.data.frame(fread("mass_loss.csv"))

# Set up treatments: 
mass$Treatment <- ifelse(grepl(pattern = "Phyllosphere", x = mass$`Bag Label`), "Phyllosphere", 
                         ifelse(grepl("RC", mass$`Bag Label`), "Closed", 
                                ifelse(grepl("RO", mass$`Bag Label`), "Open", 
                                       ifelse(grepl("RU", mass$`Bag Label`), "Overhead", 
                                              ifelse(grepl("RE", mass$`Bag Label`), "Elevated", 
                                                     ifelse(grepl("Litter", mass$`Bag Label`), "Environmental Litter", "unknown"))))))

mass$Month <- ifelse(mass$`Time Point` == 1, "May", 
                     ifelse(mass$`Time Point` == 2, "June", 
                            ifelse(mass$`Time Point` == 3, "July",
                                   ifelse(mass$`Time Point` == 4, "September", 
                                          ifelse(mass$`Time Point` == 5, "October", "error... :(")))))
mass <- mass[mass$Exclude == "No", ]


# Graph time!
not.elevated <- mass[mass$Treatment != "Elevated", ]

not.elevated$Month <- factor(not.elevated$Month, levels = c("May", "June",
                                                    "July", "September", "October"))
not.elevated$Treatment <- factor(not.elevated$Treatment, levels = c("Closed", "Overhead", "Open"))
cols <- c("Closed" = "#474748",
          "Elevated" = "#dbbd28", 
          "Overhead" = "#662f90", 
          "Open" = "#9b2a31")

ggplot(not.elevated, aes(x = Month, y = `Mass Loss (g)` * 100, fill = Treatment)) + 
  geom_boxplot(color = "black") + theme_classic() +
  ylab("Percent Mass Loss") + 
  theme(axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black")) +
  scale_fill_manual(values = cols)