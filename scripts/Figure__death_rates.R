# Figure death rate measurements
# April 7 2021
# Kendra E Walters

## Room of ()ments
require(data.table)
require(ggplot2)
require(gridExtra)

#####________ Panel A = death rate over 6 months __________ #############
death <- as.data.frame(fread("death_slides_set1_abundance_all_timepoints.csv"))

# Making a new column to represent # of days, not timepoints
days <- c(`0` = 0, `1` = 39, `2` = 60, `3` = 100, `4` = 151, `5` = 195)
death$Time_as_days <- days[as.character(death$`Time Point`)]

# Calculate intercept + slope for the graph
intercept <- mean(log(death[death$Time_as_days == 0, ]$Abundance))

lm_death <- lm(I(log(Abundance) - intercept) ~ Time_as_days + 0, death)
summary(lm_death) # estimate: -0.014771

# Graph 
A <-  ggplot() + 
        geom_point(data = death, aes(x = Time_as_days, y = log(Abundance))) + 
        theme_classic() +
        xlab("Days") + 
        ylab("ln(Abundance)") + 
        theme(axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black")) +
        geom_abline(intercept = intercept, slope = coef(lm_death), linetype = 2)


#####________ Panel B = death rate table v ground _________ #####
new_death <- as.data.frame(fread("death_slides_set2_abundance_all_timepoints_FOR_R.csv"))

# Making a new column to represent # of days, not timepoints
new.days <- c(`0` = 0, `1` = 7, `2` = 14, `3` = 21, `4` = 28, `5` = 44)
new_death$Time_as_days <- new.days[as.character(new_death$`Time Point`)]

# Calculate intercept + slope for the graph
new_death$Treatment <- ifelse(grepl(pattern = "NDU", x = new_death$`Bag Label`), "Table", 
                              ifelse(grepl(pattern = "morning", x = new_death$`Bag Label`), "Initial Abundance", 
                                     ifelse(grepl(pattern = "night", x = new_death$`Bag Label`), "Night", "Ground")))
new_table <- new_death[new_death$Treatment == "Initial Abundance", ]
new_ground <- new_death[new_death$Treatment == "Initial Abundance", ]
new_table$Treatment <- c(rep("Table", nrow(new_table)))
new_ground$Treatment <- c(rep("Ground", nrow(new_ground)))

new_table <- rbind(new_table, new_death[new_death$Treatment == "Table", ])
new_ground <- rbind(new_ground, new_death[new_death$Treatment == "Ground", ])


intercept <- mean(log(new_death[new_death$Time_as_days == 0, ]$Abundance))

lm.table <- lm(I(log(Abundance) - intercept) ~ Time_as_days + 0, new_table)
lm.ground <- lm(I(log(Abundance) - intercept) ~ Time_as_days + 0, new_ground)

cols <- c("Table" = "black", "Ground" = "grey40")

B <- ggplot() + 
      geom_point(data = new_table, aes(x = Time_as_days, y = log(Abundance), color = "Table")) + 
      geom_point(data = new_ground, aes(x = Time_as_days, y = log(Abundance), color = "Ground")) +
      scale_colour_manual(name="Treatment", values=cols) +
      geom_abline(intercept = intercept, slope = coef(lm.table), linetype = 2, color = "black") +
      geom_abline(intercept = intercept, slope = coef(lm.ground), linetype = 2, color = "grey40") + 
      xlab("Days") +
      ylab("ln(Abundance)") +
      theme_classic() + 
      theme(axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black"))




#####________ Combine panels
grid.arrange(A, B, nrow = 1)