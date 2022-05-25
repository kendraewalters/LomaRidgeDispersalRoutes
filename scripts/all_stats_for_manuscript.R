# Statistics for dispersal manuscript
# August 8th 2020
# Kendra E Walters


# Room of ()ments
require(data.table)
require(vegan)
require(car)
require(broom)
require(FSA)
require(usedist)
require(multcomp)
require(tidyr)
require(stringr)
require(emmeans)
require(ggpubr)
require(anchors)
require(tidyverse)

##### _______ subtracted abundance on glass slides ____________ ########
glass <- as.data.frame(fread("glass_slides_abundance_all_timepoints.csv"))

## Make treatment names:
glass$Treatment <- ifelse(grepl("LC", glass$`Bag Label`), "No Dispersal", 
                          ifelse(grepl("LO", glass$`Bag Label`), "All Dispersal", 
                                 ifelse(grepl("LU", glass$`Bag Label`), "Local Rain/Air", 
                                        ifelse(grepl("LE", glass$`Bag Label`), "Regional Rain/Air", "unknown"))))

glass$Month <- ifelse(glass$`Time Point` == 1, "End May", 
                      ifelse(glass$`Time Point` == 2, "End June", 
                             ifelse(glass$`Time Point` == 3, "End July",
                                    ifelse(glass$`Time Point` == 4, "Mid September", 
                                           ifelse(glass$`Time Point` == 5, "End October", "error... :(")))))

glass$Abundance_per_cm2 <- glass$`Abundance`/(2.5*7.5)
glass$Month <- factor(glass$Month, levels = c("End May", "End June", "End July", "Mid September", "End October"), ordered = TRUE)

## Remove LUT4R1 (252), LOT4R1 (257), and LOT5R4 (326), two of which were literally COVERED in soil
remove <- c("LUT4R1", "LOT4R1", "LOT5R4")
glass_red <- glass[!(glass$`Bag Label` %in% remove), ]


## Make a difference between group graph but use the paired samples within a plot
glass_red$Replicate <- gsub("L[A-Z]T[1-5]", "", glass_red$`Bag Label`)
diff.month <- list()
diff <- list()
months <- unique(glass_red$Month)

out.df <- data.frame("Abundance_per_cm2" = numeric(), "Treatment" = character(), 
                     "Replicate" = character(), "Month" = character())

for (i in 1:length(months)) {
  key <- as.character(months[[i]])
  sub <- glass_red[glass_red$Month == key, ]
  rep <- unique(sub$Replicate)
  
  for (j in 1:length(rep)) {
    reprep <- rep[[j]]
    
    subsub <- sub[sub$Replicate == reprep, ]
    regional <- subsub[subsub$Treatment == "Regional Rain/Air", ]$Abundance_per_cm2 - subsub[subsub$Treatment == "No Dispersal", ]$Abundance_per_cm2
    local <- subsub[subsub$Treatment == "Local Rain/Air", ]$Abundance_per_cm2 - subsub[subsub$Treatment == "Regional Rain/Air", ]$Abundance_per_cm2
    soil <- subsub[subsub$Treatment == "All Dispersal", ]$Abundance_per_cm2 - subsub[subsub$Treatment == "Local Rain/Air", ]$Abundance_per_cm2
    
    out <- data.frame(c("Soil" = soil, "Vegetation" = local, "Air" = regional))
    out$Treatment <- row.names(out)
    out$Replicate <- c(rep(reprep, nrow(out)))
    out$Month <- c(rep(key, nrow(out)))
    names(out) <- c("Abundance_per_cm2", "Treatment", "Replicate", "Month")
    
    out.df <- rbind(out.df, out)
  }
}

## Statistics
# Overall model statistic for Figure 2A
mod <- lm(Abundance_per_cm2 ~ Treatment, out.df) # Type III SS, ANOVA, treatment
Anova(mod, type = "III")

# Are the means greater than 0? for Figure 2A
library(ggpubr)
hist(out.df$Abundance_per_cm2)
ggdensity(out.df$Abundance_per_cm2)
ggqqplot(out.df$Abundance_per_cm2)
shapiro.test(out.df$Abundance_per_cm2)

t.test(out.df[out.df$Treatment == "Soil", ]$Abundance_per_cm2)
t.test(out.df[out.df$Treatment == "Air", ]$Abundance_per_cm2)
t.test(out.df[out.df$Treatment == "Vegetation", ]$Abundance_per_cm2)

# Overall model statistic for Figure 2B
options(contrasts = c("contr.sum","contr.poly"))

mod <- lm(Abundance_per_cm2 ~ Treatment * Month, out.df) # Type III SS, ANOVA, treatment by month interaction
Anova(mod, type = "III")
model <- Anova(mod, type = "III")

tidy_aov <- tidy(model)

sum_squares_treat_x_month <- tidy_aov[tidy_aov$term == "Treatment:Month", ]$sumsq
sum_squares_residuals <- tidy_aov[tidy_aov$term == "Residuals", ]$sumsq

R_squared <- sum_squares_treat_x_month /
  (sum_squares_treat_x_month + sum_squares_residuals)

R_squared



# Testing models per timepoint and whether the means greater than 0? for Figure 2B
t1 <- out.df[out.df$Month == "End May", ]
mod <- lm(Abundance_per_cm2 ~ Treatment, t1)
Anova(mod, type = "III") # P = 0.4062

t.test(t1[t1$Treatment == "Air", ]$Abundance_per_cm2) # 0.0005647
t.test(t1[t1$Treatment == "Vegetation", ]$Abundance_per_cm2) # 0.02806
t.test(t1[t1$Treatment == "Soil", ]$Abundance_per_cm2)


t2 <- out.df[out.df$Month == "End June", ]
mod <- lm(Abundance_per_cm2 ~ Treatment, t2)
Anova(mod, type = "III") # P = 0.1249

t.test(t2[t2$Treatment == "Air", ]$Abundance_per_cm2) # 0.0001676
t.test(t2[t2$Treatment == "Vegetation", ]$Abundance_per_cm2) # 0.04949
t.test(t2[t2$Treatment == "Soil", ]$Abundance_per_cm2)


t3 <- out.df[out.df$Month == "End July", ]
mod <- lm(Abundance_per_cm2 ~ Treatment, t3)
Anova(mod, type = "III") # P = 0.017973
TukeyHSD(aov(Abundance_per_cm2 ~ Treatment, t3))


t.test(t3[t3$Treatment == "Air", ]$Abundance_per_cm2) # 0.0001504
t.test(t3[t3$Treatment == "Vegetation", ]$Abundance_per_cm2)
t.test(t3[t3$Treatment == "Soil", ]$Abundance_per_cm2)


t4 <- out.df[out.df$Month == "Mid September", ]
mod <- lm(Abundance_per_cm2 ~ Treatment, t4)
Anova(mod, type = "III") # P = 0.02218
TukeyHSD(aov(Abundance_per_cm2 ~ Treatment, t4))

t.test(t4[t4$Treatment == "Air", ]$Abundance_per_cm2) # 0.009686
t.test(t4[t4$Treatment == "Vegetation", ]$Abundance_per_cm2)
t.test(t4[t4$Treatment == "Soil", ]$Abundance_per_cm2)


t5 <- out.df[out.df$Month == "End October", ]
mod <- lm(Abundance_per_cm2 ~ Treatment, t5)
Anova(mod, type = "III") # P = 0.0224
TukeyHSD(aov(Abundance_per_cm2 ~ Treatment, t5))

t.test(t5[t5$Treatment == "Air", ]$Abundance_per_cm2) # 0.8437
t.test(t5[t5$Treatment == "Vegetation", ]$Abundance_per_cm2) # 0.006704
t.test(t5[t5$Treatment == "Soil", ]$Abundance_per_cm2) # 0.0151


##### _______ raw abundance on glass slides _______________ #####
# Room of ()ments
require(data.table)
require(car)
require(ggpattern)
require(tidyverse)

# Load the data
remove <- c("LUT4R1", "LOT4R1", "LOT5R4") # Remove LUT4R1 (252), LOT4R1 (257), and LOT5R4 (326), two of which were literally COVERED in soil

glass <- as.data.frame(fread("glass_slides_abundance_all_timepoints.csv")) %>% 
  mutate(Treatment = ifelse(grepl("LC", glass$`Bag Label`), "Closed", 
                            ifelse(grepl("LO", glass$`Bag Label`), "Open", 
                                   ifelse(grepl("LU", glass$`Bag Label`), "Overhead", 
                                          ifelse(grepl("LE", glass$`Bag Label`), "Elevated", "unknown")))), 
         Month = ifelse(glass$`Time Point` == 1, "May", 
                        ifelse(glass$`Time Point` == 2, "June", 
                               ifelse(glass$`Time Point` == 3, "July",
                                      ifelse(glass$`Time Point` == 4, "September", 
                                             ifelse(glass$`Time Point` == 5, "October", "error... :("))))), 
         Abundance_per_cm2 = glass$`Abundance`/(2.5*7.5), 
         Month = factor(Month, levels = c("May", "June", "July", "September", "October")), 
         Treatment = factor(Treatment, levels = c("Closed", "Elevated", "Overhead", "Open"))) %>% 
  filter(!(`Bag Label` %in% remove))


## Stats time!!
# Overall model statistic
library(ggpubr)
hist(out.df$Abundance_per_cm2)
ggdensity(out.df$Abundance_per_cm2)
ggqqplot(out.df$Abundance_per_cm2)
shapiro.test(out.df$Abundance_per_cm2)



options(contrasts = c("contr.sum","contr.poly"))

mod <- lm(Abundance_per_cm2 ~ Treatment * Month, glass) # Type III SS, ANOVA, treatment by month interaction
Anova(mod, type = "III")
model <- Anova(mod, type = "III")

library(broom)
tidy_aov <- tidy(model)

sum_squares_treat <- tidy_aov[tidy_aov$term == "Treatment", ]$sumsq
sum_squares_month <- tidy_aov[tidy_aov$term == "Month", ]$sumsq

sum_squares_treat_x_month <- tidy_aov[tidy_aov$term == "Treatment:Month", ]$sumsq
sum_squares_residuals <- tidy_aov[tidy_aov$term == "Residuals", ]$sumsq

all <- sum_squares_treat + sum_squares_month + sum_squares_treat_x_month + sum_squares_residuals

sum_squares_treat / all 
sum_squares_month / all 
sum_squares_treat_x_month / all 




# Testing models per timepoint and whether the means greater than 0? 
t1 <- glass[glass$Month == "May", ] 
t1 %>% 
  lm(Abundance_per_cm2 ~ Treatment, .) %>% 
  Anova(., type = "III") 
TukeyHSD(aov(Abundance_per_cm2 ~ Treatment, t1))

t2 <- glass[glass$Month == "June", ] 
t2 %>% 
  lm(Abundance_per_cm2 ~ Treatment, .) %>% 
  Anova(., type = "III") 
TukeyHSD(aov(Abundance_per_cm2 ~ Treatment, t2))

t3 <- glass[glass$Month == "July", ] 
t3 %>% 
  lm(Abundance_per_cm2 ~ Treatment, .) %>% 
  Anova(., type = "III") 
TukeyHSD(aov(Abundance_per_cm2 ~ Treatment, t3))


t4 <- glass[glass$Month == "September", ] 
t4 %>% 
  lm(Abundance_per_cm2 ~ Treatment, .) %>% 
  Anova(., type = "III") 
TukeyHSD(aov(Abundance_per_cm2 ~ Treatment, t4))

t5 <- glass[glass$Month == "October", ] 
t5 %>% 
  lm(Abundance_per_cm2 ~ Treatment, .) %>% 
  Anova(., type = "III") 
TukeyHSD(aov(Abundance_per_cm2 ~ Treatment, t5))


##### _______ dispersal from soil over time _________ #####
glass <- as.data.frame(fread("glass_slide_abundance_all_timepoints.csv"))

## Make treatment names:
glass$Treatment <- ifelse(grepl("LC", glass$`Bag Label`), "No Dispersal", 
                          ifelse(grepl("LO", glass$`Bag Label`), "All Dispersal", 
                                 ifelse(grepl("LU", glass$`Bag Label`), "Local Rain/Air", 
                                        ifelse(grepl("LE", glass$`Bag Label`), "Regional Rain/Air", "unknown"))))

glass$Month <- ifelse(glass$`Time Point` == 1, "End May", 
                      ifelse(glass$`Time Point` == 2, "End June", 
                             ifelse(glass$`Time Point` == 3, "End July",
                                    ifelse(glass$`Time Point` == 4, "Mid September", 
                                           ifelse(glass$`Time Point` == 5, "End October", "error... :(")))))

glass$Abundance_per_cm2 <- glass$`Abundance`/(2.5*7.5)
glass$Month <- factor(glass$Month, levels = c("End May", "End June", "End July", "Mid September", "End October"), ordered = TRUE)

## Remove LUT4R1 (252), LOT4R1 (257), and LOT5R4 (326), two of which were literally COVERED in soil
remove <- c("LUT4R1", "LOT4R1", "LOT5R4")
glass_red <- glass[!(glass$`Bag Label` %in% remove), ]


## Make a difference between group graph but use the paired samples within a plot
glass_red$Replicate <- gsub("L[A-Z]T[1-5]", "", glass_red$`Bag Label`)
diff.month <- list()
diff <- list()
months <- unique(glass_red$Month)

out.df <- data.frame("Abundance_per_cm2" = numeric(), "Treatment" = character(), 
                     "Replicate" = character(), "Month" = character())

for (i in 1:length(months)) {
  key <- as.character(months[[i]])
  sub <- glass_red[glass_red$Month == key, ]
  rep <- unique(sub$Replicate)
  
  for (j in 1:length(rep)) {
    reprep <- rep[[j]]
    
    subsub <- sub[sub$Replicate == reprep, ]
    regional <- subsub[subsub$Treatment == "Regional Rain/Air", ]$Abundance_per_cm2 - subsub[subsub$Treatment == "No Dispersal", ]$Abundance_per_cm2
    local <- subsub[subsub$Treatment == "Local Rain/Air", ]$Abundance_per_cm2 - subsub[subsub$Treatment == "Regional Rain/Air", ]$Abundance_per_cm2
    soil <- subsub[subsub$Treatment == "All Dispersal", ]$Abundance_per_cm2 - subsub[subsub$Treatment == "Local Rain/Air", ]$Abundance_per_cm2
    
    out <- data.frame(c("Soil" = soil, "Vegetation" = local, "Air" = regional))
    out$Treatment <- row.names(out)
    out$Replicate <- c(rep(reprep, nrow(out)))
    out$Month <- c(rep(key, nrow(out)))
    names(out) <- c("Abundance_per_cm2", "Treatment", "Replicate", "Month")
    
    out.df <- rbind(out.df, out)
  }
}

soil <- out.df[out.df$Treatment == "Soil", ]
soil$is.T5 <- ifelse(soil$Month == "End October", "Yes", "No")

# Stats
mod <- lm(Abundance_per_cm2 ~ Month, data = soil)
summary(aov(mod))
TukeyHSD(aov(mod))

t.test(soil[soil$is.T5 == "Yes", ]$Abundance_per_cm2, 
       soil[soil$is.T5 == "No", ]$Abundance_per_cm2, 
       alternative = c("greater"))


##### _______ abundance on glass slides / abundance of litterbags _______ ######

# Incoming cells 
glass <- as.data.frame(fread("glass_slide_abundance_all_timepoints.csv"))

## Make treatment names:
glass$Treatment <- ifelse(grepl("LC", glass$`Bag Label`), "No Dispersal", 
                          ifelse(grepl("LO", glass$`Bag Label`), "All Dispersal", 
                                 ifelse(grepl("LU", glass$`Bag Label`), "Local Rain/Air", 
                                        ifelse(grepl("LE", glass$`Bag Label`), "Regional Rain/Air", "unknown"))))

glass$Month <- ifelse(glass$`Time Point` == 1, "End May", 
                      ifelse(glass$`Time Point` == 2, "End June", 
                             ifelse(glass$`Time Point` == 3, "End July",
                                    ifelse(glass$`Time Point` == 4, "Mid September", 
                                           ifelse(glass$`Time Point` == 5, "End October", "error... :(")))))

glass$Abundance_per_cm2 <- glass$`Abundance`/(2.5*7.5)
glass$Month <- factor(glass$Month, levels = c("End May", "End June", "End July", "Mid September", "End October"), ordered = TRUE)

## Remove LUT4R1 (252), LOT4R1 (257), and LOT5R4 (326), two of which were literally COVERED in soil
remove <- c("LUT4R1", "LOT4R1", "LOT5R4")
glass_red <- glass[!(glass$`Bag Label` %in% remove), ]

## Just the open treatment -- because you just want the total immigration number
open <- glass_red[glass_red$Treatment == "All Dispersal", ]

## Now the cells in the litterbags 
grass <- as.data.frame(fread("grass_litterbag_abundance_all_timepoints.csv"))

## Make treatment names:
grass$Treatment <- ifelse(grepl(pattern = "Phyllosphere", x = grass$`Bag Label`), "Phyllosphere", 
                          ifelse(grepl("RC", grass$`Bag Label`), "No Dispersal", 
                                 ifelse(grepl("RO", grass$`Bag Label`), "All Dispersal", 
                                        ifelse(grepl("RU", grass$`Bag Label`), "Local Rain/Air", 
                                               ifelse(grepl("RE", grass$`Bag Label`), "Regional Rain/Air", 
                                                      ifelse(grepl("Litter", grass$`Bag Label`), "Environmental Litter", "unknown"))))))
grass$Month <- ifelse(grass$`Time Point` == 1, "End May", 
                      ifelse(grass$`Time Point` == 2, "End June", 
                             ifelse(grass$`Time Point` == 3, "End July",
                                    ifelse(grass$`Time Point` == 4, "Mid September", 
                                           ifelse(grass$`Time Point` == 5, "End October", 
                                                  ifelse(grepl(pattern = "Phyllosphere", x = grass$`Bag Label`), "Before Dispersal", "error... :("))))))

## to be equal, just the cells in the completely open bags
open.litter <- grass[grass$Treatment == "All Dispersal", ]

## Now the MASS of those litterbags
mass <- as.data.frame(fread("dry_weight_litter_bag__for_immigration_percent_calc.txt"))

## Merge the mass of the litterbags with the abundance per gram
open.litter.mass <- merge(open.litter, mass, by.x = "Bag Label", by.y = "Sample", all = FALSE)

## Calculate the number of cells in each litterbag
open.litter.mass$Cells_in_bag <- open.litter.mass$`Abundance Per G` * open.litter.mass$Dry_Mass_g

## Calculate the number of cells in a 1 cm^2 area in litterbags
# here we are assuming a 1cm sealed rim around the bag, reducing it to 8cm x 8cm (instead of the 10cm x 10cm fabric dimension)
open.litter.mass$Litter_cells_per_cm2 <- open.litter.mass$Cells_in_bag / (8*8)

## Now we divide the means
(mean(open$Abundance_per_cm2)) / (mean(open.litter.mass$Litter_cells_per_cm2)) * 100 



##### _______ abundance of litterbags _______########
grass <- as.data.frame(fread("grass_litterbag_abundance_all_timepoints.csv"))


## Make treatment names:
grass$Treatment <- ifelse(grepl(pattern = "Phyllosphere", x = grass$`Bag Label`), "Phyllosphere", 
                          ifelse(grepl("RC", grass$`Bag Label`), "No Dispersal", 
                                 ifelse(grepl("RO", grass$`Bag Label`), "All Dispersal", 
                                        ifelse(grepl("RU", grass$`Bag Label`), "Local Rain/Air", 
                                               ifelse(grepl("RE", grass$`Bag Label`), "Regional Rain/Air", 
                                                      ifelse(grepl("Litter", grass$`Bag Label`), "Environmental Litter", "unknown"))))))
grass$Month <- ifelse(grass$`Time Point` == 1, "End May", 
                      ifelse(grass$`Time Point` == 2, "End June", 
                             ifelse(grass$`Time Point` == 3, "End July",
                                    ifelse(grass$`Time Point` == 4, "Mid September", 
                                           ifelse(grass$`Time Point` == 5, "End October", 
                                                  ifelse(grepl(pattern = "Phyllosphere", x = grass$`Bag Label`), "Before Dispersal", "error... :("))))))

## create differences between treatments at each timepoint 
just_bags <- grass[!(grass$Treatment %in% c("Environmental Litter", "Phyllosphere")), ]

just_bags$Replicate <- gsub("R[A-Z]T[1-5]", "", just_bags$`Bag Label`)
diff.month <- list()
diff <- list()
months <- unique(just_bags$Month)

out.df.grass <- data.frame("Abundance_per_g" = numeric(), "Treatment" = character(), 
                           "Replicate" = character(), "Month" = character())


for (i in 1:length(months)) {
  key <- as.character(months[[i]])
  sub <- just_bags[just_bags$Month == key, ]
  rep <- unique(sub$Replicate)
  
  for (j in 1:length(rep)) {
    reprep <- rep[[j]]
    
    subsub <- sub[sub$Replicate == reprep, ]
    regional <- subsub[subsub$Treatment == "Regional Rain/Air", ]$`Abundance Per G` - subsub[subsub$Treatment == "No Dispersal", ]$`Abundance Per G`
    local <- subsub[subsub$Treatment == "Local Rain/Air", ]$`Abundance Per G` - subsub[subsub$Treatment == "Regional Rain/Air", ]$`Abundance Per G` 
    soil <- subsub[subsub$Treatment == "All Dispersal", ]$`Abundance Per G` - subsub[subsub$Treatment == "Local Rain/Air", ]$`Abundance Per G` 
    
    out <- data.frame(c("Soil" = soil, "Vegetation" = local, "Air" = regional))
    out$Treatment <- row.names(out)
    out$Replicate <- c(rep(reprep, nrow(out)))
    out$Month <- c(rep(key, nrow(out)))
    names(out) <- c("Abundance_per_g", "Treatment", "Replicate", "Month")
    
    out.df.grass <- rbind(out.df.grass, out)
    
  }
}

## Statistics!
# Testing for normality; looks okay to me (esp. for generally robust anovas)

options(contrasts = c("contr.sum","contr.poly"))

library(ggpubr)
hist(out.df.grass$Abundance_per_g)
ggdensity(out.df.grass$Abundance_per_g)
ggqqplot(out.df.grass$Abundance_per_g)
shapiro.test(out.df.grass$Abundance_per_g)

# Overall model statistic
mod <- lm(Abundance_per_g ~ Treatment, out.df.grass) # Type III SS, ANOVA, treatment
Anova(mod, type = "III")
TukeyHSD(aov(mod))

mod <- lm(Abundance_per_g ~ Treatment * Month, out.df.grass)
Anova(mod, type = "III")


# To get the r-squared for the treatment
model <- Anova(mod, type = "III")
tidy_aov <- tidy(model) 
sum_squares_regression <- tidy_aov[tidy_aov$term == "Treatment", ]$sumsq # put the term here you are interested in
sum_squares_residuals <- tidy_aov[tidy_aov$term == "Residuals", ]$sumsq

R_squared <- sum_squares_regression /
  (sum_squares_regression + sum_squares_residuals)

R_squared


# Overall t-test statistics 
t.test(out.df.grass[out.df.grass$Treatment == "Air", ]$Abundance_per_g, alternative = "greater")
t.test(out.df.grass[out.df.grass$Treatment == "Vegetation", ]$Abundance_per_g, alternative = "greater")
t.test(out.df.grass[out.df.grass$Treatment == "Soil", ]$Abundance_per_g, alternative = "greater")

# Testing models per timepoint and whether the means greater than 0?
t1 <- out.df.grass[out.df.grass$Month == "End May", ]
mod <- lm(Abundance_per_g ~ Treatment, t1)
Anova(mod, type = "III") # P = 0.02564
TukeyHSD(aov(Abundance_per_g ~ Treatment, t1))

t.test(t1[t1$Treatment == "Air", ]$Abundance_per_g)
t.test(t1[t1$Treatment == "Vegetation", ]$Abundance_per_g)
t.test(t1[t1$Treatment == "Soil", ]$Abundance_per_g)


t2 <- out.df.grass[out.df.grass$Month == "End June", ]
mod <- lm(Abundance_per_g ~ Treatment, t2)
Anova(mod, type = "III") # P = 0.01738
TukeyHSD(aov(Abundance_per_g ~ Treatment, t2))

t.test(t2[t2$Treatment == "Air", ]$Abundance_per_g)
t.test(t2[t2$Treatment == "Vegetation", ]$Abundance_per_g)
t.test(t2[t2$Treatment == "Soil", ]$Abundance_per_g)


t3 <- out.df.grass[out.df.grass$Month == "End July", ]
mod <- lm(Abundance_per_g ~ Treatment, t3)
Anova(mod, type = "III") # P = 0.01181
TukeyHSD(aov(Abundance_per_g ~ Treatment, t3))

t.test(t3[t3$Treatment == "Air", ]$Abundance_per_g)
t.test(t3[t3$Treatment == "Vegetation", ]$Abundance_per_g)
t.test(t3[t3$Treatment == "Soil", ]$Abundance_per_g)


t4 <- out.df.grass[out.df.grass$Month == "Mid September", ]
mod <- lm(Abundance_per_g ~ Treatment, t4)
Anova(mod, type = "III") # P = 0.04005
TukeyHSD(aov(Abundance_per_g ~ Treatment, t4))

t.test(t4[t4$Treatment == "Air", ]$Abundance_per_g)
t.test(t4[t4$Treatment == "Vegetation", ]$Abundance_per_g)
t.test(t4[t4$Treatment == "Soil", ]$Abundance_per_g)


t5 <- out.df.grass[out.df.grass$Month == "End October", ]
mod <- lm(Abundance_per_g ~ Treatment, t5)
Anova(mod, type = "III") # P = 0.3556
TukeyHSD(aov(Abundance_per_g ~ Treatment, t5))

t.test(t5[t5$Treatment == "Air", ]$Abundance_per_g)
t.test(t5[t5$Treatment == "Vegetation", ]$Abundance_per_g)
t.test(t5[t5$Treatment == "Soil", ]$Abundance_per_g)


## For JUST vegetation, running a lm
out.df.grass$Time_as_days <- ifelse(out.df.grass$Month == "End May", 39, 
                                    ifelse(out.df.grass$Month == "End June", 60, 
                                           ifelse(out.df.grass$Month == "End July", 100, 
                                                  ifelse(out.df.grass$Month == "Mid September", 151, 
                                                         ifelse(out.df.grass$Month == "End October", 195, "EIIIIIEEEEEIIII!!")))))
out.df.grass$Time_as_days <- as.numeric(out.df.grass$Time_as_days)

summary(lm(Abundance_per_g ~ Time_as_days, out.df.grass[out.df.grass$Treatment == "Vegetation", ]))


##### _______ correlation abundance of litterbags and mass loss __________ #####
# Input data
mass <- as.data.frame(fread("mass_loss.csv"))

# Set up treatments: 
mass$Treatment <- ifelse(grepl(pattern = "Phyllosphere", x = mass$`Bag Label`), "Phyllosphere", 
                         ifelse(grepl("RC", mass$`Bag Label`), "No Dispersal", 
                                ifelse(grepl("RO", mass$`Bag Label`), "All Dispersal", 
                                       ifelse(grepl("RU", mass$`Bag Label`), "Local Rain/Air", 
                                              ifelse(grepl("RE", mass$`Bag Label`), "Regional Rain/Air", 
                                                     ifelse(grepl("Litter", mass$`Bag Label`), "Environmental Litter", "unknown"))))))

mass$Month <- ifelse(mass$`Time Point` == 1, "End May", 
                     ifelse(mass$`Time Point` == 2, "End June", 
                            ifelse(mass$`Time Point` == 3, "End July",
                                   ifelse(mass$`Time Point` == 4, "Mid September", 
                                          ifelse(mass$`Time Point` == 5, "End October", "error... :(")))))
mass <- mass[mass$Exclude == "No", ]

not.elevated <- mass[mass$Treatment != "Regional Rain/Air", ]



## Raw abundance on litterbags
grass <- as.data.frame(fread("grass_litterbag_abundance_all_timepoints.csv"))


grass$Treatment <- ifelse(grepl(pattern = "Phyllosphere", x = grass$`Bag Label`), "Phyllosphere", 
                          ifelse(grepl("RC", grass$`Bag Label`), "Closed", 
                                 ifelse(grepl("RO", grass$`Bag Label`), "Open", 
                                        ifelse(grepl("RU", grass$`Bag Label`), "Overhead", 
                                               ifelse(grepl("RE", grass$`Bag Label`), "Elevated", 
                                                      ifelse(grepl("Litter", grass$`Bag Label`), "Environmental Litter", "unknown"))))))
grass$Month <- ifelse(grass$`Time Point` == 1, "End May", 
                      ifelse(grass$`Time Point` == 2, "End June", 
                             ifelse(grass$`Time Point` == 3, "End July",
                                    ifelse(grass$`Time Point` == 4, "Mid September", 
                                           ifelse(grass$`Time Point` == 5, "End October", 
                                                  ifelse(grepl(pattern = "Phyllosphere", x = grass$`Bag Label`), "Before Dispersal", "error... :("))))))

# Removing those treatments that we aren't interested in: 
no <- c( "Phyllosphere", "Environmental Litter")
grass <- grass[!(grass$Treatment %in% no),  ]

grass.not.elevated <- grass[grass$Treatment != "Elevated", ]

# Merge datasets
together.all <- merge(mass, grass, by = "Bag Label")
together.not.elevated <- merge(not.elevated, grass.not.elevated, by = "Bag Label")

# Statistics time:
shapiro.test(together.not.elevated$`Abundance Per G`)
shapiro.test(together.not.elevated$Percent_Mass_Loss)

ggqqplot(together.not.elevated$`Abundance Per G`, ylab = "Abundance")
ggqqplot(together.not.elevated$Percent_Mass_Loss, ylab = "Mass Loss")

cor.test(together.not.elevated$`Abundance Per G`, together.not.elevated$Percent_Mass_Loss,  method = "spearman")
cor.test(together.not.elevated$`Abundance Per G`, together.not.elevated$Percent_Mass_Loss,  method = "pearson")


##### _______ composition of bacteria -- litterbags _______________ ####
q2.table.rarefied <- as.data.frame(fread("OTU_filtered_plus_rarefied__ROUNDED.tsv"))

# Cleaning up the dataframe 
row.names(q2.table.rarefied) <- q2.table.rarefied$V1 # set sample names as row names
q2.table.rarefied$V1 <- NULL # delete column of sample names

q2.table.rarefied$Dispersal_Route <- ifelse(grepl("^RO", row.names(q2.table.rarefied)), "Open",
                                            ifelse(grepl("^RC", row.names(q2.table.rarefied)), "Closed",
                                                   ifelse(grepl("^RU", row.names(q2.table.rarefied)), "Overhead",
                                                          ifelse(grepl("^RE", row.names(q2.table.rarefied)), "Elevated",
                                                                 ifelse(grepl("Litter", row.names(q2.table.rarefied)), "Environmental Litter", "ooops.....")))))

q2.table.rarefied$TimePoint <- str_extract(string = row.names(q2.table.rarefied), pattern = "T[1-5]")

## What relative percentage of the community is the radiation-resistant genus Hymenobacter at T5 in elevated slides?
elevated_T5 <- q2.table.rarefied[q2.table.rarefied$Dispersal_Route == "Elevated" & 
                                   q2.table.rarefied$TimePoint == "T5", ]
elevated_T5$Dispersal_Route <- NULL
elevated_T5$TimePoint <- NULL

elevated_T5 <- as.data.frame(t(elevated_T5))

q2.taxonomy <- as.data.frame(fread("15_taxonomy_filtered.tsv")) # import taxonomy file (3 columns --> OTU ID, taxonomy, and confidence)
elevated_T5.taxa <- merge(elevated_T5, q2.taxonomy[1:2], by.x = "row.names", by.y = "#OTUID") # merge taxonomy and otu table

hymeno <- colSums(elevated_T5.taxa[grepl("Hymenobacter", elevated_T5.taxa$taxonomy), 
                                   c("RET5R1", "RET5R2", "RET5R3", "RET5R4", "RET5R5", "RET5R6", "RET5R7")])

all <- colSums(elevated_T5.taxa[ , c("RET5R1", "RET5R2", "RET5R3", "RET5R4", "RET5R5", "RET5R6", "RET5R7")])

mean(hymeno/all)



##### _______ percentage of community for some genera ___________ #######

# Note: script says family throughout but I changed to look at genus level without updating every family mention

#### Script to make taxa barplots (from figures) ______________
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
q2.fam.agg.tax$family <- gsub("g__", "", q2.fam.agg.tax$family)

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




#### The statistics ___________________

# % Erwinia in the first timepoint for all litterbags
r_and_t <- c("Closed_T1", "Elevated_T1", "Overhead_T1", "Open_T1")

q2.FINAL.together[q2.FINAL.together$Genus == "Erwinia" & 
                    q2.FINAL.together$Route_AND_Time %in% r_and_t & 
                    q2.FINAL.together$TimePoint == "T1", ]$`Relative Abundance` %>% mean()


## % Janthinobacterium at the first and last timepoint for overhead and open treatment
# T1
q2.FINAL.together %>% filter(Genus == "Janthinobacterium", Dispersal_Route %in% c("Open", "Overhead"), TimePoint == "T1") %>% 
  pull("Relative Abundance") %>% mean

# T5
q2.FINAL.together %>% filter(Genus == "Janthinobacterium", Dispersal_Route %in% c("Open", "Overhead"), TimePoint == "T5") %>% 
  pull("Relative Abundance") %>% mean

## % Hymenobacter  at the first and last timepoint for Elevated treatment
# T1
q2.FINAL.together %>% filter(Genus == "Hymenobacter", Dispersal_Route %in% c("Elevated"), TimePoint == "T1") %>% 
  pull("Relative Abundance") %>% mean

# T5
q2.FINAL.together %>% filter(Genus == "Hymenobacter", Dispersal_Route %in% c("Elevated"), TimePoint == "T5") %>% 
  pull("Relative Abundance") %>% mean

## % Paenibacillus & % Sanguibacter in Closed versus other litterbags (open, overhead, and elevated) for whole experiment
# Paenibacillus in Closed
q2.FINAL.together %>% filter(Genus == "Paenibacillus", Dispersal_Route %in% c("Closed"), TimePoint == "Avg") %>% 
  pull("Relative Abundance") %>% sum # 0.08092574

# Paenibacillus in all other
q2.FINAL.together %>% filter(Genus == "Paenibacillus", Dispersal_Route %in% c("Open", "Overhead", "Elevated"), 
                             TimePoint == "Avg") %>% 
  pull("Relative Abundance") %>% sum # 0.009777722

# Sanguibacter in Closed
q2.FINAL.together %>% filter(Genus == "Sanguibacter", Dispersal_Route %in% c("Closed"), TimePoint == "Avg") %>% 
  pull("Relative Abundance") %>% sum # 0.05177081

# Sanguibacter in all other
q2.FINAL.together %>% filter(Genus == "Sanguibacter", Dispersal_Route %in% c("Open", "Overhead", "Elevated"), 
                             TimePoint == "Avg") %>% 
  pull("Relative Abundance") %>% sum # 0.01134932


q2.FINAL.together %>% filter(Dispersal_Route %in% c("Closed"), TimePoint == "Avg") %>% group_by(Genus) %>% 
  summarize(sum(`Relative Abundance`)) %>%  pull(2) %>% sum




##### _______ subtracted mass loss _______ #####
# Input data
mass <- as.data.frame(fread("mass_loss.csv"))

# Set up treatments: 
mass$Treatment <- ifelse(grepl(pattern = "Phyllosphere", x = mass$`Bag Label`), "Phyllosphere", 
                         ifelse(grepl("RC", mass$`Bag Label`), "No Dispersal", 
                                ifelse(grepl("RO", mass$`Bag Label`), "All Dispersal", 
                                       ifelse(grepl("RU", mass$`Bag Label`), "Local Rain/Air", 
                                              ifelse(grepl("RE", mass$`Bag Label`), "Regional Rain/Air", 
                                                     ifelse(grepl("Litter", mass$`Bag Label`), "Environmental Litter", "unknown"))))))

mass$Month <- ifelse(mass$`Time Point` == 1, "End May", 
                     ifelse(mass$`Time Point` == 2, "End June", 
                            ifelse(mass$`Time Point` == 3, "End July",
                                   ifelse(mass$`Time Point` == 4, "Mid September", 
                                          ifelse(mass$`Time Point` == 5, "End October", "error... :(")))))
mass <- mass[mass$Exclude == "No", ]


# Calculate the differences
mass$Replicate <- gsub("R[A-Z]T[1-5]", "", mass$`Bag Label`)
diff.month <- list()
diff <- list()
months <- unique(mass$Month)

mass.outdf <- data.frame("Percent_Mass_Loss" = numeric(), "Treatment" = character(), 
                         "Replicate" = character(), "Month" = character())

for (i in 1:length(months)) {
  key <- as.character(months[[i]])
  sub <- mass[mass$Month == key, ]
  rep <- unique(sub$Replicate)
  
  for (j in 1:length(rep)) {
    reprep <- rep[[j]]
    
    subsub <- sub[sub$Replicate == reprep, ]
    air.vegetation <- subsub[subsub$Treatment == "Local Rain/Air", ]$Percent_Mass_Loss - subsub[subsub$Treatment == "No Dispersal", ]$Percent_Mass_Loss
    soil <- subsub[subsub$Treatment == "All Dispersal", ]$Percent_Mass_Loss - subsub[subsub$Treatment == "Local Rain/Air", ]$Percent_Mass_Loss
    
    out <- data.frame(c("Soil" = soil, "Air + Vegetation" = air.vegetation))
    out$Treatment <- row.names(out)
    out$Replicate <- c(rep(reprep, nrow(out)))
    out$Month <- c(rep(key, nrow(out)))
    names(out) <- c("Percent_Mass_Loss", "Treatment", "Replicate", "Month")
    
    mass.outdf <- rbind(mass.outdf, out)
    
  }
}

## Statistics: 

# Overall model
options(contrasts = c("contr.sum","contr.poly"))

mod <- lm(Percent_Mass_Loss ~ Treatment * Month, data = mass.outdf)
Anova(mod, type = "III")

tidy_aov <- tidy(Anova(mod, type = "III"))

r_squared_interaction <- tidy_aov[tidy_aov$term == "Treatment:Month", ]$sumsq
r_squared_residuals <- tidy_aov[tidy_aov$term == "Residuals", ]$sumsq

r_squared_interaction / (r_squared_interaction + r_squared_residuals) # 0.2490881


t1 <- mass.outdf[mass.outdf$Month == "End May", ]
t2 <- mass.outdf[mass.outdf$Month == "End June", ]
t3 <- mass.outdf[mass.outdf$Month == "End July", ]
t4 <- mass.outdf[mass.outdf$Month == "Mid September", ]
t5 <- mass.outdf[mass.outdf$Month == "End October", ]


# Testing differences from zero on the graph

t.test(t1[t1$Treatment == "Air + Vegetation", ]$Percent_Mass_Loss, mu = 0) # the only one that is significant!
t.test(t1[t1$Treatment == "Soil", ]$Percent_Mass_Loss, mu = 0)

t.test(t2[t2$Treatment == "Air + Vegetation", ]$Percent_Mass_Loss, mu = 0)
t.test(t2[t2$Treatment == "Soil", ]$Percent_Mass_Loss, mu = 0)

t.test(t3[t3$Treatment == "Air + Vegetation", ]$Percent_Mass_Loss, mu = 0)
t.test(t3[t3$Treatment == "Soil", ]$Percent_Mass_Loss, mu = 0)

t.test(t4[t4$Treatment == "Air + Vegetation", ]$Percent_Mass_Loss, mu = 0)
t.test(t4[t4$Treatment == "Soil", ]$Percent_Mass_Loss, mu = 0)

t.test(t5[t5$Treatment == "Air + Vegetation", ]$Percent_Mass_Loss, mu = 0)
t.test(t5[t5$Treatment == "Soil", ]$Percent_Mass_Loss, mu = 0)


# T1 ANOVA
mod <- lm(Percent_Mass_Loss ~ Treatment, data = t1)
Anova(mod, type = "III")
TukeyHSD(aov(mod))

tidy_aov <- tidy(Anova(mod, type = "III"))

r_squared_interaction <- tidy_aov[tidy_aov$term == "Treatment", ]$sumsq
r_squared_residuals <- tidy_aov[tidy_aov$term == "Residuals", ]$sumsq

r_squared_interaction / (r_squared_interaction + r_squared_residuals) # 0.5672412


# T2 ANOVA
mod <- lm(Percent_Mass_Loss ~ Treatment, data = t2)
Anova(mod, type = "III")


# T3 ANOVA
mod <- lm(Percent_Mass_Loss ~ Treatment, data = t3)
Anova(mod, type = "III")


# T4 ANOVA
mod <- lm(Percent_Mass_Loss ~ Treatment, data = t4)
Anova(mod, type = "III")

# T5 ANOVA
mod <- lm(Percent_Mass_Loss ~ Treatment, data = t5)
Anova(mod, type = "III")



##### _______ death rate ______________________________________ #####
### ____ PANEL A ______________ 
# Input data
death <- as.data.frame(fread("death_slides_set1_abundance_all_timepoints.csv"))


death$Month <- ifelse(death$`Time Point` == 1, "End May", 
                      ifelse(death$`Time Point` == 2, "End June", 
                             ifelse(death$`Time Point` == 3, "End July",
                                    ifelse(death$`Time Point` == 4, "Mid September", 
                                           ifelse(death$`Time Point` == 5, "End October", 
                                                  ifelse(death$`Time Point` == 0, "Start", "error... :("))))))

death$Month <- factor(death$Month, levels = c("Start", "End May", "End June", "End July", "Mid September", "End October"), ordered = TRUE)

death$`Time Point` <- as.character(death$`Time Point`)
death$`Time Point` <- ifelse(death$`Time Point` == "4", "4.5", 
                             ifelse(death$`Time Point` == "5", "6", death$`Time Point`))
death$`Time Point` <- as.numeric(death$`Time Point`)

death$Time_as_days <- ifelse(grepl(pattern = "0", death$`Time Point`), 0, 
                             ifelse(grepl(pattern = "1", death$`Time Point`), 39, 
                                    ifelse(grepl(pattern = "2", death$`Time Point`), 60, 
                                           ifelse(grepl(pattern = "3", death$`Time Point`), 100, 
                                                  ifelse(grepl(pattern = "4.5", death$`Time Point`), 151, 
                                                         ifelse(grepl(pattern = "6", death$`Time Point`), 195, "EIIIIIEEEEEIIII!!"))))))
death$Time_as_days <- as.numeric(death$Time_as_days)

##### Calculating death rate for DAILY NUMBERS: 
intercept <- mean(log10(death[death$Time_as_days == 0, ]$Abundance))

lm_death <- lm(I(log10(Abundance) - intercept) ~ Time_as_days + 0, death) # estimate: -0.014771
summary(lm_death)
1 - 10^-0.014771 # 0.0334396 === 3.34%


# can do weekly
dead.weekly <- 7*-0.014771
round(100*(1 - 10^dead.weekly), 2) # 21.19% overall

# can do monthly
dead.30 <- 30*-0.014771
round(100*(1 - 10^dead.30), 2) # 63.95% overall

dead.31 <- 31*-0.014771
round(100*(1 - 10^dead.31), 2) # 65.16% overall 

### ____ PANEL B ______________ 
new_death <- as.data.frame(fread("death_slides_set2_abundance_all_timepoints.csv"))

# Entering the days since start of experiment, not just the timepoint
days.key <- c("0" = 0, "1" = 7, "2" = 14, "3" = 21, "4" = 28, "5" = 44)
new_death$`Time Point` <- as.character(new_death$`Time Point`)
new_death$Time_as_days <- days.key[new_death$`Time Point`]

# Entering treatment (sorry about all my nested ifelse statements... )
new_death$Treatment <- ifelse(grepl(pattern = "NDU", x = new_death$`Bag Label`), "Table", 
                              ifelse(grepl(pattern = "morning", x = new_death$`Bag Label`), "Initial Abundance", 
                                     ifelse(grepl(pattern = "night", x = new_death$`Bag Label`), "Night", "Ground")))

# Making sure the T0 abundance gets applied to both the table and ground samples, since it's the starting place for both
new_table <- new_death[new_death$Treatment == "Initial Abundance", ]
new_ground <- new_death[new_death$Treatment == "Initial Abundance", ]
new_table$Treatment <- c(rep("Table", nrow(new_table)))
new_ground$Treatment <- c(rep("Ground", nrow(new_ground)))

# Putting them all together (we are just duplicating the initial abundance measurements, half for table samples and half for ground samples)
small <- new_death[new_death$Treatment == "Ground" | new_death$Treatment == "Table" , ]
bi <- rbind(small, new_table)
big <- rbind(bi, new_ground)

# STATS: do table and ground differ in their rates:
big$Treatment <- as.factor(big$Treatment)
big$Time_as_days <- as.numeric(big$Time_as_days)

library(car)
intercept.table <- mean(log10(big[big$Time_as_days == 0, ]$Abundance)) # intercept.table and intercept.ground are the same
Anova(lm(formula = I(log10(Abundance)-intercept.table) ~ Time_as_days:Treatment + Time_as_days + 0, data = big), type = "II") # independent of order

# Table
table <- big[big$Treatment == "Table", ]
intercept <- mean(log10(table[table$Time_as_days == 0, ]$Abundance))

lm_death <- lm(I(log10(Abundance) - intercept) ~ Time_as_days + 0, table) # estimate: -0.014771
summary(lm_death)

# Ground
ground <- big[big$Treatment == "Ground", ]
intercept <- mean(log10(ground[ground$Time_as_days == 0, ]$Abundance))

lm_death <- lm(I(log10(Abundance) - intercept) ~ Time_as_days + 0, ground) # estimate: -0.014771
summary(lm_death)

##### _______ SourceTracker _______________ #####

# Different superscripts letters in the tables indicate significant difference of pairwise comparisons 
# to date using Dunn’s post hoc test with Bonferroni correction.

# Loading output from using Source Tracker in the terminal
proportions <- as.data.frame(fread("mixing_proportions.txt"))


# doing the thing
proportions$Dispersal_Route <- ifelse(grepl("^LO", proportions$SampleID), "Open", 
                                      ifelse(grepl("^LC", proportions$SampleID), "Closed", 
                                             ifelse(grepl("^LU", proportions$SampleID), "Up", 
                                                    ifelse(grepl("^LE", proportions$SampleID), "Elevated",
                                                           ifelse(grepl("Litter", proportions$SampleID), "Environmental Litter", 
                                                                  ifelse(grepl("LD", proportions$SampleID), "Death",
                                                                         ifelse(grepl("Soil", proportions$SampleID), "Soil",
                                                                                ifelse(grepl("Air", proportions$SampleID), "Air", "ooops....."))))))))

proportions.long <- gather(proportions, Source, Proportion, Air:Unknown, factor_key=TRUE)

# Overall model -- 
model <- lm(Proportion ~ Dispersal_Route * Source, data = proportions.long)
Anova(model, type = "III")


# Get medians of percentages
up_open <- proportions[proportions$Dispersal_Route == "Up" | proportions$Dispersal_Route == "Open", ]
median(up_open$Litter) # 0.40
mean(up_open$Litter) # 0.39

median(up_open$Soil)
median(up_open$Air)


# Dunn test -- get multiple comparisons at each dispersal route 
elevated <- proportions.long[proportions.long$Dispersal_Route == "Elevated", ]
overhead <- proportions.long[proportions.long$Dispersal_Route == "Up", ]
open <- proportions.long[proportions.long$Dispersal_Route == "Open", ]

dunnTest(Proportion ~ Source, data = elevated, method = "bonferroni")
dunnTest(Proportion ~ Source, data = overhead, method = "bonferroni")
dunnTest(Proportion ~ Source, data = open, method = "bonferroni")

model <- lm(Proportion ~ Source, data = elevated)
Anova(model, type = "III")
kruskal.test(Proportion ~ Source, data = elevated)

model <- lm(Proportion ~ Source, data = overhead)
Anova(model, type = "III")
kruskal.test(Proportion ~ Source, data = overhead)

model <- lm(Proportion ~ Source, data = open)
Anova(model, type = "III")
kruskal.test(Proportion ~ Source, data = open)







##### _______ alpha- and beta-diversity for glass + grass ____#######
### PANEL A ==== Glass Shannon 
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


# Statistics
just.bags <- merged_alpha[merged_alpha$Dispersal_Route_new_names %in% c("Elevated", "Overhead", "Open"), ]

options(contrasts = c("contr.sum","contr.poly"))
model <- lm(Shannon ~ Dispersal_Route_new_names, data = just.bags)
my.anova <- Anova(model, type = "III")
my.anova
TukeyHSD(aov(model))

# Calculating R2
tidy_aov <- tidy(my.anova) 
sum_squares_treatment <- tidy_aov[tidy_aov$term == "Dispersal_Route_new_names", ]$sumsq
sum_squares_residuals <- tidy_aov[tidy_aov$term == "Residuals", ]$sumsq

R_squared <- sum_squares_treatment /
  (sum_squares_treatment + sum_squares_residuals)

R_squared


### PANEL B ==== Litterbag Bacteria Shannon 
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

# Statistics
just.bags <- merged_alpha2[merged_alpha2$Dispersal_Route %in% c("Closed", "Elevated", "Overhead", "Open"), ]

options(contrasts = c("contr.sum","contr.poly"))
model <- lm(Shannon ~ Dispersal_Route, data = just.bags)
my.anova <- Anova(model, type = "III")
my.anova
TukeyHSD(aov(model))

# Calculating R2
tidy_aov <- tidy(my.anova) 
sum_squares_treatment <- tidy_aov[tidy_aov$term == "Dispersal_Route", ]$sumsq
sum_squares_residuals <- tidy_aov[tidy_aov$term == "Residuals", ]$sumsq

R_squared <- sum_squares_treatment /
  (sum_squares_treatment + sum_squares_residuals)

R_squared



### PANEL C ==== Litterbag Fungi Shannon 
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

# LOOKING AT DISPERSAL ROUTE:
routes_keep <- c("Environmental Litter", "Closed", "Elevated", "Up", "Open")
merged_alpha2 <- merged_alpha2[merged_alpha2$Dispersal_Route %in% routes_keep, ]
merged_alpha2$Dispersal_Route <- str_replace(merged_alpha2$Dispersal_Route, "Up", "Overhead")
merged_alpha2$Dispersal_Route <- factor(merged_alpha2$Dispersal_Route, levels = c("Phyllosphere", "Environmental Litter", "Closed", "Elevated", "Overhead", "Open"))


# Statistics
just.bags <- merged_alpha2[merged_alpha2$Dispersal_Route %in% c("Closed", "Elevated", "Overhead", "Open"), ]

options(contrasts = c("contr.sum","contr.poly"))
model <- lm(Shannon ~ Dispersal_Route, data = just.bags)
my.anova <- Anova(model, type = "III")
my.anova
TukeyHSD(aov(model))

ggplot(just.bags) + geom_boxplot(aes(x = Dispersal_Route, y = Shannon))

# Calculating R2
tidy_aov <- tidy(my.anova) 
sum_squares_treatment <- tidy_aov[tidy_aov$term == "Dispersal_Route", ]$sumsq
sum_squares_residuals <- tidy_aov[tidy_aov$term == "Residuals", ]$sumsq

R_squared <- sum_squares_treatment /
  (sum_squares_treatment + sum_squares_residuals)

R_squared

### PANEL D ==== Glass Slide Bacteria Beta-Diversity 
## Load data
bray.dist.df <- as.data.frame(fread("median_bray_curtis_sqrt_q1_single_rarefied_to_1000.txt"))
bray.dist.df <- bray.dist.df[!(grepl("LD", bray.dist.df$V1)), !(grepl("LD", names(bray.dist.df)))]
row.names(bray.dist.df) <- bray.dist.df$V1
bray.dist.df$V1 <- NULL

# Just bags
just.bags.df <- bray.dist.df[grepl("^L[UOE]", row.names(bray.dist.df)), 
                             grepl("^L[UOE]", names(bray.dist.df))]

bray.dist <- as.dist(as.matrix(just.bags.df)) # turning into a distance object

## make our groups
groups <- ifelse(grepl("Air", names(just.bags.df)), "Air", 
                 ifelse(grepl("Litter", names(just.bags.df)), "Environmental Litter", 
                        ifelse(grepl("Soil", names(just.bags.df)), "Soil", 
                               ifelse(grepl("LE", names(just.bags.df)), "Elevated", 
                                      ifelse(grepl("LU", names(just.bags.df)), "Overhead", 
                                             ifelse(grepl("LO", names(just.bags.df)), "Open",
                                                    ifelse(grepl("LC", names(just.bags.df)), "Closed", "there was a problem..")))))))

## do the test! 
disper <- betadisper(bray.dist, group = groups, type = c("centroid"))
permutest(disper, pairwise = TRUE)
0.044657 / (0.044657 + 0.235410) # R2



### PANEL E ==== Litterbag Bacteria Beta-Diversity 
## Load data
bray.dist.df <- as.data.frame(fread("median_bray_curtis_sqrt_q2_single_rarefied_to_1000.txt"))
row.names(bray.dist.df) <- bray.dist.df$V1

bray.dist <- as.dist(as.matrix(bray.dist.df))

# Just bags
just.bags.df <- bray.dist.df[grepl("^R[CUOE]", row.names(bray.dist.df)), 
                             grepl("^R[CUOE]", names(bray.dist.df))]

bray.dist <- as.dist(as.matrix(just.bags.df)) # turning into a distance object

## make our groups
groups <- ifelse(grepl("Air", names(just.bags.df)), "Air", 
                 ifelse(grepl("Litter", names(just.bags.df)), "Environmental Litter", 
                        ifelse(grepl("Soil", names(just.bags.df)), "Soil", 
                               ifelse(grepl("RE", names(just.bags.df)), "Elevated", 
                                      ifelse(grepl("RU", names(just.bags.df)), "Overhead", 
                                             ifelse(grepl("RO", names(just.bags.df)), "Open",
                                                    ifelse(grepl("RC", names(just.bags.df)), "Closed", "there was a problem..")))))))

## do the test! 
disper <- betadisper(bray.dist, group = groups, type = c("centroid"))
permutest(disper, pairwise = TRUE)

0.20534 / (0.20534 + 0.77329) # R2




### PANEL F ==== Litterbag Fungi Beta-Diversity 

## Load data
bray.dist.df <- as.data.frame(fread("median_bray_curtis_sqrt_q2_ITS_single_rarefied_to_3500.txt"))
row.names(bray.dist.df) <- bray.dist.df$V1 
bray.dist.df$V1 <- NULL


# Just bags
just.bags.df <- bray.dist.df[grepl("^R[CUOE]", row.names(bray.dist.df)), 
                             grepl("^R[CUOE]", names(bray.dist.df))]

bray.dist <- as.dist(as.matrix(just.bags.df)) # turning into a distance object

## make our groups
groups <- ifelse(grepl("Litter", names(just.bags.df)), "Environmental Litter", 
                 ifelse(grepl("RE", names(just.bags.df)), "Elevated", 
                        ifelse(grepl("RU", names(just.bags.df)), "Overhead", 
                               ifelse(grepl("RO", names(just.bags.df)), "Open",
                                      ifelse(grepl("RC", names(just.bags.df)), "Closed", "there was a problem..")))))

## do the test! 
disper.fun <- betadisper(bray.dist, group = groups, type = c("centroid"))
permutest(disper.fun, permutations = 50,  pairwise = TRUE)

0.25457 / (0.25457 + 0.43325) # R2


##### _______ raw abundance in litterbags _______ #########
grass <- as.data.frame(fread("grass_litterbag_abundance_all_timepoints.csv"))


## Make treatment names:
grass$Treatment <- ifelse(grepl(pattern = "Phyllosphere", x = grass$`Bag Label`), "Phyllosphere", 
                          ifelse(grepl("RC", grass$`Bag Label`), "No Dispersal", 
                                 ifelse(grepl("RO", grass$`Bag Label`), "All Dispersal", 
                                        ifelse(grepl("RU", grass$`Bag Label`), "Local Rain/Air", 
                                               ifelse(grepl("RE", grass$`Bag Label`), "Regional Rain/Air", 
                                                      ifelse(grepl("Litter", grass$`Bag Label`), "Environmental Litter", "unknown"))))))
grass$Month <- ifelse(grass$`Time Point` == 1, "End May", 
                      ifelse(grass$`Time Point` == 2, "End June", 
                             ifelse(grass$`Time Point` == 3, "End July",
                                    ifelse(grass$`Time Point` == 4, "Mid September", 
                                           ifelse(grass$`Time Point` == 5, "End October", 
                                                  ifelse(grepl(pattern = "Phyllosphere", x = grass$`Bag Label`), "Before Dispersal", "error... :("))))))

just_bags <- grass[!(grass$Treatment %in% c("Environmental Litter", "Phyllosphere")), ]

just_bags$Replicate <- gsub("R[A-Z]T[1-5]", "", just_bags$`Bag Label`)

## Statistics
# Overall model + r2 calculation for the entire figure
options(contrasts = c("contr.sum","contr.poly"))

mod <- lm(`Abundance Per G` ~ Treatment + Treatment:Month, just_bags)
Anova(mod, type = "III")

model <- Anova(mod, type = "III")
tidy_aov <- tidy(model) 
sum_squares_treatment <- tidy_aov[tidy_aov$term == "Treatment", ]$sumsq # put the term here you are interested in
sum_squares_treatment_by_month <- tidy_aov[tidy_aov$term == "Treatment:Month", ]$sumsq # put the term here you are interested in
sum_squares_residuals <- tidy_aov[tidy_aov$term == "Residuals", ]$sumsq

R_squared_treatment <- sum_squares_treatment /
  (sum_squares_treatment + sum_squares_treatment_by_month + sum_squares_residuals)

R_squared_treatment_by_month <- sum_squares_treatment_by_month /
  (sum_squares_treatment + sum_squares_treatment_by_month + sum_squares_residuals)

R_squared_treatment
R_squared_treatment_by_month

# by timepoint
t1 <- just_bags[just_bags$Month == "End May", ]
t1$Treatment <- as.factor(t1$Treatment)
mod <- lm(`Abundance Per G` ~ Treatment, t1)
Anova(mod, type = "III") # P = 0.0498053
TukeyHSD(aov(`Abundance Per G` ~ Treatment, t1))
summary(glht(mod, linfct = mcp(Treatment = "Tukey")))

t2 <- just_bags[just_bags$Month == "End June", ]
mod <- lm(`Abundance Per G` ~ Treatment, t2)
Anova(mod, type = "III") # P = 0.001289
TukeyHSD(aov(`Abundance Per G` ~ Treatment, t2))

t3 <- just_bags[just_bags$Month == "End July", ]
mod <- lm(`Abundance Per G` ~ Treatment, t3)
Anova(mod, type = "III") # P = 0.02518
TukeyHSD(aov(`Abundance Per G` ~ Treatment, t3))

t4 <- just_bags[just_bags$Month == "Mid September", ]
mod <- lm(`Abundance Per G` ~ Treatment, t4)
Anova(mod, type = "III") # P = 0.117349

t5 <- just_bags[just_bags$Month == "End October", ]
mod <- lm(`Abundance Per G` ~ Treatment, t5)
Anova(mod, type = "III") # P = 0.3556
TukeyHSD(aov(`Abundance Per G` ~ Treatment, t5))

# did the communities receiving air dispersal have HIGHER abundance than the control? 
var.test(`Abundance Per G` ~ Treatment, data = just_bags[just_bags$Treatment == "No Dispersal" | just_bags$Treatment == "Regional Rain/Air", ]) # P < 2.2e-16 so they DON'T have the same variance
t.test(just_bags[just_bags$Treatment == "No Dispersal", ]$`Abundance Per G`, just_bags[just_bags$Treatment == "Regional Rain/Air", ]$`Abundance Per G`, 
       var.equal = FALSE, alternative = "less")

# did the communities receiving soil dispersal ever DIFFER from communities receiving vegetation + air? 
var.test(`Abundance Per G` ~ Treatment, data = just_bags[just_bags$Treatment == "All Dispersal" | just_bags$Treatment == "Local Rain/Air", ]) # P = 0.9706 so they DO have the same variance
t.test(just_bags[just_bags$Treatment == "All Dispersal", ]$`Abundance Per G`, just_bags[just_bags$Treatment == "Local Rain/Air", ]$`Abundance Per G`, 
       var.equal = TRUE)

# did the communities receiving vegetation dispersal have HIGHER abundance than those just receiving air? 
var.test(`Abundance Per G` ~ Treatment, data = just_bags[just_bags$Treatment == "Local Rain/Air" | just_bags$Treatment == "Regional Rain/Air", ]) # p-value < 2.2e-16 so they DON'T have the same variance
t.test(just_bags[just_bags$Treatment == "Local Rain/Air", ]$`Abundance Per G`, just_bags[just_bags$Treatment == "Regional Rain/Air", ]$`Abundance Per G`, 
       var.equal = FALSE, alternative = "greater")


# Subset weight data before treatment
air <- just_bags[just_bags$Treatment == "Regional Rain/Air", ]
air <- air[order(air$`Time Point`, air$Replicate), ]
air <- air[air$`Bag Label` != "RET5R3", ] # to match vegetation which doesn't have this point

# subset weight data after treatment
vegetation <- just_bags[just_bags$Treatment == "Local Rain/Air", ]
vegetation <- vegetation[order(vegetation$`Time Point`, vegetation$Replicate), ]



##### _______ correlation between litter abundance and raw mass loss ________ #####
# Input data
mass <- as.data.frame(fread("mass_loss.csv"))

# Set up treatments: 
mass$Treatment <- ifelse(grepl(pattern = "Phyllosphere", x = mass$`Bag Label`), "Phyllosphere", 
                         ifelse(grepl("RC", mass$`Bag Label`), "No Dispersal", 
                                ifelse(grepl("RO", mass$`Bag Label`), "All Dispersal", 
                                       ifelse(grepl("RU", mass$`Bag Label`), "Local Rain/Air", 
                                              ifelse(grepl("RE", mass$`Bag Label`), "Regional Rain/Air", 
                                                     ifelse(grepl("Litter", mass$`Bag Label`), "Environmental Litter", "unknown"))))))

mass$Month <- ifelse(mass$`Time Point` == 1, "End May", 
                     ifelse(mass$`Time Point` == 2, "End June", 
                            ifelse(mass$`Time Point` == 3, "End July",
                                   ifelse(mass$`Time Point` == 4, "Mid September", 
                                          ifelse(mass$`Time Point` == 5, "End October", "error... :(")))))
mass <- mass[mass$Exclude == "No", ]

not.elevated <- mass[mass$Treatment != "Regional Rain/Air", ]



## Raw abundance on litterbags
grass <- as.data.frame(fread("grass_litterbags_abundance_all_timepoints.csv"))


grass$Treatment <- ifelse(grepl(pattern = "Phyllosphere", x = grass$`Bag Label`), "Phyllosphere", 
                          ifelse(grepl("RC", grass$`Bag Label`), "Closed", 
                                 ifelse(grepl("RO", grass$`Bag Label`), "Open", 
                                        ifelse(grepl("RU", grass$`Bag Label`), "Overhead", 
                                               ifelse(grepl("RE", grass$`Bag Label`), "Elevated", 
                                                      ifelse(grepl("Litter", grass$`Bag Label`), "Environmental Litter", "unknown"))))))
grass$Month <- ifelse(grass$`Time Point` == 1, "End May", 
                      ifelse(grass$`Time Point` == 2, "End June", 
                             ifelse(grass$`Time Point` == 3, "End July",
                                    ifelse(grass$`Time Point` == 4, "Mid September", 
                                           ifelse(grass$`Time Point` == 5, "End October", 
                                                  ifelse(grepl(pattern = "Phyllosphere", x = grass$`Bag Label`), "Before Dispersal", "error... :("))))))
# Removing those treatments that we aren't interested in: 
no <- c( "Phyllosphere", "Environmental Litter")
grass <- grass[!(grass$Treatment %in% no),  ]

grass.not.elevated <- grass[grass$Treatment != "Elevated", ]

# Merge
together.all <- merge(mass, grass, by = "Bag Label")
together.not.elevated <- merge(not.elevated, grass.not.elevated, by = "Bag Label")

## Correlation time
summary(lm(`Mass Loss (g)` ~ log(`Abundance Per G`), data = together.all))
summary(lm(`Mass Loss (g)` ~ log(`Abundance Per G`), data = together.not.elevated))




##### _______ similarity of litterbags to environmental litter _________ #######

# Load data + remove unwanted samples
to.centroid <- as.data.frame(fread("distances_each_sample_to_env_by_timepoint.txt"))
distances <- to.centroid[!(grepl(pattern = "Rep [1-5]|Rep1", to.centroid$Sample)), ] #take out replicates

# Set up metadata for visualization
months.key <- c("T1" = "May", "T2" = "June", "T3" = "July", "T4" = "September", "T5" = "October")
distances$Month <- months.key[str_extract(distances$Sample, "T[1-5]")]

treatment.key <- c("RC" = "Closed", "RE" = "Elevated", "RU" = "Overhead", "RO" = "Open")
distances$Treatment <- treatment.key[str_extract(distances$Sample, "R[EOUC]")]

# Stats
options(contrasts = c("contr.sum","contr.poly"))

mod <- lm(`Abundance Per G` ~ Treatment + Treatment:Month, just_bags)
Anova(mod, type = "III")

model <- Anova(mod, type = "III")
tidy_aov <- tidy(model) 
sum_squares_treatment <- tidy_aov[tidy_aov$term == "Treatment", ]$sumsq # put the term here you are interested in
sum_squares_treatment_by_month <- tidy_aov[tidy_aov$term == "Treatment:Month", ]$sumsq # put the term here you are interested in
sum_squares_residuals <- tidy_aov[tidy_aov$term == "Residuals", ]$sumsq

R_squared_treatment <- sum_squares_treatment /
  (sum_squares_treatment + sum_squares_treatment_by_month + sum_squares_residuals)

R_squared_treatment_by_month <- sum_squares_treatment_by_month /
  (sum_squares_treatment + sum_squares_treatment_by_month + sum_squares_residuals)

R_squared_treatment
R_squared_treatment_by_month

# by timepoint
distances$Treatment <- as.factor(distances$Treatment)

t1 <- distances[distances$Month == "May", ]
mod <- lm(Distance_to_Env_Centroid ~ Treatment, t1)
Anova(mod, type = "III") # not sig

t2 <- distances[distances$Month == "June", ]
mod <- lm(Distance_to_Env_Centroid ~ Treatment, t2)
Anova(mod, type = "III") # P = 1.396e-06 
TukeyHSD(aov(Distance_to_Env_Centroid ~ Treatment, t2))
summary(glht(mod, linfct = mcp(Treatment = "Tukey"))) #another way to do this

t3 <- distances[distances$Month == "July", ]
mod <- lm(Distance_to_Env_Centroid ~ Treatment, t3)
Anova(mod, type = "III") # P = 0.001995
TukeyHSD(aov(Distance_to_Env_Centroid ~ Treatment, t3))

t4 <- distances[distances$Month == "September", ]
mod <- lm(Distance_to_Env_Centroid ~ Treatment, t4)
Anova(mod, type = "III") # P = 0.04086
TukeyHSD(aov(Distance_to_Env_Centroid ~ Treatment, t4))

t5 <- distances[distances$Month == "October", ]
mod <- lm(Distance_to_Env_Centroid ~ Treatment, t5)
Anova(mod, type = "III") # P = 0.3556
TukeyHSD(aov(Distance_to_Env_Centroid ~ Treatment, t5))


##### _______ Figure S7: raw data mass loss ________ ######
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

not.elevated <- mass[mass$Treatment != "Elevated", ]

# What are the rates?
not.elevated %>% group_by(c(`Time Point`)) %>% summarize(mean = mean(`Mass Loss (g)`))

# Stats: 
t1 <- not.elevated[not.elevated$Month ==  "May", ]
mod <- lm(`Mass Loss (g)` ~ Treatment, data = t1)
summary(aov(mod))
0.05456 / (0.05456 + 0.04902) # R2
TukeyHSD(aov(mod))

t2 <- not.elevated[not.elevated$Month == "June", ]
mod <- lm(`Mass Loss (g)` ~ Treatment, data = t2)
Anova(mod)

t3 <- not.elevated[not.elevated$Month == "July", ]
mod <- lm(`Mass Loss (g)` ~ Treatment, data = t3)
Anova(mod)

t4 <- not.elevated[not.elevated$Month == "September", ]
mod <- lm(`Mass Loss (g)` ~ Treatment, data = t4)
Anova(mod)

t5 <- not.elevated[not.elevated$Month == "October", ]
mod <- lm(`Mass Loss (g)` ~ Treatment, data = t5)
Anova(mod)





##### _______ light intensity -- elevated v ground _______ ########
# Room of ()ments
require(data.table)
require(ggplot2)
require(lubridate)

# Process the data
files <- list.files(pattern = "*.csv")
data <- list()

for (i in 1:length(files)) {
  
  
  tmp <- as.data.frame(fread(files[[i]]))
  
  
  tmp <- tmp[ , 2:4]
  colnames(tmp) <- c("Date_Time", "Temp_C", "Light_Intensity_Lux")
  tmp$Date_Time <- mdy_hms(tmp$Date_Time)
  tmp$Day <- factor(strftime(x = tmp$Date_Time, format = "%m/%d", tz = "UTC"))
  
  tmp <- tmp[tmp$Light_Intensity_Lux != 0, ]
  tmp$Location <- c(rep(files[[i]], nrow(tmp)))
  
  
  data[[i]] <- tmp
}


# Rbind the dataframes and clean up the data a little bit
mega <- do.call("rbind", data)

mega$Location <- ifelse(grepl(pattern = "ground", x = mega$Location), "Ground", 
                        ifelse(grepl(pattern = "table", x = mega$Location), "Table", "whoops"))
days_no <- c("03/06", "03/07")
mega <- mega[!(mega$Day %in% days_no) & !(is.na(mega$Day)) , ]

# Subset to just daylight hours
day <- with(mega , mega[hour( Date_Time ) >= 7.10 & hour( Date_Time ) < 18 , ] ) # note: these are average civil twilight hours

ground.light <- day[day$Location == "Ground", ]$Light_Intensity_Lux
table.light <- day[day$Location == "Table", ]$Light_Intensity_Lux

# Test for assumptions for t-test
var.test(Light_Intensity_Lux ~ Location, data = day) # p-value = 0.113 means equal variance

qqnorm(ground.light)
qqline(ground.light)
hist(ground.light) # we have some serious skew so Welch t-test here we come

qqnorm(table.light)
qqline(table.light)
hist(table.light) # better but still skewed 

t.test(table.light, ground.light, alternative = c("greater"), var.equal = FALSE)



##### _______ fungi (not phyllosphere) ___________________ #######
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

q2.fam.agg.tax.agg <- aggregate(. ~ family, data = q2.fam.agg.tax[c(1:(length(names(q2.fam.agg.tax))-2), match("family", names(q2.fam.agg.tax)))], FUN = sum) # aggregate by fammily
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

#### The statistics ___________________

# % Aureobasidium on average for all Elevated, Overhead, and Open litterbags
route <- c("Elevated", "Overhead", "Open")

q2.FINAL.together[q2.FINAL.together$Genus == "Aureobasidium" & 
                    q2.FINAL.together$Dispersal_Route %in% route & 
                    q2.FINAL.together$TimePoint == "Avg", ]$`Relative Abundance` %>% sum / length(route)

# % Aureobasidium on average for Closed litterbags
route <- c("Closed")

q2.FINAL.together[q2.FINAL.together$Genus == "Aureobasidium" & 
                    q2.FINAL.together$Dispersal_Route %in% route & 
                    q2.FINAL.together$TimePoint == "Avg", ]$`Relative Abundance` %>% sum / length(route)


# % Filobasidium  on average for all Elevated, Overhead, and Open litterbags
route <- c("Elevated", "Overhead", "Open")

q2.FINAL.together[q2.FINAL.together$Genus == "Filobasidium" & 
                    q2.FINAL.together$Dispersal_Route %in% route & 
                    q2.FINAL.together$TimePoint == "Avg", ]$`Relative Abundance` %>% sum / length(route)


# % Filobasidium  on average for Closed litterbags
route <- c("Closed")

q2.FINAL.together[q2.FINAL.together$Genus == "Filobasidium" & 
                    q2.FINAL.together$Dispersal_Route %in% route & 
                    q2.FINAL.together$TimePoint == "Avg", ]$`Relative Abundance` %>% sum / length(route)


# % Paraconiothyrium on average for all Overhead and Open litterbags
route <- c("Overhead", "Open")

q2.FINAL.together[q2.FINAL.together$Genus == "Paraconiothyrium" & 
                    q2.FINAL.together$Dispersal_Route %in% route & 
                    q2.FINAL.together$TimePoint == "Avg", ]$`Relative Abundance` %>% sum / length(route)


# % Paraconiothyrium on average for all Closed and Elevated litterbags
route <- c("Elevated","Closed")

q2.FINAL.together[q2.FINAL.together$Genus == "Paraconiothyrium" & 
                    q2.FINAL.together$Dispersal_Route %in% route & 
                    q2.FINAL.together$TimePoint == "Avg", ]$`Relative Abundance` %>% sum / length(route)







##### _______ correlation between abundance on glass slides and litterbags ________ #######
glass <- as.data.frame(fread("glass_slides_abundance_all_timepoints.csv"))

## Make treatment names:
glass$Treatment <- ifelse(grepl("LC", glass$`Bag Label`), "Closed", 
                          ifelse(grepl("LO", glass$`Bag Label`), "Open", 
                                 ifelse(grepl("LU", glass$`Bag Label`), "Overhead", 
                                        ifelse(grepl("LE", glass$`Bag Label`), "Elevated", "unknown"))))

glass$Month <- ifelse(glass$`Time Point` == 1, "End May", 
                      ifelse(glass$`Time Point` == 2, "End June", 
                             ifelse(glass$`Time Point` == 3, "End July",
                                    ifelse(glass$`Time Point` == 4, "Mid September", 
                                           ifelse(glass$`Time Point` == 5, "End October", "error... :(")))))

glass$Abundance_per_cm2 <- glass$`Abundance`/(2.5*7.5)
glass$Month <- factor(glass$Month, levels = c("End May", "End June", "End July", "Mid September", "End October"), ordered = TRUE)
glass$Replicate <- gsub("L[A-Z]T[1-5]", "", glass$`Bag Label`)


## Remove LUT4R1 (252), LOT4R1 (257), and LOT5R4 (326), two of which were literally COVERED in soil
remove <- c("LUT4R1", "LOT4R1", "LOT5R4")
glass_red <- glass[!(glass$`Bag Label` %in% remove), ]

## Grass abundance now
grass <- as.data.frame(fread("grass_litterbag_abundance_all_timepoints.csv"))


## Make treatment names:
grass$Treatment <- ifelse(grepl(pattern = "Phyllosphere", x = grass$`Bag Label`), "Phyllosphere", 
                          ifelse(grepl("RC", grass$`Bag Label`), "Closed", 
                                 ifelse(grepl("RO", grass$`Bag Label`), "Open", 
                                        ifelse(grepl("RU", grass$`Bag Label`), "Overhead", 
                                               ifelse(grepl("RE", grass$`Bag Label`), "Elevated", 
                                                      ifelse(grepl("Litter", grass$`Bag Label`), "Environmental Litter", "unknown"))))))
grass$Month <- ifelse(grass$`Time Point` == 1, "End May", 
                      ifelse(grass$`Time Point` == 2, "End June", 
                             ifelse(grass$`Time Point` == 3, "End July",
                                    ifelse(grass$`Time Point` == 4, "Mid September", 
                                           ifelse(grass$`Time Point` == 5, "End October", 
                                                  ifelse(grepl(pattern = "Phyllosphere", x = grass$`Bag Label`), "Before Dispersal", "error... :("))))))

just_bags <- grass[!(grass$Treatment %in% c("Environmental Litter", "Phyllosphere")), ]

just_bags$Replicate <- gsub("R[A-Z]T[1-5]", "", just_bags$`Bag Label`)

## Merge the two together
glass_red$Treat_Time_Plot <- paste0(glass_red$Treatment, "_T", glass_red$`Time Point`, "_", glass_red$Replicate)
just_bags$Treat_Time_Plot <- paste0(just_bags$Treatment, "_T", just_bags$`Time Point`, "_", just_bags$Replicate)

together <- merge(x = glass_red, y = just_bags, by = "Treat_Time_Plot")

## Statistics: correlate the two!
summary(lm(`Abundance Per G` ~ Abundance_per_cm2, together)) # P = 0.07377



##### _______ immigration rate ________________________________________ ######

# Load data
glass <- as.data.frame(fread("glass_slides_abundance_all_timepoints.csv"))
glass <- glass[, 1:3]
death <- as.data.frame(fread("death_slides_set1_abundance_all_timepoints.csv"))

death$`Time Point` <- as.character(death$`Time Point`)
days.key <- c("0" = 0, "1" = 39, "2" = 60, "3" = 100, "4" = 151, "5" = 195)
death$Time_as_days <- days.key[death$`Time Point`]



##### To calculate immigration rate, we use the equation: I = r * N(^)
# we use this eq because N is at equilibrium so we can assume that death rate equals immigration rate
# r = death rate, which is the slope of the linear model when we natural log the abundance on the death slides
# because we are using the continuous rate equation: P(t) = P(o) * e^(rt) where r = death rate
# I have a sheet explaining these calculations, if needed

# calculating n.hat
open <- glass[grepl("LO", glass$`Bag Label`), ]
aggregate(open$Abundance, by = list(open$`Time Point`), FUN = mean) # checking how accurately the mean Nhat matches the individual time point Nhats (seems *fairly* consistent over time so I'm moving foward with the calculation)
n.hat <- mean(open$Abundance)


# calculating death rate
intercept <- mean(log(death[death$Time_as_days == 0, ]$Abundance)) # because we know for sure what the abundance is at timepoint 0, we want to force the linear model through this intercept

lm_death <- lm(I(log(Abundance) - intercept) ~ Time_as_days + 0, death) # calculates "r" (see eq. above)
r <- -coef(lm_death) # THIS IS THE DEATH RATE!! (we do negative because this is actually describing the continuous rate at which cells survive, not the continuous rate at which cells die)

1 - exp(coef(lm_death)) # 3.34% dies per day; note: this is NOT the death rate value we input into our equilibrium model!! 
1 - exp((coef(lm_death)*7)) # 21.19% dies per week

# calculating immigration rate for the entire slide
I <- n.hat * r
I # 148,327.2 per day 
# note: we can now multiply or divide this number to get cells/week or cells/hour, etc. whatever you are interested in

# now we divide by the area of the glass slide (2.5cm x 7.5 cm) to get immigration rate / day /cm^2
I / (2.5 * 7.5) # 7,910.782 == about 7,900 cells / day / cm^2
(I*7) / (2.5 * 7.5) # 55375.47 == about 55,000 cells / week / cm^2



### calculating the % of total community immigrating each day
grass <- as.data.frame(fread("grass_litterbag_abundance_all_timepoints.csv"))
dry.weight <- as.data.frame(fread("dry_weight_litter_bag__for_immigration_percent_calc.txt"))
combo <- merge(grass, dry.weight, by.x = "Bag Label", by.y = "Sample", all = FALSE)

combo$Abundance_per_litterbag <- combo$`Abundance Per G` * combo$Dry_Mass_g # because abundance is per gram dry weight
combo$Abundance_per_cm2 <- combo$Abundance_per_litterbag / (8*9) # assuming the 10cm x 10cm litterbags have a border of 1cm on three sides where we heat sealed the folded over windowscreen to each other

open.litter <- combo[grepl("RO", combo$`Bag Label`), ]
aggregate(open.litter$Abundance_per_cm2, by = list(open.litter$`Time Point`), FUN = mean) # checking variation, much more variable than environmental litter



I / mean(open.litter$Abundance_per_cm2) # 0.47% of community


##### _______ litter chemistry analysis ___________________ #####

# Step 1: read in data (note: baseline correction already performed on this dataset)
chemdata <- read.csv("litter_chemistry.CSV", row.names=1, comment.char="#", check.names = FALSE) %>% 
  select(!.row) %>% 
  filter(grepl("ro|ru|rc", sample.name)) %>%  # remove environmental , phyllosphere, and Elevated samples
  column_to_rownames("sample.name")


# Step 2: Find the largest negative number + add to each observance
add.value <- -min(chemdata)
chemdata.pos <- chemdata + add.value


# Step 3: standardized across all samples and check to see if it worked
chem.stand <- as.data.frame(t(apply(chemdata.pos, 1, function(x) {x}/sum(x)))) # within each sample adjust peaks proportionally (sum = 1)
apply(chem.stand, 1, sum, na.rm = T) # make sure all your samples sum to "1" 


# Step 4: calculate dissimilarity measure
chem.dist <- dist(chem.stand, method = "euclidean", diag = FALSE, upper = FALSE) # calculate euclidean distance matrix
chem.dist.df <- as.data.frame(as.matrix(chem.dist)) # save to dataframe (you will need this later)

# And do it again but scale them
chem.dist.scale <- dist(scale(chem.stand), method = "euclidean", diag = FALSE, upper = FALSE) # calculate euclidean distance matrix
chem.dist.scale.df <- as.data.frame(as.matrix(chem.dist.scale)) # save to dataframe (you will need this later)


# Step 5: create metadata (I could also read in a file here, but this is easier for me)
meta <- data.frame("SampleID" = row.names(chem.dist.df))
meta$TimePoint <- str_extract(meta$SampleID, "t[0-5]")
meta$Treatment <- str_extract(meta$SampleID, "^r[a-z]")
meta$Treatment <- ifelse(grepl("phy", meta$SampleID), "Phyllosphere", 
                         ifelse(grepl("litter", meta$SampleID), "Environmental Litter", meta$Treatment))
meta$Replicate <- str_extract(meta$SampleID, "r[1-9]")


# Step 6: do some stats!
set.seed(42)
adonis(chem.dist.scale.df ~ as.factor(Treatment), permutations = 10000, 
       data = meta, na.rm= TRUE)


nmds <- metaMDS(chem.stand.less, k = 2)
coords <- nmds$

  
# Step 7: 
pca.obj <- prcomp(chem.stand, center = TRUE, scale = TRUE)
#plot(pca.obj)
pca.meta <- merge(meta, pca.obj$x, by.x = "SampleID", by.y = "row.names")

ggplot(data = pca.meta) + geom_point(aes(x = PC1, y = PC2, color = Treatment), size = 3) + 
  theme_classic() + scale_colour_brewer(type = "qual", palette = 6) 

pca.obj <- prcomp(chem.stand.less, center = TRUE, scale = TRUE)
#plot(pca.obj)
pca.meta <- merge(meta, pca.obj$x, by.x = "SampleID", by.y = "row.names")
