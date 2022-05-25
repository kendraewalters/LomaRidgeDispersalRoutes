# Abundance on glass slides averaged for each dispersal route AND then split up by dispersal route and time
# Kendra E. Walters
# April 6th 2020

# Room of ()ments
require(data.table)
require(ggplot2)
require(tidyverse)


##### _____________ SETTING UP FUNCTIONS _________________ #####


## From: http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}


#####


# Load data
glass <- as.data.frame(fread("glass_slides_abundance_all_timepoints.csv"))

# Make treatment names:
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



########## __________________ Make a difference between group graph but use the paired samples within a plot
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


## Calculating variance stats + assigning days for x-axis
days.key <- c("End May" = 39, "End June" = 60, "End July" = 100, "Mid September" = 151, "End October" = 195)
data.to.plot <- summarySE(data = out.df, 
                          measurevar = "Abundance_per_cm2", 
                          groupvars = c("Treatment", "Month")) %>% 
  mutate(Treatment = factor(Treatment, levels = c("Air", "Vegetation", "Soil")),
         Days = days.key[Month])

t_not_0$Treatment <- factor(t_not_0$Treatment, levels = c("Environmental Litter", "Closed", 
                                                          "Elevated", "Overhead", 
                                                          "Open"))
## To make the line graph
cols <- c("Vegetation"="#50AB47","Air"="#1B51A1","Soil"="#E68D24")

ggplot(data.to.plot, aes(x = Days, y = Abundance_per_cm2, group = Treatment, color = Treatment)) + 
  geom_hline(yintercept = 0, color = "grey50", linetype = "longdash", size = 0.75) +
  geom_line(position=pd, size = 0.75) + 
  geom_point(position=pd, size = 3) +
  geom_errorbar(aes(ymin = Abundance_per_cm2 - ci, ymax = Abundance_per_cm2 + ci), 
                position = pd) +
  theme_classic() +
  ylab("Difference in Bacterial Abundance / cm2 (relative to Closed treatment)") + xlab("Days") +
  theme(axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black")) + 
  scale_color_manual(values = cols, name = "Dispersal Route")