# Script to make Figure bacterial abundance in litterbags

# Kendra E Walters
# December 3rd, 2020

# Room of ()ments
require(data.table)
require(ggplot2)
require(tidyverse)

# Read in the data
grass <- as.data.frame(fread("grass_litterbag_abundance_all_timepoints_FOR_R.csv"))


##################################_________________________________ Make treatment names:
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
# Set up color scheme: 

cols <- c("All Dispersal"="#340951","Local Rain/Air"="#612C86","Regional Rain/Air"="#A970D0", "No Dispersal" = "#DFB1FE", "Environmental Litter" = "#E9AC44")
grass$Month <- factor(grass$Month, levels = c("End May", "End June", "End July", "Mid September", "End October"), ordered = TRUE)

# Making a dataset of DIFFERENCES
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

## Calculating variance stats + assigning days for x-axis
days.key <- c("End May" = 39, "End June" = 60, "End July" = 100, "Mid September" = 151, "End October" = 195)
data.to.plot <- summarySE(data = out.df.grass, 
                          measurevar = "Abundance_per_g", 
                          groupvars = c("Treatment", "Month")) %>% 
  mutate(Treatment = factor(Treatment),
         Days = days.key[Month])


## Graph
pd <- position_dodge(6) # move them 6 to the left and right
cols <- c("Vegetation"="#50AB47","Air"="#1B51A1","Soil"="#E68D24")

ggplot(data.to.plot, aes(x = Days, y = Abundance_per_g, group = Treatment, color = Treatment)) + 
  geom_line(position=pd, size = 0.75) + 
  geom_point(position=pd, size = 3) +
  geom_errorbar(aes(ymin = Abundance_per_g - ci, ymax = Abundance_per_g + ci), 
                position = pd) +
  geom_hline(yintercept = 0, color = "grey50", linetype = "longdash") +
  theme_classic() +
  ylab("Difference in Bacterial Abundance / g (relative to Closed treatment") + xlab("Days") +
  theme(axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black")) + 
  scale_color_manual(values = cols)