# Script specifically for line graph of the DIFFERENCE between treatments (getting at dispersal routes) of MASS LOSS
# Kendra E. Walters
# June 1st 2020

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


## Calculating variance stats + assigning days for x-axis
days.key <- c("End May" = 39, "End June" = 60, "End July" = 100, "Mid September" = 151, "End October" = 195)
data.to.plot <- summarySE(data = mass.outdf, 
                          measurevar = "Percent_Mass_Loss", 
                          groupvars = c("Treatment", "Month")) %>% 
  mutate(Treatment = factor(Treatment),
         Days = days.key[Month])



## To make the line graph
cols <- c("Air + Vegetation"="#50AB47","Air"="#1B51A1","Soil"="#E68D24")
ggplot(data.to.plot, aes(x = Days, y = Percent_Mass_Loss, group = Treatment, color = Treatment)) + 
  geom_hline(yintercept = 0, color = "grey50", linetype = "longdash", size = 0.75) +
  geom_line(position=pd, size = 0.75) + 
  geom_point(position=pd, size = 3) +
  geom_errorbar(aes(ymin = Percent_Mass_Loss - ci, ymax = Percent_Mass_Loss + ci), 
                position = pd) +
  theme_classic() +
  ylab("Difference in Percent Mass Loss (relative to Closed treatment)") + xlab("Days") +
  theme(axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black")) + 
  scale_color_manual(values = cols, name = "Dispersal Route")