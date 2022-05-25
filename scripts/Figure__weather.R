# Analyze weather data from weather station near LR
# Feb 20 2019
# Kendra E. Walters

############_______________ Load data and libraries
require(ggplot2)
require(data.table)
require(tidyr)

weather <- as.data.frame(fread("East_Loma_weather_data_2018.csv"))




################_______________________ Manipulate the data frame______________#################
## Make the data frame readable and subset columns I need
weather <- weather[, c("Timestamp", " \"\"Temp F\"\"", " \"\"Wind_spd\"\"", " \"\"Wind_dir\"\"", " \"\"Rain_mm\"\"")]
colnames(weather) <- c("Timestamp", "Temp_F", "Wind_Speed", "Wind_Direction", "Rain_mm")

## Subset the dates I need: 
weather.2 <- separate(data = weather, col = Timestamp, into = c("Date", "Time"), sep = ", ")
weather.2$Date <- as.Date(weather.2$Date, "%m/%d/%Y")

use <- weather.2[weather.2$Date >= "2018-04-14" & weather.2$Date <= "2018-10-26",]

ggplot(data = use, color = "black") +
  geom_line(aes(x = Date, y = Wind_Direction))+
  theme_classic() + ylab("Wind direction") + xlab("")  +
  scale_x_date(date_breaks = "months" , date_labels = "%b-%y")

## Average or sum values per day, then merge the data
rain.sum <- aggregate(use$Rain_mm, by=list(Category=use$Date), FUN=sum)
colnames(rain.sum) <- c("Date", "Rain.Sum.mm")

other.avg <- aggregate(list(use$Temp_F, use$Wind_Speed), by=list(Category=use$Date), FUN=mean)
colnames(other.avg) <- c("Date", "Temp.F.avg", "Wind.Speed.avg")

nope <- use[!(use$Wind_Direction == 0), ]
direction <- aggregate(list(nope$Wind_Direction), by=list(Category=nope$Date), FUN=mean)
colnames(direction) <- c("Date", "Wind.Dir.avg")

max <- aggregate(list(use$Temp_F, use$Wind_Speed), by=list(Category=use$Date), FUN=max)
colnames(max) <- c("Date", "Temp.F.max", "Wind.Speed.max")

together <- merge(x = rain.sum, y = other.avg)
together <- merge(together, max)
together <- merge(together, direction)

together$Month <- ifelse(together$Date >= "2018-04-14" & together$Date < "2018-05-23", "End May", 
                         ifelse(together$Date >= "2018-05-23" & together$Date < "2018-06-13", "End June", 
                                ifelse(together$Date >= "2018-06-13" & together$Date < "2018-07-23", "End July", 
                                       ifelse(together$Date >= "2018-07-23" & together$Date < "2018-09-12", "Mid September", 
                                              ifelse(together$Date >= "2018-09-12" & together$Date < "2018-10-26", "End October", "OOOOoooooh... Noooo....")))))
together <- together[!(together$Month == "OOOOoooooh... Noooo...."), ]
together$Temp.C.avg <- (together$Temp.F.avg - 32) * (5/9)
together$Temp.C.max <- (together$Temp.F.max - 32) * (5/9)

together$Wind.Speed.km.avg <- (together$Wind.Speed.avg) * (5/3.1)
together$Wind.Speed.km.max <- (together$Wind.Speed.max) * (5/3.1)



## Get rain sum per timepoint: 
t1 <- together[together$Date >= "2018-04-14" & together$Date < "2018-05-23", ]
t2 <- together[together$Date >= "2018-05-23" & together$Date < "2018-06-13", ]
t3 <- together[together$Date >= "2018-06-13" & together$Date < "2018-07-23", ]
t4 <- together[together$Date >= "2018-07-23" & together$Date < "2018-09-12", ]
t5 <- together[together$Date >= "2018-09-12" & together$Date < "2018-10-26", ]

rain.v <- c(sum(t1$Rain.Sum.mm), sum(t2$Rain.Sum.mm), sum(t3$Rain.Sum.mm), sum(t4$Rain.Sum.mm), sum(t5$Rain.Sum.mm))


## ###############################______________________PLOT THINGS________________________###################
ggplot(data = together, color = "black") + 
  geom_bar(aes(x = Date, y = Rain.Sum.mm), stat = "identity", fill = "steelblue") +
  theme_classic() + ylab("Precipitation (mm)") + xlab("") + 
  scale_x_date(date_breaks = "months" , date_labels = "%b-%y") + 
  scale_y_continuous(expand = c(0,0))

ggplot(data = together, color = "black") + 
  geom_line(aes(x = Date, y = Wind.Speed.km.max), color = "grey60", size = 1) +
  theme_classic() + ylab("Maximum Wind Speed Per Day (km/h)") + xlab("")  + 
  scale_x_date(date_breaks = "months" , date_labels = "%b-%y") + 
  coord_cartesian(xlim = as.Date(c("2018-04-14", "2018-10-25")))

ggplot(data = together, color = "black") + 
  geom_line(aes(x = Date, y = Temp.C.avg), color = "firebrick", size = 1) +
  theme_classic() + ylab("Temperature (Degrees C)") + xlab("")  + 
  scale_x_date(date_breaks = "months" , date_labels = "%b-%y")
