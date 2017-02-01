# ------------------------------------------------------------------------------
# Borrelia Vaccine gypsy moth pupae script
# 10 January 2017
# Modified: 
# Input:  GCPupae2015.csv
#         HCPupae2015.csv
#         TCPupae2015.csv
#         GXPupae2015.csv
#         HXPupae2015.csv
#         TXPupae2015.csv
#         PupMDensity2015.csv
# Output: Graphs
#
# ------------------------------------------------------------------------------
# Reset R's brain
rm(list=ls())

# Find where R is currently working
getwd()

# Reset to desired working directory
setwd("c:/Users/Felicia Keesing/Documents/Data & Figures/Vaccination")

# Find where R is currently working
getwd()


library("ggplot2")

####################################################################
# PART ONE: 2015 data organization and combination
####################################################################

# read in the gypsy moth pupae datafiles
GC<-read.csv("GCPupae2015.csv")
HC<-read.csv("HCPupae2015.csv")
TC<-read.csv("TCPupae2015.csv")
GX<-read.csv("GXPupae2015.csv")
HX<-read.csv("HXPupae2015.csv")
TX<-read.csv("TXPupae2015.csv")
v1 <- rbind(GC, HC)
v2 <- rbind(v1,TC)
v3 <- rbind(GX,v2)
v4 <- rbind(v3,HX)
v5 <- rbind(v4,TX)

head(v5)
dim(v5)

# get the dates sorted out
v5$NewDate <- as.character(v5$Date,format="%m/%d/%y") # convert imported date into characters
v5$date <- as.Date(v5$NewDate, format="%m/%d/%Y") # convert characters into date recognized by R
v5$NewDate <- NULL
v5$Date <- NULL
head(v5)

# split "Site" into two variables: "Location" and "Treatment"
v5[,16] <- data.frame(Location =   substr(v5$Site, start = 1, stop = 1))
v5[,17] <- data.frame(Treatment = substr(v5$Site, start = 2, stop = 2) )
#head(v5)
p <- v5

#write.csv(p,"AllGridsPup2015.csv") # used this to find some missing values in original files
# some data entry omitted 0 values after a pupa was found missing at an earlier check

# aggregate by grid, adding up all the remaining pupae to get Number Remaining
p.grid <- aggregate(p$Day.7, by=list(p$Site), FUN=sum)
names(p.grid) <- c("Site", "NumRemaining")
p.grid <- aggregate(p$Day.7, by=list(p$Location, p$Treatment, p$Site), FUN=sum)
names(p.grid) <- c("Location", "Treatment", "Site", "NumRemaining")
p.grid$PropEaten <- 1-p.grid$NumRemaining/81 #calculate the proportion eaten

#read in mouse densities on each grid
d <- read.csv("PupMDensity2015.csv")

# merge to get master result
p.dens <- merge(p.grid, d, by="Site")


# THIS SECTION TREATS EACH GRID AS UNIT OF REPLICATION
mean.Man <- tapply(p.grid$PropEaten, list(p.grid$Treatment), mean)
# Calculate standard deviation of each group using tapply
# which returns a matrix--it will be the same dimension as mean.burde 
sd.Man <- tapply(p.grid$PropEaten, list(p.grid$Treatment), sd)
# Calculate sample size in each group
n.Man <- tapply(p.grid$PropEaten, list(p.grid$Treatment), length)
# Calculate the standard error 
se.Man <- sd.Man/sqrt(n.Man)
par(mfrow = c(1,2))
# Now make the barplot
# beside=TRUE makes it a grouped barplot
mids<- barplot(mean.Man,
               xlab = "Grid treatment", 
               ylab = "Mean proportion of pupae eaten",
               ylim = c(0,1),
               # col=colors,
               beside=TRUE)
# use arrows to put error bars on the plot
arrows(mids, mean.Man - se.Man, mids, mean.Man+se.Man, code = 3, angle = 90, length = 0.1)
text(mids, .1, paste ("n = ", n.Man), cex=.7, col="white")

#plot the proportion eaten versus # of animals on that grid
plot(p.dens$PropEaten ~ p.dens$MNA, xlim=c(0,85), ylim=c(0,.9),
     ylab="Proportion Pupae Eaten", xlab="MNA")
#abline(lm(p.dens$PropEaten ~ p.dens$MNA))
anova(lm(p.dens$PropEaten ~ p.dens$MNA))
text(p.dens$MNA, p.dens$PropEaten-.05,labels=p.dens$Site)

# DATA ANALYSIS HERE
# Need to incorporate density of mice on plots as covariate
# Need to incorporate paired design of grids 












shapiro.test(p.grid$PropEaten)
hist(p.grid$PropEaten)

# analysis of covariance with density, but loses paired design
n <- aov(p.dens$PropEaten ~ p.dens$Treatment.x+p.dens$MNA)
plot(lm(p.grid$PropEaten ~ p.grid$Treatment))


################################################################
# JUST LOOK AT THE ONES WITH RODENT SIGN
################################################################

# subset the ones with rodent sign
p.sign <- subset(p, Sign == "yes")
# aggregate by grid, count number with sign on each grid
p.signcount <- aggregate(p.sign$Day.7, by=list(p.sign$Location, p.sign$Treatment, p.sign$Site), FUN=length)
names(p.signcount) <- c("Location", "Treatment", "Site", "NumWithSign")
#  
p.signsum <- aggregate(p.sign$Day.7, by=list(p.sign$Location, p.sign$Treatment, p.sign$Site), FUN=sum)
names(p.signsum) <- c("Location", "Treatment", "Site", "NumRemaining")
# merge the two to get the count and sum in one data file
p.sign.grid <- merge(p.signsum, p.signcount, by="Site")
p.sign.grid[,5:6] <- NULL
p.sign.grid$PropEaten <- 1-p.sign.grid$NumRemaining/p.sign.grid$NumWithSign #calculate the proportion eaten

#read in mouse densities on each grid
d <- read.csv("PupMDensity2015.csv")

p.sign.dens <- merge(p.sign.grid, d, by="Site")


# THIS SECTION TREATS EACH GRID AS UNIT OF REPLICATION, but only for pupae with rodent sign
mean.Man <- tapply(p.sign.grid$PropEaten, list(p.sign.grid$Treatment), mean)
# Calculate standard deviation of each group using tapply
# which returns a matrix--it will be the same dimension as mean.burde 
sd.Man <- tapply(p.sign.grid$PropEaten, list(p.sign.grid$Treatment), sd)
# Calculate sample size in each group
n.Man <- tapply(p.sign.grid$PropEaten, list(p.sign.grid$Treatment), length)
# Calculate the standard error 
se.Man <- sd.Man/sqrt(n.Man)
par(mfrow = c(1,2))
# Now make the barplot
# beside=TRUE makes it a grouped barplot
mids<- barplot(mean.Man,
               xlab = "Grid treatment", 
               ylab = "Mean proportion of pupae eaten",
               ylim = c(0,1),
               # col=colors,
               beside=TRUE)
# use arrows to put error bars on the plot
arrows(mids, mean.Man - se.Man, mids, mean.Man+se.Man, code = 3, angle = 90, length = 0.1)
text(mids, .1, paste ("n = ", n.Man), cex=.7, col="white")



# THIS SECTION TREATS EACH GRID AS UNIT OF REPLICATION, but only for pupae with rodent sign
mean.Man <- tapply(p.sign.grid$NumWithSign, list(p.sign.grid$Treatment), mean)
# Calculate standard deviation of each group using tapply
# which returns a matrix--it will be the same dimension as mean.burde 
sd.Man <- tapply(p.sign.grid$NumWithSign, list(p.sign.grid$Treatment), sd)
# Calculate sample size in each group
n.Man <- tapply(p.sign.grid$NumWithSign, list(p.sign.grid$Treatment), length)
# Calculate the standard error 
se.Man <- sd.Man/sqrt(n.Man)
par(mfrow = c(1,2))
# Now make the barplot
# beside=TRUE makes it a grouped barplot
mids<- barplot(mean.Man,
               xlab = "Grid treatment", 
               ylab = "Mean number of pupae with rodent sign",
               ylim = c(0,70),
               # col=colors,
               beside=TRUE)
# use arrows to put error bars on the plot
arrows(mids, mean.Man - se.Man, mids, mean.Man+se.Man, code = 3, angle = 90, length = 0.1)
text(mids, 5, paste ("n = ", n.Man), cex=.7, col="white")


#plot the proportion eaten versus # of animals on that grid
plot(p.sign.dens$NumWithSign ~ p.sign.dens$MNA, xlim=c(0,85), ylim=c(0,85),
     ylab="Number of pupae with rodent sign", xlab="MNA", pch=16)
text(p.sign.dens$MNA, p.sign.dens$NumWithSign-5,labels=p.sign.dens$Site)
#abline(lm(p.dens$PropEaten ~ p.dens$MNA))
anova(lm(p.dens$PropEaten ~ p.dens$MNA))

# DATA ANALYSIS HERE
# Need to incorporate density of mice on plots as covariate
# Need to incorporate paired design of grids 













shapiro.test(p.grid$PropEaten)
hist(p.grid$PropEaten)

# analysis of covariance with density, but loses paired design
n <- aov(p.dens$PropEaten ~ p.dens$Treatment.x+p.dens$MNA)
plot(lm(p.grid$PropEaten ~ p.grid$Treatment))



