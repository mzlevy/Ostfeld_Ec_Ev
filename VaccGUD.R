# ------------------------------------------------------------------------------
# Borrelia Vaccine GUD Trays script
# 15 January 2017
# Modified: 26 January 2017
# Input: GUD data total 20161101.csv
# Output: 
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

##################################################################
# DEFINE THE SE function
#
# Function to calculate the mean and the standard deviation
# for each group
#+++++++++++++++++++++++++
# data : a data frame
# varname : the name of a column containing the variable
#to be summariezed
# groupnames : vector of column names to be used as
# grouping variables
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}



####################################################################
# GUD TRAYS 2015
####################################################################

g<-read.csv("GUD data total 20161101.csv")
head(g)
dim(g)

# remove trays that flipped over, etc. and focus just on 2015
g1 <- subset(g,Problem=="N",select=c(1:15)) # get rid of any data with problems (e.g. tray flipped over)
g2 <- subset(g1, Year=="2015", select=c(1:15)) # select just data from 2015

# this aggregates by date and strips down to the key column -- CORRECTED seeds remaining
# Note that the *corrected* seeds remaining is a maximum of 4 g so that we have a conservative
# estimate of the weight of seeds removed.
# The resulting number for CorrSeedRem is the mean of the several days of measurements at each tray.
g3 <- aggregate(g2[,10], by=list(g2$Grid, g2$Site, g2$Treatment, g2$Cover, g2$Trap), FUN=mean, ra.nm=FALSE)
colnames(g3)[1:6]<- c("Grid", "Site", "Treatment", "Cover", "Trap", "CorrSeedRem")

################################################################################
# CALCULATES PLOT AVERAGES WITHOUT USING PAIRED DESIGN OF SHADED/UNSHADED TRAYS

# this section aggregates by treatment, date, grid, but without using the paired design of shaded/unshaded at each trap location
# aggregate by site-treatment, e.g. the mean of GX covered trays
g4 <- aggregate(g3[,6], by=list(g3$Site, g3$Treatment, g3$Cover), FUN=mean, ra.nm=FALSE)
colnames(g4)[1:4]<- c("Site", "Treatment", "Cover", "CorrSeedRem")
head(g4)


g4summ <- data_summary(g4, varname="CorrSeedRem", 
                       groupnames=c("Treatment", "Cover"))
g4summ$se <- g4summ$sd/sqrt(3)
# Convert dose to a factor variable
g4summ$Treatment=as.factor(g4summ$Treatment)
head(g4summ)

ggplot(g4summ, aes(x=Treatment, y=CorrSeedRem, fill=Cover)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=CorrSeedRem-se, ymax=CorrSeedRem+se), width=.2,
                position=position_dodge(.9))+
  title("Ignoring paired tray design")



################################################################################
# CALCULATES PLOT AVERAGES USING PAIRED DESIGN OF SHADED/UNSHADED TRAYS


r<-read.csv("GUD data total 20161101.csv")
head(r)
dim(r)

# remove trays that flipped over, etc. and focus just on 2015
r1 <- subset(g,Problem=="N",select=c(1:15)) # get rid of any data with problems (e.g. tray flipped over)
r2 <- subset(g1, Year=="2015", select=c(1:15)) # select just data from 2015
r2$SeedEaten <- NULL # This variable is wrong -- makes assumptions about starting weight of seeds.
r2$PitReader <- NULL
r2$GUDTray <- NULL
r2$SeedRemaining <- NULL
r2$Problem <- NULL
r2$Notes <- NULL

# this aggregates by date and strips down to the key column -- CORRECTED seeds remaining
# Note that the *corrected* seeds remaining is a maximum of 4 g so that we have a conservative
# estimate of the weight of seeds removed.
# The resulting number for CorrSeedRem is the mean of the several days of measurements at each tray.
#g3 <- aggregate(g2[,10], by=list(g2$Grid, g2$Site, g2$Treatment, g2$Cover, g2$Trap), FUN=mean, ra.nm=FALSE)
#colnames(g3)[1:6]<- c("Grid", "Site", "Treatment", "Cover", "Trap", "CorrSeedRem")

# this section takes advantage of the paired shaded-unshaded trays at each location
gS <- subset(r2,Cover=="S") # subset shaded trays
gU <- subset(r2,Cover=="U") # subset unshaded trays
# merge these two to get paired design by traplocation
gPaired <- merge(gS, gU, by=c("Grid", "Trap", "Site", "Treatment", "Date"))
gPaired[,6:7] <- NULL # remove categorical coding as "S" for shaded
gPaired[,8:9] <- NULL # remove categorical coding as "U" for uncovered
colnames(gPaired)[7] <- c("SeedRemS")
colnames(gPaired)[9] <- c("SeedRemU")
# Calculate how much more there was in uncovered tray
gPaired$SeedDiff <- gPaired$SeedRemU - gPaired$SeedRemS 
# Set the Dist variable, which notes whether either of the pairs was disturbed on a particular night. 
# If it is >0, then at least one was disturbed.
gPaired$Dist <- gPaired$Disturbed.x+gPaired$Disturbed.y
# now subset just those nights where at least one of the two trays was disturbed
gPaired2 <- subset(gPaired, Dist != 0)
head(gPaired2)

# Take the mean of the SeedDiff across dates at each tray pair
g3 <- aggregate(gPaired2[,10], by=list(gPaired2$Grid, gPaired2$Site, 
                                      gPaired2$Treatment), FUN=mean, ra.nm=FALSE)
colnames(g3)[1:4] <- c("Grid", "Site", "Treatment", "SeedDiff")

# STATISTICAL ANALYSIS HERE






# GRAPH MEAN OF DIFFERENCE BETWEEN SHADED AND UNSHADED, ACROSS ALL GRIDS
mean.Man <- tapply(gPaired2$SeedDiff, list(gPaired2$Treatment), mean)
# Calculate standard deviation of each group using tapply
# which returns a matrix--it will be the same dimension as mean.burde 
sd.Man <- tapply(gPaired2$SeedDiff, list(gPaired2$Treatment), sd)
# Calculate sample size in each group
n.Man <- tapply(gPaired2$SeedDiff, list(gPaired2$Treatment), length)
# Calculate the standard error 
se.Man <- sd.Man/sqrt(n.Man)
par(mfrow = c(1,2))
# Now make the barplot
# beside=TRUE makes it a grouped barplot
mids<- barplot(mean.Man,
               xlab = "Grid treatment", 
               ylab = "Mean difference in wt. of seeds remaining (g)",
               ylim = c(0,.5))
# use arrows to put error bars on the plot
arrows(mids, mean.Man - se.Man, mids, mean.Man+se.Man, code = 3, angle = 90, length = 0.1)
text(mids, .02, paste ("n = ", n.Man), cex=.7, col="white")



# GRAPH EACH GRID MEAN OF DIFFERENCE BETWEEN SHADED AND UNSHADED
mean.Man <- tapply(gPaired2$SeedDiff, list(gPaired2$Treatment, gPaired2$Site), mean)
# Calculate standard deviation of each group using tapply
# which returns a matrix--it will be the same dimension as mean.burde 
sd.Man <- tapply(gPaired2$SeedDiff, list(gPaired2$Treatment, gPaired2$Site), sd)
# Calculate sample size in each group
n.Man <- tapply(gPaired2$SeedDiff, list(gPaired2$Treatment, gPaired2$Site), length)
# Calculate the standard error 
se.Man <- sd.Man/sqrt(n.Man)
par(mfrow = c(1,2))
# Now make the barplot
# beside=TRUE makes it a grouped barplot
mids<- barplot(mean.Man,
               xlab = "Grid treatment", 
               ylab = "Mean difference in wt. of seeds remaining (g)",
               ylim = c(-1,1),
               # col=colors,
               beside=TRUE, legend=TRUE)
# use arrows to put error bars on the plot
arrows(mids, mean.Man - se.Man, mids, mean.Man+se.Man, code = 3, angle = 90, length = 0.1)
text(mids, .1, paste ("n = ", n.Man), cex=.7, col="white")







# EFFECTS OF TREATMENT AND SITE ON *DIFFERENCE* BETWEEN COVERED & UNCOVERED TRAYS
gPaired2Summ <- data_summary(gPaired2, varname=c("SeedDiff"), 
                        groupnames=c("Treatment"))
gPaired2Summ$se <- gPaired2Summ$sd/sqrt(3)
# Convert dose to a factor variable
gPSumm3$Treatment=as.factor(gPSumm3$Treatment)
head(gP3Summ)

p <- ggplot(gP3Summ, aes(x=Treatment, y=SeedDiff)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=SeedDiff-se, ymax=SeedDiff+se), width=.2,
                position=position_dodge(.9))


p <- ggplot(gPSumm, aes(x=Site, y=SeedDiff, fill=Treatment)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=SeedDiff-se, ymax=SeedDiff+se), width=.2,
                position=position_dodge(.9))
p+ylab("Mean difference in weight of seeds remaining (g)")+
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=12))+
  theme(legend.title=element_blank())
  
  
  
#####PULL OUT JUST DISTURBED GUD TRAYS #############################################################
#####We define disturbed as either one of the trays at a trap location begin disturbed on a particular
#####night. If neither tray is disturbed on a particular night, we do not include it.
head(g)
gD <- subset(g,Year =="2015" & Problem=="N")
gDS <- subset(gD,Cover=="S",select=c(1:11))
gDU <- subset(gD,Cover=="U",select=c(1:11))
gDPaired <- merge(gDS, gDU, by=c("Grid", "Trap", "Site", "Treatment"))
#remove unnecessary columns
gDPaired[,12:16] <- NULL
gDPaired[,8:9] <- NULL # remove categorical coding as "S" for shaded
gDPaired[,5:6] <- NULL # remove categorical coding as "U" for uncovered
colnames(gDPaired)[5:9] <- c("Date", "DistS", "SeedRemS", "DistU", "SeedRemU")
# calculate difference between shaded and unshaded
gDPaired$SeedDiff <- gDPaired$SeedRemU - gDPaired$SeedRemS # how much more there was in uncovered tray

# now pull out only those dates where at least one of the two paired trays was disturbed that night
gD1 <- subset(gDPaired, DistS =="1" | DistU == "1")

# aggregate taking the mean of the date at each tray location
gD2 <- aggregate(gD1[,7:10], by=list(gD1$Grid, gD1$Site, gD1$Treatment, gD1$Trap), FUN=mean, ra.nm=FALSe)
colnames(gD2)[1:4]<- c("Grid", "Site", "Treatment", "Trap")

# ANOVA on SeedDiff blocked by Site
shapiro.test(gD2$SeedDiff)
gD2$transSeedDiff <- sqrt(gD2$SeedDiff)
shapiro.test(gD2$transSeedDiff)
block <- aov(lm(gD2$transSeedDiff ~ gD2$Treatment + gD2$Site))
summary(block)

# now to graph it
gD2summ <- data_summary(gD2, varname="SeedDiff", 
                        groupnames=c("Treatment", "Site"))
gD2summ$se <- gD2summ$sd/sqrt(3)
# Convert dose to a factor variable
gD2summ$Treatment=as.factor(gD2summ$Treatment)
head(gD2summ)

ggplot(gD2summ, aes(x=Site, y=SeedDiff, fill=Treatment)) + 
  geom_bar(stat="identity", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=SeedDiff-se, ymax=SeedDiff+se), width=.2,
                position=position_dodge(.9))

# ANOVA on SeedRemU with Treatment*Site interaction
shapiro.test(gD2$SeedRemU)
gD2$transSeedRemU <- asin(sqrt(gD2$SeedRemU))
shapiro.test(gD2$transSeedRemU)
hist(gD2$transSeedRemU)
block <- aov(lm(gD2$transSeedRemU ~ gD2$Treatment + gD2$Site))
summary(block)

# ANOVA on SeedRemS with Treatment*Site interaction
shapiro.test(gD2$SeedRemS)
gD2$transSeedRemS <- asin(sqrt(gD2$SeedRemS))
shapiro.test(gD2$transSeedRemS)
hist(gD2$transSeedRemS)
block <- aov(lm(gD2$transSeedRemS ~ gD2$Treatment + gD2$Site))
summary(block)


#####HERE FOCUS ON **COUNT** OF DISTURBED GUD TRAYS ####################################################
#####
#####
gd <- aggregate(g2[,10], by=list(g2$Grid, g2$Site, g2$Treatment, g2$Cover, g2$Trap, g2$Date), FUN=sum, ra.nm=FALSE)
colnames(gd)[1:7]<- c("Grid", "Site", "Treatment", "Cover", "Trap", "Date", "CountDist")

#aggregate across dates for each trap location, so that TotalDist is total number of visits to that tray 
gd4 <- aggregate(gd[,7], by=list(gd$Grid, gd$Site, gd$Treatment, gd$Cover, gd$Trap), FUN=sum, ra.nm=FALSE)
colnames(gd4)[1:6]<- c("Grid", "Site", "Treatment", "Cover", "Trap", "TotalDist")


# Here we are pairing S/U tray pairs to get the total number of disturbances at each location
gd4S <- subset(gd4,Cover=="S") # subset shaded trays
gd4U <- subset(gd4,Cover=="U") # subset unshaded trays
# merge these two to get paired design by trap location
gd4Paired <- merge(gd4S, gd4U, by=c("Grid", "Trap", "Site", "Treatment"))
gd4Paired[,5] <- NULL # remove categorical coding as "S" for shaded
gd4Paired[,6] <- NULL # remove categorical coding as "U" for uncovered
colnames(gd4Paired)[5] <- c("TotalDistS")
colnames(gd4Paired)[6] <- c("TotalDistU")
# Calculate total number of disturbance-nights at each tray pair
gd4Paired$TotalDist <- gd4Paired$TotalDistS + gd4Paired$TotalDistU 
head(gd4Paired)


# GRAPH MEAN # DISTURBED TRAY-NIGHTS FOR EACH TRAY PAIR, BY GRID, 
mean.Man <- tapply(gd4Paired$TotalDist, list(gd4Paired$Treatment, gd4Paired$Site), mean)
# Calculate standard deviation of each group using tapply
# which returns a matrix--it will be the same dimension as mean.burde 
sd.Man <- tapply(gd4Paired$TotalDist, list(gd4Paired$Treatment, gd4Paired$Site), sd)
# Calculate sample size in each group
n.Man <- tapply(gd4Paired$TotalDist, list(gd4Paired$Treatment, gd4Paired$Site), length)
# Calculate the standard error 
se.Man <- sd.Man/sqrt(n.Man)
par(mfrow = c(1,2))
# Now make the barplot
# beside=TRUE makes it a grouped barplot
mids<- barplot(mean.Man,
               xlab = "Grid treatment", 
               ylab = "Mean # of trays disturbed",
               ylim = c(0,10),
               # col=colors,
               beside=TRUE, legend=TRUE)
# use arrows to put error bars on the plot
arrows(mids, mean.Man - se.Man, mids, mean.Man+se.Man, code = 3, angle = 90, length = 0.1)
text(mids, .1, paste ("n = ", n.Man), cex=.7, col="white")


# GRAPH MEAN # DISTURBED TRAY-NIGHTS FOR EACH TRAY PAIR, BY TREATMENT, 
mean.Man <- tapply(gd4Paired$TotalDist, list(gd4Paired$Treatment), mean)
# Calculate standard deviation of each group using tapply
# which returns a matrix--it will be the same dimension as mean.burde 
sd.Man <- tapply(gd4Paired$TotalDist, list(gd4Paired$Treatment), sd)
# Calculate sample size in each group
n.Man <- tapply(gd4Paired$TotalDist, list(gd4Paired$Treatment), length)
# Calculate the standard error 
se.Man <- sd.Man/sqrt(n.Man)
par(mfrow = c(1,2))
# Now make the barplot
# beside=TRUE makes it a grouped barplot
mids<- barplot(mean.Man,
               xlab = "Grid treatment", 
               ylab = "Mean # of tray-pair-nights with disturbance",
               ylim = c(0,6),
               # col=colors,
               beside=TRUE)
# use arrows to put error bars on the plot
arrows(mids, mean.Man - se.Man, mids, mean.Man+se.Man, code = 3, angle = 90, length = 0.1)
text(mids, .3, paste ("n = ", n.Man), cex=.7, col="white")




