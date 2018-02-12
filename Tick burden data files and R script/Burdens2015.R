# ------------------------------------------------------------------------------
# Borrelia Vaccine tick burdens script
# 28 October 2016
# Modified: 15 December 2016
# Input: TickBurdens2015.csv, which contains the data on Peromyscus leucopus from experimental and controlgrids at Cary.
#       Xenodiagnosis Jul2015.csv
#       2014MMTrapping.csv
#       2014VaccTX.csv
#       2014VaccHX.csv
#       2014VaccGX.csv
#       GUD data total 20161101.csv
# Output: 
#
# ------------------------------------------------------------------------------
# Reset R's brain
rm(list=ls())

# Find where R is currently working
getwd()

# Reset to desired working directory
#setwd("c:/Users/Felicia Keesing/Documents/Data & Figures/Vaccination")
setwd("~/Vaccine/Tick burden data files and R script")
# Find where R is currently working
getwd()


library("ggplot2")

####################################################################
# PART ONE: TICK BURDENS
####################################################################

# read in the tick burdens datafile
t<-read.csv("TickBurdens2015.csv")
head(t)
dim(t)

# replace "NA" with 0 throughout the dataframe
t[is.na(t)]<- 0

# calculate total tick burden for nymphs and adults
t$LSUM <- t$LRE+t$LLE+t$LHD
t$NSUM <- t$NRE+t$NLE+t$NHD

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

##################################################################
# aggregate data by grid by date; each animal's tick numbers for each period are on one line

t1 <- aggregate(t[,22:31], by=list(t$Period,t$Grid,t$Sex, t$Treatment, t$Location, t$Tag), FUN=mean, ra.nm=FALSE)
colnames(t1)[1:6]<- c("Period","Grid", "Sex", "Treatment", "Location", "Tag")

# subset out just Period 2, when animals had been more fully vaccinated
t2 <- subset(t1, Period == 2,
                  select=c(Treatment, Location, Grid, Period, Tag, Sex, LSUM, NSUM))

# average burdens across all individuals on a grid, differentiating by sex
t3 <- aggregate(t2[,7:8], by=list(t2$Grid, t2$Sex, t2$Treatment, t2$Location), FUN=mean, ra.nm=FALSE)
colnames(t3)[1:4]<- c("Grid","Sex", "Treatment", "Location")

# average across males and females to get overall average tick burdens, regardless of sex
t5 <- aggregate(t3[,5:6], by=list(t3$Grid, t3$Treatment, t3$Location), FUN=mean, ra.nm=FALSE)
colnames(t5)[1]<- "Grid"
colnames(t5)[2]<- "Treatment"
colnames(t5)[3]<- "Location"

# THIS SECTION TREATS EACH MOUSE AS UNIT OF REPLICATION
# FIGURE 5 BOXPLOT
#
# Overlaid histograms
ggplot(t2, aes(x=LSUM, fill=Treatment)) +
  geom_density(alpha=.3)

t2$Sex <-as.character(t2$Sex)

a <- ggplot(t2, aes(x=Sex, y=LSUM, fill=Treatment)) + 
  geom_boxplot(position = position_dodge()) +
  guides(fill=FALSE)+
 # stat_summary(fun.y=mean, geom="point", shape=19, size=2)+
  labs(x="Sex",y="Larval burden per mouse")
plot_grid(a, ncol=2, nrow =1) # puts into one column

mean.Man <- tapply(t2$LSUM, list(t2$Sex), mean)
# Calculate standard deviation of each group using tapply
# which returns a matrix--it will be the same dimension as mean.burde 
sd.Man <- tapply(t2$LSUM, list(t2$Sex), sd)
# Calculate sample size in each group
n.Man <- tapply(t2$LSUM, list(t2$Sex), length)
# Calculate the standard error 
se.Man <- sd.Man/sqrt(n.Man)
par(mfrow = c(1,2))
# Now make the barplot
# beside=TRUE makes it a grouped barplot
mids<- barplot(mean.Man,
               xlab = "Sex", 
               ylab = "Mean larval burden per mouse",
               ylim = c(0,15),
              # col=colors,
               beside=TRUE)
# use arrows to put error bars on the plot
arrows(mids, mean.Man - se.Man, mids, mean.Man+se.Man, code = 3, angle = 90, length = 0.1)
text(mids, .5, paste ("n = ", n.Man), cex=.7, col="white")

shapiro.test(t2$LSUM)
hist(t2$LSUM)
t2$transLSUM <- log(t2$LSUM+1)
shapiro.test(t2$transLSUM)
hist(t2$transLSUM)
n <- anova(lm(t2$transLSUM ~ t2$Treatment * t2$Sex))
n
plot(lm(t2$transLSUM ~ t2$Treatment * t2$Sex))

# export t2 to be used with xeno results in VaccXeno2015.R
write.csv(t2,"2015TickSumm.csv")

  
#####################################################################################3
#   NOW SUBSET OUT JUST MALES TO TEST FOR BLOCK EFFECT AND TREATMENT EFFECT
#
t2m <- subset(t2,Sex=="1")
# THIS SECTION TREATS EACH MOUSE AS UNIT OF REPLICATION
mean.Man <- tapply(t2m$LSUM, list(t2m$Treatment, t2m$Location), mean)
# Calculate standard deviation of each group using tapply
# which returns a matrix--it will be the same dimension as mean.burde 
sd.Man <- tapply(t2m$LSUM, list(t2m$Treatment, t2m$Location), sd)
# Calculate sample size in each group
n.Man <- tapply(t2m$LSUM, list(t2m$Treatment, t2m$Location), length)
# Calculate the standard error 
se.Man <- sd.Man/sqrt(n.Man)
par(mfrow = c(1,2))
# Now make the barplot
# beside=TRUE makes it a grouped barplot
mids<- barplot(mean.Man,
               xlab = "Location", 
               ylab = "Mean larval burden per male mouse",
               ylim = c(0,20),
               # col=colors, 
               legend=TRUE, beside=TRUE, args.legend = list(x="topleft"))
# use arrows to put error bars on the plot
arrows(mids, mean.Man - se.Man, mids, mean.Man+se.Man, code = 3, angle = 90, length = 0.1)
text(mids, .5, paste ("n = ", n.Man), cex=.7, col="white")

shapiro.test(t2m$LSUM)
hist(t2m$LSUM)
t2m$transLSUM <- log(t2m$LSUM+1)
shapiro.test(t2m$transLSUM)
hist(t2m$transLSUM)
n <- anova(lm(t2m$transLSUM ~ t2m$Treatment * t2m$Location))
plot(lm(t2m$transLSUM ~ t2m$Treatment * t2m$Sex))

#####################################################################################3
#   NOW SUBSET OUT JUST FEMALES TO TEST FOR BLOCK EFFECT AND TREATMENT EFFECT
#
t2f <- subset(t2,Sex=="2")
# THIS SECTION TREATS EACH MOUSE AS UNIT OF REPLICATION
mean.Man <- tapply(t2f$LSUM, list(t2f$Treatment, t2f$Location), mean)
# Calculate standard deviation of each group using tapply
# which returns a matrix--it will be the same dimension as mean.burde 
sd.Man <- tapply(t2f$LSUM, list(t2f$Treatment, t2f$Location), sd)
# Calculate sample size in each group
n.Man <- tapply(t2f$LSUM, list(t2f$Treatment, t2f$Location), length)
# Calculate the standard error 
se.Man <- sd.Man/sqrt(n.Man)
par(mfrow = c(1,2))
# Now make the barplot
# beside=TRUE makes it a grouped barplot
mids<- barplot(mean.Man,
               xlab = "Location", 
               ylab = "Mean larval burden per female mouse",
               ylim = c(0,20),
               # col=colors, 
               legend=TRUE, beside=TRUE, args.legend = list(x="topleft"))
# use arrows to put error bars on the plot
arrows(mids, mean.Man - se.Man, mids, mean.Man+se.Man, code = 3, angle = 90, length = 0.1)
text(mids, .5, paste ("n = ", n.Man), cex=.7, col="white")

shapiro.test(t2f$LSUM)
hist(t2f$LSUM)
t2f$transLSUM <- log(t2f$LSUM+1)
shapiro.test(t2f$transLSUM)
hist(t2f$transLSUM)
n <- anova(lm(t2f$transLSUM ~ t2f$Treatment * t2f$Location))
plot(lm(t2f$transLSUM ~ t2f$Treatment * t2f$Sex))





# statistics on LSUM by male/female * vaccinated/sham
shapiro.test(t2$LSUM) # testing for need for transformation
hist(t2$LSUM)
t2$transLSUM <- log(t2$LSUM+1)
shapiro.test(t2$transLSUM)
hist(t2$transLSUM)
# a priori prediction of males<females and alpha = 0.1
int <- anova(lm(t2$transLSUM ~ t2$Treatment * t2$Sex))
plot(lm(t2$transLSUM ~ t2$Treatment * t2$Sex))
int
noint <- anova(lm(t2$transLSUM ~ t2$Treatment + t2$Sex))
plot(lm(t2$transLSUM ~ t2$Treatment + t2$Sex))
noint
# sig difference of sex and marginally of treatment with a priori predictions



#################################################################################################
#LOOKING JUST AT WHETHER MALE/FEMALE LARVAL BURDENS DIFFER BETWEEN GRIDS, REMOVING VACCINATED MICE

# Now subset out treatment and just look at LSUM vs grid and male/female
m <- subset(t2, Treatment=="C") # only mice on control grids
msumm <- data_summary(m, varname="LSUM", 
                       groupnames=c("Location", "Sex"))
#msumm$se <- msumm$sd/sqrt(3)
# Convert dose to a factor variable
#msumm$Treatment=as.factor(msumm$Treatment)
head(msumm)

p<- ggplot(msumm, aes(x=Sex, y=LSUM, fill=Location)) + 
  geom_bar(stat="identity",  
           position=position_dodge()) +
  geom_errorbar(aes(ymin=LSUM-se, ymax=LSUM+se), width=.2,
                position=position_dodge(.9)) 
print(p)

# statistics on LSUM by male/female * vaccinated/sham
shapiro.test(t2$LSUM) # testing for need for transformation
hist(t2$LSUM)
t2$transLSUM <- log(t2$LSUM+1)
shapiro.test(t2$transLSUM)
hist(t2$transLSUM)
# a priori prediction of males<females and alpha = 0.1
int <- anova(lm(t2$transLSUM ~ t2$Treatment * t2$Sex))
plot(lm(t2$transLSUM ~ t2$Treatment * t2$Sex))
int
noint <- anova(lm(t2$transLSUM ~ t2$Treatment + t2$Sex))
plot(lm(t2$transLSUM ~ t2$Treatment + t2$Sex))
noint





###NOW BY SEX AND TREATMENT
t3summ2 <- data_summary(t3, varname="LSUM", 
                        groupnames=c("Sex", "Treatment"))
t3summ2$se <- t3summ2$sd/sqrt(3)
# Convert dose to a factor variable
t3summ2$Treatment=as.factor(t3summ2$Treatment)
head(t3summ2)

z<- ggplot(t3summ2, aes(x=Sex, y=LSUM, fill=Treatment)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=LSUM-se, ymax=LSUM+se), width=.2,
                position=position_dodge(.9)) 
print(z)

######################################################################################3
# ATTEMPTING TO MERGE DATA FROM ANIMALS BROUGHT INTO THE LAB WITH DATA ON TICK BURDENS
# The code works but there were only two animals with lab results from July 2015 that
# were subsequently captured in Aug/Sept 2015.
######################################################################################
xeno <- read.csv("Xenodiagnosis Jul2015.csv")
dim(xeno)
head(xeno)
xeno <- xeno[,-14]
xeno <- xeno[,-13]
colnames(xeno)[1]<- "Tag" #rename "mouse.Individual" to "Tag" for compatability in merge

tx <- subset(t,Period==2)
write.csv(tx, "output.csv") # send to Excel to try to find the doubled tag record
tx <- tx[-374,]
tx
tx <- subset(tx,Tag != "E7625") # delete the ambiguous record with two captures on diff grids


all <- merge(tx, xeno, by=c("Tag"))
head(all)

######################################################################################3
# PART TWO. 2014 TICK BURDENS WHEN ANIMALS WERE TREATED AT INDIVIDUAL LEVEL NOT GRID LEVEL
#
######################################################################################
b <- read.csv("2014MMTrapping.csv")
dim(b)
head(b)
# replace "NA" with 0 throughout the dataframe
b[is.na(b)]<- 0

# calculate total tick burden for nymphs and adults
b$LSUM <- b$RE+b$LE+b$HD
b$NSUM <- b$RE.1+b$LE.1+b$HD.1

# pull out just those records of P. leucopus and only on experimental grids (Treatment = X)
b <- subset(b,Spp=="PL" & Treatment == "E")
b$Date <- NULL

b$NewDate <- as.character(b$Date.1,format="%y/%m/%d") # convert imported date into characters
b$Date <- as.Date(b$NewDate, format="%Y/%m/%d") # convert characters into date recognized by R
b$NewDate <- NULL
b$Date.1 <- NULL
head(b)

# here make a subset of the data just for larval season from end of July to end of September
#install.packages("dplyr")
library("dplyr")
bs <- filter(b, between(Date, as.Date("2014-07-25"), as.Date("2014-09-30")))
colnames(bs)[10]<- c("Tag")
bs$LSUM <- bs$RE+bs$LE+bs$HD
bs$NSUM <- bs$RE.1+bs$LE.1+bs$HD.1

# read in the three files of vaccination histories of mice in 2014
TX <- read.csv("2014VaccTX.csv")
HX <- read.csv("2014VaccHX.csv")
GX <- read.csv("2014VaccGX.csv")
# bind the three dataframes together
new <- rbind(TX, HX)
newer <- rbind(new,GX)
summary(newer)

# get the dates sorted out
newer$NewDate <- as.character(newer$Date,format="%y/%m/%d") # convert imported date into characters
newer$date <- as.Date(newer$NewDate, format="%Y/%m/%d") # convert characters into date recognized by R
newer$NewDate <- NULL
newer$Date <- NULL
head(newer)

# aggregate by VaccDose to get total number of doses for each animal
vacc <- aggregate(newer[11],by=list(newer$Tag, newer$Location, newer$VaccType, newer$Sex), FUN=sum, ra.nm=TRUE)
vacc[is.na(vacc)]<- 0
colnames(vacc)[1:4]<- c("Tag", "Location", "VaccType", "Sex")

tot <- merge(bs,vacc, by="Tag")
tot[,31:32] <- NULL
tot[,32] <- NULL
# determine the week of highest tick burden, which turns out to be 8-11-14, by finding the highest
# mean tick burden. 
f <- aggregate(tot[30], by=list(tot$Date, tot$Location.x, tot$Sex.x), FUN=mean, ra.nm=TRUE)
colnames(f)[1:3]<- c("Date", "Location", "Sex")
f

# now subset the full dataset for just those two dates of peak larvae
f1 <- subset(tot,Date=="2014-08-11" | Date=="2014-08-18")
f1 <- f1[,c(1,2,3,5:29,30,4,31,32)]
# and aggregate by Location and Sex and VaccType
f2 <- aggregate(f1[29], by=list(f1$Tag, f1$Sex.x, f1$VaccType, f1$Location.x), FUN=mean, ra.nm=TRUE)
colnames(f2)[1:4] <- c("Tag", "Sex", "VaccType", "Location")

shapiro.test(f2$LSUM)
hist(f2$LSUM, breaks=20)
f2$transLSUM <- sqrt(f2$LSUM)
shapiro.test(f2$transLSUM)
hist(f2$transLSUM, breaks=20)

int <- anova(lm(f2$transLSUM ~ f2$Sex * f2$VaccType))
noint <- anova(lm(f2$transLSUM ~ f2$Sex + f2$VaccType))

f2m <- subset(f2,Sex=="1")
blockm <- anova(lm(f2m$transLSUM ~ f2m$VaccType + f2m$Location))
f2f <- subset(f2,Sex=="2")
blockf <- anova(lm(f2f$transLSUM ~ f2f$VaccType + f2f$Location))

# GRAPH
mean.Man <- tapply(f2$LSUM, list(f2$VaccType, f2$Sex), mean)
# Calculate standard deviation of each group using tapply
# which returns a matrix--it will be the same dimension as mean.burde 
sd.Man <- tapply(f2$LSUM, list(f2$VaccType, f2$Sex), sd)
# Calculate sample size in each group
n.Man <- tapply(f2$LSUM, list(f2$VaccType, f2$Sex), length)
# Calculate the standard error 
se.Man <- sd.Man/sqrt(n.Man)
par(mfrow = c(1,2))
# Now make the barplot
# beside=TRUE makes it a grouped barplot
mids<- barplot(mean.Man,
               xlab = "Sex", 
               ylab = "Mean larval burden per mouse at August peak",
              ylim = c(0,25),
               # col=colors,
               beside=TRUE, legend=TRUE)
# use arrows to put error bars on the plot
arrows(mids, mean.Man - se.Man, mids, mean.Man+se.Man, code = 3, angle = 90, length = 0.1)
text(mids, .5, paste ("n = ", n.Man), cex=.7, col="white")


##NOW  JUST THOSE THAT GOT THREE OR MORE DOSES
f3 <- subset(f1,f1$VaccDose > 2)
# and aggregate by Location and Sex and VaccType
f4 <- aggregate(f3[29], by=list(f3$Tag, f3$Sex.x, f3$VaccType, f3$Location.x), FUN=mean, ra.nm=TRUE)
colnames(f4)[1:4] <- c("Tag", "Sex", "VaccType", "Location")

shapiro.test(f4$LSUM)
hist(f4$LSUM, breaks=20)
f4$transLSUM <- log(f4$LSUM+1)
shapiro.test(f4$transLSUM)
hist(f4$transLSUM, breaks=20)

int <- anova(lm(f4$transLSUM ~ f4$Sex * f4$VaccType))
noint <- anova(lm(f4$transLSUM ~ f4$Sex + f4$VaccType))

f4m <- subset(f4,Sex=="1")
blockm <- anova(lm(f4m$transLSUM ~ f4m$VaccType + f4m$Location))
f4f <- subset(f4,Sex=="2")
blockf <- anova(lm(f2f$transLSUM ~ f2f$VaccType + f2f$Location))

# GRAPH
mean.Man <- tapply(f4$LSUM, list(f4$VaccType, f4$Sex), mean)
# Calculate standard deviation of each group using tapply
# which returns a matrix--it will be the same dimension as mean.burde 
sd.Man <- tapply(f4$LSUM, list(f4$VaccType, f2$Sex), sd)
# Calculate sample size in each group
n.Man <- tapply(f4$LSUM, list(f4$VaccType, f4$Sex), length)
# Calculate the standard error 
se.Man <- sd.Man/sqrt(n.Man)
par(mfrow = c(1,2))
# Now make the barplot
# beside=TRUE makes it a grouped barplot
mids<- barplot(mean.Man,
               xlab = "Sex", 
               ylab = "Mean larval burden per mouse, 3+ doses",
               ylim = c(0,25),
               # col=colors,
               beside=TRUE, legend=TRUE)
# use arrows to put error bars on the plot
arrows(mids, mean.Man - se.Man, mids, mean.Man+se.Man, code = 3, angle = 90, length = 0.1)
text(mids, .5, paste ("n = ", n.Man), cex=.7, col="white")

####################################################################
# PART TWO: GUD TRAYS
####################################################################

g<-read.csv("GUD data total 20161101.csv")
head(g)
dim(g)



# this section aggregates by treatment, date, grid, but without using the paired design of shaded/unshaded at each trap location
g1 <- subset(g,Problem=="N",select=c(1:14)) # get rid of any data with problems (e.g. tray flipped over)
g2 <- subset(g1, Year=="2015", select=c(1:14)) # select just data from 2015

# this aggregates by date strips down to the key column seeds remaining
g3 <- aggregate(g2[,11], by=list(g2$Grid, g2$Site, g2$Treatment, g2$Cover, g2$Trap), FUN=mean, ra.nm=FALSE)
colnames(g3)[1:6]<- c("Grid", "Site", "Treatment", "Cover", "Trap", "SeedRemaining")

# graph with each trap location as a 




# aggregate by date
g4 <- aggregate(g3[,6], by=list(g3$Site, g3$Treatment, g3$Cover), FUN=mean, ra.nm=FALSE)
colnames(g4)[1:4]<- c("Site", "Treatment", "Cover", "SeedRemaining")
head(g4)

g4summ <- data_summary(g4, varname="SeedRemaining", 
                       groupnames=c("Treatment", "Cover"))
g4summ$se <- g4summ$sd/sqrt(3)
# Convert dose to a factor variable
g4summ$Treatment=as.factor(g4summ$Treatment)
head(g4summ)

ggplot(g4summ, aes(x=Treatment, y=SeedRemaining, fill=Cover)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=SeedRemaining-se, ymax=SeedRemaining+se), width=.2,
                position=position_dodge(.9))

##PAIRED
# this section takes advantage of the paired shaded-unshaded trays at each location
gS <- subset(g3,Cover=="S",select=c(1:6))
gU <- subset(g3,Cover=="U",select=c(1:6))
gPaired <- merge(gS, gU, by=c("Grid", "Trap", "Site", "Treatment"))
gPaired[,7] <- NULL # remove categorical coding as "S" for shaded
gPaired[,5] <- NULL # remove categorical coding as "U" for uncovered
colnames(gPaired)[5:6] <- c("SeedRemS", "SeedRemU")
gPaired$SeedDiff <- gPaired$SeedRemU - gPaired$SeedRemS # how much more there was in uncovered tray


block <- aov(lm(gPaired$SeedRemS ~ gPaired$Treatment + gPaired$Site))

gPSumm <- data_summary(gPaired, varname=c("SeedDiff"), 
                        groupnames=c("Treatment", "Site"))
gPSumm$se <- gPSumm$sd/sqrt(3)
# Convert dose to a factor variable
gPSumm$Treatment=as.factor(gPSumm$Treatment)
head(gPSumm)

p <- ggplot(gPSumm, aes(x=Treatment, y=SeedDiff, fill = "Site")) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=SeedDiff-se, ymax=SeedDiff+se), width=.2,
                position=position_dodge(.9))




gP3 <- aggregate(gPaired[,6:8], by=list(gPaired$Site, gPaired$Treatment, gPaired$Grid, gPaired$Trap), FUN=mean, ra.nm=FALSE)
colnames(gP3)[1:6]<- c("Site", "Treatment", "Grid", "Trap", "SeedRemS", "SeedRemU")

gP4 <- aggregate(gP3[,5:7], by=list(gP3$Site, gP3$Treatment, gP3$Grid), FUN=mean, ra.nm=FALSE)
colnames(gP4)[1:5]<- c("Site", "Treatment", "Grid", "SeedRemS", "SeedRemU")

shapiro.test(gP3$SeedDiff) #test for normality
hist(gP3$SeedDiff) #plot histogram

# taking square root makes the distribution normal
gP3$transSeedDiff <- sqrt(gP3$SeedDiff)
shapiro.test(gP3$transSeedDiff)



# GRAPH AND STATS OF EFFECTS OF TREATMENT AND SITE ON *DIFFERENCE* BETWEEN COVERED & UNCOVERED TRAYS
int <- aov(gP3$transSeedDiff~gP3$Treatment * gP3$Site)
summary(int)
noint <- aov(gP3$transSeedDiff~gP3$Treatment + gP3$Site)
summary(noint)

gP3Summ <- data_summary(gP3, varname=c("SeedDiff"), 
                       groupnames=c("Treatment"))
gP3Summ$se <- gP3Summ$sd/sqrt(3)
# Convert dose to a factor variable
gPSumm3$Treatment=as.factor(gPSumm3$Treatment)
head(gP3Summ)

p <- ggplot(gP3Summ, aes(x=Treatment, y=SeedDiff)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=SeedDiff-se, ymax=SeedDiff+se), width=.2,
                position=position_dodge(.9))






# POWER ANALYSIS AT GRID LEVEL

library("pwr")
pwr.anova.test(k=2, n=3, f=0.5, sig.level=0.1, power=NULL)
pwr.t.test(n=3, d=0.5, sig.level=0.1, power=NULL)




library(reshape2)


p <- ggplot(gPSumm, aes(x=Site, y=SeedDiff, fill=Treatment)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=SeedDiff-se, ymax=SeedDiff+se), width=.2,
                position=position_dodge(.9))
p+ylab("Mean difference in weight of seeds remaining (g)")+
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=14,face="bold"))+
  theme(legend.title=element_blank())+

  
  
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


##NOW FOR COUNT OF DISTURBED GUD TRAYS
gd <- aggregate(g2[,10], by=list(g2$Grid, g2$Site, g2$Treatment, g2$Cover, g2$Trap, g2$Date), FUN=sum, ra.nm=FALSE)
colnames(gd)[1:7]<- c("Grid", "Site", "Treatment", "Cover", "Trap", "Date", "CountDist")

#aggregate across dates for each trap location, so that TotalDist is total number of visits to that tray 
gd4 <- aggregate(gd[,7], by=list(gd$Grid, gd$Site, gd$Treatment, gd$Cover, gd$Trap), FUN=sum, ra.nm=FALSE)
colnames(gd4)[1:6]<- c("Grid", "Site", "Treatment", "Cover", "Trap", "TotalDist")



##Statistical tests of disturbed trays done two ways
# this way treating each tray pair as a replicate, not blocked by site
gd4$sqrtDist <- sqrt(gd4$TotalDist)
shapiro.test(gd4$TotalDist)
shapiro.test(gd4$sqrtDist)
int <- anova(lm(gd4$TotalDist~gd4$Treatment * gd4$Cover))
int
noint <- anova(lm(gd4$TotalDist~gd4$Treatment + gd4$Cover))


gd4summ <- data_summary(gd4, varname="TotalDist", 
                        groupnames=c("Treatment", "Cover"))
gd4summ$se <- gd4summ$sd/sqrt(3)
# Convert dose to a factor variable
gd4summ$Treatment=as.factor(gd4summ$Treatment)
head(gd4summ)

ggplot(gd4summ, aes(x=Treatment, y=TotalDist, fill=Cover)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=TotalDist-se, ymax=TotalDist+se), width=.2,
                position=position_dodge(.9))


# this way treating each grid as a replicate
shapiro.test(gd5$MeanDist)
int <- aov(gd5$MeanDist~gd5$Treatment * gd5$Cover)
noint <- aov(gd5$MeanDist~gd5$Treatment + gd5$Cover)

gd5summ <- data_summary(gd5, varname="Mean # Trays Disturbed", 
                       groupnames=c("Treatment", "Cover"))
gd5summ$se <- gd5summ$sd/sqrt(3)
# Convert dose to a factor variable
gd5summ$Treatment=as.factor(gd5summ$Treatment)
head(gd5summ)

ggplot(gd5summ, aes(x=Treatment, y=MeanDist, fill=Cover)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=MeanDist-se, ymax=MeanDist+se), width=.2,
                position=position_dodge(.9))

#Mike
kruskal.test(t2$LSUM~as.numeric(t2$Sex))
kruskal.test(t2$LSUM~t2$Treatment)
glm1<-glm(t2$LSUM~as.numeric(t2$Sex), poisson(link=log))
summary(glm1)
glm2<-glm(t2$LSUM~as.numeric(t2$Sex)+t2$Treatment, family=quasipoisson(link=log))
glm3<-glm(t2$LSUM~as.numeric(t2$Sex)+t2$Treatment + s.numeric(t2$Sex)*t2$Treatment, family=quasipoisson(link=log))

library(MASS)
LSUM<-t2$LSUM
Sex<-as.numeric(t2$Sex)
Treatment<-t2$Treatment
kruskal.test(LSUM ~ Treatment)
NB0<-glm.nb(LSUM ~ Treatment)
NB1<-glm.nb(LSUM ~ Sex + Treatment)
NB2<-glm.nb(LSUM ~ Sex + Treatment + Sex*Treatment)

mean(LSUM[Treatment=="X"])  #6.568
mean(LSUM[Treatment=="C"])  #8.231959




