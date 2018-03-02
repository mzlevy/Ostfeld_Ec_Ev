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

rm(list=ls())
library("ggplot2")

# Reset to desired working directory
setwd("~/Ostfeld_Ec_Ev/Tick burden data files and R script")





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
#write.csv(t2,"2015TickSumm.csv")

  
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

# aggregate by date
g4 <- aggregate(g3[,6], by=list(g3$Site, g3$Treatment, g3$Cover), FUN=mean, ra.nm=FALSE)
colnames(g4)[1:4]<- c("Site", "Treatment", "Cover", "SeedRemaining")
head(g4)


kruskal.test(t2$LSUM~as.numeric(t2$Sex))
kruskal.test(t2$LSUM~t2$Treatment)
glm1<-glm(t2$LSUM~as.numeric(t2$Sex), poisson(link=log))
summary(glm1)
glm2<-glm(t2$LSUM~as.numeric(t2$Sex)+t2$Treatment, family=quasipoisson(link=log))

library(MASS)
LSUM<-t2$LSUM
Sex<-as.numeric(t2$Sex)
Treatment<-t2$Treatment
kruskal.test(LSUM ~ Treatment)

#Regressions on tick burden
NB0<-glm.nb(LSUM ~ Treatment)
NB1<-glm.nb(LSUM ~ Sex + Treatment)
NB2<-glm.nb(LSUM ~ Sex + Treatment + Sex*Treatment)

mean(LSUM[Treatment=="X"])  #6.568
mean(LSUM[Treatment=="C"])  #8.231959




