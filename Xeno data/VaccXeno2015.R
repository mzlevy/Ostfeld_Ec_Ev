# ------------------------------------------------------------------------------
# Borrelia Vaccine Project -- Lab Xenodiagnosis script
# 6 February 2017
# Modified: 
# Input: 2015_all_inds_xenodiagnosisdata_year2.csv
#        2014_xenodiagnosisdata_year1.csv
#        2015TickSumm.csv
#        2015VaccSum.csv
# Output: 
#
# ------------------------------------------------------------------------------
# Reset R's brain
rm(list=ls())

# Find where R is currently working
getwd()

# Reset to desired working directory
#setwd("c:/Users/Felicia Keesing/Documents/Data & Figures/Vaccination")
setwd("~/Vaccine/Xeno data")
# Find where R is currently working
getwd()



install.packages("gridExtra")
install.packages("cowplot")
install.packages("dynlm")
library("ggplot2")
library("gridExtra")
library("cowplot")

####################################################################
# XENO 2015
####################################################################

x<-read.csv("2015_all_inds_xenodiagnosisdata_year2.csv")
head(x)
dim(x)

x1 <- aggregate(x[12],by=list(x$site, x$vac_type), FUN=mean)
colnames(x1)[1:2]<- c("Grid", "VaccType")

# split "Grid" into two variables: "Location" and "Treatment"
x1[,4] <- data.frame(Location =   substr(x1$Grid, start = 1, stop = 1))
x1[,5] <- data.frame(Treatment = substr(x1$Grid, start = 2, stop = 2) )
head(x1)

plot(x$vac_type, x$overall_perc_inf)

# CALCULATE # and % INFECTED/UNINFECTED OVERALL AND FOR EACH VACCINATION CLASS
# First, find # of mice that infected a tick
TotNumUninf <- sum(x$overall_num_pos == 0) # total number of mice that produced no infected ticks
TotNumInf <- sum(x$overall_num_pos > 0) # total # of mice that produced at least one infected tick
xSham <- subset(x, x$vac_type=="sham") #subset just sham mice
xTrue <- subset(x, x$vac_type=="TRUE") #subset just TRUE mice
xShamUninf <- sum(xSham$overall_num_pos == 0) # total number of sham mice that produced no infected ticks
xSShamInf <- sum(xSham$overall_num_pos > 0) # total number of sham mice that produced infected ticks
xTrueUninf <- sum(xTrue$overall_num_pos == 0) # total number of true mice that produced no infected ticks
xTrueInf <- sum(xTrue$overall_num_pos > 0) # total number of true mice that produced infected ticks

TotNumInf

#mike added:
#chisquare test for whether they infect at least 1 tick
matrix(c(23,37, 32,35), nrow=2, ncol=2)
twobytwo<-matrix(c(xShamUninf,xSShamInf, xTrueUninf, xTrueInf), nrow=2, ncol=2)
chisq.test(twobytwo)

#poisson model with total processed as an offset
glm1<-glm(formula = overall_num_pos ~ vac_type, family = poisson(link = log),
data = x, offset = overall_num_ticks_proc)

summary(glm1)



# CALCULATE % INFECTED/UNINFECTED and RES COMP FOR THOSE INFECTED
mean.Man <- tapply(x$overall_perc_inf, list(x1$VaccType, x1$Location), mean)
# Calculate standard deviation of each group using tapply
# which returns a matrix--it will be the same dimension as mean.burde 
sd.Man <- tapply(x1$overall_perc_inf, list(x1$VaccType, x1$Location), sd)
# Calculate sample size in each group
n.Man <- tapply(x1$overall_perc_inf, list(x1$VaccType, x1$Location), length)
# Calculate the standard error 
se.Man <- sd.Man/sqrt(n.Man)


# PLOT FOR EACH GRID
mean.Man <- tapply(x1$overall_perc_inf, list(x1$VaccType, x1$Location), mean)
# Calculate standard deviation of each group using tapply
# which returns a matrix--it will be the same dimension as mean.burde 
sd.Man <- tapply(x1$overall_perc_inf, list(x1$VaccType, x1$Location), sd)
# Calculate sample size in each group
n.Man <- tapply(x1$overall_perc_inf, list(x1$VaccType, x1$Location), length)
# Calculate the standard error 
se.Man <- sd.Man/sqrt(n.Man)
par(mfrow = c(1,2))
# Now make the barplot
# beside=TRUE makes it a grouped barplot
mids<- barplot(mean.Man,
               xlab = "Location", 
               ylab = "Mean percentage of ticks infected",
               ylim = c(0,70),
               # col=colors,
               beside=TRUE, legend=TRUE)
# use arrows to put error bars on the plot
arrows(mids, mean.Man - se.Man, mids, mean.Man+se.Man, code = 3, angle = 90, length = 0.1)
text(mids, 5, paste ("n = ", n.Man), cex=.7, col="white")


# PLOT FOR AVERAGE OF GRIDS
mean.Man <- tapply(x1$overall_perc_inf, list(x1$VaccType), mean)
# Calculate standard deviation of each group using tapply
# which returns a matrix--it will be the same dimension as mean.burde 
sd.Man <- tapply(x1$overall_perc_inf, list(x1$VaccType), sd)
# Calculate sample size in each group
n.Man <- tapply(x1$overall_perc_inf, list(x1$VaccType), length)
# Calculate the standard error 
se.Man <- sd.Man/sqrt(n.Man)
par(mfrow = c(1,2))
# Now make the barplot
# beside=TRUE makes it a grouped barplot
mids<- barplot(mean.Man,
               xlab = "Treatment", 
               ylab = "Mean percentage of ticks infected",
               ylim = c(0,50))
# use arrows to put error bars on the plot
arrows(mids, mean.Man - se.Man, mids, mean.Man+se.Man, code = 3, angle = 90, length = 0.1)
text(mids, 5, paste ("n = ", n.Man), cex=.7, col="white")



# VIOLIN PLOT
par(mfrow = c(1,2))
p <- ggplot(x, aes(x=vac_type, y=overall_perc_inf)) +
  geom_violin(trim=FALSE)+labs(title="Xenodiagnosis 2015")
p
# violin plot with mean points
p + stat_summary(fun.y=mean, geom="point", shape=23, size=2)

# violin plot with dot plot
p + geom_dotplot(binaxis='y', stackdir='center', dotsize=1)


# BOX PLOT
# The above adds a redundant legend. With the legend removed:
a <- ggplot(x, aes(x=vac_type, y=overall_perc_inf, fill=vac_type)) + geom_boxplot() +
  guides(fill=FALSE)+
  stat_summary(fun.y=mean, geom="point", shape=19, size=2)+
  labs(x="Vaccination type",y="Overall percentage of ticks infected")+
  theme_bw()
  
plot_grid(a, ncol=2, nrow =1) # puts into one column


####################################################################
# XENO 2014
####################################################################

y<-read.csv("2014_xenodiagnosisdata_year1.csv")
head(y)
dim(y)

#y1 <- aggregate(y[12],by=list(y$site, y$vac_type), FUN=mean)
#colnames(y1)[1:2]<- c("Grid", "VaccType")

# split "Grid" into two variables: "Location" and "Treatment"
y[,14] <- data.frame(Location =   substr(y$site, start = 1, stop = 1))
y[,15] <- data.frame(Treatment = substr(y$site, start = 2, stop = 2) )
head(y)


# PLOT FOR EACH GRID
mean.Man <- tapply(y$overall_perc_inf, list(y$vac_type, y$Location), mean)
# Calculate standard deviation of each group using tapply
# which returns a matrix--it will be the same dimension as mean.burde 
sd.Man <- tapply(y$overall_perc_inf, list(y$vac_type, y$Location), sd)
# Calculate sample size in each group
n.Man <- tapply(y$overall_perc_inf, list(y$vac_type, y$Location), length)
# Calculate the standard error 
se.Man <- sd.Man/sqrt(n.Man)
par(mfrow = c(1,1))
# Now make the barplot
# beside=TRUE makes it a grouped barplot
mids<- barplot(mean.Man,
               xlab = "Location", 
               ylab = "Mean percentage of ticks infected",
               ylim = c(0,100),
               # col=colors,
               beside=TRUE, legend=TRUE)
# use arrows to put error bars on the plot
arrows(mids, mean.Man - se.Man, mids, mean.Man+se.Man, code = 3, angle = 90, length = 0.1)
text(mids, 5, paste ("n = ", n.Man), cex=.7, col="white")




###################################################################################
#
# PLOT INDIVIDUAL MOUSE XENO VS TICK BURDEN
#
###################################################################################
t <- read.csv("2015TickSumm.csv")
colnames(x)[1]<- c("Tag")

xt <- merge(t,x, by="Tag")

# subset males and females
xtm <- subset(xt, xt$Sex == "1" & xt$Treatment=="C")
xtf <- subset(xt, xt$Sex == "2")

plot(xtm$LSUM ~ xtm$overall_perc_inf)
abline(lm(xtm$LSUM ~ xtm$overall_perc_inf))

plot(xtf$LSUM ~ xtf$overall_perc_inf)
abline(lm(xtf$LSUM ~ xtf$overall_perc_inf))

# import number of doses received by each mouse
all <- read.csv("2015VaccSum.csv")
all$date <- as.Date(all$date, format="%m/%d/%Y") # convert characters into date recognized by R

vacc <- aggregate(all[12],by=list(all$Tag, all$Site, all$Treatment, all$VaccType, all$Sex), FUN=sum)
colnames(vacc)[1:6]<- c("Tag", "Site", "Treatment", "VaccType", "Sex", "Dose")
vacc$VaccType <- sub("^$", "None", vacc$VaccType) # replace blanks with "None" for VaccType

# merge to compare how # of doses affects xeno
tot <- merge(x, vacc, by="Tag")

# PLOT 
mean.Man <- tapply(tot$overall_perc_inf, list(tot$vac_type, tot$Dose), mean)
# Calculate standard deviation of each group using tapply
# which returns a matrix--it will be the same dimension as mean.burde 
sd.Man <- tapply(tot$overall_perc_inf, list(tot$vac_type, tot$Dose), sd)
# Calculate sample size in each group
n.Man <- tapply(tot$overall_perc_inf, list(tot$vac_type, tot$Dose), length)
# Calculate the standard error 
se.Man <- sd.Man/sqrt(n.Man)
par(mfrow = c(1,1))
# Now make the barplot
# beside=TRUE makes it a grouped barplot
mids<- barplot(mean.Man,
               xlab = "Doses", 
               ylab = "Mean percentage of ticks infected",
               ylim = c(0,100),
               # col=colors,
               beside=TRUE, legend=TRUE)
# use arrows to put error bars on the plot
arrows(mids, mean.Man - se.Man, mids, mean.Man+se.Man, code = 3, angle = 90, length = 0.1)
text(mids, 5, paste ("n = ", n.Man), cex=.7, col="white")

head(tot)

