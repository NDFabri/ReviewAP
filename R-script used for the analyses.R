rm(list=ls()) #Empty the global environment

##Installing and loading libraries ####

# Install and load package tidyverse
if (!requireNamespace("tidyverse"))
  install.packages("tidyverse")
library("tidyverse") #Use of version 1.2.1

# Install and load package ggpubr
if (!requireNamespace("ggpubr"))
  install.packages("ggpubr")
library("ggpubr") #Use of version 0.4.0

# Install and load package grid
if (!requireNamespace("grid"))
  install.packages("grid")
library("grid") #Use of version 3.6.0

# Install and load package gridExtra
if (!requireNamespace("gridExtra"))
  install.packages("gridExtra")
library(gridExtra) #Use of version 2.3

##Loading datafiles from GitHub ####
T1 <- read_csv(url("https://raw.githubusercontent.com/NDFabri/ReviewAP/master/data/Tick burden.csv"))
AH1 <- read_csv(url("https://raw.githubusercontent.com/NDFabri/ReviewAP/master/data/AP in hosts.csv"))
AF1 <- read_csv(url("https://raw.githubusercontent.com/NDFabri/ReviewAP/master/data/AP in ticks.csv"))

##General data preparation ####
#Making column Row numerical
T1$Row <- as.numeric(gsub(",", ".", gsub("\\.", "", T1$Row)))
AH1$Row <- as.numeric(gsub(",", ".", gsub("\\.", "", AH1$Row)))
AF1$Row <- as.numeric(gsub(",", ".", gsub("\\.", "", AF1$Row)))

#####
##Tick burden - Data preparation and data selection ####
#Data preparation especially for the tick burden dataset

T2 <- subset(T1, `Multiple species`!="Yes") #Exclusion of data with multiple species from different genus
T2 <- subset(T2, Birds!="Migratory birds included") #Exclusion of data where migratory birds are included

#Exclusion of host species with less than X individuals checked for ticks and 
#exclusion of data per stage where less than X individuals are checked for that stage
X.T <- 40 #Defining X 
Excl.spec <- function(x, y, data=data, X){
  g <- aggregate(x~y, data=data, FUN=sum)
  data$Tot <- g$x[match(data$`Species scientific`, g$y)]
  data <- subset(data, Tot>X)
  data <- subset(data, select = -c(Tot))
  return(data)
} #Function to exclude host species with less than X individuals
T3 <- Excl.spec(T2$`Hosts studied`, T2$`Species scientific`, data = T2, X.T) #Exclusion of host species with less than X individuals checked for ticks
Excl.stage <- function(x, X){
  h <- subset(T3, x!="NA")
  h2 <- aggregate(`Hosts studied`~`Species scientific`, data=h, FUN=sum)
  T3$Tot <- h2$`Hosts studied`[match(T3$`Species scientific`, h2$`Species scientific`)]
  r <- subset(T3, select = c(Tot))
  r$x <- x
  r$x[r$Tot <= (X-1)] <- NA
  T3 <- subset(T3, select = -c(Tot))
  return(r$x)
} #Function to exclude data with less than X indivduals checked for a stage
T3$`Larvae on hosts` <- Excl.stage(T3$`Larvae on hosts`, X.T) #Exclusion of data from species with less than X individuals checked for larvae
T3$`Nymphs on hosts` <- Excl.stage(T3$`Nymphs on hosts`, X.T) #Exclusion of data from species with less than X individuals checked for nymphs
T3$`Adults on hosts` <- Excl.stage(T3$`Adults on hosts`, X.T) #Exclusion of data from species with less than X individuals checked for adults

T4 <- subset(T3, `Multiple species`!="Half") #Exclusion of data with multiple species within one genus

##Tick burden - Calculation tick burden per stage ####

T5 <- as.data.frame(unique(T4$`Species scientific`)) #Preparing a data frame to fill in the tick burden data
names(T5)[names(T5) == "unique(T4$`Species scientific`)"] <- "Species.scientific"

#Larval tick burden
y <- subset(T4, `Larvae on hosts`!="NA") #Subset including only individuals checked for larvae
x <- aggregate(`Larvae on hosts` ~ `Species scientific`, data = y, FUN = sum, na.action = na.omit)
T5$Larvae.host <- x$`Larvae on hosts`[match(T5$Species.scientific, x$`Species scientific`)] #Number of larvae found per host species
x <- aggregate(`Hosts studied` ~ `Species scientific`, data = y, FUN = sum, na.action = na.omit)
T5$HStud.L <- x$`Hosts studied`[match(T5$Species.scientific, x$`Species scientific`)] #Number of individuals checked for larvae per host species
T5$TB.L <- T5$Larvae.host/T5$HStud.L #Mean larval tick burden per host species

#Nymphal tick burden
y <- subset(T4, `Nymphs on hosts`!="NA") #Subset including only individuals checked for nymphs
x <- aggregate(`Nymphs on hosts` ~ `Species scientific`, data = y, FUN = sum, na.action = na.omit)
T5$Nymphs.host <- x$`Nymphs on hosts`[match(T5$Species.scientific, x$`Species scientific`)] #Number of nymphs found per host species
x <- aggregate(`Hosts studied` ~ `Species scientific`, data = y, FUN = sum, na.action = na.omit)
T5$HStud.N <- x$`Hosts studied`[match(T5$Species.scientific, x$`Species scientific`)] #Number of individuals checked for nymphs per host species
T5$TB.N <- T5$Nymphs.host/T5$HStud.N #Mean nymphal tick burden per host species

#Adult tick burden
y <- subset(T4, `Adults on hosts`!="NA") #Subset including only individuals checked for adults
x <- aggregate(`Adults on hosts` ~ `Species scientific`, data = y, FUN = sum, na.action = na.omit)
T5$Adults.host <- x$`Adults on hosts`[match(T5$Species.scientific, x$`Species scientific`)] #Number of adults found per host species
x <- aggregate(`Hosts studied` ~ `Species scientific`, data = y, FUN = sum, na.action = na.omit)
T5$HStud.A <- x$`Hosts studied`[match(T5$Species.scientific, x$`Species scientific`)] #Number of individuals checked for adults per host species
T5$TB.A <- T5$Adults.host/T5$HStud.A #Mean adult tick burden per host species

##Anaplasma phagocytophilum in hosts - Data preparation and data selection ####
#Data preparation especially for the AP in hosts dataset

AH2 <- subset(AH1, `Multiple species`!="Yes") #Exclusion of data with multiple species from different genus
AH2 <- subset(AH2, Birds!="Migratory birds included") #Exclusion of data where migratory birds are included
AH2 <- subset(AH2, Sampletype!="Faeces") #Exclusion of individuals tested in faeces
AH2 <- subset(AH2, Sampletype!="Pellet") #Exclusion of individuals tested in pellets

#Exclusion of host species with less than X individuals tested for the presence of Anaplasma phagocytophilum
X.AH <- 30 #Defining X
AH3 <- Excl.spec(AH2$`Hosts studied`, AH2$`Species scientific`, data = AH2, X.AH) #Exclusion of host species with less than X individuals tested

AH4 <- subset(AH3, `Multiple species`!="Half") #Exclusion of data with multiple species within one genus

##Anaplasma phagocytophilum in hosts - Calculation infection prevalence ####

AH5 <- as.data.frame(unique(AH4$`Species scientific`)) #Preparing a data frame to fill in the infection prevalence
names(AH5)[names(AH5) == "unique(AH4$`Species scientific`)"] <- "Species.scientific"

x <- aggregate(`Hosts studied` ~ `Species scientific`, data = AH4, FUN = sum, na.action = na.omit)
AH5$Htes <- x$`Hosts studied`[match(AH5$Species.scientific, x$`Species scientific`)] #Number of individuals tested per host species
x <- aggregate(`Host positive` ~ `Species scientific`, data = AH4, FUN = sum, na.action = na.omit)
AH5$Hpos <- x$`Host positive`[match(AH5$Species.scientific, x$`Species scientific`)] #Number of individuals positive per host species
AH5$Hneg <- AH5$Htes-AH5$Hpos #Number of individuals negative per host species
AH5$IP.H <- AH5$Hpos/AH5$Htes #Infection prevalence in hosts per host species

##Anaplasma phagocytophilum in feeding ticks - Data preparation and data selection ####
#Data preparation especially for the AP in ticks dataset

AF2 <- subset(AF1, `Multiple species`!="Yes") #Exclusion of data with multiple species from different genus
AF2 <- subset(AF2, Birds!="Migratory birds included") #Exclusion of data where migratory birds are included

#Exclusion of data where the tested ticks came from less than X animals
X.AF <- 5 #Defining X
AF3 <- subset(AF2, `Tick hosts` > (X.AF-1))

#Exclusion of data where less than Z ticks per stage were tested
Z.AF <- 10 #Defining Z
AF3$`Larvae tested`[AF3$`Larvae tested` <= (Z.AF-1)] <- NA #Exclusion of data on larvae where less than 10 larvae were tested
AF3$`Nymphs tested`[AF3$`Nymphs tested` <= (Z.AF-1)] <- NA #Exclusion of data on nymphs where less than 10 nymphs were tested
AF3$`Adults tested`[AF3$`Adults tested` <= (Z.AF-1)] <- NA #Exclusion of data on adults where less than 10 adults were tested
AF3$`Ticks tested` <- rowSums(AF3[,c("Larvae tested", "Nymphs tested", "Adults tested")], na.rm=TRUE) #Total number of ticks tested per datarow
AF3 <- subset(AF3, `Ticks tested`!=0) #Excluding datarows where due to exclusion of data on tickstages no data is included anymore

AF4 <- subset(AF3, `Multiple species`!="Half") #Exclusion of data with multiple species within one genus

##Anaplasma phagocytophilum in feeding ticks - Calculation infection prevalence per stage ####

AF5 <- as.data.frame(unique(AF4$`Species scientific`)) #Preparing a data frame to fill in the infection prevalence
names(AF5)[names(AF5) == "unique(AF4$`Species scientific`)"] <- "Species.scientific"

#Infection prevalence in feeding larvae
y <- subset(AF4, `Larvae tested`!="NA") #Subset including individuals from which larvae were tested
x <- aggregate(`Larvae tested` ~ `Species scientific`, data = y, FUN = sum, na.action = na.omit)
AF5$Ltes <- x$`Larvae tested`[match(AF5$Species.scientific, x$`Species scientific`)] #Number of feeding larvae tested per host species
x <- aggregate(`Larvae positive` ~ `Species scientific`, data = y, FUN = sum, na.action = na.omit)
AF5$Lpos <- x$`Larvae positive`[match(AF5$Species.scientific, x$`Species scientific`)] #Number of feeding larvae positive per host species
AF5$Lneg <- AF5$Ltes-AF5$Lpos #Number of feeding larvae negative per host species
AF5$IP.L <- AF5$Lpos/AF5$Ltes #Infection prevalence in feeding larvae per host species

#Infection prevalence in feeding nymphs
y <- subset(AF4, `Nymphs tested`!="NA") #Subset including individuals from which nymphs were tested
x <- aggregate(`Nymphs tested` ~ `Species scientific`, data = y, FUN = sum, na.action = na.omit)
AF5$Ntes <- x$`Nymphs tested`[match(AF5$Species.scientific, x$`Species scientific`)] #Number of feeding nymphs tested per host species
x <- aggregate(`Nymphs positive` ~ `Species scientific`, data = y, FUN = sum, na.action = na.omit)
AF5$Npos <- x$`Nymphs positive`[match(AF5$Species.scientific, x$`Species scientific`)] #Number of feeding nymphs positive per host species
AF5$Nneg <- AF5$Ntes-AF5$Npos #Number of feeding nymphs negative per host species
AF5$IP.N <- AF5$Npos/AF5$Ntes #Infection prevalence in feeding nymphs per host species

#Infection prevalence in feeding adults
y <- subset(AF4, `Adults tested`!="NA") #Subset including individuals from which adults were tested
x <- aggregate(`Adults tested` ~ `Species scientific`, data = y, FUN = sum, na.action = na.omit)
AF5$Ates <- x$`Adults tested`[match(AF5$Species.scientific, x$`Species scientific`)] #Number of feeding adults tested per host species
x <- aggregate(`Adults positive` ~ `Species scientific`, data = y, FUN = sum, na.action = na.omit)
AF5$Apos <- x$`Adults positive`[match(AF5$Species.scientific, x$`Species scientific`)] #Number of feeding adults positive per host species
AF5$Aneg <- AF5$Ates-AF5$Apos #Number of feeding adults negative per host species
AF5$IP.A <- AF5$Apos/AF5$Ates #Infection prevalence in feeding adults per host species

##Dataset with all calculations ####
Data <- unique(rbind(T5[1],AH5[1],AF5[1])) #Preparing a dataframe for all calculations
Data <- merge(Data, T5, by="Species.scientific", all = T) #Adding tick burden (and additional data)
Data <- merge(Data, AH5, by="Species.scientific", all = T) #Adding infection prevalence in host (and additional data)
Data <- merge(Data, AF5, by="Species.scientific", all = T) #Adding infection prevalence in feeding ticks (and additional data)

#Dataset with only species where tick burden is known
DataTB <- Data[!is.na(Data$Larvae.host),] #49 species
DataTB <- DataTB[!is.na(DataTB$Nymphs.host),] #48 species
DataTB <- DataTB[!is.na(DataTB$Adults.host),] #48 species

#Dataset with only species where tick burden AND infection prevalence in ticks is known
DataTBAL <- DataTB[!is.na(DataTB$IP.L),] #For 8 species the tick burden AND IP in larvae is known
DataTBAN <- DataTB[!is.na(DataTB$IP.N),] #For 8 species the tick burden AND IP in nymphs is known
DataTBAA <- DataTB[!is.na(DataTB$IP.A),] #For 6 species the tick burden AND IP in adults is known
DataTBAT <- DataTBAL[!is.na(DataTBAL$IP.N),] #For 5 species the tick burden AND IP in larvae AND nymphs is known
DataTBAT <- DataTBAT[!is.na(DataTBAT$IP.A),] #For 1 specie the tick burden AND IP in all stages is known

#Dataset with only species where tick burden AND infection prevalence in hosts is known
DataTBAH <- DataTB[!is.na(DataTB$IP.H),] #18 species

#Dataset with only species where tick burden AND infection prevalence in hosts OR ticks is known
DataTBATAH <- rbind(DataTBAT, DataTBAH) #19 species

#####
##Preparations for associations####
make_ci <- function(pred, data){
  fit <- pred$fit
  lower <- fit - 1.96*pred$se.fit
  upper <- fit + 1.96*pred$se.fit
  return(data.frame(fit, lower, upper, data))
} #Function to get the 95% CI

##Infection prevalence adults ~ Infection prevalence hosts ####
S.IPH <- Data[!is.na(Data$IP.H),] #Dataset of host species where the infection prevalence in hosts is known
S.IPH.IPA <- S.IPH[!is.na(S.IPH$IP.A),] #Dataset of host species where the infection prevalence in hosts and feeding adults is known

m.IPH.IPA <- glm(cbind(Apos, Aneg) ~ IP.H, data = S.IPH.IPA, family = binomial, na.action = na.omit) #Model for association between infection prevalence in feeding adults and hosts
summary(m.IPH.IPA)
anova(m.IPH.IPA, test = "Chisq") #p<0.001 -> there is an association
confint(m.IPH.IPA)

#Graph of association
p.IPH.IPA <- ggplot(S.IPH.IPA, aes(y = IP.A, x = IP.H)) +
  geom_point(shape = 1, size = log10(S.IPH.IPA$Htes)*2, colour = "black") +
  scale_y_continuous(limits = c(0,1)) +
  scale_x_continuous(limits = c(0,1)) +
  xlab("Infection prevalence in hosts") + ylab("Infection prevalence in feeding adults") +
  theme_bw()
IPHA <- seq(min(S.IPH.IPA$IP.H), max(S.IPH.IPA$IP.H), length.out = 100)
pred.p.IPH.IPA <- data.frame(IPHA = IPHA)
names(pred.p.IPH.IPA)[names(pred.p.IPH.IPA) == "IPHA"] <- "IP.H"
predA <- predict(m.IPH.IPA, newdata = pred.p.IPH.IPA, se.fit = TRUE, type = "response")
pred.p.IPH.IPA <- make_ci(predA, pred.p.IPH.IPA)
names(pred.p.IPH.IPA)[names(pred.p.IPH.IPA) == "fit"] <- "IP.A"
p.IPH.IPA <- p.IPH.IPA + geom_line(data = pred.p.IPH.IPA, size = 1, colour = "black") +
  geom_ribbon(data = pred.p.IPH.IPA, aes(ymin = lower, ymax = upper), alpha=0.5,
              linetype = 0, fill = "grey70")
p.IPH.IPA

##Infection prevalence nymphs ~ Infection prevalence hosts ####
S.IPH.IPN <- S.IPH[!is.na(S.IPH$IP.N),] #Dataset of host species where the infection prevalence in hosts and feeding nymphs is known

m.IPH.IPN <- glm(cbind(Npos, Nneg) ~ IP.H, data = S.IPH.IPN, family = binomial, na.action = na.omit) #Model for association between infection prevalence in feeding nymphs and hosts
summary(m.IPH.IPN)
anova(m.IPH.IPN, test = "Chisq") #p<0.001 -> there is an association
confint(m.IPH.IPN)

#Graph of association
p.IPH.IPN <- ggplot(S.IPH.IPN, aes(y = IP.N, x = IP.H)) +
  geom_point(shape = 1, size = log10(S.IPH.IPN$Htes)*2, colour = "black") +
  scale_y_continuous(limits = c(0,1)) +
  scale_x_continuous(limits = c(0,1)) +
  xlab("Infection prevalence in hosts") + ylab("Infection prevalence in feeding nymphs") +
  theme_bw()
IPHN <- seq(min(S.IPH.IPN$IP.H), max(S.IPH.IPN$IP.H), length.out = 100)
pred.p.IPH.IPN <- data.frame(IPHN = IPHN)
names(pred.p.IPH.IPN)[names(pred.p.IPH.IPN) == "IPHN"] <- "IP.H"
predN <- predict(m.IPH.IPN, newdata = pred.p.IPH.IPN, se.fit = TRUE, type = "response")
pred.p.IPH.IPN <- make_ci(predN, pred.p.IPH.IPN)
names(pred.p.IPH.IPN)[names(pred.p.IPH.IPN) == "fit"] <- "IP.N"
p.IPH.IPN <- p.IPH.IPN + geom_line(data = pred.p.IPH.IPN, size = 1, colour = "black") +
  geom_ribbon(data = pred.p.IPH.IPN, aes(ymin = lower, ymax = upper), alpha = 0.5,
              linetype = 0, fill = "grey70")
p.IPH.IPN

##Infection prevalence larve ~ Infection prevalence hosts (no association) ####
S.IPH.IPL <- S.IPH[!is.na(S.IPH$IP.L),] #Dataset of host species where the infection prevalence in hosts and feeding larvae is known

m.IPH.IPL <- glm(cbind(Lpos, Lneg) ~ IP.H, data = S.IPH.IPL, family = binomial, na.action = na.omit) #Model for association between infection prevalence in feeding larvae and hosts
summary(m.IPH.IPL)
anova(m.IPH.IPL, test = "Chisq") #p=0.671 -> there is no association (low number of species with both the infection prevalence in hosts and in feeding larvae known)
confint(m.IPH.IPL)

##Infection prevalence host ~ Adult tick burden ####
S.IPH.TBA <- S.IPH[!is.na(S.IPH$TB.A),] #Dataset of host species where the infection prevalence in hosts and the adult tick burden is known
S.IPH.TBA["Adults.host"][S.IPH.TBA["Adults.host"] == 0] <- 1
S.IPH.TBA$TB.A <- S.IPH.TBA$Adults.host / S.IPH.TBA$HStud.A

m.IPH.TBA <- glm(cbind(Hpos, Hneg) ~ log10(TB.A), data = S.IPH.TBA, family = binomial, na.action = na.omit) #Model for association between infection prevalence in host and adult tick burden
summary(m.IPH.TBA)
anova(m.IPH.TBA, test="Chisq") #p<0.001 -> there is an association
confint(m.IPH.TBA)

#Graph of association
p.IPH.TBA <- ggplot(S.IPH.TBA, aes(y = IP.H, x = TB.A)) +
  geom_point(shape = 1, size = log10(S.IPH.TBA$HStud.A)*2, colour = "black") +
  scale_x_log10(limits = c(0.0001,120), breaks = c(0.0001,0.001,0.01,0.10,1,10,100),
                label = c("0.0001","0.001","0.01","0.1","1", "10", "100")) +
  scale_y_continuous(limits = c(0.,1)) +
  xlab("Mean adult burden") + ylab("Infection prevalence in hosts") +
  theme_bw()
TB.A <- seq(min(S.IPH.TBA$TB.A), max(S.IPH.TBA$TB.A), length.out = 1000)
pred.p.IPH.TBA <- data.frame(TB.A = TB.A)
pred.TBA <- predict(m.IPH.TBA, newdata = pred.p.IPH.TBA, se.fit = TRUE, type = "response")
pred.p.IPH.TBA <- make_ci(pred.TBA, pred.p.IPH.TBA)
names(pred.p.IPH.TBA)[names(pred.p.IPH.TBA) == "fit"] <- "IP.H"
p.IPH.TBA <- p.IPH.TBA + geom_line(data = pred.p.IPH.TBA, size = 1, colour = "black") +
  geom_ribbon(data = pred.p.IPH.TBA, aes(ymin = lower, ymax = upper), alpha = 0.5,
              linetype = 0, fill = "grey70")
p.IPH.TBA

##Infection prevalence host ~ Nymphal tick burden ####
S.IPH.TBN <- S.IPH[!is.na(S.IPH$TB.N),] #Dataset of host species where the infection prevalence in hosts and the nymphal tick burden is known
S.IPH.TBN["Nymphs.host"][S.IPH.TBN["Nymphs.host"] == 0] <- 1
S.IPH.TBN$TB.N <- S.IPH.TBN$Nymphs.host / S.IPH.TBN$HStud.N

m.IPH.TBN <- glm(cbind(Hpos, Hneg) ~ log10(TB.N), data = S.IPH.TBN, family = binomial, na.action = na.omit) #Model for association between infection prevalence in host and nymphal tick burden
summary(m.IPH.TBN)
anova(m.IPH.TBN, test = "Chisq") #p<0.001 -> there is an association
confint(m.IPH.TBN)

#Graph of association
p.IPH.TBN <- ggplot(S.IPH.TBN, aes(y = IP.H, x = TB.N)) +
  geom_point(shape = 1, size = log10(S.IPH.TBN$HStud.N)*2, colour = "black") +
  scale_x_log10(limits = c(0.0001,120), breaks = c(0.0001,0.001,0.01,0.10,1,10,100),
                label = c("0.0001","0.001","0.01","0.1","1", "10", "100")) +
  scale_y_continuous(limits = c(0.,1)) +
  xlab("Mean nymphal burden") + ylab("Infection prevalence in hosts") +
  theme_bw()
TB.N <- seq(min(S.IPH.TBN$TB.N), max(S.IPH.TBN$TB.N), length.out = 1000)
pred.p.IPH.TBN <- data.frame(TB.N = TB.N)
pred.TBN <- predict(m.IPH.TBN, newdata = pred.p.IPH.TBN, se.fit = TRUE, type = "response")
pred.p.IPH.TBN <- make_ci(pred.TBN, pred.p.IPH.TBN)
names(pred.p.IPH.TBN)[names(pred.p.IPH.TBN) == "fit"] <- "IP.H"
p.IPH.TBN <- p.IPH.TBN + geom_line(data = pred.p.IPH.TBN, size = 1, colour = "black") +
  geom_ribbon(data = pred.p.IPH.TBN, aes(ymin = lower, ymax = upper), alpha = 0.5,
              linetype = 0, fill = "grey70")
p.IPH.TBN

##Infection prevalence host ~ Larval tick burden ####
S.IPH.TBL <- S.IPH[!is.na(S.IPH$TB.L),] #Dataset of host specie where the infection prevalence in hosts and the larval tick burden is known
S.IPH.TBL["Larvae.host"][S.IPH.TBL["Larvae.host"] == 0] <- 1
S.IPH.TBL$TB.N <- S.IPH.TBL$Larvae.host / S.IPH.TBL$HStud.L

m.IPH.TBL <- glm(cbind(Hpos, Hneg) ~ log10(TB.L), data = S.IPH.TBL, family = binomial, na.action = na.omit) #Model for association between infection prevalence in host and larval tick burden
summary(m.IPH.TBL)
anova(m.IPH.TBL, test = "Chisq") #p<0.001 -> there is an association
confint(m.IPH.TBL)

#Graph of association
p.IPH.TBL <- ggplot(S.IPH.TBL, aes(y = IP.H, x = TB.L)) +
  geom_point(shape = 1, size = log10(S.IPH.TBL$HStud.L)*2, colour = "black") +
  scale_x_log10(limits = c(0.0001,120), breaks = c(0.0001,0.001,0.01,0.10,1,10,100),
                label = c("0.0001","0.001","0.01","0.1","1", "10", "100")) +
  scale_y_continuous(limits = c(0.,1)) +
  xlab("Mean larval burden") + ylab("Infection prevalence in hosts") +
  theme_bw()
TB.L <- seq(min(S.IPH.TBL$TB.L), max(S.IPH.TBL$TB.L), length.out = 1000)
pred.p.IPH.TBL <- data.frame(TB.L = TB.L)
pred.TBL <- predict(m.IPH.TBL, newdata = pred.p.IPH.TBL, se.fit = TRUE, type = "response")
pred.p.IPH.TBL <- make_ci(pred.TBL, pred.p.IPH.TBL)
names(pred.p.IPH.TBL)[names(pred.p.IPH.TBL) == "fit"] <- "IP.H"
p.IPH.TBL <- p.IPH.TBL + geom_line(data = pred.p.IPH.TBL, size = 1, colour = "black") +
  geom_ribbon(data = pred.p.IPH.TBL, aes(ymin = lower, ymax = upper), alpha = 0.5,
              linetype = 0, fill = "grey70")
p.IPH.TBL

##Prediction IPH based tick burden for species where the IPL is known ####
S.IPH.IPL2 <- Data[!is.na(Data$IP.L),]
S.IPH.IPL2a <- S.IPH.IPL2[is.na(S.IPH.IPL2$IP.H),] #Subset where the infection prevalence in feeding larvae is known, but the infection prevalence in hosts is not

S.IPH.IPL2a$Species.scientific
#Since the species in S.IPH.IPL2a are three birds and a medium sized mammal, we predict the infection prevalence in hosts based on the nymphal burden
S.IPH.IPL2a$IP.H <- predict(m.IPH.TBN, S.IPH.IPL2a, type = "response")

##Infection prevalence larve ~ Infection prevalence hosts (with predicted values) ####
S.IPH.IPL2a$Shape <- 4 #column to indicate the shape for predicted values in the graph
S.IPH.IPL2b <- S.IPH.IPL2[!is.na(S.IPH.IPL2$IP.H),] #Subset where the infection prevalence in feeding larvae and in hosts is known
S.IPH.IPL2b$Shape <- 1 #column to indicate the shape for measured values in the graph
S.IPH.IPL2 <- rbind(S.IPH.IPL2a, S.IPH.IPL2b) #Subset where the infection prevalence in feeding alrvae is known, and the infection prevalence in hosts is either known or predicted

m.IPH.IPL2 <- glm(cbind(Lpos, Lneg) ~ IP.H, data = S.IPH.IPL2, family = binomial, na.action = na.omit)
summary(m.IPH.IPL2)
anova(m.IPH.IPL2, test="Chisq") #p=0.035 -> there is an association
confint(m.IPH.IPL2)

#Graph of association
S.IPH.IPL2$Htes[is.na(S.IPH.IPL2$Htes)] <- 50 #To indicate the size of the predicted values in the graph
p.IPH.IPL <- ggplot(S.IPH.IPL2, aes(y = IP.L, x = IP.H)) +
  geom_point(shape = S.IPH.IPL2$Shape, size = log10(S.IPH.IPL2$Htes)*2, colour = "black") +
  scale_y_continuous(limits = c(0,1)) +
  scale_x_continuous(limits = c(0,1)) +
  xlab("Infection prevalence in hosts") + ylab("Infection prevalence in feeding larvae") +
  theme_bw()
IPHL <- seq(min(S.IPH.IPL2$IP.H), max(S.IPH.IPL2$IP.H), length.out = 100)
pred.p.IPH.IPL <- data.frame(IPHL = IPHL)
names(pred.p.IPH.IPL)[names(pred.p.IPH.IPL) == "IPHL"] <- "IP.H"
predL <- predict(m.IPH.IPL2, newdata = pred.p.IPH.IPL, se.fit = TRUE, type = "response")
pred.p.IPH.IPL <- make_ci(predL, pred.p.IPH.IPL)
names(pred.p.IPH.IPL)[names(pred.p.IPH.IPL) == "fit"] <- "IP.L"
p.IPH.IPL <- p.IPH.IPL + geom_line(data = pred.p.IPH.IPL, size = 1, colour = "black") +
  geom_ribbon(data = pred.p.IPH.IPL, aes(ymin = lower, ymax = upper), alpha=0.5,
              linetype = 0, fill = "grey70")
p.IPH.IPL

##Graphs of all associations ####
Plot.ass <-
  ggarrange(arrangeGrob(p.IPH.TBL + xlab("Larval burden") + theme(axis.title=element_text(size=12)),
                        top = textGrob("A", x=unit(0.02,"npc"), y=unit(0.95,"npc"), just=c("left","top"), gp=gpar(col="black", fontsize=25))),
            arrangeGrob(p.IPH.TBN + xlab("Nymphal burden") + theme(axis.title=element_text(size=12)),
                        top = textGrob("B", x=unit(0.02,"npc"), y=unit(0.95,"npc"), just=c("left","top"), gp=gpar(col="black", fontsize=25))),
            arrangeGrob(p.IPH.TBA + xlab("Adult burden") + theme(axis.title=element_text(size=12)),
                        top = textGrob("C", x=unit(0.02,"npc"), y=unit(0.95,"npc"), just=c("left","top"), gp=gpar(col="black", fontsize=25))),
            arrangeGrob(p.IPH.IPL + ylab("Infection prevalence in feeding larvae") + theme(axis.title=element_text(size=12)),
                        top = textGrob("D", x=unit(0.02,"npc"), y=unit(0.95,"npc"), just=c("left","top"), gp=gpar(col="black", fontsize=25))),
            arrangeGrob(p.IPH.IPN + ylab("Infection prevalence in feeding nymphs") + theme(axis.title=element_text(size=12)),
                        top = textGrob("E", x=unit(0.02,"npc"), y=unit(0.95,"npc"), just=c("left","top"), gp=gpar(col="black", fontsize=25))),
            arrangeGrob(p.IPH.IPA + ylab("Infection prevalence in feeding adults") + theme(axis.title=element_text(size=12)), 
                        top = textGrob("F", x=unit(0.02,"npc"), y=unit(0.95,"npc"), just=c("left","top"), gp=gpar(col="black", fontsize=25))),
            nrow=2, ncol=3) 
Plot.ass

#####
##Compose host communities ####
#Host community 1 (With C. elaphus)
Species.scientific <- c("Alces alces",
                        "Capreolus capreolus",
                        "Cervus elaphus",
                        "Sus scrofa",
                        "Vulpes vulpes",
                        "Apodemus sylvaticus",
                        "Microtus agrestis",
                        "Myodes glareolus",
                        "Sorex araneus",
                        "Turdus merula")
HC1 <- data.frame(Species.scientific = Species.scientific)
HC1 $Host.tax <- c("Ungulate",
                   "Ungulate",
                   "Ungulate",
                   "Ungulate",
                   "Medium sized mammal",
                   "Small mammal",
                   "Small mammal",
                   "Small mammal",
                   "Small mammal",
                   "Medium sized bird")
HC1$DensityL <- c(1, 10, 1, 1, 1, 100, 100, 100, 100, 100) #Density during low phase year
HC1$DensityH <- c(1, 10, 1, 1, 1, 1000, 1000, 1000, 1000, 100) #Density during peak year
HC1$Ecotype1 <- c(0.65, 0.15, 0.94, 0.94, 1.00, 0.07, 0.07, 0.07, 0.07, 0.56) #Proportion of ecotype 1 in AP positive samples (according to Jaarsma et al., 2019)
HC1$Ecotype2 <- c(0.35, 0.85, 0.06, 0.06, 0.00, 0.04, 0.04, 0.04, 0.04, 0.05) #Proportion of ecotype 2 in AP positive samples (according to Jaarsma et al., 2019)
HC1$Ecotype3 <- c(0.00, 0.00, 0.00, 0.00, 0.00, 0.88, 0.88, 0.88, 0.88, 0.00) #Proportion of ecotype 3 in AP positive samples (according to Jaarsma et al., 2019)
HC1$Ecotype4 <- c(0.00, 0.00, 0.00, 0.00, 0.00, 0.02, 0.02, 0.02, 0.02, 0.39) #Proportion of ecotype 4 in AP positive samples (according to Jaarsma et al., 2019)

#Host community 2 (Without C. elaphus)
Species.scientific <- c("Alces alces",
                        "Capreolus capreolus",
                        "Sus scrofa",
                        "Vulpes vulpes",
                        "Apodemus sylvaticus",
                        "Microtus agrestis",
                        "Myodes glareolus",
                        "Sorex araneus",
                        "Turdus merula")
HC2 <- data.frame(Species.scientific = Species.scientific)
HC2 $Host.tax <- c("Ungulate",
                   "Ungulate",
                   "Ungulate",
                   "Medium sized mammal",
                   "Small mammal",
                   "Small mammal",
                   "Small mammal",
                   "Small mammal",
                   "Medium sized bird")
HC2$DensityL <- c(1, 10, 1, 1, 100, 100, 100, 100, 100) #Density during low phase year
HC2$DensityH <- c(1, 10, 1, 1, 1000, 1000, 1000, 1000, 100) #Density during peak year
HC2$Ecotype1 <- c(0.65, 0.15, 0.94, 1.00, 0.07, 0.07, 0.07, 0.07, 0.56) #Proportion of ecotype 1 in AP positive samples (according to Jaarsma et al., 2019)
HC2$Ecotype2 <- c(0.35, 0.85, 0.06, 0.00, 0.04, 0.04, 0.04, 0.04, 0.05) #Proportion of ecotype 2 in AP positive samples (according to Jaarsma et al., 2019)
HC2$Ecotype3 <- c(0.00, 0.00, 0.00, 0.00, 0.88, 0.88, 0.88, 0.88, 0.00) #Proportion of ecotype 3 in AP positive samples (according to Jaarsma et al., 2019)
HC2$Ecotype4 <- c(0.00, 0.00, 0.00, 0.00, 0.02, 0.02, 0.02, 0.02, 0.39) #Proportion of ecotype 4 in AP positive samples (according to Jaarsma et al., 2019)

#####
##Adding tick burden to host community datasets ####

#Adding larval tick burden
HC1$TB.L <- Data$TB.L[match(HC1$Species.scientific, Data$Species.scientific)]
HC2$TB.L <- Data$TB.L[match(HC2$Species.scientific, Data$Species.scientific)]

#Adding nymphal tick burden
HC1$TB.N <- Data$TB.N[match(HC1$Species.scientific, Data$Species.scientific)]
HC2$TB.N <- Data$TB.N[match(HC2$Species.scientific, Data$Species.scientific)]

#Adding adult tick burden
HC1$TB.A <- Data$TB.A[match(HC1$Species.scientific, Data$Species.scientific)]
HC2$TB.A <- Data$TB.A[match(HC2$Species.scientific, Data$Species.scientific)]

##Adding host infection prevalence to host community datasets ####

#Host community 1
HC1$IP.H <- Data$IP.H[match(HC1$Species.scientific, Data$Species.scientific)]
#Host community 2
HC2$IP.H <- Data$IP.H[match(HC2$Species.scientific, Data$Species.scientific)]

##Adding tick infection prevalence to host community datasets ####
#tick infection prevalence is predicted based on infection prevalence

#Larval infection prevalence
HC1$IP.L <- predict(m.IPH.IPL2, HC1, type = "response")
HC2$IP.L <- predict(m.IPH.IPL2, HC2, type = "response")

#Nymphal infection prevalence
HC1$IP.N <- predict(m.IPH.IPN, HC1, type = "response")
HC2$IP.N <- predict(m.IPH.IPN, HC2, type = "response")

#Adult infection prevalence
HC1$IP.A <- predict(m.IPH.IPA, HC1, type = "response")
HC2$IP.A <- predict(m.IPH.IPA, HC2, type = "response")

#####
##Preparations for plotting relative importances ####
plot.Imp <- function(s,z){
  HC$stage <- s
  ggplot(HC, aes(fill=Host.tax, y=stage, x=HC)) + 
    geom_bar(position="stack", stat="identity") +
    ylab(z) +
    xlab("") + theme(legend.position = "bottom")+ theme(legend.title = element_blank())
} #Function to make plots of relative importance
my_colors <- RColorBrewer::brewer.pal(8, "Set1")[2:6] #Colours for the graphs

##Relative importance in producing A.phagocytophilum positive ticks ####
ImpAPL <- function(data, s){
  x <- s
  data$Totaladults <- data$TB.A * data$DensityL #Total number of adults feeding per species
  Totaladults.data <- sum(data$Totaladults) #Total number of adults feeding
  data$Adultsinf1 <- data$Totaladults*data$IP.A*x #Total number of feeding adults infected per species
  Totaladultsinf.data <- sum(data$Adultsinf1) #Total number of feeding adults infected
  data$Totalnymphs <- data$TB.N * data$DensityL #Total number of nymphs feeding per species
  Totalnymphs.data <- sum(data$Totalnymphs) #Total number of nymphs feeding
  data$Nymphsinf1 <- data$Totalnymphs*data$IP.N*x #Total number of feeding nymphs infected per species
  Totalnymphsinf.data <- sum(data$Nymphsinf1) #Total number of feeding nymphs infected
  data$Totallarvae <- data$TB.L * data$DensityL #Total number of larvae feeding per species
  Totallarvae.data <- sum(data$Totallarvae) #Total number of larvae feeding
  data$Larvaeinf1 <- data$Totallarvae*data$IP.L*x #Total number of feeding larvae infected per species
  Totallarvaeinf.data <- sum(data$Larvaeinf1) #Total number of larvae infected
  return(data)
} #Function to calculate relative importance in producing AP ticks per life stage in low phase
ImpAPH <- function(data, s){
  x <- s
  data$Totaladults <- data$TB.A * data$DensityH
  Totaladults.data <- sum(data$Totaladults)
  data$Adultsinf1 <- data$Totaladults*data$IP.A*x
  Totaladultsinf.data <- sum(data$Adultsinf1)
  data$Totalnymphs <- data$TB.N * data$DensityH
  Totalnymphs.data <- sum(data$Totalnymphs)
  data$Nymphsinf1 <- data$Totalnymphs*data$IP.N*x
  Totalnymphsinf.data <- sum(data$Nymphsinf1)
  data$Totallarvae <- data$TB.L * data$DensityH
  Totallarvae.data <- sum(data$Totallarvae)
  data$Larvaeinf1 <- data$Totallarvae*data$IP.L*x
  Totallarvaeinf.data <- sum(data$Larvaeinf1)
  return(data)
} #Function to calculate relative importance in producing AP ticks per life stage in peak phase
max <- 100 #Maximum limit for the y-axis of the graphs of relative importance in infection prevalence
labsize <- 12 #Size of the text of label axes
txtsize <- 8 #Size of the text of axes
grobsize <- 14 #Size of the annotation above the graph

# # Ecotype 1 ####
HC1La1 <- ImpAPL(HC1, HC1$Ecotype1)
HC2La1 <- ImpAPL(HC2, HC2$Ecotype1)
HC1Ha1 <- ImpAPH(HC1, HC1$Ecotype1)
HC2Ha1 <- ImpAPH(HC2, HC2$Ecotype1)

#Making graphs for number of ticks
HC1La1$HC <- "CE+ Low"
HC2La1$HC <- "CE- Low"
HC1Ha1$HC <- "CE+ Peak"
HC2Ha1$HC <- "CE- Peak"
HC <- rbind(HC1La1, HC2La1, HC1Ha1, HC2Ha1) #This should always be called HC, otherwise the function plot.Imp does not work
HC$Host.tax <- factor(HC$Host.tax, levels = c("Small bird", "Medium sized bird", "Small mammal", "Medium sized mammal", "Ungulate")) #To reorder the classes
HC$HC <- factor(HC$HC, levels = c("CE+ Low", "CE+ Peak", "CE- Low", "CE- Peak")) #To reorder the host communities

plot.eco1 <-
  annotate_figure(
    ggarrange(annotate_figure(plot.Imp(HC$Larvaeinf1, "Number of infected feeding ticks") + ylim(0,max) + scale_fill_manual(values = my_colors) +
                                theme(legend.position = "none", axis.title=element_text(size=labsize), axis.text=element_text(size=txtsize)), top = text_grob("Larvae", size = labsize)),
              annotate_figure(plot.Imp(HC$Nymphsinf1, "") + ylim(0,max) + scale_fill_manual(values = my_colors) +
                                theme(legend.position = "none", axis.title=element_text(size=labsize), axis.text=element_text(size=txtsize)), top = text_grob("Nymphs", size = labsize)),
              annotate_figure(plot.Imp(HC$Adultsinf1, "") + ylim(0,max) + scale_fill_manual(values = my_colors) +
                                theme(legend.position = "none", axis.title=element_text(size=labsize), axis.text=element_text(size=txtsize)), top = text_grob("Adults", size = labsize)),
              nrow=1, ncol=3), 
    top = text_grob("Ecotype 1", size = grobsize))
plot.eco1

#Percentage of positive ticks of all feeding ticks
sum(HC1La1$Larvaeinf1) / sum(HC1La1$Totallarvae) * 100 #Percentage of positive larvae in CE+ Low
sum(HC2La1$Larvaeinf1) / sum(HC2La1$Totallarvae) * 100 #Percentage of positive larvae in CE- Low
sum(HC1Ha1$Larvaeinf1) / sum(HC1Ha1$Totallarvae) * 100 #Percentage of positive larvae in CE+ Peak
sum(HC2Ha1$Larvaeinf1) / sum(HC2Ha1$Totallarvae) * 100 #Percentage of positive larvae in CE- Peak

sum(HC1La1$Nymphsinf1) / sum(HC1La1$Totalnymphs) * 100 #Percentage of positive nymphs in CE+ Low
sum(HC2La1$Nymphsinf1) / sum(HC2La1$Totalnymphs) * 100 #Percentage of positive nymphs in CE- Low
sum(HC1Ha1$Nymphsinf1) / sum(HC1Ha1$Totalnymphs) * 100 #Percentage of positive nymphs in CE+ Peak
sum(HC2Ha1$Nymphsinf1) / sum(HC2Ha1$Totalnymphs) * 100 #Percentage of positive nymphs in CE- Peak

sum(HC1La1$Adultsinf1) / sum(HC1La1$Totaladults) * 100 #Percentage of positive adults in CE+ Low
sum(HC2La1$Adultsinf1) / sum(HC2La1$Totaladults) * 100 #Percentage of positive adults in CE- Low
sum(HC1Ha1$Adultsinf1) / sum(HC1Ha1$Totaladults) * 100 #Percentage of positive adults in CE+ Peak
sum(HC2Ha1$Adultsinf1) / sum(HC2Ha1$Totaladults) * 100 #Percentage of positive adults in CE- Peak

#Number of infected ticks per host community
sum(HC1La1$Larvaeinf1) #Number of infected larvae in CE+ Low
sum(HC2La1$Larvaeinf1) #Number of infected larvae in CE- Low
sum(HC1Ha1$Larvaeinf1) #Number of infected larvae in CE+ Peak
sum(HC2Ha1$Larvaeinf1) #Number of infected larvae in CE- Peak

sum(HC1La1$Nymphsinf1) #Number of infected nymphs in CE+ Low
sum(HC2La1$Nymphsinf1) #Number of infected nymphs in CE- Low
sum(HC1Ha1$Nymphsinf1) #Number of infected nymphs in CE+ Peak
sum(HC2Ha1$Nymphsinf1) #Number of infected nymphs in CE- Peak

sum(HC1La1$Adultsinf1) #Number of infected adults in CE+ Low
sum(HC2La1$Adultsinf1) #Number of infected adults in CE- Low
sum(HC1Ha1$Adultsinf1) #Number of infected adults in CE+ Peak
sum(HC2Ha1$Adultsinf1) #Number of infected adults in CE- Peak

#Percentage of taxomomic host species in their relative contribution to the production of AP positive feeding ticks
#Order = Medium sized bird, Medium sized mammal, Small mammal, Ungulate
data.frame(aggregate(Larvaeinf1~Host.tax, HC1La1, FUN = sum))$Larvaeinf1 / sum(HC1La1$Larvaeinf1) * 100 #Percentage of infected larvae per host taxonomic group in CE+ Low
data.frame(aggregate(Larvaeinf1~Host.tax, HC2La1, FUN = sum))$Larvaeinf1 / sum(HC2La1$Larvaeinf1) * 100 #Percentage of infected larvae per host taxonomic group in CE- Low
data.frame(aggregate(Larvaeinf1~Host.tax, HC1Ha1, FUN = sum))$Larvaeinf1 / sum(HC1Ha1$Larvaeinf1) * 100 #Percentage of infected larvae per host taxonomic group in CE+ Peak
data.frame(aggregate(Larvaeinf1~Host.tax, HC2Ha1, FUN = sum))$Larvaeinf1 / sum(HC2Ha1$Larvaeinf1) * 100 #Percentage of infected larvae per host taxonomic group in CE- Peak

data.frame(aggregate(Nymphsinf1~Host.tax, HC1La1, FUN = sum))$Nymphsinf1 / sum(HC1La1$Nymphsinf1) * 100 #Percentage of infected nymphs per host taxonomic group in CE+ Low
data.frame(aggregate(Nymphsinf1~Host.tax, HC2La1, FUN = sum))$Nymphsinf1 / sum(HC2La1$Nymphsinf1) * 100 #Percentage of infected nymphs per host taxonomic group in CE- Low
data.frame(aggregate(Nymphsinf1~Host.tax, HC1Ha1, FUN = sum))$Nymphsinf1 / sum(HC1Ha1$Nymphsinf1) * 100 #Percentage of infected nymphs per host taxonomic group in CE+ Peak
data.frame(aggregate(Nymphsinf1~Host.tax, HC2Ha1, FUN = sum))$Nymphsinf1 / sum(HC2Ha1$Nymphsinf1) * 100 #Percentage of infected nymphs per host taxonomic group in CE- Peak

data.frame(aggregate(Adultsinf1~Host.tax, HC1La1, FUN = sum))$Adultsinf1 / sum(HC1La1$Adultsinf1) * 100 #Percentage of infected adults per host taxonomic group in CE+ Low
data.frame(aggregate(Adultsinf1~Host.tax, HC2La1, FUN = sum))$Adultsinf1 / sum(HC2La1$Adultsinf1) * 100 #Percentage of infected adults per host taxonomic group in CE- Low
data.frame(aggregate(Adultsinf1~Host.tax, HC1Ha1, FUN = sum))$Adultsinf1 / sum(HC1Ha1$Adultsinf1) * 100 #Percentage of infected adults per host taxonomic group in CE+ Peak
data.frame(aggregate(Adultsinf1~Host.tax, HC2Ha1, FUN = sum))$Adultsinf1 / sum(HC2Ha1$Adultsinf1) * 100 #Percentage of infected adults per host taxonomic group in CE- Peak

# # Ecotype 2 ####
HC1La2 <- ImpAPL(HC1, HC1$Ecotype2)
HC2La2 <- ImpAPL(HC2, HC2$Ecotype2)
HC1Ha2 <- ImpAPH(HC1, HC1$Ecotype2)
HC2Ha2 <- ImpAPH(HC2, HC2$Ecotype2)

#Making graphs for number of ticks
HC1La2$HC <- "CE+ Low"
HC2La2$HC <- "CE- Low"
HC1Ha2$HC <- "CE+ Peak"
HC2Ha2$HC <- "CE- Peak"
HC <- rbind(HC1La2, HC2La2, HC1Ha2, HC2Ha2) #This should always be called HC, otherwise the function plot.Imp does not work
HC$Host.tax <- factor(HC$Host.tax, levels = c("Small bird", "Medium sized bird", "Small mammal", "Medium sized mammal", "Ungulate")) #To reorder the classes
HC$HC <- factor(HC$HC, levels = c("CE+ Low","CE+ Peak", "CE- Low", "CE- Peak")) #To reorder the host communities

plot.eco2 <-
  annotate_figure(
    ggarrange(annotate_figure(plot.Imp(HC$Larvaeinf1, "Number of infected feeding ticks") + ylim(0,max) + scale_fill_manual(values = my_colors) +
                                theme(legend.position = "none", axis.title=element_text(size=labsize), axis.text=element_text(size=txtsize)), top = text_grob("Larvae", size = labsize)),
              annotate_figure(plot.Imp(HC$Nymphsinf1, "") + ylim(0,max) + scale_fill_manual(values = my_colors) +
                                theme(legend.position = "none", axis.title=element_text(size=labsize), axis.text=element_text(size=txtsize)), top = text_grob("Nymphs", size = labsize)),
              annotate_figure(plot.Imp(HC$Adultsinf1, "") + ylim(0,max) + scale_fill_manual(values = my_colors) +
                                theme(legend.position = "none", axis.title=element_text(size=labsize), axis.text=element_text(size=txtsize)), top = text_grob("Adults", size = labsize)),
              nrow=1, ncol=3), 
    top = text_grob("Ecotype 2", size = grobsize))
plot.eco2

#Percentage of positive ticks of all feeding ticks
sum(HC1La2$Larvaeinf1) / sum(HC1La2$Totallarvae) * 100 #Percentage of positive larvae in CE+ Low
sum(HC2La2$Larvaeinf1) / sum(HC2La2$Totallarvae) * 100 #Percentage of positive larvae in CE- Low
sum(HC1Ha2$Larvaeinf1) / sum(HC1Ha2$Totallarvae) * 100 #Percentage of positive larvae in CE+ Peak
sum(HC2Ha2$Larvaeinf1) / sum(HC2Ha2$Totallarvae) * 100 #Percentage of positive larvae in CE- Peak

sum(HC1La2$Nymphsinf1) / sum(HC1La2$Totalnymphs) * 100 #Percentage of positive nymphs in CE+ Low
sum(HC2La2$Nymphsinf1) / sum(HC2La2$Totalnymphs) * 100 #Percentage of positive nymphs in CE- Low
sum(HC1Ha2$Nymphsinf1) / sum(HC1Ha2$Totalnymphs) * 100 #Percentage of positive nymphs in CE+ Peak
sum(HC2Ha2$Nymphsinf1) / sum(HC2Ha2$Totalnymphs) * 100 #Percentage of positive nymphs in CE- Peak

sum(HC1La2$Adultsinf1) / sum(HC1La2$Totaladults) * 100 #Percentage of positive adults in CE+ Low
sum(HC2La2$Adultsinf1) / sum(HC2La2$Totaladults) * 100 #Percentage of positive adults in CE- Low
sum(HC1Ha2$Adultsinf1) / sum(HC1Ha2$Totaladults) * 100 #Percentage of positive adults in CE+ Peak
sum(HC2Ha2$Adultsinf1) / sum(HC2Ha2$Totaladults) * 100 #Percentage of positive adults in CE- Peak

#Number of infected ticks per host community
sum(HC1La2$Larvaeinf1) #Number of infected larvae in CE+ Low
sum(HC2La2$Larvaeinf1) #Number of infected larvae in CE- Low
sum(HC1Ha2$Larvaeinf1) #Number of infected larvae in CE+ Peak
sum(HC2Ha2$Larvaeinf1) #Number of infected larvae in CE- Peak

sum(HC1La2$Nymphsinf1) #Number of infected nymphs in CE+ Low
sum(HC2La2$Nymphsinf1) #Number of infected nymphs in CE- Low
sum(HC1Ha2$Nymphsinf1) #Number of infected nymphs in CE+ Peak
sum(HC2Ha2$Nymphsinf1) #Number of infected nymphs in CE- Peak

sum(HC1La2$Adultsinf1) #Number of infected adults in CE+ Low
sum(HC2La2$Adultsinf1) #Number of infected adults in CE- Low
sum(HC1Ha2$Adultsinf1) #Number of infected adults in CE+ Peak
sum(HC2Ha2$Adultsinf1) #Number of infected adults in CE- Peak

#Percentage of taxomomic host species in their relative contribution to the production of AP positive feeding ticks
#Order = Medium sized bird, Medium sized mammal, Small mammal, Ungulate
data.frame(aggregate(Larvaeinf1~Host.tax, HC1La2, FUN = sum))$Larvaeinf1 / sum(HC1La2$Larvaeinf1) * 100 #Percentage of infected larvae per host taxonomic group in CE+ Low
data.frame(aggregate(Larvaeinf1~Host.tax, HC2La2, FUN = sum))$Larvaeinf1 / sum(HC2La2$Larvaeinf1) * 100 #Percentage of infected larvae per host taxonomic group in CE- Low
data.frame(aggregate(Larvaeinf1~Host.tax, HC1Ha2, FUN = sum))$Larvaeinf1 / sum(HC1Ha2$Larvaeinf1) * 100 #Percentage of infected larvae per host taxonomic group in CE+ Peak
data.frame(aggregate(Larvaeinf1~Host.tax, HC2Ha2, FUN = sum))$Larvaeinf1 / sum(HC2Ha2$Larvaeinf1) * 100 #Percentage of infected larvae per host taxonomic group in CE- Peak

data.frame(aggregate(Nymphsinf1~Host.tax, HC1La2, FUN = sum))$Nymphsinf1 / sum(HC1La2$Nymphsinf1) * 100 #Percentage of infected nymphs per host taxonomic group in CE+ Low
data.frame(aggregate(Nymphsinf1~Host.tax, HC2La2, FUN = sum))$Nymphsinf1 / sum(HC2La2$Nymphsinf1) * 100 #Percentage of infected nymphs per host taxonomic group in CE- Low
data.frame(aggregate(Nymphsinf1~Host.tax, HC1Ha2, FUN = sum))$Nymphsinf1 / sum(HC1Ha2$Nymphsinf1) * 100 #Percentage of infected nymphs per host taxonomic group in CE+ Peak
data.frame(aggregate(Nymphsinf1~Host.tax, HC2Ha2, FUN = sum))$Nymphsinf1 / sum(HC2Ha2$Nymphsinf1) * 100 #Percentage of infected nymphs per host taxonomic group in CE- Peak

data.frame(aggregate(Adultsinf1~Host.tax, HC1La2, FUN = sum))$Adultsinf1 / sum(HC1La2$Adultsinf1) * 100 #Percentage of infected adults per host taxonomic group in CE+ Low
data.frame(aggregate(Adultsinf1~Host.tax, HC2La2, FUN = sum))$Adultsinf1 / sum(HC2La2$Adultsinf1) * 100 #Percentage of infected adults per host taxonomic group in CE- Low
data.frame(aggregate(Adultsinf1~Host.tax, HC1Ha2, FUN = sum))$Adultsinf1 / sum(HC1Ha2$Adultsinf1) * 100 #Percentage of infected adults per host taxonomic group in CE+ Peak
data.frame(aggregate(Adultsinf1~Host.tax, HC2Ha2, FUN = sum))$Adultsinf1 / sum(HC2Ha2$Adultsinf1) * 100 #Percentage of infected adults per host taxonomic group in CE- Peak

# # Ecotype 3 ####
HC1La3 <- ImpAPL(HC1, HC1$Ecotype3)
HC2La3 <- ImpAPL(HC2, HC2$Ecotype3)
HC1Ha3 <- ImpAPH(HC1, HC1$Ecotype3)
HC2Ha3 <- ImpAPH(HC2, HC2$Ecotype3)

#Making graphs for number of ticks
HC1La3$HC <- "CE+ Low"
HC2La3$HC <- "CE- Low"
HC1Ha3$HC <- "CE+ Peak"
HC2Ha3$HC <- "CE- Peak"
HC <- rbind(HC1La3, HC2La3, HC1Ha3, HC2Ha3) #This should always be called HC, otherwise the function plot.Imp does not work
HC$Host.tax <- factor(HC$Host.tax, levels = c("Small bird", "Medium sized bird", "Small mammal", "Medium sized mammal", "Ungulate")) #To reorder the classes
HC$HC <- factor(HC$HC, levels = c("CE+ Low", "CE+ Peak", "CE- Low", "CE- Peak")) #To reorder the host communities

plot.eco3 <-
  annotate_figure(
    ggarrange(annotate_figure(plot.Imp(HC$Larvaeinf1, "Number of infected feeding ticks") + ylim(0,max) + scale_fill_manual(values = my_colors) +
                                theme(legend.position = "none", axis.title=element_text(size=labsize), axis.text=element_text(size=txtsize)), top = text_grob("Larvae", size = labsize)),
              annotate_figure(plot.Imp(HC$Nymphsinf1, "") + ylim(0,max) + scale_fill_manual(values = my_colors) +
                                theme(legend.position = "none", axis.title=element_text(size=labsize), axis.text=element_text(size=txtsize)), top = text_grob("Nymphs", size = labsize)),
              annotate_figure(plot.Imp(HC$Adultsinf1, "") + ylim(0,max) + scale_fill_manual(values = my_colors) +
                                theme(legend.position = "none", axis.title=element_text(size=labsize), axis.text=element_text(size=txtsize)), top = text_grob("Adults", size = labsize)),
              nrow=1, ncol=3), 
    top = text_grob("Ecotype 3", size = grobsize))
plot.eco3

#Percentage of positive ticks of all feeding ticks
sum(HC1La3$Larvaeinf1) / sum(HC1La3$Totallarvae) * 100 #Percentage of positive larvae in CE+ Low
sum(HC2La3$Larvaeinf1) / sum(HC2La3$Totallarvae) * 100 #Percentage of positive larvae in CE- Low
sum(HC1Ha3$Larvaeinf1) / sum(HC1Ha3$Totallarvae) * 100 #Percentage of positive larvae in CE+ Peak
sum(HC2Ha3$Larvaeinf1) / sum(HC2Ha3$Totallarvae) * 100 #Percentage of positive larvae in CE- Peak

sum(HC1La3$Nymphsinf1) / sum(HC1La3$Totalnymphs) * 100 #Percentage of positive nymphs in CE+ Low
sum(HC2La3$Nymphsinf1) / sum(HC2La3$Totalnymphs) * 100 #Percentage of positive nymphs in CE- Low
sum(HC1Ha3$Nymphsinf1) / sum(HC1Ha3$Totalnymphs) * 100 #Percentage of positive nymphs in CE+ Peak
sum(HC2Ha3$Nymphsinf1) / sum(HC2Ha3$Totalnymphs) * 100 #Percentage of positive nymphs in CE- Peak

sum(HC1La3$Adultsinf1) / sum(HC1La3$Totaladults) * 100 #Percentage of positive adults in CE+ Low
sum(HC2La3$Adultsinf1) / sum(HC2La3$Totaladults) * 100 #Percentage of positive adults in CE- Low
sum(HC1Ha3$Adultsinf1) / sum(HC1Ha3$Totaladults) * 100 #Percentage of positive adults in CE+ Peak
sum(HC2Ha3$Adultsinf1) / sum(HC2Ha3$Totaladults) * 100 #Percentage of positive adults in CE- Peak

#Number of infected ticks per host community
sum(HC1La3$Larvaeinf1) #Number of infected larvae in CE+ Low
sum(HC2La3$Larvaeinf1) #Number of infected larvae in CE- Low
sum(HC1Ha3$Larvaeinf1) #Number of infected larvae in CE+ Peak
sum(HC2Ha3$Larvaeinf1) #Number of infected larvae in CE- Peak

sum(HC1La3$Nymphsinf1) #Number of infected nymphs in CE+ Low
sum(HC2La3$Nymphsinf1) #Number of infected nymphs in CE- Low
sum(HC1Ha3$Nymphsinf1) #Number of infected nymphs in CE+ Peak
sum(HC2Ha3$Nymphsinf1) #Number of infected nymphs in CE- Peak

sum(HC1La3$Adultsinf1) #Number of infected adults in CE+ Low
sum(HC2La3$Adultsinf1) #Number of infected adults in CE- Low
sum(HC1Ha3$Adultsinf1) #Number of infected adults in CE+ Peak
sum(HC2Ha3$Adultsinf1) #Number of infected adults in CE- Peak

#Percentage of taxomomic host species in their relative contribution to the production of AP positive feeding ticks
#Order = Medium sized bird, Medium sized mammal, Small mammal, Ungulate
data.frame(aggregate(Larvaeinf1~Host.tax, HC1La3, FUN = sum))$Larvaeinf1 / sum(HC1La3$Larvaeinf1) * 100 #Percentage of infected larvae per host taxonomic group in CE+ Low
data.frame(aggregate(Larvaeinf1~Host.tax, HC2La3, FUN = sum))$Larvaeinf1 / sum(HC2La3$Larvaeinf1) * 100 #Percentage of infected larvae per host taxonomic group in CE- Low
data.frame(aggregate(Larvaeinf1~Host.tax, HC1Ha3, FUN = sum))$Larvaeinf1 / sum(HC1Ha3$Larvaeinf1) * 100 #Percentage of infected larvae per host taxonomic group in CE+ Peak
data.frame(aggregate(Larvaeinf1~Host.tax, HC2Ha3, FUN = sum))$Larvaeinf1 / sum(HC2Ha3$Larvaeinf1) * 100 #Percentage of infected larvae per host taxonomic group in CE- Peak

data.frame(aggregate(Nymphsinf1~Host.tax, HC1La3, FUN = sum))$Nymphsinf1 / sum(HC1La3$Nymphsinf1) * 100 #Percentage of infected nymphs per host taxonomic group in CE+ Low
data.frame(aggregate(Nymphsinf1~Host.tax, HC2La3, FUN = sum))$Nymphsinf1 / sum(HC2La3$Nymphsinf1) * 100 #Percentage of infected nymphs per host taxonomic group in CE- Low
data.frame(aggregate(Nymphsinf1~Host.tax, HC1Ha3, FUN = sum))$Nymphsinf1 / sum(HC1Ha3$Nymphsinf1) * 100 #Percentage of infected nymphs per host taxonomic group in CE+ Peak
data.frame(aggregate(Nymphsinf1~Host.tax, HC2Ha3, FUN = sum))$Nymphsinf1 / sum(HC2Ha3$Nymphsinf1) * 100 #Percentage of infected nymphs per host taxonomic group in CE- Peak

data.frame(aggregate(Adultsinf1~Host.tax, HC1La3, FUN = sum))$Adultsinf1 / sum(HC1La3$Adultsinf1) * 100 #Percentage of infected adults per host taxonomic group in CE+ Low
data.frame(aggregate(Adultsinf1~Host.tax, HC2La3, FUN = sum))$Adultsinf1 / sum(HC2La3$Adultsinf1) * 100 #Percentage of infected adults per host taxonomic group in CE- Low
data.frame(aggregate(Adultsinf1~Host.tax, HC1Ha3, FUN = sum))$Adultsinf1 / sum(HC1Ha3$Adultsinf1) * 100 #Percentage of infected adults per host taxonomic group in CE+ Peak
data.frame(aggregate(Adultsinf1~Host.tax, HC2Ha3, FUN = sum))$Adultsinf1 / sum(HC2Ha3$Adultsinf1) * 100 #Percentage of infected adults per host taxonomic group in CE- Peak

# # Ecotype 4 ####
HC1La4 <- ImpAPL(HC1, HC1$Ecotype4)
HC2La4 <- ImpAPL(HC2, HC2$Ecotype4)
HC1Ha4 <- ImpAPH(HC1, HC1$Ecotype4)
HC2Ha4 <- ImpAPH(HC2, HC2$Ecotype4)

#Making graphs for number of ticks
HC1La4$HC <- "CE+ Low"
HC2La4$HC <- "CE- Low"
HC1Ha4$HC <- "CE+ Peak"
HC2Ha4$HC <- "CE- Peak"
HC <- rbind(HC1La4, HC2La4, HC1Ha4, HC2Ha4) #This should always be called HC, otherwise the function plot.Imp does not work
HC$Host.tax <- factor(HC$Host.tax, levels = c("Small bird", "Medium sized bird", "Small mammal", "Medium sized mammal", "Ungulate")) #To reorder the classes
HC$HC <- factor(HC$HC, levels = c("CE+ Low", "CE+ Peak", "CE- Low", "CE- Peak")) #To reorder the host communities

plot.eco4 <-
  annotate_figure(
    ggarrange(annotate_figure(plot.Imp(HC$Larvaeinf1, "Number of infected feeding ticks") + ylim(0,max) + scale_fill_manual(values = my_colors) +
                                theme(legend.position = "none", axis.title=element_text(size=labsize), axis.text=element_text(size=txtsize)), top = text_grob("Larvae", size = labsize)),
              annotate_figure(plot.Imp(HC$Nymphsinf1, "") + ylim(0,max) + scale_fill_manual(values = my_colors) +
                                theme(legend.position = "none", axis.title=element_text(size=labsize), axis.text=element_text(size=txtsize)), top = text_grob("Nymphs", size = labsize)),
              annotate_figure(plot.Imp(HC$Adultsinf1, "") + ylim(0,max) + scale_fill_manual(values = my_colors) +
                                theme(legend.position = "none", axis.title=element_text(size=labsize), axis.text=element_text(size=txtsize)), top = text_grob("Adults", size = labsize)),
              nrow=1, ncol=3), 
    top = text_grob("Ecotype 4", size = grobsize))
plot.eco4

#Percentage of positive ticks of all feeding ticks
sum(HC1La4$Larvaeinf1) / sum(HC1La4$Totallarvae) * 100 #Percentage of positive larvae in CE+ Low
sum(HC2La4$Larvaeinf1) / sum(HC2La4$Totallarvae) * 100 #Percentage of positive larvae in CE- Low
sum(HC1Ha4$Larvaeinf1) / sum(HC1Ha4$Totallarvae) * 100 #Percentage of positive larvae in CE+ Peak
sum(HC2Ha4$Larvaeinf1) / sum(HC2Ha4$Totallarvae) * 100 #Percentage of positive larvae in CE- Peak

sum(HC1La4$Nymphsinf1) / sum(HC1La4$Totalnymphs) * 100 #Percentage of positive nymphs in CE+ Low
sum(HC2La4$Nymphsinf1) / sum(HC2La4$Totalnymphs) * 100 #Percentage of positive nymphs in CE- Low
sum(HC1Ha4$Nymphsinf1) / sum(HC1Ha4$Totalnymphs) * 100 #Percentage of positive nymphs in CE+ Peak
sum(HC2Ha4$Nymphsinf1) / sum(HC2Ha4$Totalnymphs) * 100 #Percentage of positive nymphs in CE- Peak

sum(HC1La4$Adultsinf1) / sum(HC1La4$Totaladults) * 100 #Percentage of positive adults in CE+ Low
sum(HC2La4$Adultsinf1) / sum(HC2La4$Totaladults) * 100 #Percentage of positive adults in CE- Low
sum(HC1Ha4$Adultsinf1) / sum(HC1Ha4$Totaladults) * 100 #Percentage of positive adults in CE+ Peak
sum(HC2Ha4$Adultsinf1) / sum(HC2Ha4$Totaladults) * 100 #Percentage of positive adults in CE- Peak

#Number of infected ticks per host community
sum(HC1La4$Larvaeinf1) #Number of infected larvae in CE+ Low
sum(HC2La4$Larvaeinf1) #Number of infected larvae in CE- Low
sum(HC1Ha4$Larvaeinf1) #Number of infected larvae in CE+ Peak
sum(HC2Ha4$Larvaeinf1) #Number of infected larvae in CE- Peak

sum(HC1La4$Nymphsinf1) #Number of infected nymphs in CE+ Low
sum(HC2La4$Nymphsinf1) #Number of infected nymphs in CE- Low
sum(HC1Ha4$Nymphsinf1) #Number of infected nymphs in CE+ Peak
sum(HC2Ha4$Nymphsinf1) #Number of infected nymphs in CE- Peak

sum(HC1La4$Adultsinf1) #Number of infected adults in CE+ Low
sum(HC2La4$Adultsinf1) #Number of infected adults in CE- Low
sum(HC1Ha4$Adultsinf1) #Number of infected adults in CE+ Peak
sum(HC2Ha4$Adultsinf1) #Number of infected adults in CE- Peak

#Percentage of taxomomic host species in their relative contribution to the production of AP positive feeding ticks
#Order = Medium sized bird, Medium sized mammal, Small mammal, Ungulate
data.frame(aggregate(Larvaeinf1~Host.tax, HC1La4, FUN = sum))$Larvaeinf1 / sum(HC1La4$Larvaeinf1) * 100 #Percentage of infected larvae per host taxonomic group in CE+ Low
data.frame(aggregate(Larvaeinf1~Host.tax, HC2La4, FUN = sum))$Larvaeinf1 / sum(HC2La4$Larvaeinf1) * 100 #Percentage of infected larvae per host taxonomic group in CE- Low
data.frame(aggregate(Larvaeinf1~Host.tax, HC1Ha4, FUN = sum))$Larvaeinf1 / sum(HC1Ha4$Larvaeinf1) * 100 #Percentage of infected larvae per host taxonomic group in CE+ Peak
data.frame(aggregate(Larvaeinf1~Host.tax, HC2Ha4, FUN = sum))$Larvaeinf1 / sum(HC2Ha4$Larvaeinf1) * 100 #Percentage of infected larvae per host taxonomic group in CE- Peak

data.frame(aggregate(Nymphsinf1~Host.tax, HC1La4, FUN = sum))$Nymphsinf1 / sum(HC1La4$Nymphsinf1) * 100 #Percentage of infected nymphs per host taxonomic group in CE+ Low
data.frame(aggregate(Nymphsinf1~Host.tax, HC2La4, FUN = sum))$Nymphsinf1 / sum(HC2La4$Nymphsinf1) * 100 #Percentage of infected nymphs per host taxonomic group in CE- Low
data.frame(aggregate(Nymphsinf1~Host.tax, HC1Ha4, FUN = sum))$Nymphsinf1 / sum(HC1Ha4$Nymphsinf1) * 100 #Percentage of infected nymphs per host taxonomic group in CE+ Peak
data.frame(aggregate(Nymphsinf1~Host.tax, HC2Ha4, FUN = sum))$Nymphsinf1 / sum(HC2Ha4$Nymphsinf1) * 100 #Percentage of infected nymphs per host taxonomic group in CE- Peak

data.frame(aggregate(Adultsinf1~Host.tax, HC1La4, FUN = sum))$Adultsinf1 / sum(HC1La4$Adultsinf1) * 100 #Percentage of infected adults per host taxonomic group in CE+ Low
data.frame(aggregate(Adultsinf1~Host.tax, HC2La4, FUN = sum))$Adultsinf1 / sum(HC2La4$Adultsinf1) * 100 #Percentage of infected adults per host taxonomic group in CE- Low
data.frame(aggregate(Adultsinf1~Host.tax, HC1Ha4, FUN = sum))$Adultsinf1 / sum(HC1Ha4$Adultsinf1) * 100 #Percentage of infected adults per host taxonomic group in CE+ Peak
data.frame(aggregate(Adultsinf1~Host.tax, HC2Ha4, FUN = sum))$Adultsinf1 / sum(HC2Ha4$Adultsinf1) * 100 #Percentage of infected adults per host taxonomic group in CE- Peak

# # Graph of all ecotypes together ####
Plot.RelImpAP <-
ggarrange(
  ggarrange(
    plot.eco1, "",
    plot.eco2, "","","",
    plot.eco3, "",
    plot.eco4,
    nrow=3, ncol=3, widths = c(11,0.5,11), heights = c(11,0.5,11)),
  get_legend(plot.Imp(HC$Adultsinf1, "") +
               scale_fill_manual(values = my_colors,
                                 labels = c("Medium sized bird (22-42cm)", "Small mammal (<1kg)", "Medium sized mammal (1-20kg)", "Ungulate")) +
               theme(legend.text = element_text(size = 16))), 
  nrow=2, ncol=1, heights = c(11,1)) #The warning message "Cannot convert object of class character into a grob" can be ignored
Plot.RelImpAP