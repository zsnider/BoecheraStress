### Boechera Stress Experiments ###

setwd("~/Desktop/StressExperiments")

#Loading in the packages
library(ggplot2)
library(tidyr)
library(hrbrthemes)
library(plotly)
library(tibble)
library(forcats)
library(scales)
library(nlme)
library(lme4)

### ION LEAKAGE ANALYSES ###

#Loading the Stress and Recovery Ion Leakage Data
stress <- read.csv("jrpam_all.csv")
recovery <- read.csv("jrpam_5day.csv")
colnames(stress) <- c('Individual','PopID','Region','Species','Heat.Temp','Cotyledon.vs.TL',
                           'rel.ms','Date.Stress','Time.Stress','Fo.prime','No.',
                      'F','Fm.','Temp','YII','Fo','Fm','Fv.Fm')
colnames(recovery) <- c('Individual','PopID','Region','Species','Heat.Temp','Cotyledon.vs.TL',
                      'rel.ms','Date.Stress','Time.Stress','Fo.prime','No.',
                      'F','Fm.','Temp','YII','Fo','Fm','Fv.Fm')
recovery$Individual <- as.factor(recovery$Individual)
stress$Individual <- as.factor(stress$Individual)
stress$YII <- as.numeric(stress$YII)
recovery$YII <- as.numeric(recovery$YII)
summary(stress$YII)
stress <- stress %>% 
  filter(Heat.Temp == "44")
recovery <- recovery %>%
  filter(Heat.Temp == "44")
head(stress)
head(recovery)

###Stress plot by Individual
stressplot1 <- ggplot(data=stress, aes(x=as.factor(Individual), y=YII, fill=as.factor(Region))) +
  geom_boxplot() +
  theme_bw() +
  ylab("YII") +
  xlab("Individual ID") +
  labs(title = "Initial Stress Measurements") +
  labs(fill="Heat Stress (C°)") +
  scale_fill_brewer(palette="Set2")
stressplot1
###Recovery by Individual
recoveryplot1 <- ggplot(data=recovery, aes(x=as.factor(Individual), y=as.numeric(YII), fill=as.factor(Region))) +
  geom_boxplot() +
  theme_bw() +
  ylab("YII") +
  xlab("Individual") +
  labs(title = "5 Day Recovery Measurements") +
  labs(fill="Heat Stress (C°)") +
  scale_fill_brewer(palette="Set2")
recoveryplot1

#Stress averaged for each subpop
stressplot2 <-  ggplot(data=stress, aes(x=as.factor(PopID), y=YII, fill=as.factor(Region))) +
  geom_boxplot() +
  theme_bw() +
  ylab("YII") +
  xlab("Population") +
  labs(title = "Initial Stress Measurements") +
  labs(fill="Heat Stress (C°)") +
  scale_fill_brewer(palette="Set2")
stressplot2

#Recovery averaged for each subpop
recoveryplot2 <- ggplot(data=recovery, aes(x=as.factor(PopID), y=as.numeric(YII), fill=as.factor(Region))) +
  geom_boxplot() +
  theme_bw() +
  ylab("YII") +
  xlab("Population") +
  labs(title = "5 Day Recovery Measurements") +
  labs(fill="Heat Stress (C°)") +
  scale_fill_brewer(palette="Set2")
recoveryplot2

###PHENOTYPE DATA###

#Load the data
pheno <- read.csv("PhenotypeScores.csv")
pheno1 <- pivot_longer(pheno, cols = 2:4, names_to = "Time", values_to = "Score")
head(pheno1)
##Phenotype by Individual
pheno2 <- aggregate(Score ~ Individual + Time, 
                 data=pheno1, 
                 function(x) { 
                   c(avg=mean(x)) 
                 })

pheno2$Individual <- as.factor(pheno2$Individual)
pheno2$Time <- factor(pheno2$Time, levels=c('Pre.Stress', 'Post.Stress', 'Post.Recovery'))
phenoplot1 <- ggplot(data=pheno2, aes(x=Time, y=Individual,fill=Score)) +
  geom_tile(color="black", lwd=0.08) +
  scale_fill_gradient2(low="#318e1e", mid = "#ffff55", high="#ffffe9", 
                       midpoint = 1.5, na.value = "grey50") +
  theme_ipsum() +
  ylab("Individual") +
  xlab("Time") +
  labs(title = "Phenotype Score") +
  scale_x_discrete(labels=c("Pre.Stress" = "Pre Stress", "Post.Stress" = "After Stress",
                            "Post.Recovery" = "After 5 Day Recovery"))
phenoplot1
ggplotly(phenoplot1)

## Phenotype Scores by Subpopulation
head(pheno1)
pheno3 <- aggregate(Score ~ PopID + Time, 
                    data=pheno1, 
                    function(x) { 
                      c(avg=mean(x)) 
                    })
head(pheno3)
pheno3$PopID <- as.factor(pheno3$PopID)
pheno3$Time <- factor(pheno3$Time, levels=c('Pre.Stress', 'Post.Stress', 'Post.Recovery'))
phenoplot2 <- ggplot(data=pheno3, aes(x=Time, y=PopID,fill=Score)) +
  geom_tile(color="black", lwd=0.08) +
  scale_fill_gradient2(low="#318e1e", mid = "#ffff55", high="#ffffe9", 
                       midpoint = 1.5, na.value = "gray") +
  theme_ipsum() +
  ylab("Subpopulation") +
  xlab("Time") +
  labs(title = "Phenotype Score") +
  scale_x_discrete(labels=c("Pre.Stress" = "Pre Stress", "Post.Stress" = "After Stress",
                            "Post.Recovery" = "After 5 Day Recovery"))
phenoplot2
ggplotly(phenoplot2)



###ION LEAKAGE PLOT ####
ion <- read.csv("IonLeakage.csv")
ion$Percent.Leaked <- as.numeric(ion$Percent.Leaked)
ion$Individual <- as.factor(ion$Individual)
ion$PopID <- as.factor(ion$PopID)

ionplot <- ggplot(ion, aes(x=PopID, y=Percent.Leaked, fill = Region)) +
  geom_boxplot() +
  theme_bw() +
  xlab("Subpopulation") +
  ylab("Percent of Ions Leaked After 4 Hours")
ionplot



#Without cotyledons, Only True Leaves
stressTL <- stress %>% 
  filter(Cotyledon.vs.TL == "True Leaf")
recoveryTL <- recovery %>%
  filter(Cotyledon.vs.TL == "True Leaf")

stressplotTL <- ggplot(data=stressTL, aes(x=as.factor(PopID), y=YII, fill=as.factor(Region))) +
  geom_boxplot() +
  theme_bw() +
  ylab("YII") +
  xlab("Individual ID") +
  labs(title = "Initial Stress Measurements (True Leaves)") +
  labs(fill="Heat Stress (C°)") +
  scale_fill_brewer(palette="Set2")
stressplotTL

recoveryplotTL <- ggplot(data=recoveryTL, aes(x=as.factor(PopID), y=as.numeric(YII), fill=as.factor(Region))) +
  geom_boxplot() +
  theme_bw() +
  ylab("YII") +
  xlab("Individual") +
  labs(title = "5 Day Recovery Measurements (True Leaves)") +
  labs(fill="Heat Stress (C°)") +
  scale_fill_brewer(palette="Set2")
recoveryplotTL

#Only Cotyledon
stressCTY <- stress %>% 
  filter(Cotyledon.vs.TL == "Cotyledon")
recoveryCTY <- recovery %>%
  filter(Cotyledon.vs.TL == "Cotyledon")

stressplotCTY <- ggplot(data=stressCTY, aes(x=as.factor(PopID), y=YII, fill=as.factor(Region))) +
  geom_boxplot() +
  theme_bw() +
  ylab("YII") +
  xlab("Individual ID") +
  labs(title = "Initial Stress Measurements (Cotyledons)") +
  labs(fill="Heat Stress (C°)") +
  scale_fill_brewer(palette="Set2")
stressplotCTY

recoveryplotCTY <- ggplot(data=recoveryCTY, aes(x=as.factor(PopID), y=as.numeric(YII), fill=as.factor(Region))) +
  geom_boxplot() +
  theme_bw() +
  ylab("YII") +
  xlab("Individual") +
  labs(title = "5 Day Recovery Measurements (Cotyledons)") +
  labs(fill="Heat Stress (C°)") +
  scale_fill_brewer(palette="Set2")
recoveryplotCTY

#True Leaf vs Cotyledon comparison
stressplot3 <- stress %>%
  drop_na(Cotyledon.vs.TL) %>%
  ggplot(aes(x=as.factor(Region), y=YII, fill=as.factor(Cotyledon.vs.TL))) +
  geom_boxplot() +
  theme_bw() +
  ylab("YII") +
  xlab("Individual ID") +
  labs(title = "Initial Stress Measurements (Cotyledons vs True Leaves)") +
  labs(fill="Type of Leaf") +
  scale_fill_manual(values = c("#8DA0CB","#66C2A5","#FC8D62"))
stressplot3

recoveryplot3 <- ggplot(data=recovery, aes(x=as.factor(Region), y=as.numeric(YII),
                                           fill=as.factor(Cotyledon.vs.TL))) +
  geom_boxplot() +
  theme_bw() +
  ylab("YII") +
  xlab("Individual") +
  labs(title = "5 Day Recovery Measurements (Cotyledons vs True Leaves)") +
  labs(fill="Type of Leaf") +
  scale_fill_brewer(palette="Set2")
recoveryplot3


##Running a multi way ANOVA

#First checking for normality
shapiro.test(recovery$YII)
#shapiro.test(recovery$YII[Region=="San Diego"])
recoverySD <- recovery %>%
  filter(Region == "San Diego")
shapiro.test(recoverySD$YII)
qqnorm(recoverySD$YII);qqline(recoverySD$YII)
recoveryLA <- recovery %>%
  filter(Region == "Los Angeles")
shapiro.test(recoveryLA$YII)
qqnorm(recoveryLA$YII);qqline(recoveryLA$YII)
recoverySJ <- recovery %>%
  filter(Region == "San Jacinto")
shapiro.test(recoverySJ$YII)
qqnorm(recoverySJ$YII);qqline(recoverySJ$YII)
recoverySB <- recovery %>%
  filter(Region == "San Bernadino")
shapiro.test(recoverySB$YII)
qqnorm(recoverySB$YII);qqline(recoverySB$YII)
#Test for homoscedascity
bartlett.test(recovery$YII~recovery$Region)



m1<-aov(recovery$YII~recovery$Region)
summary(m1)
TukeyHSD(m1)
pairwise.t.test(recovery$YII,recovery$Region,p.adj="bonf") #bonferroni adjustment
pairwise.t.test(recovery$YII,recovery$Region,p.adj="holm") #holm adjustment, a conservative adjustment that might increase type 1 error rate slightly
pairwise.t.test(recovery$YII,recovery$Region,p.adj="none")

recovery.ca <- recovery %>%
  filter(Species == "californica")
#Test for normality
rcSD <- recovery.ca %>%
  filter(Region == "San Diego")
shapiro.test(rcSD$YII)
qqnorm(rcSD$YII);qqline(rcSD$YII)
rcLA <- recovery.ca %>%
  filter(Region == "Los Angeles")
shapiro.test(rcLA$YII)
qqnorm(rcLA$YII);qqline(rcLA$YII)
rcSJ <- recovery.ca %>%
  filter(Region == "San Jacinto")
shapiro.test(rcSJ$YII)
qqnorm(rcSJ$YII);qqline(rcSJ$YII)
#Test for homoscedascity
bartlett.test(recovery.ca$YII~recovery.ca$Region)

m2<-aov(recovery.ca$YII~recovery.ca$Region)
summary(m2)
TukeyHSD(m2)
pairwise.t.test(recovery.ca$YII,recovery.ca$Region,p.adj="bonf") #bonferroni adjustment
pairwise.t.test(recovery.ca$YII,recovery.ca$Region,p.adj="holm") #holm adjustment, a conservative adjustment that might increase type 1 error rate slightly
pairwise.t.test(recovery.ca$YII,recovery.ca$Region,p.adj="none")






#Perform test of variance of individual's YII
#Each individual should have 1 number associated with it, it's variance
#Test of variances
#Nested anova region
#Region + Pop(Region) + Individual(Pop(region))

#Yield
#Linear mixed effects model
#Each individual is a random effect
#multiple measurements on a random effect
#measurements on an  individual are not independent
#average of measurements on one seedling
#Look up LMER with identity
#Can probably just do an LME



head(recovery2)
recovery2 <- aggregate(YII ~ Individual, 
                    data=recovery, 
                    function(x) { 
                      c(avg=var(x)) 
                    })
colnames(recovery2) <- c('Individual','YII.var')
recovery3 <- left_join(recovery2, recovery, by=c("Individual"))
head(recovery3)
recovery3[!duplicated(recovery3$Individual), ]

nest <- aov(recovery3$YII.var ~ recovery3$Region / factor(recovery3$PopID) / factor(recovery3$Individual))

summary(nest)
TukeyHSD(nest)








###MANOVA
library(car)
# s2 <- read.csv("jrpam2_10.03.2022.csv") #S2 only contains pops at 45
climate <- read.csv("boechera_climate_data.csv")
stress.clim <- left_join(stress, climate, by=c("Individual","Region","PopID"))

#scatterplotMatrix(~Elevation+tmean1+PPT1+tdmean1|Pop,
                  #data = s2, smooth=FALSE, regLine = FALSE, ellipse=TRUE, 
                  #by.groups=TRUE, diagonal=FALSE, legend = list(coords= 'bottomleft'), 
                  #cex.axis = 1.5, las = 1, col = carPalette('car')) 
man1 = manova(cbind(Elevation, tmean1, PPT1, tdmean1)~X1.Y..II., data = stress.clim)
summary(man1) #error due to residuals
summary.aov(man1)

man2 = manova(cbind(Elevation, tmean1, tmax_1, tmin1,  PPT1, tdmean1, vpdmax1, 
                    vpdmin1, solclear1, soltrans1)~X1.Y..II., data = stress.clim)
summary(man2) #error due to residuals
summary.aov(man2)


m1 <- (lm(X1.Y..II. ~ Elevation * tmean1 * PPT1 * tdmean1, data = stress.clim))
anova(m1)
summary(m1)
qqnorm(residuals(m1))

anova(lm(X1.Y..II. ~ Elevation + tmean1 + tmax_1 + tmin1 + PPT1 + tdmean1 + vpdmax1 +
         vpdmin1 + solclear1 + soltrans1, data = stress.clim))

anova(lm(X1.Y..II. ~ Elevation * tmean1 * tmax_1 * tmin1 * PPT1 * tdmean1 * vpdmax1 *
           vpdmin1 * solclear1 * soltrans1, data = stress.clim))
###Multivariate multiple regression
M3 <- lm(cbind(Elevation, tmean1, PPT1, tdmean1)~ X1.Y..II., data = stress.clim)
summary(Anova(M3))

###One Way ANOVA
s1anova <- aov(stress.clim$PopID~stress.clim$X1.Y..II.) 
summary(s1anova)

###Factorial ANOVA
anova(lm(stress.clim$PopID~stress.clim$X1.Y..II. + stress.clim$X1.Fo + stress.clim$X1.Fm + stress.clim$X1.Fv.Fm), data = stress.clim)















