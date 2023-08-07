## Analysis of snail temperature experiment ##


library(readr)
library(stats)
library(lme4)
library(lmerTest)
library(lsmeans)
library(emmeans)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(broom)
library(AICcmodavg)
library(multcomp)
library(MASS)
library(MCMCglmm)
library(dplyr)

#Load data
F0<-read.csv("F0stacked.csv")
F1<-read.csv("F1stacked.csv")


F0$Volume_conesum <- as.numeric(F0$Volume_conesum)
F0$p.treat <- as.factor(F0$p.treat)

#############################
# plots for data exploration #
#############################

#Scatter plot by group
#total no. of eggs x total mass area
massegg <- F0 %>%
ggplot(aes(x = totalegg, y = totalmass, color = p.treat)) +
  geom_point()
massegg

#Size difference between parents and clutch number
sizextotalegg <- F0 %>%
  ggplot(aes(x = r.sizediff, y = totalegg, color = p.treat)) +
  geom_point()
sizextotalegg

#clutch number and total eggs laid
clutchxtotalegg <- F0 %>%
  ggplot(aes(x = clutches, y = totalegg, color = p.treat)) +
  geom_point()
clutchxtotalegg

# Length and width F0 snail sizes
snailsize <- F0 %>%
  ggplot(aes(x = clutches, y = totalegg, color = p.treat)) +
  geom_point()
snailsize

eggsizevnumber <- F0 %>%
  ggplot(aes(x = totalegg, y = ave.eggsize, color = p.treat)) +
  geom_point() + theme_bw()
eggsizevnumber
  
simple <- lm(ave.clutch ~ ave.eggsize, data = F0)
summary(simple)

## even without taking parental temperature into account, there is no relationship between egg size and total number of eggs
  

#Make a ladder plot to show how egg production goes over number of clutches
#load data
cdata<-read.csv("~/Dropbox/PHD/Snails/Analysis/clutchd.csv")

ggplot(data=cdata, aes(x=clutch, y=number, group=P.ID)) +
  geom_line(aes(color=treat))+ geom_point() +
  scale_colour_manual(values=c("#3399FF", "#5ECD69", "#EB6060"), labels=c("Cold", "Fluctuating", "Hot"),
                      name= "Parental treatment") +
  ylab(bquote('Egg number')) +
  labs(x="Clutch number", title= "Number of eggs per clutch") +
  theme_bw()



########################
##### Run the models #####
########################

#But first z transform the volume of snail pair within treatment groups

# Function to perform Z-transform within groups
z_transform_within_group <- function(x) {
  mean_x <- mean(x)
  sd_x <- sd(x)
  z_transformed_x <- (x - mean_x) / sd_x
  return(z_transformed_x)
}

# Apply Z-transform within groups using dplyr
snail_volumez <- F0 %>%
  group_by(p.treat) %>%
  mutate(Volume_conesum = z_transform_within_group(Volume_conesum))

###############################
## Parental fecundity models ##
###############################

#Subset data for some variables that have been entered across multiple lines
#when needed for covariates

TEdata<- snail_volumez %>% distinct(p.ID, .keep_all = TRUE)


# Total number of eggs
totaleggmod <- lm(totalegg ~ Volume_conesum + p.treat, data = TEdata)
summary(totaleggmod)
anova(totaleggmod)
#nothing

## Use the full dataset
# Egg size
eggsizem <- lmer(eggsize ~ Volume_conesum + p.treat + (1|p.ID), data = snail_volumez)
summary(eggsizem)
anova(eggsizem)
#nothing


# Egg size standard deviation
eggsizesdm <- lmer(eggSD ~ Volume_conesum + p.treat + (1|tank) + (1|p.ID), data = snail_volumez)
summary(eggsizesdm)
anova(eggsizesdm)
#nothing

snail_volumez$clutch.size <- as.numeric(snail_volumez$clutch.size)

#subset the data 
CSdata <- snail_volumez %>%
  distinct(p.ID, clutch.no, .keep_all = TRUE)

# Clutch size
clutchm <- lmer(clutch.size ~ Volume_conesum + p.treat + (1|p.ID), data = CSdata)
summary(clutchm)
anova(clutchm)
# Parental treatment has a signficant p-value (P=0.02)

#Post hoc test
emmeans(clutchm, list(pairwise ~ p.treat), adjust = "tukey")


residuals <- residuals(ave.clutchm)
plot(fitted(ave.clutchm), residuals,
     xlab = "Fitted Values",
     ylab = "Residuals",
     main = "Residual Plot")


# Latency to lay eggs
latencym <- lm(egglatency ~ Volume_conesum + p.treat, data = TEdata)
summary(latencym)
anova(latencym)
#nothing






## Models excluded from the analysis for now

# Egg mass size
eggmassm <- lmer(mass.area ~ Volume_conesum + p.treat + (1|p.ID), data = CSdata)
summary(eggmassm)
anova(eggmassm)
#nothing

# Number of clutches
clutchesm <- glm(clutches ~ Volume_conesum + p.treat, family=poisson(link="log"), data = snail_volumez)
summary(clutchesm)
#nothing





############################
##Analysis for F1 snails ####
############################
F1$parent <- as.factor(F1$parent)
F1$survival <- as.factor(F1$survival)
hist(F1$survival)

#Survival
survivalmod <- glmer(survival ~ p.treat * o.treat + (1 | parent) + (1 | o.ID) , data = F1, family = binomial)

summary(survivalmod)
anova(survivalmod)
#nothing
#singular fit


## Incubation duration model
IDmod <- lmer(ID ~ p.treat*o.treat + (1|parent), data = F1)
summary(IDmod)
anova(IDmod)
#Offspring treatment significant, pvalue =tiny

emmeans(IDmod, list(pairwise ~ o.treat), adjust = "tukey")


## Offspring size model
sizemod <- lmer(offspring.size ~ p.treat*o.treat + eggsize + (1|parent), data = F1)
summary(sizemod)
anova(sizemod)

#offspring treatment significant, egg size significant



#post hoc test
emmeans(sizemod, list(pairwise ~ o.treat), adjust = "tukey")

## AICc 
aictreat <- lmer(size.ave ~ o.treat + (1|parent), data = F1)
aicsize <- lmer(size.ave ~ ave.eggsize + (1|parent), data = F1)

AIC(aictreat)
AIC(aicsize)

#Size variation (SD)

#subset the data 
SDdata <- F1 %>%
  distinct(o.ID, offspring.sizesd, .keep_all = TRUE)

sizevmod <- lmer(offspring.sizesd ~ p.treat*o.treat + (1|parent), data = SDdata)
summary(sizevmod)
anova(sizevmod)
# nothing



#######################################
########       Graphics     ##########
#######################################

#Making some grouped box plots

##Offspring size - offspring x parental treatment

osize <- F1 %>%
  ggplot(aes(x=o.treat, y=size.ave, fill=factor(p.treat))) +
  ylab(bquote('Average offspring size'~(mm))) +
  labs(x="Developmental treatment", title= "", fill="Parental treatment") +
  scale_fill_manual(values=c("#4e5154", "#ced1d6", "white"), labels = c("Cold", "Fluctuating", "Hot")) + 
  scale_y_continuous(expand=c(0.0,0.0), limits=c(0.6, 1)) +
  scale_x_discrete(labels = c("Cold", "Fluctuating", "Hot")) +
  theme_bw() +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank()) +
  theme(axis.text = element_text(size = 10)) +
  theme(axis.title = element_text(size = 10)) +
  theme(plot.title = element_text(size = 12)) +
  theme(legend.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 12)) +     
  geom_boxplot() + theme_bw()
osize

sizesd <- F1 %>%
  ggplot(aes(x=o.treat, y=size.sd, fill=factor(p.treat))) +
  ylab(bquote('Average offspring SD'~(mm))) +
  labs(x="Developmental treatment", title= "", fill="Parental treatment") +
  scale_fill_manual(values=c("#4e5154", "#ced1d6", "white"), labels = c("Cold", "Fluctuating", "Hot")) + 
  scale_y_continuous(expand=c(0.0,0.0), limits=c(0, 0.18)) +
  scale_x_discrete(labels = c("Cold", "Fluctuating", "Hot")) +
  theme_bw() +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank()) +
  theme(axis.text = element_text(size = 10)) +
  theme(axis.title = element_text(size = 10)) +
  theme(plot.title = element_text(size = 12)) +
  theme(legend.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 12)) +     
  geom_boxplot() + theme_bw()
sizesd



## These were my predictions when investigating the interaction between parental x offspring environment,
## so arrange these into a panel plot

ggarrange(osize, sizesd, 
          labels = c("A", "B"),
          common.legend = TRUE, legend = "right", align = "v",
          ncol = 2, nrow = 1)



## Parental investment plot

clutchplot <- F0 %>%
  ggplot(aes(x = p.treat, y = ave.clutch, fill = factor(p.treat))) +
  geom_boxplot() +
  ylab(bquote("Average clutch size")) +
  labs(x = "Parental treatment") +
  scale_x_discrete(labels = c("Cold", "Fluctuating", "Hot")) +
  scale_fill_manual(values = c("#4e5154", "#ced1d6", "white"), labels = c("Cold", "Fluctuating", "Hot")) + 
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "none",  # Hide the legend
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        plot.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12)) 
clutchplot

eggsizeplot <- F0 %>%
  ggplot(aes(x = p.treat, y = ave.eggsize, fill = factor(p.treat))) +
  geom_boxplot() +
  ylab(bquote("Average egg size")) +
  labs(x = "Parental treatment") +
  scale_x_discrete(labels = c("Cold", "Fluctuating", "Hot")) +
  scale_fill_manual(values = c("#4e5154", "#ced1d6", "white"), labels = c("Cold", "Fluctuating", "Hot")) + 
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "none",  # Hide the legend
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        plot.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12)) 
eggsizeplot



## Show the relationship we found between offspring treatment, egg size and offspring size

#Grouped scatterplot?

scatterbw <- F1 %>%
  ggplot(aes(x = ave.eggsize, y = size.ave)) +
  geom_point(aes(shape = o.treat, color= o.treat), size = 3.5) +
  scale_shape_manual(values = c(16,19,1), name= "Offspring treatment",
                     labels=c("Cold", "Fluctuating", "Hot")) +
  scale_color_manual(values = c("black", "grey", "black"), name= "Offspring treatment",
                     labels=c("Cold", "Fluctuating", "Hot"))+
  ylab(bquote('Average offspring size'~(mm))) +
  labs(x="Average egg size"~(mm)) +
  theme_bw() +
  theme(legend.title = element_text(size=15), legend.text = element_text(size=12), axis.text=element_text(size=12),axis.title=element_text(size=12))
scatterbw

## a colour plot

scattercol <- F1 %>%
  ggplot(aes(x = size.ave, y = ave.eggsize)) +
  geom_point(aes(color= o.treat), size = 3.5) +
  ylab(bquote('Average egg size'~(mm))) +
  labs(x="Average offspring size"~(mm)) +
  scale_color_manual(values = c("#2EB5F3", "#C885E8", "#F06423"), name= "Parental treatment",
                     labels=c("Cold", "Fluctuating", "Hot"))+
  theme_bw() +
  theme(legend.title = element_text(size=17), legend.text = element_text(size=15), axis.text=element_text(size=12),axis.title=element_text(size=12))
scattercol






## Trying to figure out what is going with the size measurements

lengthT <- F0 %>%
  ggplot(aes(x = l.size, y = ave.clutch, color = p.treat)) +
  geom_point() + theme_bw()
lengthT

volumeT <- F0 %>%
  ggplot(aes(x = Volume_conesum, y = ave.clutch, color = p.treat)) +
  geom_point() + theme_bw()
volumeT

ggarrange(lengthT, volumeT, 
          labels = c("A", "B"),
          common.legend = TRUE, legend = "right", align = "v",
          ncol = 2, nrow = 1)

F0 %>%
  ggplot(aes(x = Volume_conesum, y = l.size, color = p.treat)) +
  geom_point() + theme_bw()


group_tbl <- F0 %>% 
  group_by(p.treat) %>% 
  summarise(total_count = sum(Volume_conesum),
            .groups = 'drop')

snailvol <- lm(Volume_conesum~p.treat, data= F0)
summary(snailvol)
anova(snailvol)

snailvol2 <- aov(snailvol)
TukeyHSD(snailvol2, "p.treat")


##################################################
####### Not being used at the moment ###############
###################################################



## Run a PCA on the parental fecundity variables
library(psych)

#load data
F0PCA<-read.csv("~/Dropbox/PHD/Snails/Analysis/F0PCA.csv")

## Look at the correlations
pairs.panels(F0PCA[,-5],
             gap = 0,
             bg = c("red", "yellow", "blue"), 
             pch=21)


pca2 <-princomp(F0PCA)


PCAload <-loadings(pca2)
print(PCAload, digits = 3, cutoff = 0, sort = FALSE)

### Using the prcomp function as the proportion variance looks weird
PCA3 <-prcomp(F0PCA, scale= TRUE)
loadingsPCA3 <- PCA3$rotation

print(PCA3)
variance <- PCA3$sdev^2
prop_variance <- variance / sum(variance)
print(prop_variance)




## Calculate confidence intervals for the survival parental effects plot

# Group by p.treat and o.treat and calculate mean and standard error
summary_data <- dfs.mean %>%
  group_by(p.treat, o.treat) %>%
  summarise(mean_proportion = mean(survival),
            n = n())

# Calculate standard error
summary_data <- summary_data %>%
  group_by(p.treat) %>%
  mutate(se = sqrt(mean_proportion * (1 - mean_proportion) / n))

# Calculate confidence intervals (e.g., using 95% confidence level)
z_value <- qnorm(0.975)  # For 95% confidence level (two-tailed)
summary_data <- summary_data %>%
  group_by(p.treat) %>%
  mutate(lower_ci = mean_proportion - z_value * se,
         upper_ci = mean_proportion + z_value * se)

# Print the summary data with confidence intervals
print(summary_data)

# Merge summary_data with dfs.mean based on p.treat and o.treat
merged_data <- merge(dfs.mean, summary_data, by = c("p.treat", "o.treat"), all.x = TRUE)

# Print the merged data
print(merged_data)


ggplot(data=merged_data, aes(x=p.treat, y=survival, group=o.treat, shape=o.treat)) +
  geom_line() +
  geom_point(aes(colour=o.treat), size=3.5) +
  geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci), width=0.2, colour="black", size=1) +
  scale_colour_manual(values=c("black", "grey", "black"), labels=c("Cold", "Fluctuating", "Hot"),
                      name= "Developmental treatment") +
  scale_shape_manual(values=c(16, 19, 1), name= "Developmental treatment",
                     labels=c("Cold", "Fluctuating", "Hot")) +
  ylab(bquote('Proportion survived')) +
  labs(x="Parental treatment", title= "") +
  scale_x_discrete(labels=c("Cold", "Fluctuating", "Hot")) +
  theme_bw() +
  theme(legend.title = element_text(size=17),
        legend.text = element_text(size=15),
        axis.text = element_text(size=12),
        axis.title = element_text(size=12))



#try grouping by catagories within columns
#Group by mean of multiple columns
dfs.mean <- F1 %>% group_by(p.treat, o.treat) %>% 
  summarise(across(c(survival),mean),
            .groups = 'drop') %>%
  as.data.frame()

#IT WORKS

## Survival line plot

ggplot(data=merged_data, aes(x=p.treat, y=survival, group=o.treat, shape=o.treat)) +
  geom_line() + geom_point(aes(colour=o.treat), size=3.5) + 
  geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci), width=0.2, colour="black", size=1) +
  scale_colour_manual(values=c("black", "grey", "black"), labels=c("Cold", "Fluctuating", "Hot"),
                      name= "Developmental treatment") +
  scale_shape_manual(values = c(16,19,1), name= "Developmental treatment",
                     labels=c("Cold", "Fluctuating", "Hot")) + 
  ylab(bquote('Proportion survived')) +
  labs(x="Parental treatment", title= "") +
  scale_x_discrete(labels = c("Cold", "Fluctuating", "Hot")) +
  theme_bw() +
  theme(legend.title = element_text(size=17), legend.text = element_text(size=15), axis.text=element_text(size=12),axis.title=element_text(size=12))

