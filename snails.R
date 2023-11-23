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

#calculate the coefficient of variance from the SD and means

F0coef <- snail_volumez %>%
  group_by(tank) %>%
  summarise(mean_egg_size = mean(eggsize),
            cv = eggSD / mean_egg_size * 100)

#subset the data 
F0coefU <- F0coef %>%
  distinct(tank, cv, .keep_all = TRUE)

#add coefficient of variance to the data
data_merged <- snail_volumez %>%
  left_join(F0coefU, by = "tank")

#subset the data 
F0SD <- data_merged %>%
  distinct(tank, cv, .keep_all = TRUE)



# Egg size coef var model
eggsizesdm <- lmer(cv ~ Volume_conesum + p.treat + (1|p.ID), data = F0SD)
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








############################
##Analysis for F1 snails ####
############################
F1$parent <- as.factor(F1$parent)
F1$tank.ID <- as.factor(F1$tank.ID)

#We will be doing AIC selection for these models, so need an intercept to compare them

survivalint <- glmer(cbind(hatched, unhatched) ~ (1|parent), 
                     data = F1, 
                     family = binomial)

summary(survivalint)

##   AIC      BIC   logLik 
## 848.0    857.4   -422.0 


#Survival
survivalmod <- glmer(cbind(hatched, unhatched) ~ p.treat * o.treat + (1|parent), 
                     data = F1, 
                     family = binomial)
                     
summary(survivalmod)
anova(survivalmod)

#AIC = 801.7, BIC= 848.7

# test the interaction
survivalmod2 <- glmer(cbind(hatched, unhatched) ~ p.treat + o.treat + (1|parent), 
                     data = F1, 
                     family = binomial)

summary(survivalmod2)

#AIC = 813.3, BIC= 841.5

#Test the moderators
survivalmod3 <- glmer(cbind(hatched, unhatched) ~ o.treat + (1|parent), 
                      data = F1, 
                      family = binomial)

summary(survivalmod3)

# AIC= 813.9, BIC= 932.7

#Test the moderators
survivalmod4 <- glmer(cbind(hatched, unhatched) ~ p.treat + (1|parent), 
                      data = F1, 
                      family = binomial)

summary(survivalmod4)

# AIC= 849.1, BIC=867.9 

## Interesting, it likes the model with the interaction the most.

#Work out the sample sizes for each o.treat

observation_counts <- F1 %>%
  group_by(o.treat) %>%
  summarise(total_observations = sum(hatched, na.rm = TRUE) + sum(unhatched, na.rm = TRUE))

observation_counts



## Incubation duration model
IDmod <- lmer(ID ~ p.treat*o.treat + (1|parent), data = F1)
summary(IDmod)
anova(IDmod)
#Offspring treatment significant, pvalue =tiny

emmeans(IDmod, list(pairwise ~ o.treat), adjust = "tukey")

## There were too many tanks with 3 or fewer  offspring size measurements, so we could not include
## tank as a random effect. (18 tanks) - Or exclude these


Sizedata <- F1 %>%
  filter(!is.na(offspring.size)) %>%
  distinct(tank.ID, offspring.size, .keep_all = TRUE)

## Offspring size model
sizemod <- lmer(offspring.size ~ p.treat*o.treat + eggsize + (1|parent), data = Sizedata)
summary(sizemod)
anova(sizemod)

offspring_count <- Sizedata %>%
  group_by(o.treat) %>%
  summarise(count = n())

offspring_count

offspring_count


#offspring treatment significant, egg size significant



#post hoc test
emmeans(sizemod, list(pairwise ~ o.treat), adjust = "tukey")


#calculate the coefficient of variance from the SD and means
## First get rid of the double ups
F1SD <- F1 %>%
  distinct(tank.ID, offspring.sizesd, .keep_all = TRUE)%>%
  filter(!is.na(offspring.sizesd))


F1coef <- F1SD %>%
  group_by(tank.ID) %>%
  summarise(cv = offspring.sizesd / offspring.size * 100)


#add coefficient of variance to the data
data_merged <- F1 %>%
  left_join(F1coef, by = "tank.ID")

#subset the data 
F1SDU <- data_merged %>%
  distinct(tank.ID, cv, .keep_all = TRUE)%>%
  filter(!is.na(cv))


#Size variation (coefficient of variation)


sizevmod <- lmer(cv ~ p.treat*o.treat + (1|parent), data = F1SDU)
summary(sizevmod)
anova(sizevmod)
# nothing



#######################################
########       Graphics     ##########
#######################################

#Making some grouped box plots

##Offspring size - offspring x parental treatment

osize <- F1 %>%
  ggplot(aes(x=p.treat, y=offspring.size, fill=factor(o.treat))) +
  ylab(bquote('Offspring size'~(mm))) +
  labs(x="Parental Treatment", title= "", fill="Developmental treatment") +
  scale_fill_manual(values=c("#4e5154", "#ced1d6", "white"), labels = c("Cold", "Fluctuating", "Hot")) + 
  scale_y_continuous(expand=c(0.0,0.0), limits=c(0.6, 1)) +
  scale_x_discrete(labels = c("Cold", "Fluctuating", "Hot")) +
  theme_bw() +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank()) +
  theme(axis.text = element_text(size = 10)) +
  theme(axis.title = element_text(size = 10)) +
  theme(plot.title = element_text(size = 12)) +
  theme(legend.text = element_text(size = 18)) +
  theme(legend.title = element_text(size = 20)) +     
  geom_boxplot() + theme_bw()
osize

sizesd <- F1SDU %>%
  ggplot(aes(x=p.treat, y=cv, fill=factor(o.treat))) +
  ylab(bquote('Offspring Coefficient of variation')) +
  labs(x="Parental treatment", title= "", fill="Developmental treatment") +
  scale_fill_manual(values=c("#4e5154", "#ced1d6", "white"), labels = c("Cold", "Fluctuating", "Hot")) + 
  scale_y_continuous(expand=c(0.0,0.0), limits=c()) +
  scale_x_discrete(labels = c("Cold", "Fluctuating", "Hot")) +
  theme_bw() +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank()) +
  theme(axis.text = element_text(size = 10)) +
  theme(axis.title = element_text(size = 10)) +
  theme(plot.title = element_text(size = 12)) +
  theme(legend.text = element_text(size = 18)) +
  theme(legend.title = element_text(size = 20)) +     
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
  ggplot(aes(x = p.treat, y = clutch.size, fill = factor(p.treat))) +
  geom_boxplot() +
  ylab(bquote("Clutch size")) +
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
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 17)) 
clutchplot

eggsizeplot <- F0 %>%
  ggplot(aes(x = p.treat, y = eggsize, fill = factor(p.treat))) +
  geom_boxplot() +
  ylab(bquote("Egg size")) +
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

# GGarrange, I think while we found nothing for the trade off between egg size and clutch size
#its a major theme throughout.

ggarrange(eggsizeplot, clutchplot, 
          labels = c("A", "B"),
          ncol = 2, nrow = 1)

## Show the relationship we found between offspring treatment, egg size and offspring size

#Grouped scatterplot?



scatterbw <- F1 %>%
  ggplot(aes(x = eggsize, y = offspring.size)) +
  geom_point(aes(shape = o.treat, color= o.treat), size = 3.5) +
  scale_shape_manual(values = c(16,19,1), name= "Offspring treatment",
                     labels=c("Cold", "Fluctuating", "Hot")) +
  scale_color_manual(values = c("black", "grey", "black"), name= "Offspring treatment",
                     labels=c("Cold", "Fluctuating", "Hot"))+
  ylab(bquote('offspring size'~(mm))) +
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











#try grouping by catagories within columns
#Group by mean of multiple columns
# Group by mean of multiple columns
dfs.mean <- F1 %>%
  group_by(p.treat, o.treat) %>%
  summarise(survival_rate = mean(hatched, na.rm = TRUE),
            n = n(),
            .groups = 'drop') %>%
  as.data.frame()

dfs.mean <- dfs.mean %>%
  filter(!is.na(n) & n != 2)


# Calculate standard error
summary_data <- dfs.mean %>%
  group_by(p.treat) %>%
  mutate(se = sqrt(survival_rate * (1 - survival_rate) / n))

# Calculate confidence intervals (e.g., using 95% confidence level)
z_value <- qnorm(0.975)  # For 95% confidence level (two-tailed)
summary_data <- summary_data %>%
  group_by(p.treat) %>%
  mutate(lower_ci = survival_rate - z_value * se,
         upper_ci = survival_rate + z_value * se)


# Merge summary_data with dfs.mean based on p.treat and o.treat
merged_data <- merge(dfs.mean, summary_data, by = c("p.treat", "o.treat"), all.x = TRUE)


ggplot(data=merged_data, aes(x=p.treat, y=survival_rate.x, group=o.treat, shape=o.treat)) +
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

## Trying with raw values


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


