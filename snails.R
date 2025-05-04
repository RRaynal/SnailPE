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
library(glmmTMB)
library(dplyr)
library(car)
library(parameters)

#Load data
F0<-read.csv("F0stacked.csv")
F1<-read.csv("F1stacked.csv")
F1egg <-read.csv("F1egg.csv")

F0$Volume_conesum <- as.numeric(F0$Volume_conesum)
F0$p.treat <- as.factor(F0$p.treat)
F0$p.ID <- as.factor(F0$p.ID)
F0$clutch.no <- as.factor(F0$clutch.no)

F1$p.treat <- as.factor(F1$p.treat)
F1$parent <- as.factor(F1$parent)
F1$o.treat <- as.factor(F1$o.treat)




#############################
# plots for data exploration #
#############################

#Scatter plot by group
#total no. of eggs x total mass area
massegg <- F0 %>%
ggplot(aes(x = totalegg, y = mass.area, color = p.treat)) +
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
  ggplot(aes(x = clutch.size, y = eggsize, color = p.treat)) +
  geom_point() + theme_bw()
eggsizevnumber


## are egg size and egg number correlated (regardless of temp treatment)?
F0 <- na.omit(F0)
simplePP <- glmmTMB(clutch.size ~ eggsize,  data = F0)
summary(simplePP)

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
## average snail volume per treatment
average_volume <- F0 %>%
  group_by(p.treat) %>%
  summarise(
    avg_volume = mean(Volume_conesum, na.rm = TRUE),
    sd_volume = sd(Volume_conesum, na.rm = TRUE),
    .groups = 'drop'
  )


# z transform the volume of snail pair within treatment groups

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

F0$Volume_conesum <- as.numeric(F0$Volume_conesum)

###############################
## Parental fecundity models ##
###############################

#Subset data for some variables that have been entered across multiple lines
#when needed for covariates

TEdata<- snail_volumez %>% distinct(p.ID, .keep_all = TRUE)


#### Total number of eggs #######
totaleggmod <- glmmTMB(totalegg ~ Volume_conesum + p.treat + (1|p.ID), 
                       data = TEdata, 
                       family = poisson)
summary(totaleggmod)
Anova(totaleggmod, type = "II")
emmeans(totaleggmod, ~ p.treat)

## work out percentage of random effects for a poisson model
mean(TEdata$totalegg, na.rm = TRUE)

# 1. Extract variance of the random effect (from your model summary)
var_ID <- 0.2362

# 2. Calculate the mean of totalegg
mu <- mean(TEdata$totalegg, na.rm = TRUE)

# 3. Approximate the residual variance for Poisson-log model
var_resid <- log(1 + 1 / mu)

# 4. Calculate ICC
ICC <- var_ID / (var_ID + var_resid)

# 5. Print result
cat("Mean totalegg =", round(mu, 2), "\n")
cat("Residual variance ≈", round(var_resid, 4), "\n")
cat("ICC for p.ID ≈", round(ICC, 4), "or", round(ICC * 100, 1), "%\n")

#P.ID = 87.1%

# Extract fitted values and residuals
fitted_vals <- fitted(totaleggmod)
resid_vals <- residuals(totaleggmod, type = "pearson")

# Residuals vs. Fitted plot
ggplot(data = data.frame(Fitted = fitted_vals, Residuals = resid_vals), 
       aes(x = Fitted, y = Residuals)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Residuals vs. Fitted Plot", x = "Fitted Values", y = "Residuals") +
  theme_minimal()

## Use the full dataset
############ Egg size ##############
eggsizem <- glmmTMB(eggsize ~ Volume_conesum + p.treat + (1|p.ID) + (1|clutch.no), data = snail_volumez, family = gaussian(link = "identity"))
summary(eggsizem)
Anova(eggsizem, type = "II")


# Extract fitted values and residuals
fitted_vals <- fitted(eggsizem)
resid_vals <- residuals(eggsizem, type = "pearson")

# Residuals vs. Fitted plot
ggplot(data = data.frame(Fitted = fitted_vals, Residuals = resid_vals), 
       aes(x = Fitted, y = Residuals)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Residuals vs. Fitted Plot", x = "Fitted Values", y = "Residuals") +
  theme_minimal()
## good

emmeans(eggsizem, ~ p.treat)

## Percentages of random effects ##
p.ID    = 0.002097
clutch  = 0.00008871
resid   = 0.008804

icc_pID = 0.002097 / (0.002097 + 0.00008871 + 0.008804)
icc_clutch = 0.00008871 / (0.002097 + 0.00008871 + 0.008804)


######### Egg size coefficient of variation #########

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
eggsizesdm <- glmmTMB(cv ~ Volume_conesum + p.treat + (1|p.ID), data = F0SD, family = gaussian(link = "identity"))
summary(eggsizesdm)
Anova(eggsizesdm, type = "II")
#nothing
diagnose(eggsizesdm)

# Extract fitted values and residuals
fitted_vals <- fitted(eggsizesdm)
resid_vals <- residuals(eggsizesdm, type = "pearson")

# Residuals vs. Fitted plot
ggplot(data = data.frame(Fitted = fitted_vals, Residuals = resid_vals), 
       aes(x = Fitted, y = Residuals)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Residuals vs. Fitted Plot", x = "Fitted Values", y = "Residuals") +
  theme_minimal()

#good

## Percent variance of random effects
# Variance components
var_ID <- 1.798
var_resid <- 10.866

# Total variance
var_total <- var_ID + var_resid

# ICC / proportion of variance due to p.ID
ICC <- var_ID / var_total

# result
cat("Total variance:", round(var_total, 3), "\n")
cat("ICC for p.ID =", round(ICC, 4), "or", round(ICC * 100, 1), "%\n")



################ Clutch size###########################################

snail_volumez$clutch.size <- as.numeric(snail_volumez$clutch.size)

#subset the data 
CSdata <- snail_volumez %>%
  distinct(p.ID, clutch.no, .keep_all = TRUE)

# Clutch size
clutchm <- glmmTMB(clutch.size ~ Volume_conesum + p.treat + (1|p.ID), data = CSdata, family = gaussian(link = "identity"))
summary(clutchm)

Anova(clutchm, type = "II")
# Parental treatment has a signficant p-value (P=0.02)

# Extract fitted values and residuals
fitted_vals <- fitted(clutchm)
resid_vals <- residuals(clutchm, type = "pearson")

# Residuals vs. Fitted plot
ggplot(data = data.frame(Fitted = fitted_vals, Residuals = resid_vals), 
       aes(x = Fitted, y = Residuals)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Residuals vs. Fitted Plot", x = "Fitted Values", y = "Residuals") +
  theme_minimal()

# good

#Post hoc test
emmeans(clutchm, list(pairwise ~ p.treat), adjust = "tukey")
emmeans(clutchm, ~ p.treat)

## random effects as percentage variance ##

# Store variance components
var_pID <- 2.816e-11
var_clutch <- 0.003794
var_residual <- 0.006636

# Total variance
total_var <- var_pID + var_clutch + var_residual

# Percentages
percent_pID <- 100 * var_pID / total_var
percent_clutch <- 100 * var_clutch / total_var
percent_residual <- 100 * var_residual / total_var

# Print
cat("Percent variance explained:\n")
cat(sprintf("p.ID: %.4f%%\n", percent_pID))
cat(sprintf("clutch.no: %.2f%%\n", percent_clutch))
cat(sprintf("Residual: %.2f%%\n", percent_residual))

TEdata <- na.omit(TEdata)

# Latency to lay eggs
latencym <- glmmTMB(egglatency ~ Volume_conesum + p.treat + (1|p.ID), data = TEdata, family = gaussian(link = "identity"))
summary(latencym)
Anova(latencym, type = "II")
#nothing

# Extract fitted values and residuals
fitted_vals <- fitted(latencym)
resid_vals <- residuals(latencym, type = "pearson")

# Residuals vs. Fitted plot
ggplot(data = data.frame(Fitted = fitted_vals, Residuals = resid_vals), 
       aes(x = Fitted, y = Residuals)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Residuals vs. Fitted Plot", x = "Fitted Values", y = "Residuals") +
  theme_minimal()



############################
##Analysis for F1 snails ####
############################
F1$parent <- as.factor(F1$parent)
F1$tank.ID <- as.factor(F1$tank.ID)


#Survival
survivalmod <- glmmTMB(cbind(hatched, unhatched) ~ p.treat * o.treat + (1| parent) + (1|tank.ID), 
                       data = F1, 
                       family = binomial)
                     
summary(survivalmod)
Anova(survivalmod, type = "II")

# Extract fitted values and residuals
fitted_vals <- fitted(survivalmod)
resid_vals <- residuals(survivalmod, type = "pearson")

# Residuals vs. Fitted plot
ggplot(data = data.frame(Fitted = fitted_vals, Residuals = resid_vals), 
       aes(x = Fitted, y = Residuals)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Residuals vs. Fitted Plot", x = "Fitted Values", y = "Residuals") +
  theme_minimal()

#Post hoc test
emmeans(survivalmod, pairwise ~ o.treat, adjust = "tukey")

#Calculate percent change between o.treat levels
emmeans(survivalmod, ~ o.treat, type = "response")

## calculate percent variance explained by random effects

# Variance components
var_parent <- 0.5184
var_tank <- 4.0061

# Fixed residual variance on the latent scale for logistic models
var_residual <- (pi^2) / 3  # ≈ 3.29

# Total variance on the latent scale
var_total <- var_parent + var_tank + var_residual

# Percent variance explained
percent_parent <- (var_parent / var_total) * 100
percent_tank <- (var_tank / var_total) * 100
percent_residual <- (var_residual / var_total) * 100

# Combine into a data frame for display
variance_summary <- data.frame(
  Component = c("parent", "tank.ID", "Residual"),
  Variance = c(var_parent, var_tank, var_residual),
  Percent = round(c(percent_parent, percent_tank, percent_residual), 2)
)

print(variance_summary)

#Work out the sample sizes for each o.treat

observation_counts <- F1 %>%
  group_by(o.treat) %>%
  summarise(total_observations = sum(hatched, na.rm = TRUE) + sum(unhatched, na.rm = TRUE))

observation_counts



## Incubation duration model

F1_ID <- F1 %>%
  drop_na(ID) %>%
  droplevels()



IDmod <- glmmTMB(ID ~ p.treat * o.treat + (1 | parent), 
                 data = F1_ID, 
                 family = gaussian(link = "identity"))
summary(IDmod)
Anova(IDmod, type = "II")  

#calculate percent change
emmeans(IDmod, ~ o.treat)

#calculate the percentage variance explained by random effects

# Variance components 
var_parent <- 0.2177
var_residual <- 3.9428

# Total variance
var_total <- var_parent + var_residual

# Percent variance explained
percent_parent <- (var_parent / var_total) * 100
percent_residual <- (var_residual / var_total) * 100

# Combine into a data frame
variance_summary <- data.frame(
  Component = c("parent", "Residual"),
  Variance = c(var_parent, var_residual),
  Percent = round(c(percent_parent, percent_residual), 2)
)

print(variance_summary)


# Extract fitted values and residuals
fitted_vals <- fitted(IDmod)
resid_vals <- residuals(IDmod, type = "pearson")

# Residuals vs. Fitted plot
ggplot(data = data.frame(Fitted = fitted_vals, Residuals = resid_vals), 
       aes(x = Fitted, y = Residuals)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Residuals vs. Fitted Plot", x = "Fitted Values", y = "Residuals") +
  theme_minimal()

#Offsprintg treatment significant, pvalue =tiny

emmeans(IDmod, list(pairwise ~ o.treat), adjust = "tukey")

emmeans(IDmod, list(pairwise ~ p.treat), adjust = "tukey")



Sizedata <- F1 %>%
  filter(!is.na(offspring.size)) %>%
  distinct(tank.ID, offspring.size, .keep_all = TRUE)

#### Offspring size model ####

sizemod <- glmmTMB(offspring.size ~ p.treat * o.treat + eggsize + (1 | parent) + (1 | tank.ID), 
                   data = Sizedata, 
                   family = gaussian(link = "identity"))

summary(sizemod)
Anova(sizemod, type = "II")

#calculate percent change
emmeans(sizemod, ~ o.treat)

# calculate percentage variance of random effects
# Variance components from the model
var_parent <- 1.437e-11
var_tank <- 1.327e-03
var_residual <- 3.988e-03

# Total variance
var_total <- var_parent + var_tank + var_residual

# Percent variance explained
percent_parent <- (var_parent / var_total) * 100
percent_tank <- (var_tank / var_total) * 100
percent_residual <- (var_residual / var_total) * 100

# Combine into a table
variance_summary <- data.frame(
  Component = c("parent", "tank.ID", "Residual"),
  Variance = c(var_parent, var_tank, var_residual),
  Percent = round(c(percent_parent, percent_tank, percent_residual), 4)
)

print(variance_summary)

# Extract fitted values and residuals
fitted_vals <- fitted(sizemod)
resid_vals <- residuals(sizemod, type = "pearson")

# Residuals vs. Fitted plot
ggplot(data = data.frame(Fitted = fitted_vals, Residuals = resid_vals), 
       aes(x = Fitted, y = Residuals)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Residuals vs. Fitted Plot", x = "Fitted Values", y = "Residuals") +
  theme_minimal()


offspring_count <- Sizedata %>%
  group_by(o.treat) %>%
  summarise(count = n())

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


sizevmod <- glmmTMB(cv ~ p.treat * o.treat + (1 | parent),
                    data = F1SDU, 
                    family = gaussian(link = "identity"))
summary(sizevmod)
Anova(sizevmod, type = "II")
# nothing
# Extract fitted values and residuals
fitted_vals <- fitted(sizevmod)
resid_vals <- residuals(sizvemod, type = "pearson")

# Residuals vs. Fitted plot
ggplot(data = data.frame(Fitted = fitted_vals, Residuals = resid_vals), 
       aes(x = Fitted, y = Residuals)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Residuals vs. Fitted Plot", x = "Fitted Values", y = "Residuals") +
  theme_minimal()

#calculating the percentage of variance explained by random effects
# Variance components from the model
var_parent <- 2.508
var_residual <- 9.551

# Total variance
var_total <- var_parent + var_residual

# Percent variance explained
percent_parent <- (var_parent / var_total) * 100
percent_residual <- (var_residual / var_total) * 100

# Combine into a table
variance_summary <- data.frame(
  Component = c("parent", "Residual"),
  Variance = c(var_parent, var_residual),
  Percent = round(c(percent_parent, percent_residual), 2)
)

print(variance_summary)

#######################################
########       Graphics     ##########
#######################################

#Making some grouped box plots

##Offspring size - offspring x parental treatment

osize <- F1 %>%
  ggplot(aes(x=p.treat, y=offspring.size, fill=factor(o.treat))) +
  ylab(bquote('Offspring size'~(mm))) +
  labs(x="Parental Treatment", title= "", fill="Developmental\ntreatment") +
  scale_fill_manual(values=c("#4e5154", "#ced1d6", "white"), labels = c("Cold", "Fluctuating", "Hot")) + 
  scale_y_continuous(expand=c(0.0,0.0), limits=c(0.6, 1)) +
  scale_x_discrete(labels = c("Cold", "Fluctuating", "Hot")) +
  geom_boxplot() +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    plot.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18)
  )
osize

#ID plot

IDplot <- F1_ID %>%
  ggplot(aes(x = p.treat, y = ID, fill = factor(o.treat))) +
  ylab(bquote('Incubation duration'~(days))) +
  labs(x = "Parental Treatment", title = "", fill = "Developmental\ntreatment") +
  scale_fill_manual(values = c("#4e5154", "#ced1d6", "white"), labels = c("Cold", "Fluctuating", "Hot")) + 
  scale_x_discrete(labels = c("Cold", "Fluctuating", "Hot")) +
  geom_boxplot() +
  scale_y_continuous(breaks = seq(floor(min(F1_ID$ID)), ceiling(max(F1_ID$ID)), by = 2)) + # Adjust as needed
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    plot.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18)
  )

IDplot



sizesd <- F1SDU %>%
  ggplot(aes(x = p.treat, y = cv, fill = factor(o.treat))) +
  ylab(bquote('Offspring Coefficient of variation')) +
  labs(x = "Parental treatment", title = "", fill = "Developmental\ntreatment") +
  scale_fill_manual(values = c("#4e5154", "#ced1d6", "white"), labels = c("Cold", "Fluctuating", "Hot")) + 
  scale_y_continuous(expand = c(0.0, 0.0), limits = c(0, 25)) +
  scale_x_discrete(labels = c("Cold", "Fluctuating", "Hot")) +
  geom_boxplot() +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    plot.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18)
  )
sizesd



## These were my predictions when investigating the interaction between parental x offspring environment,
## so arrange these into a panel plot

ggarrange(IDplot, osize, 
          labels = c("A", "B"),
          common.legend = TRUE, legend = "right", align = "v",
          ncol = 2, nrow = 1)



## Parental investment plot

clutchplot <- F0 %>%
  ggplot(aes(x = p.treat, y = clutch.size, fill = factor(p.treat))) +
  geom_boxplot() +
  ylab(bquote("Clutch size (no. of eggs)")) +
  labs(x = "Parental treatment") +
  scale_x_discrete(labels = c("Cold", "Fluctuating", "Hot")) +
  scale_fill_manual(values = c("#4e5154", "#ced1d6", "white"), labels = c("Cold", "Fluctuating", "Hot")) + 
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "none",  # Hide the legend
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        plot.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18)) 
clutchplot

eggsizeplot <- F0 %>%
  ggplot(aes(x = p.treat, y = eggsize, fill = factor(p.treat))) +
  geom_boxplot() +
  ylab(bquote("Egg size (mm)")) +
  labs(x = "Parental treatment") +
  scale_x_discrete(labels = c("Cold", "Fluctuating", "Hot")) +
  scale_fill_manual(values = c("#4e5154", "#ced1d6", "white"), labels = c("Cold", "Fluctuating", "Hot")) + 
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "none",  # Hide the legend
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        plot.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18)) 
eggsizeplot

# GGarrange, while we found nothing for the trade off between egg size and clutch size
#its a major theme throughout.

ggarrange(eggsizeplot, clutchplot, 
          labels = c("A", "B"),
          ncol = 2, nrow = 1)

## Show the relationship we found between offspring treatment, egg size and offspring size

#Grouped scatterplot?

#average offspring size per tank
avg_offspring_size <- F1 %>%
  group_by(tank.ID) %>%
  summarise(mean_offspring_size = mean(offspring.size, na.rm = TRUE))

F1 <- F1 %>%
  group_by(tank.ID) %>%
  mutate(mean_offspring_size = mean(offspring.size, na.rm = TRUE))


scatterbw <- ggplot(F1_clean, aes(x = eggsize, y = mean_offspring_size)) +
  geom_point(aes(shape = o.treat, color = o.treat), size = 3.5) +
  scale_shape_manual(values = c(16, 19, 1), name = "Developmental\ntreatment",
                     labels = c("Cold", "Fluctuating", "Hot")) +
  scale_color_manual(values = c("black", "grey", "black"), name = "Developmental\ntreatment",
                     labels = c("Cold", "Fluctuating", "Hot")) +
  ylab(bquote('Average offspring size' ~ (mm))) +
  labs(x = "Average egg size" ~ (mm)) +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        plot.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18)) +
  geom_smooth(method = "lm", color = "black", se = FALSE) +  # Add regression line
  annotate("text", x = max(F1_clean$eggsize) * 0.96, y = max(F1_clean$mean_offspring_size) * 0.95,  # Move label to right
           label = paste("R² =", R2_value, "\n", p_label), size = 4, hjust = 0)
scatterbw


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
  geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci), width=0.2, colour="black", 
                size=1, position = position_dodge(width = 0.1)) +
  geom_point(aes(colour=o.treat, fill=o.treat), size=5, stroke=1, 
             position = position_dodge(width = 0.1)) + 
  scale_shape_manual(values=c(21, 21, 21),
                     name="Developmental\ntreatment",
                     labels=c("Cold", "Fluctuating", "Hot")) +
  scale_fill_manual(values=c("black", "grey", "white"),  
                    name="Developmental\ntreatment",
                    labels=c("Cold", "Fluctuating", "Hot")) +
  scale_colour_manual(values=c("black", "black", "black"),  # Border colors
                      name="Developmental\ntreatment",
                      labels=c("Cold", "Fluctuating", "Hot")) +
  ylab(bquote('Hatching succcess')) +
  labs(x="Parental treatment", title="") +
  scale_x_discrete(labels=c("Cold", "Fluctuating", "Hot")) +
  theme_bw() +
  theme(legend.title = element_text(size=18),
        legend.text = element_text(size=16),
        axis.text = element_text(size=14),
        axis.title = element_text(size=16))
 



ggplot(data = merged_data, aes(x = p.treat, y = survival_rate.x, group = o.treat, shape = o.treat)) + 
  geom_line() + 
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), 
                width = 0.2, 
                colour = "black", 
                size = 1, 
                position = position_jitter(width = 0.1)) +  # Add jitter to CI
  geom_jitter(aes(colour = o.treat, fill = o.treat), 
              size = 5, 
              stroke = 1, 
              width = 0.1) +  # Add jitter to points
  scale_shape_manual(values = c(21, 21, 21),
                     name = "Developmental\ntreatment", 
                     labels = c("Cold", "Fluctuating", "Hot")) +
  scale_fill_manual(values = c("black", "grey", "white"),   
                    name = "Developmental\ntreatment", 
                    labels = c("Cold", "Fluctuating", "Hot")) +
  scale_colour_manual(values = c("black", "black", "black"),  # Border colors
                      name = "Developmental\ntreatment", 
                      labels = c("Cold", "Fluctuating", "Hot")) +
  ylab(bquote('Proportion hatched')) + 
  labs(x = "Parental treatment", title = "") + 
  scale_x_discrete(labels = c("Cold", "Fluctuating", "Hot")) + 
  theme_bw() + 
  theme(legend.title = element_text(size = 18), 
        legend.text = element_text(size = 16), 
        axis.text = element_text(size = 14), 
        axis.title = element_text(size = 16))









## are the treatments significantly different from one another?

#Load data
Rooms <-read.csv("treatmenttest_long.csv")

Treat_diff <- aov(Value ~ Treatment, data = Rooms)
summary(Treat_diff)

emmeans(Treat_diff, list(pairwise ~ Treatment), adjust = "tukey")

#number of offspring tanks per combination

no.offspring.T <- F1 %>%
  group_by(tank.ID) %>%
  slice(1) %>%
  ungroup()

no.offspring.T %>%
  count(p.treat, o.treat)

# no. of eggs in every combination of treatments

F1egg %>%
  group_by(p.treat, o.treat) %>%
  summarize(total_egg_no = sum(egg_no, na.rm = TRUE), .groups = "drop")