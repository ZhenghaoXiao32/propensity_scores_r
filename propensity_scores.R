library(tidyverse)
library(readxl)
library(mice)
library(party)
library(mitools)
library(twang)
library(lattice)
library(car)
library(broom)
library(treatSens)
library(Matching)
library(MatchIt)
library(optmatch)
library(rgenoud)
library(rbounds)
setwd("~/Documents/GitHub/propensity_scores_r")
patient <- read_excel("patientdata.xlsx")
str(patient)
head(patient, n = 10)
#====================Data imputation=====================================#
# all the data looks good,  
summary(patient$kg)
ggplot(patient, aes(kg)) +
   geom_histogram(bins = 50, color = "white")
# but some of the kg is 0, which is weird,
# as the plot shows, the majority of the kg column is between 50 to 100 kg,
# although I have no knowledge about how this data is collected, according to the 
# density plot, I can make a guess that those 0s are of missing data instead of real 
# world data. 
# we should replace it with NA first:

patient$kg <- na_if(patient$kg, 0)

# check how many NA we have in kg column
round(sum(is.na(patient$kg)) / length(patient$kg), 3)

# we have 9% missing data in kg column
# we need to create a na indicator for this column, 
# because the reason of causing missing data can be a true confounder
kg_na <- data.frame(is.na(patient$kg))
names(kg_na) <- c("kg_na")
patient <- cbind(patient, kg_na)

# scale and factorize our data
colnames(patient)
factor_cols <- c("severe", "cognitivedecline", "depression", "cancer", "autoimmune", "transferedin", 
                 "sex", "surgerytype", "nolifesupportorder", "insurancetype", "resp", 
                 "infection", "trauma", "race", "kg_na")
scale_cols <- c("age", "yearseducation", "bloodpressure", "temperature", "creatinelevels", 
                "sodiumlevels", "urineweight", "kg", "income")
patient_adj <- patient %>% 
      mutate_at(factor_cols, factor) %>%
      mutate_at(scale_cols, scale) %>% 
      mutate_at(scale_cols, as.numeric)
str(patient_adj)


# use multiple imputation to replace NAs in kg column
# impute variables by from least missing to most missing
# Using multiple imputation by chained equations
# with predictive mean matching ("pmm") as the univariate imputation method
# m is the amount of imputation, set a seed for reproducibility
imputed_data <- mice(patient_adj, m = 5, maxit = 50, method = "pmm", seed = 500)

summary(imputed_data)
# have a look of the imputed data:
imputed_data$imp$kg

imputed_data1 <- complete(imputed_data, 1)
imputed_data2 <- complete(imputed_data, 2)
imputed_data3 <- complete(imputed_data, 3)
imputed_data4 <- complete(imputed_data, 4)
imputed_data5 <- complete(imputed_data, 5)

#create a list of all imputed datasets
#this object was specifically designed to be analyzed with the survey package
all_imputations <- imputationList(list(imputed_data1, imputed_data2, imputed_data3,
                                       imputed_data4, imputed_data5))

#=========================Get ropensity cores==============================================#
### after imputation finished, we can first use logistic regression to get propensity scores
# let's get our formula first:
# this is just a test formula with all possible covariates
cov_names <- subset(colnames(patient_adj), !(colnames(patient_adj) %in% c("died", "surgerytype")))
ps_formula <- paste(cov_names, collapse = "+")
ps_formula <- formula(paste("surgerytype ~ ", ps_formula, sep = ""))
ps_formula

# since our data has no strata no clusters, we can just use the glm function
# let's start with 1st imputation data:
ps_model1 <- glm(ps_formula, data = imputed_data1, family = "binomial")
# have a look on the ps model summary 
(sum_1 <- summary(ps_model1))

# some of the predictors are not significant, we may decide to remove them.
# yeadeducaiton is just at verge, we can test it with other imputation data
ps_model2 <- glm(ps_formula, data = imputed_data2, family = "binomial")
(sum_2 <- summary(ps_model2))

ps_model3 <- glm(ps_formula, data = imputed_data3, family = "binomial")
(sum_3 <- summary(ps_model3))

ps_model4 <- glm(ps_formula, data = imputed_data4, family = "binomial")
(sum_4 <- summary(ps_model4))

ps_model5 <- glm(ps_formula, data = imputed_data5, family = "binomial")
(sum_5 <- summary(ps_model5))

# as a result, yearseducation is significant in 2 of the 5 models
# and the p-value of it is very close to 0.05
# so I decide to put it in our final ps model

# let's check multicolliearity of our ps model 
vif(ps_model1)
# all the vif values if close to 1
# we have no problem in multicollinearity
# find the significant predictors for treatment variable now
names(subset(coef(sum_1)[, 4], coef(sum_1)[, 4] < 0.05))
# we need to add yearseducation into this list:
cov_names_final <- c("severe", "cognitivedecline", "depression",  
                      "cancer", "transferedin", "bloodpressure", "yearseducation", 
                      "creatinelevels", "sodiumlevels", "urineweight", "kg", 
                      "nolifesupportorder", "insurancetype", "resp", "infection", 
                      "trauma", "kg_na")
ps_formula_final <- paste(cov_names_final, collapse = "+")
ps_formula_final <- formula(paste("surgerytype ~ ", ps_formula_final, sep = ""))
ps_formula_final

ps_model_final1 <- glm(ps_formula_final, data = imputed_data1, family = "binomial")
summary(ps_model_final1)
### get propensity scores from logistic regression
p_score1 <- fitted(ps_model_final1)
imputed_data1$p_scores <- p_score1


#now let's do it with all imputaion data:
svy_design_all <- svydesign(ids = ~0, data = all_imputations)
ps_model_all <- with(svy_design_all, svyglm(ps_formula_final, family = binomial()))

# we use the average of propensity scores from 5 imputation data as our final p score
p_score_all <- sapply(ps_model_all, fitted)
p_scores <- apply(p_score_all, 1, mean)

# bind the p scores into our data frame:
all_imputations <- update(all_imputations, p_scores = p_scores)
str(all_imputations)


## let's try recursive partitioning algorithms to get p scores

# first we need to 
# check the distributions of treatment variable and outcome variable
# if they are balance, the partitioning will be fine
table(patient_adj$died)
table(patient_adj$surgerytype)

my_ctree <- ctree(ps_formula_final, data = imputed_data1)
plot(my_ctree)


## also need to check for appropriate mtry value:
## recommended number of covariates per tree is the square root of the number of predictors
sqrt(ncol(patient) - 1)


# random forest
set.seed(2019)
my_controls <- cforest_unbiased(ntree = 1000, mtry = 5)
my_cforest <- cforest(ps_formula_final, data = imputed_data1, controls = my_controls)

p_score_rf <- predict(my_cforest, type = "prob")
imputed_data1$p_score_rf <- matrix(unlist(p_score_rf), , 2, byrow = TRUE)[, 2]
str(imputed_data1)

## generalize boosted modeling 
# for complete data 1
# first we need to convert our treatment variable back to numeric for ps function
set.seed(2019)
imputed_data1$surgerytype <- as.numeric(imputed_data1$surgerytype == 1)
my_gbm <- ps(ps_formula_final, data = imputed_data1, n.trees = 10000, interaction.depth = 4,
             stop.method = c("es.max"), estimand = "ATT", verbose = TRUE)
summary(my_gbm)
plot(my_gbm)
p_score_gbm <- my_gbm$ps
names(p_score_gbm) = "p_score_gbm"
imputed_data1$p_score_gbm <- unlist(p_score_gbm)

str(imputed_data1[, 26:28])
# I didn't convert the propensity scores to linear propensity scores,
# because the advantage in using linear propensity scores is for matching,
# but I'm only using weighting and stratification.


## factor treatment variable back
imputed_data1$surgerytype <- factor(imputed_data1$surgerytype)
str(imputed_data1)
### evaluation of common support:
# for the impuatation data 1:

# descriptive statistics:

with(imputed_data1, by(p_scores, surgerytype, summary))
with(imputed_data1, by(p_score_rf, surgerytype, summary))
with(imputed_data1, by(p_score_gbm, surgerytype, summary))


by(imputed_data1[, 26:28], imputed_data1$surgerytype, summary)
#From the statistics, we can find that propensity scores getting 
# from logistic regression has best common support.

# graphical checks:
# for logistic regression
bwplot(p_scores ~ surgerytype, data = imputed_data1, 
       ylab = "Propensity Scores by logistic regression",
       xlab = "Surgery Type", auto.key = TRUE)

imputed_data1 %>%
      ggplot(aes(x = p_scores)) +
      geom_histogram(color = "white") +
      facet_wrap(~ surgerytype) +
      xlab("Probability of surgery type: PS")

# for random forest
bwplot(p_score_rf ~ surgerytype, data = imputed_data1, 
       ylab = "Propensity Scores by random forest",
       xlab = "Surgery Type", auto.key = TRUE)

imputed_data1 %>%
   ggplot(aes(x = p_score_rf)) +
   geom_histogram(color = "white") +
   facet_wrap(~ surgerytype) +
   xlab("Probability of surgery type: PS")

# for generalize boosted modeling
bwplot(p_score_gbm ~ surgerytype, data = imputed_data1, 
       ylab = "Propensity Scores by generalized boosted modeling",
       xlab = "Surgery Type", auto.key = TRUE)

imputed_data1 %>%
   ggplot(aes(x = p_score_gbm)) +
   geom_histogram(color = "white") +
   facet_wrap(~ surgerytype) +
   xlab("Probability of surgery type: PS")

# compare them
densityplot( ~ p_scores, groups = surgerytype, plot.points = FALSE, 
             xlim = c(0,1), lwd = 2,
             data = imputed_data1,  
             ylab = "Propensity Scores by logistic regression", 
             xlab = "Treatment",auto.key = TRUE)

densityplot( ~ p_score_rf, groups = surgerytype, plot.points = FALSE, 
             xlim = c(0,1), lwd = 2,
             data = imputed_data1,  
             ylab = "Propensity Scores by random forest", 
             xlab = "Treatment",auto.key = TRUE)

densityplot( ~ p_score_gbm, groups = surgerytype, plot.points = FALSE, 
             xlim = c(0,1), lwd = 2,
             data = imputed_data1,  
             ylab = "Propensity Scores by generalized boosted modeling", 
             xlab = "Treatment",auto.key = TRUE)
## from the comparison we can conclude that 
## propensity scores get from logistic regression has most common support


### for all imputed data sets 
all_imputations_stacked <- data.frame() 
for (imp in 1:5) { temp <- cbind(all_imputations$imputations[[imp]],imputation = imp)
all_imputations_stacked = rbind(all_imputations_stacked, temp) }
all_imputations_stacked$surgerytype <- factor(all_imputations_stacked$surgerytype)
all_imputations_stacked$imputation <- factor(all_imputations_stacked$imputation,
                                           labels=paste("Imputation",1:5))

densityplot( ~ p_scores | imputation, data = all_imputations_stacked, 
             plot.points = FALSE, lwd = 2,
             groups = surgerytype, xlab = "Propensity Scores by logistic regression",
             auto.key = TRUE)
bwplot(p_scores ~ surgerytype | imputation, data = all_imputations_stacked, lwd = 2,
        ylab = "Propensity Scores by logistic regression", auto.key = TRUE)

#=================================Propensity scores weighting=======================#
#### weighting

## for single imputation data:
# for propensity scores from logistic regression 
# ATT weights
imputed_data1$weight_att <- with(imputed_data1, 
                                  ifelse(surgerytype == 1, 1, p_scores/(1 - p_scores)))

with(imputed_data1, by(weight_att, surgerytype, summary))

# ATE weights
imputed_data1$weight_ate <- with(imputed_data1, 
                                  ifelse(surgerytype == 1, 1/p_scores, 1/(1 - p_scores)))

with(imputed_data1, by (weight_ate, surgerytype, summary))

## balance table for ATT weight
balance_table1 <- bal.stat(imputed_data1, vars = cov_names_final,
                          treat.var = "surgerytype", w.all = imputed_data1$weight_att,
                          get.ks = FALSE, sampw = 1, estimand = "ATT", multinom = FALSE)
# I cut off test statistics from the balance table result 
# because it is not recommended: first because covariate balance is a property of 
# the sample, and hypothesis tests refer to the population,
# second, inferential measures depend on sample size, as our data set is big,
# the test statistics will always be significant.

round(balance_table1$results[, 1:5], 3)
# we can see that the tx.sd and ct.sd are very close, 
# and the absolute values of standard effect size are lower than 0.1,
# which indicating we have achieved adequate covariate balance.
summary(abs(balance_table1$results[, 5]))


## balance table for ATE weight
balance_table2 <- bal.stat(imputed_data1, vars = cov_names_final,
                          treat.var = "surgerytype", w.all = imputed_data1$weight_ate,
                          get.ks = FALSE, sampw = 1, estimand = "ATE", multinom = FALSE)

round(balance_table2$results[, 1:5], 3)
summary(abs(balance_table2$results[, 5]))
## for ps from random forest
## ATT weights
imputed_data1$weight_att_rf <- with(imputed_data1, 
                                  ifelse(surgerytype == 1, 1, p_score_rf/(1 - p_score_rf)))

with(imputed_data1, by(weight_att_rf, surgerytype, summary))

# ATE weights
imputed_data1$weight_ate_rf <- with(imputed_data1, 
                                  ifelse(surgerytype == 1, 1/p_score_rf, 1/(1 - p_score_rf)))

with(imputed_data1, by (weight_ate_rf, surgerytype, summary))

# balance table for ATT rf weights
balance_table3 <- bal.stat(imputed_data1, vars = cov_names_final,
                           treat.var = "surgerytype", w.all = imputed_data1$weight_att_rf,
                           get.ks = FALSE, sampw = 1, estimand = "ATT", multinom = FALSE)

round(balance_table3$results[, 1:5], 3)
summary(abs(balance_table3$results[, 5]))
# reach balance in ATT 

# balance table for ATE rf weights
balance_table4 <- bal.stat(imputed_data1, vars = cov_names_final,
                           treat.var = "surgerytype", w.all = imputed_data1$weight_ate_rf,
                           get.ks = FALSE, sampw = 1, estimand = "ATE", multinom = FALSE)

round(balance_table4$results[, 1:5], 3)
summary(abs(balance_table4$results[, 5]))
# failed to reach balance in ATE

### for ps from gbm
## ATT weights
imputed_data1$weight_att_gbm <- with(imputed_data1, 
                                     ifelse(surgerytype == 1, 1, p_score_gbm/(1 - p_score_gbm)))

with(imputed_data1, by(weight_att_gbm, surgerytype, summary))

# ATE weights
imputed_data1$weight_ate_gbm <- with(imputed_data1, 
                                     ifelse(surgerytype == 1, 1/p_score_gbm, 1/(1 - p_score_gbm)))

with(imputed_data1, by (weight_ate_gbm, surgerytype, summary))

# balance table for ATT gbm weights
balance_table5 <- bal.stat(imputed_data1, vars = cov_names_final,
                           treat.var = "surgerytype", w.all = imputed_data1$weight_att_gbm,
                           get.ks = FALSE, sampw = 1, estimand = "ATT", multinom = FALSE)

round(balance_table5$results[, 1:5], 3)
summary(abs(balance_table5$results[, 5]))
# reached balance in ATT 

# balance table for ATE gbm weights
balance_table6 <- bal.stat(imputed_data1, vars = cov_names_final,
                           treat.var = "surgerytype", w.all = imputed_data1$weight_ate_gbm,
                           get.ks = FALSE, sampw = 1, estimand = "ATE", multinom = FALSE)

round(balance_table6$results[, 1:5], 3)
summary(abs(balance_table6$results[, 5]))
# reached balance in ATE
# according to the standard effect size
# weighting using propensity scores get from logistic regression 
# achieved best performance in covariate balance
# so our final estimation of treatment effects will based on ps from logistic regression

### estimation of treatment effects
## ATT 
svy_design_att <- svydesign(ids = ~0, weights = imputed_data1$weight_att,
                        data = imputed_data1)

svy_design_boot_att <- as.svrepdesign(svy_design_att, type = c("bootstrap"), replicates = 2000)

# svyglm model 
model_att <- svyglm(died ~ surgerytype + severe + cognitivedecline + depression + cancer + 
                    autoimmune + transferedin + age + sex + yearseducation + 
                    bloodpressure + temperature + creatinelevels + sodiumlevels + 
                    urineweight + kg + nolifesupportorder + insurancetype + resp + 
                    infection + trauma + race + income + kg_na, svy_design_att, family = binomial())
summary(model_att)


outcome_att <- tidy(model_att) %>%
   filter(p.value < 0.05) %>%
   mutate(odd_ratio = exp(estimate))

outcome_att

# bootstrap model

model_att_boot <- svyglm(died ~ surgerytype + severe + cognitivedecline + depression + cancer + 
                            autoimmune + transferedin + age + sex + yearseducation + 
                            bloodpressure + temperature + creatinelevels + sodiumlevels + 
                            urineweight + kg + nolifesupportorder + insurancetype + resp + 
                            infection + trauma + race + income + kg_na, svy_design_boot_att,
                         family = binomial())
summary(model_att_boot)

outcome_att_boot <- tidy(model_att_boot) %>%
   filter(p.value < 0.05) %>%
   mutate(odd_ratio = exp(estimate))
outcome_att_boot

# we can see that the results are very close
## ATE 
svy_design_ate <- svydesign(ids = ~0, weights = imputed_data1$weight_ate,
                            data = imputed_data1)

svy_design_boot_ate <- as.svrepdesign(svy_design_ate, type = c("bootstrap"), replicates = 2000)

# svyglm model 
model_ate <- svyglm(died ~ surgerytype + severe + cognitivedecline + depression + cancer + 
                       autoimmune + transferedin + age + sex + yearseducation + 
                       bloodpressure + temperature + creatinelevels + sodiumlevels + 
                       urineweight + kg + nolifesupportorder + insurancetype + resp + 
                       infection + trauma + race + income + kg_na, svy_design_ate, family = binomial())
summary(model_ate)


outcome_ate <- tidy(model_ate) %>%
   filter(p.value < 0.05) %>%
   mutate(odd_ratio = exp(estimate))

outcome_ate

# bootstrap model

model_ate_boot <- svyglm(died ~ surgerytype + severe + cognitivedecline + depression + cancer + 
                            autoimmune + transferedin + age + sex + yearseducation + 
                            bloodpressure + temperature + creatinelevels + sodiumlevels + 
                            urineweight + kg + nolifesupportorder + insurancetype + resp + 
                            infection + trauma + race + income + kg_na, svy_design_boot_ate,
                         family = binomial())
summary(model_ate_boot)

outcome_ate_boot <- tidy(model_ate_boot) %>%
   filter(p.value < 0.05) %>%
   mutate(odd_ratio = exp(estimate))
outcome_ate_boot
# results are still very similar

### weighting with multiple imputed data sets
# weight ATT
all_imputations <- update(all_imputations, 
                          weigh_att = ifelse(surgerytype == 1, 1, p_scores/(1 - p_scores)))
# weight ATE 
all_imputations <- update(all_imputations, 
                          weigh_ate = ifelse(surgerytype == 1, 1/p_scores, 1/(1 - p_scores)))
# model with weight att
svy_design_mi_att <- svydesign(ids = ~ 0, weights = ~ weigh_att, data = all_imputations)
model_mi_att <- with(svy_design_mi_att, svyglm(died ~ surgerytype, family = binomial()))
result_model_mi_att <- MIcombine(model_mi_att)
summary(result_model_mi_att)

# model with weight ate
svy_design_mi_ate <- svydesign(ids = ~ 0, weights = ~ weigh_ate, data = all_imputations)
model_mi_ate <- with(svy_design_mi_ate, svyglm(died ~ surgerytype, family = binomial()))
result_model_mi_ate <- MIcombine(model_mi_ate)
summary(result_model_mi_ate)

# for both att and ate, we lost 0% percent data, 
# which means the final estimate and standard error are very similar
# to those obtained with a single imputed data set
# so we are using the single imputed data result as our final result

### sensitivity analysis

sens_att <- treatSens(formula = died ~ surgerytype + p_scores + 
                     I(p_scores ^ 2) + I(p_scores ^ 3), resp.family = binomial(),
                  trt.family = binomial(link = "probit"), grid.dim = c(5, 5), nsim = 20, 
                  weights = imputed_data1$weight_att, data = imputed_data1)
?treatSens
# I tried to build a sensitivity analysis, 
# but for now, the treatSens function can only deal with continuous ourcome variable
# so I gave it up.

#=========================Propensity score stratification====================#
### stratification

# Stratification of ATT with logistic regression propensity scores

# we need to coerce the treatment variable to numeric again
imputed_data1$surgerytype <- as.numeric(imputed_data1$surgerytype == 1)
stratification <- matchit(ps_formula_final, distance = imputed_data1$p_scores,
                          data = imputed_data1, method = "subclass", sub.by = "treat",
                          subclass = 5)

stratification

balance_stratification <- summary(stratification, standardize = TRUE)
strat_diff <- data.frame(balance_stratification$q.table[, 3, ])
summary(strat_diff)
# It can be concluded that adequate covariate balance was only achieved in subclass 4
# with the criterion that standardized mean differences should be less than 0.1
# we can use marginal mean weighting through stratification to mitigate this problem

# Stratification of ATT with random forest propensity scores
stratification_rf <- matchit(ps_formula_final, distance = imputed_data1$p_score_rf,
                          data = imputed_data1, method = "subclass", sub.by = "treat",
                          subclass = 5)

stratification_rf

balance_stratification_rf <- summary(stratification_rf, standardize = TRUE)
strat_diff_rf <- data.frame(balance_stratification_rf$q.table[, 3, ])
summary(strat_diff_rf)

# Well the performance is even worser compared with logistic regression propensity scores

# Stratification of ATT with generalize boosted modeling propensity scores:
stratification_gbm <- matchit(ps_formula_final, distance = imputed_data1$p_score_gbm,
                             data = imputed_data1, method = "subclass", sub.by = "treat",
                             subclass = 5)

stratification_gbm
# Covariate balance check
balance_stratification_gbm <- summary(stratification_gbm, standardize = TRUE)
strat_diff_gbm <- data.frame(balance_stratification_gbm$q.table[, 3, ])
summary(strat_diff_gbm)
# seems like the gbm method generates the best result in covariate balance
# we can just try to improve the subclass numbers to see if it can reach the loose
# criterion in all subclasses which is 0.25
stratification_gbm1 <- matchit(ps_formula_final, distance = imputed_data1$p_score_gbm,
                              data = imputed_data1, method = "subclass", sub.by = "treat",
                              subclass = 8)
balance_stratification_gbm1 <- summary(stratification_gbm1, standardize = TRUE)
strat_diff_gbm1 <- data.frame(balance_stratification_gbm1$q.table[, 3, ])
summary(strat_diff_gbm1)

# yes, we reached the cut off of 0.25 by using 8 subclasses
# we can still try with the marginal mean weighting through stratification 
# for better covariate balance

# marginal mean weighting through stratification for ATT 
stratum <- match.data(stratification)
stratum_rf <- match.data(stratification_rf)
stratum_gbm <- match.data(stratification_gbm)

# for logistic regression propensity scores
design <- svydesign(ids = ~0, data = stratum)

ntreat <- data.frame(table(stratum$subclass[stratum$surgerytype == 1]))
names(ntreat) <- c("subclass", "N.1s")
ncontrol <- data.frame(table(stratum$subclass[stratum$surgerytype == 0]))
names(ncontrol) <- c("subclass", "N.0s")

scounts <- merge(ntreat, ncontrol)
stratum <-merge(stratum, scounts)
propt <- svymean(~factor(surgerytype), design)

stratum$w <- with(stratum, ifelse(surgerytype == 1, 1, 
                                  stratum$N.1s * propt[1] / stratum$N.0s * propt[2]))

xtabs(~w + subclass, stratum)

# covariate balance check:
stratum$norm_weight <- stratum$w / mean(stratum$w)
norm_table <- bal.stat(stratum, estimand = "ATT", w.all = stratum$norm_weight,
                       vars = cov_names_final, sampw = 1, get.ks = FALSE,
                       treat.var = "surgerytype", multinom = FALSE)
round(norm_table$results[, 1:5], 3)
summary(abs(norm_table$results[, 5]))
# Now the result are perfect 


# for random forest propensity scores
design_rf <- svydesign(ids = ~0, data = stratum_rf)

ntreat_rf <- data.frame(table(stratum_rf$subclass[stratum_rf$surgerytype == 1]))
names(ntreat_rf) <- c("subclass", "N.1s")
ncontrol_rf <- data.frame(table(stratum_rf$subclass[stratum_rf$surgerytype == 0]))
names(ncontrol_rf) <- c("subclass", "N.0s")

scounts_rf <- merge(ntreat_rf, ncontrol_rf)
stratum_rf <-merge(stratum_rf, scounts_rf)
propt_rf <- svymean(~factor(surgerytype), design_rf)

stratum_rf$w <- with(stratum_rf, ifelse(surgerytype == 1, 1, 
                                        stratum_rf$N.1s * propt_rf[1] / stratum_rf$N.0s * propt_rf[2]))

xtabs(~w + subclass, stratum_rf)
# Covariate balance check
stratum_rf$rf_weight <- stratum_rf$w / mean(stratum_rf$w)
rf_table <- bal.stat(stratum_rf, estimand = "ATT", w.all = stratum_rf$rf_weight,
                         vars = cov_names_final, sampw = 1, get.ks = FALSE,
                         treat.var = "surgerytype", multinom = FALSE)

summary(abs(rf_table$results[, 5]))
# this one is still not good 

# for generalize boosted modeling propensity scores:
design_gbm <- svydesign(ids = ~0, data = stratum_gbm)

ntreat_gbm <- data.frame(table(stratum_gbm$subclass[stratum_gbm$surgerytype == 1]))
names(ntreat_gbm) <- c("subclass", "N.1s")
ncontrol_gbm <- data.frame(table(stratum_gbm$subclass[stratum_gbm$surgerytype == 0]))
names(ncontrol_gbm) <- c("subclass", "N.0s")

scounts_gbm <- merge(ntreat_gbm, ncontrol_gbm)
stratum_gbm <-merge(stratum_gbm, scounts_gbm)
propt_gbm <- svymean(~factor(surgerytype), design_gbm)

stratum_gbm$w <- with(stratum_gbm, 
                                  ifelse(surgerytype == 1, 1, 
                                         stratum_gbm$N.1s * propt_gbm[1] / stratum_gbm$N.0s * propt_gbm[2]))

xtabs(~w + subclass, stratum_gbm)

# Covariate balance check:
stratum_gbm$gbm_weight <- stratum_gbm$w / mean(stratum_gbm$w)
gbm_table <- bal.stat(stratum_gbm, estimand = "ATT", w.all = stratum_gbm$gbm_weight,
                     vars = cov_names_final, sampw = 1, get.ks = FALSE,
                     treat.var = "surgerytype", multinom = FALSE)

summary(abs(gbm_table$results[, 5]))

# the result is good and very close to logistic regression propensity scores
# as the mean of standardized effect size in logistic regression propensity scores 
# is slightly lower, we will choose logistic regression propensity scores to run our final model

### Estimation of treatment effect

final_design <- svydesign(ids = ~0, weights = stratum$norm_weight, data = stratum)
final_design_boot <- as.svrepdesign(final_design, type = c("bootstrap"), replicates = 2000)
model_boot <- svyglm(died ~ factor(surgerytype) + severe + cognitivedecline + depression + cancer + 
                            autoimmune + transferedin + age + sex + yearseducation + 
                            bloodpressure + temperature + creatinelevels + sodiumlevels + 
                            urineweight + kg + nolifesupportorder + insurancetype + resp + 
                            infection + trauma + race + income + kg_na, final_design_boot,
                         family = binomial())
summary(model_boot)
outcome_strat <- tidy(model_boot) %>%
   filter(p.value < 0.05) %>%
   mutate(odd_ratio = exp(estimate))

outcome_strat
### sensitivity analysis
sens_strat <- treatSens(formula = died ~ surgerytype + p_scores + 
                         I(p_scores ^ 2) + I(p_scores ^ 3), resp.family = binomial(),
                      trt.family = binomial(link = "probit"), grid.dim = c(5, 5), nsim = 20, 
                      weights = stratum$norm_weight, data = stratum)


#================================Propensity scores matching==========================#

