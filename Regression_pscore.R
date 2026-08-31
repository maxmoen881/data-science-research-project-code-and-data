#Regression Modelling for Partnership Score

#Install Packages and load libraries

library(knitr)
library(dplyr)
library(survival)
library(ggplot2)
library(tibble)
library(lubridate)
library(ggsurvfit)
library(gtsummary)
library(tidycmprsk)
library(lmtest)
library(condsurv)
library(tidyr)
library(ranger)
library(ggfortify)
library(moments)
library(tidyverse)
library(broom)
library(flexsurv)
library(SurvRegCensCov)
library(survminer)
library(sm)
library(splines)
library(survcomp)
library(pscl)
library(AER)
library(DHARMa)
library(e1071)
library(fitdistrplus)
library(MASS)
library(ggeffects)

#Starting with some EDA

#aggregate to partnership level

score_data <- clean_cricket_data %>%
  dplyr::select(match_id, innings, partnership_number,
                partnership_score, partnership_length, partnership_start,
                venue, batting_team, bowling_team,
                season, gender, run_rate = partnership_rr_start, required_run_rate = partnership_required_rr_start) %>%
  # Remove rows where partnership_score is NA
  dplyr::filter(!is.na(partnership_score)) %>%
  filter(season %in% c("2021", "2022", "2023", "2024"))

#Assess the distribution of scores

#summary statistics
summary(score_data$partnership_score)
table(score_data$partnership_score)
mean(score_data$partnership_score)
sd(score_data$partnership_score)
median(score_data$partnership_score)
IQR(score_data$partnership_score)

range(score_data$partnership_score)

quantile(score_data$partnership_score,
         probs = c(0,0.1,0.25,0.5,0.75,0.9,0.95,0.99,1))

#Basic summaries
cat("Total partnerships:", nrow(score_data), "\n")
cat("Zero partnerships:", sum(score_data$partnership_score == 0), "\n")
cat("Proportion zeros:", mean(score_data$partnership_score == 0), "\n\n")


#Visualisations

# Histogram
hist(score_data$partnership_score, 
     main = "Partnership Score Distribution", 
     xlab = "Runs", ylab = "Frequency", breaks = 30)

mean(score_data$partnership_score)
var(score_data$partnership_score)

#Boxplot by season
boxplot(partnership_score ~ season, data = score_data,
        main = "Partnership Scores by Season")

#Boxplot by gender
boxplot(partnership_score ~ gender, data = score_data,
        main = "Partnership Scores by Gender")

#checking the skewness
skewness(score_data$partnership_score)

#Check for overdispersion (for Poisson baseline)
#If var >> mean, overdispersion is present
cat("Mean:", mean(score_data$partnership_score), "\n")
cat("Variance:", var(score_data$partnership_score), "\n")
cat("Variance/Mean ratio:", var(score_data$partnership_score) / mean(score_data$partnership_score), "\n")

#therefore poisson is not appropriate.

#Trying different models

# Define model formulas
formula_null <- partnership_score ~ 1
formula_full <- partnership_score ~ innings + gender + venue + partnership_number + partnership_start + batting_team + bowling_team

#Poisson
poisson_null <- glm(formula_null, data = score_data, family = poisson)
poisson_full <- glm(formula_full, data = score_data, family = poisson)

#AIC
cat("\nNull Model AIC:", AIC(poisson_null), "\n")
cat("Full Model AIC:", AIC(poisson_full), "\n")

#Overdispersion test
dispersiontest(poisson_full)
#the AIC is much lower for the full model, yet there is overdispersion so I wont use poisson.

#Negative binomial
nb_null <- glm.nb(formula_null, data = score_data)
nb_full <- glm.nb(formula_full, data = score_data)

cat("\nNull Model AIC:", AIC(nb_null), "\n")
cat("Full Model AIC:", AIC(nb_full), "\n")
cat("Full Model Summary:\n")
summary(nb_full)

exp(coef(nb_full))

#zero-inflated negative binomial

#create zinb formulas
formula_null_zi <- partnership_score ~ 1 | 1
formula_full_zi <- partnership_score ~ 
  innings + venue + partnership_number + partnership_start + batting_team | 1

#having problems with optimisation error - probs because teams are nested within gender.
#i have removed gender/bowling team now as gender is nested within batting_team anyway.

zinb_null <- zeroinfl(formula_null_zi, data = score_data, dist = "negbin")
zinb_full <- zeroinfl(formula_full_zi, data = score_data, dist = "negbin")

cat("\nNull Model AIC:", AIC(zinb_null), "\n")
cat("Full Model AIC:", AIC(zinb_full), "\n")
cat("Full Model Summary:\n")
summary(zinb_full)

#hurdle negative binomial (zero-inflated)
hurdnb_null <- hurdle(formula_null_zi, data = score_data, dist = "negbin")
hurdnb_full <- hurdle(formula_full_zi, data = score_data, dist = "negbin")

cat("\nNull Model AIC:", AIC(hurdnb_null), "\n")
cat("Full Model AIC:", AIC(hurdnb_full), "\n")
cat("Full Model Summary:\n")
summary(hurdnb_full)

##covariate selection
#the full model currently has ~ innings + venue + partnership_number + partnership_start + batting_team
AIC(nb_full, zinb_full, hurdnb_full)
#compare the 3 full models first
cat("\nNB Model AIC:", AIC(nb_full), "\n")
cat("Zero Inflated Model AIC:", AIC(zinb_full), "\n")
cat("Hurdle Model AIC:", AIC(hurdnb_full), "\n")

BIC(nb_full, zinb_full, hurdnb_full)

cat("\nNB Model BIC:", BIC(nb_full), "\n")
cat("Zero Inflated Model BIC:", BIC(zinb_full), "\n")
cat("Hurdle Model BIC:", BIC(hurdnb_full), "\n")

#can clearly see that zero-inflated and hurdle are the best so focus on these

#zero-inflated covariate selection
#can see that partnership_number and partnership_start are the most significant variables
#so i will keep these in for all models apart from 1 model each removing each
#also need to try with just part_number and just part_start

zinb_null <- zeroinfl(formula_null_zi, data = score_data, dist = "negbin")
zinb_full <- zeroinfl(formula_full_zi, data = score_data, dist = "negbin")

lrtest(zinb_null, zinb_full)

zinb_1_1 <- zeroinfl(partnership_score ~ 
                       partnership_number | 1,
                     data = score_data, dist = "negbin")
zinb_1_2 <- zeroinfl(partnership_score ~ 
                       partnership_start | 1,
                     data = score_data, dist = "negbin")
zinb_1_3 <- zeroinfl(partnership_score ~ 
                       innings | 1,
                     data = score_data, dist = "negbin")
zinb_2_1 <- zeroinfl(partnership_score ~ 
                       partnership_number + partnership_start | 1,
                     data = score_data, dist = "negbin")
zinb_2_2 <- zeroinfl(partnership_score ~ 
                       partnership_number + innings | 1,
                     data = score_data, dist = "negbin")
zinb_2_3 <- zeroinfl(partnership_score ~ 
                       innings + partnership_start | 1,
                     data = score_data, dist = "negbin")
zinb_3_1 <- zeroinfl(partnership_score ~ 
                       partnership_number + partnership_start + innings | 1,
                     data = score_data, dist = "negbin")
zinb_3_2 <- zeroinfl(partnership_score ~ 
                       partnership_number + partnership_start + venue | 1,
                     data = score_data, dist = "negbin")
zinb_3_3 <- zeroinfl(partnership_score ~ 
                       partnership_number + partnership_start + batting_team | 1,
                     data = score_data, dist = "negbin")

zinb_4_1 <- zeroinfl(partnership_score ~ 
                       partnership_number + partnership_start + innings + venue | 1,
                     data = score_data, dist = "negbin")
zinb_4_2 <- zeroinfl(partnership_score ~ 
                       partnership_number + partnership_start + innings + batting_team | 1,
                     data = score_data, dist = "negbin")
zinb_4_3 <- zeroinfl(partnership_score ~ 
                       partnership_number + partnership_start + venue + batting_team | 1,
                     data = score_data, dist = "negbin")

zinb_no_partnum <- zeroinfl(partnership_score ~ 
                              partnership_start + innings + venue + batting_team | 1,
                            data = score_data, dist = "negbin")

zinb_no_partstart <- zeroinfl(partnership_score ~ 
                                partnership_number + innings + venue + batting_team | 1,
                              data = score_data, dist = "negbin")

#Create comparison table for ZINB
zinb_comparison <- data.frame(
  Model = c("Null model (intercept only)",
            "Partnership number",
            "Partnership start",
            "Innings",
            "Partnership only (number + starting ball)",
            "Partnership number + innings",
            "Partnership start + innings",
            "Partnership + innings",
            "Partnership + venue",
            "Partnership + batting team",
            "Full - batting team",
            "Full - venue",
            "Full - innings",
            "Full - part number",
            "Full - part start",
            "Full model"),
  AIC = c(AIC(zinb_null), AIC(zinb_1_1), AIC(zinb_1_2), AIC(zinb_1_3),
          AIC(zinb_2_1), AIC(zinb_2_2), AIC(zinb_2_3),
          AIC(zinb_3_1), AIC(zinb_3_2), AIC(zinb_3_3),
          AIC(zinb_4_1), AIC(zinb_4_2), AIC(zinb_4_3),
          AIC(zinb_no_partnum), AIC(zinb_no_partstart),
          AIC(zinb_full)),
  BIC = c(BIC(zinb_null), BIC(zinb_1_1), BIC(zinb_1_2), BIC(zinb_1_3),
          BIC(zinb_2_1), BIC(zinb_2_2), BIC(zinb_2_3),
          BIC(zinb_3_1), BIC(zinb_3_2), BIC(zinb_3_3),
          BIC(zinb_4_1), BIC(zinb_4_2), BIC(zinb_4_3),
          BIC(zinb_no_partnum), BIC(zinb_no_partstart),
          BIC(zinb_full)),
  nParams = c(0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 4, 4, 5)
)

print(zinb_comparison)

zinb_comparison <- zinb_comparison %>%
  arrange(AIC) %>%
  mutate(DeltaAIC = AIC - min(AIC),
         DeltaBIC = BIC - min(BIC),
         AICwt = round(exp(-DeltaAIC/2) / sum(exp(-DeltaAIC/2)), 4))

cat("ZERO-INFLATED NEGATIVE BINOMIAL - Variable Selection:\n")
print(zinb_comparison)

#can clearly see that there are some outstanding models.
#part num + part start
#OR part num + part start + innings
#as a result i need to test with innings too

lrt_zinb <- lrtest(zinb_2_1, zinb_3_1)
print(lrt_zinb)

#doing likelihood ratio test between the two suggests just stick with the two predictors.

#Final chosen model is zinb_2_1. Just part start and number.
chosen_zinb <- zinb_2_1

#Now time to do the same process but for hurdle nb to compare
hurdnb_null <- hurdle(formula_null_zi, data = score_data, dist = "negbin")
hurdnb_full <- hurdle(formula_full_zi, data = score_data, dist = "negbin")
hurdnb_1_1 <- hurdle(partnership_score ~ 
                       partnership_number | 1,
                     data = score_data, dist = "negbin")
hurdnb_1_2 <- hurdle(partnership_score ~ 
                       partnership_start | 1,
                     data = score_data, dist = "negbin")
hurdnb_2 <- hurdle(partnership_score ~ 
                     partnership_number + partnership_start | 1,
                   data = score_data, dist = "negbin")
hurdnb_3_1 <- hurdle(partnership_score ~ 
                       partnership_number + partnership_start + innings | 1,
                     data = score_data, dist = "negbin")
hurdnb_3_2 <- hurdle(partnership_score ~ 
                       partnership_number + partnership_start + venue | 1,
                     data = score_data, dist = "negbin")
hurdnb_3_3 <- hurdle(partnership_score ~ 
                       partnership_number + partnership_start + batting_team | 1,
                     data = score_data, dist = "negbin")

hurdnb_4_1 <- hurdle(partnership_score ~ 
                       partnership_number + partnership_start + innings + venue | 1,
                     data = score_data, dist = "negbin")
hurdnb_4_2 <- hurdle(partnership_score ~ 
                       partnership_number + partnership_start + innings + batting_team | 1,
                     data = score_data, dist = "negbin")
hurdnb_4_3 <- hurdle(partnership_score ~ 
                       partnership_number + partnership_start + venue + batting_team | 1,
                     data = score_data, dist = "negbin")

hurdnb_no_partnum <- hurdle(partnership_score ~ 
                              partnership_start + innings + venue + batting_team | 1,
                            data = score_data, dist = "negbin")

hurdnb_no_partstart <- hurdle(partnership_score ~ 
                                partnership_number + innings + venue + batting_team | 1,
                              data = score_data, dist = "negbin")

#Create comparison table for hurdnb
hurdnb_comparison <- data.frame(
  Model = c("Null model (intercept only)",
            "Partnership number",
            "Partnership start",
            "Partnership only (number + starting ball)",
            "Partnership + innings",
            "Partnership + venue",
            "Partnership + batting team",
            "Full - batting team",
            "Full - venue",
            "Full - innings",
            "Full - part number",
            "Full - part start",
            "Full model"),
  AIC = c(AIC(hurdnb_null), AIC(hurdnb_1_1), AIC(hurdnb_1_2), 
          AIC(hurdnb_2), AIC(hurdnb_3_1), AIC(hurdnb_3_2), AIC(hurdnb_3_3),
          AIC(hurdnb_4_1), AIC(hurdnb_4_2), AIC(hurdnb_4_3),
          AIC(hurdnb_no_partnum), AIC(hurdnb_no_partstart),
          AIC(hurdnb_full)),
  BIC = c(BIC(hurdnb_null), BIC(hurdnb_1_1), BIC(hurdnb_1_2),
          BIC(hurdnb_2), BIC(hurdnb_3_1), BIC(hurdnb_3_2), BIC(hurdnb_3_3),
          BIC(hurdnb_4_1), BIC(hurdnb_4_2), BIC(hurdnb_4_3),
          BIC(hurdnb_no_partnum), BIC(hurdnb_no_partstart),
          BIC(hurdnb_full)),
  nParams = c(0, 1, 1, 2, 3, 3, 3, 4, 4, 4, 4, 4, 5)
)

print(hurdnb_comparison)

hurdnb_comparison <- hurdnb_comparison %>%
  arrange(AIC) %>%
  mutate(DeltaAIC = AIC - min(AIC),
         DeltaBIC = BIC - min(BIC),
         AICwt = round(exp(-DeltaAIC/2) / sum(exp(-DeltaAIC/2)), 4))

cat("HURDLE NEGATIVE BINOMIAL - Variable Selection:\n")
print(hurdnb_comparison)

#can see same predictor choices as for zero-inflated

#likelihood ratio test again
lrt_hurdnb <- lrtest(hurdnb_2, hurdnb_3_1)

#not significant so therefore choose partnum + partstart as model.
chosen_hurdlenb <- hurdnb_2

#now i need to compare the two chosen models
chosen_model_comparison <- data.frame(
  Model = c("ZINB - Best", "Hurdle NB - Best"),
  Formula = c("Partnership only (number + start)", "Partnership only (number + start)"),
  AIC = c(AIC(chosen_zinb), AIC(chosen_hurdlenb)),
  BIC = c(BIC(chosen_zinb), BIC(chosen_hurdlenb)),
  LogLik = c(logLik(chosen_zinb), logLik(chosen_hurdlenb))
)

chosen_model_comparison <- chosen_model_comparison %>%
  arrange(AIC) %>%
  mutate(DeltaAIC = AIC - min(AIC),
         AICwt = round(exp(-DeltaAIC/2) / sum(exp(-DeltaAIC/2)), 4))

print(chosen_model_comparison)

#this will give us overall comparison

#Identify overall best model
if(chosen_model_comparison$AIC[1] == AIC(chosen_zinb)) {
  best_model <- chosen_zinb
  best_model_name <- paste("ZINB -", "Partnership only (number + start)")
} else {
  best_model <- chosen_hurdlenb
  best_model_name <- paste("Hurdle NB -", "Partnership only (number + start)")
}

cat("\n*** BEST MODEL (by AIC):", best_model_name, "***\n")
cat("AIC:", round(AIC(best_model), 2), "\n")
cat("\nModel Summary:\n")
print(summary(best_model))

score_model_final <- zeroinfl(partnership_score ~ ns(partnership_start, df = 4) + partnership_number | 1,
                              data = score_data, dist = "negbin")
#SO we have found our BEST MODEL - Zero-inflated negative binomial, with 2 predictors -
#partnership_number, and partnership_start

#spline comparison

zinb_spline_start <- zeroinfl(partnership_score ~ 
                                partnership_number + ns(partnership_start, df = 3) | 1,
                              data = score_data, dist = "negbin")

zinb_spline_number <- zeroinfl(partnership_score ~ 
                                 ns(partnership_number, df = 3) + partnership_start | 1,
                               data = score_data, dist = "negbin")

zinb_splines <- zeroinfl(partnership_score ~ 
                           ns(partnership_number, df = 3) + ns(partnership_start, df = 3) | 1,
                         data = score_data, dist = "negbin")


AIC(zinb_spline_start, zinb_spline_number, zinb_splines)
BIC(zinb_spline_start, zinb_spline_number, zinb_splines)

#spline for partnership start only is the best
#test different df
zinb_spline2 <- zeroinfl(partnership_score ~ 
                           partnership_number + ns(partnership_start, df = 2) | 1,
                         data = score_data, dist = "negbin")

zinb_spline3 <- zeroinfl(partnership_score ~ 
                           partnership_number + ns(partnership_start, df = 3) | 1,
                         data = score_data, dist = "negbin")

zinb_spline4 <- zeroinfl(partnership_score ~ 
                           partnership_number + ns(partnership_start, df = 4) | 1,
                         data = score_data, dist = "negbin")

zinb_spline5 <- zeroinfl(partnership_score ~ 
                           partnership_number + ns(partnership_start, df = 5) | 1,
                         data = score_data, dist = "negbin")

zinb_spline6 <- zeroinfl(partnership_score ~ 
                           partnership_number + ns(partnership_start, df = 6) | 1,
                         data = score_data, dist = "negbin")

AIC(zinb_spline2, zinb_spline3, zinb_spline4, zinb_spline5, zinb_spline6)
BIC(zinb_spline2, zinb_spline3, zinb_spline4, zinb_spline5, zinb_spline6)

#spline df = 4 is the best.

#testing the logistic model for the zero-inflated component
#Create binary outcome: 0 if partnership_score == 0, 1 if > 0
score_data$score_binary <- as.numeric(score_data$partnership_score > 0)

#Fit logistic regression: does partnership_start predict zero vs non-zero
log_model_ps <- glm(score_binary ~ partnership_start, 
                    data = score_data, 
                    family = binomial)
summary(log_model_ps)

#does partnership_number predict zero vs non-zero
log_model_pn <- glm(score_binary ~ partnership_number, 
                    data = score_data, 
                    family = binomial)
summary(log_model_pn)

#do both predict zero vs non-zero
log_model_both <- glm(score_binary ~ partnership_start + partnership_number, 
                      data = score_data, 
                      family = binomial)
summary(log_model_both)

#Are partnership_start and partnership_number correlated?
cor.test(score_data$partnership_start, score_data$partnership_number)

#Visualise it
plot(score_data$partnership_start, score_data$partnership_number)

#Current model (intercept only in zero component)
current_intercept <- zeroinfl(partnership_score ~ ns(partnership_start, df = 4) + partnership_number | 1,
                              data = score_data, dist = "negbin")

#Using partnership_start in zero component
current_pstart <- zeroinfl(partnership_score ~ ns(partnership_start, df = 4) + partnership_number | partnership_start,
                           data = score_data, dist = "negbin")

current_pstart_spline <- zeroinfl(partnership_score ~ ns(partnership_start, df=4) + partnership_number | ns(partnership_start, df=4),
                                  data = score_data, dist = "negbin")

#Using partnership_number in zero component
current_pnum <- zeroinfl(partnership_score ~ ns(partnership_start, df = 4) + partnership_number | partnership_number,
                         data = score_data, dist = "negbin")

#Using both in zero component
current_both <- zeroinfl(partnership_score ~ ns(partnership_start, df = 4) + partnership_number | partnership_start + partnership_number,
                         data = score_data, dist = "negbin")

current_innings <- zeroinfl(partnership_score ~ ns(partnership_start, df = 4) + partnership_number | innings,
                            data = score_data, dist = "negbin")

current_innings_pstart <- zeroinfl(partnership_score ~ ns(partnership_start, df = 4) + partnership_number | innings + partnership_start,
                                   data = score_data, dist = "negbin")

current_innings_pstart_pnum <- zeroinfl(partnership_score ~ ns(partnership_start, df = 4) + partnership_number | innings + partnership_start + partnership_number,
                                        data = score_data, dist = "negbin")

current_bteam <- zeroinfl(partnership_score ~ ns(partnership_start, df = 4) + partnership_number | batting_team,
                          data = score_data, dist = "negbin")

current_venue <- zeroinfl(partnership_score ~ ns(partnership_start, df = 4) + partnership_number | venue,
                          data = score_data, dist = "negbin")

#Compare ll models
AIC(current_intercept, current_pstart, current_pstart_spline, current_pnum, current_both, current_innings, current_innings_pstart, current_innings_pstart_pnum, current_bteam, current_venue)
BIC(current_intercept, current_pstart, current_pstart_spline, current_pnum, current_both, current_innings, current_innings_pstart, current_innings_pstart_pnum, current_bteam, current_venue)

lrtest(current_intercept, current_innings)
#it seems as if adding a spline does not improve the model.
#i think to prevent overfitting we want to leave the logit model as intercept anyway.
#lets just test quickly between intercept and part start
#Model A: Spline in count + intercept in zero
model_A <- zeroinfl(partnership_score ~ ns(partnership_start, df=4) + partnership_number | 1,
                    data = score_data, dist = "negbin")

#Model B: Spline in count + partnership_start in zero
model_B <- zeroinfl(partnership_score ~ ns(partnership_start, df=4) + partnership_number | partnership_start,
                    data = score_data, dist = "negbin")

# Compare
AIC(model_A, model_B)
BIC(model_A, model_B)

lrtest(model_A, model_B)

#Intercept only is better.

#to plot observed vs predicted (training data)
summary(pred_zeroinfl_spline)
max(pred_zeroinfl_spline)
quantile(pred_zeroinfl_spline, c(.90, .95, .99, 1))
summary(score_data$partnership_score)
max(score_data$partnership_score)
quantile(score_data$partnership_score, c(.90, .95, .99, 1))

#try required_run_rate
score_innings2 <- score_data %>%
  filter(innings == 2)
m_current <- zeroinfl(partnership_score ~ ns(partnership_start, df=4) + partnership_number | 1,
                      data = score_innings2, dist = "negbin")

m_score_rrr <- zeroinfl(partnership_score ~ ns(partnership_start, df=4) + partnership_number + required_run_rate | 1,
                        data = score_innings2, dist = "negbin")

m_score_rrr_spline <- zeroinfl(partnership_score ~ ns(partnership_start, df=4) + partnership_number + ns(required_run_rate, df = 3) | 1,
                               data = score_innings2, dist = "negbin")

AIC(m_current, m_score_rrr, m_score_rrr_spline)
BIC(m_current, m_score_rrr, m_score_rrr_spline)

#rrr with spline seriously improves AIC of model, but not BIC
#lets try it in the logistic model

m_score_rrr_both <- zeroinfl(partnership_score ~ ns(partnership_start, df=4) + partnership_number + ns(required_run_rate, df = 3) | ns(required_run_rate, df = 3),
                             data = score_innings2, dist = "negbin")

#how about no spline in logistic
m_score_rrr_both_nospline <- zeroinfl(partnership_score ~ ns(partnership_start, df=5) + partnership_number + ns(required_run_rate, df = 3) | required_run_rate,
                                      data = score_innings2, dist = "negbin")

m_score_rrr_both_nosplines <- zeroinfl(partnership_score ~ ns(partnership_start, df=5) + partnership_number + required_run_rate | required_run_rate,
                                       data = score_innings2, dist = "negbin")

m_score_rrr_zerospline <- zeroinfl(partnership_score ~ ns(partnership_start, df=5) + partnership_number + required_run_rate | ns(required_run_rate, df = 3),
                                   data = score_innings2, dist = "negbin")

m_score_rrr_justzero <- zeroinfl(partnership_score ~ ns(partnership_start, df=5) + partnership_number  | required_run_rate,
                                 data = score_innings2, dist = "negbin")

m_score_rrr_justzerospline <- zeroinfl(partnership_score ~ ns(partnership_start, df=5) + partnership_number | ns(required_run_rate, df = 3),
                                       data = score_innings2, dist = "negbin")

AIC(m_current, m_score_rrr, m_score_rrr_spline, m_score_rrr_justzero, m_score_rrr_both_nospline, m_score_rrr_both)
BIC(m_current, m_score_rrr, m_score_rrr_spline, m_score_rrr_justzero, m_score_rrr_both_nospline, m_score_rrr_both)

#looks like RRR is best just in the count component, with a spline.
#need to try other splines now and compare to m_current.
#do not want it in the zero-inflated component.

m_score_rrr_spline2 <- zeroinfl(partnership_score ~ ns(partnership_start, df=4) + partnership_number + ns(required_run_rate, df = 2) | 1,
                                data = score_innings2, dist = "negbin")

m_score_rrr_spline3 <- zeroinfl(partnership_score ~ ns(partnership_start, df=4) + partnership_number + ns(required_run_rate, df = 3) | 1,
                                data = score_innings2, dist = "negbin")

m_score_rrr_spline4 <- zeroinfl(partnership_score ~ ns(partnership_start, df=4) + partnership_number + ns(required_run_rate, df = 4) | 1,
                                data = score_innings2, dist = "negbin")

m_score_rrr_spline5 <- zeroinfl(partnership_score ~ ns(partnership_start, df=4) + partnership_number + ns(required_run_rate, df = 5) | 1,
                                data = score_innings2, dist = "negbin")

AIC(m_current, m_score_rrr_spline2, m_score_rrr_spline3, m_score_rrr_spline4, m_score_rrr_spline5)
BIC(m_current, m_score_rrr_spline2, m_score_rrr_spline3, m_score_rrr_spline4, m_score_rrr_spline5)

#looking at m_current vs spline3
#do an lrtest
lrtest(m_current, m_score_rrr_spline3)

##m_score_rrr_spline3 is better as we are doing a predictive model.

#FINAL MODEL: m_score_rrr_spline3 <- 
#                   zeroinfl(partnership_score ~ 
#                   ns(partnership_start, df=4) + partnership_number + ns(required_run_rate, df = 3)
#                   | 1,
#                   data = score_innings2, dist = "negbin")


#Validation
#Internal Validation
score_model_final <- zeroinfl(partnership_score ~ ns(partnership_start, df = 4) + partnership_number | 1,
                              data = score_data, dist = "negbin")


innings2_zinbmodel <- zeroinfl(partnership_score ~ ns(partnership_start, df=4) + partnership_number + ns(required_run_rate, df = 3) | 1,
                               data = score_innings2, dist = "negbin")

#Residual diagnostics first
#for the best model
score_model_final
score_model_name <- "Partnership only (number + start)"

# Get residuals
residuals_pearson_score <- residuals(score_model_final, type = "pearson")

#Q-Q plot
qq_plot_pscore <- ggplot(data.frame(residuals = residuals_pearson_score), aes(sample = residuals)) +
  geom_qq(size = 1.8, alpha = 0.5, color = "black") +
  geom_qq_line(color = "red", linewidth = 1) +
  labs(
    title = "Q-Q Plot of Pearson Residuals for Partnership Score Model",
    subtitle = "Assessing normality assumption; right-skewed distribution expected for count data",
    x = "Theoretical Quantiles",
    y = "Sample Quantiles"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11, color = "gray30"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  )

ggsave("Plots/QQ_Plot_partnershipscore.png", 
       plot = qq_plot_pscore, width = 8, height = 6, dpi = 300)

#Scale-location
scale_location_pscore <- ggplot(data.frame(
  fitted = fitted(score_model_final),
  residuals = sqrt(abs(residuals_pearson_score))
), aes(x = fitted, y = residuals)) +
  geom_point(alpha = 0.5, size = 2, color = "black") +
  geom_smooth(method = "loess", se = TRUE, color = "blue", linewidth = 1.2) +
  labs(
    title = "Scale-Location Plot (Partnership Score Model)",
    x = "Fitted Values",
    y = "√|Pearson Residuals|"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  )

ggsave("Plots/Scale-location_partnershipscore.png", 
       plot = scale_location_pscore, width = 8, height = 6, dpi = 300)

#Histogram of residuals
hist_residuals_pscore <- ggplot(data.frame(residuals = residuals_pearson_score), aes(x = residuals)) +
  geom_histogram(bins = 20, fill = "gray50", color = "black", alpha = 0.7) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", linewidth = 1.2) +
  labs(
    title = "Histogram of Pearson Residuals (Partnership Score Model)",
    x = "Pearson Residuals",
    y = "Frequency"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  )

ggsave("Plots/Histogram_residuals_partnershipscore.png", 
       plot = hist_residuals_pscore, width = 8, height = 6, dpi = 300)

#Predicted vs Observed plot
pred_zeroinfl <- predict(score_model_final, type = "response")
mean(pred_zeroinfl)

pred_obs_pscore <- ggplot(plot_data, aes(x = observed, y = predicted)) +
  geom_point(alpha = 0.4, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", linewidth = 1.2) +
  labs(
    title = "Predicted vs Observed Partnership Score (Training Data)",
    x = "Observed Partnership Score (runs)",
    y = "Predicted Partnership Score (runs)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  )

ggsave("Plots/PredvsObs_training_partnershipscore.png", 
       plot = pred_obs_pscore, width = 8, height = 6, dpi = 300)

#Predicted scores represent expected values from the zero-inflated model.
#Despite the zero-inflation component, expected values remain >0 due to the non-zero tail of the conditional distribution.

obs_zeros       <- sum(score_data$partnership_score == 0)
pred_zeros_zinb <- sum(predict(score_model_final, type = "prob")[, "0"])
cat("Observed zeros:   ", obs_zeros, "\n")
cat("ZIP-predicted zeros:", round(pred_zeros_zinb, 1), "\n")

#hurdle as a comparison
pred_hurdle <- predict(hurdnb_full, type = "response")
plot(score_data$partnership_score,
     pred_hurdle,
     xlab = "Observed partnership score",
     ylab = "Predicted partnership score")

abline(0, 1)

par(mfrow = c(1, 1))

#final validation
cor(score_data$partnership_score, pred_zeroinfl)
MAE_score <- mean(abs(score_data$partnership_score - pred_zeroinfl))
RMSE_score <- sqrt(mean((score_data$partnership_score - pred_zeroinfl)^2))

MAE_score
RMSE_score

#Heatmap: Predicted Score by RRR × Partnership Number

# Create grid
rrr_grid <- expand.grid(
  required_run_rate = seq(0, 20, by=1),
  partnership_number = 1:10,
  partnership_start = median(score_innings2$partnership_start, na.rm=TRUE)  # Hold constant
)

# Get predictions from regression model
rrr_grid$predicted_score <- predict(innings2_zinbmodel, newdata = rrr_grid, type = "response")

inn2_pscore_heatmap <- ggplot(rrr_grid, aes(x = required_run_rate, y = partnership_number, fill = predicted_score)) +
  geom_tile(color = "white", linewidth = 0.1) +
  scale_fill_viridis_c(name = "Predicted Score\n(runs)", option = "viridis", breaks = seq(8, 26, by=2)) +
  scale_y_continuous(breaks = 1:10) +
  scale_x_continuous(breaks = seq(0, 20, by=2)) +
  labs(
    title = "Predicted Partnership Score by Required Run Rate and Partnership Number",
    subtitle = "Innings 2 Chasing Context Only",
    x = "Required Run Rate (runs per over)",
    y = "Partnership Number"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11, color = "gray30"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10)
  )

ggsave("Plots/innings2_heatmap_partnershipscore.png", 
       plot = inn2_pscore_heatmap, width = 8, height = 6, dpi = 300)


#Validation (external) on 2025 data.
#validation dataset is currently named validation_cricket_data.

validation_score_data <- validation_cricket_data %>%
  dplyr::select(
    partnership_length,
    partnership_score,
    wicket,
    partnership_number,
    partnership_start,
    gender,
    innings,
    batting_team,
    bowling_team,
    venue,
    run_rate,
    required_run_rate
  ) %>%
  # Remove rows where partnership_score is NA
  dplyr::filter(!is.na(partnership_score), !is.na(wicket))

cat("Number of 2025 validation partnerships:", nrow(validation_score_data), "\n\n")
range(validation_score_data$partnership_score)
#this data is now ready to be analysed.

#get predictions
pred_mean_score2025 <- predict(score_model_final, newdata = validation_score_data, type = "response")

plot_data_score_test <- data.frame(
  observed = validation_score_data$partnership_score,
  predicted = pred_mean_score2025
)

#plot predictions vs observed
pred_obs_2025_pscore <- ggplot(plot_data_score_test, aes(x = observed, y = predicted)) +
  geom_point(alpha = 0.4, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", linewidth = 1.2) +
  labs(
    title = "Predicted vs Observed Partnership Score (2025 Test Data)",
    x = "Observed Partnership Score (runs)",
    y = "Predicted Partnership Score (runs)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  )

ggsave("Plots/PredvsObs_test_partnershipscore.png", 
       plot = pred_obs_2025_pscore, width = 8, height = 6, dpi = 300)

#Residual plot
residuals_val <- validation_score_data$partnership_score - pred_mean_score2025

residuals_val_test <- data.frame(
  predicted = pred_mean_score2025,
  residuals = validation_score_data$partnership_score - pred_mean_score2025
)

resid_pred_2025_pscore <- ggplot(residuals_val_test, aes(x = predicted, y = residuals)) +
  geom_point(alpha = 0.4, size = 2, color = "black") +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 1.2) +
  geom_smooth(method = "loess", se = TRUE, color = "blue", fill = "lightblue", alpha = 0.2, linewidth = 1.2) +
  labs(
    title = "Residuals vs Predicted Score (2025 Test Data)",
    x = "Predicted Partnership Score (runs)",
    y = "Residuals (Observed - Predicted)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  )

ggsave("Plots/ResidVsPred_test_partnershipscore.png", 
       plot = resid_pred_2025_pscore, width = 8, height = 6, dpi = 300)

#How many observed zeros are actually in the 2025 test data
obs_zeros_2025 <- sum(validation_cricket_data$partnership_score == 0)
cat("2025 observed zeros:", obs_zeros_2025, "\n")
pred_zeros_2025 <- sum(predict(score_model_final, newdata = validation_cricket_data, type = "prob")[, "0"])
cat("2025-predicted zeros:", round(pred_zeros_2025, 1), "\n")

#MAE/RMSE
mae_score2025 <- mean(abs(validation_score_data$partnership_score - pred_mean_score2025))
rmse_score2025 <- sqrt(mean((validation_score_data$partnership_score - pred_mean_score2025)^2))
mape_score2025 <- mean(abs((validation_score_data$partnership_score - pred_mean_score2025) / (validation_score_data$partnership_score + 1))) * 100

cat("Mean Absolute Error (MAE):", round(mae_score2025, 3), "\n")
cat("Root Mean Squared Error (RMSE):", round(rmse_score2025, 3), "\n")
cat("Mean Absolute Percentage Error (MAPE):", round(mape_score2025, 2), "%\n\n")

#Calibration by score ranges
validation_score_data$pred <- pred_mean_score2025
validation_score_data$bin <- cut(validation_score_data$partnership_score, 
                                 breaks = c(0, 5, 15, 30, 100))
calibration <- validation_score_data %>%
  group_by(bin) %>%
  summarise(mean_actual = mean(partnership_score),
            mean_pred = mean(pred),
            n = n(),
            .groups = 'drop')
print(calibration)


##innings 2 model validation
#get predictions
zeroinfl_model_inn2 <- zeroinfl(partnership_score ~ 
                                  ns(partnership_start, df=4) + partnership_number + ns(required_run_rate, df = 3)| 1,
                                data = score_innings2, dist = "negbin")
zeroinfl_model_inn2

validation_score_innings2 <- validation_score_data %>%
  filter(innings == 2)

nrow(validation_score_innings2)

pred_mean_innings2 <- predict(zeroinfl_model_inn2, newdata = validation_innings2, type = "response")
print(pred_mean_innings2)

mae_inn2 <- mean(abs(validation_innings2$partnership_score - pred_mean_innings2))
rmse_inn2 <- sqrt(mean((validation_innings2$partnership_score - pred_mean_innings2)^2))
mape_inn2 <- mean(abs((validation_innings2$partnership_score - pred_mean_innings2) / (validation_innings2$partnership_score + 1))) * 100

cat("Mean Absolute Error (MAE):", round(mae_inn2, 3), "\n")
cat("Root Mean Squared Error (RMSE):", round(rmse_inn2, 3), "\n")
cat("Mean Absolute Percentage Error (MAPE):", round(mape_inn2, 2), "%\n\n")

#Get predictions
survival_model_final <- flexsurvreg(Surv(partnership_length, wicket) ~ ns(partnership_start, df=3),
                                    data = parametric_cricket_data,
                                    dist = "weibull")
zeroinfl_model_final <- zeroinfl(partnership_score ~ ns(partnership_start, df=4) + partnership_number | 1,
                                 data = score_data, dist = "negbin")

survival_model_inn2 <- flexsurvreg(Surv(partnership_length, wicket) ~ ns(partnership_start, df=3) + ns(required_run_rate, df = 4),
                                   data = parametric_innings2,
                                   dist = "weibull")

zeroinfl_model_inn2 <- zeroinfl(partnership_score ~ ns(partnership_start, df=4) + partnership_number + ns(required_run_rate, df = 3) | 1,
                                data = score_innings2, dist = "negbin")


pred_grid$predicted_score <- predict(zeroinfl_model_final, newdata = pred_grid, type = "response")

#PARTNERSHIP SCORE vs REQUIRED RUN RATE

score_rrr <- data.frame(
  partnership_start = 0,
  partnership_number = 1,
  required_run_rate = seq(0, 20, by=1)
)

predicted_score <- predict(
  zeroinfl_model_inn2,
  newdata = score_rrr,
  type = "response"
)

score_rrr$predicted_score <- predict(zeroinfl_model_inn2, newdata = score_rrr, type = "response")

pscore_vs_RRR_inn2 <- ggplot(score_rrr, aes(x = required_run_rate, y = predicted_score)) +
  geom_line(linewidth = 1.5, color = "darkgreen") +
  geom_point(size = 4, color = "darkgreen") +
  labs(title = "Partnership Score vs Required Run Rate",
       x = "Required Run Rate (runs per over)",
       y = "Predicted Partnership Score (runs)") +
  theme_gray() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  )

ggsave("Plots/pscore_vs_RRR_inn2.png", 
       plot = pscore_vs_RRR_inn2, width = 8, height = 6, dpi = 300)


#Partnership score vs start, by RRR, faceted by Partnership Number
score_grid <- expand.grid(
  partnership_start = seq(0, 95, by=5),
  partnership_number = c(1, 3, 5, 7, 10),
  required_run_rate = c(3, 6, 10, 15)
)

#Get predictions
score_grid$predicted_score <- predict(zeroinfl_model_inn2, newdata = score_grid, type = "response")

faceted_plot_inn2scores <- ggplot(score_grid, aes(x = partnership_start, y = predicted_score, color = factor(required_run_rate))) +
  geom_line(linewidth = 1) +
  facet_wrap(~partnership_number, labeller = labeller(partnership_number = c("1" = "Partnership 1", "3" = "Partnership 3", "5" = "Partnership 5", "7" = "Partnership 7", "10" = "Partnership 10"))) +
  scale_x_continuous(breaks = seq(0, 95, by=10)) +
  scale_y_continuous(breaks = seq(0, 30, by=5)) +
  scale_color_manual(name = "Required Run Rate\n(runs/over)", 
                     values = c("3" = "darkgreen", "6" = "yellow", "10" = "orange", "15" = "red")) +
  labs(
    title = "Partnership Score vs Innings Progression Under Different Chase Pressures",
    x = "Partnership Start (ball number)",
    y = "Predicted Partnership Score (runs)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11),
    legend.position = "top",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 11, face = "bold")
  )

ggsave("Plots/Inn2Score_vs_RRR_Faceted.png", 
       plot = faceted_plot_inn2scores, width = 12, height = 8, dpi = 300)
