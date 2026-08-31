#Survival Analysis of Partnership Length

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

#lets do some survival analysis
#we want to work with 2021-2024 data
sa_cricket_data <- clean_cricket_data %>%
  filter(season %in% c("2021", "2022", "2023", "2024"))

survival_null <- survfit(partnerships ~ 1)
str(survival_null)

#basic kaplan-meier curve for all data - excluding labels etc.
survfit2(partnerships ~ 1) |> 
  ggsurvfit() +
  labs(
    x = "Balls",
    y = "Overall partnership survival probability"
  ) + 
  add_confidence_interval() +
  add_risktable()

#this shows a KM curve for length of each partnership

str(sa_cricket_data$gender)
levels(sa_cricket_data$gender)

#define survival object
sa_obj <- Surv(time = sa_cricket_data$partnership_length,
               event = sa_cricket_data$wicket)

#we already have the partnerships object no?

#null cox model
cox_null <- coxph(
  sa_obj ~ 1,
  data = sa_cricket_data
)

summary(cox_null)

##trying it with each predictor
cox_gender <- coxph(
  sa_obj ~ gender,
  data = sa_cricket_data
)

cox_innings <- coxph(
  sa_obj ~ innings,
  data = sa_cricket_data
)

cox_pstart <- coxph(
  sa_obj ~ partnership_start,
  data = sa_cricket_data
)

summary(cox_gender)
summary(cox_innings)
summary(cox_pstart)

#try a full model
cox_full <- coxph(
  sa_obj ~
    gender +
    innings +
    partnership_start +
    partnership_number +
    venue +
    batting_team +
    bowling_team,
  data = sa_cricket_data
)

summary(cox_full)

#i now have 5 cox models. lets check the PH assumption first.
cox.zph(cox_full)
cox.zph(cox_gender)
cox.zph(cox_innings)
cox.zph(cox_pstart)

plot(cox.zph(cox_full))
plot(cox.zph(cox_gender))
plot(cox.zph(cox_innings))
plot(cox.zph(cox_pstart))

print(sa_obj)

#this shows many violations of the PH assumption.

#could stratify the offending coariate 0 although it seems to be all of them mostly
#or could add a time-varying coefficient (tt())

#remove unused factor levels - e.g. teams from 2026 Hundred
sa_cricket_data <- sa_cricket_data %>%
  droplevels()
str(sa_cricket_data)




#Trying some parametric models

str(clean_cricket_data)

#data cleaning and preparation

#set data frame
parametric_cricket_data <- sa_cricket_data %>%
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
    run_rate = partnership_rr_start,
    required_run_rate = partnership_required_rr_start
  ) %>%
  filter(!is.na(partnership_length), !is.na(wicket))

parametric_cricket_data$gender <- as.factor(parametric_cricket_data$gender)

head(parametric_cricket_data)

#create survival object
#this tells us how long the partnership lasted and whether it ended or not
surv_obj <- Surv(time = parametric_cricket_data$partnership_length,
                 event = parametric_cricket_data$wicket)

#create Kaplan Meier curve for all data
km_fit <- survfit(surv_obj ~ 1, data = parametric_cricket_data)
km_data <- data.frame(
  time = km_fit$time,
  surv = km_fit$surv,
  lower = km_fit$lower,
  upper = km_fit$upper
)

km_gender <- survfit(surv_obj ~ gender, data = parametric_cricket_data)

Km_plot_bygender <- ggsurvplot(km_gender, data = parametric_cricket_data, conf.int = TRUE,
                               censor = FALSE,
                               size = 1, palette = c("Blue", "Red"),
                               conf.int.alpha = 0.15,
                               surv.scale = "percent",
                               xlab = "Partnership Length (balls)",
                               ylab = "Wicket Survival Probability",
                               title = "Kaplan-Meier: Overall Partnership Survival in The Hundred By Gender",
                               ggtheme = theme_bw() +
                                 theme(
                                   plot.title = element_text(
                                     face = "bold",
                                     size = 14
                                   )
                                 ),
                               pval = TRUE,
                               break.time.by = 10,
)

ggsave("Plots/Km_bygender.png", 
       plot = Km_plot_bygender$plot, width = 8, height = 6, dpi = 300)

#p = 0.071, therefore not statistically significant in terms of difference.

km_innings <- survfit(surv_obj ~ innings, data = parametric_cricket_data)

Km_plot_byinnings <- ggsurvplot(km_innings, data = parametric_cricket_data, conf.int = TRUE,
                                censor = FALSE,
                                size = 1, palette = c("darkgreen", "purple"),
                                conf.int.alpha = 0.15,
                                surv.scale = "percent",
                                xlab = "Partnership Length (balls)",
                                ylab = "Wicket Survival Probability",
                                title = "Kaplan-Meier: Overall Partnership Survival in The Hundred By Innings",
                                ggtheme = theme_bw() +
                                  theme(
                                    plot.title = element_text(
                                      face = "bold",
                                      size = 14
                                    )
                                  ),
                                pval = TRUE,
                                break.time.by = 10,
                                
)

ggsave("Plots/Km_byinnings.png", 
       plot = Km_plot_byinnings$plot, width = 8, height = 6, dpi = 300)

##this is for innings stratification

###parametric models
#fit the models with all the covariates - these are AFT models

weibull_aft_model <- flexsurvreg(surv_obj ~ partnership_number + partnership_start
                                 + innings + batting_team + venue,
                                 data = parametric_cricket_data, dist = "weibull")

exp_aft_model <- flexsurvreg(surv_obj ~ partnership_number + partnership_start
                             + innings + batting_team + venue,
                             data = parametric_cricket_data, dist = "exp")

lnorm_aft_model <- flexsurvreg(surv_obj ~ partnership_number + partnership_start
                               + innings + batting_team + venue,
                               data = parametric_cricket_data, dist = "lnorm")

llogis_aft_model <- flexsurvreg(surv_obj ~ partnership_number + partnership_start
                                + innings + batting_team + venue,
                                data = parametric_cricket_data, dist = "llogis")

#compare using AIC and BIC
model_comparison_aft <- data.frame(
  Model = c("Weibull", "Exponential", "Log-normal", "Log-logistic"),
  AIC = c(weibull_aft_model$AIC, exp_aft_model$AIC, lnorm_aft_model$AIC, llogis_aft_model$AIC),
  BIC = c(weibull_aft_model$BIC, exp_aft_model$BIC, lnorm_aft_model$BIC, llogis_aft_model$BIC)
)
model_comparison_aft <- model_comparison_aft[order(model_comparison_aft$AIC), ]
print(model_comparison_aft)

#this shows that the Weibull is the best fit with lowest AIC/BIC

#Trying some PH models to compare
weibull_ph_model <- flexsurvreg(surv_obj ~ partnership_number + partnership_start
                                + innings + batting_team + venue,
                                data = parametric_cricket_data, dist = "weibullPH")

gompertz_ph_model <- flexsurvreg(surv_obj ~ partnership_number + partnership_start
                                 + innings + batting_team + venue,
                                 data = parametric_cricket_data, dist = "gompertz")

#compare using AIC and BIC
model_comparison_ph <- data.frame(
  Model = c("Weibull AFT", "Weibull PH", "Gompertz"),
  AIC = c(weibull_aft_model$AIC, weibull_ph_model$AIC, gompertz_ph_model$AIC),
  BIC = c(weibull_aft_model$BIC, weibull_ph_model$BIC, gompertz_ph_model$BIC)
)
model_comparison_ph <- model_comparison_ph[order(model_comparison_ph$AIC), ]
print(model_comparison_ph)

#we can clearly see that Weibull AFT is the best here too - and is the same as Weibull PH
#bad idea to compare with cox in this comparison, as they are based on different likelihoods

#testing likelihood ratio - only want to do this between nested models (i.e. for covariate selection)
lrtest(lnorm_aft_model, weibull_aft_model)
#weibull is far better fit

#Now I have decided that the Weibull AFT model is the best - need to do covariate selection

#Covariate selection
#starting with full model and null model
weibull_full <- flexsurvreg(surv_obj ~ partnership_number + partnership_start
                            + innings + batting_team + venue,
                            data = parametric_cricket_data, dist = "weibull")

weibull_null <- flexsurvreg(surv_obj ~ 1,
                            data = parametric_cricket_data,
                            dist = "weibull")

AIC(weibull_full, weibull_null)
BIC(weibull_null)
BIC(weibull_full)

#excellent - the full model is considerably better than the null model on AIC/BIC.

#check the summary of the full model
summary(weibull_full)

tidy(weibull_full, conf.int = TRUE)

#fitting different combinations of covariates

#four covariates
weibull_novenue <- flexsurvreg(surv_obj ~ partnership_number + partnership_start
                               + innings + batting_team,
                               data = parametric_cricket_data, dist = "weibull")
weibull_noteam <- flexsurvreg(surv_obj ~ partnership_number + partnership_start
                              + innings + venue,
                              data = parametric_cricket_data, dist = "weibull")
weibull_noinn <- flexsurvreg(surv_obj ~ partnership_number + partnership_start
                             + batting_team + venue,
                             data = parametric_cricket_data, dist = "weibull")
weibull_nopstart <- flexsurvreg(surv_obj ~ partnership_number
                                + innings + batting_team + venue,
                                data = parametric_cricket_data, dist = "weibull")
weibull_nopnum <- flexsurvreg(surv_obj ~ partnership_start
                              + innings + batting_team + venue,
                              data = parametric_cricket_data, dist = "weibull")

AIC(weibull_full,
    weibull_novenue,
    weibull_noteam,
    weibull_noinn,
    weibull_nopstart,
    weibull_nopnum
)

#the AIC is similar for all of these. the one without batting_team is the best, then all the others are similar.
#without partnership_start is definitely the worst fit.

#lets use the one without batting_team and remove one more variable
weibull_nobat_1 <- flexsurvreg(surv_obj ~ partnership_number + partnership_start
                               + innings + venue,
                               data = parametric_cricket_data, dist = "weibull")
weibull_nobat_2 <- flexsurvreg(surv_obj ~ partnership_number + partnership_start
                               + innings,
                               data = parametric_cricket_data, dist = "weibull")
weibull_nobat_3 <- flexsurvreg(surv_obj ~ partnership_number + partnership_start
                               + venue,
                               data = parametric_cricket_data, dist = "weibull")
weibull_nobat_4 <- flexsurvreg(surv_obj ~ partnership_number
                               + innings + venue,
                               data = parametric_cricket_data, dist = "weibull")
weibull_nobat_5 <- flexsurvreg(surv_obj ~ partnership_start
                               + innings + venue,
                               data = parametric_cricket_data, dist = "weibull")

AIC(weibull_noteam,
    weibull_nobat_1,
    weibull_nobat_2,
    weibull_nobat_3,
    weibull_nobat_4,
    weibull_nobat_5
)


lrtest(weibull_nobat_1, weibull_nobat_2)

#at this point we can see that it is slightly lower AIC without venue or innings. But it is very similar.
#clearly removing partnership_start makes the model much worse

#lets do a systematic loop to test every model combination. this will give us the best AIC.
covariates <- c(
  "partnership_start",
  "partnership_number",
  "innings",
  "venue",
  "batting_team"
)

results <- data.frame(
  Model = character(),
  Formula = character(),
  Parameters = integer(),
  LogLik = numeric(),
  AIC = numeric(),
  BIC = numeric(),
  stringsAsFactors = FALSE
)

# Loop over model sizes
for (k in 1:length(covariates)) {
  
  # All combinations of size k
  combs <- combn(covariates, k, simplify = FALSE)
  
  for (vars in combs) {
    
    formula_text <- paste(
      "Surv(partnership_length, wicket) ~",
      paste(vars, collapse = " + ")
    )
    
    formula <- as.formula(formula_text)
    
    fit <- flexsurvreg(
      formula,
      data = parametric_cricket_data,
      dist = "weibull"
    )
    
    results <- rbind(
      results,
      data.frame(
        Model = paste0("M", nrow(results) + 1),
        Formula = formula_text,
        Parameters = length(fit$res.t),
        LogLik = fit$loglik,
        AIC = fit$AIC,
        BIC = fit$BIC
      )
    )
  }
}

results <- results[order(results$AIC), ]
results$DeltaAIC <- results$AIC - min(results$AIC)

results

#from this table we can see 8 good models:
#partnership_start
#partnership_start + partnership_number
#partnership_start + venue
#partnership_start + innings
#partnership_start + partnership_number + innings
#partnership_start + partnership_number + venue
#partnership_start + innings + venue
#partnership_start + partnership_number + innings + venue
model_partstart <- flexsurvreg(surv_obj ~ partnership_start,
                               data = parametric_cricket_data, dist = "weibull")
model_partstart_num <- flexsurvreg(surv_obj ~ partnership_start + partnership_number,
                                   data = parametric_cricket_data, dist = "weibull")
model_partstart_num_venue <- flexsurvreg(surv_obj ~ partnership_start + partnership_number + venue,
                                         data = parametric_cricket_data, dist = "weibull")
model_partstart_num_inn_venue <- flexsurvreg(surv_obj ~ partnership_start + partnership_number + innings + venue,
                                             data = parametric_cricket_data, dist = "weibull")

lrtest(model_partstart, model_partstart_num)

#similar AIC values. just partnership_start still just edges it, by doing lrtests.


#compare to spline model

model_partstart_spline <- flexsurvreg(
  surv_obj ~ ns(partnership_start, df = 3),
  data = parametric_cricket_data,
  dist = "weibull"
)

lrtest(model_partstart, model_partstart_spline)
lrtest(model_partstart_spline, model_partstart)
#higher loglik is better - spline is significantly better

AIC(model_partstart, model_partstart_spline)
BIC(model_partstart)
BIC(model_partstart_spline)
#spline is far better.
#lets compare: partstart and partnum with splines for both, splines for start and splines for num

model_partstart_partnum_splines <- flexsurvreg(
  surv_obj ~ ns(partnership_start, df = 3) + ns(partnership_number, df = 3),
  data = parametric_cricket_data,
  dist = "weibull"
)

model_partstart_spline_partnum <- flexsurvreg(
  surv_obj ~ ns(partnership_start, df = 3) + partnership_number,
  data = parametric_cricket_data,
  dist = "weibull"
)

model_partstart_partnum_spline <- flexsurvreg(
  surv_obj ~ ns(partnership_number, df = 3) + partnership_start,
  data = parametric_cricket_data,
  dist = "weibull"
)

AIC(model_partstart_spline, model_partstart_partnum_splines, model_partstart_spline_partnum, model_partstart_partnum_spline)

#it looks like the spline model with just partnership_start is the best
#lets just compare with different numbers of splines

model_spline2 <- flexsurvreg(
  surv_obj ~ ns(partnership_start, df = 2),
  data = parametric_cricket_data,
  dist = "weibull"
)

model_spline3 <- flexsurvreg(
  surv_obj ~ ns(partnership_start, df = 3),
  data = parametric_cricket_data,
  dist = "weibull"
)

model_spline4 <- flexsurvreg(
  surv_obj ~ ns(partnership_start, df = 4),
  data = parametric_cricket_data,
  dist = "weibull"
)

model_spline5 <- flexsurvreg(
  surv_obj ~ ns(partnership_start, df = 5),
  data = parametric_cricket_data,
  dist = "weibull"
)

AIC(model_spline2, model_spline3, model_spline4, model_spline5)
BIC(model_spline2)
BIC(model_spline3)
BIC(model_spline4)
BIC(model_spline5)

lrtest(model_spline3, model_spline4)

#AIC is better for 4 splines, but only slightly. BIC better for 3 and p-value = 0.051 so not significant (just) so we choose 3 to prevent overfitting.


##INNINGS 2 MODEL
#lets try required run rate for an innings 2 model
#Filter to innings 2 only
parametric_innings2 <- parametric_cricket_data %>%
  filter(innings == 2)

nrow(parametric_innings2)

#set innings 2 survival object
surv_obj_innings2 <- Surv(time = parametric_innings2$partnership_length,
                          event = parametric_innings2$wicket)

#Fit models with required_run_rate
#Model 1: partnership_start only (baseline)
model_inn2_ps <- flexsurvreg(surv_obj_innings2 ~ ns(partnership_start, df=3),
                             data = parametric_innings2,
                             dist = "weibull")

#Model 2: partnership_start + required_run_rate
model_inn2_ps_rrr <- flexsurvreg(surv_obj_innings2 ~ ns(partnership_start, df=3)  + required_run_rate,
                                 data = parametric_innings2,
                                 dist = "weibull")

#Model 3: partnership_start + partnership_number + required_run_rate
model_inn2_both <- flexsurvreg(surv_obj_innings2 ~ ns(partnership_start, df=3) + partnership_number + required_run_rate,
                               data = parametric_innings2,
                               dist = "weibull")


#Model 4: all predictors + rrr
model_inn2_full <- flexsurvreg(surv_obj_innings2 ~ ns(partnership_start, df=3) + partnership_number + venue + batting_team + required_run_rate,
                               data = parametric_innings2,
                               dist = "weibull")
#cant include innings as a predictor here.

#Compare
AIC(model_inn2_ps, model_inn2_ps_rrr, model_inn2_both, model_inn2_full)
BIC(model_inn2_ps)
BIC(model_inn2_ps_rrr)
BIC(model_inn2_both)
BIC(model_inn2_full)

summary(model_inn2_full)

lrtest(model_inn2_both, model_inn2_ps_rrr)

#so adding required run rate does slightly improve the model
#As a result lets try with a spline.
#Model 1: partnership_start only (baseline)
model_inn2_ps <- flexsurvreg(surv_obj_innings2 ~ ns(partnership_start, df=3),
                             data = parametric_innings2,
                             dist = "weibull")

#Model 2: partnership_start + required_run_rate(spline)
model_inn2_ps_rrr_spline <- flexsurvreg(surv_obj_innings2 ~ ns(partnership_start, df=3)  + ns(required_run_rate, df = 3),
                                        data = parametric_innings2,
                                        dist = "weibull")

#Model 3: partnership_start + partnership_number + required_run_rate (spline)
model_inn2_both_spline <- flexsurvreg(surv_obj_innings2 ~ ns(partnership_start, df=3) + partnership_number + ns(required_run_rate, df =3),
                                      data = parametric_innings2,
                                      dist = "weibull")

AIC(model_inn2_ps, model_inn2_ps_rrr, model_inn2_both, model_inn2_full, model_inn2_ps_rrr_spline, model_inn2_both_spline)
BIC(model_inn2_ps)
BIC(model_inn2_ps_rrr)
BIC(model_inn2_both)
BIC(model_inn2_full)
BIC(model_inn2_ps_rrr_spline)
BIC(model_inn2_both_spline)

#we now have a second innings model - part_start(spline) + required_run_rate(spline)

#lets check which spline is best
model_2spline <- flexsurvreg(Surv(partnership_length, wicket) ~ ns(partnership_start, df=3) + ns(required_run_rate, df = 2),
                             data = parametric_innings2,
                             dist = "weibull")

model_3spline <- flexsurvreg(Surv(partnership_length, wicket) ~ ns(partnership_start, df=3) + ns(required_run_rate, df = 3),
                             data = parametric_innings2,
                             dist = "weibull")

model_4spline <- flexsurvreg(Surv(partnership_length, wicket) ~ ns(partnership_start, df=3) + ns(required_run_rate, df = 4),
                             data = parametric_innings2,
                             dist = "weibull")

model_5spline <- flexsurvreg(Surv(partnership_length, wicket) ~ ns(partnership_start, df=3) + ns(required_run_rate, df = 5),
                             data = parametric_innings2,
                             dist = "weibull")

model_6spline <- flexsurvreg(Surv(partnership_length, wicket) ~ ns(partnership_start, df=3) + ns(required_run_rate, df = 6),
                             data = parametric_innings2,
                             dist = "weibull")


AIC(model_2spline, model_3spline, model_4spline, model_5spline, model_6spline)
BIC(model_2spline)
BIC(model_3spline)
BIC(model_4spline)
BIC(model_5spline)
BIC(model_6spline)



#Validation
setwd("C:/R_Projects/Dissertation")
dir.create("Plots", showWarnings = FALSE)

#Model choices
#Survival model on full data
survival_model_final <- flexsurvreg(Surv(partnership_length, wicket) ~ ns(partnership_start, df=3),
                                    data = parametric_cricket_data,
                                    dist = "weibull")

#coefficients
coef(survival_model_final)

#Survival model for innings 2
survival_model_inn2 <- flexsurvreg(Surv(partnership_length, wicket) ~ ns(partnership_start, df=3) + ns(required_run_rate, df = 4),
                                   data = parametric_innings2,
                                   dist = "weibull")

#Weibull Model vs Kaplan Meier

#Kaplan-Meier model
km_fit_clean <- survfit(surv_obj ~ 1, data = parametric_cricket_data)

#Extract KM estimates
km_data <- data.frame(
  time = km_fit_clean$time,
  survival = km_fit_clean$surv
)

#Weibull predictions
weibull_data <- summary(
  survival_model_final,
  type = "survival",
  ci = FALSE
)

weibull_data <- data.frame(
  time = weibull_data[[1]]$time,
  survival = weibull_data[[1]]$est
)

#Plot
figure_results4 <- ggplot() +
  geom_step(data = km_data, aes(x = time, y = survival, colour = "Kaplan-Meier"), linewidth = 1) +
  geom_line(data = weibull_data, aes(x = time, y = survival, colour = "Survival Model"), linewidth = 1) +
  scale_colour_manual(values = c("Kaplan-Meier" = "black", "Survival Model" = "red")) +
  labs(
    x = "Partnership Length (balls)",
    y = "Wicket Survival Probability",
    title = "Training Data: Observed Kaplan-Meier vs Fitted Weibull Model Survival Curve",
    colour = NULL
  ) +
  theme_gray() +
  theme(plot.title = element_text(
    face = "bold",
    margin = margin(10, 0, 10, 0),
    size = 14
  ),
  legend.position = "top"
  )

ggsave("Plots/KM_vs_Weibull_Training.png", 
       plot = figure_results4, width = 8, height = 6, dpi = 300)

###some validation plots - diagnostics
#final model is survival_model_final

#Using surv_reg for validation - have checked i can use either function for same model
survival_model_survreg <- survreg(Surv(partnership_length, wicket) ~ ns(partnership_start, df=3),
                                  data = parametric_cricket_data,
                                  dist = "weibull")

#Extract residuals
deviance_resid <- residuals(survival_model_survreg, type = "deviance")
response_resid <- residuals(survival_model_survreg, type = "response")

#Deviance residuals QQ plot
fitted_vals <- predict(survival_model_survreg, type = "response")
linear_pred <- predict(survival_model_survreg, type = "linear")

figure_results5 <- ggplot(data.frame(residuals = deviance_resid), aes(sample = residuals)) +
  geom_qq(size = 1.8, alpha = 0.5, color = "black") +
  geom_qq_line(color = "red", linewidth = 1) +
  labs(title = "Q-Q Plot of Deviance Residuals",
       subtitle = "Assessing normality assumption for Weibull model",
       x = "Theoretical Quantiles",
       y = "Sample Quantiles") +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11, color = "gray30"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  )

ggsave("Plots/QQ_Plot.png", 
       plot = figure_results5, width = 8, height = 6, dpi = 300)

#Log-log plot
time_points <- km_fit_clean$time[km_fit_clean$n.event > 0]
survival_probs <- km_fit_clean$surv[km_fit_clean$n.event > 0]

valid_idx <- (survival_probs > 0.01 & survival_probs < 0.99)
time_points <- time_points[valid_idx]
survival_probs <- survival_probs[valid_idx]

loglog_data <- data.frame(
  log_time = log(time_points),
  log_neg_log_surv = log(-log(survival_probs))
)

figure_resultsloglog <- ggplot(loglog_data, aes(x = log_time, y = log_neg_log_surv)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "lightblue", alpha = 0.2, linewidth = 1.2) +
  labs(title = "Log-Log Survival Plot: Weibull Linearity Check",
       x = "log(Partnership Length)",
       y = "log(-log(Survival))") +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11, color = "gray30"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  )

ggsave("Plots/LogLog_Plot.png", plot = figure_resultsloglog, width = 8, height = 6, dpi = 300)

#deviance residuals vs fitted values
resid_data <- data.frame(
  fitted = fitted_vals,
  residuals = deviance_resid
)

figure_results_devresid <- ggplot(resid_data, aes(x = fitted, y = residuals)) +
  geom_point(alpha = 0.5, size = 2, color = "black") +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 1) +
  geom_smooth(method = "loess", se = TRUE, color = "blue", alpha = 0.2, linewidth = 1.2) +
  labs(title = "Training Data: Deviance Residuals vs Fitted Values",
       x = "Fitted Partnership Length (balls)",
       y = "Deviance Residuals") +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  )

ggsave("Plots/Deviance_vs_Residuals_Plot.png", plot = figure_results_devresid, width = 8, height = 6, dpi = 300)


##External Validation

#I have my training data, now we need the test data.
#Filtering the 2025 data from the original full dataset
validation_cricket_data <- clean_cricket_data %>%
  filter(season == 2025)

#Calculating number of matches in test data
validation_cricket_data %>%
  group_by(gender) %>%
  summarise(
    total_matches = n_distinct(match_id)
  )

#Calculating number of matches in training data
sa_cricket_data %>%
  group_by(gender) %>%
  summarise(
    total_matches = n_distinct(match_id)
  )

#Number of innings in test data including percentages
validation_cricket_data %>%
  group_by(innings) %>%
  summarise(
    total_innings = n_distinct(match_id, innings)
  ) %>%
  mutate(
    percentage = total_innings / sum(total_innings) * 100
  )

#Number of innings in training data
sa_cricket_data %>%
  group_by(innings) %>%
  summarise(
    total_innings = n_distinct(match_id, innings)
  )

#we can see from this that there are 33 men's matches and 34 women's matches

#To check how many partnerships in total
training_matches <- nrow(sa_cricket_data)
testing_matches <- nrow(validation_cricket_data)

#Calculating percentages
total_matches <- training_matches + testing_matches
training_percentage <- training_matches / total_matches * 100
test_percentage <- testing_matches / total_matches * 100
print(training_percentage)
print(test_percentage)

#Censored partnerships calculate
validation_cricket_data %>%
  group_by(gender) %>%
  summarise(
    total_partnerships = n(),
    not_out_partnerships = sum(wicket == 0, na.rm = TRUE),
    wicket_ended_partnerships = sum(wicket == 1, na.rm = TRUE)
  )

#Ensure dataset is processed
validation_cricket_data <- validation_cricket_data %>%
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
    run_rate = partnership_rr_start,
    required_run_rate = partnership_required_rr_start
  ) %>%
  filter(!is.na(partnership_length), !is.na(wicket))

validation_cricket_data$gender <- as.factor(validation_cricket_data$gender)
#we now have our testing data

str(validation_cricket_data)

#set validation data
#survival_model_final is the model still

#Get predictions on 2025 data
pred_surv_time <- predict(survival_model_final, 
                          newdata = validation_cricket_data, 
                          type = "response")

#Concordance-index (measures rank correlation)

#Extract predicted partnership lengths
pred_time <- pred_surv_time$.pred_time

cindex_survival_model <- concordance(
  Surv(
    validation_cricket_data$partnership_length,
    validation_cricket_data$wicket
  ) ~ pred_time
)

cindex_survival_model
#Extract C-index
cindex_survival_model$concordance

#This gives me a C-index value of 0.609, moderate discriminatory ability.

#MAE and RMSE
mae_survival <- mean(abs(validation_cricket_data$partnership_length - pred_surv_time$.pred_time), na.rm = TRUE)
rmse_survival <- sqrt(mean((validation_cricket_data$partnership_length - pred_surv_time$.pred_time)^2, na.rm = TRUE))

cat("MAE (partnership length):", round(mae_survival, 3), "balls\n")
cat("RMSE (partnership length):", round(rmse_survival, 3), "balls\n")

#Calibration by bins
validation_cricket_data$pred <- pred_surv_time$.pred_time
validation_cricket_data$bin <- cut(validation_cricket_data$partnership_length, 
                                   breaks = quantile(validation_cricket_data$partnership_length, 
                                                     probs = seq(0, 1, 0.25), na.rm = TRUE))
calib <- validation_cricket_data %>%
  group_by(bin) %>%
  summarise(median_actual = median(partnership_length),
            median_pred = median(pred),
            n = n())
print(calib)

#Weibull vs KM

#Fit KM for 2025 data
km_test <- survfit(
  Surv(partnership_length, wicket) ~ 1,
  data = validation_cricket_data
)

km_test_data <- data.frame(
  time = km_test$time,
  survival = km_test$surv
)

#create prediction partnership lengths
times <- seq(0, max(validation_cricket_data$partnership_length), by = 1)

#using the chosen model
pred_testing <- summary(
  survival_model_final,
  newdata = validation_cricket_data,
  t = times,
  type = "survival"
)

pred_matrix <- sapply(pred_testing, function(x) x$est)
mean_surv <- rowMeans(pred_matrix)

weibull_test_data <- data.frame(
  time = times,
  survival = mean_surv
)

figure_results_test_WeibullKM <- ggplot() +
  geom_step(data = km_test_data, aes(x = time, y = survival, colour = "Kaplan-Meier"), linewidth = 1) +
  geom_line(data = weibull_test_data, aes(x = time, y = survival, colour = "Survival Model"), linewidth = 1) +
  scale_colour_manual(values = c("Kaplan-Meier" = "black", "Survival Model" = "red")) +
  labs(
    title = "Test Data (2025): Observed Kaplan-Meier vs Fitted Weibull Model",
    x = "Partnership Length (balls)",
    y = "Wicket Survival Probability",
    colour = NULL
  ) +
  theme_gray() +
  theme(plot.title = element_text(
    face = "bold",
    margin = margin(10, 0, 10, 0),
    size = 14
  ),
  legend.position = "top"
  )

ggsave("Plots/KM_vs_Weibull_Test.png", 
       plot = figure_results_test_WeibullKM, width = 8, height = 6, dpi = 300)


#lets validate the required_run_rate model

#Model 3: (3 spline)partnership_start + (4 spline)required_run_rate
innings2_survmodel <- flexsurvreg(Surv(partnership_length, wicket) ~ ns(partnership_start, df=3) + ns(required_run_rate, df = 4),
                                  data = parametric_innings2,
                                  dist = "weibull")

#C-index
# Filter validation data to innings 2
validation_innings2 <- validation_cricket_data %>%
  filter(innings == 2)
nrow(validation_innings2)

pred_surv_time_rrr <- predict(innings2_survmodel, 
                              newdata = validation_innings2, 
                              type = "response")

cindex_rrr <- concordance(
  Surv(
    validation_innings2$partnership_length,
    validation_innings2$wicket
  ) ~ pred_surv_time_rrr$.pred_time
)

cindex_rrr
cindex_rrr$concordance

#0.6336

#Calibration for second innings only
calib_innings2 <- validation_cricket_data %>%
  filter(innings == 2) %>%
  mutate(bin = cut(partnership_length, 
                   breaks = quantile(partnership_length, 
                                     probs = seq(0, 1, 0.25), na.rm = TRUE))) %>%
  group_by(bin) %>%
  summarise(median_actual = median(partnership_length),
            median_pred = median(pred),
            n = n())

calib_innings2


#Partnership Length vs Required Run Rate

length_rrr <- data.frame(
  partnership_start = 0,
  partnership_number = 1,
  required_run_rate = seq(0, 20, by=1)
)

pred_length <- predict(
  survival_model_inn2,
  newdata = length_rrr,
  type = "response"
)

length_rrr$predicted_length <- pred_length$.pred_time

plength_vs_RRR_inn2 <- ggplot(length_rrr, aes(x = required_run_rate, y = predicted_length)) +
  geom_line(linewidth = 1.5, color = "steelblue") +
  geom_point(size = 4, color = "steelblue") +
  labs(
    title = "Partnership Length vs Required Run Rate",
    x = "Required Run Rate (runs per over)",
    y = "Predicted Partnership Length (balls)"
  ) +
  theme_gray() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  )

ggsave("Plots/plength_vs_RRR_inn2.png", 
       plot = plength_vs_RRR_inn2, width = 8, height = 6, dpi = 300)


#Trying negative binomial for partnership length

#data is parametric_cricket_data

nb_length_null <- glm.nb(partnership_length ~ 1, data = parametric_cricket_data)
nb_length_full <- glm.nb(partnership_length ~ partnership_start + partnership_number + venue + innings + batting_team + bowling_team,
                         data = parametric_cricket_data)

AIC(nb_length_null, nb_length_full)
BIC(nb_length_null, nb_length_full)

#full model way better

nb_length_model <- glm.nb(partnership_length ~ ns(partnership_start, df = 3),
                          data = parametric_cricket_data,
                          link = "log")
AIC(nb_length_model, nb_length_full)
BIC(nb_length_model, nb_length_full)

#compare predictive performance
#parametric_cricket data
#validation_cricket_data

#Fit Negative Binomial on training data
#nb_length_model

#Predict on test data
pred_nb_length <- predict(nb_length_model, newdata = validation_cricket_data, type = "response")

#Calculate MAE and RMSE
mae_nb_length <- mean(abs(pred_nb_length - validation_cricket_data$partnership_length))
rmse_nb_length <- sqrt(mean((pred_nb_length - validation_cricket_data$partnership_length)^2))

#Compare with Weibull AFT
cat("Negative Binomial - MAE:", round(mae_nb_length, 3), "RMSE:", round(rmse_nb_length, 3), "\n")
cat("Weibull AFT       - MAE: 9.076          RMSE: 12.02\n")

#Also calculate C-index
cindex_nb <- concordance(
  validation_cricket_data$partnership_length ~ pred_nb_length
)

cindex_nb

#c-index 0.654
