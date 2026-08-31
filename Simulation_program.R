##Match Simulator

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

#define survival and regression model
survival_model <- flexsurvreg(
  surv_obj ~ ns(partnership_start, df = 3),
  data = parametric_cricket_data,
  dist = "weibull"
)

regression_model <- zeroinfl(partnership_score ~ ns(partnership_start, df=4) + partnership_number | 1,
                             data = score_data, dist = "negbin")

#how many innings are there
total_innings <- sa_cricket_data %>%
  filter(!is.na(innings_runs))

no_innings <- nrow(total_innings)

#How many of those ended in censoring
censored_innings <- total_innings %>%
  filter(wicket == 0) %>%
  nrow()

#Proportion of final partnerships that are censored
censoring_prob <- censored_innings/no_innings
final_partnership_censoring_rate <- (censored_innings / no_innings) * 100

cat("Total innings:", no_innings, "\n")
cat("Censored (innings end):", censored_innings, "\n")
cat("Censoring rate for final partnerships:", final_partnership_censoring_rate, "%", "\n")

#Simulate one innings 1 match
simulate_innings_1 <- function(survival_model, regression_model, censoring_prob = censoring_prob) {
  
  partnership <- 1
  total_balls <- 0
  total_runs <- 0
  total_wickets <- 0
  partnerships_list <- list()
  innings_ended_not_on_wicket <- FALSE
  
  while (partnership <= 10 & total_balls < 100) {
    
    #Create data for this partnership
    newdata <- data.frame(
      partnership_start = total_balls,
      partnership_number = partnership
    )
    
    #Predict partnership length (from survival model)
    pred_length_raw <- predict(survival_model, newdata = newdata, type = "response")
    #Extract numeric value if it's in a list or data frame
    if (is.list(pred_length_raw)) {
      pred_length <- as.numeric(pred_length_raw[[1]])
    } else {
      pred_length <- as.numeric(pred_length_raw[1])
    }
    #Add some random noise around the prediction
    partnership_length <- max(1, round(pred_length + rnorm(1, 0, 3)))
    
    #Cap at 100 balls if necessary
    if (total_balls + partnership_length > 100) {
      partnership_length <- 100 - total_balls
      #Randomly decide: does innings end by censoring or wicket
      innings_ended_not_on_wicket <- rbinom(1, 1, censoring_prob) == 1
      is_censored <- innings_ended_not_on_wicket
    } else {
      is_censored <- FALSE
    }
    
    #Predict partnership score (from regression model)
    pred_score_raw <- predict(regression_model, newdata = newdata, type = "response")
    if (is.list(pred_score_raw)) {
      partnership_score <- max(0, round(as.numeric(pred_score_raw[[1]])))
    } else {
      partnership_score <- max(0, round(as.numeric(pred_score_raw[1])))
    }
    
    #Update totals
    total_balls <- total_balls + partnership_length
    total_runs <- total_runs + partnership_score
    
    #Count wickets only if not censored
    if (!is_censored) {
      total_wickets <- total_wickets + 1
    }
    
    #Store partnership info
    partnerships_list[[partnership]] <- data.frame(
      partnership = partnership,
      start_ball = total_balls - partnership_length + 1,
      balls_faced = partnership_length,
      runs_scored = partnership_score,
      cumulative_balls = total_balls,
      cumulative_runs = total_runs,
      ended_by = ifelse(is_censored, "Innings End - No Wicket", "Wicket")
    )
    
    partnership <- partnership + 1
    
    #Stop if we've completed all 10 partnerships or reached 100 balls
    if (total_balls >= 100 | partnership > 10) break
  }
  
  #Combine results
  partnerships_df <- do.call(rbind, partnerships_list)
  rownames(partnerships_df) <- NULL
  
  return(list(
    final_score = total_runs,
    final_balls = total_balls,
    final_wickets = total_wickets,
    partnerships = partnerships_df
  ))
}

#Run the simulation
result <- simulate_innings_1(survival_model, regression_model, censoring_prob = censoring_prob )

#Print results
cat("SIMULATED INNINGS 1\n")
cat("Final Score:", result$final_score, "\n")
cat("Final Balls:", result$final_balls, "\n\n")
cat("Final Wickets:", result$final_wickets, "\n\n")
cat("Partnership Breakdown:\n")
print(result$partnerships)

#time to loop the simulation
#run simulation 1000 times
num_simulations <- 1000
results_list <- list()

for (i in 1:num_simulations) {
  result <- simulate_innings_1(survival_model, regression_model, censoring_prob = censoring_prob)
  results_list[[i]] <- result
  
  cat("SIMULATION", i, "\n")
  cat("Final Score:", result$final_score, "\n")
  cat("Final Balls:", result$final_balls, "\n")
  cat("Final Wickets:", result$final_wickets, "\n\n")
  cat("Partnership Breakdown:\n")
  print(result$partnerships)
}

#Summary statistics across all simulations
scores <- sapply(results_list, function(x) x$final_score)
cat("SUMMARY OF ALL 1000 SIMULATIONS\n")
cat("Mean score:", round(mean(scores), 2), "\n")
cat("SD score:", round(sd(scores), 2), "\n")
cat("Min score:", min(scores), "\n")
cat("Max score:", max(scores), "\n")

#Need to simulate innings 2
#then need to simulate a match to get an outcome

innings2_survmodel <- flexsurvreg(Surv(partnership_length, wicket) ~ ns(partnership_start, df=3) + ns(required_run_rate, df = 4),
                                  data = parametric_innings2,
                                  dist = "weibull")

innings2_zinbmodel <- zeroinfl(partnership_score ~ ns(partnership_start, df=4) + partnership_number + ns(required_run_rate, df = 3) | 1,
                               data = score_innings2, dist = "negbin")


simulate_innings_2 <- function(innings2_survmodel, innings2_zinbmodel, target, censoring_prob = censoring_prob) {
  
  partnership <- 1
  total_balls <- 0
  total_runs <- 0
  total_wickets <- 0
  partnerships_list <- list()
  
  while (partnership <= 10 & total_balls < 100) {
    
    #Calculate required run rate based on current situation
    runs_required <- target - total_runs
    balls_remaining <- 100 - total_balls
    required_run_rate <- (runs_required / balls_remaining) * 6  # Convert to per-over
    
    #Cap at 0 if target already reached
    required_run_rate <- max(0, required_run_rate)
    
    #Create data for this partnership
    newdata <- data.frame(
      partnership_start = total_balls,
      partnership_number = partnership,
      required_run_rate = required_run_rate
    )
    
    #Predict partnership length (from innings 2 survival model)
    pred_length_raw <- predict(innings2_survmodel, newdata = newdata, type = "response")
    if (is.list(pred_length_raw)) {
      pred_length <- as.numeric(pred_length_raw[[1]])
    } else {
      pred_length <- as.numeric(pred_length_raw[1])
    }
    
    partnership_length <- max(1, round(pred_length + rnorm(1, 0, 3)))
    
    #Cap at 100 balls if necessary
    if (total_balls + partnership_length > 100) {
      partnership_length <- 100 - total_balls
      #Randomly decide: does innings end by censoring or wicket
      is_censored <- rbinom(1, 1, censoring_prob) == 1
    } else {
      is_censored <- FALSE
    }
    
    #Predict partnership score (from innings 2 regression model)
    pred_score_raw <- predict(innings2_zinbmodel, newdata = newdata, type = "response")
    if (is.list(pred_score_raw)) {
      partnership_score <- max(0, round(as.numeric(pred_score_raw[[1]])))
    } else {
      partnership_score <- max(0, round(as.numeric(pred_score_raw[1])))
    }
    
    #Update totals
    total_balls <- total_balls + partnership_length
    total_runs <- total_runs + partnership_score
    
    #Count wickets only if not censored
    if (!is_censored) {
      total_wickets <- total_wickets + 1
    }
    
    #Check if target reached
    target_reached <- total_runs >= target
    
    #Store partnership info
    partnerships_list[[partnership]] <- data.frame(
      partnership = partnership,
      start_ball = total_balls - partnership_length + 1,
      balls_faced = partnership_length,
      runs_scored = partnership_score,
      cumulative_balls = total_balls,
      cumulative_runs = total_runs,
      wickets_down = total_wickets,
      ended_by = ifelse(is_censored, "Innings End", "Wicket"),
      target_reached = target_reached
    )
    
    partnership <- partnership + 1
    
    #Stop if we've completed all 10 partnerships, reached 100 balls, or target reached
    if (total_balls >= 100 | partnership > 10 | target_reached | total_wickets >= 10) break
  }
  
  partnerships_df <- do.call(rbind, partnerships_list)
  rownames(partnerships_df) <- NULL
  
  return(list(
    final_score = total_runs,
    final_balls = total_balls,
    final_wickets = total_wickets,
    target = target,
    target_reached = total_runs >= target,
    partnerships = partnerships_df
  ))
}

#Run full match simulation inside a function
simulate_full_match <- function(survival_model, regression_model, 
                                innings2_survmodel, innings2_zinbmodel, 
                                censoring_prob = censoring_prob) {
  
  inn1_result <- simulate_innings_1(survival_model, regression_model, censoring_prob = censoring_prob)
  cat("INNINGS 1\n")
  cat("Final Score:", inn1_result$final_score, "\n")
  cat("Balls:", inn1_result$final_balls, "\n")
  cat("Wickets:", inn1_result$final_wickets, "\n\n")
  
  #Use innings 1 score + 1 as target
  target <- inn1_result$final_score + 1
  
  #Simulate innings 2
  inn2_result <- simulate_innings_2(innings2_survmodel, innings2_zinbmodel, target, censoring_prob = censoring_prob)
  cat("INNINGS 2 (Target:", target, ")\n")
  cat("Final Score:", inn2_result$final_score, "\n")
  cat("Balls:", inn2_result$final_balls, "\n")
  cat("Wickets:", inn2_result$final_wickets, "\n")
  cat("Target Reached:", inn2_result$target_reached, "\n\n")
  
  #Match result
  if (inn2_result$final_score == inn1_result$final_score) {
    cat("MATCH RESULT: TIE\n")
    cat("Team 1:", paste(inn1_result$final_score, sep = "/", inn1_result$final_wickets), "from", inn1_result$final_balls, "balls\n")
    cat("Team 2:", paste(inn2_result$final_score, sep = "/", inn2_result$final_wickets), "from", inn2_result$final_balls, "balls (target:", target, "). Match Tied.\n")
  } else if (inn2_result$target_reached) {
    cat("MATCH RESULT: TEAM 2 WON\n")
    cat("Team 1:", paste(inn1_result$final_score, sep = "/", inn1_result$final_wickets), "from", inn1_result$final_balls, "balls\n")
    cat("Team 2:", paste(inn2_result$final_score, sep = "/", inn2_result$final_wickets), "from", inn2_result$final_balls, "balls (target:", target, "). Won by", 10 - inn2_result$final_wickets, "wickets (with", 100 - inn2_result$final_balls, "balls remaining).\n")
    cat("TEAM 2 WINS BY", 10 - inn2_result$final_wickets, "WICKETS.\n")
  } else {
    cat("MATCH RESULT: TEAM 1 WON\n")
    cat("Team 1:", paste(inn1_result$final_score, sep = "/", inn1_result$final_wickets), "from", inn1_result$final_balls, "balls\n")
    cat("Team 2:", paste(inn2_result$final_score, sep = "/", inn2_result$final_wickets), "from", inn2_result$final_balls, "balls (target:", target, "). Lost by", target - inn2_result$final_score - 1, "runs.\n")
    cat("TEAM 1 WINS BY", target - inn2_result$final_score - 1, "RUNS.\n")
  }
  
  invisible(list(innings1 = inn1_result, innings2 = inn2_result))
}

#Use it:
simulate_full_match(survival_model, regression_model, 
                    innings2_survmodel, innings2_zinbmodel,
                    censoring_prob = censoring_prob)


#time to loop the whole match simulation
#Run simulation 10 times
match_num_simulations <- 1000
match_results <- vector("list", match_num_simulations)

for (i in 1:match_num_simulations) {
  
  match_results[[i]] <- simulate_full_match(
    survival_model,
    regression_model,
    innings2_survmodel,
    innings2_zinbmodel,
    censoring_prob = censoring_prob
  )
}

#Summary statistics across all simulations
inn1_scores <- sapply(
  match_results,
  function(x) x$innings1$final_score
)
inn2_scores <- sapply(
  match_results,
  function(x) x$innings2$final_score
)
inn1_wickets <- sapply(
  match_results,
  function(x) x$innings1$final_wickets
)
inn2_wickets <- sapply(
  match_results,
  function(x) x$innings2$final_wickets
)
inn1_balls <- sapply(
  match_results,
  function(x) x$innings1$final_balls
)
inn2_balls <- sapply(
  match_results,
  function(x) x$innings2$final_balls
)
target_reached <- sapply(
  match_results,
  function(x) x$innings2$target_reached
)

cat("1000 SIMULATED MATCHES\n\n")

cat("Mean Innings 1 score:",
    round(mean(inn1_scores), 2), "\n")
cat("Mean Innings 2 score:",
    round(mean(inn2_scores), 2), "\n")
cat("Mean Innings 1 wickets:",
    round(mean(inn1_wickets), 2), "\n")
cat("Mean Innings 2 wickets:",
    round(mean(inn2_wickets), 2), "\n")
cat("Mean Innings 1 balls:",
    round(mean(inn1_balls), 2), "\n")
cat("Mean Innings 2 balls:",
    round(mean(inn2_balls), 2), "\n")
cat("Percentage of chases successful:",
    round(mean(target_reached) * 100, 2), "%\n")

cat("Innings 1:\n")
cat("Minimum:", min(inn1_scores), "\n")
cat("Median:", median(inn1_scores), "\n")
cat("Maximum:", max(inn1_scores), "\n\n")

cat("Innings 2:\n")
cat("Minimum:", min(inn2_scores), "\n")
cat("Median:", median(inn2_scores), "\n")
cat("Maximum:", max(inn2_scores), "\n")
