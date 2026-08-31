#Data Collection, Cleaning, and EDA

#Set working directory to help save plots
setwd("C:/R_Projects/Dissertation")
dir.create("Plots", showWarnings = FALSE)

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

##Import the datasets

data_file_men <- ("C:\\R_Projects\\Dissertation\\Datasets\\21-26_matches_men.csv")
cricket_men <- read.csv(data_file_men, header = TRUE)

data_file_women <- ("C:\\R_Projects\\Dissertation\\Datasets\\21-26_matches_women.csv")
cricket_women <- read.csv(data_file_women, header = TRUE)

#remove all N/A values from the men partnership survival times
cricket_clean_men <- cricket_men %>%
  drop_na(partnership_length)

#same for women data
cricket_clean_women <- cricket_women %>%
  drop_na(partnership_length)

#create a partnerships survival object for men data
partnerships_men <- Surv(cricket_clean_men$partnership_length, cricket_clean_men$wicket)

survival_men <- survfit(partnerships_men ~ 1)
str(survival_men)

#KM curv of men data
survfit2(partnerships_men ~ 1) |> 
  ggsurvfit() +
  labs(
    x = "Balls",
    y = "Overall partnership survival probability"
  ) + 
  add_confidence_interval() +
  add_risktable()

#same for women
partnerships_women <- Surv(cricket_clean_women$partnership_length, cricket_clean_women$wicket)

survival_women <- survfit(partnerships_women ~ 1)
str(survival_women)

#KM for womens data
survfit2(partnerships_women ~ 1) |> 
  ggsurvfit() +
  labs(
    x = "Balls",
    y = "Overall partnership survival probability"
  ) + 
  add_confidence_interval() +
  add_risktable()

##Combine the data
full_cricket_data <- rbind(cricket_clean_men, cricket_clean_women)

length(unique(full_cricket_data))
#47 different variables

#changing variables to factors
full_cricket_data <- full_cricket_data %>%
  mutate(
    season = factor(season),
    venue = factor(venue),
    innings = factor(innings),
    batting_team = factor(batting_team),
    bowling_team = factor(bowling_team),
    wicket_sum = factor(wicket_sum)
  )

#provide an overview of the data
summary(full_cricket_data)
sapply(full_cricket_data, class)

#matches by year
full_cricket_data %>%
  group_by(season) %>%
  summarise(
    matches = n_distinct(match_id),
    .groups = "drop"
  )

#as the data was obtained halfway through the 2026 season, we have 35 matches from this.
#Will not be using these, so remove.

full_cricket_data <- full_cricket_data %>%
  filter(season != 2026)

#matches by gender
summary(full_cricket_data$innings)
full_cricket_data %>%
  group_by(gender) %>%
  summarise(
    matches = n_distinct(match_id),
    .groups = "drop"
  )

#total partnerships
nrow(full_cricket_data)
#4274

#remove innings to only be 1/2 - due to super over
full_cricket_data <- full_cricket_data %>%
  filter(innings %in% c(1,2)) %>%
  mutate(innings = factor(innings))

nrow(full_cricket_data)

#this removes 2 observations for innings 3, and 1 for innings 4. so 3 in total.

nrow(full_cricket_data$gender)
full_cricket_data %>%
  count(gender)
full_cricket_data %>%
  count(season)

#we need to see if there are any other data points to remove

#check histogram of innings_balls
ggplot(full_cricket_data, aes(x = innings_balls)) +
  geom_histogram(binwidth = 5, boundary = 0) +
  labs(
    title = "Distribution of Balls per Innings",
    x = "Balls per innings",
    y = "Number of innings"
  ) +
  theme_minimal()

ggplot(full_cricket_data, aes(x = innings_balls)) +
  geom_histogram(binwidth = 5, boundary = 0) +
  facet_wrap(~ gender) +
  labs(
    title = "Distribution of Balls per Innings by Gender",
    x = "Balls per innings",
    y = "Number of innings"
  ) +
  theme_minimal()
#we can see some low values. this needs investigating.

#checking how many innings below 50 balls
low_innings <- full_cricket_data %>%
  filter(innings_balls < 50) %>%
  dplyr::select(match_id, gender, season, innings, innings_balls)

low_innings

full_cricket_data %>%
  filter(innings_balls < 40) %>%
  dplyr::select(match_id, gender, season, innings, innings_balls, batting_team, bowling_team)


full_cricket_data %>%
  filter(innings_balls < 40) %>%
  nrow()

#Lets remove any innings where it is not complete in any form - i.e. 100 balls 10 wickets or target reached
#this is an idea to check, or i could do it manually as above.

clean_cricket_data <- full_cricket_data %>%
  filter(
    !is.na(innings_balls),
    (
      (innings == 1 & innings_balls < 100 & innings_wickets < 10) |
        (innings == 2 & innings_balls < 100 & innings_wickets < 10 & innings_runs < target)
    )
  ) %>%
  dplyr::select(
    match_id, gender, season, innings,
    innings_balls, innings_wickets, innings_runs,
    target, batting_team, bowling_team
  ) %>%
  arrange(match_id)

clean_cricket_data

#create an object to represent these innings.
reduced_innings <- clean_cricket_data

#there are 26 innings that are incomplete.
#from here, I want to remove any MATCH from an incomplete innings 1, and innings from any incomplete innings 2.


#I want to remove the matches that these 26 innings are involved from the data. can use anti_join()
#create a backup of the data
full_cricket_data_original <- full_cricket_data

#make a function for incomplete first innings matches
incomplete_first_innings <- reduced_innings %>%
  filter(as.character(innings) == "1") %>%
  distinct(match_id)

#add the second innings of these matches
innings_to_remove <- bind_rows(
  reduced_innings %>%
    mutate(innings = as.character(innings)) %>%
    dplyr::select(match_id, innings),
  
  incomplete_first_innings %>%
    mutate(innings = "2") %>%
    dplyr::select(match_id, innings)
) %>%
  distinct(match_id, innings)

clean_cricket_data <- full_cricket_data %>%
  mutate(innings = as.character(innings)) %>%
  anti_join(
    innings_to_remove,
    by = c("match_id", "innings")
  )

nrow(innings_to_remove)

#this identifies 36 innings to remove.
#lets remove them and get our FINAL dataset.

#check if these have been removed
clean_cricket_data %>%
  semi_join(
    innings_to_remove,
    by = c("match_id", "innings")
  ) %>%
  distinct(match_id, innings) %>%
  nrow()

full_cricket_data %>%
  distinct(match_id, innings) %>%
  nrow()

clean_cricket_data %>%
  distinct(match_id, innings) %>%
  nrow()

innings_to_remove %>%
  anti_join(
    full_cricket_data %>%
      distinct(match_id, innings),
    by = c("match_id", "innings")
  )
#this confirms that.
#continue with clean_cricket_data as your dataset.

#total partnerships now
nrow(clean_cricket_data)

#percentage reduction
percentage_reduction <- (nrow(full_cricket_data) - nrow(clean_cricket_data)) /
  nrow(full_cricket_data) * 100

percentage_reduction

###Data description
summary(clean_cricket_data)

clean_cricket_data %>%
  count(gender) %>%
  mutate(
    percentage = n / sum(n) * 100
  )
clean_cricket_data %>%
  count(season) %>%
  mutate(
    percentage = n / sum(n) * 100
  )
clean_cricket_data %>%
  count(innings) %>%
  mutate(
    percentage = n / sum(n) * 100
  )

#partnership length details
partnership_length_stats <- clean_cricket_data %>%
  summarise(
    n = sum(!is.na(partnership_length)),
    mean = mean(partnership_length, na.rm = TRUE),
    SD = sd(partnership_length, na.rm = TRUE),
    median = median(partnership_length, na.rm = TRUE),
    Q1 = quantile(partnership_length, 0.25, na.rm = TRUE),
    Q3 = quantile(partnership_length, 0.75, na.rm = TRUE),
    IQR = IQR(partnership_length, na.rm = TRUE),
    minimum = min(partnership_length, na.rm = TRUE),
    maximum = max(partnership_length, na.rm = TRUE),
    censored = sum(wicket == 0, na.rm = TRUE),
    percent_censored = mean(wicket == 0, na.rm = TRUE) * 100
  )

partnership_length_stats

#partnership score details
partnership_score_stats <- clean_cricket_data %>%
  summarise(
    n = sum(!is.na(partnership_score)),
    mean = mean(partnership_score, na.rm = TRUE),
    SD = sd(partnership_score, na.rm = TRUE),
    median = median(partnership_score, na.rm = TRUE),
    Q1 = quantile(partnership_score, 0.25, na.rm = TRUE),
    Q3 = quantile(partnership_score, 0.75, na.rm = TRUE),
    IQR = IQR(partnership_score, na.rm = TRUE),
    minimum = min(partnership_score, na.rm = TRUE),
    maximum = max(partnership_score, na.rm = TRUE),
    zeros = sum(partnership_score == 0, na.rm = TRUE),
    percent_zeros = mean(partnership_score == 0, na.rm = TRUE) * 100
  )

partnership_score_stats


#predictor variables descriptions

#partnership start
#histogram
partnership_start_hist <- ggplot(clean_cricket_data, aes(x = partnership_start)) +
  geom_histogram(binwidth = 5, boundary = 0, color = "black", fill = "#0099F8", alpha = 0.8) +
  scale_x_continuous(breaks = seq(0, 120, by=10)) +
  labs(
    title = "Distribution of Partnership Start",
    x = "Ball Number at Partnership Start",
    y = "Number of Partnerships"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  )

ggsave("Plots/EDA_1_Partnership_Start_Distribution.png", 
       plot = partnership_start_hist, width = 8, height = 6, dpi = 300)

#heavily skewed due to lots of values of 1. first partnership per innings is always 1 thats why!
#the rest is surprisingly uniform, maybe a slight increase towards the end which is expected.

#partnership number stats
clean_cricket_data %>%
  count(partnership_number, name = "n_partnerships")

#Bar Chart
partnership_start_bar <- ggplot(clean_cricket_data, aes(x = factor(partnership_number))) +
  geom_bar(fill = "#2ca02c", color = "black", alpha = 0.8) +
  scale_y_continuous(breaks = seq(0, 700, by=100)) +
  labs(
    title = "Bar Chart: Distribution of Partnership Number",
    x = "Partnership Number",
    y = "Number of Partnerships"
  ) +
  theme_gray() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  )

ggsave("Plots/EDA_2_Partnership_number_barchart.png", 
       plot = partnership_start_bar, width = 8, height = 6, dpi = 300)

#compare between innings
partnership_start_byinnings_bar <- ggplot(clean_cricket_data, aes(x = factor(partnership_number),
                                                                  fill = factor(innings))) +
  geom_bar(position = "dodge", alpha = 0.8) +
  scale_fill_manual(name = "Innings",
                    values = c("1" = "#1f77b4", "2" = "#ffd700"),
                    labels = c("1" = "Innings 1 (Batting)", "2" = "Innings 2 (Chasing)")) +
  scale_y_continuous(breaks = seq(0, 400, by=50)) +
  labs(
    title = "Partnership Number Distribution: Innings 1 vs Innings 2",
    x = "Partnership Number",
    y = "Number of Partnerships"
  ) +
  theme_gray() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11),
    legend.position = "top",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10)
  )

ggsave("Plots/EDA_3_Partnership_number_byinnings_barchart.png", 
       plot = partnership_start_byinnings_bar, width = 8, height = 6, dpi = 300)

#innings 2 seems to be lower to start but higher towards end.

#venue
ggplot(clean_cricket_data, aes(x = venue, y = innings_runs)) +
  geom_boxplot() +
  labs(
    title = "Distribution of Innings Scores by Venue",
    x = "Venue",
    y = "Innings Score (runs)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

#by gender side by side
inningscores_venue_gender <- ggplot(
  clean_cricket_data,
  aes(x = venue, y = innings_runs, fill = factor(gender))
) +
  geom_boxplot(position = position_dodge(), alpha = 0.7, color = "black") +
  scale_fill_manual(name = "Gender",
                    values = c("men" = "blue", "women" = "red"),
                    labels = c("men" = "Men", "women" = "Women")) +
  labs(
    title = "Distribution of Innings Scores by Venue and Gender",
    x = "Venue",
    y = "Innings Score (runs)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 11),
    legend.position = "top",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10)
  )

ggsave("Plots/EDA_4_Inningscores_venue_gender.png", 
       plot = inningscores_venue_gender, width = 8, height = 6, dpi = 300)

##compare required run rate with partnership length
#Filter to innings 2 only (where RRR exists)
innings_2_clean_data <- clean_cricket_data %>%
  filter(innings == 2)

#run rate + required run rate
#Run rate at start (innings 1 and 2)
clean_cricket_data %>%
  summarise(
    n = sum(!is.na(partnership_rr_start)),
    mean_rr_start = mean(partnership_rr_start, na.rm = TRUE),
    sd_rr_start = sd(partnership_rr_start, na.rm = TRUE),
    median_rr_start = median(partnership_rr_start, na.rm = TRUE),
    min_rr_start = min(partnership_rr_start, na.rm = TRUE),
    max_rr_start = max(partnership_rr_start, na.rm = TRUE)
  )

#Required run rate (innings 2 only)
innings_2_clean_data %>%
  summarise(
    n = sum(!is.na(partnership_required_rr_start)),
    mean_rrr = mean(partnership_required_rr_start, na.rm = TRUE),
    sd_rrr = sd(partnership_required_rr_start, na.rm = TRUE),
    median_rrr = median(partnership_required_rr_start, na.rm = TRUE),
    min_rrr = min(partnership_required_rr_start, na.rm = TRUE),
    max_rrr = max(partnership_required_rr_start, na.rm = TRUE)
  )



#Correlation between RRR and partnership length
cor_test <- cor.test(innings_2_clean_data$partnership_required_rr_start, 
                     innings_2_clean_data$partnership_length, 
                     method = "pearson")
cor_test

#-0.19 correlation test suggests a negative relationship as expected
#visualisation
partnership_length_vs_RRR <- ggplot(innings_2_clean_data, 
                                    aes(x = partnership_required_rr_start, 
                                        y = partnership_length)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "gam",
              formula = y ~ s(x, bs = "cs"), linewidth = 1.5,
              se = TRUE) +
  coord_cartesian(xlim = c(0, 30)) +
  scale_x_continuous(breaks = seq(0, 40, by=5)) +
  scale_y_continuous(breaks = seq(0, 100, by=25)) +
  labs(
    title = "Partnership Length vs Required Run Rate",
    x = "Required Run Rate at partnership start (runs per over)",
    y = "Partnership length (balls)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  )

ggsave("Plots/EDA_7_partnership_length_vs_RRR.png", 
       plot = partnership_length_vs_RRR, width = 8, height = 6, dpi = 300)

#coord_cartesian zooms in on the important part.


##key variables
#wickets
wickets_censored_stackbar <- ggplot(clean_cricket_data, aes(x = factor(innings), fill = factor(wicket))) +
  geom_bar() +
  labs(
    title = "Stacked Bar Plot showing Proportion of Not Out Partnerships per Innings",
    x = "Innings",
    y = "Number of partnerships",
    fill = "Partnership outcome"
  ) +
  scale_fill_manual(
    values = c("0" = "forestgreen",
               "1" = "indianred3"),
    labels = c("0" = "Censored",
               "1" = "Non-censored")
  ) +
  scale_y_continuous(breaks = seq(0, 2500, 500)) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  )

ggsave("Plots/EDA_8_wickets_censored.png", 
       plot = wickets_censored_stackbar, width = 8, height = 6, dpi = 300)

#distribution of partnership_length
partnership_length_distribution <- ggplot(clean_cricket_data, aes(x = partnership_length)) +
  geom_histogram(binwidth = 1, boundary = 0.5, fill = "darkblue", color = "white", alpha = 0.8) +
  scale_x_continuous(breaks = seq(0, 90, by=5)) +
  scale_y_continuous(breaks = seq(0, 600, by=50)) +
  labs(
    title = "Distribution of Partnership Length",
    x = "Partnership Length (balls)",
    y = "Number of Partnerships"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  )

ggsave("Plots/EDA_5_Partnership_length_distribution.png", 
       plot = partnership_length_distribution, width = 8, height = 6, dpi = 300)

#this shows a clear decreasing slope

#partnership score histogram
n_total   <- nrow(clean_cricket_data)
n_zero    <- sum(clean_cricket_data$partnership_score == 0, na.rm = TRUE)
pct_zero  <- round(100 * n_zero / n_total, 1)

partnership_score_distribution <- ggplot(clean_cricket_data, aes(x = partnership_score)) +
  geom_histogram(binwidth = 1, boundary = 0.5, aes(fill = partnership_score == 0), colour = "white") +
  scale_x_continuous(breaks = seq(0, 150, by=10)) +
  scale_y_continuous(breaks = seq(0, 400, by=50)) +
  scale_fill_manual(
    values = c("FALSE" = "steelblue", "TRUE" = "red"),
    guide = "none"
  ) +
  annotate("text", x = 0, y = n_zero, label = paste0(pct_zero, "%"), vjust = -0.5, colour = "red", fontface = "bold") +
  labs(
    title = "Distribution of Partnership Score",
    x = "Partnership Score (runs)",
    y = "Frequency"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  )

ggsave("Plots/EDA_5_Partnership_score_distribution.png", 
       plot = partnership_score_distribution, width = 8, height = 6, dpi = 300)

clean_cricket_data %>%
  summarise(
    skewness(partnership_score, na.rm = TRUE),
    kurtosis(partnership_score, na.rm = TRUE))


clean_cricket_data %>%
  arrange(desc(partnership_length)) %>%
  dplyr::select(match_id, gender, season, innings, partnership_length) %>%
  head(20)


#partnerships per innings
partnerships_per_innings <- clean_cricket_data %>%
  filter(!is.na(partnership_length)) %>%
  group_by(match_id, gender, season, innings) %>%
  summarise(
    partnerships = n(),
    .groups = "drop"
  )


#EDA
###Relationships between variables
##comparing length vs score
partnership_length_vs_score <- ggplot(clean_cricket_data, aes(x = partnership_length, y = partnership_score)) +
  geom_point(alpha = 0.3) +
  geom_smooth(
    method = "gam",
    formula = y ~ s(x, bs = "cs"), linewidth = 1.5,
    se = TRUE
  ) +
  scale_x_continuous(breaks = seq(0, 100, by=25)) +
  scale_y_continuous(breaks = seq(0, 150, by=25)) +
  labs(
    x = "Partnership length (balls)",
    y = "Partnership score (runs)",
    title = "Partnership score vs partnership length"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  )

ggsave("Plots/EDA_6_Partnership_length_vs_score.png", 
       plot = partnership_length_vs_score, width = 8, height = 6, dpi = 300)

cor.test(
  clean_cricket_data$partnership_length,
  clean_cricket_data$partnership_score,
  method = "pearson"
)
cor.test(
  clean_cricket_data$partnership_length,
  clean_cricket_data$partnership_score,
  method = "spearman"
)

#Partnership score vs start (better)
partnership_score_vs_start <- ggplot(clean_cricket_data,
                                     aes(x = partnership_start,
                                         y = partnership_score,
                                         color = innings)) +
  geom_point(alpha = 0.4, size = 1.5) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), 
              size = 1.2, se = TRUE) +
  scale_color_manual(values = c("1" = "#1f77b4", "2" = "#ffd700"),
                     labels = c("1st Innings (Batting first)", "2nd Innings (Chasing)")) +
  scale_x_continuous(breaks = seq(0, 100, 10)) +
  labs(
    title = "Partnership Score vs Ball number at Partnership Start by Innings",
    x = "Ball Number in Innings",
    y = "Runs Scored in Partnership",
    color = "Innings"
  ) +
  theme_minimal() +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  )

ggsave("Plots/EDA_6_Partnership_score_vs_start.png", 
       plot = partnership_score_vs_start, width = 8, height = 6, dpi = 300)

