install.packages("tidyverse")
install.packages("readxl")
install.packages("dplyr")
library(tidyverse)
library(readxl)
library(dplyr)


############### descriptive #####################

numeric_stats <- data %>% 
  summarise(
    Mean_Age_HC = mean(Age_HC, na.rm = TRUE),
    Median_Age_HC = median(Age_HC, na.rm = TRUE),
    IQR_Age_HC = IQR(Age_HC, na.rm = TRUE),
  )
print(numeric_stats)
View(numeric_stats)

numeric_stats <- data %>% 
  group_by(Group2) %>%
  summarise(
    Mean_Age = mean(Age, na.rm = TRUE),
    Median_Age = median(Age, na.rm = TRUE),
    IQR_Age = IQR(Age, na.rm = TRUE),
    Median_CCI = median(CCI, na.rm = TRUE),
    IQR_CCI = IQR(CCI, na.rm = TRUE),
    Median_Enrollment_duration = median(Enrollment_duration, na.rm = TRUE),
    IQR_Enrollment_duration = IQR(Enrollment_duration, na.rm = TRUE),
    Median_GCS = median(GCS, na.rm = TRUE),
    Min_GCS = min(GCS, na.rm = TRUE),
    Max_GCS = max(GCS, na.rm = TRUE),
    Median_Four = median(Four, na.rm = TRUE),
    Min_Four = min(Four, na.rm = TRUE),
    Max_Four = max(Four, na.rm = TRUE)
  )
print(numeric_stats)
View(numeric_stats)


############### posthoc vs realtime #####################

kappa_result <- kappa2(data.frame(rater1 = data$fNIRS_posthoc, rater2 = data$fNIRS_realtime), "unweighted")
print(kappa_result)

percent_agreement <- mean(data$fNIRS_posthoc == data$fNIRS_realtime, na.rm = TRUE) * 100
cat("Percent Agreement: ", percent_agreement, "%\n")

confusion_matrix <- table(Reference = data$fNIRS_posthoc, Test = data$fNIRS_realtime)
print(confusion_matrix)

sensitivity <- confusion_matrix[2, 2] / sum(confusion_matrix[2, ])  # True P / All P
specificity <- confusion_matrix[1, 1] / sum(confusion_matrix[1, ])  # True N / All N

cat("Sensitivity: ", round(sensitivity * 100, 2), "%\n")
cat("Specificity: ", round(specificity * 100, 2), "%\n")


### over or underestimate

fp <- confusion_matrix[1, 2]  # False Positives
fn <- confusion_matrix[2, 1]  # False Negatives

total_cases <- sum(confusion_matrix)

fp_rate <- fp / total_cases * 100
fn_rate <- fn / total_cases * 100

cat("False Positives (Overestimation):", fp, "(", round(fp_rate, 2), "% of cases)\n")
cat("False Negatives (Underestimation):", fn, "(", round(fn_rate, 2), "% of cases)\n")


############### clinical outcomes #####################


model_oneweek <- glm(oneweek ~ fNIRS_posthoc, data = data, family = binomial)

summary(model_oneweek)
odds_ratios_oneweek <- exp(coef(model_oneweek))
conf_intervals_oneweek <- exp(confint(model_oneweek))
results_oneweek <- data.frame(
  Predictor = names(odds_ratios_oneweek),
  Odds_Ratio = odds_ratios_oneweek,
  Lower_CI = conf_intervals_oneweek[, 1],
  Upper_CI = conf_intervals_oneweek[, 2]
)
print(results_oneweek)

