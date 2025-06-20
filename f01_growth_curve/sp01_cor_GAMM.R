# Install and load the necessary packages
library(itsadug)
library(mgcv)

library(dplyr)
library(zoo)


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
dd <- read.csv("../raw_data/full_res_for_gamm00_use.csv",header = F)

for (i in 2:48) {
  
df <- data.frame(
  brain_volume = dd[[paste0("V", i)]], 
  PMA = dd$V86,
  Group = factor(dd$V1, levels = c(0, 1))
)

dt_pt <- df %>% filter(Group == 0)  # Preterm data
df_ft <- df %>% filter(Group == 1)  # Full-term data


model <- gamm(
  brain_volume ~ Group + s(PMA, by = Group,k=3),
  data = df,
  method = "REML"
)


summary(model)
#plot(model, pages = 1)


library(ggplot2)

# Predict values for plotting
#df$pred <- predict.gam(model$gam, type = "link", se.fit = TRUE)



# compute for MATLAB produce.
# Get confidence intervals for predictions (95% CI)
#conf_int <- predict(model$gam, type = "link", se.fit = TRUE)
conf_int <- predict.gam(model$gam, type = "link", se.fit = TRUE)

# Calculate 95% CI bounds (2.5% and 97.5% quantiles)
upper <- conf_int$fit + 1.96 * conf_int$se.fit
lower <- conf_int$fit - 1.96 * conf_int$se.fit
# Calculate 5% and 95% CI bounds (90% CI)
lower_5 <- conf_int$fit + 1.645 * conf_int$se.fit  # 5% CI
upper_95 <- conf_int$fit - 1.645 * conf_int$se.fit  # 95% CI

# save the data for GAMM,.
GAMM_data <- data.frame(
  PMA = df$PMA,
  Group = df$Group,
  brain_volume = conf_int$fit,  # Predicted values
  lower_5 = lower_5,
  upper_95 = upper_95,
  lower = lower,
  upper = upper,
  raw_y = df$brain_volume,
  raw_x = df$PMA
)


# part for compute percentile values,
window_size <- 2

# Function to calculate percentiles within a moving PMA range
moving_percentile <- function(pma, y, target_pma, p) {
  in_window <- which(abs(pma - target_pma) <= window_size)
  if (length(in_window) < 2) return(NA) # Need at least 2 points
  return(quantile(y[in_window], probs = p, na.rm = TRUE))
}

# Ensure data is arranged by PMA within each group
df <- df %>%
  arrange(Group, PMA) 

# Separate calculation for each group
percentile_df <- df %>%
  group_by(Group) %>%
  mutate(
    p_5  = sapply(PMA, function(p) moving_percentile(PMA, brain_volume, p, 0.05)),
    p_25 = sapply(PMA, function(p) moving_percentile(PMA, brain_volume, p, 0.25)),
    p_75 = sapply(PMA, function(p) moving_percentile(PMA, brain_volume, p, 0.75)),
    p_95 = sapply(PMA, function(p) moving_percentile(PMA, brain_volume, p, 0.95))
  ) %>%
  ungroup()

###
# another way to produce percentiles, from 2014 holland paper,
# Predict the mean (mu) at each PMA value
df$mu <- predict.gam(model$gam, newdata = df, type = "response")

# Estimate the standard deviation (sigma) from the residuals (using residuals of the GAMM model)
# Step 2: Get residuals and predicted mu
df$mu <- predict(model$gam, newdata = df)
df$residuals <- df$brain_volume - df$mu
df$abs_residual <- abs(df$residuals)  # optional: use absolute residuals

# Step 3: Fit a GAM to model local sigma as a function of PMA

#  weighting compensate for the small data points at the end. Or the sigma woul go crazy.
df$weights <- ifelse(
  df$PMA <= 46,
  1,  # keep weights constant up to 44w
  1 / (1 + ((df$PMA - 46)^1.1))  # decay after 44w using a power function
)#df$weights <- exp(-((df$PMA - 40)^2) / (2 * 1.5^2)) 
sigma_model <- gam(
  abs_residual ~ s(PMA, by = Group, k = 1, sp =20),
  weights = df$weights,
  data = df
)

#sigma_model <- gam(abs_residual ~ s(PMA, by = Group,k=1, sp = 20), data = df)

# Step 4: Predict local sigma
df$sigma <- predict(sigma_model, newdata = df)



# Define the percentiles (5%, 25%, 75%, 95%)
percentiles <- c(0.05, 0.25, 0.75, 0.95)

# Calculate the percentiles using the inverse normal CDF
# Calculate the percentiles using the inverse normal CDF
df$percentiles_5 <- df$mu + df$sigma * qnorm(0.05)
df$percentiles_25 <- df$mu + df$sigma * qnorm(0.25)
df$percentiles_50 <- df$mu + df$sigma * qnorm(0.5)
df$percentiles_75 <- df$mu + df$sigma * qnorm(0.75)
df$percentiles_95 <- df$mu + df$sigma * qnorm(0.95)

smoothed_p5 <- smooth.spline(df$PMA, df$percentiles_5)
smoothed_p25 <- smooth.spline(df$PMA, df$percentiles_25)
smoothed_p50 <- smooth.spline(df$PMA, df$percentiles_50)
smoothed_p75 <- smooth.spline(df$PMA, df$percentiles_75)
smoothed_p95 <- smooth.spline(df$PMA, df$percentiles_95)

df$smoothed_p5 <- predict(smoothed_p5, df$PMA)$y
df$smoothed_p25 <- predict(smoothed_p25, df$PMA)$y
df$smoothed_p50 <- predict(smoothed_p50, df$PMA)$y
df$smoothed_p75 <- predict(smoothed_p75, df$PMA)$y
df$smoothed_p95 <- predict(smoothed_p95, df$PMA)$y
# Plot the smoothed percentiles



# smooth curve fit saparately for pt and term,'
# Create separate data frames for preterm and full-term groups
df_pt <- df[df$Group == 0, ]  # Preterm group (Group == 0)
df_ft <- df[df$Group == 1, ]  # Full-term group (Group == 1)

# Compute smoothed percentiles for preterm group
smoothed_p5_pt <- smooth.spline(df_pt$PMA, df_pt$percentiles_5)
smoothed_p25_pt <- smooth.spline(df_pt$PMA, df_pt$percentiles_25)
smoothed_p50_pt <- smooth.spline(df_pt$PMA, df_pt$percentiles_50)
smoothed_p75_pt <- smooth.spline(df_pt$PMA, df_pt$percentiles_75)
smoothed_p95_pt <- smooth.spline(df_pt$PMA, df_pt$percentiles_95)

# Extract the smoothed percentiles for preterm group
df_pt$smoothed_p5 <- predict(smoothed_p5_pt, df_pt$PMA)$y
df_pt$smoothed_p25 <- predict(smoothed_p25_pt, df_pt$PMA)$y
df_pt$smoothed_p50 <- predict(smoothed_p50_pt, df_pt$PMA)$y
df_pt$smoothed_p75 <- predict(smoothed_p75_pt, df_pt$PMA)$y
df_pt$smoothed_p95 <- predict(smoothed_p95_pt, df_pt$PMA)$y

# Compute smoothed percentiles for full-term group
smoothed_p5_ft <- smooth.spline(df_ft$PMA, df_ft$percentiles_5)
smoothed_p25_ft <- smooth.spline(df_ft$PMA, df_ft$percentiles_25)
smoothed_p50_ft <- smooth.spline(df_ft$PMA, df_ft$percentiles_50)
smoothed_p75_ft <- smooth.spline(df_ft$PMA, df_ft$percentiles_75)
smoothed_p95_ft <- smooth.spline(df_ft$PMA, df_ft$percentiles_95)

# Extract the smoothed percentiles for full-term group
df_ft$smoothed_p5 <- predict(smoothed_p5_ft, df_ft$PMA)$y
df_ft$smoothed_p25 <- predict(smoothed_p25_ft, df_ft$PMA)$y
df_ft$smoothed_p50 <- predict(smoothed_p50_ft, df_ft$PMA)$y
df_ft$smoothed_p75 <- predict(smoothed_p75_ft, df_ft$PMA)$y
df_ft$smoothed_p95 <- predict(smoothed_p95_ft, df_ft$PMA)$y

# Combine the two data frames (preterm and full-term) into one
smoothed_percentiles_combined <- rbind(df_pt, df_ft)

# Save the combined data to CSV
#write.csv(smoothed_percentiles_combined, "smoothed_percentiles_combined.csv", row.names = FALSE)


save_dir <- paste0("./saved_GAMMs_K3/V", i)

dir.create(save_dir)
# Optionally, save preterm and full-term groups separately as well
write.csv(df_pt, file.path(save_dir,"smoothed_percentiles_pt.csv"), row.names = FALSE)  # Preterm group
write.csv(df_ft, file.path(save_dir,"smoothed_percentiles_ft.csv"), row.names = FALSE)  # Full-term group

write.csv(GAMM_data, file.path(save_dir,"GAMM_example.csv"), row.names = FALSE)

# save p-value also
model_summary <- summary(model$gam)
parametric_p_values <- model_summary$p.table[, "Pr(>|t|)"]
smooth_p_values <- model_summary$s.table[, "p-value"]
# Create a data frame with the p-value
p_values_df <- data.frame(
  Term = c("Intercept", "Group1", "s(PMA):Group0",  "s(PMA)"),
  p_value = c(parametric_p_values[1], parametric_p_values[2], smooth_p_values[1], smooth_p_values[2])
)
# Export the p-value to a CSV file
write.csv(p_values_df, file.path(save_dir,"GAMM_group_pvalue.csv"),, row.names = FALSE)

}

