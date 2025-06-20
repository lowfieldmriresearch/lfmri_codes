# Install and load the necessary packages
library(itsadug)
library(mgcv)

library(dplyr)
library(zoo)


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
dd <- read.csv("../raw_data/full_res_for_gamm00_use.csv",header = F)
region_names <- read.csv("V_var_names.csv", header = TRUE, stringsAsFactors = FALSE)[[1]]
results <- data.frame(Region = character(), P_value = numeric(), stringsAsFactors = FALSE)
p_val_res <- data.frame(Region = character(), P_value = numeric(), stringsAsFactors = FALSE)

for (i in 2:48) {
  
  region_name <- region_names[i-1]
  df <- data.frame(
    brain_volume = dd[[paste0("V", i)]], 
    PMA = dd$V86,
    Group = factor(dd$V1, levels = c(0, 1))
  )
  
  dt_pt <- df %>% filter(Group == 0)  # Preterm data
  df_ft <- df %>% filter(Group == 1)  # Full-term data
  
  
  mod_reduced <- gam(brain_volume ~ Group + s(PMA,k=3), data = df, method = "ML")
  mod_full <- gam(brain_volume ~ Group + s(PMA, by = Group,k=3), data = df, method = "ML")
  
  # Likelihood ratio test comparing the two models
  anova_res <- anova(mod_reduced, mod_full, test = "Chisq")
  
  p_val <- anova_res[2, "Pr(>Chi)"]
  
  # Determine significance level
  sig_level <- case_when(
    p_val < 0.001 ~ 4,
    p_val < 0.005 ~ 3,
    p_val < 0.01  ~ 2,
    p_val < 0.05  ~ 1,
    TRUE          ~ 0
  )
  results <- rbind(results, data.frame(Region = region_name, P_value = p_val, Sig_Level = sig_level))
  p_val_res <-rbind(p_val_res,data.frame(Region = region_name, P_value = p_val))
}

p_bh <- p.adjust(p_val_res$P_value, method = "BH")
p_bonf <- p.adjust(p_val_res$P_value, method = "bonferroni")
p_holm <- p.adjust(p_val_res$P_value, method = "holm")
p_by <- p.adjust(p_val_res$P_value, method = "BY")

results$p_bh <- p_bh
results$p_bonf <- p_bonf
results$p_holm <- p_holm
results$p_by <- p_by

write.csv(results, "ANOVA_results.csv", row.names = FALSE)
