# Install and load the necessary packages
library(itsadug)
library(tidyverse)


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
df <- read.csv("creatinine_vs_aveBV_for_GAMM.csv",header = F)
df <- subset(df, V2<65)

df_pt <- df %>% filter(V1 == 0)  # Preterm data
df_ft <- df %>% filter(V1 == 1)  # Full-term data


#mri_res <- dd[,c("V2")]
#lab_res <- dd[,c("V1")]


ft_gam <- gam(V3 ~ s(V2,k=3), data = df_ft)
pt_gam <- gam(V3 ~ s(V2,k=3), data = df_pt)


#visualized the fitted non-parametric curve by GAMM:
par(mfrow = c(1, 2))
pt <- plot_smooth(pt_gam, view = "V2", rm.ranef = TRUE, main = "Brain Volume vs PMA",n.grid = 100)
ft <- plot_smooth(ft_gam, view = "V2", rm.ranef = TRUE, main = "Brain Volume vs PMA",n.grid = 100)


# save the fitted curve and range.
write.csv(ft$fv, "creatinine_fitted_GAMM_curve_for_ft.csv", row.names = FALSE)
write.csv(pt$fv, "creatinine_fitted_GAMM_curve_for_pt.csv", row.names = FALSE)


