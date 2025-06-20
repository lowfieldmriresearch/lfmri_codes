# Define SEM model (simplified for clarity)
library(lavaan)
library(dplyr)


# load prepared data,
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
raw_data_root = "../03_serum_features/mat_for_nutr_analysis_rawvolume.csv"
rawd <- read.csv(raw_data_root)
# rearrange data parts,
mri_res <- rawd[, c(2:48)]
lab_res <- rawd[, c(49:84)]

# Sum up all GrayMatters and all Ventricles,
rawd$GrayMatter <- rowSums(rawd[, c(15:48)])
rawd <- rawd %>%
  mutate(AllVentricle = Lateral.Ventricle + Third.Ventricle + Fourth.Ventricle)  # Summing the columns

sem_model <- '
  # Mediators regressed on exposure
  GA   ~ a0*pt.0
  V.D  ~ a1*pt.0
  Albumin ~ a2*pt.0
  Creatinine ~ a3*pt.0
  CysC  ~ a4*pt.0
  

  # Each brain region regressed on mediators and exposure
  AllVentricle    ~ c1*pt.0 + b01*GA + b11*V.D + b12*Albumin + b13*Creatinine + b14*CysC
  Corpus.Callosum  ~ c2*pt.0 + b02*GA + b21*V.D + b22*Albumin + b23*Creatinine + b24*CysC
  Amyglada        ~ c3*pt.0 + b03*GA + b31*V.D + b32*Albumin + b33*Creatinine + b34*CysC
  White.Matter     ~ c4*pt.0 + b04*GA + b41*V.D + b42*Albumin + b43*Creatinine + b44*CysC
  Hippocampus     ~ c5*pt.0 + b05*GA + b51*V.D + b52*Albumin + b53*Creatinine + b54*CysC
  Thalamus        ~ c6*pt.0 + b06*GA + b61*V.D + b62*Albumin + b63*Creatinine + b64*CysC
  Caudate         ~ c7*pt.0 + b07*GA + b71*V.D + b72*Albumin + b73*Creatinine + b74*CysC
  Putamen         ~ c8*pt.0 + b08*GA + b81*V.D + b82*Albumin + b83*Creatinine + b84*CysC
  Pallidum        ~ c9*pt.0 + b09*GA + b91*V.D + b92*Albumin + b93*Creatinine + b94*CysC
  Accumbens       ~ c10*pt.0 + b010*GA + b101*V.D + b102*Albumin + b103*Creatinine + b104*CysC
  Brainstem       ~ c11*pt.0 + b011*GA + b111*V.D + b112*Albumin + b113*Creatinine + b114*CysC
  GrayMatter      ~ c12*pt.0 + b012*GA + b121*V.D + b122*Albumin + b123*Creatinine + b124*CysC



# Indirect effects
  ind_AllVentricle_GA         := a0 * b01
  ind_AllVentricle_VD         := a1 * b11
  ind_AllVentricle_Albumin    := a2 * b12
  ind_AllVentricle_Creatinine := a3 * b13
  ind_AllVentricle_CysC       := a4 * b14

  ind_CorpusCallosum_GA         := a0 * b02
  ind_CorpusCallosum_VD         := a1 * b21
  ind_CorpusCallosum_Albumin    := a2 * b22
  ind_CorpusCallosum_Creatinine := a3 * b23
  ind_CorpusCallosum_CysC       := a4 * b24

  ind_Amyglada_GA         := a0 * b03
  ind_Amyglada_VD         := a1 * b31
  ind_Amyglada_Albumin    := a2 * b32
  ind_Amyglada_Creatinine := a3 * b33
  ind_Amyglada_CysC       := a4 * b34

  ind_WhiteMatter_GA         := a0 * b04
  ind_WhiteMatter_VD         := a1 * b41
  ind_WhiteMatter_Albumin    := a2 * b42
  ind_WhiteMatter_Creatinine := a3 * b43
  ind_WhiteMatter_CysC       := a4 * b44

  ind_Hippocampus_GA         := a0 * b05
  ind_Hippocampus_VD         := a1 * b51
  ind_Hippocampus_Albumin    := a2 * b52
  ind_Hippocampus_Creatinine := a3 * b53
  ind_Hippocampus_CysC       := a4 * b54

  ind_Thalamus_GA         := a0 * b06
  ind_Thalamus_VD         := a1 * b61
  ind_Thalamus_Albumin    := a2 * b62
  ind_Thalamus_Creatinine := a3 * b63
  ind_Thalamus_CysC       := a4 * b64

  ind_Caudate_GA         := a0 * b07
  ind_Caudate_VD         := a1 * b71
  ind_Caudate_Albumin    := a2 * b72
  ind_Caudate_Creatinine := a3 * b73
  ind_Caudate_CysC       := a4 * b74

  ind_Putamen_GA         := a0 * b08
  ind_Putamen_VD         := a1 * b81
  ind_Putamen_Albumin    := a2 * b82
  ind_Putamen_Creatinine := a3 * b83
  ind_Putamen_CysC       := a4 * b84

  ind_Pallidum_GA         := a0 * b09
  ind_Pallidum_VD         := a1 * b91
  ind_Pallidum_Albumin    := a2 * b92
  ind_Pallidum_Creatinine := a3 * b93
  ind_Pallidum_CysC       := a4 * b94

  ind_Accumbens_GA         := a0 * b010
  ind_Accumbens_VD         := a1 * b101
  ind_Accumbens_Albumin    := a2 * b102
  ind_Accumbens_Creatinine := a3 * b103
  ind_Accumbens_CysC       := a4 * b104
  
  ind_Brainstem_GA         := a0 * b011
  ind_Brainstem_VD         := a1 * b111
  ind_Brainstem_Albumin    := a2 * b112
  ind_Brainstem_Creatinine := a3 * b113
  ind_Brainstem_CysC       := a4 * b114

  ind_GrayMatter_GA         := a0 * b012
  ind_GrayMatter_VD         := a1 * b121
  ind_GrayMatter_Albumin    := a2 * b122
  ind_GrayMatter_Creatinine := a3 * b123
  ind_GrayMatter_CysC       := a4 * b124
'



fit <-lavaan::sem(sem_model,data = rawd,std.ov=T,std.lv=T)
summary(fit, fit.measures = TRUE)

# visualize method1
library(semPlot)
# 可视化路径图
semPaths(
  fit,
  what = "std",             # standardized estimates
  edge.label.cex = 0.8,
  layout = "tree",          # or "circle", "spring", etc.
  cut = 0.05,                # only show paths with |estimate| > 0.2
  nCharNodes = 0,           # full variable names
  sizeMan = 6, sizeLat = 8, # control node sizes
  mar = c(5, 5, 5, 5)
)


# Get standardized parameter estimates
est_std <- parameterEstimates(fit, standardized = TRUE)

# Keep all regression paths (op == "~")
all_effects <- est_std %>%
  filter(op == "~") %>%
  select(rhs, lhs, est, std.all, pvalue) %>%
  arrange(desc(abs(std.all)))  # Sort by magnitude of standardized effect

# Optionally rename columns
colnames(all_effects) <- c("Predictor", "Outcome", "Estimate", "Standardized", "p")

# View as table
print(all_effects)
# Save data frame 'data' to a CSV file
write.csv(all_effects, "full_effect.csv", row.names = FALSE)  # Set row.names to FALSE if you don't want row indices

