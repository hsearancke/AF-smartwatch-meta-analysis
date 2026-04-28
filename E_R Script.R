# ===============================
# 1. Packages
# ===============================
required_packages <- c("mada", "meta", "metafor", "dplyr", "ggplot2", "writexl")
for (pkg in required_packages) 
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
# ===============================
# 2. Package versions and session info
# ===============================

packageVersion("mada")
packageVersion("meta")
packageVersion("metafor")
packageVersion("dplyr")
packageVersion("ggplot2")
packageVersion("writexl")

sessionInfo()

# ===============================
# 3. Hardcoded extracted 2x2 data
# ===============================
data <- data.frame(
  StudyID = c("AbuAlrub2022","AbuAlrub2022","AbuAlrub2022",
              "Mannhardt2023","Mannhardt2023","Mannhardt2023","Mannhardt2023",
              "Muller2024","Briosa2025","OscaAsensi2021","Wouters2025","Racine2022"),
  Device = c("Apple","Samsung","Withings",
             "Apple","Fitbit","Samsung","Withings",
             "Apple","Apple","RITHMI","Apple","Apple"),
  TP = c(87,88,78,52,52,35,40,13,134,92,30,120),
  FN = c(13,12,22,9,9,26,21,5,60,6,0,34),
  FP = c(14,19,20,35,35,35,29,6,79,1,2,98),
  TN = c(86,81,80,105,105,105,111,69,211,34,90,482)
)

#===============================
# SCRIPT: 1_Bivariate_Model.R
# ===============================
# Step 1: Install package (only once)
install.packages("mada")
# Step 2: Load package
library(mada)
# Step 3: Create your dataset
data <- data.frame(
  StudyID = c("AbuAlrub2022","AbuAlrub2022","AbuAlrub2022",
              "Mannhardt2023","Mannhardt2023","Mannhardt2023","Mannhardt2023",
              "Muller2024","BrisessionInfo()
osa2025","OscaAsensi2021","Wouters2025","Racine2022"),
  Device = c("Apple","Samsung","Withings",
             "Apple","Fitbit","Samsung","Withings",
             "Apple","Apple","RITHMI","Apple","Apple"),
  TP = c(87,88,78,52,52,35,40,13,134,92,30,120),
  FN = c(13,12,22,9,9,26,21,5,60,6,0,34),
  FP = c(14,19,20,35,35,35,29,6,79,1,2,98),
  TN = c(86,81,80,105,105,105,111,69,211,34,90,482)
)
# Step 4: Run the model
model <- reitsma(data)
# Step 5: PRINT RESULTS 
summary(model)
packageVersion("mada")
fpr <- 0.175
fpr_lower <- 0.133
fpr_upper <- 0.227
specificity <- 1 - fpr
spec_lower <- 1 - fpr_upper
spec_upper <- 1 - fpr_lower
specificity
spec_lower
spec_upper
# Step 6: Save SROC Plot at EMF 
# Install package (only once)
install.packages("mada")
# Install if needed
install.packages("devEMF")
# Load libraries
library(mada)
library(devEMF)
# Open EMF device
emf("SROC_plot.emf", width = 7, height = 7)
# Fix margins (prevents "figure margins too large")
par(mar = c(4, 4, 3, 1))
# Plot SROC curve
plot(model, sroclwd = 2)
# Add study points
points(
1 - (data$TN / (data$TN + data$FP)),  # FPR
data$TP / (data$TP + data$FN),        # Sensitivity
pch = 16,
col = "blue"
)
# Close device
dev.off()

# ===============================
# SCRIPT: 2_HSROC_Model.R
# ===============================
library(mada)
library(dplyr)
packageVersion("mada")
packageVersion("dplyr")
data <- data.frame(
StudyID = c("AbuAlrub2022","AbuAlrub2022","AbuAlrub2022",
"Mannhardt2023","Mannhardt2023","Mannhardt2023","Mannhardt2023",
"Muller2024","Briosa2025","OscaAsensi2021","Wouters2025","Racine2022"),
TP = c(87,88,78,52,52,35,40,13,134,92,30,120),
FN = c(13,12,22,9,9,26,21,5,60,6,0,34),
FP = c(14,19,20,35,35,35,29,6,79,1,2,98),
TN = c(86,81,80,105,105,105,111,69,211,34,90,482)
)
data_study <- data %>%
group_by(StudyID) %>%
summarise(
TP = sum(TP),
FN = sum(FN),
FP = sum(FP),
TN = sum(TN),
.groups = "drop"
)
# Apply 0.5 continuity correction after study-level aggregation
data_study_cc <- data_study
data_study_cc$TP <- data_study_cc$TP + 0.5
data_study_cc$FN <- data_study_cc$FN + 0.5
data_study_cc$FP <- data_study_cc$FP + 0.5
data_study_cc$TN <- data_study_cc$TN + 0.5
# Create mada matrix using corrected data
mada_matrix <- as.matrix(data_study_cc[, c("TP", "FN", "FP", "TN")])
# Fit HSROC/Reitsma model
hsroc_model <- reitsma(mada_matrix)
summary(hsroc_model)
plot(hsroc_model)
png("hsroc_plot.png", width = 1200, height = 1200, res = 150)
par(mar = c(4, 4, 3, 1))
plot(hsroc_model, sroclwd = 2)
points(
1 - (data_study$TN / (data_study$TN + data_study$FP)),
data_study$TP / (data_study$TP + data_study$FN),
pch = 16,
col = "blue"
)
dev.off()
# Look at coefficient structure
coef(hsroc_model)
# Extract values by column position/name
sens_logit <- coef(hsroc_model)[1, "tsens"]
fpr_logit  <- coef(hsroc_model)[1, "tfpr"]
# Convert from logit scale
sensitivity <- plogis(sens_logit)
fpr <- plogis(fpr_logit)
# Specificity
specificity <- 1 - fpr
sensitivity
specificity
s <- summary(hsroc_model)
# Sensitivity CI directly from summary
sens <- s$coefficients["sensitivity", "Estimate"]
sens_lower <- s$coefficients["sensitivity", "95%ci.lb"]
sens_upper <- s$coefficients["sensitivity", "95%ci.ub"]
# FPR CI from summary
fpr <- s$coefficients["false pos. rate", "Estimate"]
fpr_lower <- s$coefficients["false pos. rate", "95%ci.lb"]
fpr_upper <- s$coefficients["false pos. rate", "95%ci.ub"]
# Convert FPR to specificity
specificity <- 1 - fpr
spec_lower <- 1 - fpr_upper
spec_upper <- 1 - fpr_lower
sens
sens_lower
sens_upper
specificity
spec_lower
spec_upper


# ===============================
# SCRIPT: 3_DOR.R
# ===============================
# Install packages if needed
install.packages("meta")
# Load libraries
library(meta)
# Create dataset
data <- data.frame(
StudyID = c("AbuAlrub2022","AbuAlrub2022","AbuAlrub2022",
"Mannhardt2023","Mannhardt2023","Mannhardt2023","Mannhardt2023",
"Muller2024","Briosa2025","OscaAsensi2021","Wouters2025","Racine2022"),
Device = c("Apple","Samsung","Withings",
"Apple","Fitbit","Samsung","Withings",
"Apple","Apple","RITHMI","Apple","Apple"),
TP = c(87,88,78,52,52,35,40,13,134,92,30,120),
FN = c(13,12,22,9,9,26,21,5,60,6,0,34),
FP = c(14,19,20,35,35,35,29,6,79,1,2,98),
TN = c(86,81,80,105,105,105,111,69,211,34,90,482)
)
# Apply continuity correction of 0.5 to all cells
data_cc <- data
data_cc$TP <- data$TP + 0.5
data_cc$FN <- data$FN + 0.5
data_cc$FP <- data$FP + 0.5
data_cc$TN <- data$TN + 0.5
# Calculate DOR
data_cc$DOR <- (data_cc$TP * data_cc$TN) / (data_cc$FP * data_cc$FN)
# Log-transform DOR
data_cc$logDOR <- log(data_cc$DOR)
# Standard error of log(DOR)
data_cc$SE_logDOR <- sqrt(
(1 / data_cc$TP) +
(1 / data_cc$FN) +
(1 / data_cc$FP) +
(1 / data_cc$TN)
)
# Run random-effects meta-analysis using DerSimonian-Laird estimator
dor_meta <- metagen(
TE = data_cc$logDOR,
seTE = data_cc$SE_logDOR,
studlab = data_cc$StudyID,
sm = "OR",
method.tau = "DL"
)
# View results
summary(dor_meta)
# Optional forest plot
forest(dor_meta)

# ===============================
# SCRIPT: 4_PPV_MetaAnlayisis.R
# ===============================
install.packages("meta")
packageVersion("meta")
# Load library
library(meta)
# Create dataset
data <- data.frame(
StudyID = c("AbuAlrub2022","AbuAlrub2022","AbuAlrub2022",
"Mannhardt2023","Mannhardt2023","Mannhardt2023","Mannhardt2023",
"Muller2024","Briosa2025","OscaAsensi2021","Wouters2025","Racine2022"),
Device = c("Apple","Samsung","Withings",
"Apple","Fitbit","Samsung","Withings",
"Apple","Apple","RITHMI","Apple","Apple"),
TP = c(87,88,78,52,52,35,40,13,134,92,30,120),
FN = c(13,12,22,9,9,26,21,5,60,6,0,34),
FP = c(14,19,20,35,35,35,29,6,79,1,2,98),
TN = c(86,81,80,105,105,105,111,69,211,34,90,482)
)
# Calculate PPV per study/device
data$PPV <- data$TP / (data$TP + data$FP)
# Run random-effects meta-analysis of PPV
ppv_meta <- metaprop(
event = TP,
n = TP + FP,
studlab = StudyID,
data = data,
sm = "PLOGIT",
method.tau = "ML"
)
# View results
summary(ppv_meta)
# Optional forest plot
forest(ppv_meta)

# ===============================
# SCRIPT: 5_NPV_MetaAnalysis.R
# ===============================
install.packages("meta")
packageVersion("meta")
# Load library
library(meta)
# Create dataset
data <- data.frame(
StudyID = c("AbuAlrub2022","AbuAlrub2022","AbuAlrub2022",
"Mannhardt2023","Mannhardt2023","Mannhardt2023","Mannhardt2023",
"Muller2024","Briosa2025","OscaAsensi2021","Wouters2025","Racine2022"),
Device = c("Apple","Samsung","Withings",
"Apple","Fitbit","Samsung","Withings",
"Apple","Apple","RITHMI","Apple","Apple"),
TP = c(87,88,78,52,52,35,40,13,134,92,30,120),
FN = c(13,12,22,9,9,26,21,5,60,6,0,34),
FP = c(14,19,20,35,35,35,29,6,79,1,2,98),
TN = c(86,81,80,105,105,105,111,69,211,34,90,482)
)
# Calculate NPV per study/device
data$NPV <- data$TN / (data$TN + data$FN)
# Run random-effects meta-analysis of NPV
npv_meta <- metaprop(
event = TN,
n = TN + FN,
studlab = StudyID,
data = data,
sm = "PLOGIT",
method.tau = "ML"
)
# View results
summary(npv_meta)
# Optional forest plot
forest(npv_meta)

# ===============================
# SCRIPT: 6_DOR_MetaAnalysis.R
# ===============================
install.packages("metafor")
# Load library
library(metafor)
packageVersion("metafor")
# Create dataset
data <- data.frame(
StudyID = c("AbuAlrub2022","AbuAlrub2022","AbuAlrub2022",
"Mannhardt2023","Mannhardt2023","Mannhardt2023","Mannhardt2023",
"Muller2024","Briosa2025","OscaAsensi2021","Wouters2025","Racine2022"),
Device = c("Apple","Samsung","Withings",
"Apple","Fitbit","Samsung","Withings",
"Apple","Apple","RITHMI","Apple","Apple"),
TP = c(87,88,78,52,52,35,40,13,134,92,30,120),
FN = c(13,12,22,9,9,26,21,5,60,6,0,34),
FP = c(14,19,20,35,35,35,29,6,79,1,2,98),
TN = c(86,81,80,105,105,105,111,69,211,34,90,482)
)
# Apply continuity correction
data_cc <- data
data_cc$TP <- data$TP + 0.5
data_cc$FN <- data$FN + 0.5
data_cc$FP <- data$FP + 0.5
data_cc$TN <- data$TN + 0.5
# Calculate DOR and log(DOR)
data_cc$DOR <- (data_cc$TP * data_cc$TN) / (data_cc$FP * data_cc$FN)
data_cc$logDOR <- log(data_cc$DOR)
# Calculate SE and variance of log(DOR)
data_cc$SE_logDOR <- sqrt(
(1 / data_cc$TP) +
(1 / data_cc$FN) +
(1 / data_cc$FP) +
(1 / data_cc$TN)
)
data_cc$vi <- data_cc$SE_logDOR^2
data_cc$yi <- data_cc$logDOR
# Run multilevel random-effects model
res_ml <- rma.mv(
yi = yi,
V = vi,
random = ~ 1 | StudyID/Device,
data = data_cc,
method = "REML"
)
# View model results
summary(res_ml)
# Back-transform pooled log(DOR) to DOR
exp(res_ml$b)
# Back-transform 95% CI to DOR scale
exp(c(res_ml$ci.lb, res_ml$ci.ub))

# ===============================
# SCRIPT: 7_Device_Meta_Regression.R
# ===============================
# Install packages if needed
install.packages("metafor")
# Load libraries
library(metafor)
library(dplyr)
# Create dataset
data <- data.frame(
StudyID = c("AbuAlrub2022","AbuAlrub2022","AbuAlrub2022",
"Mannhardt2023","Mannhardt2023","Mannhardt2023","Mannhardt2023",
"Muller2024","Briosa2025","OscaAsensi2021","Wouters2025","Racine2022"),
Device = c("Apple","Samsung","Withings",
"Apple","Fitbit","Samsung","Withings",
"Apple","Apple","RITHMI","Apple","Apple"),
TP = c(87,88,78,52,52,35,40,13,134,92,30,120),
FN = c(13,12,22,9,9,26,21,5,60,6,0,34),
FP = c(14,19,20,35,35,35,29,6,79,1,2,98),
TN = c(86,81,80,105,105,105,111,69,211,34,90,482)
)
# -------------------------------
# Step 1: Remove devices used only once
# -------------------------------
data_meta <- data %>%
filter(Device != "RITHMI", Device != "Fitbit")
# -------------------------------
# Step 2: Continuity correction
# -------------------------------
data_meta <- data_meta %>%
mutate(
TP = ifelse(TP == 0, TP + 0.5, TP),
FP = ifelse(FP == 0, FP + 0.5, FP),
FN = ifelse(FN == 0, FN + 0.5, FN),
TN = ifelse(TN == 0, TN + 0.5, TN)
)
# -------------------------------
# Step 3: Calculate DOR
# -------------------------------
data_meta <- data_meta %>%
mutate(
DOR = (TP * TN) / (FP * FN),
logDOR = log(DOR),
var_logDOR = 1/TP + 1/FN + 1/FP + 1/TN
)
# -------------------------------
# Step 4: Run meta-regression
# -------------------------------
meta_model <- rma(
yi = logDOR,
vi = var_logDOR,
mods = ~ factor(Device),
method = "REML",
data = data_meta
)
# -------------------------------
# Step 5: View results
# -------------------------------
summary(meta_model)

# ===============================
# SCRIPT: 8_Deeks_Publication_Bias.R
# ===============================
library(metafor)
# Continuity correction
data_cc <- data
data_cc$TP <- data$TP + 0.5
data_cc$FN <- data$FN + 0.5
data_cc$FP <- data$FP + 0.5
data_cc$TN <- data$TN + 0.5
# DOR and log(DOR)
data_cc$DOR <- (data_cc$TP * data_cc$TN) / (data_cc$FP * data_cc$FN)
data_cc$logDOR <- log(data_cc$DOR)
# Standard error of log(DOR)
data_cc$SE_logDOR <- sqrt(
(1 / data_cc$TP) +
(1 / data_cc$FN) +
(1 / data_cc$FP) +
(1 / data_cc$TN)
)
# Effective sample size for Deeks test
data_cc$ESS <- 4 / (
(1 / data_cc$TP) +
(1 / data_cc$FN) +
(1 / data_cc$FP) +
(1 / data_cc$TN)
)
# Predictor for Deeks test
data_cc$inv_sqrt_ESS <- 1 / sqrt(data_cc$ESS)
# Deeks meta-regression
deeks_model <- rma(
yi = logDOR,
vi = SE_logDOR^2,
mods = ~ inv_sqrt_ESS,
method = "REML",
data = data_cc
)
summary(deeks_model)
# Deeks funnel plot
# Deeks funnel plot saved as EMF
library(devEMF)

dir.create("Outputs/Figures", recursive = TRUE, showWarnings = FALSE)

emf("Outputs/Figures/deeks_funnel_plot.emf", width = 7, height = 7)

par(mar = c(5, 5, 4, 2))

plot(
  data_cc$inv_sqrt_ESS,
  data_cc$logDOR,
  xlab = "1 / sqrt(Effective Sample Size)",
  ylab = "Log Diagnostic Odds Ratio",
  pch = 16,
  main = "Deeks' Funnel Plot"
)

abline(deeks_model, col = "red", lwd = 2)

dev.off()

# ===============================
# SCRIPT: 9_Device_Subgroup_Analysis.R
# ===============================
# Function to extract estimates + CI from a reitsma model
extract_results <- function(model) {
s <- summary(model)
# Extract logit estimates
tsens <- coef(model)[1]
tfpr  <- coef(model)[2]
# Extract CI (logit scale)
tsens_lower <- s$coef["tsens.(Intercept)", "95%ci.lb"]
tsens_upper <- s$coef["tsens.(Intercept)", "95%ci.ub"]
tfpr_lower <- s$coef["tfpr.(Intercept)", "95%ci.lb"]
tfpr_upper <- s$coef["tfpr.(Intercept)", "95%ci.ub"]
# Convert to probability scale
sens <- plogis(tsens)
sens_lower <- plogis(tsens_lower)
sens_upper <- plogis(tsens_upper)
fpr <- plogis(tfpr)
spec <- 1 - fpr
spec_lower <- 1 - plogis(tfpr_upper)
spec_upper <- 1 - plogis(tfpr_lower)
return(c(sens, sens_lower, sens_upper, spec, spec_lower, spec_upper))
}
apple_vals    <- extract_results(model_apple)
samsung_vals  <- extract_results(model_samsung)
withings_vals <- extract_results(model_withings)
results_devices <- data.frame(
Device = c("Apple", "Samsung", "Withings"),
Sensitivity = c(apple_vals[1], samsung_vals[1], withings_vals[1]),
Sens_Lower  = c(apple_vals[2], samsung_vals[2], withings_vals[2]),
Sens_Upper  = c(apple_vals[3], samsung_vals[3], withings_vals[3]),
Specificity = c(apple_vals[4], samsung_vals[4], withings_vals[4]),
Spec_Lower  = c(apple_vals[5], samsung_vals[5], withings_vals[5]),
Spec_Upper  = c(apple_vals[6], samsung_vals[6], withings_vals[6])
)
print(results_devices)
results_table <- data.frame(
Device = results_devices$Device,
Sensitivity = sprintf(
"%.2f (%.2f–%.2f)",
results_devices$Sensitivity,
results_devices$Sens_Lower,
results_devices$Sens_Upper
),
Specificity = sprintf(
"%.2f (%.2f–%.2f)",
results_devices$Specificity,
results_devices$Spec_Lower,
results_devices$Spec_Upper
)
)
results_table
colnames(results_table) <- c(
"Device",
"Sensitivity (95% CI)",
"Specificity (95% CI)"
)
install.packages("writexl")
library(writexl)
write_xlsx(results_table, "device_results_table.xlsx")

#Device Specific Bar Graph
install.packages("metafor")
library(mada)
library(ggplot2)
data_apple <- subset(data, Device == "Apple")
data_samsung <- subset(data, Device == "Samsung")
data_withings <- subset(data, Device == "Withings")
model_apple <- reitsma(data_apple)
model_samsung <- reitsma(data_samsung)
model_withings <- reitsma(data_withings)
# Apple
sens_apple <- plogis(coef(model_apple)[1])
spec_apple <- 1 - plogis(coef(model_apple)[2])
# Samsung
sens_samsung <- plogis(coef(model_samsung)[1])
spec_samsung <- 1 - plogis(coef(model_samsung)[2])
# Withings
sens_withings <- plogis(coef(model_withings)[1])
spec_withings <- 1 - plogis(coef(model_withings)[2])
sens_apple; spec_apple
results_devices <- data.frame(
  Device = c("Apple", "Samsung", "Withings"),
  Sensitivity = c(sens_apple, sens_samsung, sens_withings),
  Specificity = c(spec_apple, spec_samsung, spec_withings)
)
results_devices
plot_data <- data.frame(
  Device = rep(results_devices$Device, 2),
  Metric = rep(c("Sensitivity", "Specificity"), each = 3),
  Value = c(results_devices$Sensitivity, results_devices$Specificity)
)
ggplot(plot_data, aes(x = Device, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Sensitivity" = "red", "Specificity" = "blue")) +
  ylim(0, 1) +
  labs(
    x = "Device",
    y = "Value",
    title = "Sensitivity and Specificity by Device"
  ) +
  theme_minimal()
# -------------------------------
# Device-specific HSROC/SROC plots with 95% confidence ellipses
# Saved as Windows Metafile (.wmf)
# -------------------------------
library(mada)
# Split data by device
data_apple <- subset(data, Device == "Apple")
data_samsung <- subset(data, Device == "Samsung")
data_withings <- subset(data, Device == "Withings")
# Fit Reitsma models
model_apple <- reitsma(data_apple)
model_samsung <- reitsma(data_samsung)
model_withings <- reitsma(data_withings)
# Save as Windows Metafile
win.metafile("device_hsroc_curves.wmf", width = 10, height = 4)
# Three plots side-by-side
par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
# Apple
plot(model_apple, main = "Apple", level = 0.95, sroclwd = 2)
points(
  1 - (data_apple$TN / (data_apple$TN + data_apple$FP)),
  data_apple$TP / (data_apple$TP + data_apple$FN),
  pch = 16,
  col = "blue"
)
# Samsung
plot(model_samsung, main = "Samsung", level = 0.95, sroclwd = 2)
points(
  1 - (data_samsung$TN / (data_samsung$TN + data_samsung$FP)),
  data_samsung$TP / (data_samsung$TP + data_samsung$FN),
  pch = 16,
  col = "blue"
)
# Withings
plot(model_withings, main = "Withings", level = 0.95, sroclwd = 2)
points(
  1 - (data_withings$TN / (data_withings$TN + data_withings$FP)),
  data_withings$TP / (data_withings$TP + data_withings$FN),
  pch = 16,
  col = "blue"
)
# Close file
dev.off()
getwd()
# ===============================
# SCRIPT: 10_Sensitivity_Analysis.R
# ===============================
# Worst-case vs discarded data handling analysis
# -------------------------------
# Load libraries
library(dplyr)
library(meta)
library(mada)
library(ggplot2)
# -------------------------------
# 1. Create grouping variable
# -------------------------------
data_grouped <- data %>%
  mutate(
    Handling = case_when(
      StudyID %in% c("Racine2022", "Briosa2025", "Mannhardt2023", "AbuAlrub2022") ~ "Worst-case",
      StudyID %in% c("Muller2024", "Wouters2025") ~ "Discarded",
      TRUE ~ NA_character_
    )
  )
# Exclude studies with unclear handling, e.g. OscaAsensi2021
data_grouped_clean <- data_grouped %>%
  filter(!is.na(Handling))
# Check included studies
table(data_grouped_clean$Handling)
data_grouped_clean[, c("StudyID", "Device", "Handling")]
# -------------------------------
# 2. Sensitivity meta-analysis by handling method
# -------------------------------
sens_handling <- metaprop(
  event = TP,
  n = TP + FN,
  studlab = StudyID,
  data = data_grouped_clean,
  sm = "PLOGIT",
  method.tau = "ML",
  subgroup = Handling
)
summary(sens_handling)
# -------------------------------
# 3. Specificity meta-analysis by handling method
# -------------------------------
spec_handling <- metaprop(
  event = TN,
  n = TN + FP,
  studlab = StudyID,
  data = data_grouped_clean,
  sm = "PLOGIT",
  method.tau = "ML",
  subgroup = Handling
)
summary(spec_handling)
# -------------------------------
# 4. Extract pooled subgroup estimates
# -------------------------------
logit_to_prob <- function(x) exp(x) / (1 + exp(x))
handling_results <- data.frame(
  Handling = levels(factor(data_grouped_clean$Handling)),
  Sensitivity = logit_to_prob(sens_handling$TE.random.w),
  Sens_Lower = logit_to_prob(sens_handling$lower.random.w),
  Sens_Upper = logit_to_prob(sens_handling$upper.random.w),
  Specificity = logit_to_prob(spec_handling$TE.random.w),
  Spec_Lower = logit_to_prob(spec_handling$lower.random.w),
  Spec_Upper = logit_to_prob(spec_handling$upper.random.w)
)
print(handling_results)
# -------------------------------
# 5. Bar graph with 95% CI
# -------------------------------
plot_data <- data.frame(
  Handling = rep(handling_results$Handling, 2),
  Metric = rep(c("Sensitivity", "Specificity"), each = nrow(handling_results)),
  Value = c(handling_results$Sensitivity, handling_results$Specificity),
  Lower = c(handling_results$Sens_Lower, handling_results$Spec_Lower),
  Upper = c(handling_results$Sens_Upper, handling_results$Spec_Upper)
)
ggplot(plot_data, aes(x = Handling, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_errorbar(
    aes(ymin = Lower, ymax = Upper),
    position = position_dodge(width = 0.9),
    width = 0.2
  ) +
  scale_fill_manual(values = c("Sensitivity" = "red", "Specificity" = "blue")) +
  ylim(0, 1) +
  labs(
    x = "Data handling method",
    y = "Value",
    title = "Sensitivity and Specificity by Data Handling Method"
  ) +
  theme_minimal()
win.metafile("handling_bar_plot.wmf", width = 7, height = 5)
ggplot(plot_data, aes(x = Handling, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_errorbar(
    aes(ymin = Lower, ymax = Upper),
    position = position_dodge(width = 0.9),
    width = 0.2
  ) +
  scale_fill_manual(values = c("Sensitivity" = "red", "Specificity" = "blue")) +
  ylim(0, 1) +
  labs(
    x = "Data handling method",
    y = "Value",
    title = "Sensitivity and Specificity by Data Handling Method"
  ) +
  theme_minimal()
dev.off()
# -------------------------------
# HSROC / SROC curve with study points
# Save directly to file
# -------------------------------

library(mada)
library(devEMF)

handling_model <- reitsma(data_grouped_clean)

# Make sure FPR and Sensitivity columns exist
data_grouped_clean$FPR <- data_grouped_clean$FP /
  (data_grouped_clean$FP + data_grouped_clean$TN)

data_grouped_clean$Sensitivity <- data_grouped_clean$TP /
  (data_grouped_clean$TP + data_grouped_clean$FN)

dir.create("Outputs/Figures", recursive = TRUE, showWarnings = FALSE)

emf("Outputs/Figures/handling_hsroc_plot.emf", width = 7, height = 7)

par(mar = c(5, 5, 4, 2))

plot(
  handling_model,
  sroclwd = 2,
  main = "HSROC Curve by Data Handling Method"
)

points(
  data_grouped_clean$FPR[data_grouped_clean$Handling == "Worst-case"],
  data_grouped_clean$Sensitivity[data_grouped_clean$Handling == "Worst-case"],
  pch = 16,
  col = "red"
)

points(
  data_grouped_clean$FPR[data_grouped_clean$Handling == "Discarded"],
  data_grouped_clean$Sensitivity[data_grouped_clean$Handling == "Discarded"],
  pch = 16,
  col = "blue"
)

legend(
  "bottomright",
  legend = c("Worst-case", "Discarded"),
  col = c("red", "blue"),
  pch = 16
)

dev.off()
# ===============================
# Printed summary outputs
# ===============================
if (exists("model")) {
cat("\n--- Bivariate Reitsma Model ---\n")
print(summary(model))
}
if (exists("hsroc_model")) {
cat("\n--- HSROC Model ---\n")
print(summary(hsroc_model))
}
if (exists("dor_meta")) {
cat("\n--- DOR Meta-analysis ---\n")
print(summary(dor_meta))
}
if (exists("ppv_meta")) {
cat("\n--- PPV Meta-analysis ---\n")
print(summary(ppv_meta))
}
if (exists("npv_meta")) {
cat("\n--- NPV Meta-analysis ---\n")
print(summary(npv_meta))
}
if (exists("res_ml")) {
cat("\n--- Multilevel Model ---\n")
print(summary(res_ml))
cat("\nPooled multilevel DOR:\n")
print(exp(res_ml$b))
}
if (exists("meta_model")) {
cat("\n--- Device Type Meta-regression ---\n")
print(summary(meta_model))
}
if (exists("deeks_model")) {
cat("\n--- Deeks' Funnel Plot Asymmetry Test ---\n")
print(summary(deeks_model))
}
