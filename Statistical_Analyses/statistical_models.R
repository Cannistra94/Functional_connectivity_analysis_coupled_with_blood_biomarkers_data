###############################################
#  Statistical Analysis
# Looping over multiple FC and cognition columns
###############################################

# Load libraries
library(ggplot2)
library(dplyr)
library(lmtest)
library(sandwich)
library(psych)
library(mediation)
library(ggpubr)

#ispect data
data <- read.csv('data.csv', header=0)

tau_cols <- c("pTau217", "pTau181")
fc_cols <- c("functional_connectivity_change_1", "functional_connectivity_change_2", "functional_connectivity_change_3")
cog_cols <- c("cognitive_1", "cognitive_2", "cognitive_3")

cat("Data structure:\n")
str(data)
cat("\nSummary statistics:\n")
summary(data)

cat("\nDescriptive statistics for numeric variables:\n")
numeric_vars <- c(tau_cols, fc_cols, cog_cols, "age")
describe(data[, numeric_vars])

# Histograms for distributions
for(var in numeric_vars){
  p <- ggplot(data, aes_string(x = var)) +
    geom_histogram(binwidth = 0.3, fill = "skyblue", color = "black") +
    labs(title = paste("Distribution of", var)) +
    theme_minimal()
  print(p)
}


# Linear models 

lm_results <- list()

for(tau_var in tau_cols){
  for(fc_var in fc_cols){
    for(cog_var in cog_cols){
      
      cat("\n-----------------------------------\n")
      cat("Linear models for Tau:", tau_var, "FC:", fc_var, "Cognition:", cog_var, "\n")
      
      # Tau ~ group
      lm_tau <- lm(as.formula(paste(tau_var, "~ group + age + sex")), data = data)
      print(coeftest(lm_tau, vcov = vcovHC(lm_tau, type = "HC3")))
      
      # FC ~ group
      lm_fc <- lm(as.formula(paste(fc_var, "~ group + age + sex")), data = data)
      print(coeftest(lm_fc, vcov = vcovHC(lm_fc, type = "HC3")))
      
      # Cognition ~ group
      lm_cog <- lm(as.formula(paste(cog_var, "~ group + age + sex")), data = data)
      print(coeftest(lm_cog, vcov = vcovHC(lm_cog, type = "HC3")))
      
      # Optional correlation plot
      p <- ggscatter(data, x = fc_var, y = cog_var,
                     add = "reg.line", conf.int = TRUE,
                     cor.coef = TRUE, cor.method = "pearson") +
        labs(title = paste("Correlation:", fc_var, "vs", cog_var))
      print(p)
      
      # Store results
      lm_results[[paste(tau_var, fc_var, cog_var, sep = "_")]] <- list(
        lm_tau = lm_tau,
        lm_fc = lm_fc,
        lm_cog = lm_cog
      )
    }
  }
}


# Mediation analyses

# fc1 -> pTau217 -> cognitive_1
cat("\nMediation Analysis 1: functional_connectivity_change_1 -> pTau217 -> cognitive_1\n")
med1_fit <- lm(pTau217 ~ functional_connectivity_change_1 + age + sex, data = data)
out1_fit <- lm(cognitive_1 ~ functional_connectivity_change_1 + pTau217 + age + sex, data = data)
set.seed(123)
med1 <- mediate(med1_fit, out1_fit, treat = "functional_connectivity_change_1",
                mediator = "pTau217", boot = TRUE, sims = 2000)
summary(med1)
plot(med1)

# Example mediation 2: fc2 -> pTau181 -> cognitive_2
cat("\nMediation Analysis 2: functional_connectivity_change_2 -> pTau181 -> cognitive_2\n")
med2_fit <- lm(pTau181 ~ functional_connectivity_change_2 + age + sex, data = data)
out2_fit <- lm(cognitive_2 ~ functional_connectivity_change_2 + pTau181 + age + sex, data = data)
set.seed(123)
med2 <- mediate(med2_fit, out2_fit, treat = "functional_connectivity_change_2",
                mediator = "pTau181", boot = TRUE, sims = 2000)
summary(med2)
plot(med2)


# diagnostics

par(mfrow = c(2,2))
plot(lm_tau)
plot(lm_fc)
plot(lm_cog)
shapiro.test(residuals(lm_tau))
shapiro.test(residuals(lm_fc))
shapiro.test(residuals(lm_cog))
