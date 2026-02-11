# ==============================================================================
# Sufficiency Evidence Analysis: 4-Year-Old Children Data
# ==============================================================================
# 
# Description:
#   This script analyzes peeking behavior in consistent and inconsistent reward-size conditions
#   for 4-year-old children using mixed-effects logistic regression models
#   and creates visualizations.
#
# Author: Calder Hilde-Jones
# Last Modified: Jan 15, 2025
# 
# Requirements:
#   - R packages: tidyverse, lme4
#   - Data file: sufev_4s.csv
#   - Function file: boot_glmm.r
#   - Workspace file (optional): model_plot_sufevkids_Nov25_2025.RData
#
# ==============================================================================

# Load Required Libraries ------------------------------------------------------
library(tidyverse)
library(lme4)

# Load Previous Workspace (Optional) -------------------------------------------
load("images/sufev_4_Feb10.RData")


# ==============================================================================
# PART 1: DATA LOADING AND PREPROCESSING
# ==============================================================================

# Load Raw Data ----------------------------------------------------------------
data_suf_4 <- read.csv("data/sufev_4_data.csv")

# Data Preprocessing -----------------------------------------------------------
suf_4 <- data_suf_4 %>%
  mutate(
    # Convert Peek to binary (hp = 1, other = 0)
    Peek = if_else(Peek == "hp", 1, 0),
    
    # Convert Choice to binary (n = 0, other = 1)
    Choice = if_else(Correct == "no", 0, 1),
    
    # Convert Condition to binary factor (test = 1, control = 0)
    Condition.Binary = factor(
      if_else(Condition == "test", 1, 0),
      levels = c(1, 0)
    ),
    
    # Create scaled trial predictor
    z.trial = as.numeric(scale(Trial)),
    
    # Center condition for random effects
    Condition.test = as.numeric(Condition.Binary) - mean(as.numeric(Condition.Binary)),
    
    # Rename conditions for interpretability
    Condition = factor(
      if_else(Condition == "test", "inconsistent", "consistent"),
      levels = c("consistent", "inconsistent")
    )
  )


suf_4_trial_1 <- subset(suf_4, suf_4$Trial == 1)



# ==============================================================================
# PART 2: STATISTICAL MODELS - PEEKING BEHAVIOR
# ==============================================================================

# Set Optimizer Control Parameters ---------------------------------------------
contr <- glmerControl(
  optimizer = "bobyqa",
  optCtrl = list(maxfun = 100000000)
)

# Model 1: Full Model (Condition × Trial Interaction) -------------------------
suf_4_full <- glmer(
  Peek ~ Condition * z.trial + (1 + z.trial | ID),
  data = suf_4,
  control = contr,
  family = binomial(link = "logit")
)

summary(suf_4_full)

# Model 2: Trial 1 GLM ------------------------------

suf_4_full_glm <- glm(
  Peek ~ Condition,
  data = suf_4_trial_1,
  family = binomial(link = "logit")
)

drop1(suf_4_full_glm, test = "Chisq")

# Model 3: Reduced Model (Additive Effects Only) ------------------------------
suf_4_red <- glmer(
  Peek ~ Condition + z.trial + (1 + z.trial | ID),
  data = suf_4,
  control = contr,
  family = binomial(link = "logit")
)

summary(suf_4_red)

# Model 4: Null Model (Random Effects Only) -----------------------------------
suf_4_null <- glmer(
  Peek ~ 1 + (1 + z.trial | ID),
  data = suf_4,
  control = contr,
  family = binomial(link = "logit")
)


# Extract Variance Components --------------------------------------------------
VarCorr(suf_4_red)


# ==============================================================================
# PART 3: MODEL SUMMARIES AND COMPARISONS
# ==============================================================================
cat("\n=== Model Summaries ===\n")
summary(suf_4_full)
summary(suf_4_red)

cat("\n=== Term Significance (Full Model) ===\n")
drop1(suf_4_full, test = "Chisq")

cat("\n=== Term Significance (Reduced Model) ===\n")
drop1(suf_4_red, test = "Chisq")

# Cross-Tabulation of Peek by Condition ----------------------------------------
table(suf_4$Peek, suf_4$Condition)


# ==============================================================================
# PART 4: BOOTSTRAP CONFIDENCE INTERVALS
# ==============================================================================

# Source Bootstrap Function ----------------------------------------------------
source("scripts/boot_glmm.r")

# Generate Bootstrap Predictions (Condition × Trial) ---------------------------
boot_full_suf_kids <- boot.glmm.pred(
  model.res = suf_4_red,
  excl.warnings = TRUE,
  nboots = 1000,
  para = FALSE,
  level = 0.95,
  use = c("Condition", "z.trial")
)

# Generate Bootstrap Predictions (Condition Only) ------------------------------
boot_plot_suf_kids <- boot.glmm.pred(
  model.res = suf_4_red,
  excl.warnings = TRUE,
  nboots = 1000,
  para = FALSE,
  level = 0.95,
  use = c("Condition")
)


# ==============================================================================
# PART 5: PREPARE DATA FOR VISUALIZATION
# ==============================================================================

# Bootstrap Prediction Summaries by Condition ----------------------------------
ci_predicted_summary <- boot_plot_suf_kids$ci.predicted %>%
  group_by(Condition) %>%
  summarize(
    fitted = mean(fitted),
    lower.cl = mean(lower.cl),
    upper.cl = mean(upper.cl),
    .groups = "drop"
  )

# Individual Participant Means by Condition ------------------------------------
suf_individual_agg <- suf_4 %>%
  group_by(ID, Condition) %>%
  summarize(ind_mean = mean(Peek, na.rm = TRUE), .groups = "drop") %>%
  ungroup()


# ==============================================================================
# PART 6: VISUALIZATION 2 - VIOLIN PLOT WITH INDIVIDUAL DATA POINTS
# ==============================================================================

# Create Jittered Positions for Visualization ----------------------------------
set.seed(123)

# Generate one jitter value per participant
jitter_values_suf <- data.frame(
  ID = unique(suf_individual_agg$ID),
  jitter_x = runif(length(unique(suf_individual_agg$ID)), -0.10, 0.10),
  jitter_y = runif(length(unique(suf_individual_agg$ID)), -0.015, 0.015)
)

# Add jitter to suf_individual_agg
suf_individual_agg_jittered <- suf_individual_agg %>%
  left_join(jitter_values_suf, by = "ID") %>%
  mutate(
    x_num = ifelse(Condition == "inconsistent", 2, 1),  # adjust based on which condition you want on which side
    x_jittered = x_num + jitter_x,
    y_jittered = ind_mean + jitter_y
  )

# Add Numeric Positions to Prediction Summary ----------------------------------
ci_predicted_summary_suf <- ci_predicted_summary %>%
  mutate(
    x_num = ifelse(Condition == "inconsistent", 2, 1)  # match the above
  )

# Create Violin Plot -----------------------------------------------------------
violin_suf_kids <- ggplot(
  suf_4,
  aes(x = Condition, y = Peek, fill = Condition)
) +
  # Violin plot showing distribution
  geom_violin(
    alpha = 0.25,
    color = NA
  ) +
  
  # Lines connecting individual participants
  geom_line(
    data = suf_individual_agg_jittered,
    aes(x = x_jittered, y = y_jittered, group = ID),
    color = "grey70",
    alpha = 0.4,
    linewidth = 0.4
  ) +
  
  # Individual participant means
  geom_point(
    data = suf_individual_agg_jittered,
    aes(x = x_jittered, y = y_jittered, color = Condition),
    size = 1.5,
    alpha = 0.45
  ) +
  
  # Model predictions with confidence intervals (commented out)
  geom_errorbar(
   data = ci_predicted_summary_suf,
  aes(
   x = x_num,
  y = fitted,
     ymin = lower.cl,
     ymax = upper.cl,
     color = Condition
   ),
   width = 0.1,
   linewidth = 0.5
  ) +
  geom_point(
   data = ci_predicted_summary_suf,
   aes(x = x_num, y = fitted),
   size = 2
  ) +
  
  # Styling
  scale_fill_viridis_d(begin = 0.75, end = 0.25, option = "magma") +
  scale_color_viridis_d(begin = 0.75, end = 0.25, option = "magma") +
  scale_x_discrete(labels = str_to_title) +
  scale_y_continuous(
    name = "Rate of High-Peeks",
    breaks = c(0, 0.25, 0.5, 0.75, 1.0),
    labels = scales::percent
  ) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 13),
    strip.text.x = element_text(size = 13),
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    axis.title.y.left = element_text(vjust = 3),
    plot.title = element_text(color = "black", size = 15, face = "bold"),
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
    legend.position = "none"
  )

violin_suf_kids


# ==============================================================================
# PART 7: VISUALIZATION 3 - INTERACTION PLOT (CONDITION × TRIAL)
# ==============================================================================

# Calculate Trial Break Points for X-Axis --------------------------------------
trial_breaks <- sapply(1:12, function(x) {
  (x - mean(suf_4$Trial)) / sd(suf_4$Trial)
})

# Create Interaction Plot ------------------------------------------------------
interaction_plot_suf_4 <- ggplot() +
  # Raw data points (commented out)
  #geom_point(
  #  data = xdata,
  #  aes(x = z.trial, y = Peek, fill = Condition),
  #  size = 2.5,
  #  alpha = 0.8,
  #  shape = 21,
  #  color = "black",
  #  position = position_jitter(width = 0.1, height = 0.05, seed = 81293)
  #) +
  
  # Confidence ribbons
  geom_ribbon(
    data = boot_full_suf_kids$ci.predicted,
    aes(
      ymin = lower.cl,
      ymax = upper.cl,
      x = z.trial,
      group = Condition,
      fill = Condition
    ),
    alpha = 0.2
  ) +
  
  # Fitted lines
  geom_line(
    data = boot_full_suf_kids$ci.predicted,
    aes(x = z.trial, y = fitted, group = Condition, color = Condition),
    linewidth = 1.3
  ) +
  
  # Styling
  scale_y_continuous(
    name = "Rate of High-Peeks",
    breaks = c(0, .25, .50, .75, 1.0),
    labels = scales::percent
  ) +
  scale_x_continuous(
    name = "Trial",
    breaks = trial_breaks,
    labels = 1:12
  ) +
  scale_fill_viridis_d(begin = 0.75, end = 0.25, option = "magma", labels = str_to_title) +
  scale_color_viridis_d(begin = 0.75, end = 0.25, option = "magma", labels = str_to_title) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 13),
    strip.text.x = element_text(size = 13),
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    axis.title.y.left = element_text(vjust = 3),
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
    legend.position = "right"
  )

print(interaction_plot_suf_4)


# ==============================================================================
# SAVE WORKSPACE
# ==============================================================================

save.image("images/sufev_4_Feb10.RData")

# ==============================================================================
# END OF SCRIPT
# ==============================================================================