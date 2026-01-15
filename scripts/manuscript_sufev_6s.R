# ==============================================================================
# Sufficiency Evidence Analysis: 6-Year-Old Children Data
# ==============================================================================
# 
# Description:
#   This script analyzes peeking behavior in consistent and inconsistent reward-size conditions
#   for 6-year-old children using mixed-effects logistic regression models
#   and creates visualizations.
#
# Author: Calder Hilde-Jones
# Last Modified: Jan 15, 2025
# 
# Requirements:
#   - R packages: tidyverse, lme4, PairedData
#   - Data file: sufev_6_data.csv
#   - Function file: boot_glmm.r
#   - Workspace file (optional): model_plot_sufev6_prelim_Nov24_2025.RData
#
# ==============================================================================

# Load Required Libraries ------------------------------------------------------
library(tidyverse)
library(lme4)

# Load Previous Workspace (Optional) -------------------------------------------
load("model_plot_sufev6_prelim_Jan15.RData")


# ==============================================================================
# PART 1: DATA LOADING AND PREPROCESSING
# ==============================================================================

# Load Raw Data ----------------------------------------------------------------
sufev6_data <- read.csv("sufev_6_data.csv")

# Data Preprocessing -----------------------------------------------------------
suf_6 <- sufev6_data %>%
  mutate(
    # Convert Peek to binary (hp = 1, other = 0)
    Peek = if_else(Peek == "hp", 1, 0),
    
    # Convert Choice to binary (n = 0, other = 1)
    #Choice = if_else(Correct == "no", 0, 1),
    
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


# Dataset with only test condition to test whether they were more likely to choose correctly when peeking
suf_6_test <- suf_6 %>%
  mutate(choice = if_else(Correct == "no", 0, 1)) %>%
  filter(Condition.Binary == 1)

# Dataset with only high-peeking in test condition to see whether they ever pick wrong

suf_6_hp <- suf_6_test %>%
  filter(Peek == 1)




# ==============================================================================
# PART 2a: STATISTICAL MODELS - PEEKING BEHAVIOR
# ==============================================================================

# Set Optimizer Control Parameters ---------------------------------------------
contr <- glmerControl(
  optimizer = "bobyqa",
  optCtrl = list(maxfun = 1000000000)
)

# Model 1: Full Model (Condition × Trial Interaction) -------------------------
suf_6_full <- glmer(
  Peek ~ Condition * z.trial + (1 + Condition.test * z.trial | ID),
  data = suf_6,
  control = contr,
  family = binomial(link = "logit")
)

summary(suf_6_full)

# Model 2: Reduced Model (Additive Effects Only) ------------------------------
suf_6_red <- glmer(
  Peek ~ Condition + z.trial + (1 + Condition.test * z.trial | ID),
  data = suf_6,
  control = contr,
  family = binomial(link = "logit")
)

summary(suf_6_red)

# Model 3: Null Model (Random Effects Only) -----------------------------------
suf_6_null <- glmer(
  Peek ~ 1 + (1 + z.trial * Condition.test | ID),
  data = suf_6,
  control = contr,
  family = binomial(link = "logit")
)



# ==============================================================================
# PART 2b: STATISTICAL MODELS - CHOICE
# ==============================================================================


## Just test

suf_choice_test <- glmer(
  choice ~ Peek * z.trial + (1 + z.trial * Condition.test | ID),
  data = suf_6_test,
  control = contr,
  family = binomial(link = "logit")
)

# Reduced Model

suf_choice_test_red <- glmer(
  choice ~ Peek + z.trial + (1 + z.trial * Condition.test | ID),
  data = suf_6_test,
  control = contr,
  family = binomial(link = "logit"))



# ==============================================================================
# PART 3: MODEL SUMMARIES AND COMPARISONS
# ==============================================================================

#Model Summaries
cat("\n=== Full and Reduced Model Summaries ===\n")
summary(suf_6_full)
summary(suf_6_red)

# Test Term Significance -------------------------------------------------------
cat("\n=== Drop1 Reduced Model ===\n")
drop1(suf_6_red, test = 'Chisq')

#Choice Model Summaries
cat("\n=== Choice Test Only ===\n")
summary(suf_choice_test)

cat("\n=== Choice Test Only Red ===\n")
summary(suf_choice_test_red)


# ==============================================================================
# PART 4: VISUALIZATION 1 - PAIRED PLOT
# ==============================================================================

# Prepare Data for Paired Plot -------------------------------------------------
x_summary <- suf_6 %>%
  group_by(ID, Condition.Binary) %>%
  summarise(average_peek = mean(Peek, na.rm = TRUE), .groups = "drop") %>%
  ungroup()

# Add Jitter for Visualization (Reproducible) ----------------------------------
set.seed(012)
x_summary <- x_summary %>%
  mutate(
    jittered_condition = as.numeric(Condition.Binary) + runif(n(), -0.1, 0.1),
    subject = as.factor(ID)
  )

# Create Paired Plot -----------------------------------------------------------
paired_plot <- ggplot(
  x_summary,
  aes(x = jittered_condition, y = average_peek, group = ID)
) +
  # Lines connecting paired observations
  geom_line(aes(color = Condition.Binary)) +
  
  # Points for each participant in each condition
  geom_point(aes(color = Condition.Binary), size = 5) +
  
  # Styling
  scale_x_continuous(
    breaks = c(2, 1),
    labels = levels(suf_6$Condition),
    expand = c(0.25, 0.25)
  ) +
  labs(
    x = "Condition",
    y = "Average Number of Peeks"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.ticks.y = element_line(),
    panel.border = element_rect(fill = NA, color = "black")
  )

print(paired_plot)


# ==============================================================================
# PART 5: BOOTSTRAP CONFIDENCE INTERVALS
# ==============================================================================

# Source Bootstrap Function ----------------------------------------------------
source("scripts/boot_glmm.r")

# Generate Bootstrap Predictions (Condition Only) ------------------------------
boot_plot_suf <- boot.glmm.pred(
  model.res = suf_6_red,
  excl.warnings = TRUE,
  nboots = 1000,
  para = FALSE,
  level = 0.95,
  use = c("Condition")
)

# Generate Bootstrap Predictions (Condition × Trial) ---------------------------
boot_full_suf <- boot.glmm.pred(
  model.res = suf_6_red,
  excl.warnings = TRUE,
  nboots = 1000,
  para = FALSE,
  level = 0.95,
  use = c("Condition", "z.trial")
)


# ==============================================================================
# PART 6: PREPARE DATA FOR VISUALIZATION
# ==============================================================================

# Bootstrap Prediction Summaries by Condition ----------------------------------
ci_predicted_summary <- boot_plot_suf$ci.predicted %>%
  group_by(Condition) %>%
  summarize(
    fitted = mean(fitted),
    lower.cl = mean(lower.cl),
    upper.cl = mean(upper.cl),
    .groups = "drop"
  )

# Individual Participant Means by Condition ------------------------------------
suf_individual_agg <- suf_6 %>%
  group_by(ID, Condition) %>%
  summarize(ind_mean = mean(Peek, na.rm = TRUE), .groups = "drop") %>%
  ungroup()


# ==============================================================================
# PART 7: VISUALIZATION 2 - VIOLIN PLOT WITH INDIVIDUAL DATA POINTS
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
violin_suf <- ggplot(
  suf_6,
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
  
  # Model predictions with confidence intervals
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

violin_suf


# ==============================================================================
# PART 8: VISUALIZATION 3 - INTERACTION PLOT (CONDITION × TRIAL)
# ==============================================================================

# Calculate Trial Break Points for X-Axis --------------------------------------
trial_breaks <- sapply(1:12, function(x) {
  (x - mean(suf_6$Trial)) / sd(suf_6$Trial)
})

# Create Interaction Plot ------------------------------------------------------
interaction_plot_suf_6 <- ggplot() +
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
    data = boot_full_suf$ci.predicted,
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
    data = boot_full_suf$ci.predicted,
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

print(interaction_plot_suf_6)


# ==============================================================================
# SAVE WORKSPACE
# ==============================================================================

save.image("model_plot_sufev6_prelim_Jan15.RData")

# ==============================================================================
# END OF SCRIPT
# ==============================================================================