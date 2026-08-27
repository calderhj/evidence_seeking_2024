# ==============================================================================
# Sufficiency Evidence Analysis: Consistent vs Inconsistent Conditions
# ==============================================================================
# 
# Description:
#   This script analyzes peeking behavior in consistent and inconsistent reward-size conditions 
#   using mixed-effects logistic regression models and creates visualizations.Chimpanzees.
#
# Author: Calder Hilde-Jones
# Last Modified: Jan 15, 2025
# 
# Requirements:
#   - R packages: tidyverse, lme4
#   - Data file: sufev_data.csv
#   - Function file: boot_glmm.r
#   - Workspace file (optional): model_plot_sufev_Oct2_2025.RData
#
# ==============================================================================

# Load Required Libraries ------------------------------------------------------
library(tidyverse)
library(lme4)

# Load Previous Workspace (Optional) -------------------------------------------
load("images/sufev_chimps_Jan15.RData")


# ==============================================================================
# PART 1: DATA LOADING AND PREPROCESSING
# ==============================================================================

# Load Raw Data ----------------------------------------------------------------
mydata_suf <- read.csv("data/sufev_chimp_data.csv")

# Data Preprocessing -----------------------------------------------------------
suf_data <- mydata_suf %>%
  mutate(
    # Convert Peek to binary (hp = 1, other = 0)
    Peek = if_else(Peek == "hp", 1, 0),
    
    # Convert Choice to binary (n = 0, other = 1)
    choice = if_else(Choice == "no", 0, 1),
    
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

suf_data_test <- suf_data %>%
  filter(Condition.Binary == 1)


# ==============================================================================
# PART 2a: STATISTICAL MODELS - PEEKING BEHAVIOR
# ==============================================================================

# Set Optimizer Control Parameters ---------------------------------------------
contr <- glmerControl(
  optimizer = "bobyqa",
  optCtrl = list(maxfun = 50000000)
)

# Model 1: Full Model (Condition × Trial Interaction) -------------------------
suf_model_peek_full <- glmer(
  Peek ~ Condition * z.trial + (1 + z.trial * Condition.test | ID),
  data = suf_data,
  control = contr,
  family = binomial(link = "logit")
)

# Model 2: Reduced Model (Additive Effects Only) ------------------------------
suf_model_peek_red <- glmer(
  Peek ~ Condition + z.trial + (1 + z.trial * Condition.test | ID),
  data = suf_data,
  control = contr,
  family = binomial(link = "logit")
)

# Model 3: Null Model (Random Effects Only) -----------------------------------
suf_model_peek_null <- glmer(
  Peek ~ 1 + (1 + z.trial * Condition.test | ID),
  data = suf_data,
  control = contr,
  family = binomial(link = "logit")
)

# Extract Variance Components --------------------------------------------------
VarCorr(suf_model_peek_red)

# ==============================================================================
# PART 2b: STATISTICAL MODELS - CHOICE
# ==============================================================================

suf_choice_model <- glmer(
  choice ~ Condition * z.trial + (1 + z.trial * Condition.test | ID),
  data = suf_data,
  control = contr,
  family = binomial(link = "logit")
)

suf_choice_model_red <- glmer(
  choice ~ Condition + z.trial + (1 + z.trial * Condition.test | ID),
  data = suf_data,
  control = contr,
  family = binomial(link = "logit")
)


## Just test

suf_choice_test <- glmer(
  choice ~ Peek * z.trial + (1 + z.trial | ID),
  data = suf_data_test,
  control = contr,
  family = binomial(link = "logit")
)

# Reduced Model

suf_choice_test_red <- glmer(
  choice ~ Peek + z.trial + (1 + z.trial | ID),
  data = suf_data_test,
  control = contr,
  family = binomial(link = "logit")
)

# ==============================================================================
# PART 3: MODEL SUMMARIES AND COMPARISONS
# ==============================================================================

cat("\n=== Reduced Model Summary ===\n")
summary(suf_model_peek_red)

cat("\n=== Term Significance (Full Model) ===\n")
drop1(suf_model_peek_full, test = "Chisq")

cat("\n=== Term Significance (Reduced Model) ===\n")
drop1(suf_model_peek_red, test = "Chisq")

cat("\n=== Choice Full Model ===\n")
summary(suf_choice_model)

cat("\n=== Choice Red Model ===\n")
summary(suf_choice_model_red)

cat("\n=== Choice Test Only ===\n")
summary(suf_choice_test)

cat("\n=== Choice Test Only Red ===\n")
summary(suf_choice_test_red)

drop1(suf_choice_test_red, test = "Chisq")




# ==============================================================================
# PART 4: BOOTSTRAP CONFIDENCE INTERVALS
# ==============================================================================

# Source Bootstrap Function ----------------------------------------------------
source("scripts/boot_glmm.r")

# Generate Bootstrap Predictions (Condition × Trial) ---------------------------
boot_full_suf <- boot.glmm.pred(
  model.res = suf_model_peek_red,
  excl.warnings = TRUE,
  nboots = 1000,
  para = FALSE,
  level = 0.95,
  use = c("Condition", "z.trial")
)

# Generate Bootstrap Predictions (Condition Only) ------------------------------
boot_plot_suf <- boot.glmm.pred(
  model.res = suf_model_peek_red,
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
ci_predicted_summary <- boot_plot_suf$ci.predicted %>%
  group_by(Condition) %>%
  summarize(
    fitted = mean(fitted),
    lower.cl = mean(lower.cl),
    upper.cl = mean(upper.cl),
    .groups = "drop"
  )

# Individual Participant Means by Condition ------------------------------------
suf_individual_agg <- suf_data %>%
  group_by(ID, Condition) %>%
  summarize(ind_mean = mean(Peek, na.rm = TRUE), .groups = "drop") %>%
  ungroup()


# ==============================================================================
# PART 6: VISUALIZATION 1 - VIOLIN PLOT WITH INDIVIDUAL DATA POINTS
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
  suf_data,
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
    plot.margin = unit(c(.05, 0, .2, .3), "cm"),
    axis.title.y.left = element_text(vjust = 3),
    plot.title = element_text(color = "black", size = 15, face = "bold"),
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
    legend.position = "none"
  )

print(violin_suf)


# ==============================================================================
# PART 7: VISUALIZATION 2 - INTERACTION PLOT (CONDITION × TRIAL)
# ==============================================================================

# Calculate Trial Break Points for X-Axis --------------------------------------
trial_breaks <- sapply(1:12, function(x) {
  (x - mean(suf_data$Trial)) / sd(suf_data$Trial)
})

# Create Interaction Plot ------------------------------------------------------
interaction_plot_suf_chimps <- ggplot() +
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
    plot.margin = unit(c(.05, .05, .2, .3), "cm"),
    axis.title.y.left = element_text(vjust = 3),
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
    legend.position = "right"
  )

print(interaction_plot_suf_chimps)



# ==============================================================================
# For Revisions
# ==============================================================================

### Do they pick above chance when high-peeking in the size-inconsistent condition

suf_data_test_hp <- suf_data_test %>%
  filter(Peek == 1)

table(suf_data_test_hp$Choice)

suf_choice_test2 <- glmer(
  choice ~ Peek * z.trial + (1 + z.trial * Condition.test | ID),
  data = suf_data_test,
  control = contr,
  family = binomial(link = "logit")
)

summary(suf_choice_test2)

suf_choice_test2_red <- glmer(
  choice ~ Peek + z.trial + (1 + z.trial | ID),
  data = suf_data_test,
  control = contr,
  family = binomial(link = "logit")
)

summary(suf_choice_test2_red)


suf_choice_test2_null <- glmer(
  choice ~ 1 + (1 + z.trial | ID),
  data = suf_data_test,
  control = contr,
  family = binomial(link = "logit")
)

summary(suf_choice_test2_null)

anova(suf_choice_test2_red, suf_choice_test2_null)


###only hp and test condition

suf_hp_test <- suf_data_test %>%
  filter(Peek == 1)

suf_hp_test_choice <- glmer(
  choice ~ 1 + (1 + z.trial | ID),
  data = suf_hp_test,
  control = contr,
  family = binomial(link = "logit")
)

summary(suf_hp_test_choice)

## Is there a difference in the main comparison if we account for repeated subjects across studies?

suf_data$repeated <- ifelse(suf_data$ID %in% c(4, 5, 9, 12), 0, 1)

### Three way interaction
suf_model_repeat_full <- glmer(
  Peek ~ Condition * z.trial * repeated + (1 + z.trial * Condition.test | ID),
  data = suf_data,
  control = contr,
  family = binomial(link = "logit")
)

summary(suf_model_repeat_full)
drop1(suf_model_repeat_full, test = "Chisq")


### All two way interactions

suf_model_repeat_red <- glmer(
  Peek ~ (Condition + repeated + z.trial)^2 + (1 + z.trial * Condition.test | ID),
  data = suf_data,
  control = contr,
  family = binomial(link = "logit")
)


summary(suf_model_repeat_red)
drop1(suf_model_repeat_red, test = "Chisq")

## Full reduced
suf_model_repeat_red2 <- glmer(
  Peek ~ Condition * repeated + z.trial + (1 + z.trial * Condition.test | ID),
  data = suf_data,
  control = contr,
  family = binomial(link = "logit")
)

summary(suf_model_repeat_red2)
drop1(suf_model_repeat_red2, test = "Chisq")

### The model from which we get the effects we report
drop1(suf_model_peek_red, test = "Chisq")

plot(Effect(c("repeated"), suf_model_repeat_red2))

### ---- Effect of familiarization length

suf_data$fam_length <- as.numeric(scale(suf_data$fam_length))

model_fam_full_full <- glmer(
  Peek ~ Condition * fam_length * z.trial + (1 + z.trial * Condition.test | ID),
  data = suf_data,
  control = contr,
  family = binomial(link = "logit")
)

summary(model_fam_full_full)

drop1(model_fam_full_full, test = "Chisq")


model_fam_full_two_way <- glmer(
  Peek ~ (Condition + fam_length + z.trial)^2 + (1 + z.trial * Condition.test | ID),
  data = suf_data,
  control = contr,
  family = binomial(link = "logit")
)

summary(model_fam_full_two_way)

drop1(model_fam_full_two_way, test = "Chisq")


model_fam_full_red <- glmer(
  Peek ~ (Condition * fam_length) + z.trial + (1 + z.trial * Condition.test | ID),
  data = suf_data,
  control = contr,
  family = binomial(link = "logit")
)

summary(model_fam_full_red)

drop1(model_fam_full_red, test = "Chisq")



# ==============================================================================
# SAVE WORKSPACE
# ==============================================================================

save.image("images/sufev_chimps_Jan15.RData")

cat("\n=== Analysis Complete ===\n")
cat("Workspace saved to: model_plot_sufev_Oct2_2025.RData\n")

# ==============================================================================
# END OF SCRIPT
# ==============================================================================
