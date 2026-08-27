# ==============================================================================
# Evidence Seeking Study: Information Seeking Behavior Analysis
# ==============================================================================
# 
# Description:
#   This script analyzes peeking behavior in ambiguous vs unambiguous conditions
#   using mixed-effects models and visualizations.Chimpanzees.
#
# Requirements:
#   - R packages: tidyverse, lme4, lmerTest
#   - Data file: evsk_data.csv
#   - Function file: boot_glmm.r
# 
#
# ==============================================================================

# Load Required Libraries ------------------------------------------------------
library(tidyverse)
library(lme4)
library(lmerTest)

# R Settings -------------------------------------------------------------------
options(scipen = 999)

# Load Previous Workspace (Optional) -------------------------------------------
load("images/evsk_chimps_image.RData")


# ==============================================================================
# PART 1: DATA LOADING AND PREPROCESSING
# ==============================================================================

# Load Raw Data ----------------------------------------------------------------
mydata <- read.csv("data/evsk_chimp_data.csv")

# Data Preprocessing -----------------------------------------------------------
xdata <- mydata %>%
  mutate(
    # Convert Peek to binary (p = 1, other = 0)
    Peek = if_else(Peek == "p", 1, 0),
    
    # Convert Choice to binary (n = 0, other = 1)
    choice = if_else(Choice == "n", 0, 1),
    
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
      if_else(Condition == "test", "Ambiguous", "Unambiguous"),
      levels = c("Unambiguous", "Ambiguous")
    ),
    
    # Create factor version of ID
    ID.factor = as.factor(ID)
  )


# ==============================================================================
# PART 2: STATISTICAL MODELS - PEEKING BEHAVIOR
# ==============================================================================

# Set Optimizer Control Parameters ---------------------------------------------
contr <- glmerControl(
  optimizer = "bobyqa",
  optCtrl = list(maxfun = 50000000)
)
# Alternative optimizer (if needed):

# Model 1: Full Model (Condition × Trial Interaction) -------------------------
model_peek_full <- glmer(
  Peek ~ Condition * z.trial + (1 + z.trial * Condition.test | ID),
  data = xdata,
  control = contr,
  family = binomial(link = "logit")
)

# Model 2: Reduced Model (Additive Effects Only) ------------------------------
model_peek_red <- glmer(
  Peek ~ Condition + z.trial + (1 + z.trial * Condition.test | ID),
  data = xdata,
  control = contr,
  family = binomial(link = "logit")
)

# Model 3: Null Model (Random Effects Only) -----------------------------------
model_peek_null <- glmer(
  Peek ~ 1 + (1 + z.trial * Condition.test | ID),
  data = xdata,
  control = contr,
  family = binomial(link = "logit")
)

# Model Summaries and Comparisons ----------------------------------------------
cat("\n=== Full Model Summary ===\n")
summary(model_peek_full)

cat("\n=== Term Significance (Full Model) ===\n")
drop1(model_peek_full, test = "Chisq")
drop1(model_peek_red, test = "Chisq")


# ==============================================================================
# PART 3: BOOTSTRAP CONFIDENCE INTERVALS
# ==============================================================================

# Source Bootstrap Function ----------------------------------------------------
source("scripts/boot_glmm_submission.r")

# Generate Bootstrap Predictions -----------------------------------------------
boot_full <- boot.glmm.pred(
  model.res = model_peek_full,
  excl.warnings = TRUE,
  nboots = 1000,
  para = FALSE,
  level = 0.95,
  use = c("Condition", "z.trial")
)

cat("\n=== Bootstrap Predictions ===\n")
print(head(boot_full$ci.predicted))


# ==============================================================================
# PART 4: VISUALIZATION 2 - INTERACTION PLOT (CONDITION × TRIAL)
# ==============================================================================

# Calculate Trial Break Points for X-Axis --------------------------------------
trial_breaks <- sapply(1:12, function(x) {
  (x - mean(xdata$Trial)) / sd(xdata$Trial)
})

# Create Interaction Plot ------------------------------------------------------
interaction_plot <- ggplot() +
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
    data = boot_full$ci.predicted,
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
    data = boot_full$ci.predicted,
    aes(x = z.trial, y = fitted, group = Condition, color = Condition),
    linewidth = 1.7
  ) +
  
  # Styling
  scale_y_continuous(
    name = "Rate of Peeking",
    breaks = c(0, 0.25, 0.5, 0.75, 1.0),
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
    plot.margin = unit(c(-.4, .22, .22, .3), "cm"),
    axis.title.y.left = element_text(vjust = 3),
    plot.title = element_text(color = "black", size = 15, face = "bold"),
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
    legend.position = "right"
  )

print(interaction_plot)


# ==============================================================================
# PART 5: VISUALIZATION 3 - VIOLIN PAIRED DATA PLOT
# ==============================================================================

# Prepare Summary Statistics ---------------------------------------------------
ci_predicted_summary <- boot_full$ci.predicted %>%
  group_by(Condition) %>%
  summarize(
    fitted = mean(fitted),
    lower.cl = mean(lower.cl),
    upper.cl = mean(upper.cl),
    .groups = "drop"
  )

# Individual Participant Means by Condition ------------------------------------
individual_agg <- xdata %>%
  group_by(ID, Condition) %>%
  summarize(ind_mean = mean(Peek, na.rm = TRUE), .groups = "drop") %>%
  ungroup()

# Create Jittered Positions Manually -------------------------------------------
set.seed(123)

# Generate one jitter value per participant
jitter_values <- data.frame(
  ID = unique(individual_agg$ID),
  jitter = runif(length(unique(individual_agg$ID)), -0.05, 0.05)
)


# Add jitter to individual_agg
individual_agg_jittered <- individual_agg %>%
  left_join(jitter_values, by = "ID") %>%
  mutate(
    x_num = ifelse(Condition == "Unambiguous", 1, 2),  # reversed to match limits = rev
    x_jittered = x_num + jitter
  )

# Add x_num to ci_predicted_summary to match
ci_predicted_summary <- ci_predicted_summary %>%
  mutate(x_num = ifelse(Condition == "Unambiguous", 1, 2))


# Create Violin Plot -----------------------------------------------------------
violin_evsk <- ggplot(
  xdata,
  aes(x = Condition, y = Peek, fill = Condition)
) +
  geom_violin(
    alpha = 0.25,
    color = NA
  ) +
  geom_line(
    data = individual_agg_jittered,
    aes(x = x_jittered, y = ind_mean, group = ID),
    color = "grey70",
    alpha = 0.4,
    linewidth = 0.4
  ) +
  geom_point(
    data = individual_agg_jittered,
    aes(x = x_jittered, y = ind_mean, color = Condition),
    size = 1.5,
    alpha = 0.45
  ) +
  geom_errorbar(
    data = ci_predicted_summary,
    aes(x = Condition, y = fitted, ymin = lower.cl, ymax = upper.cl, color = Condition),
    width = 0.1,
    linewidth = 0.5
  ) +
  geom_point(
    data = ci_predicted_summary,
    aes(x = Condition, y = fitted),
    size = 2
  ) +
  scale_fill_viridis_d(begin = 0.75, end = 0.25, option = "magma") +
  scale_color_viridis_d(begin = 0.75, end = 0.25, option = "magma") +
  scale_x_discrete(labels = str_to_title) +
  scale_y_continuous(
    name = "Peeking Rate",
    breaks = c(0, 0.25, 0.5, 0.75, 1.0),
    labels = scales::percent
  ) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 13),
    strip.text.x = element_text(size = 13),
    plot.margin = unit(c(.22, .22, .22, .3), "cm"),
    axis.title.y.left = element_text(vjust = 3),
    plot.title = element_text(color = "black", size = 15, face = "bold"),
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
    legend.position = "none"
  )

violin_evsk


# ==============================================================================
# PART 6: SIDE BIAS ANALYSES
# ==============================================================================

# Prepare Data for Side Bias Analysis -----------------------------------------
xdata_bias <- xdata %>%
  mutate(
    choice_side = if_else(choice == 1, Location, if_else(Location == "R", "L", "R")),
    choice_binary = if_else(choice_side == "R", 1, 0),
  )

# Calculate Bias Summary -------------------------------------------------------
bias_summary <- xdata_bias %>%
  group_by(ID, Name, Condition, Condition.test) %>%
  summarise(average_side = mean(choice_binary, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(Bias = abs(average_side - .50))


# Linear Mixed Model for Bias --------------------------------------------------
model_bias <- lmer(Bias ~ Condition.test + (1 | ID),
                   data = bias_summary)
summary(model_bias)

drop1(model_bias, test = 'Chisq')

# Bootstrap Predictions for Bias -----------------------------------------------

boot_full_bias <-
  boot.glmm.pred(
    model.res = model_bias, excl.warnings = TRUE,
    nboots = 1000, para = FALSE, level = 0.95,
    use = c("Condition")
  )

bias_pred_summary <- boot_full_bias$ci.predicted %>%
  group_by(Condition.test) %>%
  summarize(
    fitted = mean(fitted),
    lower.cl = mean(lower.cl),
    upper.cl = mean(upper.cl)
  )


# ==============================================================================
# PART 7: VISUALIZATION 4 - BIAS VIOLIN PLOT
# ==============================================================================

bias_summary$Condition.test <- as.factor(bias_summary$Condition.test)
bias_pred_summary$Condition.test <- as.factor(bias_pred_summary$Condition.test)

# Quick visualization for side biasing. -.5 = test condition
violin_bias <- ggplot(bias_summary, aes(x = Condition.test, y = Bias, fill = Condition.test)) + 
  geom_violin(alpha = .25, color = NA) +
  geom_point(
    aes(color = Condition.test), 
    size = 1.5, alpha = .45, 
    position = position_jitter(w = 0.015, h = 0.0)
  ) +
  geom_line(
    aes(group = ID),  # Replace 'participant_id' with your actual column name
    alpha = 0.3, 
    color = "gray50",
    linewidth = 0.3
  ) +
  # geom_errorbar(
  #   data = bias_pred_summary,
  #   aes(x = Condition, y = fitted, ymin = lower.cl, ymax = upper.cl, color = Condition),
  #   width = 0.1, linewidth = 0.5
  # ) +
  # geom_point(
  #   data = bias_pred_summary,
  #   aes(x = Condition, y = fitted),
  #   size = 2
  # ) +
  scale_fill_viridis_d(begin = .25, end = .79, option = "magma") +
  scale_color_viridis_d(begin = .25, end = .79, option = "magma") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 13),
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    axis.title.y = element_text(vjust = 3),
    axis.title.x = element_text(margin = margin(t = 10)),
    legend.position = "none"
  )

violin_bias

# ==============================================================================
# PART 8: BOX CHOICE ANALYSES
# ==============================================================================

# Prepare Choice Data ----------------------------------------------------------
choice_dat <- xdata %>%
  mutate(Choice = ifelse(Choice == "y", 1, ifelse(Choice == "n", 0, Choice))) %>%
  mutate(Choice = ifelse(Choice == "yes", 1, ifelse(Choice == "no", 0, Choice)))

choice_dat$Choice <- as.numeric(choice_dat$Choice)

# Model 1: Full Choice Model ---------------------------------------------------
model_choice_full <- glmer(Choice ~ Condition * z.trial + (1 + Condition.test * z.trial | ID), 
                           data = choice_dat,
                           control = contr,
                           family = binomial(link = "logit"))

summary(model_choice_full)

# Model 2: Reduced Choice Model ------------------------------------------------
model_choice_red <- glmer(Choice ~ Condition + z.trial + (1 + Condition.test * z.trial | ID), 
                          data = choice_dat,
                          control = contr,
                          family = binomial(link = "logit"))

# Model 3: Null Choice Model ---------------------------------------------------
model_choice_null <- glmer(Choice ~ 1 + (1 + Condition.test * z.trial | ID), 
                           data = choice_dat,
                           control = contr,
                           family = binomial(link = "logit"))

# Model Summaries --------------------------------------------------------------
drop1(model_choice_full, test = "Chisq")
summary(model_choice_red)
drop1(model_choice_red, test = "Chisq")


boot_full_choice <- boot.glmm.pred(
  model.res = model_choice_full,
  excl.warnings = TRUE,
  nboots = 1000,
  para = FALSE,
  level = 0.95,
  use = c("Condition", "z.trial")
)



### Only test 

test_choice <- choice_dat %>%
  filter(Condition == "Ambiguous")

model_choice_test_full <- glmer(Choice ~ Peek * z.trial + (1 + z.trial | ID), 
                           data = test_choice,
                           control = contr,
                           family = binomial(link = "logit"))

summary(model_choice_test_full)

drop1(model_choice_test_full, test = "Chisq")

model_choice_test <- glmer(Choice ~ Peek + z.trial + (1 + z.trial | ID), 
                           data = test_choice,
                           control = contr,
                           family = binomial(link = "logit"))

drop1(model_choice_test, test = "Chisq")

boot_full_choice_test <- boot.glmm.pred(
  model.res = model_choice_test_full,
  excl.warnings = TRUE,
  nboots = 1000,
  para = FALSE,
  level = 0.95,
  use = c("Condition", "z.trial")
)

# ==============================================================================
# For revisions: Familiarization length analyses (In supplements)
# ==============================================================================

### Add familiarization lengths to data

lookup <- c("1"=9, "2"=10, "3"=13, "4"=15, "5"=10, 
            "6"=12, "7"=8, "8"=14, "9"=12, "10"=16)

xdata$familiarization <- as.numeric(scale(lookup[as.character(xdata$ID)]))

### The full full model

model_fam_full_full <- glmer(
  Peek ~ Condition * familiarization * z.trial + (1 + z.trial * Condition.test | ID),
  data = xdata,
  control = contr,
  family = binomial(link = "logit")
)

summary(model_fam_full_full)

drop1(model_fam_full_full, test = "Chisq")


model_fam_full_two_way <- glmer(
  Peek ~ (Condition + familiarization + z.trial)^2 + (1 + z.trial * Condition.test | ID),
  data = xdata,
  control = contr,
  family = binomial(link = "logit")
)

summary(model_fam_full_two_way)

drop1(model_fam_full_two_way, test = "Chisq")


model_fam_full_red <- glmer(
  Peek ~ (Condition * z.trial) + familiarization + (1 + z.trial * Condition.test | ID),
  data = xdata,
  control = contr,
  family = binomial(link = "logit")
)

summary(model_fam_full_red)

drop1(model_fam_full_red, test = "Chisq")

# ==============================================================================
# SAVE WORKSPACE
# ==============================================================================

save.image("images/evsk_chimps_Jan15.RData")

cat("\n=== Analysis Complete ===\n")
cat("Workspace saved to: model_plot_evsk_Oct2_2025.RData\n")

# ==============================================================================
# END OF SCRIPT
# ==============================================================================