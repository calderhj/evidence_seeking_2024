# ==============================================================================
# Peeking Behavior Analysis: Ambiguous vs Unambiguous Conditions (Human Data)
# ==============================================================================
# 
# Description:
#   This script analyzes peeking behavior across ambiguous and unambiguous evidence conditions using
#   mixed-effects logistic regression models and creates visualizations for
#   human participant data.
#
# 
# Requirements:
#   - R packages: tidyverse, PairedData, lme4
#   - Data file: evsk_humans.csv
#   - Function file: boot_glmm.r
#   
#
# ==============================================================================

# Load Required Libraries ------------------------------------------------------
library(tidyverse)
library(lme4)

# Load Previous Workspace (Optional) -------------------------------------------
load("images/evsk_4_image.RData")


# ==============================================================================
# PART 1: DATA LOADING AND PREPROCESSING
# ==============================================================================

# Load Raw Data ----------------------------------------------------------------
xdata <- read.csv("data/evsk_4_data.csv")

# Data Preprocessing -----------------------------------------------------------
xdata <- xdata %>%
  # Convert Peek to binary (p = 1, np = 0)
  mutate(
    Peek = case_when(
      Peek == "p" ~ 1,
      Peek == "np" ~ 0,
      TRUE ~ as.numeric(Peek)
    ),
    
    # Convert Condition to binary (ambiguous = 1, unambiguous = 0)
    Condition.Binary = case_when(
      Condition == "ambiguous" ~ 1,
      Condition == "unambiguous" ~ 0,
      TRUE ~ as.numeric(Condition)
    ),
    
    # Ensure numeric types
    Peek = as.numeric(Peek),
    Condition.Binary = as.numeric(Condition.Binary),
    
    # Create scaled trial predictor
    z.trial = as.numeric(scale(Trial)),
    
    # Center condition for random effects
    Condition.test = Condition.Binary - mean(Condition.Binary),
    
    # Set factor levels for plotting
    Condition = factor(Condition, levels = c("unambiguous", "ambiguous"))
  )


# ==============================================================================
# PART 2: STATISTICAL MODELS - PEEKING BEHAVIOR
# ==============================================================================

# Set Optimizer Control Parameters ---------------------------------------------
contr <- glmerControl(
  optimizer = "bobyqa",
  optCtrl = list(maxfun = 50000000)
)

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



# ==============================================================================
# PART 3: MODEL SUMMARIES AND COMPARISONS
# ==============================================================================

cat("\n=== Full Model Summary ===\n")
summary(model_peek_full)

cat("\n=== Reduced Model Summary ===\n")
summary(model_peek_red)

cat("\n=== Term Significance (Full Model) ===\n")
drop1(model_peek_full, test = "Chisq")

cat("\n=== Term Significance (Reduced Model) ===\n")
drop1(model_peek_red, test = "Chisq")


# ==============================================================================
# PART 4: BOOTSTRAP CONFIDENCE INTERVALS
# ==============================================================================

# Source Bootstrap Function ----------------------------------------------------
source("scripts/boot_glmm_submission.r")

# Generate Bootstrap Predictions (Condition × Trial) ---------------------------
boot_full <- boot.glmm.pred(
  model.res = model_peek_full,
  excl.warnings = TRUE,
  nboots = 1000,
  para = FALSE,
  level = 0.95,
  use = c("Condition", "z.trial")
)

# Generate Bootstrap Predictions (Condition Only) ------------------------------
boot_plot <- boot.glmm.pred(
  model.res = model_peek_red,
  excl.warnings = TRUE,
  nboots = 1000,
  para = FALSE,
  level = 0.95,
  use = c("Condition")
)


# ==============================================================================
# PART 5: PREPARE DATA FOR VISUALIZATION
# ==============================================================================

# Overall Means by Condition ---------------------------------------------------
violin_agg <- xdata %>%
  group_by(Condition) %>%
  summarise(peek.mean = mean(Peek, na.rm = TRUE)) %>%
  ungroup()

# Individual Participant Means -------------------------------------------------
violin_individual_agg <- xdata %>%
  group_by(ID, Condition) %>%
  summarize(ind_peek_mean = mean(Peek, na.rm = TRUE), .groups = "drop") %>%
  ungroup()

# Bootstrap Prediction Summaries -----------------------------------------------
ci_predicted_summary <- boot_full$ci.predicted %>%
  group_by(Condition) %>%
  summarize(
    fitted = mean(fitted),
    lower.cl = mean(lower.cl),
    upper.cl = mean(upper.cl),
    .groups = "drop"
  )

# Ensure Consistent Factor Levels for Plotting ---------------------------------
boot_full$ci.predicted$Condition <- factor(
  boot_full$ci.predicted$Condition,
  levels = c("unambiguous", "ambiguous")
)


# ==============================================================================
# PART 6: VISUALIZATION 1 - VIOLIN PLOT WITH INDIVIDUAL DATA POINTS
# ==============================================================================

# Create Jittered Positions for Visualization ----------------------------------
set.seed(123)

# Generate one jitter value per participant
jitter_values_hum <- data.frame(
  ID = unique(violin_individual_agg$ID),
  jitter_x = runif(length(unique(violin_individual_agg$ID)), -0.03, 0.03),
  jitter_y = runif(length(unique(violin_individual_agg$ID)), -0.030, 0.030)
)

# Add jitter to violin_individual_agg
violin_individual_agg_jittered <- violin_individual_agg %>%
  left_join(jitter_values_hum, by = "ID") %>%
  mutate(
    x_num = ifelse(Condition == "ambiguous", 2, 1),
    x_jittered = x_num + jitter_x,
    y_jittered = ind_peek_mean + jitter_y
  )

# Add Numeric Positions to Prediction Summary ----------------------------------
ci_predicted_summary_positioned <- ci_predicted_summary %>%
  mutate(
    x_num = ifelse(Condition == "ambiguous", 2, 1)
  )

# Create Violin Plot -----------------------------------------------------------
violin_hum_evsk <- ggplot(
  violin_individual_agg_jittered,
  aes(x = Condition, y = ind_peek_mean, fill = Condition)
) +
  # Violin plot showing distribution
  geom_violin(
    data = xdata,
    aes(x = Condition, y = Peek, fill = Condition),
    alpha = 0.2,
    color = NA
  ) +
  
  # Lines connecting individual participants
  geom_line(
    aes(x = x_jittered, y = y_jittered, group = ID),
    color = "grey70",
    alpha = 0.4,
    linewidth = 0.4
  ) +
  
  # Individual participant means
  geom_point(
    aes(x = x_jittered, y = y_jittered, color = Condition),
    size = 1.5,
    alpha = 0.4
  ) +
  
  # Model predictions with confidence intervals
  geom_errorbar(
    data = ci_predicted_summary_positioned,
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
    data = ci_predicted_summary_positioned,
    aes(x = x_num, y = fitted),
    size = 2
  ) +
  
  # Styling
  scale_fill_viridis_d(begin = 0.75, end = 0.25, option = "magma") +
  scale_color_viridis_d(begin = 0.75, end = 0.25, option = "magma") +
  scale_x_discrete(labels = str_to_title) +
  scale_y_continuous(
    name = "Rate of Peeking",
    breaks = c(0, 0.25, 0.5, 0.75, 1.0),
    labels = scales::percent
  ) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 13),
    strip.text.x = element_text(size = 13),
    plot.margin = unit(c(0.05, 0.05, 0.22, 0.30), "cm"),
    axis.title.y.left = element_text(vjust = 3),
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
    legend.position = "none"
  )

violin_hum_evsk


# ==============================================================================
# PART 7: VISUALIZATION 2 - INTERACTION PLOT (CONDITION × TRIAL)
# ==============================================================================

# Calculate Trial Break Points for X-Axis --------------------------------------
trial_breaks <- sapply(1:12, function(x) {
  (x - mean(xdata$Trial)) / sd(xdata$Trial)
})

# Create Interaction Plot ------------------------------------------

interaction_plot_hum <- ggplot() +
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
    linewidth = 1.3
  ) +
  
  # Styling
  scale_y_continuous(
    name = "Rate of Peeking",
    breaks = c(0, .25, .50, .75, 1.0),
    labels = scales::percent
  ) +
  scale_x_continuous(
    name = "Trial",
    breaks = trial_breaks,
    labels = 1:12
  ) +
  scale_fill_viridis_d(
    begin = 0.75, 
    end = 0.25, 
    option = "magma", 
    labels = str_to_title,
    breaks = c("ambiguous", "unambiguous")
  ) +
  scale_color_viridis_d(
    begin = 0.75, 
    end = 0.25, 
    option = "magma", 
    labels = str_to_title,
    breaks = c("ambiguous", "unambiguous")
  ) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 13),
    strip.text.x = element_text(size = 13),
    plot.margin = unit(c(.22, .22, .22, .30), "cm"),
    axis.title.y.left = element_text(vjust = 3),
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
    legend.position = "right"
  )

print(interaction_plot_hum)


# ==============================================================================
# PART 8: BOX CHOICE ANALYSES
# ==============================================================================

# Prepare Choice Data ----------------------------------------------------------
choice_dat <- xdata %>%
  mutate(Choice = ifelse(Choice == "y", 1, ifelse(Choice == "n", 0, Choice))) %>%
  mutate(Choice = ifelse(Choice == "yes", 1, ifelse(Choice == "no", 0, Choice)))

choice_dat$Choice <- as.numeric(choice_dat$Choice)

# Model 1: Full Choice Model (Condition × Trial Interaction) ------------------
model_choice_full <- glmer(Choice ~ Condition * z.trial + (1 + Condition.test * z.trial| ID), 
                           data = choice_dat,
                           control = contr,
                           family = binomial(link = "logit"))

summary(model_choice_full)

# Model 2: Reduced Choice Model ------------------------------------------------
model_choice_red <- glmer(Choice ~ Condition + z.trial + (1 + Condition.test * z.trial | ID), 
                          data = choice_dat,
                          control = contr,
                          family = binomial(link = "logit"))

summary(model_choice_red)

# Model 3: Null Choice Model ---------------------------------------------------
model_choice_null <- glmer(Choice ~ 1 + (1 + Condition.test | ID), 
                           data = choice_dat,
                           control = contr,
                           family = binomial(link = "logit"))

# Test Term Significance -------------------------------------------------------
drop1(model_choice_red, test = 'Chisq')
drop1(model_choice_full, test = 'Chisq')
# Filter Data by Condition -----------------------------------------------------
choice_dat_ambiguous <- choice_dat %>%
  filter(Condition == 'ambiguous') %>%
  filter(Peek == 1)

choice_dat_unambiguous <- choice_dat %>%
  filter(Condition == 'unambiguous')


# Tabulate Choices in Unambiguous Condition ------------------------------------
table(choice_dat_unambiguous$Choice)


# ==============================================================================
# PART 9: BOOTSTRAP FOR CHOICE ANALYSIS
# ==============================================================================

# Generate Bootstrap Predictions for Choice Model -----------------------------
boot_choice_full <- boot.glmm.pred(
  model.res = model_choice_full,
  excl.warnings = TRUE,
  nboots = 1000,
  para = FALSE,
  level = 0.95,
  use = c("Condition", "z.trial")
)


# ==============================================================================
# PART 10: VISUALIZATION 3 - CHOICE INTERACTION PLOT
# ==============================================================================

# Calculate Trial Break Points (reuse from earlier) ----------------------------
trial_breaks <- sapply(1:12, function(x) {
  (x - mean(xdata$Trial)) / sd(xdata$Trial)
})

# Create Choice Interaction Plot -----------------------------------------------
interaction_plot_choice <- ggplot() +
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
    data = boot_choice_full$ci.predicted,
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
    data = boot_choice_full$ci.predicted,
    aes(x = z.trial, y = fitted, group = Condition, color = Condition),
    size = 1.3
  ) +
  
  # Styling
  scale_y_continuous(
    name = "correct choice rate",
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

interaction_plot_choice


# ==============================================================================
# SAVE WORKSPACE
# ==============================================================================

save.image("images/evsk_4_Jan15.RData")



# ==============================================================================
# END OF SCRIPT
# ==============================================================================