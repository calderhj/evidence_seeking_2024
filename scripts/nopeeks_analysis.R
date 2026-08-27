library(ggplot2)
library(tidyr)
library(dplyr)


suf_nopeek_raw <- read.csv("data/suf_chimps_nopeeks.csv")

### Make dataset with nps and lps represented as 0. nps counted as lps
suf_np_with_lp <- suf_nopeek_raw %>%
  mutate(
    # Convert Peek to binary (hp = 1, other = 0)
    Peek = if_else(Peek == "hp", 1, 0),
    
    
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


### make dataset where nps and hps are represented by 1. nps counted as hps
suf_np_with_hp <- suf_nopeek_raw %>%
  mutate(
    # Convert Peek to binary (hp = 1, other = 0)
    Peek = if_else(Peek == "lp", 1, 0),
    
    
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

suf_lp_with_hp <- suf_nopeek_raw %>%
  mutate(
    # Convert Peek to binary (hp = 1, other = 0)
    Peek = if_else(Peek == "np", 1, 0),
    
    
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



### Analysis 1 ------ hps vs lps and nps
contr <- glmerControl(
  optimizer = "bobyqa",
  optCtrl = list(maxfun = 50000000)
)
# Full Model

lp_np_full <- glmer(
  Peek ~ Condition * z.trial + (1 + z.trial * Condition.test | ID),
  data = suf_np_with_lp,
  control = contr,
  family = binomial(link = "logit")
)

summary(lp_np_full)

# Reduced Model (Additive Effects Only) ------------------------------
lp_np_red <- glmer(
  Peek ~ Condition + z.trial + (1 + z.trial * Condition.test | ID),
  data = suf_np_with_lp,
  control = contr,
  family = binomial(link = "logit")
)

summary(lp_np_red)

drop1(lp_np_red, test = "Chisq")


### Analysis 2 ------ hps and nps vs lps


hp_np_full <- glmer(
  Peek ~ Condition * z.trial + (1 + z.trial * Condition.test | ID),
  data = suf_np_with_hp,
  control = contr,
  family = binomial(link = "logit")
)

summary(hp_np_full)


hp_np_red <- glmer(
  Peek ~ Condition + z.trial + (1 + z.trial * Condition.test | ID),
  data = suf_np_with_hp,
  control = contr,
  family = binomial(link = "logit")
)

summary(hp_np_red)

table(suf_np_with_hp$Peek, suf_np_with_hp$Condition)


### Analysis 3 ---- np vs lp and hp

np_full <- glmer(
  Peek ~ Condition * z.trial + (1 + z.trial * Condition.test | ID),
  data = suf_lp_with_hp,
  control = contr,
  family = binomial(link = "logit")
)

summary(hp_np_full)


np_red <- glmer(
  Peek ~ Condition + z.trial + (1 + z.trial * Condition.test | ID),
  data = suf_np_with_hp,
  control = contr,
  family = binomial(link = "logit")
)

summary(np_red)



##### Stacked bar plot

data_long <- suf_nopeek_raw %>%
  group_by(Condition, Peek) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Condition) %>%
  mutate(percentage = n / sum(n) * 100,
         Peek = factor(Peek,
                       levels = c("np", "lp", "hp"),
                       labels = c("No-peek", "Low-peek", "High-peek")),
         Condition = factor(Condition,
                            levels = c("test", "control"),
                            labels = c("Inconsistent", "Consistent")))


ggplot(data_long, aes(x = Condition, y = percentage, fill = Peek)) +
  geom_bar(stat = "identity", width = 0.5) +
  scale_fill_manual(values = c("High-peek" = "#534AB7",
                               "Low-peek"  = "#1D9E75",
                               "No-peek"   = "#888780")) +
  scale_y_continuous(labels = scales::percent_format(scale = 1),
                     expand = c(0, 0), limits = c(0, 101)) +
  labs(x = NULL, y = "Percentage of trials", fill = NULL) +
  theme_classic() +
  theme(legend.position = "top",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))

#### Analysis 4 ------ Mixed Multinomial

library(mclogit)


suf_nopeek <- suf_nopeek_raw %>%
  mutate(z.trial = scale(Trial),
         Peek = factor(Peek, levels = c("lp", "hp", "np")),
         Condition = factor(Condition, levels = c("control", "test")))


# Model 1: Full model
multinom_np_full <- mblogit(
  Peek ~ Condition * z.trial,
  random = ~ 1 | ID,
  data = suf_nopeek
)



# Model 2: Reduced model
multinom_np_red <- mblogit(
  Peek ~ Condition + z.trial,
  random = ~ 1 + z.trial + Condition | ID,
  data = suf_nopeek
)


summary(multinom_np_red)






#### Analysis 5 ------- Wilcoxon for consistent condition np vs lp


np_no_repeat_raw <- read.csv("data/sufev_np_norepeat.csv")

# Summarise per subject in consistent condition only
consistent_summary <- np_no_repeat_raw %>%
  filter(Condition == "control") %>%
  group_by(ID) %>%
  summarise(
    lp_rate = mean(Peek == "lp"),
    np_rate = mean(Peek == "np")
  )

# Wilcoxon signed-rank test: is lp rate > np rate within subjects?
wilcox.test(consistent_summary$lp_rate, 
            consistent_summary$np_rate, 
            paired = TRUE, 
            alternative = "greater",
            exact = FALSE)



