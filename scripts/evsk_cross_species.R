
### Load data from both populations
xdata <- read.csv("data/cross_species_evsk.csv")

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

### Models and analysis

contr <- glmerControl(
  optimizer = "bobyqa",
  optCtrl = list(maxfun = 50000000)
)

# Model 1: Full Model (Condition × Trial Interaction) -------------------------
model_comp_full <- glmer(
  Peek ~ Condition * z.trial * species + (1 + z.trial * Condition.test | ID),
  data = xdata,
  control = contr,
  family = binomial(link = "logit")
)

summary(model_comp_full)

# Model 2: Reduced Model (two-way) ------------------------------
model_comp_red <- glmer(
  Peek ~ (Condition + z.trial + species)^2 + (1 + z.trial * Condition.test | ID),
  data = xdata,
  control = contr,
  family = binomial(link = "logit")
)

summary(model_comp_red)
drop1(model_comp_red, test = "Chisq")

# Model 3: Reduced Model (interaction I expected) ------------------------------

model_comp_red2 <- glmer(
  Peek ~ Condition + species + z.trial + (1 + z.trial * Condition.test | ID),
  data = xdata,
  control = contr,
  family = binomial(link = "logit")
)

summary(model_comp_red2)
drop1(model_comp_red2, test = "Chisq")
plot(Effect(c("Condition", "species"), model_comp_red2))

### test just first 4 chimp trials

xdata_chimp <- subset(xdata, species == "chimp")

model_chimp_first4 <- glmer(
  Peek ~ Condition + z.trial + (1 + z.trial * Condition.test | ID),
  data = xdata_chimp,
  control = contr,
  family = binomial(link = "logit")
)

summary(model_chimp_first4)


plot(Effect(c("species", "Condition"), model_comp_red2))
