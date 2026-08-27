# ==============================================================================
# Inter-Rater Reliability Analysis: Cohen's Kappa for All Studies
# ==============================================================================
# 
# Description:
#   This script calculates inter-rater reliability using Cohen's Kappa for
#   peeking and choice coding across all studies (chimpanzees and children).
#
# Requirements:
#   - R package: irr
#   - Data files: Multiple CSV files for each study/age group
#
# ==============================================================================

# Load Required Libraries ------------------------------------------------------
library(irr)


# ==============================================================================
# PART 1: EVIDENCE SEEKING - CHIMPANZEES
# ==============================================================================


# Load  Reliability Data ---------------------------------------------
rel_evsk <- read.csv("data/reliability_evsk_chimps.csv")

# Calculate Kappa for Peek Coding (Alternative Dataset) -----------------------
evsk_chimp_peek_kappa <- kappa2(
  rel_evsk[, c("Peek", "Peek.OG.")], 
  weight = "unweighted"
)
evsk_chimp_peek_kappa

# Calculate Kappa for Choice Coding --------------------------------------------
evsk_chimp_choice_kappa <- kappa2(
  rel_evsk[, c("Choice", "Choice.OG.")], 
  weight = "unweighted"
)
evsk_chimp_choice_kappa


# ==============================================================================
# PART 2: SUFFICIENCY EVIDENCE - CHIMPANZEES
# ==============================================================================

# Load Reliability Data (Two Coders) -------------------------------------------
rel_sufev_1 <- read.csv("data/reliability_sufev_chimps1.csv")
rel_sufev_2 <- read.csv("data/reliability_sufev_chimps2.csv")

# Prepare Data Frame for Analysis ----------------------------------------------
sufev_reliability_peek <- data.frame(
  peek1 = rel_sufev_1$Peek, 
  peek2 = rel_sufev_2$Peek
)

# Calculate Kappa for Peek Coding ----------------------------------------------
sufev_rel_peek <- kappa2(sufev_reliability_peek[, c("peek1", "peek2")], weight = "unweighted")

# Load Reliability Data for Choice Coding --------------------------------------
rel_suf <- read.csv("data/reliability_sufev_chimps_choicedata.csv")

# Calculate Kappa for Choice Coding --------------------------------------------
sufev_rel_choice <- kappa2(
  rel_suf[, c("Choice.lil.", "Choice.og.")], 
  weight = "unweighted"
)
sufev_rel_choice


# ==============================================================================
# PART 3: EVIDENCE SEEKING - CHILDREN (HUMANS)
# ==============================================================================

# Load Reliability Data --------------------------------------------------------
rel_evsk_kids <- read.csv("data/reliability_evsk_4.csv")

# Calculate Kappa for Peek Coding ----------------------------------------------
evsk_kid_kappa_peek <- kappa2(
  rel_evsk_kids[, c("Peek_rel", "Peek_og")], 
  weight = "unweighted"
)
evsk_kid_kappa_peek

# Calculate Kappa for Choice Coding --------------------------------------------
evsk_kid_kappa_choice <- kappa2(
  rel_evsk_kids[, c("Choice", "Choice_og")], 
  weight = "unweighted"
)
evsk_kid_kappa_choice


# ==============================================================================
# PART 4: SUFFICIENCY EVIDENCE - 4-YEAR-OLD CHILDREN
# ==============================================================================

# Load Reliability Data --------------------------------------------------------
rel_sufev_4 <- read.csv("data/reliability_evsk_4.csv")

# Calculate Kappa for Peek Coding ----------------------------------------------
sufev_4_kappa_peek <- kappa2(
  rel_sufev_4[, c("Peek_lil", "Peek")], 
  weight = "unweighted"
)
sufev_4_kappa_peek

# Calculate Kappa for Choice Coding --------------------------------------------
sufev_4_kappa_choice <- kappa2(
  rel_sufev_4[, c("Correct_lil", "Correct")], 
  weight = "unweighted"
)
sufev_4_kappa_choice


# ==============================================================================
# PART 5: SUFFICIENCY EVIDENCE - 6-YEAR-OLD CHILDREN
# ==============================================================================

# Load Reliability Data --------------------------------------------------------
rel_sufev_6 <- read.csv("data/reliability_evsk_6.csv")

# Calculate Kappa for Peek Coding ----------------------------------------------
sufev_6_kappa_peek <- kappa2(
  rel_sufev_6[, c("Peek_lil", "Peek")], 
  weight = "unweighted"
)
sufev_6_kappa_peek

# Calculate Kappa for Choice Coding --------------------------------------------
sufev_6_kappa_choice <- kappa2(
  rel_sufev_6[, c("Correct_lil", "Correct")], 
  weight = "unweighted"
)
sufev_6_kappa_choice


# ==============================================================================
# SUMMARY OF RELIABILITY RESULTS
# ==============================================================================

cat("\n=== Inter-Rater Reliability Summary ===\n")
cat("\nEvidence Seeking - Chimpanzees:\n")
cat("  Peek Kappa:", round(evsk_chimp_peek_kappa$value, 3), "\n")
cat("  Choice Kappa:", round(evsk_chimp_choice_kappa$value, 3), "\n")

cat("\nSufficiency Evidence - Chimpanzees:\n")
cat("  Choice Kappa:", round(sufev_rel_choice$value, 3), "\n")

cat("\nEvidence Seeking - Children:\n")
cat("  Peek Kappa:", round(evsk_kid_kappa_peek$value, 3), "\n")
cat("  Choice Kappa:", round(evsk_kid_kappa_choice$value, 3), "\n")

cat("\nSufficiency Evidence - 4-Year-Olds:\n")
cat("  Peek Kappa:", round(sufev_4_kappa_peek$value, 3), "\n")
cat("  Choice Kappa:", round(sufev_4_kappa_choice$value, 3), "\n")

cat("\nSufficiency Evidence - 6-Year-Olds:\n")
cat("  Peek Kappa:", round(sufev_6_kappa_peek$value, 3), "\n")
cat("  Choice Kappa:", round(sufev_6_kappa_choice$value, 3), "\n")

# ==============================================================================
# END OF SCRIPT
# ==============================================================================