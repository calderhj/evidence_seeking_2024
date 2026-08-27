# Evidence Seeking 2024-2026
Chimpanzees and Young Children Resolve Ambiguity Through Evidence Search

Experiment 1: Chimpanzees and 4-year-old children are presented with ambiguous or unambiguous evidence as to the location of a reward. After, they have the opportunity
to peek over an occluder at indirect trace evidence. In this study peeking is defined as looking over the occluder. The "peek" variable is based on whether or not they peek.
This analysis explores whether their likelihood to peek is predicted condition, and whether they choose correctly at different rates between the two conditions. The chimpanzees 
also have two analyses of whether peeking influenced their likelihood of a correct choice and whether they were more side biased in the ambiguous evidence condition.

Experiment 2: Chimpanzees, 4-year-old children and 6-year-old children are able to seek evidence at two different levels. The lower level (low-peek) provides partial evidence
while the upper level (high-peek) provides conclusive, discriminating evidence. The rewards they must find are either the same size (not requiring a high-peek to make the right choice) or 
two different sizes (high peek required to tell the difference). The peek variable is defined by whether or not the participant high-peeked on a given trial. The analyses explore whether chimpanzees and human
children were more likely to high-peek in the size-inconsistent condition.

### Getting Started
Each script refers to the folder "submission_code" as the working directory. Data, scripts and images can be found in their labeled folders.

### Script Overview
All analyses were run on R version 4.5.1. Each script requires the packages tidyverse and lme4 in addition to an included bootstrap function.

Within the 'scripts' folder, there are two scripts for Experiment 1, one for chimpanzee data and another for data with 4-year-old human children.
There are three scripts for Experiment 2, one for chimpanzee data, one for 4-year-old data and one for 6-year-old data.
These five scripts all include data cleaning, preregistered models, bootstrapping and data visualization.
Another script has the code for all inter-rater reliability tests performed across all instances of Experiment 1 and 2.

Scripts, 'evsk_4s_submission.R', 'evsk_chimps_submission.R', 'sufev_4s_submission.R', 'sufev_6s_submission.R', and 'sufev_chimps_submission.R' each take one .csv data file from the 'data' folder. Each data file is titled
with the experiment (evsk = experiment 1, sufev = experiment 2), and the population tested (chimpanzees, four-year-olds or six-year-olds).

  Experiment 1: Chimpanzees
  Load the data file "evsk_chimp_data.csv"

  Experiment 2: Chimpanzees
  Load the data file "sufev_chimp_data.csv"
  
  Experiment 1: 4s
  Load the data file "evsk_4_data.csv"
  
  Experiment 2: 4s:
  Load the data file "sufev_4_data.csv"
  
  Experiment 2: 6s
  Load the data file "sufev_6_data.csv"
  
Once the data is loaded, the script will take the user through data cleaning, modelling, bootstrapping and figure creation for the main analyses.
Experiment 1 with chimpanzees and 4-year-olds includes extra analyses of box choice and side biasing (chimpanzees only). These analyses use the same data loaded at the beginning of the script.


