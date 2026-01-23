# Evidence Seeking 2024-2026
Ngamba summer 2024 and Berkeley 2024-2025 study. Chimpanzees and young children search for indirect evidence and understand the sufficiency of partial evidence.

### Quick Start
To get started, clone this repository, e.g., execute the following in the terminal:

```{bash}
git clone https://github.com/calderhj/evidence_seeking_2024.git
```

### Script Overview
All analyses were run on R version 4.5.1. Each script requires the packages tidyverse and lme4 in addition to an included bootstrap function.

Within the 'scripts' folder, there are two scripts for Experiment 1, one for chimpanzee data and another for data with 4-year-old human children.
There are three scripts for Experiment 2, one for chimpanzee data, one for 4-year-old data and one for 6-year-old data.
These five scripts all include data cleaning, preregistered models, bootstrapping and data visualization.
Another script has the code for all inter-rater reliability tests performed across all instances of Experiment 1 and 2.

Scripts, 'evsk_4s.R', 'evsk_chimps.R', 'sufev_4s.R', 'sufev_6s.R', and 'sufev_chimps.R' each take one .csv data file from the 'data' folder. Each data file is titled
with the experiment (evsk = experiment 1, sufev = experiment 2), and the population tested (chimpanzees, four-year-olds or six-year-olds).



