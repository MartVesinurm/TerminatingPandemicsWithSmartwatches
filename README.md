Welcome to the  GitHub repository for the paper 'Terminating Pandemics with Smartwatches', currently in peer-review. This repository contains the scripts, and additional information to facilitate replication, evaluation, and further experimentation based on the findings of our research.

In the paper in question we developed a multiscale modeling framework that integrates within-host viral dynamics and between-host interactions to estimate the risk of viral disease outbreaks within a given population. We used the model to evaluate the population-level effectiveness of smartwatch detection in reducing the transmission of three COVID-19 variants and seasonal and pandemic influenza. This repository contains the R scripts to replicate each of the individual results / figures. The results presented in the paper were generated using RStudio (2024.04.2 Build 764), Posit Software, PBC.


Publication Reference:

Authors: Märt Vesinurm, Martial Ndeffo-Mbah, Dan Yamin, Margaret L. Brandeau.

Published in: Currently in peer-review.

DOI: TBA

*************************************
All script are mostly self-contained, however, you might have to install packages to have access to all the used libraries. Please refer to the manuscript for ellaboration on baseline assumptions and arguments behind them.

To generate the main results under baseline assumptions as they are presented in Figure 2 and Table 2 in the manuscript, run (source) the R script "Figure 2". Then, to generate Figure 2, you must simply run the ggplot object "Figure2". The numerical results are stored in the following dataframes: **Covids_66**(For Ancestral, Delta, and Omircron strains of COVID-19 assuming 66% withdrawal rate), **Covids_75** (For Ancestral, Delta, and Omircron strains of COVID-19 assuming 75% withdrawal rate), **Influenzas_66** (For Pandemic and Seasonal strains of Influenza assuming 66% withdrawal rate), and **Influenzas_75** (For Pandemic and Seasonal strains of Influenza assuming 75% withdrawal rate).

To generate the results presenting sensitivity analyses of probability of achieving an effective reproduction number below 1 over smartwatch detection days, which are presented in Figure 3 of the manuscript, run (source) the R script "Figure 3". Then, to generate Figure 3, you must simply run the ggplot object "Figure3". The numerical results are stored in the following dataframes: **Uncertainty_covid**for Ancestral, Delta, and Omircron strains of COVID-19 and **Uncertainty_influenza** for the Pandemic and Seasonal strains of Influenza.

To generate the results presenting the required reduction in social contacts for 50% and 75% probability of achieving an effective reproduction number below 1 over smartwatch detection days, which are presented in Figure 4 of the manuscript, run (source) the R script "Figure 4". Note that this is a very slow and computationally heavy script. Then, to generate Figure 4, you must simply run the ggplot object "Figure4". The numerical results are stored in the following dataframes: 
**withdrawal_requirement_covids_50** for Ancestral, Delta, and Omircron strains of COVID-19 showing the withdrawal requirement to achieve 50% probablity of disease elimination.
**withdrawal_requirement_covids_75** for Ancestral, Delta, and Omircron strains of COVID-19 showing the withdrawal requirement to achieve 75% probablity of disease elimination.
**withdrawal_requirement_influenzas_50** for the Pandemic and Seasonal strains of Influenza showing the withdrawal requirement to achieve 50% probablity of disease elimination.
**withdrawal_requirement_influenzas_75** for the Pandemic and Seasonal strains of Influenza showing the withdrawal requirement to achieve 75% probablity of disease elimination.

Finally, to generate the results presenting the percentage reduction in the reproduction number given variation in day of detecion and rate of withdrawal, which are presented in Figure 5 of the manuscript, run (source) the R script "Figure 5". Then, to generate Figure 5, you must simply run the ggplot object "Figure5". The numerical results are stored in the following dataframes: 
**ancestral_combinations**, which shows the %-reduction in R for the Ancestral strain of COVID-19 given different rates of withdrawal and days of detection.
**delta_combinations**, which shows the %-reduction in R for the Delta strain of COVID-19 given different rates of withdrawal and days of detection.
**omicron_combinations**, which shows the %-reduction in R for the Omicron strain of COVID-19 given different rates of withdrawal and days of detection.
**pandemic_combinations**, which shows the %-reduction in R for the Pandemic strain of Influenza given different rates of withdrawal and days of detection.
**seasonal_combinations**, which shows the %-reduction in R for the Seasonal strain of Influenza given different rates of withdrawal and days of detection.
