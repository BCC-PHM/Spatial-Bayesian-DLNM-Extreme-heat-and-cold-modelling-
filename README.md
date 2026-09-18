# **Spatial Bayesian Distributed Lag Non-linear Modelling of Temperature-Related Mortality in Birmingham**

An application of the spatial Bayesian Distributed Lag Non-linear models (SB-DLNMs) introduced by **Quijal-Zamorano et al. (2024)**.

The SB-DLNMs estimate temperature-mortality relationships at fine geographical level with limited data, allowing the possibility for local authorities to implement the same method without relying only on national level estimates.

## Method

Analysis is conducted at **2022 electoral ward level across Birmingham**.

1.  **Temperature exposure**\
    Area-weighted average at Birmingham 2022 ward boundaries using the mean of daily maximum and minimum values of air temperature on 1 km grid across UK obtained from the HadUK-Grid database.
2.  **Death register**\
    Daily register of deaths held by Birmingham City Council (BCC) is not publicly available.

Statistics employed:

1.  Time-stratified case-crossover design using Poisson regression
2.  Distributed lag non-linear model (DLNM) to capture the delayed and non-linear temperature-mortality relationship
3.  Bayesian hierarchical spatial modelling using the BYM2 model, allowing ward's exposure-response curve to borrow strength from adjacent neighbours
4.  Integrated Nested Laplace Approximation (INLA) for model fitting
5.  Ward-specific minimum mortality temperature (MMT) estimated as a posterior distribution
6.  Backward attributable fraction framework to calculate attributable numbers and excess mortality
7.  Posterior summaries reported as medians with 95% credible intervals and exceedance probabilities

## Data

| Measure                         | Dataset / source         |      Year | Role                     |
|------------------|------------------|-----------------:|------------------|
| Daily mean temperature          | HadUK-Grid at 1km        | 2005-2025 | Exposure                 |
| Daily deaths                    | BCC death register       | 2005–2025 | Outcome                  |
| Ward-level population estimates | ONS                      |  mid-2022 | Rate denominator         |
| Spatial boundaries              | ONS ward 2022 boundaries |      2022 | Spatial unit of analysis |

## Main code

-   `01_raster_to_time_series_ward.R` — converts raster data into an area-weighted daily time series at 2022 ward level
-   `02_all_cause_mortality_process_ward.R` — processes the BCC all-cause death register into daily ward-level counts and builds the case-crossover strata
-   `03_bayesain_spatial_dnlm.R` — fits the spatial Bayesian DLNM (BYM2 on the cross-basis coefficients) in INLA and draws posterior coefficients
-   `04_plot_inla_res.R` — plots the ward-specific exposure-response curves and derives the posterior MMT draws
-   `05_excess_mort_attr.R` — calculates attributable fractions and annual excess deaths by ward
-   `06_DPH_mmt_draws_0_3.R` — MMT draws specifically for the DPH report
-   `06a_DPH_raster_to_time_series_ward.R` — ward temperature series for the DPH report period using provisional HadUK-Grid
-   `07_DPH_efficient_excess_mort_attr.R` — excess mortality attribution for the DPH report
-   `08_baseline_mort_ward.R` — ward-level baseline mortality
-   `09_DPH_heatwave_em2026.R` — excess mortality plots for the June and July 2026 heatwave
-   `10_DPH_report_static_graphs.R` — non-interactive plots made for DPH report

## Output

-   All plots are in the `figs` folder

## Requirements

The analysis was conducted in **R**. Relevant packages include `tidyverse`, `sf`, `exactextractr`, `dlnm`, `INLA`, `spdep`, `tmap`, `readxl`, `doParallel` and `foreach`.

## Reference

Quijal-Zamorano M, Martinez-Beneito MA, Ballester J, Marí-Dell’Olmo M. Spatial bayesian distributed lag non-linear models (SB-DLNM) for small-area exposure-lag-response epidemiological modelling. International Journal of Epidemiology 2024;53:dyae061. <https://doi.org/10.1093/ije/dyae061>.
