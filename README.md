# Alcohol-attributable and avoidable injury deaths in Chile (PAF / PIF)

R code to estimate **alcohol-attributable fractions (PAF)** and **potential impact fractions (PIF)** for injury mortality in Chile, and to translate them into **attributable and avoidable deaths** and **years of potential life lost (YPLL)** by year, sex and age group.

The analysis combines:

- **Alcohol consumption** from Chilean national drug surveys (SENDA), modelled as gamma distributions of daily grams of pure alcohol, split into drinkers with and without heavy episodic drinking (HED / non-HED).
- **Relative-risk functions** for three injury groups: road injuries, unintentional injuries and intentional injuries.
- **Death records** from DEIS (Ministry of Health, Chile), classified into injury causes using the GBD ICD-10 mapping (Shield et al.).

Two counterfactual scenarios are evaluated, each with 10 %, 20 % and 30 % reductions:

| Scenario | Description |
|---|---|
| **HED shift** (`cf_type = "hed"`) | A fraction of HED drinkers consuming ≤ 60 g/day move to the non-HED group. HED drinkers above 60 g/day stay unchanged. |
| **Volume shift** (`cf_type = "volume"`) | All current drinkers reduce their daily volume by the given percentage. HED/non-HED membership doesn't change. |

Uncertainty is propagated by Monte Carlo simulation (10,000 draws) over the RR coefficients and the exposure prevalences. The results are 95 % uncertainty intervals.

---

## Repository structure

```
PIF addiction/
├── functions.R                 # All statistical functions (no data, no execution)
├── icd10_codes_inj.R           # ICD-10 code lists for injury causes (GBD classification)
├── Mortality injuries.R        # Step 1 – builds the injury mortality dataset from DEIS
├── run_comparison.R            # Step 2 – computes PAF and PIF (main + sensitivity datasets)
├── PIF results.R               # Step 3 – attributable/avoidable deaths, YPLL, tables, figures
├── Data/                       # Input and intermediate data (not versioned, see below)
└── results/                    # Outputs (tables .xlsx, figures .tiff, results .rds)
```

---

## Requirements

- **R ≥ 4.1** (the code uses the native pipe `|>`)
- Packages (installed automatically through `pacman::p_load()` where possible):

```r
install.packages("pacman")
pacman::p_load(
  tidyverse, rio, janitor, readr,          # data handling / IO
  fitdistrplus, MASS, purrr,               # distribution fitting (functions.R)
  RColorBrewer, scales, patchwork,         # figures
  Hmisc                                    # weighted quantiles (in-text statistics)
)
```

---

## Input data

### Mortality data (included in processed form)

The raw DEIS death files are too large for GitHub (~870 MB and ~90 MB), so they are **not** included. The repository ships the processed dataset instead:

| File | Location | Description |
|---|---|---|
| `mortality_injuries_cleanV2.rds` | `PIF addiction/Data/` | Injury deaths, age ≥ 15, 2008–2024, already classified by cause (`year`, `gender`, `comuna`, `region`, `age`, `age_group`, `disease`, `disease_group`). **Included.** |

This file is all that `PIF results.R` needs, so **you can skip Step 1**.

To rebuild it from scratch, download the death records from the [DEIS open data portal](https://deis.minsal.cl/#datosabiertos) (Ministry of Health, Chile) and place them in `PIF addiction/Data/` with these names:

| File | Description |
|---|---|
| `DEFUNCIONES_FUENTE_DEIS_1990_2023_CIFRAS_OFICIALES.csv` | Official death records, 1990–2023 |
| `DEFUNCIONES_FUENTE_DEIS_2024_2026_02062026.csv` | Death records, 2024–2026 (release of 02-06-2026) |

Then run Step 1.

### Other inputs

| File | Location | Source / description |
|---|---|---|
| `data_binge_original.rds` | `PIF addiction/Data/` | Processed SENDA survey data, 2008–2022 (original alcohol definition) |
| `data_binge_sensitivity.rds` | `PIF addiction/Data/` | Processed SENDA survey data, 2012–2024 (**main analysis**) |
| `ine_proyecciones.xlsx` | **repository root** (one level above `PIF addiction/`) | INE population projections by single age; sheet 1 = men, sheet 2 = women; columns `edad`, `ano_YYYY` |

The two `data_binge_*.rds` files must contain at least the following variables: `year`, `sexo` (`"Hombre"`/`"Mujer"`), `edad_tramo` (1–4), `volajohdia` (g of alcohol/day), `volajohdiams`, `hed` (0/1), `cvolaj` (`"ltabs"` = lifetime abstainer, `"fd"` = former drinker), `exp` (survey expansion weight). The sensitivity file also needs `oh1` and `oh2`.

Age groups used throughout: `1` = 15–29, `2` = 30–44, `3` = 45–59, `4` = 60+.

---

## How to run

> **Working directory:** all paths in the scripts are written relative to the **repository root** (the folder that *contains* `PIF addiction/`), e.g. `"PIF addiction/Data/..."`. Open the R project at the root, or run `setwd("<path-to-repo-root>")` before running anything.

Run the scripts in this order, each one in a **fresh R session**:

### 1. Build the injury mortality dataset (optional)

> Skip this step if you use the included `Data/mortality_injuries_cleanV2.rds`. You only need it to rebuild that file from the raw DEIS CSVs (see *Input data*).

```r
source("PIF addiction/Mortality injuries.R")
```

- Reads both DEIS files and keeps deaths aged ≥ 15 from 2008 onward (plus 2024).
- Classifies each death by external cause (`diag2`) with `icd10_codes_inj.R`.
- **Outputs:**
  - `Data/mortality_data_injuriesV2.rds` (intermediate, all deaths)
  - `Data/mortality_injuries_cleanV2.rds` (injury deaths only, with `disease` and `disease_group`)

### 2. Compute PAF and PIF

```r
source("PIF addiction/run_comparison.R")
```

- Sources `functions.R`.
- Runs the full pipeline twice: on `data_binge_original.rds` (survey years 2008–2022) and on `data_binge_sensitivity.rds` (survey years 2012–2024).
- For each sex × age group × survey year it fits gamma distributions (HED / non-HED), computes the PAF and the PIF under both counterfactuals, and draws 10,000 Monte Carlo simulations.
- ⏱ This is the computationally heavy step (can take several hours).
- **Outputs** (in `results/`):
  - `paf_comparisonV2.rds` / `.xlsx`
  - `pif_comparisonV2.rds` / `.xlsx`

### 3. Attributable and avoidable deaths, YPLL, tables and figures

```r
source("PIF addiction/PIF results.R")
```

- Joins PAF/PIF with observed deaths and YPLL (life expectancy: 76 years for men, 82 for women).
- **Main analysis:** `dataset == "sensitivity"` (2012–2024).
- **Sensitivity analysis:** original vs. sensitivity datasets (overlapping years 2012–2022).
- **Outputs** (in `results/`):

| File | Content |
|---|---|
| `tab1_descriptive.xlsx` | Descriptive statistics of current drinkers (HED / non-HED / overall) |
| `tab1_att_main.xlsx` | Median attributable deaths by cause, sex and age group |
| `tab2_avoid_main.xlsx` | Median avoidable deaths by cause, scenario and reduction |
| `fig1_att_total.tiff` | Total vs. alcohol-attributable injury deaths over time |
| `fig2_att_avoided_MAIN_SEX.tiff` | Avoided deaths by sex – HED shift (main) |
| `fig2_att_avoided_SUPP_SEX.tiff` | Avoided deaths by sex – volume shift (supplementary) |
| `fig3_avoidxage.tiff`, `fig3_avoidxageSUPP.tiff` | Avoided deaths by age group over time (HED / volume) |
| `fig3_avoidxage_MAINFOREST.tiff`, `fig3_avoidxage_SUPPFOREST.tiff` | Same, as forest plots (median 2012–2024) |
| `fig4_att_ypll.tiff` | Total vs. attributable YPLL |
| `fig5_ypll_avoid.tiff` | % of YPLL avoided over time (HED shift) |
| `fig5_ypll_avoid_MAIN.tiff`, `fig5_ypll_avoid_SUPP.tiff` | % of YPLL avoided, forest plots (HED / volume) |
| `tab_sens_paf_diff.xlsx`, `tab_sens_avoid_diff.xlsx`, `tab_sens_summary.xlsx` | Sensitivity analysis tables |
| `figS1_paf_og_vs_sens.tiff`, `figS2_avoid_og_vs_sens.tiff` | Sensitivity analysis figures |

The script also prints the in-text statistics (consumption, HED prevalence, totals and percentages) to the console.

---

## Main functions (`functions.R`)

| Function | Purpose |
|---|---|
| `build_cd_hed_list()` | Consumption vectors by year × age group × HED status |
| `fit_gamma_by_tramo()` | Gamma fits of consumption per year and age group |
| `build_s_hed_list_weighted()` | Weighted share of HED among current drinkers |
| `build_abs_list_weighted()` / `build_form_list_weighted()` | Weighted prevalence of lifetime abstainers / former drinkers |
| `paf_hed_function()` / `confint_paf_hed()` / `compute_paf()` | PAF (point, Monte Carlo CI, loop over years and age groups) |
| `pif_hed_function()` / `confint_pif_hed()` / `compute_inj()` | PIF under the HED-shift counterfactual |
| `pif_volume_function()` / `confint_vol_hed()` / `compute_vol()` | PIF under the volume-shift counterfactual |
| `tidy_paf()` / `tidy_pif()` | Reshape wide outputs to long format |

**Relative-risk specification** (`run_comparison.R`):

- Non-HED: `RR(x) = exp(b1 · x)`
- HED: `RR(x) = exp(b1 · x + b2)`

Coefficients used: road injuries `b1 = 0.0029955`, `b2 = 0.95935`. Unintentional and intentional injuries share one set: `b1 = 0.0019980`, `b2 = 0.70310`. Standard errors are in the script.

---

## Notes

- The DEIS CSV files are very large (~870 MB and ~90 MB) and exceed GitHub's 100 MB file limit. Keep `Data/` and heavy `results/` files out of version control, for example with this `.gitignore`:

  ```gitignore
  PIF addiction/Data/
  PIF addiction/results/*.rds
  *.RData
  .Rhistory
  .Rproj.user/
  ```

- Random seeds are fixed (`seed = 145`), so results are reproducible.
- Variable and level names in the data are in Spanish (`sexo`, `Hombre`/`Mujer`, `edad_tramo`). They are recoded to English (`Men`/`Women`) in `PIF results.R`.
