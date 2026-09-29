# ============================================================ #
#  run_comparison.R                                            #
#  Pipeline de cálculo PAF/PIF para data_og y data_sens.      #
#  Exporta resultados comparables en RDS y Excel.             #
#  Punto de entrada único: correr este script completo.       #
# ============================================================ #

rm(list = ls()); gc()

# ── 0. FUNCIONES Y LIBRERÍAS ──────────────────────────────────
source("PIF addiction/functions.R")      # todas las funciones estadísticas
library(rio)

# ── 1. DATOS ─────────────────────────────────────────────────
data_og <- import("PIF addiction/Data/data_binge_original.rds") %>%
  filter(volajohdia <= 150)

data_sens <- import("PIF addiction/Data/data_binge_sensitivity.rds") %>%
  mutate(
    volajohdia   = ifelse(oh1 == "No" | oh2 != "30 dias", 0, volajohdia),
    volajohdiams = ifelse(oh1 == "No" | oh2 != "30 dias", 0, volajohdiams)
  ) %>%
  filter(volajohdia <= 150)

# ── 2. CONSTANTES ─────────────────────────────────────────────
# data_og:   encuestas 2008–2022 (paso 2) → 8 puntos
# data_sens: encuestas 2012–2024 (paso 2) → 7 puntos
years_vec_og   <- c(2008L, 2010L, 2012L, 2014L, 2016L, 2018L, 2020L, 2022L)
years_vec_sens <- c(2012L, 2014L, 2016L, 2018L, 2020L, 2022L, 2024L)

# Parámetros RR — Road Injuries
b1_ri      <- 0.00299550897979837
var_b1_ri  <- 0.000508678216036836^2
b2_ri      <- 0.959350221334602
vcov_hed   <- matrix(c(0.000508678222^2, 0, 0, 0.227875857649849^2), nrow = 2)

# Parámetros RR — Unintentional + Intentional Injuries (misma RR)
b1_inj     <- 0.00199800266267306
var_b1_inj <- 0.000509185879264017^2
b2_inj     <- 0.703097511413113
var_b2_inj <- 0.198243000172775^2

# Funciones RR
rr_linear  <- function(x, b)    exp(b * x)
ri_hed_fun <- function(x, beta) exp(beta[1] * x + beta[2])

# ── 3. PIPELINE ───────────────────────────────────────────────

run_pipeline <- function(data, years_vec, n_sim = 10000, seed = 145) {

  # 3a. Distribuciones de consumo y proporciones ---
  cd_fem  <- build_cd_hed_list(data, "Mujer",  years_vec)
  g_fem   <- fit_gamma_by_tramo(cd_fem)
  p_hed_f <- build_s_hed_list_weighted(data, "Mujer",  years_vec)
  p_abs_f <- build_abs_list_weighted(data,   "Mujer",  years_vec)
  p_frm_f <- build_form_list_weighted(data,  "Mujer",  years_vec)

  cd_male <- build_cd_hed_list(data, "Hombre", years_vec)
  g_male  <- fit_gamma_by_tramo(cd_male)
  p_hed_m <- build_s_hed_list_weighted(data, "Hombre", years_vec)
  p_abs_m <- build_abs_list_weighted(data,   "Hombre", years_vec)
  p_frm_m <- build_form_list_weighted(data,  "Hombre", years_vec)

  # 3b. PAF ------------------------------------------------
  paf_out <- tidy_paf(
    df_fem = bind_rows(
      compute_paf("Mujer",  g_fem,  p_abs_f, p_frm_f, p_hed_f, rr_linear, ri_hed_fun,
                  b1_ri,  var_b1_ri,  b2_ri,  vcov_hed[2,2], n_sim = n_sim, seed = seed,
                  disease = "Road Injuries"),
      compute_paf("Mujer",  g_fem,  p_abs_f, p_frm_f, p_hed_f, rr_linear, ri_hed_fun,
                  b1_inj, var_b1_inj, b2_inj, var_b2_inj,     n_sim = n_sim, seed = seed,
                  disease = "Unintentional Injuries"),
      compute_paf("Mujer",  g_fem,  p_abs_f, p_frm_f, p_hed_f, rr_linear, ri_hed_fun,
                  b1_inj, var_b1_inj, b2_inj, var_b2_inj,     n_sim = n_sim, seed = seed,
                  disease = "Unintentional Injuries") %>% mutate(disease = "Intentional Injuries")
    ),
    df_male = bind_rows(
      compute_paf("Hombre", g_male, p_abs_m, p_frm_m, p_hed_m, rr_linear, ri_hed_fun,
                  b1_ri,  var_b1_ri,  b2_ri,  vcov_hed[2,2], n_sim = n_sim, seed = seed,
                  disease = "Road Injuries"),
      compute_paf("Hombre", g_male, p_abs_m, p_frm_m, p_hed_m, rr_linear, ri_hed_fun,
                  b1_inj, var_b1_inj, b2_inj, var_b2_inj,     n_sim = n_sim, seed = seed,
                  disease = "Unintentional Injuries"),
      compute_paf("Hombre", g_male, p_abs_m, p_frm_m, p_hed_m, rr_linear, ri_hed_fun,
                  b1_inj, var_b1_inj, b2_inj, var_b2_inj,     n_sim = n_sim, seed = seed,
                  disease = "Unintentional Injuries") %>% mutate(disease = "Intentional Injuries")
    )
  )

  # 3c. Helpers internos para los loops de PIF ---
  run_hed <- function(sexo, g, p_hed, b1, vb1, b2, vb2, dis)
    lapply(c(`10%` = 0.9, `20%` = 0.8, `30%` = 0.7), function(sv)
      compute_inj(sexo, g, p_hed, sv, rr_linear, ri_hed_fun,
                  b1, vb1, b2, vb2, n_sim, seed, dis)) %>% bind_rows()

  run_vol <- function(sexo, g, p_hed, b1, vb1, b2, vb2, dis)
    lapply(c(`10%` = 0.9, `20%` = 0.8, `30%` = 0.7), function(vv)
      compute_vol(sexo, g, p_hed, vv, rr_linear, ri_hed_fun,
                  b1, vb1, b2, vb2, n_sim, seed, dis)) %>% bind_rows()

  rename_fem  <- function(df)
    rename(df, PIF_ag1 = Fem1_point, LL_ag1 = Fem1_lower, UL_ag1 = Fem1_upper,
               PIF_ag2 = Fem2_point, LL_ag2 = Fem2_lower, UL_ag2 = Fem2_upper,
               PIF_ag3 = Fem3_point, LL_ag3 = Fem3_lower, UL_ag3 = Fem3_upper,
               PIF_ag4 = Fem4_point, LL_ag4 = Fem4_lower, UL_ag4 = Fem4_upper)

  rename_male <- function(df)
    rename(df, PIF_ag1 = Male1_point, LL_ag1 = Male1_lower, UL_ag1 = Male1_upper,
               PIF_ag2 = Male2_point, LL_ag2 = Male2_lower, UL_ag2 = Male2_upper,
               PIF_ag3 = Male3_point, LL_ag3 = Male3_lower, UL_ag3 = Male3_upper,
               PIF_ag4 = Male4_point, LL_ag4 = Male4_lower, UL_ag4 = Male4_upper)

  build_pif <- function(run_fn, rename_f, rename_m) {
    bind_rows(
      tidy_pif(rename_f(bind_rows(
        run_fn("Mujer",  g_fem,  p_hed_f, b1_ri,  var_b1_ri,  b2_ri,  vcov_hed[2,2], "Road Injuries"),
        run_fn("Mujer",  g_fem,  p_hed_f, b1_inj, var_b1_inj, b2_inj, var_b2_inj,    "Unintentional Injuries"),
        run_fn("Mujer",  g_fem,  p_hed_f, b1_inj, var_b1_inj, b2_inj, var_b2_inj,    "Unintentional Injuries") %>%
          mutate(disease = "Intentional Injuries")
      )), "Mujer"),
      tidy_pif(rename_m(bind_rows(
        run_fn("Hombre", g_male, p_hed_m, b1_ri,  var_b1_ri,  b2_ri,  vcov_hed[2,2], "Road Injuries"),
        run_fn("Hombre", g_male, p_hed_m, b1_inj, var_b1_inj, b2_inj, var_b2_inj,    "Unintentional Injuries"),
        run_fn("Hombre", g_male, p_hed_m, b1_inj, var_b1_inj, b2_inj, var_b2_inj,    "Unintentional Injuries") %>%
          mutate(disease = "Intentional Injuries")
      )), "Hombre")
    )
  }

  # 3d. PIF HED shift y PIF volume shift ---
  pif_hed_out <- build_pif(run_hed, rename_fem, rename_male)
  pif_vol_out <- build_pif(run_vol, rename_fem, rename_male)

  list(paf = paf_out, pif_hed = pif_hed_out, pif_vol = pif_vol_out)
}

# ── 4. EJECUTAR PARA AMBOS DATASETS ──────────────────────────
message(">>> Pipeline: data_og")
res_og   <- run_pipeline(data_og,   years_vec = years_vec_og,   n_sim = 10000, seed = 145)

message(">>> Pipeline: data_sens")
res_sens <- run_pipeline(data_sens, years_vec = years_vec_sens, n_sim = 10000, seed = 145)

# ── 5. COMBINAR ───────────────────────────────────────────────
paf_comp <- bind_rows(
  res_og$paf   %>% mutate(dataset = "original"),
  res_sens$paf %>% mutate(dataset = "sensitivity")
)

pif_comp <- bind_rows(
  res_og$pif_hed   %>% mutate(dataset = "original"),
  res_og$pif_vol   %>% mutate(dataset = "original"),
  res_sens$pif_hed %>% mutate(dataset = "sensitivity"),
  res_sens$pif_vol %>% mutate(dataset = "sensitivity")
)

# ── 6. EXPORTAR ───────────────────────────────────────────────
dir.create("PIF addiction/results", showWarnings = FALSE, recursive = TRUE)

write_rds(paf_comp, "PIF addiction/results/paf_comparisonV2.rds")
write_rds(pif_comp, "PIF addiction/results/pif_comparisonV2.rds")

export(
  paf_comp %>%
    mutate(value = paste0(round(point, 3), " (", round(lower, 3), ", ", round(upper, 3), ")")) %>%
    dplyr::select(dataset, year, gender, age_group, disease, value) %>%
    arrange(dataset, disease, gender, age_group, year),
  "PIF addiction/results/paf_comparisonV2.xlsx"
)

export(
  pif_comp %>%
    mutate(value = paste0(round(pif, 3), " (", round(pif_low, 3), ", ", round(pif_up, 3), ")")) %>%
    dplyr::select(dataset, cf_type, shift, year, gender, age_group, disease, value) %>%
    arrange(dataset, cf_type, shift, disease, gender, age_group, year),
  "PIF addiction/results/pif_comparisonV2.xlsx"
)

# ── 7. VERIFICACIÓN RÁPIDA ────────────────────────────────────
paf_diff <- paf_comp %>%
  dplyr::select(dataset, year, gender, age_group, disease, point) %>%
  pivot_wider(names_from = dataset, values_from = point) %>%
  mutate(diff = original - sensitivity)

pif_diff <- pif_comp %>%
  dplyr::select(dataset, cf_type, shift, year, gender, age_group, disease, pif) %>%
  pivot_wider(names_from = dataset, values_from = pif) %>%
  mutate(diff = original - sensitivity)

cat("\n--- Diferencias PAF (original - sensitivity) ---\n")
print(paf_diff %>% summarise(max = max(abs(diff), na.rm = TRUE),
                             mean = mean(abs(diff), na.rm = TRUE)))

cat("\n--- Diferencias PIF por cf_type y shift ---\n")
print(pif_diff %>% group_by(cf_type, shift) %>%
        summarise(max = max(abs(diff), na.rm = TRUE),
                  mean = mean(abs(diff), na.rm = TRUE), .groups = "drop"))

message(">>> Listo. Resultados en PIF addiction/results/")
