# ============================================================ #
#  run_comparison_reclassify.R                                 #
#  Prueba de concepto: volume shift con reclasificación        #
#  HED → NHED cuando el consumo desplazado cae bajo el umbral. #
#                                                              #
#  Solo corre data_sens. Compara pif_vol original vs           #
#  pif_vol_reclassify. NO modifica run_comparison.R ni         #
#  functions.R.                                                #
#                                                              #
#  Si el resultado convence → reemplazar en run_comparison.R   #
# ============================================================ #

rm(list = ls()); gc()
source("PIF addiction/functions.R")
library(rio)

# ── 1. DATOS ─────────────────────────────────────────────────
data_sens <- import("PIF addiction/Data/data_binge_sensitivity.rds") %>%
  mutate(
    volajohdia   = ifelse(oh1 == "No" | oh2 != "30 dias", 0, volajohdia),
    volajohdiams = ifelse(oh1 == "No" | oh2 != "30 dias", 0, volajohdiams)
  ) %>%
  filter(volajohdia <= 150)

years_vec_sens <- c(2012L, 2014L, 2016L, 2018L, 2020L, 2022L, 2024L)

# ── 2. CONSTANTES ─────────────────────────────────────────────
b1_ri      <- 0.00299550897979837
var_b1_ri  <- 0.000508678216036836^2
b2_ri      <- 0.959350221334602
vcov_hed   <- matrix(c(0.000508678222^2, 0, 0, 0.227875857649849^2), nrow = 2)

b1_inj     <- 0.00199800266267306
var_b1_inj <- 0.000509185879264017^2
b2_inj     <- 0.703097511413113
var_b2_inj <- 0.198243000172775^2

rr_linear  <- function(x, b)    exp(b * x)
ri_hed_fun <- function(x, beta) exp(beta[1] * x + beta[2])

# ── 3. NUEVAS FUNCIONES ───────────────────────────────────────

# Contrafáctico: reduce consumo x vol_shift.
# Los bebedores HED cuyo consumo desplazado (x * vol_shift) cae bajo
# el umbral 'cut' se reclasifican como NHED → pierden el riesgo b2.
pif_volume_reclassify_function <- function(x, fit_nhed, fit_hed,
                                            rr_fun_nhed, rr_fun_hed,
                                            b, s_hed, vol_shift = 1, cut = 60) {
  trap_int <- function(x, f) {
    if (length(x) < 2L) return(0)
    dx <- diff(x); sum((f[-1] + f[-length(f)]) * dx / 2)
  }
  dens_from_fit <- function(fit, x) {
    est <- fit$estimate
    if ("rate"  %in% names(est)) return(dgamma(x, shape = est[["shape"]], rate  = est[["rate"]]))
    if ("scale" %in% names(est)) return(dgamma(x, shape = est[["shape"]], scale = est[["scale"]]))
    stop("fit$estimate debe tener 'rate' o 'scale'.")
  }
  y_nhed   <- dens_from_fit(fit_nhed, x)
  y_hed    <- dens_from_fit(fit_hed,  x)
  par_nhed <- if (is.list(b)) b$nhed else b
  par_hed  <- if (is.list(b)) b$hed  else b

  idx_nhed <- x <= cut
  Z_nhed   <- trap_int(x[idx_nhed], y_nhed[idx_nhed])
  Z_hed    <- trap_int(x, y_hed)
  if (!is.finite(Z_nhed) || Z_nhed <= 0 || !is.finite(Z_hed) || Z_hed <= 0) return(NA_real_)

  # Riesgo observado
  R_nhed_obs <- trap_int(x[idx_nhed], (y_nhed[idx_nhed] / Z_nhed) * rr_fun_nhed(x[idx_nhed], par_nhed))
  R_hed_obs  <- trap_int(x,           (y_hed  / Z_hed)            * rr_fun_hed( x,            par_hed))

  # Contrafáctico NHED: consumo desplazado, misma función RR NHED
  R_nhed_cf  <- trap_int(x[idx_nhed],
                         (y_nhed[idx_nhed] / Z_nhed) * rr_fun_nhed(x[idx_nhed] * vol_shift, par_nhed))

  # Contrafáctico HED: los que caen bajo el umbral usan RR NHED en su consumo desplazado;
  # los que permanecen sobre el umbral usan RR HED en su consumo desplazado.
  rr_hed_cf <- ifelse(x * vol_shift <= cut,
                      rr_fun_nhed(x * vol_shift, par_nhed),   # reclasificados → RR NHED
                      rr_fun_hed( x * vol_shift, par_hed))    # permanecen HED → RR HED
  R_hed_cf   <- trap_int(x, (y_hed / Z_hed) * rr_hed_cf)

  s_nhed <- 1 - s_hed
  R_obs  <- s_nhed * R_nhed_obs + s_hed * R_hed_obs
  R_cf   <- s_nhed * R_nhed_cf  + s_hed * R_hed_cf

  if (!is.finite(R_obs) || R_obs <= 0) return(NA_real_)
  round(1 - R_cf / R_obs, 3)
}

confint_vol_reclassify <- function(x, fit_nhed, fit_hed,
                                    rr_fun_nhed, rr_fun_hed,
                                    beta_mean_nhed, beta_var_nhed,
                                    beta_mean_b2,   beta_var_b2,
                                    s_hed, vol_shift = 1, n_sim = 5000, seed = 145) {
  set.seed(seed)
  sims <- numeric(n_sim)
  for (i in seq_len(n_sim)) {
    b1_i <- if (beta_var_nhed > 0) rnorm(1, beta_mean_nhed, sqrt(beta_var_nhed)) else beta_mean_nhed
    b2_i <- if (beta_var_b2   > 0) rnorm(1, beta_mean_b2,   sqrt(beta_var_b2))   else beta_mean_b2
    sims[i] <- tryCatch(
      pif_volume_reclassify_function(
        x, fit_nhed, fit_hed, rr_fun_nhed, rr_fun_hed,
        b = list(nhed = b1_i, hed = c(b1_i, b2_i)),
        s_hed = s_hed, vol_shift = vol_shift
      ),
      error = function(e) NA_real_
    )
  }
  sims <- sims[is.finite(sims)]
  if (!length(sims)) stop("Todas las simulaciones fallaron.")
  list(point_estimate = round(mean(sims), 3),
       lower_ci       = round(max(0, quantile(sims, 0.025)), 3),
       upper_ci       = round(min(1, quantile(sims, 0.975)), 3),
       sims           = sims)
}

compute_vol_reclassify <- function(sexo = c("Mujer", "Hombre"),
                                    g_hed_list, p_hed_list, vol_shift_val,
                                    rr_fun_nhed, rr_fun_hed,
                                    beta_mean_nhed, beta_var_nhed,
                                    beta_mean_b2,   beta_var_b2,
                                    n_sim = 10000, seed = 145,
                                    disease = "Road Injuries") {
  sexo      <- match.arg(sexo)
  prefix    <- if (sexo == "Mujer") "Fem" else "Male"
  x         <- seq(0.1, 150, length.out = 1500)
  years_vec <- as.integer(names(g_hed_list))
  df <- data.frame(
    Year = years_vec,
    setNames(as.list(rep(NA_real_, 12)),
             c(paste0(prefix, 1:4, rep(c("_point","_lower","_upper"), each = 4))))
  )
  for (i in seq_along(df$Year)) {
    year_i <- as.character(df$Year[i])
    for (j in 1:4) {
      fit_nhed <- g_hed_list[[year_i]][[j]]$nhed
      fit_hed  <- g_hed_list[[year_i]][[j]]$hed
      s_hed    <- p_hed_list[[j]][i]
      if (is.null(fit_nhed) || is.null(fit_hed) || is.na(s_hed)) {
        cat("Skipping PIF vol reclassify — Year:", df$Year[i], sexo, "tramo", j, "\n"); next
      }
      res <- tryCatch(
        confint_vol_reclassify(x, fit_nhed, fit_hed, rr_fun_nhed, rr_fun_hed,
                               beta_mean_nhed, beta_var_nhed, beta_mean_b2, beta_var_b2,
                               s_hed, vol_shift_val, n_sim, seed),
        error = function(e) {
          cat("Error PIF vol reclassify:", df$Year[i], sexo, j, "->", e$message, "\n"); NULL
        }
      )
      if (!is.null(res)) {
        df[i, paste0(prefix, j, "_point")] <- res$point_estimate
        df[i, paste0(prefix, j, "_lower")] <- res$lower_ci
        df[i, paste0(prefix, j, "_upper")] <- res$upper_ci
      }
    }
  }
  df$disease <- disease
  df$shift   <- paste0(((1 - vol_shift_val) * 100), "%")
  df$cf_type <- "volume_reclassify"
  df
}

# ── 4. DATOS DE CONSUMO (solo data_sens) ─────────────────────
message(">>> Construyendo distribuciones — data_sens")

cd_fem  <- build_cd_hed_list(data_sens, "Mujer",  years_vec_sens)
g_fem   <- fit_gamma_by_tramo(cd_fem)
p_hed_f <- build_s_hed_list_weighted(data_sens, "Mujer",  years_vec_sens)

cd_male <- build_cd_hed_list(data_sens, "Hombre", years_vec_sens)
g_male  <- fit_gamma_by_tramo(cd_male)
p_hed_m <- build_s_hed_list_weighted(data_sens, "Hombre", years_vec_sens)

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

# ── 5. PIF VOLUME ORIGINAL (referencia) ──────────────────────
message(">>> Calculando PIF volume original (referencia)...")

run_vol_orig <- function(sexo, g, p_hed, b1, vb1, b2, vb2, dis)
  lapply(c(`10%` = 0.9, `20%` = 0.8, `30%` = 0.7), function(vv)
    compute_vol(sexo, g, p_hed, vv, rr_linear, ri_hed_fun,
                b1, vb1, b2, vb2, n_sim = 5000, seed = 145, disease = dis)) %>%
  bind_rows()

pif_vol_orig <- bind_rows(
  tidy_pif(rename_fem(bind_rows(
    run_vol_orig("Mujer",  g_fem,  p_hed_f, b1_ri,  var_b1_ri,  b2_ri,  vcov_hed[2,2], "Road Injuries"),
    run_vol_orig("Mujer",  g_fem,  p_hed_f, b1_inj, var_b1_inj, b2_inj, var_b2_inj,    "Unintentional Injuries"),
    run_vol_orig("Mujer",  g_fem,  p_hed_f, b1_inj, var_b1_inj, b2_inj, var_b2_inj,    "Unintentional Injuries") %>%
      mutate(disease = "Intentional Injuries")
  )), "Mujer"),
  tidy_pif(rename_male(bind_rows(
    run_vol_orig("Hombre", g_male, p_hed_m, b1_ri,  var_b1_ri,  b2_ri,  vcov_hed[2,2], "Road Injuries"),
    run_vol_orig("Hombre", g_male, p_hed_m, b1_inj, var_b1_inj, b2_inj, var_b2_inj,    "Unintentional Injuries"),
    run_vol_orig("Hombre", g_male, p_hed_m, b1_inj, var_b1_inj, b2_inj, var_b2_inj,    "Unintentional Injuries") %>%
      mutate(disease = "Intentional Injuries")
  )), "Hombre")
) %>% mutate(gender = recode(gender, "Mujer" = "Women", "Hombre" = "Men"),
             cf_type = "volume_original")

# ── 6. PIF VOLUME CON RECLASIFICACIÓN ─────────────────────────
message(">>> Calculando PIF volume con reclasificación HED→NHED...")

run_vol_recl <- function(sexo, g, p_hed, b1, vb1, b2, vb2, dis)
  lapply(c(`10%` = 0.9, `20%` = 0.8, `30%` = 0.7), function(vv)
    compute_vol_reclassify(sexo, g, p_hed, vv, rr_linear, ri_hed_fun,
                           b1, vb1, b2, vb2, n_sim = 5000, seed = 145, disease = dis)) %>%
  bind_rows()

pif_vol_recl <- bind_rows(
  tidy_pif(rename_fem(bind_rows(
    run_vol_recl("Mujer",  g_fem,  p_hed_f, b1_ri,  var_b1_ri,  b2_ri,  vcov_hed[2,2], "Road Injuries"),
    run_vol_recl("Mujer",  g_fem,  p_hed_f, b1_inj, var_b1_inj, b2_inj, var_b2_inj,    "Unintentional Injuries"),
    run_vol_recl("Mujer",  g_fem,  p_hed_f, b1_inj, var_b1_inj, b2_inj, var_b2_inj,    "Unintentional Injuries") %>%
      mutate(disease = "Intentional Injuries")
  )), "Mujer"),
  tidy_pif(rename_male(bind_rows(
    run_vol_recl("Hombre", g_male, p_hed_m, b1_ri,  var_b1_ri,  b2_ri,  vcov_hed[2,2], "Road Injuries"),
    run_vol_recl("Hombre", g_male, p_hed_m, b1_inj, var_b1_inj, b2_inj, var_b2_inj,    "Unintentional Injuries"),
    run_vol_recl("Hombre", g_male, p_hed_m, b1_inj, var_b1_inj, b2_inj, var_b2_inj,    "Unintentional Injuries") %>%
      mutate(disease = "Intentional Injuries")
  )), "Hombre")
) %>% mutate(gender = recode(gender, "Mujer" = "Women", "Hombre" = "Men"))

# ── 7. COMPARACIÓN ────────────────────────────────────────────
comp <- bind_rows(pif_vol_orig, pif_vol_recl) %>%
  group_by(cf_type, shift, disease) %>%
  summarise(
    pif_median = round(median(pif, na.rm = TRUE), 4),
    pif_min    = round(min(pif,    na.rm = TRUE), 4),
    pif_max    = round(max(pif,    na.rm = TRUE), 4),
    .groups = "drop"
  )

message("\n>>> Comparación PIF: volume original vs volume con reclasificación\n")
print(comp, n = Inf)

# Guardar resultados de la prueba
dir.create("PIF addiction/results", showWarnings = FALSE, recursive = TRUE)
write_rds(pif_vol_recl,  "PIF addiction/results/pif_vol_reclassify_test.rds")
export(comp, "PIF addiction/results/pif_vol_reclassify_comparison.xlsx")

message("\n>>> Listo. Resultados en PIF addiction/results/")
