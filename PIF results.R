# ============================================================ #
#  PIF results.R                                               #
#  Análisis de muertes atribuibles y evitables por alcohol.   #
#                                                              #
#  MAIN ANALYSIS   : data_sens, compara HED shift vs Volume   #
#  SENSITIVITY     : data_og vs data_sens                     #
#                                                              #
#  Requiere haber corrido run_comparison.R previamente.       #
# ============================================================ #

rm(list = ls()); gc()
pacman::p_load(tidyverse, rio, RColorBrewer, scales, patchwork, janitor)

# ── 0. DATOS DE ENTRADA ───────────────────────────────────────

paf_comp <- read_rds("PIF addiction/results/paf_comparisonV2.rds") %>%
  mutate(gender = recode(gender, "Mujer" = "Women", "Hombre" = "Men"))

pif_comp <- read_rds("PIF addiction/results/pif_comparisonV2.rds") %>%
  mutate(gender = recode(gender, "Mujer" = "Women", "Hombre" = "Men"))

# Mortalidad individual
def <- import("PIF addiction/Data/mortality_injuries_cleanV2.rds") %>%
  mutate(
    gender = recode(gender, "Hombre" = "Men", "Mujer" = "Women")
  )

# YPLL y muertes por grupo
ypll_mort <- def %>%
  mutate(life_exp = ifelse(gender == "Men", 76, 82),
         ypll_ind = pmax(life_exp - age, 0)) %>%
  group_by(year, gender, age_group, disease) %>%
  summarise(deaths = n(), ypll = sum(ypll_ind), .groups = "drop") %>%
  mutate(disease = case_when(
    disease %in% c("Interpersonal Violence", "Self Harm") ~ "Intentional Injuries",
    disease == "Road Injuries"                             ~ "Road Injuries",
    TRUE                                                   ~ "Unintentional Injuries"
  )) %>%
  group_by(year, gender, age_group, disease) %>%
  summarise(deaths = sum(deaths), ypll = sum(ypll), .groups = "drop")

# Proyecciones INE (población ≥ 15 años)
load_ine <- function(sheet, gdr) {
  import("ine_proyecciones.xlsx", sheet = sheet) %>%
    mutate(edad = as.numeric(ifelse(edad == "100+", 100, edad))) %>%
    filter(edad >= 15) %>%
    pivot_longer(starts_with("ano_"), names_to = "year", values_to = "pop") %>%
    group_by(year) %>%
    summarise(pop = sum(pop), .groups = "drop") %>%
    mutate(year   = as.numeric(gsub("ano_", "", year)),
           gender = gdr)
}
pop <- bind_rows(load_ine(1, "Men"), load_ine(2, "Women")) %>%
  filter(year != 2017)

# ── HELPER: construir att a partir de un PAF filtrado ─────────
build_att <- function(paf_df) {
  ypll_mort %>%
    left_join(paf_df, by = c("year", "gender", "age_group", "disease")) %>%
    mutate(
      att         = deaths * point,
      att_low     = deaths * lower,
      att_up      = deaths * upper,
      ypll_att    = ypll   * point,
      ypll_att_low = ypll  * lower,
      ypll_att_up  = ypll  * upper
    )
}

# ── HELPER: construir muertes y YPLL evitables ────────────────
build_avoid <- function(att_df, pif_df) {
  att_df %>%
    dplyr::select(year, gender, age_group, disease, deaths, ypll,
                  att, att_low, att_up, ypll_att, ypll_att_low, ypll_att_up) %>%
    # inner_join: solo años presentes en ambos (survey years)
    inner_join(pif_df, by = c("year", "gender", "age_group", "disease")) %>%
    mutate(
      avoid      = att      * pif,
      avoid_low  = att      * pif_low,
      avoid_up   = att      * pif_up,
      avoid_ypll     = ypll_att     * pif,
      avoid_ypll_low = ypll_att_low * pif,
      avoid_ypll_up  = ypll_att_up  * pif
    )
}

# ── PALETAS ───────────────────────────────────────────────────
# Sexo / Total — azul, rojo, gris oscuro
col_group  <- c("Men" = "#2166ac", "Women" = "#b2182b", "Overall" = "#555555")
fill_group <- c("Men"     = adjustcolor("#2166ac", 0.25),
                "Women"   = adjustcolor("#b2182b", 0.25),
                "Overall" = adjustcolor("#555555", 0.15))
lt_group   <- c("Men" = "solid", "Women" = "solid", "Overall" = "dashed")

# Grupos etarios — secuencial cálido (naranja→rojo oscuro)
pal_age <- c("#fdbe85", "#fd8d3c", "#e6550d", "#a63603")

# Escenarios — púrpura / verde (claramente distinto de sexo y grupos etarios)
col_cf  <- c("hed" = "#7b2d8b", "volume" = "#1a9850")
fill_cf <- c("hed"    = adjustcolor("#7b2d8b", 0.2),
             "volume" = adjustcolor("#1a9850", 0.2))

# Dataset (sensibilidad) — azul petróleo / rojo
col_dataset  <- c("original" = "#0571b0", "sensitivity" = "#ca0020")
fill_dataset <- c("original"    = adjustcolor("#0571b0", 0.2),
                  "sensitivity" = adjustcolor("#ca0020", 0.2))

label_cf   <- c("hed" = "HED shift", "volume" = "Volume shift")
label_shift <- c("10%" = "10 % reduction", "20%" = "20 % reduction", "30%" = "30 % reduction")


# ╔══════════════════════════════════════════════════════════╗
# ║              MAIN ANALYSIS — data_sens                  ║
# ╚══════════════════════════════════════════════════════════╝

paf_main <- paf_comp %>% filter(dataset == "sensitivity")
pif_main <- pif_comp %>% filter(dataset == "sensitivity")

att_main   <- build_att(paf_main)
avoid_main <- build_avoid(att_main, pif_main)

# ─── TABLA 1: Estadísticas descriptivas — bebedores actuales ──
wtd_med <- function(x, w) {
  ord  <- order(x); x <- x[ord]; w <- w[ord]
  cumw <- cumsum(w) / sum(w)
  x[which(cumw >= 0.5)[1]]
}

d_binge <- readRDS("PIF addiction/Data/data_binge_sensitivity.rds") |>
  filter(!is.na(hed), !is.na(volajohdia), volajohdia <= 150) |>
  mutate(
    hed_grp = factor(hed, 0:1, c("No HED", "HED")),
    sex_lab = ifelse(sexo == "Hombre", "Male", "Female"),
    age_lab = factor(edad_tramo, 1:4, c("15-29", "30-44", "45-59", ">=60"))
  )

d_binge |>
  group_by(year, hed) |>
  summarise(n = sum(exp), .groups = "drop_last") |>
  mutate(pct = n / sum(n) * 100)

fmt_cont <- function(vals)
  paste0(round(median(vals), 1), " (", round(min(vals), 1), "-", round(max(vals), 1), ")")
fmt_pct  <- function(vals)
  paste0(round(mean(vals), 1), "% (", round(min(vals), 1), "-", round(max(vals), 1), "%)")

# Consumo diario — mediana ponderada anual
cons_ann <- d_binge |>
  group_by(year, hed_grp) |>
  summarise(med = wtd_med(volajohdiams, exp), .groups = "drop")

cons_tot <- d_binge |>
  group_by(year) |>
  summarise(med = wtd_med(volajohdia, exp), .groups = "drop") |>
  mutate(hed_grp = factor("Overall"))

cons_cell <- bind_rows(cons_ann, cons_tot) |>
  group_by(hed_grp) |>
  summarise(cell = fmt_cont(med), .groups = "drop") |>
  mutate(variable = "Daily alcohol use (g)", subgroup = "")

# Sexo — % ponderado anual
sex_ann <- d_binge |>
  group_by(year, hed_grp) |>
  summarise(
    Male   = sum(exp[sex_lab == "Male"])   / sum(exp) * 100,
    Female = sum(exp[sex_lab == "Female"]) / sum(exp) * 100,
    .groups = "drop"
  )

sex_tot <- d_binge |>
  group_by(year) |>
  summarise(
    Male   = sum(exp[sex_lab == "Male"])   / sum(exp) * 100,
    Female = sum(exp[sex_lab == "Female"]) / sum(exp) * 100,
    .groups = "drop"
  ) |>
  mutate(hed_grp = factor("Overall"))

sex_cell <- bind_rows(sex_ann, sex_tot) |>
  group_by(hed_grp) |>
  summarise(Male = fmt_pct(Male), Female = fmt_pct(Female), .groups = "drop") |>
  pivot_longer(c(Male, Female), names_to = "subgroup", values_to = "cell") |>
  mutate(variable = "Sex")

# Grupo etario — % ponderado anual
age_labels <- c(a1 = "15-29", a2 = "30-44", a3 = "45-59", a4 = ">=60")

age_ann <- d_binge |>
  group_by(year, hed_grp) |>
  summarise(
    a1 = sum(exp[edad_tramo == 1]) / sum(exp) * 100,
    a2 = sum(exp[edad_tramo == 2]) / sum(exp) * 100,
    a3 = sum(exp[edad_tramo == 3]) / sum(exp) * 100,
    a4 = sum(exp[edad_tramo == 4]) / sum(exp) * 100,
    .groups = "drop"
  )

age_tot <- d_binge |>
  group_by(year) |>
  summarise(
    a1 = sum(exp[edad_tramo == 1]) / sum(exp) * 100,
    a2 = sum(exp[edad_tramo == 2]) / sum(exp) * 100,
    a3 = sum(exp[edad_tramo == 3]) / sum(exp) * 100,
    a4 = sum(exp[edad_tramo == 4]) / sum(exp) * 100,
    .groups = "drop"
  ) |>
  mutate(hed_grp = factor("Overall"))

age_cell <- bind_rows(age_ann, age_tot) |>
  group_by(hed_grp) |>
  summarise(across(c(a1, a2, a3, a4), fmt_pct), .groups = "drop") |>
  pivot_longer(c(a1, a2, a3, a4), names_to = "subgroup", values_to = "cell") |>
  mutate(subgroup = age_labels[subgroup],
         variable = "Age group")

tab1_desc <- bind_rows(cons_cell, sex_cell, age_cell) |>
  pivot_wider(names_from = hed_grp, values_from = cell) |>
  dplyr::select(variable, subgroup, HED, `No HED`, Overall)

export(tab1_desc, "PIF addiction/results/tab1_descriptive.xlsx")

# ─── TABLA S1: Muertes y YPLL atribuibles (mediana entre años) ─
tab1_att <- att_main %>%
  filter(!is.na(att)) %>%
  group_by(gender, age_group, disease) %>%
  summarise(across(c(att, att_low, att_up, ypll_att, ypll_att_low, ypll_att_up),
                   ~ median(.x, na.rm = TRUE)),
            .groups = "drop") %>%
  mutate(
    age_lab  = recode(as.character(age_group),
                      `1` = "15–29", `2` = "30–44", `3` = "45–59", `4` = "60+"),
    cell_att  = paste0(round(att),  " (", round(att_low),  "–", round(att_up),  ")"),
    cell_ypll = paste0(round(ypll_att), " (", round(ypll_att_low), "–", round(ypll_att_up), ")")
  )

export(tab1_att %>%
         pivot_wider(id_cols = c(disease, gender),
                     names_from = age_lab,
                     values_from = cell_att),
       "PIF addiction/results/tab1_att_main.xlsx")

# ─── TABLA 2: Muertes evitables — HED vs Volume, por shift ───
tab2_avoid <- avoid_main %>%
  group_by(year, disease, shift, cf_type) %>%
  summarise(avoid     = sum(avoid,     na.rm = TRUE),
            avoid_low = sum(avoid_low, na.rm = TRUE),
            avoid_up  = sum(avoid_up,  na.rm = TRUE),
            att       = sum(att,       na.rm = TRUE),
            .groups = "drop") %>%
  group_by(disease, shift, cf_type) %>%
  summarise(avoid     = round(median(avoid),     0),
            avoid_low = round(median(avoid_low), 0),
            avoid_up  = round(median(avoid_up),  0),
            att       = round(median(att),       0),
            .groups = "drop") %>%
  mutate(
    p_avoid = round(avoid / att * 100, 1),
    cell    = paste0(avoid, " (", avoid_low, "–", avoid_up, ")")
  )

export(tab2_avoid %>%
         dplyr::select(disease, cf_type, shift, cell, p_avoid) %>%
         pivot_wider(names_from = cf_type, values_from = c(cell, p_avoid)),
       "PIF addiction/results/tab2_avoid_main.xlsx")

# ─── FIGURA 1: Muertes evitables en el tiempo — HED 


# ── Build the three-layer data ─────────────────────────────────
fig1b_data <- avoid_main |>
  group_by(year, disease, cf_type) %>% 
  filter(cf_type == "hed") %>% 
  summarise(
    deaths = sum(deaths, na.rm = TRUE),
    att    = sum(att,    na.rm = TRUE),
    avoid  = sum(avoid,  na.rm = TRUE),
    .groups = "drop"
  ) |>
  pivot_longer(c(deaths, att, avoid),
               names_to  = "layer",
               values_to = "value") |>
  mutate(layer = factor(layer,
                        levels = c("deaths", "att", "avoid"),
                        labels = c("Total deaths",
                                   "Attributable deaths")))

# ── Palette ───────────────────────────────────────────────────
layer_fills <- c(
  "Total deaths"                       = "#BDBDBD",
  "Attributable deaths"                = "#6BAED6",
  "Avoidable deaths (counterfactual)"  = "#08519C"
)
layer_cols <- c(
  "Total deaths"                       = "#969696",
  "Attributable deaths"                = "#3182BD",
  "Avoidable deaths (counterfactual)"  = "#08306B"
)

# ── Plot ──────────────────────────────────────────────────────

panel_a_data <- avoid_main |>
  distinct(year, disease, gender, age_group, deaths, att) |>
  group_by(year, disease) |>
  summarise(deaths = sum(deaths, na.rm = TRUE),
            att    = sum(att,    na.rm = TRUE),
            .groups = "drop") |>
  pivot_longer(c(deaths, att),
               names_to = "layer", values_to = "value") |>
  mutate(layer = factor(layer,
                        levels = c("deaths", "att"),
                        labels = c("Total deaths", "Alcohol-attributable deaths")))

pa <- ggplot(panel_a_data,
             aes(x = year, y = value, fill = layer, color = layer)) +
  geom_area(alpha = 0.65, position = "identity") +
  facet_wrap(~ disease, nrow = 1,
             labeller = labeller(disease = label_wrap_gen(18))) +
  scale_fill_manual(
    values = c("Total deaths" = "#BDBDBD", "Alcohol-attributable deaths" = "#6BAED6"),
    name = NULL) +
  scale_color_manual(
    values = c("Total deaths" = "#969696", "Alcohol-attributable deaths" = "#3182BD"),
    name = NULL) +
  scale_x_continuous(breaks = seq(2012, 2024,2)) +
  scale_y_continuous(labels = label_number(big.mark = ".", decimal.mark = ",")) +
  labs(x = NULL, y = "Number of deaths") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position   = "bottom",
    axis.text.x       = element_text(angle = 45, hjust = 1),
    strip.background  = element_rect(fill = "grey95", color = NA),
    panel.grid.minor  = element_blank()
  )

ggsave("PIF addiction/results/fig1_att_total.tiff",
       plot = pa, device = "tiff", dpi = 300,
       width = 11, height = 8, units = "in", compression = "lzw")

# ─── FIGURA 2: Muertes evitables ─────────────

# líneas × escenario, coloreadas por cf_type, trazo por shift
# Paleta secuencial de azules — coherente con el azul de panel A
pal_shift <- c(
  "10 % reduction" = "#6BAED6",   # azul claro (igual al att en panel A)
  "20 % reduction" = "#2171B5",   # azul medio
  "30 % reduction" = "#08306B"    # azul oscuro
)


pb2_summary_sex <- avoid_main |>
  filter(cf_type == "hed") |>
  group_by(year, disease, gender, shift) |>
  summarise(
    avoid     = sum(avoid,     na.rm = TRUE),
    avoid_low = sum(avoid_low, na.rm = TRUE),
    avoid_up  = sum(avoid_up,  na.rm = TRUE),
    .groups   = "drop"
  ) |>
  mutate(shift = factor(shift, levels = c("10%","20%","30%"), labels = label_shift)) |>
  group_by(disease, gender, shift) |>
  summarise(
    avoid_med = median(avoid),
    low_med   = median(avoid_low),
    up_med    = median(avoid_up),
    .groups   = "drop"
  )

pb2_summary_sex_supp <- avoid_main |>
  filter(cf_type == "volume") |>
  group_by(year, disease, gender, shift) |>
  summarise(
    avoid     = sum(avoid,     na.rm = TRUE),
    avoid_low = sum(avoid_low, na.rm = TRUE),
    avoid_up  = sum(avoid_up,  na.rm = TRUE),
    .groups   = "drop"
  ) |>
  mutate(shift = factor(shift, levels = c("10%","20%","30%"), labels = label_shift)) |>
  group_by(disease, gender, shift) |>
  summarise(
    avoid_med = median(avoid),
    low_med   = median(avoid_low),
    up_med    = median(avoid_up),
    .groups   = "drop"
  )

pb2_main_sex <- ggplot(pb2_summary_sex,
                            aes(x = avoid_med, y = disease,
                                color = shift, shape = shift)) +
  geom_linerange(aes(xmin = low_med, xmax = up_med),
                 linewidth = 0.8, alpha = 0.7,
                 position = position_dodge(width = 0.6)) +
  geom_point(size = 3.5,
             position = position_dodge(width = 0.6)) +
  facet_wrap(~ gender) +
  scale_color_manual(values = pal_shift, name = "Reduction scenario") +
  scale_shape_manual(
    values = c("10 % reduction" = 16,
               "20 % reduction" = 17,
               "30 % reduction" = 15),
    name = "Reduction scenario") +
  scale_x_continuous(
    labels = label_number(big.mark = ".", decimal.mark = ","),
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(x = "Alcohol-attributable injury deaths avoided (median 2012–2024)",
       y = NULL) +
  theme_minimal(base_size = 13) +
  theme(legend.position    = "bottom",
        panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        strip.background   = element_rect(fill = "grey90", color = NA))+geom_vline(
          xintercept = 0,
          color = "grey40",
          linewidth = 0.8
        )

pb2_main_sex_supp <- ggplot(pb2_summary_sex_supp,
       aes(x = avoid_med, y = disease,
           color = shift, shape = shift)) +
  geom_linerange(aes(xmin = low_med, xmax = up_med),
                 linewidth = 0.8, alpha = 0.7,
                 position = position_dodge(width = 0.6)) +
  geom_point(size = 3.5,
             position = position_dodge(width = 0.6)) +
  facet_wrap(~ gender) +
  scale_color_manual(values = pal_shift, name = "Reduction scenario") +
  scale_shape_manual(
    values = c("10 % reduction" = 16,
               "20 % reduction" = 17,
               "30 % reduction" = 15),
    name = "Reduction scenario") +
  scale_x_continuous(
    labels = label_number(big.mark = ".", decimal.mark = ","),
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(x = "Alcohol-attributable injury deaths avoided (median 2012–2024)",
       y = NULL) +
  theme_minimal(base_size = 13) +
  theme(legend.position    = "bottom",
        panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        strip.background   = element_rect(fill = "grey90", color = NA))+geom_vline(
          xintercept = 0,
          color = "grey40",
          linewidth = 0.8
        )


ggsave("PIF addiction/results/fig2_att_avoided_MAIN_SEX.tiff",
       plot = pb2_main_sex, device = "tiff", dpi = 300,
       width = 11, height = 8, units = "in", compression = "lzw")

ggsave("PIF addiction/results/fig2_att_avoided_SUPP_SEX.tiff",
       plot = pb2_main_sex_supp, device = "tiff", dpi = 300,
       width = 11, height = 8, units = "in", compression = "lzw")

# ─── FIGURA 3: Muertes evitables por grupo etario ─────────────
fig3_data <- avoid_main %>%
  group_by(year, disease, age_group, shift, cf_type,gender) %>%
  summarise(avoid     = sum(avoid,     na.rm = TRUE),
            avoid_low = sum(avoid_low, na.rm = TRUE),
            avoid_up  = sum(avoid_up,  na.rm = TRUE),
            .groups = "drop") %>%
  mutate(age_group = factor(age_group, 1:4, c("15–29","30–44","45–59","60+")),
         shift     = factor(shift,   levels = c("10%","20%","30%"), labels = label_shift),
         cf_type   = factor(cf_type, levels = c("hed","volume"),    labels = label_cf))

fig3 <- ggplot(filter(fig3_data, cf_type == "HED shift", gender == "Men"),
               aes(x = year, y = avoid, color = age_group, group = age_group)) +
  geom_hline(yintercept = 0, color = "grey35", linewidth = 0.6) +
  geom_errorbar(aes(ymin = avoid_low, ymax = avoid_up),
                width = 0.3, linewidth = 0.4, alpha = 0.6) +
  geom_line(linewidth = 1) +
  geom_point(size = 2, shape = 21, fill = "white") +
  facet_grid(rows = vars(shift), cols = vars(disease)) +
  scale_color_manual(values = pal_age, name = "Age group") +
  scale_x_continuous(breaks = seq(2012, 2024, 2)) +
  scale_y_continuous(labels = label_number(big.mark = ".", decimal.mark = ",")) +
  labs(x = "", y = "Deaths avoided") +
  theme_minimal(base_size = 13) +
  theme(legend.position   = "bottom",
        axis.text.x       = element_text(angle = 45, hjust = 1),
        strip.background  = element_rect(fill = "grey95", color = NA),
        panel.grid.minor  = element_blank())

fig3

fig3_hed_men <- fig3 %+% filter(fig3_data, cf_type == "HED shift", gender == "Men")
fig3_hed_women <- fig3 %+% filter(fig3_data, cf_type == "HED shift", gender == "Women")

fig3_main <- fig3_hed_men/fig3_hed_women

fig3_vol_men <- fig3 %+% filter(fig3_data, cf_type == "Volume shift", gender == "Men")
fig3_vol_women <- fig3 %+% filter(fig3_data, cf_type == "Volume shift", gender == "Women")

fig3_supp <- fig3_vol_men/fig3_vol_women

#---------------#
# MAIN ANALYSIS #
#---------------#

pal_shift <- c(
  "10 % reduction" = "#6BAED6",
  "20 % reduction" = "#2171B5",
  "30 % reduction" = "#08306B"
)
fig3_summary <- avoid_main |>
  filter(cf_type == "hed") |>
  group_by(year, disease, gender, age_group, shift) |>
  summarise(
    avoid     = sum(avoid,     na.rm = TRUE),
    avoid_low = sum(avoid_low, na.rm = TRUE),
    avoid_up  = sum(avoid_up,  na.rm = TRUE),
    .groups   = "drop"
  ) |>
  mutate(
    age_group = factor(age_group, 1:4, c("15–29","30–44","45–59","60+")),
    shift     = factor(shift, levels = c("10%","20%","30%"), labels = label_shift)
  ) |>
  group_by(disease, gender, age_group, shift) |>
  summarise(
    avoid_med = median(avoid),
    low_med   = median(avoid_low),
    up_med    = median(avoid_up),
    .groups   = "drop"
  )

fig3_summary_supp <- avoid_main |>
  filter(cf_type == "volume") |>
  group_by(year, disease, gender, age_group, shift) |>
  summarise(
    avoid     = sum(avoid,     na.rm = TRUE),
    avoid_low = sum(avoid_low, na.rm = TRUE),
    avoid_up  = sum(avoid_up,  na.rm = TRUE),
    .groups   = "drop"
  ) |>
  mutate(
    age_group = factor(age_group, 1:4, c("15–29","30–44","45–59","60+")),
    shift     = factor(shift, levels = c("10%","20%","30%"), labels = label_shift)
  ) |>
  group_by(disease, gender, age_group, shift) |>
  summarise(
    avoid_med = median(avoid),
    low_med   = median(avoid_low),
    up_med    = median(avoid_up),
    .groups   = "drop"
  )

fig3_main_forest <- ggplot(fig3_summary,
       aes(x = avoid_med, y = age_group,
           color = shift, shape = shift)) +
  geom_linerange(aes(xmin = low_med, xmax = up_med),
                 linewidth = 0.8, alpha = 0.7,
                 position = position_dodge(width = 0.6)) +
  geom_point(size = 3.5,
             position = position_dodge(width = 0.6)) +
  facet_grid(rows = vars(disease), cols = vars(gender),
             labeller = labeller(disease = label_wrap_gen(18))) +
  scale_color_manual(values = pal_shift, name = "Reduction scenario") +
  scale_shape_manual(
    values = c("10 % reduction" = 16,
               "20 % reduction" = 17,
               "30 % reduction" = 15),
    name = "Reduction scenario") +
  scale_x_continuous(
    labels = label_number(big.mark = ".", decimal.mark = ","),
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(x = "Alcohol-attributable injury deaths avoided (median 2012–2024)",
       y = NULL) +
  theme_minimal(base_size = 13) +
  theme(legend.position    = "bottom",
        panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        strip.background   = element_rect(fill = "grey95", color = NA))+geom_vline(
          xintercept = 0,
          color = "grey40",
          linewidth = 0.8
        )

fig3_main_forest_supp <- ggplot(fig3_summary_supp,
                           aes(x = avoid_med, y = age_group,
                               color = shift, shape = shift)) +
  geom_linerange(aes(xmin = low_med, xmax = up_med),
                 linewidth = 0.8, alpha = 0.7,
                 position = position_dodge(width = 0.6)) +
  geom_point(size = 3.5,
             position = position_dodge(width = 0.6)) +
  facet_grid(rows = vars(disease), cols = vars(gender),
             labeller = labeller(disease = label_wrap_gen(18))) +
  scale_color_manual(values = pal_shift, name = "Reduction scenario") +
  scale_shape_manual(
    values = c("10 % reduction" = 16,
               "20 % reduction" = 17,
               "30 % reduction" = 15),
    name = "Reduction scenario") +
  scale_x_continuous(
    labels = label_number(big.mark = ".", decimal.mark = ","),
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(x = "Alcohol-attributable injury deaths avoided (median 2012–2024)",
       y = NULL) +
  theme_minimal(base_size = 13) +
  theme(legend.position    = "bottom",
        panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        strip.background   = element_rect(fill = "grey95", color = NA))+geom_vline(
          xintercept = 0,
          color = "grey40",
          linewidth = 0.8
        )

ggsave("PIF addiction/results/fig3_avoidxage.tiff",
       plot = fig3_main, device = "tiff", dpi = 300,
       width = 11, height = 8, units = "in", compression = "lzw")

ggsave("PIF addiction/results/fig3_avoidxageSUPP.tiff",
       plot = fig3_supp, device = "tiff", dpi = 300,
       width = 11, height = 8, units = "in", compression = "lzw")


ggsave("PIF addiction/results/fig3_avoidxage_MAINFOREST.tiff",
       plot = fig3_main_forest, device = "tiff", dpi = 300,
       width = 11, height = 8, units = "in", compression = "lzw")

ggsave("PIF addiction/results/fig3_avoidxage_SUPPFOREST.tiff",
       plot = fig3_main_forest_supp, device = "tiff", dpi = 300,
       width = 11, height = 8, units = "in", compression = "lzw")

# ─── FIGURA 4: Total and attributable YPLL ─────
pa_ypll_data <- avoid_main |>
  filter(cf_type == "hed") |>
  distinct(year, disease, gender, age_group, ypll, ypll_att) |>
  group_by(year, disease, gender) |>
  summarise(
    ypll     = sum(ypll,     na.rm = TRUE),
    ypll_att = sum(ypll_att, na.rm = TRUE),
    .groups  = "drop"
  ) |>
  pivot_longer(c(ypll, ypll_att),
               names_to = "layer", values_to = "value") |>
  mutate(layer = factor(layer,
                        levels = c("ypll", "ypll_att"),
                        labels = c("Total YPLL", "Attributable YPLL")))

pa_ypll <- ggplot(pa_ypll_data,
                  aes(x = year, y = value, fill = layer, color = layer)) +
  geom_area(alpha = 0.65, position = "identity") +
  facet_grid(rows = vars(disease), cols = vars(gender),
             labeller = labeller(disease = label_wrap_gen(18))) +
  scale_fill_manual(
    values = c("Total YPLL" = "#BDBDBD", "Attributable YPLL" = "#6BAED6"),
    name = NULL) +
  scale_color_manual(
    values = c("Total YPLL" = "#969696", "Attributable YPLL" = "#3182BD"),
    name = NULL) +
  scale_x_continuous(breaks = seq(2012, 2024, 2)) +
  scale_y_continuous(labels = label_number(big.mark = ".", decimal.mark = ",")) +
  labs(x = NULL, y = "YPLL", tag = "A") +
  theme_minimal(base_size = 12) +
  theme(legend.position   = "bottom",
        axis.text.x       = element_text(angle = 45, hjust = 1),
        strip.background  = element_rect(fill = "grey95", color = NA),
        panel.grid.minor  = element_blank())

ggsave("PIF addiction/results/fig4_att_ypll.tiff",
       plot = pa_ypll, device = "tiff", dpi = 300,
       width = 11, height = 8, units = "in", compression = "lzw")

# ─── FIGURA 5: % de muertes atribuibles que son evitables ─────
pal_shift <- c(
  "10 % reduction" = "#6BAED6",
  "20 % reduction" = "#2171B5",
  "30 % reduction" = "#08306B"
)

fig_ypll_data <- avoid_main |>
  filter(cf_type == "hed") |>
  group_by(year, disease, gender, shift) |>
  summarise(
    ypll         = sum(ypll,         na.rm = TRUE),
    avoid_ypll   = sum(avoid_ypll,   na.rm = TRUE),
    avoid_ypll_low = sum(avoid_ypll_low, na.rm = TRUE),
    avoid_ypll_up  = sum(avoid_ypll_up,  na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    pct     = avoid_ypll     / ypll * 100,
    pct_low = avoid_ypll_low / ypll * 100,
    pct_up  = avoid_ypll_up  / ypll * 100,
    shift   = factor(shift, levels = c("10%","20%","30%"), labels = label_shift)
  )

fig_ypll_data_supp <- avoid_main |>
  filter(cf_type == "volume") |>
  group_by(year, disease, gender, shift) |>
  summarise(
    ypll         = sum(ypll,         na.rm = TRUE),
    avoid_ypll   = sum(avoid_ypll,   na.rm = TRUE),
    avoid_ypll_low = sum(avoid_ypll_low, na.rm = TRUE),
    avoid_ypll_up  = sum(avoid_ypll_up,  na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    pct     = avoid_ypll     / ypll * 100,
    pct_low = avoid_ypll_low / ypll * 100,
    pct_up  = avoid_ypll_up  / ypll * 100,
    shift   = factor(shift, levels = c("10%","20%","30%"), labels = label_shift)
  )

fig5_summary <- fig_ypll_data |>
  group_by(disease, gender, shift) |>
  summarise(
    pct_med = median(pct),
    low_med = median(pct_low),
    up_med  = median(pct_up),
    .groups = "drop"
  )

fig5_summary_supp <- fig_ypll_data_supp |>
  group_by(disease, gender, shift) |>
  summarise(
    pct_med = median(pct),
    low_med = median(pct_low),
    up_med  = median(pct_up),
    .groups = "drop"
  )

fig5_main_forest <- ggplot(fig5_summary,
                           aes(x = pct_med, y = disease,
                               color = shift, shape = shift)) +
  geom_linerange(aes(xmin = low_med, xmax = up_med),
                 linewidth = 0.8, alpha = 0.7,
                 position = position_dodge(width = 0.6)) +
  geom_point(size = 3.5,
             position = position_dodge(width = 0.6)) +
  facet_wrap(~ gender) +
  scale_color_manual(values = pal_shift, name = "Reduction scenario") +
  scale_shape_manual(
    values = c("10 % reduction" = 16,
               "20 % reduction" = 17,
               "30 % reduction" = 15),
    name = "Reduction scenario") +
  scale_x_continuous(
    labels = function(x) paste0(x, "%"),
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(x = "% of total YPLL avoided (median 2012–2024)",
       y = NULL) +
  theme_minimal(base_size = 13) +
  theme(legend.position    = "bottom",
        panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        strip.background   = element_rect(fill = "grey95", color = NA))+geom_vline(
          xintercept = 0,
          color = "grey40",
          linewidth = 0.8
        )

fig5_main_forest_supp <- ggplot(fig5_summary_supp,
                           aes(x = pct_med, y = disease,
                               color = shift, shape = shift)) +
  geom_linerange(aes(xmin = low_med, xmax = up_med),
                 linewidth = 0.8, alpha = 0.7,
                 position = position_dodge(width = 0.6)) +
  geom_point(size = 3.5,
             position = position_dodge(width = 0.6)) +
  facet_wrap(~ gender) +
  scale_color_manual(values = pal_shift, name = "Reduction scenario") +
  scale_shape_manual(
    values = c("10 % reduction" = 16,
               "20 % reduction" = 17,
               "30 % reduction" = 15),
    name = "Reduction scenario") +
  scale_x_continuous(
    labels = function(x) paste0(x, "%"),
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(x = "% of total YPLL avoided (median 2012–2024)",
       y = NULL) +
  theme_minimal(base_size = 13) +
  theme(legend.position    = "bottom",
        panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        strip.background   = element_rect(fill = "grey95", color = NA))+geom_vline(
          xintercept = 0,
          color = "grey40",
          linewidth = 0.8
        )




fig5_main <- ggplot(fig_ypll_data,
       aes(x = year, y = pct, color = shift, fill = shift, group = shift)) +
  geom_errorbar(aes(ymin = pct_low, ymax = pct_up),
                width = 0.3, linewidth = 0.4, alpha = 0.6) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2, shape = 21, fill = "white") +
  facet_grid(rows = vars(disease), cols = vars(gender),
             labeller = labeller(disease = label_wrap_gen(18))) +
  scale_color_manual(values = pal_shift, name = "Reduction scenario") +
  scale_fill_manual( values = pal_shift, name = "Reduction scenario") +
  scale_x_continuous(breaks = seq(2012, 2024, 2)) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(x = "Year",
       y = "% of total YPLL avoided") +
  theme_minimal(base_size = 13) +
  theme(legend.position   = "bottom",
        axis.text.x       = element_text(angle = 45, hjust = 1),
        strip.background  = element_rect(fill = "grey95", color = NA),
        panel.grid.minor  = element_blank())

ggsave("PIF addiction/results/fig5_ypll_avoid.tiff",
       plot = fig5_main, device = "tiff", dpi = 300,
       width = 11, height = 8, units = "in", compression = "lzw")

ggsave("PIF addiction/results/fig5_ypll_avoid_MAIN.tiff",
       plot = fig5_main_forest, device = "tiff", dpi = 300,
       width = 11, height = 8, units = "in", compression = "lzw")


ggsave("PIF addiction/results/fig5_ypll_avoid_SUPP.tiff",
       plot = fig5_main_forest_supp, device = "tiff", dpi = 300,
       width = 11, height = 8, units = "in", compression = "lzw")

# ─── ESTADÍSTICAS IN-TEXT ──────────────────────────────────────
# 1. Consumo promedio ponderado por sexo y HED/NHED (en tragos/día — unidad de volajohdia)
cat("── 1. Consumo de alcohol (tragos/día de consumo) por HED/NHED y sexo ──\n")
d_binge |>
  filter(volajohdia<=150) %>% 
  group_by(year, sex_lab, hed_grp) |>
  summarise(wmean = weighted.mean(volajohdia, exp, na.rm = TRUE), .groups = "drop") |>
  group_by(sex_lab, hed_grp) |>
  summarise(
    mean_drinks = round(mean(wmean), 1),
    .groups = "drop"
  ) |>
  pivot_wider(names_from = hed_grp, values_from = mean_drinks)



library(Hmisc)
wtd_med <- function(x, w) {
  wtd.quantile(x, weights = w, probs = 0.5, na.rm = TRUE)
}

# ── 1. Consumo (igual que tab1_desc) ──────────────────────────
cat("── 1. Consumo diario de alcohol (g/día) ──\n")
cons_ann |>
  group_by(hed_grp) |>
  summarise(cell = fmt_cont(med), .groups = "drop")

# ── 2. Prevalencia HED por sexo (% de bebedores que son HED) ──
cat("\n── 2. Prevalencia HED entre bebedores ──\n")

d_binge %>% 
  group_by(sexo, year, edad_tramo) %>% 
  count(hed) %>% 
  mutate(pct = n/sum(n)*100) %>% 
  filter(sexo == "Mujer")

d_binge |>
  group_by(year, sex_lab) |>
  summarise(prev_hed = sum(exp[hed == 1]) / sum(exp) * 100,
            .groups  = "drop") |>
  group_by(sex_lab) |>
  summarise(cell = fmt_pct(prev_hed), .groups = "drop") |>
  bind_rows(
    d_binge |>
      group_by(year) |>
      summarise(prev_hed = sum(exp[hed == 1]) / sum(exp) * 100,
                .groups  = "drop") |>
      summarise(cell = fmt_pct(prev_hed)) |>
      mutate(sex_lab = "Total")
  )
# ── 3. Total muertes por lesiones ─────────────────────────────
cat("── 3. Muertes por lesiones (promedio anual 2012–2024) ──\n")

deaths_total <- att_main |>
  distinct(year, gender, age_group, disease, deaths) |>
  group_by(year) |>
  summarise(total = sum(deaths), .groups = "drop")

cat("Total (promedio anual):", round(mean(deaths_total$total)), "\n")
cat("Rango:", min(deaths_total$total), "–", max(deaths_total$total), "\n\n")

deaths_cause <- att_main |>
  distinct(year, gender, age_group, disease, deaths) |>
  group_by(year, disease) |>
  summarise(n = sum(deaths), .groups = "drop") |>
  group_by(disease) |>
  summarise(mean_n = round(mean(n)), .groups = "drop") |>
  mutate(pct = round(mean_n / sum(mean_n) * 100, 1))

deaths_cause

# ── 4. Muertes atribuibles ────────────────────────────────────
cat("── 4. Muertes atribuibles (promedio anual, mediana entre años) ──\n")

att_year <- att_main |>
  filter(!is.na(att)) |>
  group_by(year) |>
  summarise(
    att    = sum(att,    na.rm = TRUE),
    deaths = sum(deaths, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(pct = att / deaths * 100)

cat("Mediana anual muertes atribuibles:", round(median(att_year$att)), "\n")
cat("Rango:", round(min(att_year$att)), "–", round(max(att_year$att)), "\n")
cat("% del total de muertes por lesiones (mediana):", round(median(att_year$pct), 1), "%\n\n")

# % por sexo
cat("── % muertes atribuibles por sexo ──\n")
att_main |>
  filter(!is.na(att)) |>
  group_by(year, gender) |>
  summarise(att = sum(att, na.rm = TRUE), .groups = "drop") |>
  group_by(year) |>
  mutate(pct = att / sum(att) * 100) |>
  group_by(gender) |>
  summarise(pct = round(median(pct), 1), .groups = "drop")




# ── 5 & 6. % muertes atribuibles y YPLL atribuible por grupo etario ──
cat("── 5. % muertes atribuibles por grupo etario ──\n")
att_main |>
  filter(!is.na(att)) |>
  mutate(age_lab = factor(age_group, 1:4, c("15–29","30–44","45–59","60+"))) |>
  group_by(year, age_lab) |>
  summarise(att = sum(att, na.rm = TRUE), .groups = "drop") |>
  group_by(year) |>
  mutate(pct = att / sum(att) * 100) |>
  group_by(age_lab) |>
  summarise(pct = round(median(pct), 1), .groups = "drop")

cat("\n── 6. % YPLL atribuible por grupo etario ──\n")
att_main |>
  filter(!is.na(ypll_att)) |>
  mutate(age_lab = factor(age_group, 1:4, c("15–29","30–44","45–59","60+"))) |>
  group_by(year, age_lab) |>
  summarise(ypll_att = sum(ypll_att, na.rm = TRUE), .groups = "drop") |>
  group_by(year) |>
  mutate(pct = ypll_att / sum(ypll_att) * 100) |>
  group_by(age_lab) |>
  summarise(pct = round(median(pct), 1), .groups = "drop")


# % muertes evitables promedio por cf_type y shift
avoid_main %>%
  group_by(year, shift, cf_type) %>%
  summarise(avoid = sum(avoid), att = sum(att), .groups = "drop") %>%
  mutate(p = avoid / att * 100) %>%
  group_by(shift, cf_type) %>%
  summarise(mean_p    = round(mean(p), 1),
            mean_p_lo = round(mean(avoid_main %>%
                                     filter(shift == first(shift), cf_type == first(cf_type)) %>%
                                     group_by(year) %>%
                                     summarise(p = sum(avoid_low)/sum(att)*100) %>%
                                     pull(p)), 1),
            .groups = "drop")

# Versión simple y robusta:
avoid_main %>%
  group_by(year, shift, cf_type) %>%
  summarise(avoid     = sum(avoid,     na.rm = TRUE),
            avoid_low = sum(avoid_low, na.rm = TRUE),
            avoid_up  = sum(avoid_up,  na.rm = TRUE),
            att       = sum(att,       na.rm = TRUE),
            .groups = "drop") %>%
  mutate(p     = round(avoid     / att * 100, 1),
         p_low = round(avoid_low / att * 100, 1),
         p_up  = round(avoid_up  / att * 100, 1)) %>%
  group_by(shift, cf_type) %>%
  summarise(mean_p     = round(mean(p),     1),
            mean_p_low = round(mean(p_low), 1),
            mean_p_up  = round(mean(p_up),  1),
            .groups = "drop")


# ╔══════════════════════════════════════════════════════════╗
# ║               SUPPLEMENTARY TABLES                       ║
# ╚══════════════════════════════════════════════════════════╝
# ── Extend pop to 2024 using 2022 as proxy ─────────────────────
pop_ext <- pop_ag |>
  filter(year == 2022) |>
  mutate(year = 2024) |>
  bind_rows(pop_ag)

age_labels <- c(`1` = "15–29", `2` = "30–44", `3` = "45–59", `4` = "60+")

# ══════════════════════════════════════════════════════════════
# TABLE 1: Median annual alcohol-attributable deaths
#          by disease × sex × age group
# ══════════════════════════════════════════════════════════════
tab1 <- att_main |>
  filter(!is.na(att)) |>
  group_by(gender, age_group, disease) |>
  summarise(
    med     = round(median(att),     1),
    low     = round(median(att_low), 1),
    up      = round(median(att_up),  1),
    .groups = "drop"
  ) |>
  mutate(
    age_lab = age_labels[as.character(age_group)],
    cell    = paste0(med, " (", low, "–", up, ")")
  ) |>
  select(disease, gender, age_lab, cell) |>
  pivot_wider(names_from = age_lab, values_from = cell)

# ══════════════════════════════════════════════════════════════
# TABLE 2: Median YPLL rate (per 100,000) with 95% CI
#          by sex × disease × age group
# ══════════════════════════════════════════════════════════════
tab2 <- att_main |>
  filter(!is.na(ypll_att)) |>
  left_join(pop_ext, by = c("year", "gender", "age_group")) |>
  mutate(
    rate_att     = ypll_att     / pop * 1e5,
    rate_att_low = ypll_att_low / pop * 1e5,
    rate_att_up  = ypll_att_up  / pop * 1e5
  ) |>
  group_by(gender, age_group, disease) |>
  summarise(
    med = round(median(rate_att)),
    low = round(median(rate_att_low)),
    up  = round(median(rate_att_up)),
    .groups = "drop"
  ) |>
  mutate(
    age_lab = age_labels[as.character(age_group)],
    cell    = paste0(med, " (", low, "–", up, ")")
  ) |>
  select(disease, gender, age_lab, cell) |>
  pivot_wider(names_from = age_lab, values_from = cell)

tab2


# TAB 3
# Helper: median avoided deaths table (by cf_type)
make_avoid_tab <- function(cf) {
  avoid_main |>
    filter(cf_type == cf) |>
    group_by(disease, gender, age_group, shift) |>
    summarise(
      med = round(median(avoid)),
      low = round(median(avoid_low)),
      up  = round(median(avoid_up)),
      .groups = "drop"
    ) |>
    mutate(
      age_lab   = age_labels[as.character(age_group)],
      shift_lab = recode(shift, "10%" = "10% reduction",
                                "20%" = "20% reduction",
                                "30%" = "30% reduction"),
      cell = paste0(med, " (", low, "–", up, ")")
    ) |>
    select(disease, gender, shift_lab, age_lab, cell) |>
    pivot_wider(names_from = age_lab, values_from = cell) |>
    arrange(disease, gender, shift_lab)
}

# Helper: median YPLL rate avoided (per 100,000) table (by cf_type)
make_ypll_avoid_tab <- function(cf) {
  avoid_main |>
    filter(cf_type == cf) |>
    left_join(pop_ext, by = c("year", "gender", "age_group")) |>
    mutate(
      rate     = avoid_ypll     / pop * 1e5,
      rate_low = avoid_ypll_low / pop * 1e5,
      rate_up  = avoid_ypll_up  / pop * 1e5
    ) |>
    group_by(disease, gender, age_group, shift) |>
    summarise(
      med = round(median(rate)),
      low = round(median(rate_low)),
      up  = round(median(rate_up)),
      .groups = "drop"
    ) |>
    mutate(
      age_lab   = age_labels[as.character(age_group)],
      shift_lab = recode(shift, "10%" = "10% reduction",
                                "20%" = "20% reduction",
                                "30%" = "30% reduction"),
      cell = paste0(med, " (", low, "–", up, ")")
    ) |>
    select(disease, gender, shift_lab, age_lab, cell) |>
    pivot_wider(names_from = age_lab, values_from = cell) |>
    arrange(disease, gender, shift_lab)
}

tab3 <- make_avoid_tab("hed")

tab4 <- make_avoid_tab("volume")
tab5 <- make_ypll_avoid_tab("hed")
tab6 <- make_ypll_avoid_tab("volume")

cat("Table 3 (avoided deaths – HED):\n"); print(tab3)
cat("\nTable 4 (avoided deaths – Volume):\n"); print(tab4)
# ╔══════════════════════════════════════════════════════════╗
# ║          SENSITIVITY ANALYSIS — og vs sens              ║
# ╚══════════════════════════════════════════════════════════╝

# ─── Muertes atribuibles: og vs sens ──────────────────────────
att_og   <- build_att(paf_comp %>% filter(dataset == "original"))
att_sens <- att_main  # ya construido arriba

att_comp <- bind_rows(
  att_og   %>% mutate(dataset = "original"),
  att_sens %>% mutate(dataset = "sensitivity")
)

# Diferencias en PAF por año/enfermedad/sexo/grupo etario
paf_diff_tab <- paf_comp %>%
  pivot_wider(id_cols = c(year, gender, age_group, disease),
              names_from = dataset, values_from = point) %>%
  mutate(diff_abs  = original - sensitivity,
         diff_rel  = (original - sensitivity) / sensitivity * 100) %>%
  group_by(disease, gender) %>%
  summarise(mean_diff_abs = round(mean(diff_abs, na.rm = TRUE), 3),
            mean_diff_rel = round(mean(diff_rel, na.rm = TRUE), 1),
            .groups = "drop")

export(paf_diff_tab, "PIF addiction/results/tab_sens_paf_diff.xlsx")

# ─── FIGURA S1: PAF og vs sens en el tiempo ───────────────────
fig_s1_data <- paf_comp %>%
  group_by(year, gender, disease, dataset) %>%
  summarise(point = mean(point, na.rm = TRUE),
            lower = mean(lower, na.rm = TRUE),
            upper = mean(upper, na.rm = TRUE),
            .groups = "drop")

fig_s1 <- ggplot(fig_s1_data, aes(x = year, y = point,
                                   color = dataset, fill = dataset,
                                   group = dataset)) +
  geom_hline(yintercept = 0, color = "grey35", linewidth = 0.6) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.8) +
  facet_grid(rows = vars(gender), cols = vars(disease)) +
  scale_color_manual(values = col_dataset,  name = "Dataset") +
  scale_fill_manual( values = fill_dataset, name = "Dataset") +
  scale_x_continuous(breaks = pretty_breaks(4)) +
  scale_y_continuous(labels = label_number(accuracy = 0.001)) +
  labs(x = "Year", y = "PAF (point estimate)") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(fill = "grey95", color = NA),
        panel.grid.minor = element_blank())

ggsave("PIF addiction/results/figS1_paf_og_vs_sens.tiff",
       plot = fig_s1, device = "tiff", dpi = 300,
       width = 11, height = 7, units = "in", compression = "lzw")

# ─── Muertes evitables: og vs sens (años en común) ────────────
avoid_comp <- bind_rows(
  build_avoid(att_og,   pif_comp %>% filter(dataset == "original"))   %>% mutate(dataset = "original"),
  build_avoid(att_sens, pif_comp %>% filter(dataset == "sensitivity")) %>% mutate(dataset = "sensitivity")
)

# Resumen: diferencia en muertes evitables en años solapados (2012–2022)
avoid_diff_tab <- avoid_comp %>%
  filter(year %in% 2012:2022) %>%
  group_by(year, disease, shift, cf_type, dataset) %>%
  summarise(avoid = sum(avoid, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = dataset, values_from = avoid) %>%
  mutate(diff_abs = original - sensitivity,
         diff_rel = (original - sensitivity) / sensitivity * 100) %>%
  group_by(disease, shift, cf_type) %>%
  summarise(mean_diff_abs = round(mean(diff_abs, na.rm = TRUE), 1),
            mean_diff_rel = round(mean(diff_rel, na.rm = TRUE), 1),
            .groups = "drop")

export(avoid_diff_tab, "PIF addiction/results/tab_sens_avoid_diff.xlsx")

# ─── FIGURA S2: Muertes evitables og vs sens ──────────────────
fig_s2_data <- avoid_comp %>%
  filter(year %in% 2012:2022) %>%
  group_by(year, disease, shift, cf_type, dataset) %>%
  summarise(avoid     = sum(avoid,     na.rm = TRUE),
            avoid_low = sum(avoid_low, na.rm = TRUE),
            avoid_up  = sum(avoid_up,  na.rm = TRUE),
            .groups = "drop") %>%
  mutate(shift   = factor(shift,   levels = c("10%","20%","30%"), labels = label_shift),
         cf_type = factor(cf_type, levels = c("hed","volume"),    labels = label_cf))

fig_s2 <- ggplot(fig_s2_data, aes(x = year, y = avoid,
                                   color = dataset, fill = dataset,
                                   linetype = dataset, group = dataset)) +
  geom_hline(yintercept = 0, color = "grey35", linewidth = 0.6) +
  geom_ribbon(aes(ymin = avoid_low, ymax = avoid_up), alpha = 0.2, color = NA) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.8) +
  facet_grid(rows = vars(shift), cols = vars(cf_type ~ disease)) +
  scale_color_manual(values = col_dataset,  name = "Dataset") +
  scale_fill_manual( values = fill_dataset, name = "Dataset") +
  scale_linetype_manual(values = c("original" = "dashed", "sensitivity" = "solid"),
                        name = "Dataset") +
  scale_x_continuous(breaks = pretty_breaks(4)) +
  scale_y_continuous(labels = label_number(big.mark = ".", decimal.mark = ",")) +
  labs(x = "Year", y = "Avoidable deaths") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(fill = "grey95", color = NA),
        panel.grid.minor = element_blank())

ggsave("PIF addiction/results/figS2_avoid_og_vs_sens.tiff",
       plot = fig_s2, device = "tiff", dpi = 300,
       width = 14, height = 9, units = "in", compression = "lzw")

# ─── TABLA RESUMEN sensitivity ─────────────────────────────────
tab_sens_summary <- avoid_comp %>%
  filter(year %in% 2012:2022) %>%
  group_by(year, disease, shift, cf_type, dataset) %>%
  summarise(avoid     = sum(avoid,     na.rm = TRUE),
            avoid_low = sum(avoid_low, na.rm = TRUE),
            avoid_up  = sum(avoid_up,  na.rm = TRUE),
            att       = sum(att,       na.rm = TRUE),
            .groups = "drop") %>%
  group_by(disease, shift, cf_type, dataset) %>%
  summarise(avoid     = round(median(avoid),     0),
            avoid_low = round(median(avoid_low), 0),
            avoid_up  = round(median(avoid_up),  0),
            att       = round(median(att),       0),
            .groups = "drop") %>%
  mutate(p_avoid = round(avoid / att * 100, 1),
         cell    = paste0(avoid, " (", avoid_low, "–", avoid_up, ")")) %>%
  pivot_wider(id_cols = c(disease, shift, cf_type),
              names_from = dataset,
              values_from = c(cell, p_avoid))

export(tab_sens_summary, "PIF addiction/results/tab_sens_summary.xlsx")

message(">>> PIF results.R completado. Outputs en PIF addiction/results/")
