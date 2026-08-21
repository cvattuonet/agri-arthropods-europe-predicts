#---------------------------------------------------------------------------------------------------------------------------------------#
# CODE DETAILS                                                                                                                          #
#                                                                                                                                       #
# Author: Catalina Vattuone  (cvattuonet@gmail.com)                                                                                     #        
# Date of latest update: 19-08-2026                                                                                                     #
# Main use: Study on "Assessing how agricultural practices and landscape composition shape differences in arthropod biodiversity among  #  
#           habitat types in Europe using the PREDICTS database". Internship at INRAE                                                   #
#                                                                                                                                       #
# Content: Extraction of figures of data characterization and main results.                                                             # 
#---------------------------------------------------------------------------------------------------------------------------------------#

#LIBRARIES  ------------------------------------------------------------
pacman::p_load(readr, dplyr, tidyr, stringr, forcats, purrr, lubridate,ggplot2, ggridges, patchwork, scales,sf, rnaturalearth, rnaturalearthdata,corrplot, MuMIn, glmmTMB, lme4)

#LOAD   ------------------------------------------------------------
predicts_df <- readRDS("Databases/PREDICTS/6fa1dedf-c546-41e0-a470-17c4863686b8.rds")

#datasets used for modelling
sites_df <- read_delim("Intermediate_dataset/sites_for_abundance_models.csv")
sites_df$Habitat_category <- relevel(factor(sites_df$Habitat_category), ref = "Natural vegetation")

sites_df_richness <- read_delim("Intermediate_dataset/sites_for_richness_models.csv")
sites_df_richness$Habitat_category <- relevel(factor(sites_df_richness$Habitat_category), ref = "Natural vegetation")

df_1000_viz <- read_delim("Intermediate_dataset/abundance_sites_for_vizualization_nonscale.csv")

observation_df <- read_delim("Intermediate_dataset/observations_mod.csv")

observation_df <- observation_df %>%
  mutate(SSBS = as.character(SSBS)) %>%
  left_join(sites_df %>% mutate(SSBS = as.character(SSBS), Habitat_category = as.character(Habitat_category)) %>% select(SSBS, Habitat_category), by = "SSBS")

#models
global_gaussian <- readRDS("Models_results/global_gaussian.rds")
dredge_gaussian <- readRDS("Models_results/dredge_gaussian.rds")

global_zinb_richness <- readRDS("Models_results/global_zinb_richness.rds")
dredge_zinb_richness <- readRDS("Models_results/dredge_zinb_richness.rds")

#OUTPUT DIRECTORY  ------------------------------------------------------------
out_dir <- "Figures/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

site_meta <- predicts_df %>%
  mutate(SSBS = as.character(SSBS)) %>%
  distinct(SSBS, Reference, Latitude, Longitude)

#SHARED LABELS, PALETTES AND FUNCTIONS ACROSS FIGURES ------------------------------------------------------------
#predictors labels
predictor_labels <- c(
  urban_pct_1000m_mod = "Urban cover",
  water_pct_1000m_mod = "Water cover",
  natural_pct_1000m_mod = "Natural habitat cover",
  pasture_pct_1000m_mod  = "Pasture cover",
  pesticide_1000m  = "Pesticide use",
  crop_diversity_machefer_alpha_1000m = "Crop diversity",
  pct_very_small_and_small_field_1000m = "Small field prevalence",
  pct_very_large_and_large_field_1000m = "Large field prevalence",
  climate_mean_temp_c_1000m = "Mean temperature")

interacting_predictors <- c("pesticide_1000m", "crop_diversity_machefer_alpha_1000m",
                            "pct_very_small_and_small_field_1000m", "natural_pct_1000m_mod",
                            "urban_pct_1000m_mod", "water_pct_1000m_mod")

#land use and use intensity labeld
land_use_levels      <- c("Natural vegetation", "Semi-natural vegetation", "Pasture", "Cropland")
intensity_levels     <- c("Reference", "Minimal", "Light", "Intense")
intensity_shapes     <- c("Reference" = 15, "Minimal" = 16, "Light" = 17, "Intense" = 18)
land_use_colors      <- c("Natural vegetation" = "darkgreen", "Semi-natural vegetation" = "#2166AC",
                          "Pasture" = "#E377C2", "Cropland" = "#E69F00")

land_use_short_labels <- c("Natural vegetation" = "Natural veg.", "Semi-natural vegetation" = "Semi-natural veg.",
                           "Pasture" = "Pasture", "Cropland" = "Cropland")

#shot names
habitat_code_map <- c(
  "Natural vegetation"                  = "NV",
  "Semi-natural vegetation_Minimal use" = "SNV-M",
  "Semi-natural vegetation_Light use"   = "SNV-L",
  "Semi-natural vegetation_Intense use" = "SNV-I",
  "Cropland_Minimal use"                = "C-M",
  "Cropland_Light use"                  = "C-L",
  "Cropland_Intense use"                = "C-I",
  "Pasture_Minimal use"                 = "P-M",
  "Pasture_Light use"                   = "P-L",
  "Pasture_Intense use"                 = "P-I")

hab_levels <- c("NV", "SNV-M", "SNV-L", "SNV-I", "C-M", "C-L", "C-I", "P-M", "P-L", "P-I")

habitat_colors <- c("NV"    = "#228B22", "SNV-M" = "#FFF7BC", "SNV-L" = "#FEE391", "SNV-I" = "#FEC44F",
                    "C-M"   = "#FFD59E", "C-L"   = "#FFBB78", "C-I"   = "#FF7F00",
                    "P-M"   = "#DEEBF7", "P-L"   = "#9ECAE1", "P-I"   = "#4292C6")

#habitat_category string (e.g. "Cropland_Intense use") to its short code (e.g. "C-I")
recode_habitat_codes <- function(x) {
  factor(unname(habitat_code_map[as.character(x)]), levels = hab_levels)
}

#split a Habitat_category term/level string into land_use + intensity
parse_habitat_level <- function(x) {
  data.frame(intensity = factor(case_when(
    str_detect(x, "Minimal") ~ "Minimal",
    str_detect(x, "Light")   ~ "Light",
    str_detect(x, "Intense") ~ "Intense",
    TRUE ~ "Reference"), levels = intensity_levels),
    land_use = str_remove(x, "_(Minimal|Light|Intense) use$"))
}

#pulls model-averaged coefficients + vcov out of a dredge object
get_coefs_vcov_avg <- function(dredge_obj, delta_cutoff = 2) {
  n_models <- nrow(subset(dredge_obj, delta < delta_cutoff))
  cat("Averaging", n_models, "models with delta <", delta_cutoff, "\n")
  
  avg <- model.avg(dredge_obj, subset = delta < delta_cutoff, fit = TRUE)
  b  <- coef(avg, full = TRUE)
  vc <- as.matrix(vcov(avg, full = TRUE))
  
  clean <- function(x) sub("^cond\\((.*)\\)$", "\\1", x)
  names(b) <- clean(names(b))
  dimnames(vc) <- list(clean(rownames(vc)), clean(colnames(vc)))
  
  keep <- !grepl("^zi\\(|^zi~", names(b))
  list(b = b[keep], vc = vc[keep, keep, drop = FALSE])
}

get_habitat_pctdiff <- function(cv) {
  terms <- grep("^Habitat_category[^:]+$", names(cv$b), value = TRUE)
  df <- data.frame(term = terms, estimate = as.numeric(cv$b[terms]), se = sqrt(diag(cv$vc)[terms])) %>%
    mutate(level = str_remove(term, "^Habitat_category")) %>%
    bind_cols(parse_habitat_level(.$level)) %>%
    mutate(conf.low = estimate - 1.96 * se,
           conf.high = estimate + 1.96 * se,
           pct_diff = (exp(estimate) - 1) * 100,
           pct_low  = (exp(conf.low) - 1) * 100,
           pct_high = (exp(conf.high) - 1) * 100)
  
  #add Natural vegetation itself as the 0% reference point
  bind_rows(data.frame(term = NA_character_, land_use = "Natural vegetation", intensity = "Reference",
                       pct_diff = 0, pct_low = 0, pct_high = 0),
            df %>% select(term, land_use, intensity, pct_diff, pct_low, pct_high)) %>%
    mutate(land_use = factor(land_use, levels = land_use_levels))
}

count_habitat_sites <- function(data) {
  counts <- data %>% count(Habitat_category, name = "n")
  bind_cols(counts, parse_habitat_level(as.character(counts$Habitat_category))) %>%
    select(land_use, intensity, n)
}

get_landscape_pctdiff <- function(cv) {
  terms <- setdiff(names(cv$b), c("(Intercept)", "(Int)", grep("Habitat_category", names(cv$b), value = TRUE)))
  terms <- terms[!str_detect(terms, ":")]
  
  data.frame(term = terms, estimate = as.numeric(cv$b[terms]), se = sqrt(diag(cv$vc)[terms])) %>%
    mutate(conf.low = estimate - 1.96 * se, conf.high = estimate + 1.96 * se,
           pct_diff = (exp(estimate) - 1) * 100,
           pct_low  = (exp(conf.low) - 1) * 100,
           pct_high = (exp(conf.high) - 1) * 100,
           label = recode(term, !!!predictor_labels),
           depends_on_habitat = purrr::map_lgl(term, ~ any(grepl(
             paste0("(^Habitat_category.*:", .x, "$)|(^", .x, ":Habitat_category.*$)"), names(cv$b)))),
           significant = pct_low > 0 | pct_high < 0,
           legend_group = case_when(
             depends_on_habitat & significant  ~ "Significant, but differs by habitat",
             !depends_on_habitat & significant ~ "Significant",
             TRUE                              ~ "Not significant"),
           legend_group = factor(legend_group, levels = c("Significant, but differs by habitat", "Significant", "Not significant")))
}

compute_simple_slopes <- function(cv, predictor) {
  b <- cv$b; vc <- cv$vc
  if (!(predictor %in% names(b))) return(NULL)
  int_terms <- grep(paste0("(^Habitat_category.*:", predictor, "$)|(^", predictor, ":Habitat_category.*$)"), names(b), value = TRUE)
  if (length(int_terms) == 0) return(NULL)   #this predictor's slope wasn't retained as varying by habitat
  levels_found <- str_remove(str_extract(int_terms, "Habitat_category[^:]+"), "^Habitat_category")
  ref_row <- data.frame(level = "Natural vegetation", slope = as.numeric(b[[predictor]]), se = sqrt(vc[predictor, predictor]))
  
  other_rows <- purrr::map2_dfr(int_terms, levels_found, function(term, lvl) {
    slope_val <- as.numeric(b[[predictor]] + b[[term]])
    var_val   <- vc[predictor, predictor] + vc[term, term] + 2 * vc[predictor, term]
    data.frame(level = lvl, slope = slope_val, se = sqrt(var_val))
  })
  
  bind_rows(ref_row, other_rows) %>%
    mutate(predictor = predictor,
           conf.low = slope - 1.96 * se, conf.high = slope + 1.96 * se,
           pct_diff = (exp(slope) - 1) * 100,
           pct_low  = (exp(conf.low) - 1) * 100,
           pct_high = (exp(conf.high) - 1) * 100) %>%
    bind_cols(parse_habitat_level(.$level))
}

get_interaction_slopes <- function(cv) {
  df <- purrr::map_dfr(interacting_predictors, ~compute_simple_slopes(cv, .x))
  present <- interacting_predictors[interacting_predictors %in% unique(df$predictor)]
  df %>%
    mutate(land_use = factor(land_use, levels = land_use_levels),
           label = factor(recode(predictor, !!!predictor_labels), levels = unname(predictor_labels[present])))
}

DELTA_CUTOFF <- 2
cv_abundance <- get_coefs_vcov_avg(dredge_gaussian, delta_cutoff = DELTA_CUTOFF)
cv_richness  <- get_coefs_vcov_avg(dredge_zinb_richness, delta_cutoff = DELTA_CUTOFF)

#Figure 1: Global sites distribution  ------------------------------------------------------------
world_map <- ne_countries(scale = "medium", returnclass = "sf")
coords_cleaned <- predicts_df %>%
  filter(!is.na(Latitude) & !is.na(Longitude)) %>%
  distinct(Longitude, Latitude)

sites_sf <- st_as_sf(coords_cleaned, coords = c("Longitude", "Latitude"), crs = 4326)

robinson_crs <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

lon_labels <- tibble(Lon = c(-120, -60, 60, 120), Lat = c(-40, -40, -40, -40),
                     Label = c("120°W", "60°W", "60°E", "120°E")) %>%
  st_as_sf(coords = c("Lon", "Lat"), crs = 4326)

global_sampling_map <- ggplot() +
  geom_sf(data = world_map, fill = "#F8FAFC", color = "#CBD5E1", size = 0.2) +
  geom_sf(data = sites_sf, color = "#0284C7", alpha = 0.6, size = 0.6, stroke = 0.1) +
  geom_sf_text(data = lon_labels, aes(label = Label), color = "#64748B",
               size = 2.4, fontface = "plain", vjust = -0.5) +
  coord_sf(crs = robinson_crs, datum = st_crs(4326), label_axes = list(left = "N")) +
  theme_minimal(base_size = 11) +
  labs(title = NULL, subtitle = NULL, x = NULL, y = NULL) +
  theme(axis.text.y = element_text(size = 8, color = "#334155"),
        axis.text.x = element_blank(),
        panel.grid.major = element_line(color = "#E2E8F0", size = 0.2),
        panel.background = element_rect(fill = "#FFFFFF", color = NA),
        plot.background = element_rect(fill = "#FFFFFF", color = NA),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt"))

ggsave(paste0(out_dir, "Figure1_global_sites_map.png"), global_sampling_map, width = 10, height = 5, dpi = 300)

#Figure 2: Europe sites distribution ------------------------------------------------------------
europe_lcc_crs <- 3034
coords_sites <- sites_df %>%
  mutate(SSBS = as.character(SSBS)) %>%
  left_join(site_meta %>% select(SSBS, Latitude, Longitude), by = "SSBS") %>%
  filter(!is.na(Latitude) & !is.na(Longitude)) %>%
  distinct(Longitude, Latitude)

sites_projected <- st_as_sf(coords_sites, coords = c("Longitude", "Latitude"), crs = 4326) %>%
  st_transform(crs = europe_lcc_crs)

world_projected <- st_transform(world_map, crs = europe_lcc_crs)

bbox_meters <- st_bbox(sites_projected)
x_min <- bbox_meters["xmin"] - 200000
x_max <- bbox_meters["xmax"] + 200000
y_min <- bbox_meters["ymin"] - 200000
y_max <- bbox_meters["ymax"] + 200000

europe_sampling_map <- ggplot() +
  geom_sf(data = world_projected, fill = "#F8FAFC", color = "#CBD5E1", size = 0.2) +
  geom_sf(data = sites_projected, color = "#0284C7", alpha = 0.7, size = 0.8, stroke = 0.1) +
  coord_sf(xlim = c(x_min, x_max), ylim = c(y_min, y_max), expand = FALSE, datum = st_crs(4326),
           label_axes = list(lower = "E", left = "N")) +
  theme_minimal(base_size = 11) +
  labs(title = NULL, subtitle = NULL, x = NULL, y = NULL) +
  theme(axis.text = element_text(size = 8, color = "#334155"),
        panel.grid.major = element_line(color = "#E2E8F0", size = 0.2),
        panel.background = element_rect(fill = "#FFFFFF", color = NA),
        plot.background = element_rect(fill = "#FFFFFF", color = NA),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt"))

ggsave(paste0(out_dir, "Figure2_europe_sites_map.png"), europe_sampling_map, width = 8, height = 8, dpi = 300)

#Figure 3: Main results - habitat category effects ------------------------------------------------------------
habitat_abundance <- get_habitat_pctdiff(cv_abundance) %>%
  left_join(count_habitat_sites(sites_df), by = c("land_use", "intensity")) %>%
  mutate(land_use = factor(land_use, levels = land_use_levels))

habitat_richness <- get_habitat_pctdiff(cv_richness) %>%
  left_join(count_habitat_sites(sites_df_richness), by = c("land_use", "intensity")) %>%
  mutate(land_use = factor(land_use, levels = land_use_levels))

#hand-tuned offsets so the "(n)" labels don't collide with the error bars - adjust if the
#underlying effect sizes change enough to shift where the bars land
label_y_bump_abundance <- c(
  "Natural vegetation Reference"    = 0.15,
  "Semi-natural vegetation Minimal" = -0.45,
  "Semi-natural vegetation Light"   = 0.05,
  "Semi-natural vegetation Intense" = -1.35,
  "Pasture Minimal"                 = 0.05,
  "Pasture Light"                   = -1.2,
  "Pasture Intense"                 = -0.9,
  "Cropland Minimal"                = 0.1,
  "Cropland Light"                  = -0.6,
  "Cropland Intense"                = 0.1)

label_y_bump_richness <- c(
  "Natural vegetation Reference"    = 0.2,
  "Semi-natural vegetation Minimal" = -0.25,
  "Semi-natural vegetation Light"   = 0.1,
  "Semi-natural vegetation Intense" = -1.4,
  "Pasture Minimal"                 = 0.1,
  "Pasture Light"                   = -1,
  "Pasture Intense"                 = -0.45,
  "Cropland Minimal"                = 0.2,
  "Cropland Light"                  = -0.8,
  "Cropland Intense"                = 0.2)

label_hjust <- c(
  "Natural vegetation Reference"    = 0.5,
  "Semi-natural vegetation Minimal" = 0.15,
  "Semi-natural vegetation Light"   = 0.5,
  "Semi-natural vegetation Intense" = 0.7,
  "Pasture Minimal"                 = 0.5,
  "Pasture Light"                   = 0.2,
  "Pasture Intense"                 = 0.8,
  "Cropland Minimal"                = 0.2,
  "Cropland Light"                  = 0.5,
  "Cropland Intense"                = 0.8)

add_label_positions <- function(df, bump) {
  df %>%
    group_by(land_use) %>%
    mutate(group_max = max(pct_high, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(overall_max = max(pct_high, na.rm = TRUE),
           label_y  = group_max + overall_max * (0.08 + bump[paste(land_use, intensity)]),
           hjust_val = label_hjust[paste(land_use, intensity)])
}

habitat_abundance <- add_label_positions(habitat_abundance, label_y_bump_abundance)
habitat_richness  <- add_label_positions(habitat_richness,  label_y_bump_richness)

plot_habitat_effects <- function(df, ylab, y_break) {
  ggplot(df, aes(x = land_use, y = pct_diff, color = land_use, shape = intensity, group = intensity)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    geom_errorbar(aes(ymin = pct_low, ymax = pct_high), width = 0.3, position = position_dodge(width = 0.5)) +
    geom_point(size = 5, position = position_dodge(width = 0.5)) +
    geom_text(aes(y = label_y, label = paste0("(", n, ")"), hjust = hjust_val),
              position = position_dodge(width = 0.5), size = 5, color = "black", show.legend = FALSE) +
    scale_color_manual(values = land_use_colors, name = "Land use", guide = "none") +
    scale_shape_manual(values = intensity_shapes, name = "Land-use intensity") +
    scale_x_discrete(labels = land_use_short_labels) +
    scale_y_continuous(breaks = scales::breaks_width(y_break), expand = expansion(mult = c(0.05, 0.16))) +
    theme_minimal(base_size = 22) +
    theme(panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 30, hjust = 1),
          axis.title.y = element_text(size = 17)) +
    labs(x = NULL, y = ylab)
}

p_habitat_abundance <- plot_habitat_effects(habitat_abundance, "Abundance difference (%)", 50)
p_habitat_richness  <- plot_habitat_effects(habitat_richness,  "Taxa richness difference (%)", 25)

habitat_combined <- (p_habitat_abundance / p_habitat_richness) +
  plot_layout(guides = "collect") &
  theme(legend.position = "top",
        legend.title = element_text(size = 15, face = "plain"),
        legend.text = element_text(size = 15, face = "plain"),
        legend.key.size = unit(0.7, "cm"))

ggsave(paste0(out_dir, "Figure3_habitat_effects.png"), habitat_combined, width = 12, height = 14, dpi = 300)

#Figure 4: Main results - landscape/practice main effects ------------------------------------------------------------
landscape_abundance <- get_landscape_pctdiff(cv_abundance)
landscape_richness  <- get_landscape_pctdiff(cv_richness)

#use the abundance-model ordering for both panels, so the same predictor lines up in the same position
landscape_order <- landscape_abundance %>% arrange(pct_diff) %>% pull(label) %>% as.character()
landscape_abundance <- landscape_abundance %>% mutate(label = factor(label, levels = landscape_order))
landscape_richness  <- landscape_richness  %>% mutate(label = factor(label, levels = landscape_order))

plot_landscape_effects <- function(df, ylab) {
  group_colors <- c("Significant, but differs by habitat" = "#00A3A5", "Significant" = "orange", "Not significant" = "darkgrey")
  group_alphas <- c("Significant, but differs by habitat" = 1, "Significant" = 1, "Not significant" = 1)
  
  ggplot(df, aes(x = label, y = pct_diff, color = legend_group, alpha = legend_group)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    geom_errorbar(aes(ymin = pct_low, ymax = pct_high), width = 0.15) +
    geom_point(size = 3) +
    scale_color_manual(values = group_colors, name = "Effect") +
    scale_alpha_manual(values = group_alphas, name = "Effect") +
    theme_minimal(base_size = 22) +
    theme(panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 30, hjust = 1),
          axis.title.y = element_text(size = 17)) +
    labs(x = NULL, y = ylab)
}

p_landscape_abundance <- plot_landscape_effects(landscape_abundance, "Abundance difference (%)")
p_landscape_richness  <- plot_landscape_effects(landscape_richness,  "Taxa richness difference (%)")

landscape_combined <- (p_landscape_abundance / plot_spacer() / p_landscape_richness) +
  plot_layout(heights = c(1, 0.1, 1), guides = "collect") &
  theme(legend.position = "top",
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.key.size = unit(0.35, "cm"),
        legend.box.spacing = unit(30, "pt"))

ggsave(paste0(out_dir, "Figure4_landscape_effects.png"), landscape_combined, width = 12, height = 14, dpi = 300)

#Figure 5: Main results - interaction (simple slopes by habitat) ------------------------------------------------------------
interactions_abundance <- get_interaction_slopes(cv_abundance)
interactions_richness  <- get_interaction_slopes(cv_richness)

plot_interaction_slopes <- function(df, ylab) {
  ggplot(df, aes(x = land_use, y = pct_diff, color = land_use, shape = intensity, group = intensity)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    geom_errorbar(aes(ymin = pct_low, ymax = pct_high), width = 0.17, position = position_dodge(width = 0.5)) +
    geom_point(size = 4.5, position = position_dodge(width = 0.5)) +
    facet_wrap(~label, nrow = 1, scales = "fixed") +
    scale_color_manual(values = land_use_colors, name = "Land use", guide = "none") +
    scale_shape_manual(values = intensity_shapes, name = "Land-use intensity") +
    scale_x_discrete(labels = land_use_short_labels) +
    theme_minimal(base_size = 20) +
    theme(panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 30, hjust = 1),
          axis.title.y = element_text(size = 17),
          strip.text = element_text(face = "bold", size = 13),
          strip.background = element_blank(),
          panel.spacing = unit(2, "lines"),
          legend.position = "top",
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 12),
          legend.key.size = unit(0.35, "cm")) +
    labs(x = NULL, y = ylab)
}

p_int_abundance <- plot_interaction_slopes(interactions_abundance, "Abundance change (%)")
p_int_richness  <- plot_interaction_slopes(interactions_richness,  "Taxa richness change (%)")

p_int_combined <- (p_int_abundance / plot_spacer() / p_int_richness) +
  plot_layout(heights = c(1, 0.12, 1), guides = "collect") &
  theme(legend.position = "top",
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.key.size = unit(0.7, "cm"))

ggsave(paste0(out_dir, "Figure5_interaction_slopes.png"), p_int_combined, width = 16, height = 10, dpi = 300)


# Figure 6: Taxonomic composition of the dataset------------------------------------------------------------
taxa_clean_data <- observation_df %>%
  filter(Measurement > 0) %>%
  filter(!is.na(Order), !is.na(Family)) %>%
  group_by(Order, Family) %>%
  tally(name = "obs_count") %>%
  ungroup() %>%
  mutate(percentage = (obs_count / sum(obs_count)) * 100) %>%   #percentage share of total records
  filter(percentage >= 0.5) %>%                                 #only families representing >= 0.5%, to keep it clean
  group_by(Order) %>%
  mutate(order_total = sum(percentage)) %>%
  ungroup() %>%
  mutate(Family = reorder(Family, percentage))

plot_taxa_grid <- ggplot(taxa_clean_data, aes(x = percentage, y = Family, fill = Order)) +
  geom_col(width = 0.8, alpha = 0.9, show.legend = FALSE) +
  geom_text(aes(label = paste0(round(percentage, 1), "%")), hjust = -0.1, size = 3, color = "#334155") +
  facet_wrap(~ Order, scales = "free_y", ncol = 2) +
  scale_x_continuous(labels = label_percent(scale = 1), expand = expansion(mult = c(0, 0.15))) +
  scale_fill_brewer(palette = "Set2") +
  theme_bw(base_size = 11) +
  labs(title = "Taxonomic composition of dataset",
       subtitle = "Blocks represent major Orders (>1% of total share); bars show individual Family percentage shares",
       x = "Percentage of total records", y = "") +
  theme(plot.title = element_text(face = "bold", size = 13, color = "#0F172A"),
        plot.subtitle = element_text(size = 10, color = "#475569"),
        strip.text = element_text(face = "bold", size = 11, color = "#0F172A"),
        strip.background = element_rect(fill = "#F1F5F9", color = "#CBD5E1"),
        axis.text.y = element_text(size = 9, color = "#1E293B"),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(1, "lines"))

ggsave(paste0(out_dir, "Figure6_taxonomic_composition.png"), plot_taxa_grid, width = 10, height = 14, dpi = 300)


# Figure 7: Taxonomic composition by habitat (heatmap)------------------------------------------------------------
europe_bias <- observation_df %>%
  filter(Use_intensity != "Cannot decide", !is.na(Habitat_category)) %>%
  mutate(Habitat_category = ifelse(
    Habitat_category %in% c("Natural vegetation_Light use", "Natural vegetation_Minimal use"),
    "Natural vegetation", Habitat_category))

habitat_levels_full <- c("Natural vegetation",
                         "Semi-natural vegetation_Minimal use", "Semi-natural vegetation_Light use",
                         "Semi-natural vegetation_Intense use",
                         "Cropland_Minimal use", "Cropland_Light use", "Cropland_Intense use",
                         "Pasture_Minimal use", "Pasture_Light use", "Pasture_Intense use")

europe_bias <- europe_bias %>% mutate(Habitat_category = factor(Habitat_category, levels = habitat_levels_full))

order_long <- europe_bias %>%
  filter(!is.na(Order), Order != "") %>%
  group_by(Habitat_category, Order) %>%
  summarise(n_records = n(), .groups = "drop") %>%
  mutate(Habitat_category = recode_habitat_codes(as.character(Habitat_category)))

taxonomic_bias_heatmap <- order_long %>%
  group_by(Habitat_category) %>%
  mutate(pct = 100 * n_records / sum(n_records)) %>%
  ungroup() %>%
  mutate(Order = fct_reorder(Order, pct, .fun = sum),
         label = case_when(
           pct == 0 ~ "",           #genuinely absent
           pct < 1  ~ "<1",         #present but rare
           TRUE     ~ as.character(round(pct, 0))),
         fill_val = ifelse(pct == 0, NA, pct)) %>%   #NA -> distinct colour for true zeros
  ggplot(aes(x = Habitat_category, y = Order, fill = fill_val)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = label, color = pct > 45), size = 2.7, fontface = "bold") +
  scale_fill_gradient(low = "#F2FAFB", high = "#00A3A5", na.value = "grey96", name = "% of records\nin habitat") +
  scale_color_manual(values = c("TRUE" = "white", "FALSE" = "grey15"), guide = "none") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank()) +
  labs(x = NULL, y = NULL,
       title = "Order composition of each habitat category",
       subtitle = "Each column sums to 100%. Grey = taxon absent; <1 = present but under 1% of records.")

ggsave(paste0(out_dir, "Figure7_taxonomic_bias_heatmap.png"), taxonomic_bias_heatmap, width = 9, height = 8, dpi = 300)


# Figure 8: Temporal distribution of sampling events by reference------------------------------------------------------------
reference_metadata <- observation_df %>%
  mutate(Lat_Long_mod = paste(Latitude_mod, Longitude_mod, sep = "_"),
         Year_Extracted = year(as.Date(Sample_start_earliest))) %>%
  filter(!is.na(Year_Extracted)) %>%
  group_by(Reference) %>%
  summarise(start_yr = min(Year_Extracted, na.rm = TRUE),
            end_yr   = max(Year_Extracted, na.rm = TRUE),
            unique_coor_count = n_distinct(Lat_Long_mod), .groups = "drop")

sampling_events <- observation_df %>%
  mutate(Year_Extracted = year(as.Date(Sample_start_earliest))) %>%
  filter(!is.na(Year_Extracted)) %>%
  select(Reference, Year_Extracted) %>%
  distinct()

timeline_data <- sampling_events %>%
  left_join(reference_metadata, by = "Reference") %>%
  mutate(Ref_Label = str_trunc(Reference, 40)) %>%
  mutate(Ref_Label = reorder(Ref_Label, start_yr))

plot_event_timeline_reference <- ggplot(timeline_data) +
  geom_hline(yintercept = seq_along(levels(timeline_data$Ref_Label)), color = "#E2E8F0", size = 0.4, linetype = "dashed") +
  geom_segment(data = . %>% select(Ref_Label, start_yr, end_yr) %>% distinct(),
               aes(x = start_yr, xend = end_yr, y = Ref_Label, yend = Ref_Label),
               color = "#00A3A5", size = 1.0, alpha = 0.7) +
  geom_point(aes(x = Year_Extracted, y = Ref_Label), color = "#00A3A5", size = 2.5, alpha = 0.9) +
  scale_x_continuous(breaks = seq(min(timeline_data$start_yr), max(timeline_data$end_yr), by = 2)) +
  theme_bw(base_size = 10) +
  labs(title = "Timeline covered by Reference",
       subtitle = "Lines show study duration; dots indicate active sampling years",
       x = "Year", y = "") +
  theme(plot.title = element_text(face = "bold", size = 13, color = "#0F172A"),
        plot.subtitle = element_text(size = 10, color = "#475569"),
        axis.text.y = element_text(size = 8.5, color = "#334155", family = "mono"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank())

n_refs <- n_distinct(timeline_data$Ref_Label)
ggsave(paste0(out_dir, "Figure8_temporal_distribution.png"), plot_event_timeline_reference,
       width = 10, height = max(8, n_refs * 0.28), dpi = 300, limitsize = FALSE)


# Figure 9: Response distribution by habitat, as density ------------------------------------------------------------
dens_abundance <- sites_df %>%
  transmute(Habitat_category = recode_habitat_codes(Habitat_category), Value = log_TA,
            Metric = "Log Abundance log(TA + 1)")

dens_richness <- sites_df_richness %>%
  transmute(Habitat_category = recode_habitat_codes(Habitat_category), Value = taxa_richness,
            Metric = "Taxa Richness")

fig9_hab_order <- c("P-M", "P-L", "P-I", "C-M", "C-L", "C-I", "SNV-M", "SNV-L", "SNV-I", "NV")

df_distribution <- bind_rows(dens_abundance, dens_richness) %>%
  filter(!is.na(Habitat_category), !is.na(Value)) %>%
  mutate(Habitat_category = factor(as.character(Habitat_category), levels = fig9_hab_order),
         Metric = factor(Metric, levels = c("Log Abundance log(TA + 1)", "Taxa Richness")))

n_counts <- df_distribution %>%
  group_by(Metric, Habitat_category) %>%
  summarise(n = n(), max_val = max(Value, na.rm = TRUE), .groups = "drop") %>%
  mutate(label = paste0("n=", n))

response_density_by_habitat <- ggplot(df_distribution, aes(x = Value, y = Habitat_category, fill = Habitat_category)) +
  geom_density_ridges(alpha = 0.85, scale = 1.1, color = "white", linewidth = 0.3) +
  facet_wrap(~Metric, scales = "free_x", ncol = 2) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "grey90", color = "grey80"),
        strip.text = element_text(face = "bold", size = 9),
        legend.position = "none",
        axis.title.y = element_blank()) +
  geom_text(data = n_counts, aes(x = max_val + (0.1 * max_val), y = Habitat_category, label = label),
            hjust = 0, size = 3, inherit.aes = FALSE) +
  scale_fill_manual(values = habitat_colors) +
  labs(x = "Value", y = "Habitat Category") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.2)))

ggsave(paste0(out_dir, "Figure9_response_density_by_habitat.png"), response_density_by_habitat, width = 10, height = 8, dpi = 300)

# Figure 10: abundance by reference  ------------------------------------------------------------
abundance_by_ref <- sites_df %>%
  mutate(SSBS = as.character(SSBS)) %>%
  left_join(site_meta %>% select(SSBS, Reference), by = "SSBS") %>%
  filter(!is.na(Reference))

abundance_by_reference_plot <- ggplot(abundance_by_ref, aes(x = Reference, y = log_TA)) +
  geom_boxplot(outlier.size = 1, fill = "#00A3A5", color = "black", linewidth = 0.3) +
  coord_flip() +
  theme_bw(base_size = 10) +
  labs(x = "Reference", y = "Log(abundance + 1)",
       title = "Distribution of site-level abundance by reference")

n_ref_ab <- n_distinct(abundance_by_ref$Reference)
ggsave(paste0(out_dir, "Figure10_abundance_by_reference.png"), abundance_by_reference_plot,
       width = 8, height = max(8, n_ref_ab * 0.25), dpi = 300, limitsize = FALSE)

# Figure 11: richness by reference  ------------------------------------------------------------
richness_by_ref <- sites_df_richness %>%
  mutate(SSBS = as.character(SSBS)) %>%
  left_join(site_meta %>% select(SSBS, Reference), by = "SSBS") %>%
  filter(!is.na(Reference))

richness_by_reference_plot <- ggplot(richness_by_ref, aes(x = Reference, y = taxa_richness)) +
  geom_boxplot(outlier.size = 1, fill = "#00A3A5", color = "black", linewidth = 0.3) +
  coord_flip() +
  theme_bw(base_size = 10) +
  labs(x = "Reference", y = "Taxa richness",
       title = "Distribution of site-level taxa richness by reference")

n_ref_rich <- n_distinct(richness_by_ref$Reference)
ggsave(paste0(out_dir, "Figure11_richness_by_reference.png"), richness_by_reference_plot,
       width = 8, height = max(8, n_ref_rich * 0.25), dpi = 300, limitsize = FALSE)


# Figure 12: Each landscape predictor by habitat category  ------------------------------------------------------------
df_long <- df_1000_viz %>%
  select(Habitat_category, urban_pct_1000m_mod, water_pct_1000m_mod,
         natural_pct_1000m_mod, pasture_pct_1000m_mod, crop_pct_1000m,
         pesticide_1000m, fertilizer_1000m, crop_diversity_shannon_1000m, crop_diversity_machefer_alpha_1000m,
         pct_very_small_and_small_field_1000m, pct_very_large_and_large_field_1000m,
         climate_annual_precip_mm_1000m, climate_mean_temp_c_1000m) %>%
  pivot_longer(cols = -Habitat_category, names_to = "Predictor", values_to = "Value")

df_viz <- df_long %>%
  mutate(Habitat_category = recode_habitat_codes(Habitat_category)) %>%
  filter(!is.na(Habitat_category)) %>%
  mutate(Habitat_category = factor(as.character(Habitat_category), levels = rev(hab_levels)),
         Predictor = factor(Predictor,
                            levels = c("urban_pct_1000m_mod", "water_pct_1000m_mod", "natural_pct_1000m_mod", "pasture_pct_1000m_mod", "crop_pct_1000m", "pesticide_1000m",
                                       "fertilizer_1000m", "crop_diversity_shannon_1000m", "crop_diversity_machefer_alpha_1000m", "pct_very_small_and_small_field_1000m", "pct_very_large_and_large_field_1000m",
                                       "climate_mean_temp_c_1000m", "climate_annual_precip_mm_1000m"),
                            labels = c("Urban cover %", "Water cover %", "Natural vegetation cover %", "Pasture cover %", "Cropland cover %", "Average pesticide (kg/ha-yr)",
                                       "Average N fertilizer (kg/ha-yr)", "Crops Shannon Index", "Crop Alpha Diversity", "Small fields prevalence %", "Large fields prevalence %",
                                       "Mean temperature (°C)", "Annual precipitation (mm)")))

predictors_by_habitat_plot <- ggplot(df_viz, aes(x = Habitat_category, y = Value, fill = Habitat_category)) +
  geom_boxplot(outlier.size = 0.5, linewidth = 0.2) +
  scale_fill_manual(values = habitat_colors, drop = TRUE) +
  facet_wrap(~Predictor, scales = "free", ncol = 3) +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "none") +
  labs(title = "Distribution of landscape predictors (1,000 m radius) by habitat",
       y = "Value", x = "Local Habitat Category")

ggsave(paste0(out_dir, "Figure12_predictors_by_habitat.png"), predictors_by_habitat_plot, width = 12, height = 14, dpi = 300)

# Figure 13: Correlation matrix of predictors  ------------------------------------------------------------
continuous_preds <- c("urban_pct_1000m_mod", "water_pct_1000m_mod",
                      "natural_pct_1000m_mod", "pasture_pct_1000m_mod", "crop_pct_1000m",
                      "pesticide_1000m", "fertilizer_1000m",
                      "crop_diversity_shannon_1000m", "crop_diversity_machefer_alpha_1000m",
                      "pct_very_small_and_small_field_1000m", "pct_very_large_and_large_field_1000m",
                      "climate_mean_temp_c_1000m", "climate_annual_precip_mm_1000m")

cor_data_renamed <- sites_df %>%
  distinct(Coordinate_ID, .keep_all = TRUE) %>%
  select(all_of(continuous_preds)) %>%
  mutate(across(everything(), as.numeric)) %>%
  rename("Urban cover" = "urban_pct_1000m_mod",
         "Water cover" = "water_pct_1000m_mod",
         "Natural vegetation cover" = "natural_pct_1000m_mod",
         "Pasture cover" = "pasture_pct_1000m_mod",
         "Cropland cover" = "crop_pct_1000m",
         "Pesticide" = "pesticide_1000m",
         "N Fertilizer" = "fertilizer_1000m",
         "Crops Shannon Index" = "crop_diversity_shannon_1000m",
         "Crop Alpha Diversity" = "crop_diversity_machefer_alpha_1000m",
         "Small fields prevalence" = "pct_very_small_and_small_field_1000m",
         "Large fields prevalence" = "pct_very_large_and_large_field_1000m",
         "Mean temperature" = "climate_mean_temp_c_1000m",
         "Annual precipitation" = "climate_annual_precip_mm_1000m")

cor_matrix <- cor(cor_data_renamed, use = "complete.obs", method = "pearson")
print(round(cor_matrix, 2))

png(paste0(out_dir, "Figure13_correlation_matrix.png"), width = 8, height = 8, units = "in", res = 300)
corrplot(cor_matrix, method = "color", type = "upper",
         order = "hclust",
         addCoef.col = "black",
         number.cex = 0.8,
         tl.col = "black", tl.srt = 45, tl.cex = 0.7,
         col = colorRampPalette(c("#B2182B", "white", "#2166AC"))(200),
         mar = c(0, 0, 1, 0))
dev.off()


# Figure 14: Significant predictors by study  ------------------------------------------------------------
var_labels <- c(
  urban_pct_1000m_mod                  = "Urban cover (%)",
  water_pct_1000m_mod                  = "Water cover (%)",
  natural_pct_1000m_mod                = "Natural habitat cover (%)",
  pasture_pct_1000m_mod                = "Pasture cover (%)",
  crop_pct_1000m                       = "Cropland cover (%)",
  pesticide_1000m                      = "Average pesticide application rate (kg/ha-year)",
  fertilizer_1000m                     = "Average N fertilizer application rate",
  crop_diversity_shannon_1000m         = "Crop diversity (Shannon index)",
  crop_diversity_machefer_alpha_1000m  = "Crop diversity (alpha index)",
  pct_very_small_and_small_field_1000m = "Small field prevalence (%)",
  pct_very_large_and_large_field_1000m = "Large field prevalence (%)",
  climate_mean_temp_c_1000m            = "Mean temperature (°C)",
  climate_annual_precip_mm_1000m       = "Annual precipitation (mm)")

plot_var_by_study <- function(data, var, label) {
  ggplot(data, aes(x = reorder(SS, .data[[var]], FUN = median, na.rm = TRUE), y = .data[[var]])) +
    geom_boxplot(outlier.size = 0.8, fill = "#2166AC", alpha = 0.6, na.rm = TRUE) +
    coord_flip() +
    theme_minimal(base_size = 10) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
    labs(x = "Study (SS)", y = label,
         title = paste0("Distribution of ", label, " in buffer of sampled sites by study"))
}

for (var in names(var_labels)) {
  p <- plot_var_by_study(df_1000_viz, var, var_labels[[var]])
  ggsave(filename = paste0(out_dir, "Figure14_dist_by_study_", var, ".png"), plot = p, width = 8, height = 10, dpi = 300)
}