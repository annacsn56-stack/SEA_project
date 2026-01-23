# ===============================================
# FIGURE – Cambodia 2016–2022
# (Version PRO v7 - Polices "Poster")
# ===============================================

# --- Packages ---
library(dplyr)
library(ggplot2)
library(sf)
library(readr)
library(rnaturalearth)
library(rnaturalearthdata)
library(cowplot)
library(RColorBrewer)

# --- (optionnel) ton dossier de travail ---
setwd("/media/annacosson/ANNA_DD/STAGE_M2/Projet_Asie_du_Sud_Est/Article/Map/Data/")

# +++ PÉRIODES DE TRAITEMENT (CORRIGÉES) +++
treatment_data <- data.frame(
  policy = factor(c("DHA-PPQ", "AS-MQ"), levels = c("DHA-PPQ", "AS-MQ")),
  start = c(2015.5, 2018.0), # <-- MODIFIÉ
  end = c(2018.0, 2022.5)  # <-- MODIFIÉ
)

# +++ Données d'incidence (API/1000) +++
incidence_data <- data.frame(
  Year = 2016:2022,
  API_1000 = c(1.6, 3.13, 4.18, 2.1, 0.59, 0.27, 0.25)
)

# --- 1) Lire le fichier propre ---
df <- read.csv("Final_data/metadata_merged_with_CRT_with_KEL1.csv",
               header = TRUE, sep = ",", dec = ",",
               stringsAsFactors = FALSE)

# Nettoyage basique + filtre années
df <- df %>%
  mutate(
    Year = as.integer(Year),
    Province = trimws(Province),
    K13 = toupper(trimws(K13))
  ) %>%
  filter(Year >= 2016, Year <= 2022)

# --- 1b) Calculer le N total par année ---
year_totals <- df %>%
  count(Year, name = "N_Total")

# --- 2) Agrégations par province × année ---
agg <- df %>%
  group_by(Province, Year) %>%
  summarise(
    # K13
    N_k13 = sum(!is.na(K13) & trimws(K13) != ""),
    N_count_C580Y = sum(K13 == "C580Y", na.rm = TRUE),
    N_count_Y493H = sum(K13 == "Y493H", na.rm = TRUE),
    N_count_WT    = sum(K13 == "WT",    na.rm = TRUE),
    pct_C580Y = 100 * N_count_C580Y / ifelse(N_k13 > 0, N_k13, NA),
    pct_Y493H = 100 * N_count_Y493H / ifelse(N_k13 > 0, N_k13, NA),
    pct_WT    = 100 * N_count_WT    / ifelse(N_k13 > 0, N_k13, NA),
    
    # PMII
    N_pmii = sum(!is.na(PMII) & trimws(PMII) != ""),
    N_count_PMII_multi = sum(tolower(PMII) == "multi", na.rm = TRUE),
    pct_PMII_multi = 100 * N_count_PMII_multi / ifelse(N_pmii > 0, N_pmii, NA),
    
    # MDR1
    N_mdr1 = sum(!is.na(MDR1) & trimws(MDR1) != ""),
    N_count_MDR1_multi = sum(tolower(MDR1) == "multi", na.rm = TRUE),
    pct_MDR1_multi = 100 * N_count_MDR1_multi / ifelse(N_mdr1 > 0, N_mdr1, NA),
    
    # CRT
    N_crt = sum(!is.na(CRT_PPQ_status) & trimws(CRT_PPQ_status) != ""),
    N_count_CRT_mutant = sum(tolower(CRT_PPQ_status) == "mutant", na.rm = TRUE),
    pct_CRT_mutant = 100 * N_count_CRT_mutant / ifelse(N_crt > 0, N_crt, NA),
    
    # KEL1
    N_kel1 = sum(!is.na(KEL1_PLA1_proxy) & trimws(KEL1_PLA1_proxy) != ""),
    N_count_KEL1_PLA1 = sum(toupper(KEL1_PLA1_proxy) == "KEL1/PLA1", na.rm = TRUE),
    pct_KEL1_PLA1 = 100 * N_count_KEL1_PLA1 / ifelse(N_kel1 > 0, N_kel1, NA),
    
    .groups = "drop"
  )

# --- 4) Géodonnées & normalisation des noms ---
norm_name <- function(x){
  x |>
    iconv(to = "ASCII//TRANSLIT") |>
    tolower() |>
    gsub("[^a-z0-9]+"," ", x = _) |>
    trimws()
}

provinces_kh <- rnaturalearth::ne_states(country = "Cambodia", returnclass = "sf") %>%
  mutate(name_norm = norm_name(name))

alias_map <- c(
  "kampong som"      = "krong preah sihanouk",
  "kampong speu"     = "kampong spoe",
  "kratie"           = "kracheh",
  "mondolkiri"       = "mondol kiri",
  "oddar manchey"    = "otdar mean chey",
  "pursat"           = "pouthisat",
  "rattanakiri"      = "rotanokiri",
  "stung treng"      = "stoeng treng",
  "sieam reap"       = "siemreab",
  "preah sihanouk"   = "krong preah sihanouk",
  "preah vihear"     = "preah vihear",
  "krong preah vihear" = "preah vihear",
  "pailin"           = "krong pailin"
)

agg_norm <- agg %>%
  mutate(name_norm = norm_name(Province)) %>%
  filter(!name_norm %in% c("chhoukmeas","promoy")) %>%
  mutate(name_norm = dplyr::recode(name_norm, !!!alias_map))


# --- 5) Palettes & thème ---
fill_pal_K13  <- colorRampPalette(brewer.pal(9, "YlOrRd"))
fill_pal_PMII <- colorRampPalette(brewer.pal(9, "PuBu"))
fill_pal_MDR1 <- colorRampPalette(brewer.pal(9, "YlGn"))
fill_pal_CRT  <- colorRampPalette(RColorBrewer::brewer.pal(9, "Purples"))
fill_pal_KEL1 <- colorRampPalette(c("#fde0dd", "#f768a1", "#c51b8a"))

# --- MODIFIÉ : Ajout de la taille des barres de légende ---
base_theme <- theme_minimal(base_size = 19) + # (Garde la grande taille de police)
  theme(panel.grid = element_blank(),
        strip.text = element_text(face = "bold"),
        plot.title  = element_text(face = "bold", size = 21),
        legend.position = "right",
        # Supprime les coordonnées des cartes
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        
        # +++ NOUVEAU : Agrandir les barres de légende +++
        # (Tu peux ajuster les valeurs 3.0 et 1.2 si c'est trop/pas assez)
        legend.key.height = unit(1.5, "cm"), # Hauteur de la barre
        legend.key.width = unit(1.5, "cm")  # Épaisseur de la barre
        # +++ FIN NOUVEAU +++
  )

# --- 6) Carte unitaire ---
map_one <- function(df_sf, fill_var, title,
                    limits = c(0,100),
                    legend_title = "%",
                    show_title = TRUE,
                    show_legend = TRUE) {
  
  pal <- if (fill_var == "pct_PMII_multi") {
    fill_pal_PMII
  } else if (fill_var == "pct_MDR1_multi") {
    fill_pal_MDR1
  } else if (fill_var == "pct_CRT_mutant") {
    fill_pal_CRT
  } else if (fill_var == "pct_KEL1_PLA1") {
    fill_pal_KEL1
  } else {
    fill_pal_K13
  }
  
  num_label <- if (fill_var == "pct_C580Y") {
    sum(df_sf$N_count_C580Y, na.rm = TRUE)
  } else if (fill_var == "pct_Y493H") {
    sum(df_sf$N_count_Y493H, na.rm = TRUE)
  } else if (fill_var == "pct_WT") {
    sum(df_sf$N_count_WT, na.rm = TRUE)
  } else if (fill_var == "pct_PMII_multi") {
    sum(df_sf$N_count_PMII_multi, na.rm = TRUE)
  } else if (fill_var == "pct_MDR1_multi") {
    sum(df_sf$N_count_MDR1_multi, na.rm = TRUE)
  } else if (fill_var == "pct_CRT_mutant") {
    sum(df_sf$N_count_CRT_mutant, na.rm = TRUE)
  } else if (fill_var == "pct_KEL1_PLA1") {
    sum(df_sf$N_count_KEL1_PLA1, na.rm = TRUE)
  } else NA
  
  den_label <- if (fill_var %in% c("pct_C580Y","pct_Y493H","pct_WT")) {
    sum(df_sf$N_k13,  na.rm = TRUE)
  } else if (fill_var == "pct_PMII_multi") {
    sum(df_sf$N_pmii, na.rm = TRUE)
  } else if (fill_var == "pct_MDR1_multi") {
    sum(df_sf$N_mdr1, na.rm = TRUE)
  } else if (fill_var == "pct_CRT_mutant") {
    sum(df_sf$N_crt,  na.rm = TRUE)
  } else if (fill_var == "pct_KEL1_PLA1") {
    sum(df_sf$N_kel1, na.rm = TRUE)
  } else NA
  
  N_label <- if (!is.na(num_label) && !is.na(den_label)) {
    paste0("N = ", num_label, " / ", den_label)
  } else NULL
  
  leg <- if (fill_var %in% c("pct_PMII_multi","pct_MDR1_multi")) {
    "%"
  } else {
    legend_title
  }
  
  ggplot() +
    geom_sf(data = provinces_kh, fill = NA, color = "grey70", linewidth = 0.25) +
    geom_sf(data = df_sf, aes(fill = .data[[fill_var]]), color = "grey50", linewidth = 0.2) +
    scale_fill_gradientn(name = leg, colours = pal(100),
                         limits = c(0,100), na.value = "#cccccc") +
    labs(
      title = if(show_title) title else NULL,
      subtitle = N_label
    ) +
    base_theme + # Utilise le base_theme mis à jour (base_size = 19)
    theme(
      # --- MODIFIÉ : Taille du N=... augmentée ---
      plot.subtitle = element_text(size = 27, hjust = 0, margin = margin(b = 4)), # <-- MODIFIÉ (12->17)
      plot.margin = margin(2,2,2,2), 
      legend.position = if(show_legend) "right" else "none"
    )
}


# ===============================================
# --- SECTION DES FONCTIONS D'ASSEMBLAGE ---
# ===============================================

# +++ Fonction pour la barre de traitement ---
plot_treatment_bar <- function(df_treat, year_range) {
  
  x_limits <- c(min(year_range) - 0.5, max(year_range) + 0.5)
  
  ggplot(df_treat) +
    geom_rect(aes(xmin = start, xmax = end, ymin = 0, ymax = 1, fill = policy), 
              color = "white", linewidth = 0.5) +
    # --- MODIFIÉ : Taille augmentée ---
    geom_text(aes(x = (start + end) / 2, y = 0.5, label = policy), 
              color = "white", fontface = "bold", size = 11.5) + # <-- MODIFIÉ (6.5->11.5)
    scale_x_continuous(limits = x_limits, expand = c(0,0)) +
    scale_fill_manual(values = c("DHA-PPQ" = "#1b9e77", "AS-MQ" = "#d95f02")) +
    theme_void() +
    theme(legend.position = "none")
}

# +++ Fonction pour la courbe d'incidence ---
plot_incidence_curve <- function(df_incidence) {
  
  x_limits <- c(min(df_incidence$Year) - 0.5, max(df_incidence$Year) + 0.5)
  
  ggplot(df_incidence, aes(x = Year, y = API_1000)) +
    geom_line(color = "#005a9e", linewidth = 1.2) +
    geom_point(color = "#005a9e", size = 3) +
    # --- MODIFIÉ : Taille des valeurs API augmentée ---
    geom_text(
      aes(label = format(API_1000, decimal.mark = ",")),
      vjust = -0.8, 
      size = 11.5, # <-- MODIFIÉ (6.5->11.5)
      color = "black"
    ) +
    scale_x_continuous(breaks = df_incidence$Year, limits = x_limits, expand = c(0,0)) +
    scale_y_continuous(
      limits = c(0, max(df_incidence$API_1000) * 1.15),
      expand = c(0, 0) 
    ) +
    labs(
      title = "Annual Parasite Incidence (API) per 1,000 population, Cambodia",
      x = NULL, 
      y = "API / 1,000"
    ) +
    # --- MODIFIÉ : Thème de base plus grand ---
    theme_classic(base_size = 26) + # <-- MODIFIÉ (16->21)
    theme(
      plot.title = element_text(hjust = 0.5, size = 23, face = "bold"), # <-- MODIFIÉ (18->23)
      axis.title.y = element_text(face = "bold", size = 19), # <-- MODIFIÉ (14->19)
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
      panel.grid.major.y = element_line(color = "grey85", linetype = "dashed"),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(), 
      axis.ticks.x = element_blank(),
      axis.line.x = element_blank(),
      axis.line.y = element_blank() 
    )
}

# +++ Fonction pour la rangée de labels d'année ---
plot_year_labels_row <- function(years, year_totals_df) {
  
  labels_list <- lapply(years, function(y) {
    total_n <- year_totals_df %>%
      filter(Year == y) %>%
      pull(N_Total)
    
    year_label <- paste0(y, " (N = ", total_n, ")")
    
    # --- MODIFIÉ : Taille augmentée ---
    cowplot::ggdraw() + 
      cowplot::draw_label(year_label, fontface = "bold", size = 30) + # <-- MODIFIÉ (22->28)
      theme(plot.margin = margin(t = 20, b = 5)) 
  })
  
  cowplot::plot_grid(plotlist = labels_list, nrow = 1)
}

# +++ Fonction pour la colonne de labels des marqueurs ---
plot_marker_labels <- function() {
  
  # --- MODIFIÉ : 'label_size' augmenté ---
  label_size <- 35 # <-- MODIFIÉ (16->22)
  l1 <- ggdraw() + draw_label(bquote(bold(italic("k13")~" C580Y (%)")), size = label_size, x = 0.95, hjust = 1)
  l2 <- ggdraw() + draw_label(bquote(bold(italic("k13")~" Y493H (%)")), size = label_size, x = 0.95, hjust = 1)
  l3 <- ggdraw() + draw_label(bquote(bold(italic("k13")~" WT (%)")), size = label_size, x = 0.95, hjust = 1)
  l4 <- ggdraw() + draw_label(bquote(bold(italic("pfpm2")~" multicopy (%)")), size = label_size, x = 0.95, hjust = 1)
  l5 <- ggdraw() + draw_label(bquote(bold(italic("pfmdr1")~" multicopy (%)")), size = label_size, x = 0.95, hjust = 1)
  l6 <- ggdraw() + draw_label(bquote(bold(italic("crt")~" mutant (%)")), size = label_size, x = 0.95, hjust = 1)
  l7 <- ggdraw() + draw_label(
    bquote(atop(bold("KEL1/PLA1 (%)"),
                scriptstyle(italic("pfpm2")~"multicopy +"~italic("k13")~"C580Y"))),
    size = label_size, x = 0.95, hjust = 1)
    
  plot_grid(l1, l2, l3, l4, l5, l6, l7, ncol = 1, align = "v")
}

# +++ Fonction pour la colonne de cartes (simplifiée) ---
plot_year_column <- function(year){ 
  
  d <- provinces_kh %>%
    dplyr::left_join(agg_norm %>% dplyr::filter(Year == year), by = "name_norm")
  
  p1 <- map_one(d, "pct_C580Y",      "k13 C580Y",      show_title = FALSE, show_legend = FALSE)
  p2 <- map_one(d, "pct_Y493H",      "k13 Y493H",      show_title = FALSE, show_legend = FALSE)
  p3 <- map_one(d, "pct_WT",         "k13 WT",         show_title = FALSE, show_legend = FALSE)
  p4 <- map_one(d, "pct_PMII_multi", "pfpm2 multicopy",show_title = FALSE, show_legend = FALSE)
  p5 <- map_one(d, "pct_MDR1_multi", "pfmdr1 multicopy",show_title = FALSE, show_legend = FALSE)
  p6 <- map_one(d, "pct_CRT_mutant", "crt mutant",     show_title = FALSE, show_legend = FALSE)
  p7 <- map_one(d, "pct_KEL1_PLA1",  "KEL1/PLA1",      show_title = FALSE, show_legend = FALSE)
  
  maps_column <- cowplot::plot_grid(p1, p2, p3, p4, p5, p6, p7,
                                    ncol = 1, 
                                    align = "v")
  
  return(maps_column)
}

# +++ Fonction pour la grille de cartes (simplifiée) ---
plot_all_years_grid <- function(years = sort(unique(agg_norm$Year))){
  
  columns <- lapply(years, plot_year_column) 
  
  full_map_grid <- cowplot::plot_grid(
    plotlist = columns,
    nrow = 1 
  )
  
  return(full_map_grid)
}


# +++ Fonction pour la colonne de légendes (avec alignement K13) ---
plot_legend_stack <- function(df_sf_dummy) {
  
  df_sf_dummy <- df_sf_dummy %>%
    mutate(
      pct_C580Y = NA_real_,
      pct_Y493H = NA_real_,
      pct_WT = NA_real_,
      pct_PMII_multi = NA_real_,
      pct_MDR1_multi = NA_real_,
      pct_CRT_mutant = NA_real_,
      pct_KEL1_PLA1 = NA_real_,
      N_count_C580Y = NA_integer_,
      N_count_Y493H = NA_integer_,
      N_count_WT = NA_integer_,
      N_count_PMII_multi = NA_integer_,
      N_count_MDR1_multi = NA_integer_,
      N_count_CRT_mutant = NA_integer_,
      N_count_KEL1_PLA1 = NA_integer_,
      N_k13 = NA_integer_,
      N_pmii = NA_integer_,
      N_mdr1 = NA_integer_,
      N_crt = NA_integer_,
      N_kel1 = NA_integer_
    )
  
  # map_one() utilise maintenant base_theme(base_size=19),
  # donc les légendes seront beaucoup plus grandes.
  p_k13 <- map_one(df_sf_dummy, "pct_C580Y", "K13 (%)", show_legend = TRUE)
  l_k13 <- get_legend(p_k13)
  
  p_pmii <- map_one(df_sf_dummy, "pct_PMII_multi", "pfpm2 multicopy (%)", show_legend = TRUE)
  l_pmii <- get_legend(p_pmii)
  
  p_mdr1 <- map_one(df_sf_dummy, "pct_MDR1_multi", "pfmdr1 multicopy (%)", show_legend = TRUE)
  l_mdr1 <- get_legend(p_mdr1)
  
  p_crt <- map_one(df_sf_dummy, "pct_CRT_mutant", "crt mutant (%)", show_legend = TRUE)
  l_crt <- get_legend(p_crt)
  
  p_kel1 <- map_one(df_sf_dummy, "pct_KEL1_PLA1", "KEL1/PLA1 (%)", show_legend = TRUE)
  l_kel1 <- get_legend(p_kel1)
  
  plot_grid(l_k13, l_pmii, l_mdr1, l_crt, l_kel1, 
            ncol = 1, 
            align = "v",
            rel_heights = c(4, 1.2, 1.2, 1.2, 1.2) 
  )
}


# ===============================================
# --- SECTION D'ASSEMBLAGE FINAL & EXPORT ---
# ===============================================

# +++ Fonction d'assemblage final (logique 5 blocs) +++
save_final_plot <- function(
    out_png = "prevalence_KH_multi_years_FINAL.png",
    out_pdf = "prevalence_KH_multi_years_FINAL.pdf"
){
  
  # --- 0. Données de base ---
  yrs <- sort(unique(agg_norm$Year))
  
  # --- 1. Générer TOUS les composants "bruts" ---
  treatment_bar_plot <- plot_treatment_bar(treatment_data, yrs)
  incidence_plot <- plot_incidence_curve(incidence_data)
  year_label_row <- plot_year_labels_row(yrs, year_totals) 
  map_grid_raw <- plot_all_years_grid(yrs) 
  
  # --- MODIFIÉ : Taille du caption augmentée ---
  global_caption <- ggdraw() + 
    cowplot::draw_label(
      "Grey = missing data",
      fontface = "italic", x = 0.01, hjust = 0, size = 18 # <-- MODIFIÉ (13->18)
    )
  
  labels_col_raw <- plot_marker_labels()
  legends_col_raw <- plot_legend_stack(provinces_kh) 
  
  
  # --- 2. Assemblage en 5 blocs verticaux ---
  
  # [1] Barre, [2] Courbe, [3] Labels Année, [4] Cartes, [5] Caption
  # --- MODIFIÉ : Augmentation de la hauteur des labels d'année ---
  main_heights <- c(0.7, 5, 2.0, 22, 0.5) # <-- MODIFIÉ (1.5 -> 2.0, 0.3 -> 0.5)
  
  
  # --- 2a. Colonne centrale ---
  main_content_col <- plot_grid(
    treatment_bar_plot,
    incidence_plot,
    year_label_row,     # Bloc 3
    map_grid_raw,       # Bloc 4
    global_caption,
    ncol = 1,
    align = "v",
    rel_heights = main_heights
  )
  
  # --- 2b. Colonne de gauche (Labels) ---
  labels_col_padded <- plot_grid(
    NULL,             
    NULL,             
    NULL,             
    labels_col_raw,   # Bloc 4
    NULL,             
    ncol = 1,
    rel_heights = main_heights 
  )
  
  # --- 2c. Colonne de droite (Légendes) ---
  legends_col_padded <- plot_grid(
    NULL,               
    NULL,               
    NULL,               
    legends_col_raw,    # Bloc 4
    NULL,               
    ncol = 1,
    rel_heights = main_heights 
  )
  
  # --- 3. Assemblage horizontal ---
  final_plot <- plot_grid(
    labels_col_padded,  
    main_content_col,   
    legends_col_padded, 
    nrow = 1,
    align = "h", 
    
    # --- MODIFIÉ : Plus de place pour les labels/légendes ---
    rel_widths = c(2.8, length(yrs) + 1.5, 1.5) # <-- MODIFIÉ (1.8/1.2 -> 2.8/1.5)
  )
  
  # --- 4. Exporter ---
  # --- MODIFIÉ : Augmentation de la taille TOTALE de l'image ---
  new_width <- (length(yrs) * 4) + 7  # <-- MODIFIÉ (env. 35 pouces)
  new_height <- 36 # <-- MODIFIÉ (env. 36 pouces)
  
  ggplot2::ggsave(out_png, final_plot,
                  width = new_width, height = new_height,
                  dpi = 300, bg = "white")
  
  ggplot2::ggsave(out_pdf, final_plot,
                  width = new_width, height = new_height,
                  device = cairo_pdf, bg = "white")
  
  print(paste("Figure finale (POLICES POSTER) sauvegardée en :", out_png, "et", out_pdf))
  
  # return(final_plot) 
}

# --- 10) Lancer ---
save_final_plot()
