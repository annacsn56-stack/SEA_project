# ================== IBD networks — Cambodia (PUBLICATION READY - ITALIC FIX) ==================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidygraph)
  library(ggraph)
  library(igraph)
  library(ggplot2)
  library(cowplot)
  library(scales)
  library(patchwork)
  library(tidyr)
})

setwd("/media/annacosson/ANNA_DD/STAGE_M2/Projet_Asie_du_Sud_Est/Article/IBD")

# ---------- PARAMS ----------
file_ibd    <- "/media/annacosson/ANNA_DD/STAGE_M2/Projet_Asie_du_Sud_Est/Article/IBD/hmmIBD/IBD_fract.txt"
file_kel1   <- "metadata_complet_fusionne.csv"

ibd_threshold <- 0.95
year_min <- 2016L; year_max <- 2022L

# Couleurs (clés texte standard pour le mapping)
kel_palette <- c(
  "KEL1/PLA1 (ART-R/PPQ-R)" = "#FF66CC",
  "k13-mut/MDR1-multi (ART-R/MQ-R)" = "#00C080",
  "other" = "#BDBDBD"
)
to_norm <- function(x) tolower(trimws(x))

# ---------- 1) LOAD & PREP ----------
ibd_raw <- read_delim(file_ibd, delim = "\t", show_col_types = FALSE)
kel1_raw <- read_csv(file_kel1, show_col_types = FALSE)

colnames(ibd_raw)[1:3] <- c("sample1","sample2","IBD_value")
ibd_raw <- ibd_raw %>% mutate(sample1 = to_norm(sample1), sample2 = to_norm(sample2))

kel1 <- kel1_raw %>%
  mutate(
    Sample = to_norm(Sample), 
    Year = suppressWarnings(as.integer(trimws(Year))),
    K13_raw = tolower(trimws(K13)), 
    MDR1_raw = tolower(trimws(MDR1)), 
    KEL1_PLA1_proxy_raw = trimws(KEL1_PLA1_proxy)
  ) %>%
  filter(!is.na(Year) & dplyr::between(Year, year_min, year_max)) %>%
  mutate(
    plot_genotype = case_when(
      tolower(KEL1_PLA1_proxy_raw) %in% c("kel1/pla1","kel1pla1","kel1_pla1","kel1") ~ "KEL1/PLA1 (ART-R/PPQ-R)",
      !is.na(K13_raw) & !K13_raw %in% c("wt", "wild-type", "wild type") & MDR1_raw == "multi" ~ "k13-mut/MDR1-multi (ART-R/MQ-R)",
      TRUE ~ "other"
    )
  ) %>%
  filter(!is.na(plot_genotype) & plot_genotype != "")

nodes <- kel1 %>% distinct(Sample, .keep_all = TRUE)
years <- sort(unique(nodes$Year))
edges <- ibd_raw %>%
  filter(IBD_value >= ibd_threshold) %>%
  filter(sample1 %in% nodes$Sample & sample2 %in% nodes$Sample)

if (nrow(nodes) == 0) stop("No eligible nodes.")

# ---------- 2) THEME "PUBLICATION" ----------
theme_pub <- theme_classic(base_size = 15) +
  theme(
    text = element_text(color = "black"),
    axis.text = element_text(color = "black", size = 17),
    axis.title = element_text(face = "bold", size = 18),
    plot.title = element_text(face = "bold", size = 20, hjust = 0),
    plot.subtitle = element_text(size = 14, color = "gray30"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 16),
    legend.text = element_text(size = 14),
    panel.grid.major.y = element_line(color = "gray92", linewidth = 0.4),
    panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 11)
  )

# ---------- 3) PANEL A: NETWORKS ----------
graph_all <- tbl_graph(nodes = nodes %>% rename(name = Sample), edges = edges %>% rename(from = sample1, to = sample2), directed = FALSE)
set.seed(42)
coords <- layout_with_fr(as.igraph(graph_all), niter = 3000, area = vcount(graph_all)^2) %>% as.data.frame()
colnames(coords) <- c("x","y")
coords$name <- graph_all %>% activate(nodes) %>% as_tibble() %>% pull(name)

plot_year_clean <- function(y){
  nodes_y <- nodes %>% filter(Year == y)
  edges_y <- edges %>% filter(sample1 %in% nodes_y$Sample & sample2 %in% nodes_y$Sample)
  lbl <- paste0(y, "\n(N=", nrow(nodes_y), ")")
  
  p <- ggplot() + theme_void() + 
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      plot.margin = margin(2,2,2,2),
      plot.title = element_text(hjust = 0.5, size = 19, face = "bold", margin = margin(b=3))
    ) +
    ggtitle(lbl)
  
  if (nrow(nodes_y) == 0) return(p)
  
  coords_y <- coords %>% filter(name %in% nodes_y$Sample)
  
  if (nrow(edges_y) == 0) {
    p <- p + geom_point(data = coords_y %>% left_join(nodes_y, by = c("name"="Sample")),
                        aes(x=x, y=y, fill=plot_genotype),
                        shape=21, size=2.5, color="white", stroke=0.3) 
  } else {
    g_y <- tbl_graph(nodes = nodes_y %>% rename(name = Sample), edges = edges_y %>% rename(from=sample1, to=sample2), directed=FALSE)
    coords_y <- coords %>% filter(name %in% (g_y %>% activate(nodes) %>% as_tibble() %>% pull(name)))
    p <- ggraph(g_y, layout = "manual", x = coords_y$x, y = coords_y$y) +
      geom_edge_link(color = "black", alpha = 0.6, width = 0.3) +
      geom_node_point(aes(fill = plot_genotype), shape = 21, size = 2.5, color = "black", stroke = 0.3) +
      theme_void() +
      theme(
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        plot.margin = margin(2,2,2,2),
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold", margin = margin(b=3)),
        legend.position = "none"
      ) +
      ggtitle(lbl)
  }
  
  # --- MODIF 1: Ajout des labels avec expression italique ---
  return(p + scale_fill_manual(
    values = kel_palette,
    labels = c(
      "other", 
      expression(paste(italic("k13"), "-mut/MDR1-multi (ART-R/MQ-R)")), 
      "KEL1/PLA1 (ART-R/PPQ-R)"
    )
  ))
}

plots_facets <- lapply(years, plot_year_clean)
panel_A <- wrap_plots(plots_facets, ncol = 3) + plot_layout(guides = "collect") & theme(legend.position = "none")

# ---------- 4) PANEL B: METRICS ----------
metrics <- lapply(years, function(y){
  nodes_y <- nodes %>% filter(Year == y)
  edges_y <- edges %>% filter(sample1 %in% nodes_y$Sample & sample2 %in% nodes_y$Sample)
  if (nrow(nodes_y) == 0) return(NULL)
  g <- graph_from_data_frame(edges_y, directed = FALSE, vertices = nodes_y)
  comps <- components(g)
  sizes <- sort(comps$csize, decreasing = TRUE)
  if (length(sizes) == 0) sizes <- rep(1, nrow(nodes_y))
  n_singletons <- sum(sizes == 1)
  share_clustered <- if (nrow(nodes_y) > 0) (nrow(nodes_y) - n_singletons) / nrow(nodes_y) else 0
  
  frac_kel1 <- mean(nodes_y$plot_genotype == "KEL1/PLA1 (ART-R/PPQ-R)")
  frac_mqr  <- mean(nodes_y$plot_genotype == "k13-mut/MDR1-multi (ART-R/MQ-R)")
  frac_other <- mean(nodes_y$plot_genotype == "other")
  
  tibble(Year = y, n_clusters = length(sizes), share_clustered = share_clustered,
         frac_kel1 = frac_kel1, frac_mqr = frac_mqr, frac_other = frac_other)
}) %>% bind_rows()

metrics_long <- metrics %>%
  select(Year, frac_kel1, frac_mqr, frac_other) %>%
  pivot_longer(cols = c(frac_kel1, frac_mqr, frac_other), names_to = "genotype", values_to = "fraction") %>%
  mutate(genotype = factor(genotype, levels = c("frac_other", "frac_mqr", "frac_kel1"),
                           labels = c("other", "k13-mut/MDR1-multi (ART-R/MQ-R)", "KEL1/PLA1 (ART-R/PPQ-R)")))

v_line <- geom_vline(xintercept = 2018, linetype = "dashed", color = "#D55E00", linewidth = 0.6, alpha = 0.8)

p1 <- ggplot(metrics, aes(Year, n_clusters)) +
  geom_line(color = "#0072B2", linewidth = 0.8) + geom_point(color = "#0072B2", size = 2) + v_line +
  labs(title = "Number of IBD Clusters", y = "Count", x = NULL) + theme_pub

p2 <- ggplot(metrics, aes(Year, share_clustered)) +
  geom_line(color = "#0072B2", linewidth = 0.8) + geom_point(color = "#0072B2", size = 2) + v_line +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(title = "Clustered Samples", subtitle = "(Size > 1)", x = NULL, y = "Proportion") + theme_pub

p3 <- ggplot(metrics_long, aes(x = Year, y = fraction, fill = genotype)) +
  geom_col(position = "stack", width = 0.75, color = "white", size = 0.2) +
  
  # --- MODIF 2: Ajout des labels avec expression italique ---
  scale_fill_manual(
    values = kel_palette,
    labels = c(
      "other", 
      expression(paste(italic("k13"), "-mut/MDR1-multi (ART-R/MQ-R)")), 
      "KEL1/PLA1 (ART-R/PPQ-R)"
    )
  ) +
  
  v_line +
  # --- ICI : ON FORCE TOUTES LES ANNÉES ---
  scale_x_continuous(breaks = 2016:2022) + 
  # ----------------------------------------
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0,1), expand = c(0,0)) +
  labs(title = "Genotype Prevalence", x = NULL, y = "Proportion", fill = "Genotype") +
  theme_pub + theme(legend.position = c(-0.45,0.5),
                    legend.title = element_text(face = "bold", size = 30), # Taille du titre "Genotype"
                    legend.text = element_text(size = 25)
                    )

panel_B_col <- cowplot::plot_grid(p1, p2, p3, ncol = 1, align = "v")



# ---------- 6) FINAL ASSEMBLY (MANUAL LAYOUT / GGDRAW) ----------

# 1. Rangée du haut (50/50)
top_row <- cowplot::plot_grid(panel_A, panel_B_col, ncol = 2, rel_widths = c(1, 1))

# 2. Structure complète avec légende en bas
final_structure <- cowplot::plot_grid(top_row, ncol = 1, rel_heights = c(1, 0.15))

# 3. Calques (Labels A/B + Texte)
final_canvas <- ggdraw(final_structure) +
  
  # Label A
  draw_label("A", x = 0.01, y = 0.985, size = 24, fontface = "bold", hjust = 0) +
  
  # Label B (Milieu)
  draw_label("B", x = 0.51, y = 0.985, size = 24, fontface = "bold", hjust = 0) +
  
  # Phrase
  draw_label(
    "Dashed line (2018): switch from DHA-PPQ to AS-MQ", 
    x = 0.98, y = 0.95, 
    size = 17, fontface = "italic", color = "gray30", hjust = 1
  )

print(final_canvas)

# Sauvegarde
ggsave("IBD_Cambodia_Final.png", final_canvas, width = 21, height = 13, dpi = 350)
ggsave("IBD_Cambodia_Final.pdf", final_canvas, width = 21, height = 13, dpi = 350)

message("✅ Figure sauvegardée : Layout manuel conservé, k13 en italique.")
