# ============================================================
#  04 - MAPAS DE LAS REDES (Figuras 2 y 3)
#  - Descarga la batimetría y carga el mapa base
#  - Dibuja cada red sobre la batimetría (nodos y aristas)
#
#  Salida: fig_red_PM (Magallánico), fig_red_PPA (Rockhopper)
#          y los objetos mapa_base, bathy_PM_df, bathy_PPA_df,
#          que usa 06_mpa_overlap.R
#
#  Requiere 01, 02 y 03 en la misma sesión.
# ============================================================

library(sf)
library(dplyr)
library(ggplot2)
library(ggspatial)
library(marmap)
library(scales)
library(patchwork)

if (!exists("mitre") || !exists("edges_df_PM_OI_a_winter")) {
  stop("Correr primero 01, 02 y 03")
}

options(timeout = 600)

# ------------------------------------------------------------
#  Mapa base (ruta: adaptar)
# ------------------------------------------------------------
f_mapa_base <- "C:/Users/User/Desktop/.../mapa_base.shp"
mapa_base <- st_read(f_mapa_base)

# ------------------------------------------------------------
#  Batimetrías (se descargan de NOAA; keep = TRUE las guarda
#  en el directorio de trabajo)
# ------------------------------------------------------------
bathy_PM  <- getNOAA.bathy(lon1 = -77,  lon2 = -52, lat1 = -65, lat2 = -50,
                           resolution = 1, keep = TRUE)
bathy_PPA <- getNOAA.bathy(lon1 = -100, lon2 = -52, lat1 = -66, lat2 = -48,
                           resolution = 1, keep = TRUE)
bathy_PM_df  <- fortify.bathy(bathy_PM)
bathy_PPA_df <- fortify.bathy(bathy_PPA)

# ------------------------------------------------------------
#  Función: red sobre batimetría
#   - nodes, edges_df: nodos y aristas de la red
#   - color_nodo: color de los nodos
# ------------------------------------------------------------
mapa_red_baty <- function(bathy_df, edges_df, nodes, xlim, ylim, titulo,
                          min_weight = 1, color_nodo = "firebrick") {

  edges_plot <- edges_df %>% filter(weight >= min_weight)

  ggplot() +
    geom_raster(data = bathy_df, aes(x = x, y = y, fill = z),
                interpolate = TRUE, show.legend = FALSE) +
    scale_fill_gradientn(
      colours = c("#748FC4", "#869DC8", "lightsteelblue2", "lightsteelblue1"),
      limits = c(-4000, 0), na.value = "lightsteelblue1", oob = scales::squish) +
    geom_sf(data = mapa_base, fill = "#e8e2d0", color = "black") +
    geom_sf(data = mitre,         fill = NA, color = "black", linewidth = 0.4) +
    geom_sf(data = yaganesI,      fill = NA, color = "black", linewidth = 0.4) +
    geom_sf(data = yaganesII,     fill = NA, color = "black", linewidth = 0.4) +
    geom_sf(data = burdwood,      fill = NA, color = "black", linewidth = 0.4) +
    geom_sf(data = diego_ramirez, fill = NA, color = "black", linewidth = 0.4) +
    geom_segment(data = edges_plot,
                 aes(x = x, y = y, xend = xend, yend = yend, linewidth = weight),
                 alpha = 0.35, color = "gray20") +
    geom_point(data = nodes, aes(x = lon, y = lat, size = n),
               color = color_nodo, alpha = 0.9) +
    scale_linewidth_continuous(range = c(0.15, 1.5), guide = "none") +
    scale_size_continuous(range = c(1.5, 6), guide = "none") +
    annotation_scale(location = "bl") +
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
    labs(title = titulo) +
    theme_minimal() +
    theme(panel.border = element_rect(color = "black", fill = NA),
          axis.title = element_blank(),
          panel.grid = element_blank(),
          axis.ticks = element_line(color = "black", linewidth = 0.4),
          plot.title = element_text(size = 9, face = "bold"))
}

# ------------------------------------------------------------
#  Magallánico (nodos naranja)
# ------------------------------------------------------------
xlim_PM <- c(-76, -55); ylim_PM <- c(-61, -50.3)
col_PM  <- "#FF7052"

red_PM_OI_pre <- mapa_red_baty(bathy_PM_df, edges_df_PM_OI_a,        nodes_PM_OI_a,        xlim_PM, ylim_PM, "Isla Observatorio — pre-molt", color_nodo = col_PM)
red_PM_BF_pre <- mapa_red_baty(bathy_PM_df, edges_df_PM_BF_a,        nodes_PM_BF_a,        xlim_PM, ylim_PM, "Bahía Franklin — pre-molt",    color_nodo = col_PM)
red_PM_OI_win <- mapa_red_baty(bathy_PM_df, edges_df_PM_OI_a_winter, nodes_PM_OI_a_winter, xlim_PM, ylim_PM, "Isla Observatorio — winter",   color_nodo = col_PM)
red_PM_BF_win <- mapa_red_baty(bathy_PM_df, edges_df_PM_BF_a_winter, nodes_PM_BF_a_winter, xlim_PM, ylim_PM, "Bahía Franklin — winter",      color_nodo = col_PM)

fig_red_PM <- (red_PM_BF_pre | red_PM_BF_win) / (red_PM_OI_pre | red_PM_OI_win)
fig_red_PM

# ------------------------------------------------------------
#  Rockhopper (nodos violeta); pre-muda con un recorte más angosto
# ------------------------------------------------------------
xlim_PPA     <- c(-100, -54)
xlim_PPA_pre <- c(-76, -54)
ylim_PPA     <- c(-66, -50)
col_PPA      <- "#A94BB7"

red_PPA_BF_pre <- mapa_red_baty(bathy_PPA_df, edges_df_PPA_BF_a,        nodes_PPA_BF_a,        xlim_PPA_pre, ylim_PPA, "Bahía Franklin — pre-molt", color_nodo = col_PPA)
red_PPA_SJ_pre <- mapa_red_baty(bathy_PPA_df, edges_df_PPA_SJ_a,        nodes_PPA_SJ_a,        xlim_PPA_pre, ylim_PPA, "San Juan — pre-molt",       color_nodo = col_PPA)
red_PPA_BF_win <- mapa_red_baty(bathy_PPA_df, edges_df_PPA_BF_a_winter, nodes_PPA_BF_a_winter, xlim_PPA,     ylim_PPA, "Bahía Franklin — winter",   color_nodo = col_PPA)
red_PPA_SJ_win <- mapa_red_baty(bathy_PPA_df, edges_df_PPA_SJ_a_winter, nodes_PPA_SJ_a_winter, xlim_PPA,     ylim_PPA, "San Juan — winter",         color_nodo = col_PPA)

fig_red_PPA <- (red_PPA_BF_pre | red_PPA_BF_win) / (red_PPA_SJ_pre | red_PPA_SJ_win)
fig_red_PPA
