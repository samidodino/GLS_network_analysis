# ============================================================
#  03b - MAPAS DEL ÍNDICE COMPUESTO DE CENTRALIDAD (opcional)
#  Índice = promedio de in-degree, out-degree y betweenness,
#  cada uno normalizado 0-1 dentro de su red.
#  Valores altos = nodos críticos.
#
#  Salida: all_PM_indice, all_PPA_indice,
#          all_PM_indice_winter, all_PPA_indice_winter
#
#  Requiere 01, 02 y 03 en la misma sesión.
# ============================================================

library(sf)
library(dplyr)
library(ggplot2)
library(patchwork)
library(viridis)
library(rnaturalearth)
library(rnaturalearthdata)

if (!exists("metricas_todas")) stop("Correr primero 01, 02 y 03")

# ------------------------------------------------------------
#  Capas auxiliares (rutas: adaptar)
# ------------------------------------------------------------
f_limite <- "C:/Users/User/Desktop/GLS projectos/Limite_con_Chile/linea_de_limite_FA004.shp"
f_frentes <- "C:/Users/User/Desktop/GLS projectos/Fronts_of_the_Antarctic_Circumpolar_Current_-_GIS_data/shapefile/antarctic_circumpolar_current_fronts.shp"

limite_arg_chile <- st_read(f_limite) %>% filter(grepl("Chile", fna))
polar_front      <- st_read(f_frentes) %>% filter(NAME %in% "Polar Front (PF)")
world            <- ne_countries(scale = "medium", returnclass = "sf")

# ------------------------------------------------------------
#  Función: mapa de una red con nodos coloreados por el índice
#   - capas: lista de sf (AMP, límite, frente) dibujadas como contorno
# ------------------------------------------------------------
norm01 <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

mapa_indice <- function(nodes, metricas, edges_df, titulo, xlim, ylim, capas) {

  nodos <- nodes %>%
    dplyr::select(id, lon, lat) %>%
    left_join(metricas, by = c("id" = "nodo")) %>%
    mutate(indice = (norm01(in_degree) + norm01(out_degree) +
                     norm01(betweenness)) / 3)

  ggplot() +
    geom_sf(data = world, fill = "gray90") +
    geom_segment(data = edges_df,
                 aes(x = x, y = y, xend = xend, yend = yend),
                 color = "gray20", alpha = 0.4) +
    geom_point(data = nodos,
               aes(x = lon, y = lat, size = indice, color = indice),
               alpha = 0.9) +
    geom_text(data = nodos, aes(x = lon, y = lat, label = id),
              size = 3, fontface = "bold", vjust = -1) +
    lapply(capas, function(x) geom_sf(data = x, fill = NA,
                                      color = "black", linewidth = 0.5)) +
    scale_color_viridis_c(name = "Index", limits = c(0, 1)) +
    scale_size_continuous(guide = "none") +
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
    labs(title = titulo) +
    theme_minimal() +
    theme(panel.grid.major = element_line(color = "gray90"),
          axis.title = element_blank(),
          legend.position = "bottom", legend.box = "horizontal")
}

capas_PM  <- list(mitre)
capas_PPA <- list(mitre, yaganesI, yaganesII, burdwood, diego_ramirez,
                  limite_arg_chile, polar_front)

# ------------------------------------------------------------
#  Pre-muda
# ------------------------------------------------------------
xlim_PM <- c(-70, -60); ylim_PM <- c(-58, -53)
xlim_PPA <- c(-75, -50); ylim_PPA <- c(-65, -50)

all_PM_indice <-
  mapa_indice(nodes_PM_OI_a, metricas_PM_OI, edges_df_PM_OI_a,
              "Magellanic – Observatorio (pre-molt)", xlim_PM, ylim_PM, capas_PM) +
  mapa_indice(nodes_PM_BF_a, metricas_PM_BF, edges_df_PM_BF_a,
              "Magellanic – Franklin (pre-molt)", xlim_PM, ylim_PM, capas_PM)

all_PPA_indice <-
  mapa_indice(nodes_PPA_BF_a, metricas_PPA_BF, edges_df_PPA_BF_a,
              "Rockhoppers – Franklin (pre-molt)", xlim_PPA, ylim_PPA, capas_PPA) +
  mapa_indice(nodes_PPA_SJ_a, metricas_PPA_SJ, edges_df_PPA_SJ_a,
              "Rockhoppers – San Juan (pre-molt)", xlim_PPA, ylim_PPA, capas_PPA)

# ------------------------------------------------------------
#  Invierno
# ------------------------------------------------------------
xlim_PM_w <- c(-70, -63); ylim_PM_w <- c(-56, -52)
xlim_PPA_w <- c(-100, -52); ylim_PPA_w <- c(-65, -50)

all_PM_indice_winter <-
  mapa_indice(nodes_PM_OI_a_winter, metricas_PM_OI_winter, edges_df_PM_OI_a_winter,
              "Magellanic – Observatorio (winter)", xlim_PM_w, ylim_PM_w, capas_PM) +
  mapa_indice(nodes_PM_BF_a_winter, metricas_PM_BF_winter, edges_df_PM_BF_a_winter,
              "Magellanic – Franklin (winter)", xlim_PM_w, ylim_PM_w, capas_PM)

all_PPA_indice_winter <-
  mapa_indice(nodes_PPA_BF_a_winter, metricas_PPA_BF_winter, edges_df_PPA_BF_a_winter,
              "Rockhoppers – Franklin (winter)", xlim_PPA_w, ylim_PPA_w, capas_PPA) +
  mapa_indice(nodes_PPA_SJ_a_winter, metricas_PPA_SJ_winter, edges_df_PPA_SJ_a_winter,
              "Rockhoppers – San Juan (winter)", xlim_PPA_w, ylim_PPA_w, capas_PPA)

all_PM_indice
all_PPA_indice
all_PM_indice_winter
all_PPA_indice_winter
