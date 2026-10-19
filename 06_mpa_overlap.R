# ============================================================
#  06 - SUPERPOSICIÓN DE NODOS CON AMP
#  - Clasifica cada nodo como protegido (centroide dentro de AMP)
#  - Tabla de cobertura y conectividad dentro/fuera de AMP (Tabla 2)
#  - Mapas de validación nodos dentro/fuera de AMP (S3, S4)
#
#  Requiere en la misma sesión (scripts previos):
#   - 01 y 02: nodes_<red>_a y nodes_<red>_a_winter
#   - 03: metricas_<red>[_winter] y los polígonos de AMP
#         (mitre, yaganesI, yaganesII, burdwood, diego_ramirez)
#   - 04: bathy_PM_df, bathy_PPA_df, mapa_base
#  No depende de 05 (robustez).
# ============================================================

library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggspatial)
library(patchwork)

if (!exists("bathy_PM_df") || !exists("metricas_todas")) {
  stop("Correr primero 01, 02, 03 y 04")
}

# ------------------------------------------------------------
#  1. Unir las AMP en un solo polígono
# ------------------------------------------------------------
amps <- bind_rows(
  st_transform(mitre,         4326) %>% transmute(amp = "Mitre"),
  st_transform(yaganesI,      4326) %>% transmute(amp = "Yaganes I"),
  st_transform(yaganesII,     4326) %>% transmute(amp = "Yaganes II"),
  st_transform(burdwood,      4326) %>% transmute(amp = "Burdwood"),
  st_transform(diego_ramirez, 4326) %>% transmute(amp = "Diego Ramirez"))
amps <- st_make_valid(amps)       # evita fallas por vértices duplicados
amps_union <- st_union(amps)

# ------------------------------------------------------------
#  2. Clasificar nodos: protegido = centroide dentro de alguna AMP
#     (se trabaja en UTM 19S para la intersección)
# ------------------------------------------------------------
clasificar_nodos_amp <- function(nodes_sf, amps_union, crs_metrico = 32719) {
  nodos  <- st_transform(st_as_sf(nodes_sf), crs_metrico)
  amps_m <- st_transform(amps_union, crs_metrico)
  nodos$protegido <- st_intersects(nodos, amps_m, sparse = FALSE)[, 1]
  list(n_nodes_total = nrow(nodos),
       n_protegidos  = sum(nodos$protegido),
       nodes_sf      = nodos)
}

res_amp <- lapply(list(
  PM_OI_premolt  = nodes_PM_OI_a,
  PM_BF_premolt  = nodes_PM_BF_a,
  PPA_BF_premolt = nodes_PPA_BF_a,
  PPA_SJ_premolt = nodes_PPA_SJ_a,
  PM_OI_winter   = nodes_PM_OI_a_winter,
  PM_BF_winter   = nodes_PM_BF_a_winter,
  PPA_BF_winter  = nodes_PPA_BF_a_winter,
  PPA_SJ_winter  = nodes_PPA_SJ_a_winter),
  clasificar_nodos_amp, amps_union = amps_union)

# ------------------------------------------------------------
#  3. Tabla 2: cobertura y métricas dentro (in) vs fuera (out) de AMP
# ------------------------------------------------------------
metricas <- list(
  PM_OI_premolt  = metricas_PM_OI,
  PM_BF_premolt  = metricas_PM_BF,
  PPA_BF_premolt = metricas_PPA_BF,
  PPA_SJ_premolt = metricas_PPA_SJ,
  PM_OI_winter   = metricas_PM_OI_winter,
  PM_BF_winter   = metricas_PM_BF_winter,
  PPA_BF_winter  = metricas_PPA_BF_winter,
  PPA_SJ_winter  = metricas_PPA_SJ_winter)

# cobertura: una fila por red
cobertura <- tibble(
  red         = names(res_amp),
  N_nodes     = sapply(res_amp, `[[`, "n_nodes_total"),
  N_protected = sapply(res_amp, `[[`, "n_protegidos")) %>%
  mutate(Pct_protected = round(100 * N_protected / N_nodes))

# in-degree, out-degree y betweenness medios dentro y fuera de AMP
resumen_metricas <- bind_rows(lapply(names(res_amp), function(r) {
  res_amp[[r]]$nodes_sf %>%
    st_drop_geometry() %>%
    mutate(id = as.character(id)) %>%
    left_join(mutate(metricas[[r]], nodo = as.character(nodo)),
              by = c("id" = "nodo")) %>%
    group_by(protegido) %>%
    summarise(in_deg_media  = mean(in_degree,   na.rm = TRUE),
              out_deg_media = mean(out_degree,  na.rm = TRUE),
              betw_media    = mean(betweenness, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(red = r)
}))

metricas_wide <- resumen_metricas %>%
  mutate(estado = ifelse(protegido, "in", "out")) %>%
  select(red, estado, in_deg_media, out_deg_media, betw_media) %>%
  pivot_wider(names_from = estado,
              values_from = c(in_deg_media, out_deg_media, betw_media)) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

tabla_X <- cobertura %>%
  left_join(metricas_wide, by = "red") %>%
  mutate(
    Species = ifelse(grepl("^PM", red), "Magellanic", "Southern Rockhopper"),
    Colony  = case_when(grepl("_OI_", red) ~ "Isla Observatorio",
                        grepl("_BF_", red) ~ "Bahía Franklin",
                        grepl("_SJ_", red) ~ "San Juan de Salvamento"),
    Stage   = ifelse(grepl("premolt", red), "Pre-molt", "Winter")) %>%
  select(Species, Colony, Stage, N_nodes, N_protected, Pct_protected,
         InDeg_in  = in_deg_media_in,  InDeg_out  = in_deg_media_out,
         OutDeg_in = out_deg_media_in, OutDeg_out = out_deg_media_out,
         Betw_in   = betw_media_in,    Betw_out   = betw_media_out) %>%
  arrange(Species, Colony, factor(Stage, levels = c("Pre-molt", "Winter")))

# guiones en las redes sin nodos protegidos
tabla_X_print <- tabla_X %>%
  mutate(across(where(is.numeric), ~ ifelse(is.na(.x), "—", as.character(.x))))
print(as.data.frame(tabla_X_print), row.names = FALSE)

# write.csv(tabla_X_print, "tabla_X_amp.csv", row.names = FALSE)

# ------------------------------------------------------------
#  4. Mapas de validación: nodos dentro/fuera de AMP
# ------------------------------------------------------------
mapa_validacion_amp <- function(bathy_df, res, xlim, ylim, titulo) {
  nodos <- st_transform(res$nodes_sf, 4326)
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
    geom_sf(data = nodos, aes(color = protegido), size = 3) +
    scale_color_manual(name = "Within MPA",
                       values = c(`TRUE` = "firebrick", `FALSE` = "gray20"),
                       labels = c(`TRUE` = "Inside", `FALSE` = "Outside")) +
    annotation_scale(location = "bl") +
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
    labs(title = titulo) +
    theme_minimal() +
    theme(panel.border = element_rect(color = "black", fill = NA),
          axis.title = element_blank(), panel.grid = element_blank(),
          axis.ticks = element_line(color = "black", linewidth = 0.4),
          plot.title = element_text(size = 9, face = "bold"),
          legend.position = "bottom")
}

# Magallánico
m1 <- mapa_validacion_amp(bathy_PM_df, res_amp$PM_OI_premolt, c(-76, -55), c(-61, -50.3), "PM Observatorio (pre-molt)")
m2 <- mapa_validacion_amp(bathy_PM_df, res_amp$PM_BF_premolt, c(-76, -55), c(-61, -50.3), "PM Franklin (pre-molt)")
m5 <- mapa_validacion_amp(bathy_PM_df, res_amp$PM_OI_winter,  c(-76, -55), c(-61, -50.3), "PM Observatorio (winter)")
m6 <- mapa_validacion_amp(bathy_PM_df, res_amp$PM_BF_winter,  c(-76, -55), c(-61, -50.3), "PM Franklin (winter)")

# Rockhopper
m3 <- mapa_validacion_amp(bathy_PPA_df, res_amp$PPA_BF_premolt, c(-76, -54),  c(-65, -53), "PPA Franklin (pre-molt)")
m4 <- mapa_validacion_amp(bathy_PPA_df, res_amp$PPA_SJ_premolt, c(-76, -54),  c(-65, -53), "PPA San Juan (pre-molt)")
m7 <- mapa_validacion_amp(bathy_PPA_df, res_amp$PPA_BF_winter,  c(-100, -54), c(-66, -50), "PPA Franklin (winter)")
m8 <- mapa_validacion_amp(bathy_PPA_df, res_amp$PPA_SJ_winter,  c(-100, -54), c(-66, -50), "PPA San Juan (winter)")

# Franklin arriba, Observatorio abajo; pre-molt a la izquierda, winter a la derecha
(m2 + m6) / (m1 + m5)

# Franklin arriba, San Juan abajo
(m3 + m7) / (m4 + m8)

