# ============================================================
#  03 - MÉTRICAS DE CONECTIVIDAD POR NODO + POLÍGONOS DE AMP
#  - Métricas de cada nodo de las 8 redes
#  - Carga de los polígonos de AMP (los usan 04 y 06)
#
#  Salida:
#    metricas_<red>[_winter]   una tabla por red (columna 'nodo' = id del nodo)
#    metricas_todas            las 8 tablas juntas
#    mitre, yaganesI, yaganesII, burdwood, diego_ramirez
#
#  Requiere haber corrido 01 y 02 en la misma sesión.
# ============================================================

library(sf)
library(igraph)
library(dplyr)

if (!exists("g_PM_OI_a_winter")) stop("Correr primero 01_redes_premolt.R y 02_redes_invierno.R")

# ------------------------------------------------------------
#  Rutas de los polígonos de AMP (adaptar)
# ------------------------------------------------------------
dir_amp  <- "C:/Users/User/Desktop/GLS projectos/MarineAreas"
f_mitre  <- "C:/Users/User/Desktop/GLS projectos/Poligonos_Áreas_Protegidas-20251118T185019Z-1-001/Poligonos_Áreas_Protegidas/Poligono_Peninsula_Mitre_Merge_P07F2.shp"

mitre         <- st_read(f_mitre)
yaganesI      <- st_read(file.path(dir_amp, "Yaganes I.shp"))
yaganesII     <- st_read(file.path(dir_amp, "Yaganes II.shp"))
burdwood      <- st_read(file.path(dir_amp, "BB_I-II/BB_I_II_limite_zonificacion.shp"))
diego_ramirez <- st_read(file.path(dir_amp, "Diego_Ramirez/WDPA_WDOECM_Sep2023_Public_555643508_shp-polygons.shp"))

# ------------------------------------------------------------
#  Función: métricas de conectividad de cada nodo de un grafo
# ------------------------------------------------------------
calcular_metricas_conectividad <- function(grafo, etapa, spp, colonia) {

  # asegurar que la red tenga pesos
  if (is.null(E(grafo)$weight)) E(grafo)$weight <- 1

  # degree (no ponderado) y strength (ponderado), entrada y salida
  in_degree    <- degree(grafo, mode = "in")
  out_degree   <- degree(grafo, mode = "out")
  in_strength  <- strength(grafo, mode = "in",  weights = E(grafo)$weight)
  out_strength <- strength(grafo, mode = "out", weights = E(grafo)$weight)

  # betweenness y closeness ponderados: distancia = 1/peso
  # (transición frecuente = nodos más cercanos)
  betw  <- betweenness(grafo, directed = TRUE,
                       weights = 1 / E(grafo)$weight, normalized = TRUE)
  close <- closeness(grafo, mode = "out",
                     weights = 1 / E(grafo)$weight, normalized = TRUE)

  # centralidad de autovector
  eig <- eigen_centrality(grafo, directed = TRUE,
                          weights = E(grafo)$weight)$vector

  data.frame(
    nodo = names(in_strength),
    in_degree, out_degree, in_strength, out_strength,
    betweenness = betw, closeness = close, eigen = eig,
    etapa = etapa, especie = spp, colonia = colonia)
}

# ------------------------------------------------------------
#  Aplicar a las 8 redes
# ------------------------------------------------------------
# pre-muda
metricas_PM_OI  <- calcular_metricas_conectividad(g_PM_OI_a,  "premolt", "PM",  "Observatorio")
metricas_PM_BF  <- calcular_metricas_conectividad(g_PM_BF_a,  "premolt", "PM",  "Franklin")
metricas_PPA_BF <- calcular_metricas_conectividad(g_PPA_BF_a, "premolt", "PPA", "Franklin")
metricas_PPA_SJ <- calcular_metricas_conectividad(g_PPA_SJ_a, "premolt", "PPA", "San_Juan")

# invierno
metricas_PM_OI_winter  <- calcular_metricas_conectividad(g_PM_OI_a_winter,  "winter", "PM",  "Observatorio")
metricas_PM_BF_winter  <- calcular_metricas_conectividad(g_PM_BF_a_winter,  "winter", "PM",  "Franklin")
metricas_PPA_BF_winter <- calcular_metricas_conectividad(g_PPA_BF_a_winter, "winter", "PPA", "Franklin")
metricas_PPA_SJ_winter <- calcular_metricas_conectividad(g_PPA_SJ_a_winter, "winter", "PPA", "San_Juan")

metricas_todas <- bind_rows(
  metricas_PM_OI, metricas_PM_BF, metricas_PPA_BF, metricas_PPA_SJ,
  metricas_PM_OI_winter, metricas_PM_BF_winter,
  metricas_PPA_BF_winter, metricas_PPA_SJ_winter)
