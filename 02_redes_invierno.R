# ============================================================
#  02 - REDES DE MOVIMIENTO, INVIERNO
#  Construye una red por especie-colonia (4 redes) con el mismo
#  procedimiento que la pre-muda (ver construir_red() en
#  01_redes_premolt.R). Al final arma la tabla de nodos y
#  aristas de las 8 redes.
#
#  Salida, para cada red (PM_OI, PM_BF, PPA_BF, PPA_SJ):
#    g_<red>_a_winter          grafo igraph dirigido y ponderado
#    nodes_<red>_a_winter      nodos (sf, con id, n, lon, lat)
#    edges_df_<red>_a_winter   aristas con coordenadas (para graficar)
#
#  Requiere haber corrido 01_redes_premolt.R en la misma sesión.
# ============================================================

library(sf)
library(igraph)
library(dplyr)

if (!exists("construir_red")) stop("Correr primero 01_redes_premolt.R")

# ------------------------------------------------------------
#  Datos de invierno (posiciones GLS agrupadas por colonia)
# ------------------------------------------------------------
PM_BF_winter  <- read.csv("PM_Franklin/Results/datos_agrupados_PM_franklin_winter.csv")
PM_OI_winter  <- read.csv("PM_Obser/Results/datos_agrupados_Obser_winter.csv")
PPA_BF_winter <- read.csv("PPA_Franklin/Results/datos_agrupados_PPA_franklin_winter.csv")
PPA_SJ_winter <- read.csv("PPA_San_Juan/Results/datos_agrupados_san_juan_winter.csv")

# ------------------------------------------------------------
#  Construir las 4 redes (d = 100 km en todas)
# ------------------------------------------------------------
redes_winter <- list(
  PM_OI  = construir_red(PM_OI_winter,  d = 100000),
  PM_BF  = construir_red(PM_BF_winter,  d = 100000),
  PPA_BF = construir_red(PPA_BF_winter, d = 100000),
  PPA_SJ = construir_red(PPA_SJ_winter, d = 100000))

# objetos con los nombres que usan los scripts siguientes
for (r in names(redes_winter)) {
  assign(paste0("g_",        r, "_a_winter"), redes_winter[[r]]$g)
  assign(paste0("nodes_",    r, "_a_winter"), redes_winter[[r]]$nodes)
  assign(paste0("edges_df_", r, "_a_winter"), redes_winter[[r]]$edges_df)
}

# control: nodos y aristas por red 
sapply(redes_winter, function(x) c(nodos = vcount(x$g), aristas = ecount(x$g)))

# ------------------------------------------------------------
#  Tabla de las 8 redes: número de nodos y aristas
# ------------------------------------------------------------
redes <- list(
  list("Magellanic",          "Bahía Franklin",         "pre-molt", g_PM_BF_a),
  list("Magellanic",          "Bahía Franklin",         "winter",   g_PM_BF_a_winter),
  list("Magellanic",          "Isla Observatorio",      "pre-molt", g_PM_OI_a),
  list("Magellanic",          "Isla Observatorio",      "winter",   g_PM_OI_a_winter),
  list("Southern Rockhopper", "Bahía Franklin",         "pre-molt", g_PPA_BF_a),
  list("Southern Rockhopper", "Bahía Franklin",         "winter",   g_PPA_BF_a_winter),
  list("Southern Rockhopper", "San Juan de Salvamento", "pre-molt", g_PPA_SJ_a),
  list("Southern Rockhopper", "San Juan de Salvamento", "winter",   g_PPA_SJ_a_winter))

tabla_redes <- bind_rows(lapply(redes, function(x) tibble(
  Penguin   = x[[1]],
  Colony    = x[[2]],
  Stage     = x[[3]],
  `N nodes` = vcount(x[[4]]),
  `N edges` = ecount(x[[4]]))))

tabla_redes
