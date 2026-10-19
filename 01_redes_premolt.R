# ============================================================
#  01 - REDES DE MOVIMIENTO, PRE-MUDA
#  Construye una red por especie-colonia (4 redes):
#    1. Clustering jerárquico de las posiciones (nodos candidatos)
#    2. Nodos = clusters usados por >= 50% de los individuos
#    3. Aristas = transiciones entre nodos (una por individuo)
#
#  Salida, para cada red (PM_OI, PM_BF, PPA_BF, PPA_SJ):
#    g_<red>_a          grafo igraph dirigido y ponderado
#    nodes_<red>_a      nodos (sf, con id, n, lon, lat)
#    edges_df_<red>_a   aristas con coordenadas (para graficar)
#
#  Ejecutar en una sesión nueva de R. Define construir_red(),
#  que también usa 02_redes_invierno.R.
# ============================================================

library(sf)
library(sp)
library(igraph)
library(dplyr)
library(geosphere)

# ------------------------------------------------------------
#  Función: construye la red de una especie-colonia
#   - datos: data.frame con id, dtime, lon, lat
#   - d: distancia de corte del clustering (m)
#   - prop_individuos: fracción mínima de individuos por nodo
# ------------------------------------------------------------
construir_red <- function(datos, d, prop_individuos = 0.5) {

  # clustering jerárquico (complete linkage) sobre distancias geodésicas
  xy <- SpatialPointsDataFrame(
    matrix(c(datos$lon, datos$lat), ncol = 2),
    data = datos,
    proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
  hc <- hclust(as.dist(distm(xy)), method = "complete")
  xy$clust    <- cutree(hc, h = d)
  datos$clust <- as.character(xy$clust)
  xy_sf <- st_as_sf(xy)

  # clusters válidos: usados por >= prop_individuos de los individuos
  umbral <- ceiling(prop_individuos * length(unique(xy_sf$id)))
  clusters_validos <- xy_sf %>%
    st_drop_geometry() %>%
    group_by(clust) %>%
    summarise(n_ind = n_distinct(id), .groups = "drop") %>%
    filter(n_ind >= umbral) %>%
    pull(clust)

  # nodos: centroide y número de posiciones de cada cluster válido
  nodos_resumen <- xy_sf %>%
    filter(clust %in% clusters_validos) %>%
    group_by(clust) %>%
    summarise(geometry = st_centroid(st_union(geometry)),
              n = n(), .groups = "drop") %>%
    mutate(id = as.character(clust))

  # aristas: transiciones entre clusters válidos, una por individuo;
  # el peso es el número de individuos que usan la transición
  transiciones <- datos %>%
    filter(clust %in% clusters_validos) %>%
    arrange(id, dtime) %>%
    group_by(id) %>%
    mutate(clust_from = lag(clust), clust_to = clust) %>%
    filter(!is.na(clust_from), clust_from != clust_to) %>%
    ungroup() %>%
    distinct(id, clust_from, clust_to) %>%
    count(clust_from, clust_to, name = "weight") %>%
    filter(clust_from %in% clusters_validos,
           clust_to   %in% clusters_validos)

  g <- graph_from_data_frame(transiciones, directed = TRUE,
                             vertices = nodos_resumen)

  nodes <- nodos_resumen %>%
    mutate(lon = st_coordinates(.)[, 1],
           lat = st_coordinates(.)[, 2])

  edges_df <- igraph::as_data_frame(g, what = "edges") %>%
    left_join(nodes, by = c("from" = "id")) %>%
    rename(x = lon, y = lat) %>%
    left_join(nodes, by = c("to" = "id")) %>%
    rename(xend = lon, yend = lat)

  list(g = g, nodes = nodes, edges_df = edges_df, transiciones = transiciones)
}

# ------------------------------------------------------------
#  Datos de pre-muda (posiciones GLS agrupadas por colonia)
# ------------------------------------------------------------
PM_BF_premolt  <- read.csv("PM_Franklin/Results/datos_agrupados_PM_franklin_pre-molt.csv")
PM_OI_premolt  <- read.csv("PM_Obser/Results/datos_agrupados_Obser_pre-molt.csv")
PPA_BF_premolt <- read.csv("PPA_Franklin/Results/overlaps_with_best_track/datos_agrupados_pre-molt_withSex.csv")
PPA_SJ_premolt <- read.csv("PPA_San_Juan/Results/datos_agrupados_san_juan_pre-molt.csv")

# ------------------------------------------------------------
#  Construir las 4 redes
#   Magallánico: d = 30 km (error del geolocalizador; evita crear
#                nodos más finos que la resolución del GLS)
#   Rockhopper:  d = 100 km (error espacial mayor en esta especie)
# ------------------------------------------------------------
redes_premolt <- list(
  PM_OI  = construir_red(PM_OI_premolt,  d = 30000),
  PM_BF  = construir_red(PM_BF_premolt,  d = 30000),
  PPA_BF = construir_red(PPA_BF_premolt, d = 100000),
  PPA_SJ = construir_red(PPA_SJ_premolt, d = 100000))

# objetos con los nombres que usan los scripts siguientes
for (r in names(redes_premolt)) {
  assign(paste0("g_",        r, "_a"), redes_premolt[[r]]$g)
  assign(paste0("nodes_",    r, "_a"), redes_premolt[[r]]$nodes)
  assign(paste0("edges_df_", r, "_a"), redes_premolt[[r]]$edges_df)
}

# control: nodos y aristas por red 
sapply(redes_premolt, function(x) c(nodos = vcount(x$g), aristas = ecount(x$g)))

