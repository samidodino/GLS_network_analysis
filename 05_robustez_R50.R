# ============================================================
#  07 - ROBUSTEZ DE LAS REDES FRENTE A EXTINCIONES
#  - Convierte las redes de igraph a objetos compatibles con network
#  - Simula extinciones aleatorias y por criterios dirigidos
#  - Calcula in-degree, out-degree, betweenness e índice combinado
#  - Compara R50 entre especies, colonias y etapas del ciclo anual
#  - Genera tablas resumen y figuras de robustez de las redes
#
#  Salidas principales:
#  - Tablas R50 por red y criterio de extinción
#  - tabla_R50_total_1, tabla_R50_total_2 y tabla_R50_IS050
#  - Figuras comparativas de robustez (R50)
#
#  Requiere que los objetos de red g_* estén cargados previamente.

# ============================================================

library(igraph)
library(network)
library(NetworkExtinction)
library(ggplot2)
library(dplyr)
library(ggh4x)

## 1. Preparación y análisis inicial: Magallánico, Observatorio, premolt ----
#g_simp <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
g_simp <- igraph::simplify(g_PM_OI_a, remove.multiple = TRUE, remove.loops = TRUE)
el <- igraph::as_edgelist(g_simp, names = TRUE)
net <- network::as.network(el, directed = TRUE)
network::is.network(net)

## Conversión a matriz de adyacencia para NetworkExtinction ----
matrix <- as.matrix(net)
colnames(matrix) <- NULL
rownames(matrix) <- NULL
#colnames(matrix) <- c(1:379)
Mnet <- as.network(matrix)

## Extinciones simuladas y criterios dirigidos ----
res <- SimulateExtinctions(Network = Mnet, Method = "Mostconnected")
res$sims

ExtinctionPlot(History = res[[1]], Variable = "AccSecExt")

# ExtinctionOrder
#Calcular in-degree
in_deg <- degree(g_simp, mode = "in")
#Ordenar de mayor a menor el in_degree:
orden_in_degree <- order(in_deg, decreasing = TRUE)

res_in <- ExtinctionOrder(
  Network = Mnet,
  Order = orden_in_degree,
  NetworkType = "Trophic",
  IS = 0.5)

#Calcular out-degree
out_deg <- degree(g_simp, mode = "out")
#Ordenar de mayor a menor el out_degree:
orden_out_degree <- order(out_deg, decreasing = TRUE)

res_out <- ExtinctionOrder(
  Network = Mnet,
  Order = orden_out_degree,
  NetworkType = "Trophic",
  IS = 0.5)

#Calcular betweenness
betw <- betweenness(
  g_simp,
  directed = TRUE,
  weights = 1 / E(g_simp)$weight  # más peso = más cercanía. peso alto = transición frecuente = más conectividad.
)
#Ordenar de mayor a menor el betweenness:
orden_betw <- order(betw, decreasing = TRUE)

res_betw <- ExtinctionOrder(
  Network = Mnet,
  Order = orden_betw,
  NetworkType = "Trophic",
  IS = 0.5)

# Calcular indice conjunto (centralidad de flujo?)
#Funcion auxiliar para estandarizar
norm01 <- function(x) {
  (x - min(x, na.rm = TRUE)) /
    (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

in_deg_norm <- norm01(in_deg)
out_deg_norm <- norm01(out_deg)
betw_norm <- norm01(betw)

indice <- (in_deg_norm + out_deg_norm + betw_norm) / 3
#Ordenar de mayor a menor indice:
orden_indice <- order(indice, decreasing = TRUE)

res_indice <- ExtinctionOrder(
  Network = Mnet,
  Order = orden_indice,
  NetworkType = "Trophic",
  IS = 0.5)

tabla_R50_PM_OI_premolt <- data.frame(
  criterio = c("In-degree","Out-degree",
    "Betweenness","Índice combinado"),
  R50 = c(res_in$R50, res_out$R50, res_betw$R50,res_indice$R50),
  R100 = c(res_in$R100,res_out$R100,res_betw$R100,res_indice$R100))

tabla_R50_PM_OI_premolt

tabla_metricas <- data.frame(
  nodo = seq_along(in_deg),
  in_degree = in_deg,
  out_degree = out_deg,
  betweenness = betw,
  in_deg_n = in_deg_norm,
  out_deg_n = out_deg_norm,
  betw_n = betw_norm,
  indice = indice)

data.frame(
  in_d = orden_in_degree,
  out_d = orden_out_degree,
  betw = orden_betw,
  indice = orden_indice)

#Random extinctions---
set.seed(707)
random_ext_PM_OI_premolt <- RandomExtinctions(Network= Mnet, nsim= 100,SimNum = 7)

## 2. Magallánico: Franklin, premolt ----
g_simp_PM_BF_premolt <- igraph::simplify(g_PM_BF_a, remove.multiple = TRUE, remove.loops = TRUE)
el_PM_BF_premolt <- igraph::as_edgelist(g_simp_PM_BF_premolt, names = TRUE)
net_PM_BF_premolt <- network::as.network(el_PM_BF_premolt, directed = TRUE)
network::is.network(net_PM_BF_premolt)

matrix_PM_BF_premolt <- as.matrix(net_PM_BF_premolt)
colnames(matrix_PM_BF_premolt) <- NULL
rownames(matrix_PM_BF_premolt) <- NULL

Mnet_PM_BF_premolt <- as.network(matrix_PM_BF_premolt)

res_PM_BF_premolt <- SimulateExtinctions(Network = Mnet_PM_BF_premolt, Method = "Mostconnected")

#in-degree
in_deg_PM_BF_premolt <- degree(g_simp_PM_BF_premolt, mode = "in")
#Ordenar de mayor a menor el in_degree:
orden_in_deg_PM_BF_premolt <- order(in_deg_PM_BF_premolt, decreasing = TRUE)

res_in_PM_BF_premolt <- ExtinctionOrder(
  Network = Mnet_PM_BF_premolt,
  Order = orden_in_deg_PM_BF_premolt,
  NetworkType = "Trophic",
  IS = 0.5)

#Out-degree
out_deg_PM_BF_premolt <- degree(g_simp_PM_BF_premolt, mode = "out")
#Ordenar de mayor a menor el out_degree:
orden_out_deg_PM_BF_premolt <- order(out_deg_PM_BF_premolt, decreasing = TRUE)

res_out_PM_BF_premolt <- ExtinctionOrder(
  Network = Mnet_PM_BF_premolt,
  Order = orden_out_deg_PM_BF_premolt,
  NetworkType = "Trophic",
  IS = 0.5)

#Calcular betweenness
betw_PM_BF_premolt <- betweenness(
  g_simp_PM_BF_premolt,
  directed = TRUE,
  weights = 1 / E(g_simp_PM_BF_premolt)$weight  # más peso = más cercanía. peso alto = transición frecuente = más conectividad.
)
#Ordenar de mayor a menor el betweenness:
orden_betw_PM_BF_premolt <- order(betw_PM_BF_premolt, decreasing = TRUE)

res_betw_PM_BF_premolt <- ExtinctionOrder(
  Network = Mnet_PM_BF_premolt,
  Order = orden_betw_PM_BF_premolt,
  NetworkType = "Trophic",
  IS = 0.5)

# Calcular indice conjunto (centralidad de flujo?)
#Funcion auxiliar para estandarizar
norm01 <- function(x) {
  (x - min(x, na.rm = TRUE)) /
    (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

in_deg_norm_PM_BF_premolt <- norm01(in_deg_PM_BF_premolt)
out_deg_norm_PM_BF_premolt <- norm01(out_deg_PM_BF_premolt)
betw_norm_PM_BF_premolt <- norm01(betw_PM_BF_premolt)

indice_PM_BF_premolt <- (in_deg_norm_PM_BF_premolt + out_deg_norm_PM_BF_premolt + betw_norm_PM_BF_premolt) / 3
#Ordenar de mayor a menor indice:
orden_indice_PM_BF_premolt <- order(indice_PM_BF_premolt, decreasing = TRUE)

res_indice_PM_BF_premolt <- ExtinctionOrder(
  Network = Mnet_PM_BF_premolt,
  Order = orden_indice_PM_BF_premolt,
  NetworkType = "Trophic",
  IS = 0.5)

tabla_R50_PM_BF_premolt <- data.frame(
  criterio = c(
    "In-degree",
    "Out-degree",
    "Betweenness",
    "Índice combinado"),
  R50 = c(
    res_in_PM_BF_premolt$R50,
    res_out_PM_BF_premolt$R50,
    res_betw_PM_BF_premolt$R50,
    res_indice_PM_BF_premolt$R50),
  R100 = c(
    res_in_PM_BF_premolt$R100,
    res_out_PM_BF_premolt$R100,
    res_betw_PM_BF_premolt$R100,
    res_indice_PM_BF_premolt$R100))

tabla_R50_PM_BF_premolt

#Random extinctions---
set.seed(707)
random_ext_PM_BF_premolt <- RandomExtinctions(Network= Mnet_PM_BF_premolt, nsim= 100,SimNum = 8)

## 3. Magallánico: Observatorio, winter ----
g_simp_PM_OI_winter <- igraph::simplify(g_PM_OI_a_winter, remove.multiple = TRUE, remove.loops = TRUE)
el_PM_OI_winter <- igraph::as_edgelist(g_simp_PM_OI_winter, names = TRUE)
net_PM_OI_winter <- network::as.network(el_PM_OI_winter, directed = TRUE)
network::is.network(net_PM_OI_winter)

matrix_PM_OI_winter <- as.matrix(net_PM_OI_winter)
colnames(matrix_PM_OI_winter) <- NULL
rownames(matrix_PM_OI_winter) <- NULL

Mnet_PM_OI_winter <- as.network(matrix_PM_OI_winter)

res_PM_OI_winter <- SimulateExtinctions(Network = Mnet_PM_OI_winter, Method = "Mostconnected")

# 4 criterios de extincion
#in-degree
in_deg_PM_OI_winter <- degree(g_simp_PM_OI_winter, mode = "in")
#Ordenar de mayor a menor el in_degree:
orden_in_deg_PM_OI_winter <- order(in_deg_PM_OI_winter, decreasing = TRUE)

res_in_PM_OI_winter <- ExtinctionOrder(
  Network = Mnet_PM_OI_winter,
  Order = orden_in_deg_PM_OI_winter,
  NetworkType = "Trophic",
  IS = 0.5)

#Out-degree
out_deg_PM_OI_winter <- degree(g_simp_PM_OI_winter, mode = "out")
#Ordenar de mayor a menor el out_degree:
orden_out_deg_PM_OI_winter <- order(out_deg_PM_OI_winter, decreasing = TRUE)

res_out_PM_OI_winter <- ExtinctionOrder(
  Network = Mnet_PM_OI_winter,
  Order = orden_out_deg_PM_OI_winter,
  NetworkType = "Trophic",
  IS = 0.5)

#Calcular betweenness
betw_PM_OI_winter <- betweenness(
  g_simp_PM_OI_winter,
  directed = TRUE,
  weights = 1 / E(g_simp_PM_OI_winter)$weight  # más peso = más cercanía. peso alto = transición frecuente = más conectividad.
)
#Ordenar de mayor a menor el betweenness:
orden_betw_PM_OI_winter <- order(betw_PM_OI_winter, decreasing = TRUE)

res_betw_PM_OI_winter <- ExtinctionOrder(
  Network = Mnet_PM_OI_winter,
  Order = orden_betw_PM_OI_winter,
  NetworkType = "Trophic",
  IS = 0.5)

# Calcular indice conjunto (centralidad de flujo?)
#Funcion auxiliar para estandarizar
norm01 <- function(x) {
  (x - min(x, na.rm = TRUE)) /
    (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

in_deg_norm_PM_OI_winter <- norm01(in_deg_PM_OI_winter)
out_deg_norm_PM_OI_winter <- norm01(out_deg_PM_OI_winter)
betw_norm_PM_OI_winter <- norm01(betw_PM_OI_winter)

indice_PM_OI_winter <- (in_deg_norm_PM_OI_winter + out_deg_norm_PM_OI_winter + betw_norm_PM_OI_winter) / 3
#Ordenar de mayor a menor indice:
orden_indice_PM_OI_winter <- order(indice_PM_OI_winter, decreasing = TRUE)

res_indice_PM_OI_winter <- ExtinctionOrder(
  Network = Mnet_PM_OI_winter,
  Order = orden_indice_PM_OI_winter,
  NetworkType = "Trophic",
  IS = 0.5)

tabla_R50_PM_OI_winter <- data.frame(
  criterio = c(
    "In-degree",
    "Out-degree",
    "Betweenness",
    "Índice combinado"),
  R50 = c(
    res_in_PM_OI_winter$R50,
    res_out_PM_OI_winter$R50,
    res_betw_PM_OI_winter$R50,
    res_indice_PM_OI_winter$R50),
  R100 = c(
    res_in_PM_OI_winter$R100,
    res_out_PM_OI_winter$R100,
    res_betw_PM_OI_winter$R100,
    res_indice_PM_OI_winter$R100))

tabla_R50_PM_OI_winter

#Random extinction
set.seed(707)
random_ext_PM_OI_winter <- RandomExtinctions(Network= Mnet_PM_OI_winter, nsim= 100,SimNum = 17)

## 4. Magallánico: Franklin, winter ----
g_simp_PM_BF_winter <- igraph::simplify(g_PM_BF_a_winter, remove.multiple = TRUE, remove.loops = TRUE)
el_PM_BF_winter <- igraph::as_edgelist(g_simp_PM_BF_winter, names = TRUE)
net_PM_BF_winter <- network::as.network(el_PM_BF_winter, directed = TRUE)
network::is.network(net_PM_BF_winter)

matrix_PM_BF_winter <- as.matrix(net_PM_BF_winter)
colnames(matrix_PM_BF_winter) <- NULL
rownames(matrix_PM_BF_winter) <- NULL

Mnet_PM_BF_winter <- as.network(matrix_PM_BF_winter)

res_PM_BF_winter <- SimulateExtinctions(Network = Mnet_PM_BF_winter, Method = "Mostconnected")

# 4 criterios de extincion
#in-degree
in_deg_PM_BF_winter <- degree(g_simp_PM_BF_winter, mode = "in")
#Ordenar de mayor a menor el in_degree:
orden_in_deg_PM_BF_winter <- order(in_deg_PM_BF_winter, decreasing = TRUE)

res_in_PM_BF_winter <- ExtinctionOrder(
  Network = Mnet_PM_BF_winter,
  Order = orden_in_deg_PM_BF_winter,
  NetworkType = "Trophic",
  IS = 0.5)
#Out-degree
out_deg_PM_BF_winter <- degree(g_simp_PM_BF_winter, mode = "out")
#Ordenar de mayor a menor el out_degree:
orden_out_deg_PM_BF_winter <- order(out_deg_PM_BF_winter, decreasing = TRUE)

res_out_PM_BF_winter <- ExtinctionOrder(
  Network = Mnet_PM_BF_winter,
  Order = orden_out_deg_PM_BF_winter,
  NetworkType = "Trophic",
  IS = 0.5)
#Calcular betweenness
betw_PM_BF_winter <- betweenness(
  g_simp_PM_BF_winter,
  directed = TRUE,
  weights = 1 / E(g_simp_PM_BF_winter)$weight  # más peso = más cercanía. peso alto = transición frecuente = más conectividad.
)
#Ordenar de mayor a menor el betweenness:
orden_betw_PM_BF_winter <- order(betw_PM_BF_winter, decreasing = TRUE)

res_betw_PM_BF_winter <- ExtinctionOrder(
  Network = Mnet_PM_BF_winter,
  Order = orden_betw_PM_BF_winter,
  NetworkType = "Trophic",
  IS = 0.5
)

# Calcular indice conjunto (centralidad de flujo?)
#Funcion auxiliar para estandarizar
norm01 <- function(x) {
  (x - min(x, na.rm = TRUE)) /
    (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

in_deg_norm_PM_BF_winter <- norm01(in_deg_PM_BF_winter)
out_deg_norm_PM_BF_winter <- norm01(out_deg_PM_BF_winter)
betw_norm_PM_BF_winter <- norm01(betw_PM_BF_winter)

indice_PM_BF_winter <- (in_deg_norm_PM_BF_winter + out_deg_norm_PM_BF_winter + betw_norm_PM_BF_winter) / 3
#Ordenar de mayor a menor indice:
orden_indice_PM_BF_winter <- order(indice_PM_BF_winter, decreasing = TRUE)

res_indice_PM_BF_winter <- ExtinctionOrder(
  Network = Mnet_PM_BF_winter,
  Order = orden_indice_PM_BF_winter,
  NetworkType = "Trophic",
  IS = 0.5)

tabla_R50_PM_BF_winter <- data.frame(
  criterio = c(
    "In-degree",
    "Out-degree",
    "Betweenness",
    "Índice combinado"
  ),
  R50 = c(
    res_in_PM_BF_winter$R50,
    res_out_PM_BF_winter$R50,
    res_betw_PM_BF_winter$R50,
    res_indice_PM_BF_winter$R50
  ),
  R100 = c(
    res_in_PM_BF_winter$R100,
    res_out_PM_BF_winter$R100,
    res_betw_PM_BF_winter$R100,
    res_indice_PM_BF_winter$R100
  )
)

tabla_R50_PM_BF_winter

#Random extinction
set.seed(707)
random_ext_PM_BF_winter <- RandomExtinctions(Network= Mnet_PM_BF_winter, nsim= 100,SimNum = 14)

## 5. Rockhopper: Franklin, premolt ----

g_simp_PPA_BF_premolt <- igraph::simplify(g_PPA_BF_a, remove.multiple = TRUE, remove.loops = TRUE)
el_PPA_BF_premolt <- igraph::as_edgelist(g_simp_PPA_BF_premolt, names = TRUE)
net_PPA_BF_premolt <- network::as.network(el_PPA_BF_premolt, directed = TRUE)
network::is.network(net_PPA_BF_premolt)

matrix_PPA_BF_premolt <- as.matrix(net_PPA_BF_premolt)
colnames(matrix_PPA_BF_premolt) <- NULL
rownames(matrix_PPA_BF_premolt) <- NULL

Mnet_PPA_BF_premolt <- as.network(matrix_PPA_BF_premolt)

res_PPA_BF_premolt <- SimulateExtinctions(Network = Mnet_PPA_BF_premolt, Method = "Mostconnected")

# 4 criterios de extincion
#in-degree
in_deg_PPA_BF_premolt <- degree(g_simp_PPA_BF_premolt, mode = "in")
#Ordenar de mayor a menor el in_degree:
orden_in_deg_PPA_BF_premolt <- order(in_deg_PPA_BF_premolt, decreasing = TRUE)

res_in_PPA_BF_premolt <- ExtinctionOrder(
  Network = Mnet_PPA_BF_premolt,
  Order = orden_in_deg_PPA_BF_premolt,
  NetworkType = "Trophic",
  IS = 0.5
)
#Out-degree
out_deg_PPA_BF_premolt <- degree(g_simp_PPA_BF_premolt, mode = "out")
#Ordenar de mayor a menor el out_degree:
orden_out_deg_PPA_BF_premolt <- order(out_deg_PPA_BF_premolt, decreasing = TRUE)

res_out_PPA_BF_premolt <- ExtinctionOrder(
  Network = Mnet_PPA_BF_premolt,
  Order = orden_out_deg_PPA_BF_premolt,
  NetworkType = "Trophic",
  IS = 0.5
)
#Calcular betweenness
betw_PPA_BF_premolt <- betweenness(
  g_simp_PPA_BF_premolt,
  directed = TRUE,
  weights = 1 / E(g_simp_PPA_BF_premolt)$weight  # más peso = más cercanía. peso alto = transición frecuente = más conectividad.
)
#Ordenar de mayor a menor el betweenness:
orden_betw_PPA_BF_premolt <- order(betw_PPA_BF_premolt, decreasing = TRUE)

res_betw_PPA_BF_premolt <- ExtinctionOrder(
  Network = Mnet_PPA_BF_premolt,
  Order = orden_betw_PPA_BF_premolt,
  NetworkType = "Trophic",
  IS = 0.5
)

# Calcular indice conjunto (centralidad de flujo?)
#Funcion auxiliar para estandarizar
norm01 <- function(x) {
  (x - min(x, na.rm = TRUE)) /
    (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

in_deg_norm_PPA_BF_premolt <- norm01(in_deg_PPA_BF_premolt)
out_deg_norm_PPA_BF_premolt <- norm01(out_deg_PPA_BF_premolt)
betw_norm_PPA_BF_premolt <- norm01(betw_PPA_BF_premolt)

indice_PPA_BF_premolt <- (in_deg_norm_PPA_BF_premolt + out_deg_norm_PPA_BF_premolt + betw_norm_PPA_BF_premolt) / 3
#Ordenar de mayor a menor indice:
orden_indice_PPA_BF_premolt <- order(indice_PPA_BF_premolt, decreasing = TRUE)

res_indice_PPA_BF_premolt <- ExtinctionOrder(
  Network = Mnet_PPA_BF_premolt,
  Order = orden_indice_PPA_BF_premolt,
  NetworkType = "Trophic",
  IS = 0.5
)

tabla_R50_PPA_BF_premolt <- data.frame(
  criterio = c(
    "In-degree",
    "Out-degree",
    "Betweenness",
    "Índice combinado"
  ),
  R50 = c(
    res_in_PPA_BF_premolt$R50,
    res_out_PPA_BF_premolt$R50,
    res_betw_PPA_BF_premolt$R50,
    res_indice_PPA_BF_premolt$R50
  ),
  R100 = c(
    res_in_PPA_BF_premolt$R100,
    res_out_PPA_BF_premolt$R100,
    res_betw_PPA_BF_premolt$R100,
    res_indice_PPA_BF_premolt$R100
  )
)
tabla_R50_PPA_BF_premolt

#Random extinction
set.seed(707)
random_ext_PPA_BF_premolt <- RandomExtinctions(Network= Mnet_PPA_BF_premolt, nsim= 100,SimNum = 13)

## 6. Rockhopper: San Juan, premolt ----
g_simp_PPA_SJ_premolt <- igraph::simplify(g_PPA_SJ_a, remove.multiple = TRUE, remove.loops = TRUE)
el_PPA_SJ_premolt <- igraph::as_edgelist(g_simp_PPA_SJ_premolt, names = TRUE)
net_PPA_SJ_premolt <- network::as.network(el_PPA_SJ_premolt, directed = TRUE)
network::is.network(net_PPA_SJ_premolt)

matrix_PPA_SJ_premolt <- as.matrix(net_PPA_SJ_premolt)
colnames(matrix_PPA_SJ_premolt) <- NULL
rownames(matrix_PPA_SJ_premolt) <- NULL

Mnet_PPA_SJ_premolt <- as.network(matrix_PPA_SJ_premolt)
res_PPA_SJ_premolt <- SimulateExtinctions(Network = Mnet_PPA_SJ_premolt, Method = "Mostconnected")

# 4 criterios de extincion
#in-degree
in_deg_PPA_SJ_premolt <- degree(g_simp_PPA_SJ_premolt, mode = "in")
#Ordenar de mayor a menor el in_degree:
orden_in_deg_PPA_SJ_premolt <- order(in_deg_PPA_SJ_premolt, decreasing = TRUE)

res_in_PPA_SJ_premolt <- ExtinctionOrder(
  Network = Mnet_PPA_SJ_premolt,
  Order = orden_in_deg_PPA_SJ_premolt,
  NetworkType = "Trophic",
  IS = 0.5
)
#Out-degree
out_deg_PPA_SJ_premolt <- degree(g_simp_PPA_SJ_premolt, mode = "out")
#Ordenar de mayor a menor el out_degree:
orden_out_deg_PPA_SJ_premolt <- order(out_deg_PPA_SJ_premolt, decreasing = TRUE)

res_out_PPA_SJ_premolt <- ExtinctionOrder(
  Network = Mnet_PPA_SJ_premolt,
  Order = orden_out_deg_PPA_SJ_premolt,
  NetworkType = "Trophic",
  IS = 0.5
)
#Calcular betweenness
betw_PPA_SJ_premolt <- betweenness(
  g_simp_PPA_SJ_premolt,
  directed = TRUE,
  weights = 1 / E(g_simp_PPA_SJ_premolt)$weight  # más peso = más cercanía. peso alto = transición frecuente = más conectividad.
)
#Ordenar de mayor a menor el betweenness:
orden_betw_PPA_SJ_premolt <- order(betw_PPA_SJ_premolt, decreasing = TRUE)

res_betw_PPA_SJ_premolt <- ExtinctionOrder(
  Network = Mnet_PPA_SJ_premolt,
  Order = orden_betw_PPA_SJ_premolt,
  NetworkType = "Trophic",
  IS = 0.5
)
# Calcular indice conjunto (centralidad de flujo?)
#Funcion auxiliar para estandarizar
norm01 <- function(x) {
  (x - min(x, na.rm = TRUE)) /
    (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

in_deg_norm_PPA_SJ_premolt <- norm01(in_deg_PPA_SJ_premolt)
out_deg_norm_PPA_SJ_premolt <- norm01(out_deg_PPA_SJ_premolt)
betw_norm_PPA_SJ_premolt <- norm01(betw_PPA_SJ_premolt)

indice_PPA_SJ_premolt <- (in_deg_norm_PPA_SJ_premolt + out_deg_norm_PPA_SJ_premolt + betw_norm_PPA_SJ_premolt) / 3
#Ordenar de mayor a menor indice:
orden_indice_PPA_SJ_premolt <- order(indice_PPA_SJ_premolt, decreasing = TRUE)

res_indice_PPA_SJ_premolt <- ExtinctionOrder(
  Network = Mnet_PPA_SJ_premolt,
  Order = orden_indice_PPA_SJ_premolt,
  NetworkType = "Trophic",
  IS = 0.5
)

tabla_R50_PPA_SJ_premolt <- data.frame(
  criterio = c(
    "In-degree",
    "Out-degree",
    "Betweenness",
    "Índice combinado"
  ),
  R50 = c(
    res_in_PPA_SJ_premolt$R50,
    res_out_PPA_SJ_premolt$R50,
    res_betw_PPA_SJ_premolt$R50,
    res_indice_PPA_SJ_premolt$R50
  ),
  R100 = c(
    res_in_PPA_SJ_premolt$R100,
    res_out_PPA_SJ_premolt$R100,
    res_betw_PPA_SJ_premolt$R100,
    res_indice_PPA_SJ_premolt$R100
  )
)
tabla_R50_PPA_SJ_premolt

#Random extinction
set.seed(707)
random_ext_PPA_SJ_premolt <- RandomExtinctions(Network= Mnet_PPA_SJ_premolt, nsim= 100,SimNum = 7)

## 7. Rockhopper: Franklin, winter ----

g_simp_PPA_BF_winter <- igraph::simplify(g_PPA_BF_a_winter, remove.multiple = TRUE, remove.loops = TRUE)
el_PPA_BF_winter <- igraph::as_edgelist(g_simp_PPA_BF_winter, names = TRUE)
net_PPA_BF_winter <- network::as.network(el_PPA_BF_winter, directed = TRUE)
network::is.network(net_PPA_BF_winter)

matrix_PPA_BF_winter <- as.matrix(net_PPA_BF_winter)
colnames(matrix_PPA_BF_winter) <- NULL
rownames(matrix_PPA_BF_winter) <- NULL

Mnet_PPA_BF_winter <- as.network(matrix_PPA_BF_winter)

res_PPA_BF_winter <- SimulateExtinctions(Network = Mnet_PPA_BF_winter, Method = "Mostconnected")

# 4 criterios de extincion
#in-degree
in_deg_PPA_BF_winter <- degree(g_simp_PPA_BF_winter, mode = "in")
#Ordenar de mayor a menor el in_degree:
orden_in_deg_PPA_BF_winter <- order(in_deg_PPA_BF_winter, decreasing = TRUE)

res_in_PPA_BF_winter <- ExtinctionOrder(
  Network = Mnet_PPA_BF_winter,
  Order = orden_in_deg_PPA_BF_winter,
  NetworkType = "Trophic",
  IS = 0.5
)
#Out-degree
out_deg_PPA_BF_winter <- degree(g_simp_PPA_BF_winter, mode = "out")
#Ordenar de mayor a menor el out_degree:
orden_out_deg_PPA_BF_winter <- order(out_deg_PPA_BF_winter, decreasing = TRUE)

res_out_PPA_BF_winter <- ExtinctionOrder(
  Network = Mnet_PPA_BF_winter,
  Order = orden_out_deg_PPA_BF_winter,
  NetworkType = "Trophic",
  IS = 0.5
)
#Calcular betweenness
betw_PPA_BF_winter <- betweenness(
  g_simp_PPA_BF_winter,
  directed = TRUE,
  weights = 1 / E(g_simp_PPA_BF_winter)$weight  # más peso = más cercanía. peso alto = transición frecuente = más conectividad.
)
#Ordenar de mayor a menor el betweenness:
orden_betw_PPA_BF_winter <- order(betw_PPA_BF_winter, decreasing = TRUE)

res_betw_PPA_BF_winter <- ExtinctionOrder(
  Network = Mnet_PPA_BF_winter,
  Order = orden_betw_PPA_BF_winter,
  NetworkType = "Trophic",
  IS = 0.5
)

# Calcular indice conjunto (centralidad de flujo?)
#Funcion auxiliar para estandarizar
norm01 <- function(x) {
  (x - min(x, na.rm = TRUE)) /
    (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

in_deg_norm_PPA_BF_winter <- norm01(in_deg_PPA_BF_winter)
out_deg_norm_PPA_BF_winter <- norm01(out_deg_PPA_BF_winter)
betw_norm_PPA_BF_winter <- norm01(betw_PPA_BF_winter)

indice_PPA_BF_winter <- (in_deg_norm_PPA_BF_winter + out_deg_norm_PPA_BF_winter + betw_norm_PPA_BF_winter) / 3
#Ordenar de mayor a menor indice:
orden_indice_PPA_BF_winter <- order(indice_PPA_BF_winter, decreasing = TRUE)

res_indice_PPA_BF_winter <- ExtinctionOrder(
  Network = Mnet_PPA_BF_winter,
  Order = orden_indice_PPA_BF_winter,
  NetworkType = "Trophic",
  IS = 0.5
)

tabla_R50_PPA_BF_winter <- data.frame(
  criterio = c(
    "In-degree",
    "Out-degree",
    "Betweenness",
    "Índice combinado"
  ),
  R50 = c(
    res_in_PPA_BF_winter$R50,
    res_out_PPA_BF_winter$R50,
    res_betw_PPA_BF_winter$R50,
    res_indice_PPA_BF_winter$R50
  ),
  R100 = c(
    res_in_PPA_BF_winter$R100,
    res_out_PPA_BF_winter$R100,
    res_betw_PPA_BF_winter$R100,
    res_indice_PPA_BF_winter$R100
  )
)
tabla_R50_PPA_BF_winter

#Random extinction
set.seed(707)
random_ext_PPA_BF_winter <- RandomExtinctions(Network= Mnet_PPA_BF_winter, nsim= 100,SimNum = 117)

## 8. Rockhopper: San Juan, winter ----

g_simp_PPA_SJ_winter <- igraph::simplify(g_PPA_SJ_a_winter, remove.multiple = TRUE, remove.loops = TRUE)
el_PPA_SJ_winter <- igraph::as_edgelist(g_simp_PPA_SJ_winter, names = TRUE)
net_PPA_SJ_winter <- network::as.network(el_PPA_SJ_winter, directed = TRUE)
network::is.network(net_PPA_SJ_winter)

matrix_PPA_SJ_winter <- as.matrix(net_PPA_SJ_winter)
colnames(matrix_PPA_SJ_winter) <- NULL
rownames(matrix_PPA_SJ_winter) <- NULL

Mnet_PPA_SJ_winter <- as.network(matrix_PPA_SJ_winter)

res_PPA_SJ_winter <- SimulateExtinctions(Network = Mnet_PPA_SJ_winter, Method = "Mostconnected")

# 4 criterios de extincion
#in-degree
in_deg_PPA_SJ_winter <- degree(g_simp_PPA_SJ_winter, mode = "in")
#Ordenar de mayor a menor el in_degree:
orden_in_deg_PPA_SJ_winter <- order(in_deg_PPA_SJ_winter, decreasing = TRUE)

res_in_PPA_SJ_winter <- ExtinctionOrder(
  Network = Mnet_PPA_SJ_winter,
  Order = orden_in_deg_PPA_SJ_winter,
  NetworkType = "Trophic",
  IS = 0.5
)
#Out-degree
out_deg_PPA_SJ_winter <- degree(g_simp_PPA_SJ_winter, mode = "out")
#Ordenar de mayor a menor el out_degree:
orden_out_deg_PPA_SJ_winter <- order(out_deg_PPA_SJ_winter, decreasing = TRUE)

res_out_PPA_SJ_winter <- ExtinctionOrder(
  Network = Mnet_PPA_SJ_winter,
  Order = orden_out_deg_PPA_SJ_winter,
  NetworkType = "Trophic",
  IS = 0.5
)
#Calcular betweenness
betw_PPA_SJ_winter <- betweenness(
  g_simp_PPA_SJ_winter,
  directed = TRUE,
  weights = 1 / E(g_simp_PPA_SJ_winter)$weight  # más peso = más cercanía. peso alto = transición frecuente = más conectividad.
)
#Ordenar de mayor a menor el betweenness:
orden_betw_PPA_SJ_winter <- order(betw_PPA_SJ_winter, decreasing = TRUE)

res_betw_PPA_SJ_winter <- ExtinctionOrder(
  Network = Mnet_PPA_SJ_winter,
  Order = orden_betw_PPA_SJ_winter,
  NetworkType = "Trophic",
  IS = 0.5
)

# Calcular indice conjunto (centralidad de flujo?)
#Funcion auxiliar para estandarizar
norm01 <- function(x) {
  (x - min(x, na.rm = TRUE)) /
    (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

in_deg_norm_PPA_SJ_winter <- norm01(in_deg_PPA_SJ_winter)
out_deg_norm_PPA_SJ_winter <- norm01(out_deg_PPA_SJ_winter)
betw_norm_PPA_SJ_winter <- norm01(betw_PPA_SJ_winter)

indice_PPA_SJ_winter <- (in_deg_norm_PPA_SJ_winter + out_deg_norm_PPA_SJ_winter + betw_norm_PPA_SJ_winter) / 3
#Ordenar de mayor a menor indice:
orden_indice_PPA_SJ_winter <- order(indice_PPA_SJ_winter, decreasing = TRUE)

res_indice_PPA_SJ_winter <- ExtinctionOrder(
  Network = Mnet_PPA_SJ_winter,
  Order = orden_indice_PPA_SJ_winter,
  NetworkType = "Trophic",
  IS = 0.5
)

tabla_R50_PPA_SJ_winter <- data.frame(
  criterio = c(
    "In-degree",
    "Out-degree",
    "Betweenness",
    "Índice combinado"
  ),
  R50 = c(
    res_in_PPA_SJ_winter$R50,
    res_out_PPA_SJ_winter$R50,
    res_betw_PPA_SJ_winter$R50,
    res_indice_PPA_SJ_winter$R50
  ),
  R100 = c(
    res_in_PPA_SJ_winter$R100,
    res_out_PPA_SJ_winter$R100,
    res_betw_PPA_SJ_winter$R100,
    res_indice_PPA_SJ_winter$R100
  )
)
tabla_R50_PPA_SJ_winter

#Random extinction
set.seed(707)
random_ext_PPA_SJ_winter <- RandomExtinctions(Network= Mnet_PPA_SJ_winter, nsim= 100,SimNum = 49)

## 9. Integración y comparación de resultados R50 ----
armar_tabla_R50 <- function(tabla_dirigida, random_obj,
                            especie, colonia, etapa) {

  # Extinciones dirigidas
  df_dir <- tabla_dirigida %>%
    mutate(
      especie = especie,
      colonia = colonia,
      etapa = etapa,
      tipo = "Dirigida"
    ) %>%
    dplyr::select(especie, colonia, etapa, criterio, R50, tipo)

  # Extinción aleatoria
  df_rand <- data.frame(
    especie = especie,
    colonia = colonia,
    etapa = etapa,
    criterio = "Random",
    R50 = random_obj$R50result["Mean"],
    tipo = "Aleatoria",
    CI_random = random_obj$R50result["CI"]
  )

  bind_rows(df_dir, df_rand)
}

tabla_R50_total_1 <- bind_rows(

  # Magallánicos – Observatorio
  armar_tabla_R50(tabla_R50_PM_OI_premolt, random_ext_PM_OI_premolt,
                  "Magellanic", "Observatorio", "premolt"),

  armar_tabla_R50(tabla_R50_PM_OI_winter, random_ext_PM_OI_winter,
                  "Magellanic", "Observatorio", "winter"),

  # Magallánicos – Franklin
  armar_tabla_R50(tabla_R50_PM_BF_premolt, random_ext_PM_BF_premolt,
                  "Magellanic", "Franklin", "premolt"),

  armar_tabla_R50(tabla_R50_PM_BF_winter, random_ext_PM_BF_winter,
                  "Magellanic", "Franklin", "winter"),

  # Rockhopper – Franklin
  armar_tabla_R50(tabla_R50_PPA_BF_premolt, random_ext_PPA_BF_premolt,
                  "Rockhoppers", "Franklin", "premolt"),

  armar_tabla_R50(tabla_R50_PPA_BF_winter, random_ext_PPA_BF_winter,
                  "Rockhoppers", "Franklin", "winter"),

  # Rockhopper – San Juan
  armar_tabla_R50(tabla_R50_PPA_SJ_premolt, random_ext_PPA_SJ_premolt,
                  "Rockhoppers", "San Juan", "premolt"),

  armar_tabla_R50(tabla_R50_PPA_SJ_winter, random_ext_PPA_SJ_winter,
                  "Rockhoppers", "San Juan", "winter"))

tabla_R50_total_1

tabla_resumen <- tabla_R50_total_1 %>%
  group_by(especie, colonia, etapa) %>%
  summarise(
    R50_min_dirigido = min(R50[tipo == "Dirigida"]),
    R50_random = R50[tipo == "Aleatoria"][1])

# Identificar si R50= Random R50 o si R50<Random R50
# Criterio: “= random” si R50 cae dentro del intervalo de confianza del R50_random.

tabla_R50_total_2 <- tabla_R50_total_1 %>%
  group_by(especie, colonia, etapa) %>%
  mutate(
    R50_random = R50[criterio == "Random"],
    CI_random  = CI_random[criterio == "Random"],
    similar_random = ifelse(
      tipo == "Dirigida" &
        R50 >= (R50_random - CI_random) &
        R50 <= (R50_random + CI_random),
      TRUE, FALSE
    )

  ) %>%
  ungroup()

tabla_R50_total_2 <- tabla_R50_total_1 %>%
  group_by(especie, colonia, etapa) %>%
  mutate(
    R50_random = R50[criterio == "Random"],
    CI_random  = CI_random[criterio == "Random"],
    # Limites
    CI_inf_random = R50_random - CI_random,
    CI_sup_random = R50_random + CI_random,
    # Comparar dirigida vs random
    similar_random = ifelse(
      tipo == "Dirigida" &
        R50 >= CI_inf_random &
        R50 <= CI_sup_random,
      TRUE, FALSE)) %>%
  ungroup()

tabla_R50_IS050 <- tabla_R50_total_2 %>% filter(tipo == "Dirigida") %>%
  dplyr::select(especie, colonia, etapa, criterio, R50_050 = R50)

View(tabla_R50_IS050)
df_dir <- tabla_R50_total_2 %>%
  filter(tipo == "Dirigida")

df_rand <- tabla_R50_total_2 %>%
  filter(criterio == "Random")

# quiero fijar el orden del ploteo de los criterios:

df_dir$criterio <- factor(
  df_dir$criterio,
  levels = c(
    "Random",
    "Índice combinado",
    "Betweenness",
    "Out-degree",
    "In-degree"))

ggplot(df_dir,aes(x = R50, y = criterio, fill = similar_random)) +
  geom_col(width = 0.6) +
  # Línea vertical del R50 aleatorio
  geom_vline(data = df_rand,aes(xintercept = R50, linetype = "Random extinction"),
             linetype = "dashed",linewidth = 0.8,color = "black") +
  facet_grid(especie ~ colonia + etapa) +
  scale_fill_manual(values = c("FALSE" = "gray60", "TRUE" = "firebrick"),name = "≈ Random") +
    labs(x = "R50 ",y = "Criterio de extinción",title = "Robustez (R50)",
    subtitle = "Línea punteada: R50 esperado bajo extinción aleatoria") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "gray90"),strip.text = element_text(face = "bold"),
    legend.position = "top")
# renombrar el criterio al inglés antes de plotear
df_dir$criterio <- recode(df_dir$criterio, "Índice combinado" = "Composite index")
df_rand$criterio <- recode(df_rand$criterio, "Índice combinado" = "Composite index")

ggplot() +
  # Intervalo de confianza del random
  geom_errorbarh(data = df_rand,aes(y = criterio,xmin = R50 - CI_random,xmax = R50 + CI_random),
    height = 0.25,linewidth = 0.8,color = "gray40") +
  # Línea vertical del R50 random medio
  geom_vline(data = df_rand,aes(xintercept = R50, linetype = "Random extinction"),linewidth = 0.8,
    color = "black") +
  # Puntos: R50 dirigido
  geom_point(data = df_dir,aes(x = R50,y = criterio,color = similar_random),size = 3) +
  facet_grid(especie ~ colonia + etapa) +
  scale_color_manual(values = c("FALSE" = "gray60", "TRUE" = "firebrick"),name = "≈ Random") +
  # scale_color_manual(
  #   values = c("FALSE" = "firebrick", "TRUE" = "gray70"),   # invertido: destaca el que supera random
  #   name = "More efficient\nthan random",
  #   labels = c("FALSE" = "Yes", "TRUE" = "No")) +
  #
  scale_linetype_manual(values = c("Random extinction" = "dashed"),name = "") +
  labs(x = "R50 ",y = "Extinction criterion", title = "Network robustness (R50)") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "gray90"),strip.text = element_text(face = "bold"),
    legend.position = "top")

# sin los paneles vacios
df_dir <- df_dir %>%
  mutate(panel = paste(especie, colonia, etapa, sep = " – "))
df_rand <- df_rand %>%
  mutate(panel = paste(especie, colonia, etapa, sep = " – "))

ggplot() +
  # Intervalo de confianza del random
  geom_errorbarh(data = df_rand,aes(y = criterio,xmin = R50 - CI_random,xmax = R50 + CI_random),
                 height = 0.25,linewidth = 0.8,color = "gray40") +
  # Línea vertical del R50 random medio
  geom_vline(data = df_rand,aes(xintercept = R50, linetype = "Random extinction"),linewidth = 0.8,
             color = "black") +
  # Puntos: R50 dirigido
  geom_point(data = df_dir,aes(x = R50,y = criterio,color = similar_random),size = 3) +
  #facet_grid(especie ~ colonia + etapa) +
  facet_wrap(~ panel, ncol = 3, scales = "free_y")+
  scale_color_manual(values = c("FALSE" = "gray60", "TRUE" = "firebrick"),name = "≈ Random") +
  # scale_color_manual(
  #   values = c("FALSE" = "firebrick", "TRUE" = "gray70"),   # invertido: destaca el que supera random
  #   name = "More efficient\nthan random",
  #   labels = c("FALSE" = "Yes", "TRUE" = "No")) +
  #
  scale_linetype_manual(values = c("Random extinction" = "dashed"),name = "") +
  labs(x = "R50 ",y = "Extinction criterion", title = "Network robustness (R50)") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "gray90"),strip.text = element_text(face = "bold"),
        legend.position = "top")

# etiqueta de panel y orden
df_dir <- df_dir %>%
  mutate(panel = paste(especie, colonia, etapa, sep = " – "),
         criterio = recode(criterio, "Índice combinado" = "Composite index"))
df_rand <- df_rand %>%
  mutate(panel = paste(especie, colonia, etapa, sep = " – "),
         criterio = recode(criterio, "Índice combinado" = "Composite index"))

orden_paneles <- c(
  "Magellanic – Observatorio – premolt", "Magellanic – Franklin – premolt",
  "Magellanic – Observatorio – winter",  "Magellanic – Franklin – winter",
  "Rockhoppers – San Juan – premolt",    "Rockhoppers – Franklin – premolt",
  "Rockhoppers – San Juan – winter",     "Rockhoppers – Franklin – winter")
df_dir$panel <- factor(df_dir$panel,  levels = orden_paneles)
df_rand$panel <- factor(df_rand$panel, levels = orden_paneles)

# orden de criterios en el eje y
orden_crit <- c("Composite index","Betweenness","Out-degree","In-degree")
df_dir$criterio <- factor(df_dir$criterio,  levels = orden_crit)

ggplot() +
  # geom_errorbarh(data = df_rand,
  #                aes(y = criterio, xmin = R50 - CI_random, xmax = R50 + CI_random),
  #                height = 0.25, linewidth = 0.7, color = "gray40") +
   geom_rect(data = df_rand,
            aes(xmin = R50 - CI_random, xmax = R50 + CI_random,
                ymin = -Inf, ymax = Inf),
            fill = "gray70", alpha = 0.3, inherit.aes = FALSE) +
  geom_vline(data = df_rand,
             aes(xintercept = R50), linetype = "dashed", linewidth = 0.7, color = "black") +
  geom_point(data = df_dir,
             aes(x = R50, y = criterio, color = similar_random), size = 2.8) +
  facet_wrap(~ panel, ncol = 4) +
  scale_color_manual(
    values = c("FALSE" = "firebrick", "TRUE" = "gray70"),
    name = "More efficient than random",
    labels = c("FALSE" = "Yes", "TRUE" = "No")) +
  labs(x = "R50", y = "Extinction criterion") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "gray90"),
        strip.text = element_text(face = "bold", size = 7),
        legend.position = "top",
        panel.spacing = unit(0.6, "lines"))

ggplot() +
  geom_rect(data = df_rand,
            aes(xmin = R50 - CI_random, xmax = R50 + CI_random,
                ymin = -Inf, ymax = Inf),
            fill = "gray70", alpha = 0.3, inherit.aes = FALSE) +
  geom_vline(data = df_rand,
             aes(xintercept = R50), linetype = "dashed", linewidth = 0.7, color = "black") +
  geom_point(data = df_dir,
             aes(x = R50, y = criterio, color = similar_random), size = 2.8) +
  facet_wrap(~ panel, ncol = 4) +
  scale_color_manual(
    values = c("FALSE" = "firebrick", "TRUE" = "gray70"),
    name = "More efficient than random",
    labels = c("FALSE" = "Yes", "TRUE" = "No")) +
  scale_x_continuous(
    breaks = seq(0, 1, by = 0.1),
    minor_breaks = seq(0, 1, by = 0.05),
    labels = function(x) sprintf("%.2f", x)) +
  labs(x = "R50", y = "Extinction criterion") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "gray90"),
        strip.text = element_text(face = "bold", size = 7),
        legend.position = "top",
        panel.spacing = unit(0.6, "lines"))

# recodificar criterio y fijar orden
df_dir <- df_dir %>%
  mutate(criterio = recode(criterio, "Índice combinado" = "Composite index"))
df_rand <- df_rand %>%
  mutate(criterio = recode(criterio, "Índice combinado" = "Composite index"))

orden_crit <- c("Composite index","Betweenness","Out-degree","In-degree")
df_dir$criterio <- factor(df_dir$criterio, levels = orden_crit)

# orden de etapa (premolt arriba, winter abajo) y especie
df_dir$etapa <- factor(df_dir$etapa,   levels = c("premolt","winter"))
df_rand$etapa <- factor(df_rand$etapa,  levels = c("premolt","winter"))
df_dir$especie <- factor(df_dir$especie,  levels = c("Magellanic","Rockhoppers"))
df_rand$especie <- factor(df_rand$especie, levels = c("Magellanic","Rockhoppers"))

ggplot() +
  geom_rect(data = df_rand,
            aes(xmin = R50 - CI_random, xmax = R50 + CI_random,
                ymin = -Inf, ymax = Inf),
            fill = "gray70", alpha = 0.3, inherit.aes = FALSE) +
  geom_vline(data = df_rand,
             aes(xintercept = R50), linetype = "dashed", linewidth = 0.7, color = "black") +
  geom_point(data = df_dir,
             aes(x = R50, y = criterio, color = similar_random), size = 2.8) +
  facet_nested(etapa ~ especie + colonia, scales = "free_x") +
  scale_color_manual(
    values = c("FALSE" = "firebrick", "TRUE" = "gray70"),
    name = "More efficient than random",
    labels = c("FALSE" = "Yes", "TRUE" = "No")) +
  scale_x_continuous(
    breaks = seq(0, 1, by = 0.1),
    minor_breaks = NULL,
    labels = function(x) sprintf("%.2f", x)) +
  coord_cartesian(xlim = c(0.2, 0.65))  +
  labs(x = "R50", y = "Extinction criterion") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "gray90"),
        strip.text = element_text(face = "bold", size = 10),
        legend.position = "top",
        panel.spacing = unit(0.6, "lines"))

random_ext_PPA_SJ_winter$R50result   # ¿Mean y CI?
# y mirá la dispersión real de las réplicas:
str(random_ext_PPA_SJ_winter, max.level = 1)

random_ext_PPA_SJ_premolt$R50result   # ¿CI es 0 o >0 ahora?
random_ext_PPA_SJ_winter$R50result

vcount(g_simp_PPA_SJ_winter)   # ¿sigue teniendo los mismos nodos que antes?
ecount(g_simp_PPA_SJ_winter)   # ¿las mismas aristas?

# en cuantas redes cada criterio es mejor que el random?::
df_dir %>%
  group_by(criterio) %>%
  summarise(n_redes        = n(),
    n_mejor_random = sum(similar_random == FALSE),  # FALSE = mejor que random
    prop_mejor     = mean(similar_random == FALSE)) %>%
  arrange(desc(prop_mejor))

#R50 promedio de cada criterio en las redes donde es informativo
df_dir %>%
  filter(similar_random == FALSE) %>%   # solo donde supera al azar
  group_by(criterio) %>%
  summarise(R50_medio = mean(R50), R50_min = min(R50))
