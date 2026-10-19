# ============================================================
#  05 - ROBUSTEZ DE LAS REDES (R50)
#  Para cada red se remueven nodos en cuatro órdenes dirigidos
#  (in-degree, out-degree, betweenness, índice compuesto) y se
#  compara el R50 con el de la extinción aleatoria (100 réplicas).
#
#  Salida: df_dir (R50 por red y criterio), df_rand (extinción
#          aleatoria por red) y la figura de robustez.
#
#  Requiere 01 y 02 en la misma sesión.
# ============================================================

library(igraph)
library(network)
library(NetworkExtinction)
library(dplyr)
library(ggplot2)
library(ggh4x)

if (!exists("g_PM_OI_a_winter")) stop("Correr primero 01_redes_premolt.R y 02_redes_invierno.R")

norm01 <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

# ------------------------------------------------------------
#  Función: robustez de una red
#   La red 'network' se arma desde la matriz de adyacencia de
#   igraph, así los índices de order() corresponden a los mismos
#   nodos que en el grafo.
#   - incluir_aislados = FALSE: se quitan los nodos sin ninguna
#     transición (no aportan a la conectividad)
# ------------------------------------------------------------
robustez_red <- function(g, nsim = 100, semilla = 707, incluir_aislados = FALSE) {

  g <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
  aislados <- sum(igraph::degree(g) == 0)
  if (!incluir_aislados) g <- delete_vertices(g, which(igraph::degree(g) == 0))

  M    <- as_adjacency_matrix(g, sparse = FALSE)
  Mnet <- network::as.network(unname(M), matrix.type = "adjacency", directed = TRUE)

  # métricas que definen el orden de extinción
  in_deg  <- igraph::degree(g, mode = "in")
  out_deg <- igraph::degree(g, mode = "out")
  betw    <- igraph::betweenness(g, directed = TRUE, weights = 1 / E(g)$weight)
  indice  <- (norm01(in_deg) + norm01(out_deg) + norm01(betw)) / 3

  metricas <- list(`In-degree` = in_deg, `Out-degree` = out_deg,
                   Betweenness = betw, `Composite index` = indice)

  # extinción dirigida: de mayor a menor valor de cada métrica
  res <- lapply(metricas, function(m) {
    ExtinctionOrder(Network = Mnet, Order = order(m, decreasing = TRUE),
                    NetworkType = "Trophic", IS = 0.5)
  })
  dirigida <- data.frame(criterio = names(res),
                         R50  = sapply(res, function(r) r$R50),
                         R100 = sapply(res, function(r) r$R100),
                         row.names = NULL)

  # extinción aleatoria
  set.seed(semilla)
  aleatoria <- RandomExtinctions(Network = Mnet, nsim = nsim,
                                 SimNum = network::network.size(Mnet))

  list(dirigida   = dirigida,
       random_R50 = as.numeric(aleatoria$R50result[["Mean"]]),
       random_CI  = as.numeric(aleatoria$R50result[["CI"]]),
       n_nodos    = vcount(g),
       n_aislados = aislados)
}

# ------------------------------------------------------------
#  Aplicar a las 8 redes
# ------------------------------------------------------------
redes <- list(
  PM_OI_premolt  = list(g = g_PM_OI_a,        especie = "Magellanic",  colonia = "Observatorio", etapa = "premolt"),
  PM_BF_premolt  = list(g = g_PM_BF_a,        especie = "Magellanic",  colonia = "Franklin",     etapa = "premolt"),
  PM_OI_winter   = list(g = g_PM_OI_a_winter, especie = "Magellanic",  colonia = "Observatorio", etapa = "winter"),
  PM_BF_winter   = list(g = g_PM_BF_a_winter, especie = "Magellanic",  colonia = "Franklin",     etapa = "winter"),
  PPA_BF_premolt = list(g = g_PPA_BF_a,        especie = "Rockhoppers", colonia = "Franklin",     etapa = "premolt"),
  PPA_SJ_premolt = list(g = g_PPA_SJ_a,        especie = "Rockhoppers", colonia = "San Juan",     etapa = "premolt"),
  PPA_BF_winter  = list(g = g_PPA_BF_a_winter, especie = "Rockhoppers", colonia = "Franklin",     etapa = "winter"),
  PPA_SJ_winter  = list(g = g_PPA_SJ_a_winter, especie = "Rockhoppers", colonia = "San Juan",     etapa = "winter"))

resultados <- lapply(redes, function(r) robustez_red(r$g))

# control: nodos de la red analizada y nodos aislados quitados
sapply(resultados, function(x) c(nodos = x$n_nodos, aislados = x$n_aislados))

# ------------------------------------------------------------
#  Tabla de R50: cada criterio dirigido vs el intervalo aleatorio
#  (similar_random = TRUE si R50 cae dentro de aleatorio +/- CI)
# ------------------------------------------------------------
df_dir <- bind_rows(lapply(names(redes), function(nm) {
  data.frame(especie = redes[[nm]]$especie,
             colonia = redes[[nm]]$colonia,
             etapa   = redes[[nm]]$etapa,
             resultados[[nm]]$dirigida,
             R50_random = resultados[[nm]]$random_R50,
             CI_random  = resultados[[nm]]$random_CI)
})) %>%
  mutate(similar_random = R50 >= R50_random - CI_random &
                          R50 <= R50_random + CI_random)

df_rand <- df_dir %>%
  distinct(especie, colonia, etapa, R50_random, CI_random) %>%
  rename(R50 = R50_random)

print(df_dir)

# ------------------------------------------------------------
#  Resumen: en cuántas redes cada criterio supera al azar
# ------------------------------------------------------------
df_dir %>%
  group_by(criterio) %>%
  summarise(n_redes        = n(),
            n_mejor_random = sum(!similar_random),
            prop_mejor     = mean(!similar_random)) %>%
  arrange(desc(prop_mejor))

# R50 medio de cada criterio en las redes donde supera al azar
df_dir %>%
  filter(!similar_random) %>%
  group_by(criterio) %>%
  summarise(R50_medio = mean(R50), R50_min = min(R50))

# ------------------------------------------------------------
#  Figura: R50 por criterio; la banda gris es el intervalo de la
#  extinción aleatoria y la línea punteada su media
# ------------------------------------------------------------
df_dir$criterio <- factor(df_dir$criterio,
  levels = c("Composite index", "Betweenness", "Out-degree", "In-degree"))
df_dir$etapa    <- factor(df_dir$etapa,   levels = c("premolt", "winter"))
df_rand$etapa   <- factor(df_rand$etapa,  levels = c("premolt", "winter"))
df_dir$especie  <- factor(df_dir$especie,  levels = c("Magellanic", "Rockhoppers"))
df_rand$especie <- factor(df_rand$especie, levels = c("Magellanic", "Rockhoppers"))

ggplot() +
  geom_rect(data = df_rand,
            aes(xmin = R50 - CI_random, xmax = R50 + CI_random,
                ymin = -Inf, ymax = Inf),
            fill = "gray70", alpha = 0.3, inherit.aes = FALSE) +
  geom_vline(data = df_rand, aes(xintercept = R50),
             linetype = "dashed", linewidth = 0.7, color = "black") +
  geom_point(data = df_dir,
             aes(x = R50, y = criterio, color = similar_random), size = 2.8) +
  facet_nested(etapa ~ especie + colonia, scales = "free_x") +
  scale_color_manual(
    values = c("FALSE" = "firebrick", "TRUE" = "gray70"),
    name = "More efficient than random",
    labels = c("FALSE" = "Yes", "TRUE" = "No")) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1), minor_breaks = NULL,
                     labels = function(x) sprintf("%.2f", x)) +
  coord_cartesian(xlim = c(0.2, 0.65)) +
  labs(x = "R50", y = "Extinction criterion") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "gray90"),
        strip.text = element_text(face = "bold", size = 10),
        legend.position = "top",
        panel.spacing = unit(0.6, "lines"))
