# ============================================================
#  05 - ROBUSTEZ DE LAS REDES (R50)
#  Para cada red se remueven nodos en cuatro órdenes dirigidos
#  (in-degree, out-degree, betweenness, índice combinado) y se
#  compara el R50 con el de la extinción aleatoria (100 réplicas).
#
#  Construcción de la red para NetworkExtinction: lista de aristas
#  del grafo igraph -> objeto network -> matriz de adyacencia sin
#  nombres -> objeto network.
#
#  Salida: tabla_R50_total_2 (R50 por red y criterio), df_dir,
#          df_rand y la figura de robustez.
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
#  Función: robustez de una red (mismos pasos, en el mismo orden,
#  que el análisis original)
#   - simnum: nodos a remover en cada réplica aleatoria
# ------------------------------------------------------------
robustez_red <- function(g, simnum) {
  
  # grafo simple -> objeto network
  g_simp <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
  el     <- igraph::as_edgelist(g_simp, names = TRUE)
  net    <- network::as.network(el, directed = TRUE)
  m      <- as.matrix(net)
  colnames(m) <- NULL
  rownames(m) <- NULL
  Mnet <- as.network(m)
  
  res_mc <- SimulateExtinctions(Network = Mnet, Method = "Mostconnected")
  
  # métricas que definen el orden de extinción
  in_deg  <- igraph::degree(g_simp, mode = "in")
  out_deg <- igraph::degree(g_simp, mode = "out")
  betw    <- igraph::betweenness(g_simp, directed = TRUE,
                                 weights = 1 / E(g_simp)$weight)
  indice  <- (norm01(in_deg) + norm01(out_deg) + norm01(betw)) / 3
  
  # extinción dirigida: de mayor a menor valor de cada métrica
  ext <- function(metrica) {
    ExtinctionOrder(Network = Mnet, Order = order(metrica, decreasing = TRUE),
                    NetworkType = "Trophic", IS = 0.5)
  }
  res_in     <- ext(in_deg)
  res_out    <- ext(out_deg)
  res_betw   <- ext(betw)
  res_indice <- ext(indice)
  
  tabla <- data.frame(
    criterio = c("In-degree", "Out-degree", "Betweenness", "Índice combinado"),
    R50  = c(res_in$R50,  res_out$R50,  res_betw$R50,  res_indice$R50),
    R100 = c(res_in$R100, res_out$R100, res_betw$R100, res_indice$R100))
  
  # extinción aleatoria
  set.seed(707)
  aleatoria <- RandomExtinctions(Network = Mnet, nsim = 100, SimNum = simnum)
  
  list(tabla = tabla, aleatoria = aleatoria)
}

# ------------------------------------------------------------
#  Aplicar a las 8 redes
# ------------------------------------------------------------
redes <- list(
  list(g = g_PM_OI_a,         simnum = 7,   especie = "Magellanic",  colonia = "Observatorio", etapa = "premolt"),
  list(g = g_PM_BF_a,         simnum = 8,   especie = "Magellanic",  colonia = "Franklin",     etapa = "premolt"),
  list(g = g_PM_OI_a_winter,  simnum = 17,  especie = "Magellanic",  colonia = "Observatorio", etapa = "winter"),
  list(g = g_PM_BF_a_winter,  simnum = 14,  especie = "Magellanic",  colonia = "Franklin",     etapa = "winter"),
  list(g = g_PPA_BF_a,        simnum = 13,  especie = "Rockhoppers", colonia = "Franklin",     etapa = "premolt"),
  list(g = g_PPA_SJ_a,        simnum = 7,   especie = "Rockhoppers", colonia = "San Juan",     etapa = "premolt"),
  list(g = g_PPA_BF_a_winter, simnum = 117, especie = "Rockhoppers", colonia = "Franklin",     etapa = "winter"),
  list(g = g_PPA_SJ_a_winter, simnum = 49,  especie = "Rockhoppers", colonia = "San Juan",     etapa = "winter"))

resultados <- lapply(redes, function(r) robustez_red(r$g, r$simnum))

# ------------------------------------------------------------
#  Tabla de R50: dirigidas + aleatoria por red
# ------------------------------------------------------------
tabla_R50_total_1 <- bind_rows(lapply(seq_along(redes), function(i) {
  r   <- redes[[i]]
  res <- resultados[[i]]
  bind_rows(
    res$tabla %>%
      mutate(especie = r$especie, colonia = r$colonia, etapa = r$etapa,
             tipo = "Dirigida") %>%
      dplyr::select(especie, colonia, etapa, criterio, R50, tipo),
    data.frame(especie = r$especie, colonia = r$colonia, etapa = r$etapa,
               criterio = "Random",
               R50 = res$aleatoria$R50result["Mean"],
               tipo = "Aleatoria",
               CI_random = res$aleatoria$R50result["CI"]))
}))

# similar_random = TRUE si el R50 dirigido cae dentro de aleatorio +/- CI
tabla_R50_total_2 <- tabla_R50_total_1 %>%
  group_by(especie, colonia, etapa) %>%
  mutate(R50_random    = R50[criterio == "Random"],
         CI_random     = CI_random[criterio == "Random"],
         CI_inf_random = R50_random - CI_random,
         CI_sup_random = R50_random + CI_random,
         similar_random = tipo == "Dirigida" &
           R50 >= CI_inf_random & R50 <= CI_sup_random) %>%
  ungroup()

df_dir  <- tabla_R50_total_2 %>% filter(tipo == "Dirigida") %>%
  mutate(criterio = recode(criterio, "Índice combinado" = "Composite index"))
df_rand <- tabla_R50_total_2 %>% filter(criterio == "Random")

print(df_dir %>% dplyr::select(especie, colonia, etapa, criterio, R50, similar_random))

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
df_dir$etapa    <- factor(df_dir$etapa,    levels = c("premolt", "winter"))
df_rand$etapa   <- factor(df_rand$etapa,   levels = c("premolt", "winter"))
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
