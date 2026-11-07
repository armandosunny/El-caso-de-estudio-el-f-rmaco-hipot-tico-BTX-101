# ==============================================================
# CLASE MUESTRA - BIOTECNOLOGIA / ONCOLOGIA SIMULADA
# Tres pruebas de una muestra, diagnosticos y visualizaciones
# SOLO PARA DOCENCIA: NO SON DATOS CLINICOS REALES
# ==============================================================

# 0. PAQUETES
paquetes <- c("ggplot2", "ggstatsplot", "effectsize", "outliers", "BayesFactor")
faltantes <- setdiff(paquetes, rownames(installed.packages()))
if (length(faltantes) > 0) install.packages(faltantes)
library(ggplot2)
library(ggstatsplot)
library(effectsize)
library(outliers)
library(BayesFactor)

# 1. CARPETA Y DATOS SIMULADOS (no seleccionar semillas por valor p)
dir.create("resultados_oncologia", showWarnings = FALSE)
set.seed(2026)
n <- 60
mu0_t <- 25          # % de reduccion tumoral de referencia (hipotetico)
mu0_z <- 15          # ng/mL de reduccion del biomarcador de referencia (hipotetico)
sigma_conocida <- 8 # ng/mL; conocida de validacion externa HIPOTETICA
p0 <- 0.50           # proporcion de respondedores de referencia (hipotetica)

# Reduccion del tamano tumoral (%) tras 8 semanas, positiva = mejoria
reduccion <- rnorm(n, mean = 34, sd = 11)
# Reduccion del biomarcador (ng/mL), positiva = mejoria
biomarcador <- rnorm(n, mean = 19, sd = sigma_conocida)
# Respuesta tumoral (definicion didactica: reduccion >=30%)
respuesta <- as.integer(reduccion >= 30)

datos <- data.frame(
  paciente = seq_len(n),
  reduccion_tumoral_pct = reduccion,
  reduccion_biomarcador_ng_ml = biomarcador,
  respuesta = respuesta
)
write.csv(datos, "resultados_oncologia/datos_oncologia_simulados.csv", row.names = FALSE)

# 2. ESTADISTICA DESCRIPTIVA
x <- datos$reduccion_tumoral_pct
cat("\n--- DESCRIPTIVOS: REDUCCION TUMORAL ---\n")
print(summary(x))
cat("n=",length(x)," media=",mean(x)," DE=",sd(x)," mediana=",median(x),"\n")

# 3. DETECCION EXPLORATORIA DE OUTLIERS - IQR
Q1 <- unname(quantile(x, 0.25)); Q3 <- unname(quantile(x, 0.75))
RIQ <- IQR(x)
lim_inf <- Q1 - 1.5 * RIQ
lim_sup <- Q3 + 1.5 * RIQ
datos$outlier_IQR <- x < lim_inf | x > lim_sup
cat("\n--- IQR ---\n")
cat("Q1:",Q1," Q3:",Q3," IQR:",RIQ," limites:",lim_inf,lim_sup,"\n")
print(datos[datos$outlier_IQR, c("paciente", "reduccion_tumoral_pct")])
cat("Outliers mediante boxplot.stats:", paste(boxplot.stats(x)$out, collapse=", "), "\n")

g1 <- ggplot(datos, aes(x = "Pacientes", y = reduccion_tumoral_pct)) +
  geom_boxplot(width = 0.32, fill = "#7CB8C9", color = "#24566A", outlier.shape = NA, alpha = .6) +
  geom_jitter(aes(color = outlier_IQR), width = .09, height = 0, size = 2.5) +
  scale_color_manual(values = c(`FALSE`="#28678B", `TRUE`="#CB514F"), name="Posible outlier IQR") +
  geom_hline(yintercept = mu0_t, linetype = 2, color = "#B45E36") +
  labs(title="Reduccion tumoral: boxplot y datos individuales", subtitle="Criterio exploratorio IQR; no elimina observaciones", x=NULL, y="Reduccion del tumor (%)") +
  theme_minimal(base_size = 14)
print(g1)
ggsave("resultados_oncologia/01_boxplot_IQR.png", g1, width=9, height=6, dpi=300)


# ==========================================================
# GRAFICA DEL IQR Y LIMITES DE DETECCION DE OUTLIERS
# ==========================================================

library(ggplot2)

g_IQR <- ggplot(
  datos,
  aes(x = paciente, y = reduccion_tumoral_pct)
) +

  # Región intercuartílica: Q1 a Q3
  annotate(
    "rect",
    xmin = -Inf, xmax = Inf,
    ymin = Q1, ymax = Q3,
    fill = "#7BC8B6",
    alpha = 0.15
  ) +

  # Observaciones individuales
  geom_point(
    aes(color = outlier_IQR),
    size = 2.8,
    alpha = 0.85
  ) +

  # Límite inferior IQR
  geom_hline(
    yintercept = lim_inf,
    color = "#D64F4F",
    linetype = "dashed",
    linewidth = 0.9
  ) +

  # Límite superior IQR
  geom_hline(
    yintercept = lim_sup,
    color = "#D64F4F",
    linetype = "dashed",
    linewidth = 0.9
  ) +

  # Primer cuartil
  geom_hline(
    yintercept = Q1,
    color = "#2A9D8F",
    linetype = "dotted",
    linewidth = 1
  ) +

  # Tercer cuartil
  geom_hline(
    yintercept = Q3,
    color = "#2A9D8F",
    linetype = "dotted",
    linewidth = 1
  ) +

  # Mediana
  geom_hline(
    yintercept = median(x),
    color = "#235B91",
    linewidth = 1.1
  ) +

  # Colores de los pacientes
  scale_color_manual(
    values = c(
      "FALSE" = "#28678B",
      "TRUE" = "#D64F4F"
    ),
    labels = c(
      "FALSE" = "Observación regular",
      "TRUE" = "Posible outlier"
    ),
    name = "Clasificación IQR"
  ) +

  # Etiquetas y diseño
  labs(
    title = "Detección de valores atípicos mediante IQR",
    subtitle = "Reducción tumoral en pacientes tratados con BTX-101 (datos simulados)",
    x = "Paciente",
    y = "Reducción del tamaño tumoral (%)",
    caption = "Líneas rojas: Q1 - 1.5 × IQR y Q3 + 1.5 × IQR"
  ) +

  theme_minimal(base_size = 14) +

  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "gray35"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

# Mostrar gráfica
print(g_IQR)

# Guardar gráfica
ggsave(
  "resultados_oncologia/02_limites_IQR.png",
  plot = g_IQR,
  width = 11,
  height = 6,
  dpi = 300
)



# 4. PRUEBA DE GRUBBS (UN solo posible outlier; requiere normalidad subyacente)
# No utilizar su valor p como permiso automatico de eliminacion.
cat("\n--- GRUBBS ---\n")
grubbs_resultado <- outliers::grubbs.test(x, type = 10, two.sided = TRUE)
print(grubbs_resultado)


###############################################################################
################################################################################
################################################################################



# 5. NORMALIDAD: HISTOGRAMA, Q-Q, SHAPIRO-WILK
# La normalidad se evalua para el desenlace continuo, NO para respuesta binaria.
g2 <- ggplot(datos, aes(x = reduccion_tumoral_pct)) +
  geom_histogram(aes(y = after_stat(density)), bins = 10, fill="#4D98BE", color="white", alpha=.85) +
  stat_function(fun=dnorm, args=list(mean=mean(x), sd=sd(x)), color="#C95B54", linewidth=1.2) +
  geom_vline(xintercept=mu0_t, linetype=2, color="#444444") +
  labs(title="Histograma y curva normal ajustada", subtitle="La curva es descriptiva; no prueba la normalidad", x="Reduccion tumoral (%)", y="Densidad") +
  theme_minimal(base_size=14)
print(g2)
ggsave("resultados_oncologia/02_histograma_normal.png",g2,width=9,height=6,dpi=300)

g3 <- ggplot(datos, aes(sample = reduccion_tumoral_pct)) +
  stat_qq(color="#276E98", size=2.5) +
  stat_qq_line(color="#C95B54", linewidth=1.2) +
  labs(title="Q-Q plot de reduccion tumoral", x="Cuantiles teoricos normales", y="Cuantiles observados (%)") +
  theme_minimal(base_size=14)
print(g3)
ggsave("resultados_oncologia/03_QQ_plot.png",g3,width=8,height=6,dpi=300)

cat("\n--- SHAPIRO-WILK ---\n")
shapiro_resultado <- shapiro.test(x)
print(shapiro_resultado)
cat("El valor p muestra normalidad; considerar Q-Q, extremos y diseno.\n")

# 6. INDEPENDENCIA: depende del diseno (pacientes distintos), no de Shapiro.
# Solo exploracion del orden de registro: no prueba independencia.
g4 <- ggplot(datos,aes(x=paciente,y=reduccion_tumoral_pct)) +
  geom_line(color="#93BAC7") + geom_point(color="#28678B",size=2) +
  geom_hline(yintercept=mu0_t,linetype=2,color="#B45E36") +
  labs(title="Reduccion tumoral por orden de registro", subtitle="Exploratorio: no demuestra independencia",x="Paciente (orden de registro)",y="Reduccion tumoral (%)") +
  theme_minimal(base_size=14)
print(g4)
ggsave("resultados_oncologia/04_orden_registro.png",g4,width=9,height=5,dpi=300)


# ==========================================================
# 7. PRUEBA t DE STUDENT PARA UNA MUESTRA
# ==========================================================

cat("\n--- PRUEBA t DE UNA MUESTRA ---\n")

t_resultado <- t.test(
  x,
  mu = mu0_t,
  alternative = "two.sided",
  conf.level = 0.95
)

print(t_resultado)

# Calculo manual
t_manual <- (mean(x) - mu0_t) /
  (sd(x) / sqrt(n))

p_manual <- 2 * pt(
  -abs(t_manual),
  df = n - 1
)

# Diferencias
diferencia_pp <- mean(x) - mu0_t

cambio_relativo <- (
  (mean(x) - mu0_t) / mu0_t
) * 100

# Mostrar resultados
cat("\nMedia observada:", round(mean(x), 4), "%\n")

cat("Media hipotetica:", mu0_t, "%\n")

cat("Diferencia en puntos porcentuales:",
    round(diferencia_pp, 4), "\n")

cat("Cambio porcentual relativo:",
    round(cambio_relativo, 2), "%\n")

cat("Estadistico t:",
    round(t_manual, 4), "\n")

cat("Grados de libertad:", n - 1, "\n")

cat("Valor p:",
    format.pval(p_manual, digits = 4), "\n")


# ==========================================================
# 8. TAMANO DEL EFECTO
# d DE COHEN Y g DE HEDGES
# ==========================================================

cat("\n--- TAMANO DEL EFECTO ---\n")

# Calculo con effectsize
resultado_cohen <- effectsize::cohens_d(
  x,
  mu = mu0_t,
  ci = 0.95
)

resultado_hedges <- effectsize::hedges_g(
  x,
  mu = mu0_t,
  ci = 0.95
)

# Mostrar resultados completos
print(resultado_cohen)
print(resultado_hedges)

# Extraer estimaciones e intervalos
d_cohen <- resultado_cohen$Cohens_d
g_hedges <- resultado_hedges$Hedges_g

cat("\n--- RESUMEN DE TAMANOS DEL EFECTO ---\n")

cat(
  "d de Cohen:",
  round(d_cohen, 4), "\n"
)

cat(
  "IC95% de Cohen: [",
  round(resultado_cohen$CI_low, 4), ", ",
  round(resultado_cohen$CI_high, 4), "]\n",
  sep = ""
)

cat(
  "g de Hedges:",
  round(g_hedges, 4), "\n"
)

cat(
  "IC95% de Hedges: [",
  round(resultado_hedges$CI_low, 4), ", ",
  round(resultado_hedges$CI_high, 4), "]\n",
  sep = ""
)

# Interpretacion orientativa del efecto
interpretar_efecto <- function(valor) {

  magnitud <- abs(valor)

  if (magnitud < 0.2) {
    return("Muy pequeno")
  } else if (magnitud < 0.5) {
    return("Pequeno")
  } else if (magnitud < 0.8) {
    return("Moderado")
  } else {
    return("Grande")
  }
}

cat(
  "\nMagnitud orientativa de Cohen:",
  interpretar_efecto(d_cohen), "\n"
)

cat(
  "Magnitud orientativa de Hedges:",
  interpretar_efecto(g_hedges), "\n"
)



# 9. GGHISTOSTATS: FRECUENTISTA Y BAYESIANA
# gghistostats usa x=, NO y=.
g5 <- ggstatsplot::gghistostats(
  data=datos, x=reduccion_tumoral_pct, type="parametric",
  test.value=mu0_t, alternative="two.sided", effsize.type="g",
  results.subtitle=TRUE, bf.message=FALSE,
  title="Prueba t de una muestra: reduccion tumoral",
  xlab="Reduccion tumoral (%)", ggtheme=theme_minimal(base_size=14)
) + labs(caption=NULL)
print(g5)
ggsave("resultados_oncologia/05_gghistostats_t.png",g5,width=12,height=8,dpi=300)

g6 <- ggstatsplot::gghistostats(
  data=datos, x=reduccion_tumoral_pct, type="bayes",
  test.value=mu0_t, alternative="two.sided", bf.prior=0.707,
  title="Prueba t bayesiana: reduccion tumoral",
  xlab="Reduccion tumoral (%)", ggtheme=theme_minimal(base_size=14)
) + labs(caption=NULL)
print(g6)
ggsave("resultados_oncologia/06_gghistostats_bayes.png",g6,width=12,height=8,dpi=300)
cat("\n--- BAYES FACTOR ---\n")
bf <- BayesFactor::ttestBF(x=x,mu=mu0_t,rscale=0.707)
print(bf)
cat("BF10:", BayesFactor::extractBF(bf, logbf=FALSE)$bf,"\n")










###############################################################################
################################################################################
################################################################################
###############################################################################
################################################################################
################################################################################
###############################################################################
################################################################################
################################################################################
###############################################################################
################################################################################
################################################################################



# 10. Z PARA UNA MEDIA: reduccion de biomarcador (ng/mL)
# sigma=8 es conocida por supuesto desde validacion EXTERNA hipotetica.
cat("\n--- Z PARA MEDIA DE BIOMARCADOR ---\n")
y <- datos$reduccion_biomarcador_ng_ml
z_media <- (mean(y)-mu0_z)/(sigma_conocida/sqrt(n))
p_z_media <- 2*pnorm(-abs(z_media))
IC_z_media <- mean(y)+c(-1,1)*qnorm(.975)*sigma_conocida/sqrt(n)
cat("media:",mean(y)," sigma conocida:",sigma_conocida," Z:",z_media," p:",p_z_media," IC95%:",IC_z_media,"\n")
# Si sigma no se conociera, no se justifica Z clasica; usar t.
g7 <- ggplot(datos,aes(x=reduccion_biomarcador_ng_ml)) +
  geom_histogram(bins=10,fill="#5A9A8E",color="white",alpha=.85) +
  geom_vline(xintercept=mu0_z,color="#CB574F",linetype=2,linewidth=1) +
  geom_vline(xintercept=mean(y),color="#246A94",linewidth=1) +
  labs(title="Prueba Z: reduccion del biomarcador",subtitle="Rojo: referencia; azul: media observada",x="Reduccion del biomarcador (ng/mL)",y="Pacientes") +
  theme_minimal(base_size=14)
print(g7)
ggsave("resultados_oncologia/07_Z_media_biomarcador.png",g7,width=9,height=6,dpi=300)

# 11. Z PARA UNA PROPORCION (RESPONDEDORES)
cat("\n--- Z PARA PROPORCION DE RESPUESTA ---\n")
exitos <- sum(datos$respuesta)
p_hat <- exitos/n
esperados <- c(exitos=n*p0, fracasos=n*(1-p0))
cat("Respuestas:",exitos," n:",n," proporcion:",p_hat," esperados H0:",esperados,"\n")
if (any(esperados < 10)) warning("Aproximacion normal cuestionable: consultar binom.test")
z_prop <- (p_hat-p0)/sqrt(p0*(1-p0)/n)
p_z_prop <- 2*pnorm(-abs(z_prop))
cat("Z proporcion:",z_prop," p bilateral:",p_z_prop,"\n")
cat("prop.test sin correccion de continuidad:\n")
print(prop.test(exitos,n,p=p0,alternative="two.sided",correct=FALSE))
cat("Prueba binomial exacta de contraste:\n")
print(binom.test(exitos,n,p=p0,alternative="two.sided"))
# Intervalo Wilson aproximado generado por prop.test
ic_wilson <- prop.test(exitos,n,correct=FALSE)$conf.int

g8 <- ggplot(data.frame(criterio="Respuesta", propor= p_hat),aes(x=criterio,y=propor)) +
  geom_col(width=.45,fill="#397EA1",alpha=.85) +
  geom_errorbar(aes(ymin=ic_wilson[1],ymax=ic_wilson[2]),width=.1,linewidth=1) +
  geom_hline(yintercept=p0,linetype=2,color="#CB574F",linewidth=1) +
  scale_y_continuous(labels=function(z) paste0(round(100*z),"%"),limits=c(0,1)) +
  labs(title="Proporcion de pacientes con respuesta",subtitle="Barra: proporción observada; intervalo Wilson 95%; rojo: 50%",x=NULL,y="Proporcion de respondedores") +
  theme_minimal(base_size=14)
print(g8)
ggsave("resultados_oncologia/08_Z_proporcion_respuesta.png",g8,width=8,height=6,dpi=300)

# 12. EXPORTAR TABLA GENERAL DE RESULTADOS
resumen <- data.frame(
  prueba=c("t: reduccion tumoral", "Z: biomarcador", "Z: proporcion respondedores"),
  n=c(n,n,n),
  estimacion=c(mean(x),mean(y),p_hat),
  referencia=c(mu0_t,mu0_z,p0),
  estadistico=c(unname(t_resultado$statistic),z_media,z_prop),
  p_valor=c(t_resultado$p.value,p_z_media,p_z_prop)
)
print(resumen)
write.csv(resumen,"resultados_oncologia/resumen_t_Z_Z.csv",row.names=FALSE)
cat("\nListo. Graficas, datos y resumen en: resultados_oncologia/\n")












# ==========================================================
# 12. TABLA FINAL COMPLETA
# t DE STUDENT, Z Y ANALISIS BAYESIANO
# ==========================================================

library(effectsize)
library(BayesFactor)

# ----------------------------------------------------------
# 1. TAMANOS DEL EFECTO
# ----------------------------------------------------------

cohen <- effectsize::cohens_d(
  x,
  mu = mu0_t,
  ci = 0.95
)

hedges <- effectsize::hedges_g(
  x,
  mu = mu0_t,
  ci = 0.95
)

# Extraer valores
d_cohen <- cohen$Cohens_d
g_hedges <- hedges$Hedges_g

d_IC_inf <- cohen$CI_low
d_IC_sup <- cohen$CI_high

g_IC_inf <- hedges$CI_low
g_IC_sup <- hedges$CI_high


# ----------------------------------------------------------
# 2. PRUEBA t FRECUENTISTA
# ----------------------------------------------------------

t_resultado <- t.test(
  x,
  mu = mu0_t,
  alternative = "two.sided",
  conf.level = 0.95
)


# ----------------------------------------------------------
# 3. PRUEBA t BAYESIANA
# ----------------------------------------------------------

bayes <- BayesFactor::ttestBF(
  x = x,
  mu = mu0_t,
  rscale = 0.707
)

# Extraer BF10
BF10 <- as.vector(
  extractBF(bayes)$bf
)[1]

# BF01
BF01 <- 1 / BF10

# Logaritmos naturales
logBF10 <- log(BF10)
logBF01 <- log(BF01)


# ----------------------------------------------------------
# 4. AGREGAR COLUMNAS AL RESUMEN EXISTENTE
# ----------------------------------------------------------

resumen$IC95_inf <- NA_real_
resumen$IC95_sup <- NA_real_

resumen$Cohens_d <- NA_real_
resumen$Cohens_d_IC_inf <- NA_real_
resumen$Cohens_d_IC_sup <- NA_real_

resumen$Hedges_g <- NA_real_
resumen$Hedges_g_IC_inf <- NA_real_
resumen$Hedges_g_IC_sup <- NA_real_

resumen$BF10 <- NA_real_
resumen$BF01 <- NA_real_
resumen$logBF10 <- NA_real_
resumen$logBF01 <- NA_real_

# Identificar fila de prueba t
fila_t <- which(
  resumen$prueba == "t: reduccion tumoral"
)

stopifnot(length(fila_t) == 1)

# Intervalo de confianza de la media
resumen$IC95_inf[fila_t] <- t_resultado$conf.int[1]
resumen$IC95_sup[fila_t] <- t_resultado$conf.int[2]

# Cohen
resumen$Cohens_d[fila_t] <- d_cohen
resumen$Cohens_d_IC_inf[fila_t] <- d_IC_inf
resumen$Cohens_d_IC_sup[fila_t] <- d_IC_sup

# Hedges
resumen$Hedges_g[fila_t] <- g_hedges
resumen$Hedges_g_IC_inf[fila_t] <- g_IC_inf
resumen$Hedges_g_IC_sup[fila_t] <- g_IC_sup

# Bayes
resumen$BF10[fila_t] <- BF10
resumen$BF01[fila_t] <- BF01
resumen$logBF10[fila_t] <- logBF10
resumen$logBF01[fila_t] <- logBF01


# ----------------------------------------------------------
# 5. IMPRIMIR RESULTADOS
# ----------------------------------------------------------

cat("\n========== TABLA FINAL ==========\n")

print(resumen, row.names = FALSE)

cat("\n--- TAMANOS DEL EFECTO ---\n")
cat("Cohen d:", round(d_cohen, 4), "\n")
cat("Hedges g:", round(g_hedges, 4), "\n")

cat("\n--- ANALISIS BAYESIANO ---\n")
cat("BF10:", round(BF10, 4), "\n")
cat("BF01:", round(BF01, 4), "\n")
cat("log(BF10):", round(logBF10, 4), "\n")
cat("log(BF01):", round(logBF01, 4), "\n")


# ----------------------------------------------------------
# 6. EXPORTAR CSV COMPLETO
# ----------------------------------------------------------

write.csv(
  resumen,
  "resultados_oncologia/resumen_t_Z_Z_completo.csv",
  row.names = FALSE,
  na = ""
)

cat("\nArchivo exportado correctamente.\n")











###############################################################################
################################################################################
################################################################################
###############################################################################
################################################################################
################################################################################
###############################################################################
################################################################################
################################################################################
###############################################################################
################################################################################
################################################################################
###############################################################################
################################################################################
################################################################################
###############################################################################
################################################################################
################################################################################
###############################################################################
################################################################################
################################################################################
###############################################################################
################################################################################
################################################################################
###############################################################################
################################################################################
################################################################################
###############################################################################
################################################################################
################################################################################
###############################################################################
################################################################################
################################################################################
###############################################################################
################################################################################
################################################################################
###############################################################################
################################################################################
################################################################################
###############################################################################
################################################################################
################################################################################
###############################################################################
################################################################################
################################################################################
###############################################################################
################################################################################
################################################################################














# =====================================================
# DETECCIÓN DE VALORES ATÍPICOS
# PRUEBA Z DE UNA MUESTRA
# =====================================================

# Seleccionar variable
y <- datos$reduccion_biomarcador_ng_ml

# 1. Diagrama de caja
boxplot(
  y,
  main = "Valores atípicos del biomarcador BTX-M1",
  ylab = "Reducción del biomarcador (ng/mL)",
  col = "lightblue"
)

# 2. Identificar outliers mediante IQR
boxplot.stats(y)$out

# 3. Prueba de Grubbs
library(outliers)

grubbs.test(
  y,
  type = 10,
  two.sided = TRUE
)


# =====================================================
# BOXPLOT CON JITTER
# Biomarcador BTX-M1
# =====================================================

library(ggplot2)

# Seleccionar variable
y <- datos$reduccion_biomarcador_ng_ml

# Crear base para graficar
df <- data.frame(
  Grupo = "BTX-101",
  Biomarcador = y
)

# Boxplot + jitter
ggplot(df, aes(x = Grupo, y = Biomarcador)) +
  geom_boxplot(
    width = 0.35,
    fill = "#A8DADC",
    color = "#1D3557",
    outlier.shape = NA,
    alpha = 0.7
  ) +
  geom_jitter(
    width = 0.12,
    size = 2.5,
    alpha = 0.7,
    color = "#E63946"
  ) +
  labs(
    title = "Distribución de la reducción del biomarcador",
    subtitle = "Tratamiento experimental BTX-101",
    x = "",
    y = "Reducción del biomarcador (ng/mL)"
  ) +
  theme_classic(base_size = 14)







# =====================================================
# GRÁFICA DE IQR Y LÍMITES DE OUTLIERS
# =====================================================

library(ggplot2)

# Calcular cuartiles y límites
Q1 <- quantile(y, 0.25, na.rm = TRUE)
Q3 <- quantile(y, 0.75, na.rm = TRUE)

IQR_val <- IQR(y, na.rm = TRUE)

LI <- Q1 - 1.5 * IQR_val
LS <- Q3 + 1.5 * IQR_val

# Clasificar observaciones
df$Atipico <- ifelse(
  df$Biomarcador < LI | df$Biomarcador > LS,
  "Valor atípico",
  "Valor dentro de límites"
)

# Graficar
ggplot(df, aes(x = Grupo, y = Biomarcador)) +

  # Región intercuartílica
  annotate(
    "rect",
    xmin = 0.65, xmax = 1.35,
    ymin = as.numeric(Q1),
    ymax = as.numeric(Q3),
    fill = "#A8DADC",
    alpha = 0.4
  ) +

  # Límites para detectar outliers
  geom_hline(
    yintercept = c(LI, LS),
    color = "#E63946",
    linetype = "dashed",
    linewidth = 0.9
  ) +

  # Cuartiles
  geom_hline(
    yintercept = c(Q1, Q3),
    color = "#1D3557",
    linewidth = 0.8
  ) +

  # Mediana
  geom_hline(
    yintercept = median(y, na.rm = TRUE),
    color = "#2A9D8F",
    linewidth = 1
  ) +

  # Puntos individuales
  geom_jitter(
    aes(color = Atipico),
    width = 0.12,
    size = 2.5,
    alpha = 0.8
  ) +

  scale_color_manual(
    values = c(
      "Valor atípico" = "#E63946",
      "Valor dentro de límites" = "#457B9D"
    )
  ) +

  labs(
    title = "Identificación de valores atípicos mediante IQR",
    subtitle = "Reducción del biomarcador BTX-M1",
    x = "",
    y = "Reducción del biomarcador (ng/mL)",
    color = "Clasificación"
  ) +

  theme_classic(base_size = 14)






cat("Q1 =", Q1, "\n")
cat("Q3 =", Q3, "\n")
cat("IQR =", IQR_val, "\n")
cat("Límite inferior =", LI, "\n")
cat("Límite superior =", LS, "\n")

cat("Número de outliers =",
    sum(df$Atipico == "Valor atípico", na.rm = TRUE),
    "\n")

# Mostrar los pacientes con valores atípicos
datos[
  !is.na(datos$reduccion_biomarcador_ng_ml) &
  (datos$reduccion_biomarcador_ng_ml < LI |
   datos$reduccion_biomarcador_ng_ml > LS),
]




# =====================================================
# GENERAR Y GUARDAR BOXPLOT + JITTER E IQR
# BIOMARCADOR BTX-M1
# =====================================================

library(ggplot2)

# 1. Preparar los datos
y <- datos$reduccion_biomarcador_ng_ml

df <- data.frame(
  Grupo = "BTX-101",
  Biomarcador = y
)

df <- na.omit(df)

# 2. Calcular IQR y límites
Q1 <- quantile(y, 0.25, na.rm = TRUE)
Q3 <- quantile(y, 0.75, na.rm = TRUE)

IQR_val <- IQR(y, na.rm = TRUE)

LI <- Q1 - 1.5 * IQR_val
LS <- Q3 + 1.5 * IQR_val

df$Atipico <- ifelse(
  df$Biomarcador < LI | df$Biomarcador > LS,
  "Valor atípico",
  "Valor dentro de límites"
)

# =====================================================
# GRÁFICA 1: BOXPLOT CON JITTER
# =====================================================

grafica_boxplot <- ggplot(
  df,
  aes(x = Grupo, y = Biomarcador)
) +
  geom_boxplot(
    width = 0.35,
    fill = "#A8DADC",
    color = "#1D3557",
    outlier.shape = NA,
    alpha = 0.7
  ) +
  geom_jitter(
    width = 0.12,
    size = 2.5,
    alpha = 0.7,
    color = "#E63946"
  ) +
  labs(
    title = "Distribución de la reducción del biomarcador",
    subtitle = "Tratamiento experimental BTX-101",
    x = "",
    y = "Reducción del biomarcador (ng/mL)"
  ) +
  theme_classic(base_size = 14)

# Mostrar gráfica 1
print(grafica_boxplot)


# =====================================================
# GRÁFICA 2: IQR Y LÍMITES DE OUTLIERS
# =====================================================

grafica_iqr <- ggplot(
  df,
  aes(x = Grupo, y = Biomarcador)
) +

  # Región intercuartílica
  annotate(
    "rect",
    xmin = 0.65,
    xmax = 1.35,
    ymin = as.numeric(Q1),
    ymax = as.numeric(Q3),
    fill = "#A8DADC",
    alpha = 0.4
  ) +

  # Límites inferior y superior
  geom_hline(
    yintercept = c(LI, LS),
    color = "#E63946",
    linetype = "dashed",
    linewidth = 0.9
  ) +

  # Primer y tercer cuartil
  geom_hline(
    yintercept = c(Q1, Q3),
    color = "#1D3557",
    linewidth = 0.8
  ) +

  # Mediana
  geom_hline(
    yintercept = median(y, na.rm = TRUE),
    color = "#2A9D8F",
    linewidth = 1
  ) +

  # Valores individuales
  geom_jitter(
    aes(color = Atipico),
    width = 0.12,
    size = 2.5,
    alpha = 0.8
  ) +

  scale_color_manual(
    values = c(
      "Valor atípico" = "#E63946",
      "Valor dentro de límites" = "#457B9D"
    )
  ) +

  labs(
    title = "Identificación de valores atípicos mediante IQR",
    subtitle = "Reducción del biomarcador BTX-M1",
    x = "",
    y = "Reducción del biomarcador (ng/mL)",
    color = "Clasificación"
  ) +

  theme_classic(base_size = 14)

# Mostrar gráfica 2
print(grafica_iqr)


# =====================================================
# GUARDAR AMBAS GRÁFICAS EN PNG
# =====================================================

dir.create(
  "Graficas_BTX101",
  showWarnings = FALSE
)

ggsave(
  "Graficas_BTX101/Boxplot_Jitter_BTX101_Z.png",
  plot = grafica_boxplot,
  width = 8,
  height = 6,
  dpi = 300,
  bg = "white"
)

ggsave(
  "Graficas_BTX101/IQR_Outliers_BTX101_Z.png",
  plot = grafica_iqr,
  width = 8,
  height = 6,
  dpi = 300,
  bg = "white"
)

# Mostrar ubicación de los archivos
normalizePath("Graficas_BTX101")



y <- datos$reduccion_biomarcador_ng_ml

grubbs_resultado <- outliers::grubbs.test(
  y,
  type = 10,
  two.sided = TRUE
)

print(grubbs_resultado)

# =====================================================
# EVALUACIÓN DEL SUPUESTO DE NORMALIDAD
# PRUEBA Z DE UNA MUESTRA
# Biomarcador BTX-M1
# =====================================================

library(ggplot2)
library(nortest)

# =====================================================
# 1. SELECCIONAR VARIABLE
# =====================================================

y <- datos$reduccion_biomarcador_ng_ml

# Eliminar valores faltantes
y <- y[!is.na(y)]

# Estadística descriptiva
cat("Número de pacientes =", length(y), "\n")
cat("Media =", mean(y), "\n")
cat("Mediana =", median(y), "\n")
cat("Desviación estándar =", sd(y), "\n")


# =====================================================
# 2. PRUEBA DE SHAPIRO-WILK
# =====================================================

shapiro_resultado <- shapiro.test(y)

print(shapiro_resultado)


# =====================================================
# 3. PRUEBA DE LILLIEFORS
# Kolmogorov-Smirnov con parámetros estimados
# =====================================================

lillie_resultado <- nortest::lillie.test(y)

print(lillie_resultado)


# =====================================================
# 4. PREPARAR DATOS PARA GRAFICAR
# =====================================================

df_normalidad <- data.frame(
  Biomarcador = y
)


# =====================================================
# 5. HISTOGRAMA CON CURVA NORMAL TEÓRICA
# =====================================================

grafica_histograma <- ggplot(
  df_normalidad,
  aes(x = Biomarcador)
) +

  geom_histogram(
    aes(y = after_stat(density)),
    bins = 12,
    fill = "#A8DADC",
    color = "#1D3557",
    alpha = 0.8
  ) +

  stat_function(
    fun = dnorm,
    args = list(
      mean = mean(y),
      sd = sd(y)
    ),
    color = "#E63946",
    linewidth = 1.2
  ) +

  labs(
    title = "Distribución del biomarcador BTX-M1",
    subtitle = "Histograma y curva normal de referencia",
    x = "Reducción del biomarcador (ng/mL)",
    y = "Densidad"
  ) +

  theme_classic(base_size = 14)

print(grafica_histograma)


# =====================================================
# 6. GRÁFICO Q-Q PLOT
# =====================================================

grafica_qq <- ggplot(
  df_normalidad,
  aes(sample = Biomarcador)
) +

  stat_qq(
    color = "#457B9D",
    size = 2.5,
    alpha = 0.8
  ) +

  stat_qq_line(
    color = "#E63946",
    linewidth = 1
  ) +

  labs(
    title = "Normalidad del biomarcador BTX-M1",
    subtitle = "Gráfico Q-Q normal",
    x = "Cuantiles teóricos",
    y = "Cuantiles observados (ng/mL)"
  ) +

  theme_classic(base_size = 14)

print(grafica_qq)


# =====================================================
# 7. TABLA RESUMEN DE NORMALIDAD
# =====================================================

tabla_normalidad <- data.frame(
  Prueba = c("Shapiro-Wilk", "Lilliefors"),
  Estadistico = c(
    unname(shapiro_resultado$statistic),
    unname(lillie_resultado$statistic)
  ),
  Valor_p = c(
    shapiro_resultado$p.value,
    lillie_resultado$p.value
  )
)

tabla_normalidad$Interpretacion <- ifelse(
  tabla_normalidad$Valor_p < 0.05,
  "Se rechaza normalidad",
  "No se rechaza normalidad"
)

print(tabla_normalidad)


# =====================================================
# 8. GUARDAR GRÁFICAS A 300 DPI
# =====================================================

dir.create(
  "Graficas_BTX101",
  showWarnings = FALSE
)

ggsave(
  "Graficas_BTX101/Histograma_Normalidad_Z.png",
  plot = grafica_histograma,
  width = 8,
  height = 6,
  dpi = 300,
  bg = "white"
)

ggsave(
  "Graficas_BTX101/QQplot_Normalidad_Z.png",
  plot = grafica_qq,
  width = 8,
  height = 6,
  dpi = 300,
  bg = "white"
)

# Guardar tabla
write.csv(
  tabla_normalidad,
  "Graficas_BTX101/Resultados_Normalidad_Z.csv",
  row.names = FALSE
)

cat("\nAnálisis de normalidad finalizado.\n")






# =====================================================
# PRUEBA Z DE UNA MUESTRA
# Biomarcador BTX-M1
# =====================================================

# Instalar una sola vez (si no lo tienes)
# install.packages("BSDA")

library(BSDA)

# Seleccionar variable
y <- datos$reduccion_biomarcador_ng_ml

# Valor de referencia
mu0_z <- 15

# Desviación estándar poblacional conocida
sigma <- 8

# Ejecutar prueba Z bilateral
cat("\n--- PRUEBA Z DE UNA MUESTRA ---\n")

z_resultado <- BSDA::z.test(
  y,
  mu = mu0_z,
  sigma.x = sigma,
  alternative = "two.sided",
  conf.level = 0.95
)

print(z_resultado)



# =====================================================
# TAMAÑO DEL EFECTO
# d DE COHEN Y g DE HEDGES
# Biomarcador BTX-M1
# =====================================================

library(effectsize)

# Seleccionar variable
y <- datos$reduccion_biomarcador_ng_ml
y <- na.omit(y)

# Media de referencia
mu0_z <- 15

# Desviación estándar poblacional conocida
sigma <- 8

# =====================================================
# 1. d DE COHEN (desviación estándar muestral)
# =====================================================

cohen <- effectsize::cohens_d(
  y,
  mu = mu0_z,
  ci = 0.95
)

print(cohen)

# =====================================================
# 2. g DE HEDGES
# =====================================================

hedges <- effectsize::hedges_g(
  y,
  mu = mu0_z,
  ci = 0.95
)

print(hedges)

# =====================================================
# 3. EFECTO ESTANDARIZADO CON SIGMA CONOCIDA
# =====================================================

d_sigma <- (mean(y) - mu0_z) / sigma

cat("\n--- TAMAÑOS DEL EFECTO ---\n")

cat("d de Cohen =", round(cohen$Cohens_d, 4), "\n")

cat("g de Hedges =", round(hedges$Hedges_g, 4), "\n")

cat("Efecto estandarizado con sigma conocida =",
    round(d_sigma, 4), "\n")
    
    
    
    
    
    
    
    # =====================================================
# PRUEBA Z DE UNA MUESTRA
# GRÁFICA ESTILO ggstatsplot
# Biomarcador BTX-M1
# =====================================================

library(ggplot2)

# 1. Seleccionar variable
y <- datos$reduccion_biomarcador_ng_ml
y <- y[is.finite(y)]

# 2. Parámetros
mu0 <- 15
sigma <- 8    # Valor poblacional supuesto
n <- length(y)

# 3. Estadística descriptiva
media <- mean(y)
DE <- sd(y)

# 4. Prueba Z bilateral
EE <- sigma / sqrt(n)

Z <- (media - mu0) / EE

p_valor <- 2 * pnorm(-abs(Z))

# 5. Tamaño del efecto estandarizado con sigma
d_sigma <- (media - mu0) / sigma

# 6. Crear subtítulo estadístico
subtitulo <- sprintf(
  "Z = %.2f, p = %.3g, d_sigma = %.2f, n = %d",
  Z, p_valor, d_sigma, n
)

# 7. Preparar base de datos
df_z <- data.frame(
  Biomarcador = y
)

# =====================================================
# 8. GRÁFICA ESTILO ggstatsplot
# =====================================================

grafica_z <- ggplot(
  df_z,
  aes(x = Biomarcador)
) +

  geom_histogram(
    bins = 10,
    fill = "grey65",
    color = "black",
    alpha = 0.85
  ) +

  # Línea azul: media observada
  geom_vline(
    xintercept = media,
    color = "blue",
    linetype = "dashed",
    linewidth = 1
  ) +

  # Etiqueta de la media
  annotate(
    "text",
    x = media,
    y = Inf,
    label = paste0(
      "Media = ", round(media, 2)
    ),
    vjust = 1.5,
    hjust = -0.1,
    size = 3.5
  ) +

  labs(
    title = "Prueba Z de una muestra: reducción del biomarcador",
    subtitle = subtitulo,
    x = "Reducción del biomarcador (ng/mL)",
    y = "Frecuencia",
    caption = "BTX-101 | sigma poblacional supuesta = 8 ng/mL"
  ) +

  theme_bw(base_size = 12) +

  theme(
    plot.title = element_text(size = 13),
    plot.subtitle = element_text(size = 11),
    plot.caption = element_text(size = 9),
    panel.grid.minor = element_blank()
  )

# Mostrar gráfica
print(grafica_z)

# =====================================================
# 9. GUARDAR GRÁFICA
# =====================================================

dir.create(
  "Graficas_BTX101",
  showWarnings = FALSE
)

ggsave(
  "Graficas_BTX101/Prueba_Z_BTX101_ggstatsplot_estilo.png",
  plot = grafica_z,
  width = 9,
  height = 6,
  dpi = 300,
  bg = "white"
)









# =====================================================
# ANÁLISIS BAYESIANO PARA UNA MEDIA
# EQUIVALENTE BAYESIANO DE LA PRUEBA Z
# Biomarcador BTX-M1 | BTX-101
# =====================================================

library(ggplot2)

# =====================================================
# 1. SELECCIONAR VARIABLE
# =====================================================

y <- datos$reduccion_biomarcador_ng_ml
y <- y[is.finite(y)]

n <- length(y)
media <- mean(y)

# =====================================================
# 2. DEFINIR LOS PARÁMETROS
# =====================================================

# Hipótesis nula: mu = 15
mu0 <- 15

# Desviación estándar poblacional conocida
# Supuesto didáctico
sigma <- 8

# Distribución previa bajo H1:
# mu ~ Normal(mu0, tau^2)
#
# tau = 8 ng/mL, elección didáctica
# Representa la incertidumbre previa sobre mu

tau <- 8

# =====================================================
# 3. FACTOR DE BAYES
# =====================================================

# Error estándar
EE <- sigma / sqrt(n)

# Densidad marginal de la media bajo H0
m0 <- dnorm(
  media,
  mean = mu0,
  sd = EE
)

# Densidad marginal de la media bajo H1
m1 <- dnorm(
  media,
  mean = mu0,
  sd = sqrt(EE^2 + tau^2)
)

# Factor de Bayes
BF10 <- m1 / m0
BF01 <- 1 / BF10

cat("\n--- PRUEBA BAYESIANA ---\n")

cat("n =", n, "\n")
cat("Media observada =", round(media, 4), "\n")
cat("Sigma poblacional =", sigma, "\n")
cat("Escala previa tau =", tau, "\n")

cat("\nBF10 =", round(BF10, 4), "\n")
cat("BF01 =", round(BF01, 6), "\n")

cat("ln(BF10) =", round(log(BF10), 4), "\n")

# =====================================================
# 4. DISTRIBUCIÓN POSTERIOR DE LA MEDIA
# =====================================================

# Varianza posterior
var_post <- 1 / (
  n / sigma^2 + 1 / tau^2
)

# Desviación estándar posterior
sd_post <- sqrt(var_post)

# Media posterior
media_post <- var_post * (
  n * media / sigma^2 +
  mu0 / tau^2
)

# Intervalo creíble del 95 %
IC_post <- qnorm(
  c(0.025, 0.975),
  mean = media_post,
  sd = sd_post
)

cat("\n--- ESTIMACIÓN POSTERIOR ---\n")

cat("Media posterior =",
    round(media_post, 4), "\n")

cat("DE posterior =",
    round(sd_post, 4), "\n")

cat("IC creíble 95% = [",
    round(IC_post[1], 4), ",",
    round(IC_post[2], 4), "]\n")

# Probabilidad posterior de que mu supere 15
prob_superior <- 1 - pnorm(
  mu0,
  mean = media_post,
  sd = sd_post
)

cat("P(mu > 15 | datos, H1) =",
    round(prob_superior, 6), "\n")

# =====================================================
# 5. GRÁFICA DE LA DISTRIBUCIÓN POSTERIOR
# =====================================================

x_post <- seq(
  media_post - 4 * sd_post,
  media_post + 4 * sd_post,
  length.out = 1000
)

df_post <- data.frame(
  Media = x_post,
  Densidad = dnorm(
    x_post,
    mean = media_post,
    sd = sd_post
  )
)

grafica_bayes_z <- ggplot(
  df_post,
  aes(x = Media, y = Densidad)
) +

  geom_area(
    fill = "#A8DADC",
    alpha = 0.65
  ) +

  geom_line(
    color = "#1D3557",
    linewidth = 1.2
  ) +

  geom_vline(
    xintercept = mu0,
    color = "#E63946",
    linetype = "dashed",
    linewidth = 1
  ) +

  geom_vline(
    xintercept = media_post,
    color = "#2A9D8F",
    linewidth = 1
  ) +

  labs(
    title = "Análisis bayesiano: biomarcador BTX-M1",
    subtitle = paste0(
      "BF10 = ", round(BF10, 2),
      " | Media posterior = ",
      round(media_post, 2),
      " ng/mL"
    ),
    x = "Media de reducción del biomarcador (ng/mL)",
    y = "Densidad posterior",
    caption = paste0(
      "Línea roja: H0 = 15 ng/mL",
      " | Línea verde: media posterior"
    )
  ) +

  theme_classic(base_size = 14)

print(grafica_bayes_z)

# =====================================================
# 6. GUARDAR GRÁFICA
# =====================================================

dir.create(
  "Graficas_BTX101",
  showWarnings = FALSE
)

ggsave(
  "Graficas_BTX101/Bayesiana_Z_BTX101.png",
  plot = grafica_bayes_z,
  width = 9,
  height = 6,
  dpi = 300,
  bg = "white"
)







# =====================================================
# ANÁLISIS BAYESIANO MCMC
# BIOMARCADOR BTX-M1
# 4 CADENAS + CONVERGENCIA + GRÁFICAS
# =====================================================

library(ggplot2)
library(posterior)

set.seed(2026)

# =====================================================
# 1. DATOS Y PARÁMETROS
# =====================================================

y <- datos$reduccion_biomarcador_ng_ml
y <- y[is.finite(y)]

n <- length(y)
media <- mean(y)

mu0 <- 15
sigma <- 8    # DE poblacional conocida: supuesto
tau <- 8      # DE de la distribución previa

# =====================================================
# 2. SOLUCIÓN BAYESIANA ANALÍTICA
# =====================================================

var_post <- 1 / (
  n / sigma^2 + 1 / tau^2
)

sd_post <- sqrt(var_post)

media_post <- var_post * (
  n * media / sigma^2 +
  mu0 / tau^2
)

# Factor de Bayes: cálculo analítico
EE <- sigma / sqrt(n)

log_BF10 <- dnorm(
  media,
  mean = mu0,
  sd = sqrt(EE^2 + tau^2),
  log = TRUE
) -
  dnorm(
    media,
    mean = mu0,
    sd = EE,
    log = TRUE
  )

BF10 <- exp(log_BF10)

# =====================================================
# 3. DISTRIBUCIÓN POSTERIOR (LOG)
# =====================================================

log_posterior <- function(mu) {

  # Verosimilitud
  log_lik <- -n * (media - mu)^2 /
    (2 * sigma^2)

  # Distribución previa
  log_prior <- -(mu - mu0)^2 /
    (2 * tau^2)

  log_lik + log_prior
}

# =====================================================
# 4. FUNCIÓN METROPOLIS-HASTINGS
# =====================================================

metropolis <- function(inicio, iteraciones,
                       propuesta_sd) {

  cadena <- numeric(iteraciones)
  cadena[1] <- inicio

  aceptadas <- 0

  for (i in 2:iteraciones) {

    actual <- cadena[i - 1]

    propuesta <- rnorm(
      1,
      mean = actual,
      sd = propuesta_sd
    )

    log_alpha <- log_posterior(propuesta) -
      log_posterior(actual)

    if (log(runif(1)) < log_alpha) {

      cadena[i] <- propuesta
      aceptadas <- aceptadas + 1

    } else {

      cadena[i] <- actual
    }
  }

  list(
    muestras = cadena,
    aceptacion = aceptadas / (iteraciones - 1)
  )
}

# =====================================================
# 5. EJECUTAR 4 CADENAS
# =====================================================

n_cadenas <- 4
n_iter <- 6000
warmup <- 2000

# Diferentes puntos iniciales
inicios <- c(7, 15, 23, 31)

# Escala de propuestas
propuesta_sd <- 2 * sd_post

resultados <- lapply(
  seq_len(n_cadenas),
  function(i) {
    metropolis(
      inicio = inicios[i],
      iteraciones = n_iter,
      propuesta_sd = propuesta_sd
    )
  }
)

# =====================================================
# 6. CONSERVAR MUESTRAS POSTERIORES
# =====================================================

cadenas <- sapply(
  resultados,
  function(r) r$muestras
)

# Eliminar iteraciones de calentamiento
cadenas_post <- cadenas[
  (warmup + 1):n_iter,
  ,
  drop = FALSE
]

colnames(cadenas_post) <- paste0(
  "Cadena ", 1:n_cadenas
)

# Convertir al formato del paquete posterior
array_post <- array(
  cadenas_post,
  dim = c(
    nrow(cadenas_post),
    n_cadenas,
    1
  ),
  dimnames = list(
    NULL,
    paste0("Cadena_", 1:n_cadenas),
    "mu"
  )
)

draws <- posterior::as_draws_array(array_post)

# =====================================================
# 7. DIAGNÓSTICOS DE CONVERGENCIA
# =====================================================

diagnosticos <- posterior::summarise_draws(
  draws,
  mean,
  sd,
  rhat,
  ess_bulk,
  ess_tail
)

print(diagnosticos)

# Tasas de aceptación
aceptacion <- sapply(
  resultados,
  function(r) r$aceptacion
)

cat("\n--- TASAS DE ACEPTACIÓN ---\n")
print(round(aceptacion, 3))

# =====================================================
# 8. RESUMEN POSTERIOR
# =====================================================

muestras <- as.vector(cadenas_post)

IC_MCMC <- quantile(
  muestras,
  probs = c(0.025, 0.975)
)

cat("\n--- RESULTADOS BAYESIANOS ---\n")

cat("Pacientes =", n, "\n")

cat("Muestras posteriores =",
    length(muestras), "\n")

cat("BF10 analítico =",
    round(BF10, 4), "\n")

cat("Media posterior exacta =",
    round(media_post, 4), "\n")

cat("Media posterior MCMC =",
    round(mean(muestras), 4), "\n")

cat("IC creíble MCMC 95% =",
    round(IC_MCMC[1], 4), "a",
    round(IC_MCMC[2], 4), "\n")

# =====================================================
# 9. PREPARAR DATOS PARA GRÁFICAS
# =====================================================

df_mcmc <- data.frame(
  Iteracion = rep(
    seq_len(nrow(cadenas_post)),
    times = n_cadenas
  ),
  Mu = as.vector(cadenas_post),
  Cadena = factor(
    rep(1:n_cadenas, each = nrow(cadenas_post))
  )
)

# =====================================================
# 10. GRÁFICA DE TRAZAS (TRACE PLOT)
# =====================================================

grafica_trazas <- ggplot(
  df_mcmc,
  aes(x = Iteracion, y = Mu)
) +
  geom_line(
    color = "#457B9D",
    linewidth = 0.25,
    alpha = 0.8
  ) +
  facet_wrap(
    ~ Cadena,
    ncol = 2
  ) +
  geom_hline(
    yintercept = media_post,
    color = "#E63946",
    linewidth = 0.7
  ) +
  labs(
    title = "Convergencia MCMC: biomarcador BTX-M1",
    subtitle = "Cuatro cadenas de Metropolis-Hastings",
    x = "Iteración posterior",
    y = "Media poblacional (ng/mL)"
  ) +
  theme_classic(base_size = 13)

print(grafica_trazas)

# =====================================================
# 11. DISTRIBUCIONES POR CADENA
# =====================================================

grafica_densidades <- ggplot(
  df_mcmc,
  aes(x = Mu, color = Cadena)
) +
  geom_density(
    linewidth = 1
  ) +
  geom_vline(
    xintercept = mu0,
    color = "#E63946",
    linetype = "dashed"
  ) +
  labs(
    title = "Distribuciones posteriores por cadena",
    subtitle = paste0(
      "BF10 = ", round(BF10, 2),
      " | Media posterior = ",
      round(media_post, 2), " ng/mL"
    ),
    x = "Media del biomarcador (ng/mL)",
    y = "Densidad posterior"
  ) +
  theme_classic(base_size = 13)

print(grafica_densidades)

# =====================================================
# 12. GUARDAR RESULTADOS
# =====================================================

dir.create(
  "Graficas_BTX101",
  showWarnings = FALSE
)

ggsave(
  "Graficas_BTX101/MCMC_Traceplot_BTX101.png",
  plot = grafica_trazas,
  width = 10,
  height = 7,
  dpi = 300
)

ggsave(
  "Graficas_BTX101/MCMC_Densidades_BTX101.png",
  plot = grafica_densidades,
  width = 9,
  height = 6,
  dpi = 300
)

write.csv(
  as.data.frame(diagnosticos),
  "Graficas_BTX101/MCMC_Diagnosticos.csv",
  row.names = FALSE
)














##############################################################################
################################################################################
################################################################################
###############################################################################
################################################################################
################################################################################
###############################################################################
################################################################################
################################################################################
###############################################################################
################################################################################
################################################################################
###############################################################################
################################################################################
################################################################################
###############################################################################
################################################################################
################################################################################
###############################################################################
################################################################################
################################################################################
###############################################################################
################################################################################
################################################################################
###############################################################################
################################################################################
################################################################################
###############################################################################
################################################################################
################################################################################
###############################################################################
################################################################################
################################################################################
###############################################################################
################################################################################
################################################################################
###############################################################################
################################################################################
################################################################################
###############################################################################
################################################################################
################################################################################
###############################################################################
################################################################################
################################################################################
###############################################################################
################################################################################
################################################################################

#Z proporcion


# =====================================================
# PRUEBA Z PARA UNA PROPORCIÓN
# BTX-101 | RESPUESTA AL TRATAMIENTO
# =====================================================

library(ggplot2)

# =====================================================
# 1. SELECCIONAR VARIABLE
# =====================================================

respuesta <- datos$respuesta

# Eliminar valores faltantes
respuesta <- respuesta[!is.na(respuesta)]

# Comprobar que la variable sea binaria
stopifnot(
  all(respuesta %in% c(0, 1))
)

# =====================================================
# 2. ESTADÍSTICA DESCRIPTIVA
# =====================================================

n <- length(respuesta)

exitos <- sum(respuesta == 1)

fracasos <- sum(respuesta == 0)

p_hat <- exitos / n

cat("\n--- ESTADÍSTICA DESCRIPTIVA ---\n")

cat("Número de pacientes =", n, "\n")

cat("Respondieron =", exitos, "\n")

cat("No respondieron =", fracasos, "\n")

cat("Proporción observada =", p_hat, "\n")

cat("Porcentaje de respuesta =",
    round(p_hat * 100, 2), "%\n")


# =====================================================
# 3. VERIFICAR SUPUESTOS DE LA PRUEBA Z
# =====================================================

# Proporción hipotética
p0 <- 0.50

# Éxitos y fracasos esperados bajo H0
exitos_esperados <- n * p0

fracasos_esperados <- n * (1 - p0)

cat("\n--- SUPUESTOS ---\n")

cat("Éxitos esperados =",
    exitos_esperados, "\n")

cat("Fracasos esperados =",
    fracasos_esperados, "\n")

if (exitos_esperados >= 10 &&
    fracasos_esperados >= 10) {

  cat("Aproximación normal: ADECUADA\n")

} else {

  cat("Aproximación normal: NO ADECUADA\n")

}


# =====================================================
# 4. ESTADÍSTICO Z PARA UNA PROPORCIÓN
# =====================================================

# H0: p = 0.50
# H1: p != 0.50

# Error estándar bajo H0
EE <- sqrt(
  p0 * (1 - p0) / n
)

# Estadístico Z
Z <- (p_hat - p0) / EE

# Valor p bilateral
p_valor <- 2 * pnorm(-abs(Z))

cat("\n--- PRUEBA Z PARA UNA PROPORCIÓN ---\n")

cat("Z =", round(Z, 4), "\n")

cat("Valor p =",
    round(p_valor, 6), "\n")


# =====================================================
# 5. INTERVALO DE CONFIANZA DEL 95 %
# MÉTODO DE WILSON
# =====================================================

alpha <- 0.05

z_critico <- qnorm(1 - alpha / 2)

denominador <- 1 + z_critico^2 / n

centro <- (
  p_hat + z_critico^2 / (2 * n)
) / denominador

margen <- (
  z_critico *
  sqrt(
    p_hat * (1 - p_hat) / n +
    z_critico^2 / (4 * n^2)
  )
) / denominador

IC_inf <- centro - margen

IC_sup <- centro + margen

cat("IC 95% Wilson = [",
    round(IC_inf, 4), ",",
    round(IC_sup, 4), "]\n")


# =====================================================
# 6. DECISIÓN ESTADÍSTICA
# =====================================================

if (p_valor < alpha) {

  decision <- "Se rechaza H0"

} else {

  decision <- "No se rechaza H0"

}

cat("Decisión =", decision, "\n")


# =====================================================
# 7. TABLA FINAL
# =====================================================

tabla_Z_prop <- data.frame(

  n = n,

  Respondieron = exitos,

  No_respondieron = fracasos,

  Proporcion_observada = p_hat,

  Porcentaje_respuesta = p_hat * 100,

  Proporcion_hipotetica = p0,

  Error_estandar_H0 = round(EE, 4),

  Z = round(Z, 4),

  Valor_p = signif(p_valor, 6),

  IC95_inferior = round(IC_inf, 4),

  IC95_superior = round(IC_sup, 4),

  Decision = decision

)

print(tabla_Z_prop)


# =====================================================
# 8. GRÁFICA 1: RESPONDIERON VS NO RESPONDIERON
# =====================================================

df_respuesta <- data.frame(

  Respuesta = c(
    "Respondieron",
    "No respondieron"
  ),

  Pacientes = c(
    exitos,
    fracasos
  )
)

grafica_respuesta <- ggplot(
  df_respuesta,
  aes(
    x = Respuesta,
    y = Pacientes,
    fill = Respuesta
  )
) +

  geom_col(
    width = 0.55,
    color = "#1D3557"
  ) +

  geom_text(
    aes(label = Pacientes),
    vjust = -0.5,
    size = 5
  ) +

  scale_fill_manual(
    values = c(
      "Respondieron" = "#2A9D8F",
      "No respondieron" = "#A8DADC"
    )
  ) +

  scale_y_continuous(
    limits = c(0, max(df_respuesta$Pacientes) * 1.2)
  ) +

  labs(
    title = "Respuesta al tratamiento BTX-101",
    subtitle = paste0(
      "Pacientes evaluados: n = ", n
    ),
    x = "",
    y = "Número de pacientes"
  ) +

  theme_classic(base_size = 14) +

  theme(
    legend.position = "none"
  )

print(grafica_respuesta)


# =====================================================
# 9. GRÁFICA 2: PROPORCIÓN OBSERVADA VS REFERENCIA
# =====================================================

df_prop <- data.frame(
  Grupo = "BTX-101",
  Proporcion = p_hat
)

grafica_proporcion <- ggplot(
  df_prop,
  aes(x = Grupo, y = Proporcion)
) +

  geom_col(
    width = 0.35,
    fill = "#2A9D8F",
    color = "#1D3557"
  ) +

  geom_errorbar(
    aes(
      ymin = IC_inf,
      ymax = IC_sup
    ),
    width = 0.12,
    linewidth = 0.9
  ) +

  geom_hline(
    yintercept = p0,
    color = "#E63946",
    linetype = "dashed",
    linewidth = 1
  ) +

  annotate(
    "text",
    x = 1.3,
    y = p0,
    label = "Referencia = 50%",
    color = "#E63946",
    hjust = 0,
    vjust = -0.5,
    size = 4
  ) +

  scale_y_continuous(
    labels = function(x) paste0(round(x * 100), "%"),
    limits = c(0, 1)
  ) +

  labs(
    title = "Prueba Z para una proporción",
    subtitle = paste0(
      "Z = ", round(Z, 3),
      " | p = ", round(p_valor, 4),
      " | IC 95% de Wilson"
    ),
    x = "",
    y = "Proporción de pacientes que respondieron"
  ) +

  theme_classic(base_size = 14)

print(grafica_proporcion)


# =====================================================
# 10. GUARDAR RESULTADOS
# =====================================================

dir.create(
  "Graficas_BTX101",
  showWarnings = FALSE
)

ggsave(
  "Graficas_BTX101/Z_Proporcion_Respuestas.png",
  plot = grafica_respuesta,
  width = 8,
  height = 6,
  dpi = 300,
  bg = "white"
)

ggsave(
  "Graficas_BTX101/Z_Proporcion_IC95.png",
  plot = grafica_proporcion,
  width = 8,
  height = 6,
  dpi = 300,
  bg = "white"
)

write.csv(
  tabla_Z_prop,
  "Graficas_BTX101/Resultados_Z_Proporcion.csv",
  row.names = FALSE
)

cat("\nAnálisis finalizado.\n")






# =====================================================
# VERIFICAR APROXIMACIÓN NORMAL
# PRUEBA Z PARA UNA PROPORCIÓN
# =====================================================

respuesta <- datos$respuesta
respuesta <- respuesta[!is.na(respuesta)]

n <- length(respuesta)
exitos <- sum(respuesta == 1)

p0 <- 0.50

# Éxitos y fracasos esperados
E_exitos <- n * p0
E_fracasos <- n * (1 - p0)

cat("Éxitos esperados =", E_exitos, "\n")
cat("Fracasos esperados =", E_fracasos, "\n")

if (E_exitos >= 10 && E_fracasos >= 10) {
  cat("Aproximación normal adecuada\n")
} else {
  cat("Considerar prueba binomial exacta\n")
}

# Prueba Z sin corrección de continuidad
prueba_z <- prop.test(
  x = exitos,
  n = n,
  p = p0,
  alternative = "two.sided",
  correct = FALSE
)

print(prueba_z)

# Prueba binomial exacta
prueba_exacta <- binom.test(
  x = exitos,
  n = n,
  p = p0,
  alternative = "two.sided"
)

print(prueba_exacta)




# =====================================================
# PRUEBA Z PARA UNA PROPORCIÓN
# BTX-101
# =====================================================

# Seleccionar variable
respuesta <- datos$respuesta
respuesta <- respuesta[!is.na(respuesta)]

# Estadística descriptiva
n <- length(respuesta)
exitos <- sum(respuesta == 1)
p_hat <- exitos / n

# Proporción hipotética
p0 <- 0.50

# Calcular error estándar bajo H0
EE <- sqrt(p0 * (1 - p0) / n)

# Calcular estadístico Z
Z <- (p_hat - p0) / EE

# Calcular valor p bilateral
p_valor <- 2 * pnorm(-abs(Z))

# Mostrar resultados
cat("\n--- PRUEBA Z PARA UNA PROPORCIÓN ---\n")
cat("n =", n, "\n")
cat("Respondieron =", exitos, "\n")
cat("Proporción observada =", p_hat, "\n")
cat("Estadístico Z =", round(Z, 4), "\n")
cat("Valor p =", round(p_valor, 6), "\n")





respuesta <- datos$respuesta

x <- sum(respuesta == 1, na.rm = TRUE)
n <- sum(!is.na(respuesta))

p_hat <- x / n

cat("Proporción observada =", p_hat, "\n")
cat("Porcentaje observado =", p_hat * 100, "%\n")




# =====================================================
# PRUEBA Z PARA UNA PROPORCIÓN
# GRÁFICA SIN INTERVALO DE WILSON
# Tratamiento experimental BTX-101
# =====================================================

library(ggplot2)

# =====================================================
# 1. SELECCIONAR VARIABLE
# =====================================================

respuesta <- datos$respuesta
respuesta <- respuesta[!is.na(respuesta)]

# Verificar que sea binaria
stopifnot(all(respuesta %in% c(0, 1)))

# =====================================================
# 2. ESTADÍSTICA DESCRIPTIVA
# =====================================================

n <- length(respuesta)
exitos <- sum(respuesta == 1)
p_hat <- exitos / n

# Proporción hipotética
p0 <- 0.50

# =====================================================
# 3. PRUEBA Z BILATERAL
# =====================================================

EE <- sqrt(p0 * (1 - p0) / n)

Z <- (p_hat - p0) / EE

p_valor <- 2 * pnorm(-abs(Z))

# =====================================================
# 4. PREPARAR DATOS PARA LA GRÁFICA
# =====================================================

df_prop <- data.frame(
  Grupo = "BTX-101",
  Proporcion = p_hat
)

# =====================================================
# 5. GRÁFICA SIN INTERVALO DE CONFIANZA
# =====================================================

grafica_z_prop <- ggplot(
  df_prop,
  aes(x = Grupo, y = Proporcion)
) +

  # Barra de proporción observada
  geom_col(
    width = 0.35,
    fill = "#2A9D8F",
    color = "#1D3557",
    linewidth = 0.8
  ) +

  # Línea de referencia: 50 %
  geom_hline(
    yintercept = p0,
    color = "#E63946",
    linetype = "dashed",
    linewidth = 1
  ) +

  # Etiqueta del valor observado
  geom_text(
    aes(
      label = paste0(
        round(Proporcion * 100, 1), "%"
      )
    ),
    vjust = -0.8,
    size = 5,
    fontface = "bold"
  ) +

  # Etiqueta de referencia
  annotate(
    "text",
    x = 1.27,
    y = p0,
    label = "Referencia = 50%",
    color = "#E63946",
    hjust = 0,
    vjust = -0.6,
    size = 4
  ) +

  # Eje Y en porcentajes
  scale_y_continuous(
    breaks = seq(0, 1, 0.20),
    labels = function(x) paste0(x * 100, "%"),
    limits = c(0, 1),
    expand = expansion(mult = c(0, 0.03))
  ) +

  # Títulos y estadísticas
  labs(
    title = "Prueba Z para una proporción",
    subtitle = paste0(
      "Z = ", round(Z, 3),
      " | p = ", round(p_valor, 4),
      " | n = ", n
    ),
    x = "",
    y = "Proporción de pacientes que respondieron"
  ) +

  theme_classic(base_size = 14) +

  theme(
    plot.title = element_text(size = 16),
    plot.subtitle = element_text(size = 12),
    axis.text = element_text(color = "black"),
    panel.grid.major.y = element_line(
      color = "grey90",
      linewidth = 0.4
    )
  )

# Mostrar gráfica
print(grafica_z_prop)

# =====================================================
# 6. GUARDAR GRÁFICA A 300 DPI
# =====================================================

dir.create(
  "Graficas_BTX101",
  showWarnings = FALSE
)

ggsave(
  filename = paste0(
    "Graficas_BTX101/",
    "Z_Proporcion_Sin_Wilson_BTX101.png"
  ),
  plot = grafica_z_prop,
  width = 8,
  height = 6,
  dpi = 300,
  bg = "white"
)




# =====================================================
# PRUEBA BINOMIAL BAYESIANA PARA UNA PROPORCIÓN
# TRATAMIENTO EXPERIMENTAL BTX-101
# =====================================================

library(ggplot2)

# =====================================================
# 1. SELECCIONAR VARIABLE
# =====================================================

respuesta <- datos$respuesta
respuesta <- respuesta[!is.na(respuesta)]

stopifnot(all(respuesta %in% c(0, 1)))

# Estadística descriptiva
n <- length(respuesta)
x <- sum(respuesta == 1)

p_hat <- x / n

cat("\n--- DATOS DE LA MUESTRA ---\n")
cat("Pacientes =", n, "\n")
cat("Respondieron =", x, "\n")
cat("No respondieron =", n - x, "\n")
cat("Proporción observada =", p_hat, "\n")


# =====================================================
# 2. DEFINIR HIPÓTESIS Y DISTRIBUCIÓN PREVIA
# =====================================================

# H0: p = 0.50
# H1: p != 0.50

p0 <- 0.50

# Distribución previa uniforme
# p ~ Beta(a, b)

a <- 1
b <- 1


# =====================================================
# 3. DISTRIBUCIÓN POSTERIOR
# =====================================================

# p | datos ~ Beta(a + x, b + n - x)

a_post <- a + x
b_post <- b + n - x

# Media posterior
media_post <- a_post / (a_post + b_post)

# Mediana posterior
mediana_post <- qbeta(
  0.50,
  a_post,
  b_post
)

# Moda posterior
moda_post <- (a_post - 1) /
  (a_post + b_post - 2)

# Intervalo creíble del 95 %
IC_post <- qbeta(
  c(0.025, 0.975),
  a_post,
  b_post
)

cat("\n--- DISTRIBUCIÓN POSTERIOR ---\n")

cat("Beta posterior = Beta(",
    a_post, ",", b_post, ")\n")

cat("Media posterior =",
    round(media_post, 4), "\n")

cat("Mediana posterior =",
    round(mediana_post, 4), "\n")

cat("Moda posterior =",
    round(moda_post, 4), "\n")

cat("IC creíble 95% = [",
    round(IC_post[1], 4), ",",
    round(IC_post[2], 4), "]\n")


# =====================================================
# 4. PROBABILIDAD POSTERIOR
# =====================================================

# Probabilidad de que p sea mayor a 0.50
# Condicionada al modelo H1

prob_superior <- 1 - pbeta(
  p0,
  a_post,
  b_post
)

cat("\nP(p > 0.50 | datos, H1) =",
    round(prob_superior, 6), "\n")


# =====================================================
# 5. FACTOR DE BAYES
# =====================================================

# H0: p = p0
# H1: p ~ Beta(a,b)
#
# BF10 = P(datos | H1) / P(datos | H0)
#
# Se calcula en escala logarítmica
# para evitar problemas numéricos

log_BF10 <- lbeta(
  a_post,
  b_post
) -
  lbeta(a, b) -
  x * log(p0) -
  (n - x) * log1p(-p0)

BF10 <- exp(log_BF10)
BF01 <- 1 / BF10

cat("\n--- FACTOR DE BAYES ---\n")

cat("BF10 =", round(BF10, 4), "\n")
cat("BF01 =", round(BF01, 4), "\n")


# =====================================================
# 6. PREPARAR LA DISTRIBUCIÓN POSTERIOR
# =====================================================

p_valores <- seq(
  0.001,
  0.999,
  length.out = 2000
)

df_posterior <- data.frame(
  Proporcion = p_valores,
  Densidad = dbeta(
    p_valores,
    a_post,
    b_post
  )
)


# =====================================================
# 7. GRÁFICA BAYESIANA
# =====================================================

grafica_bayes_prop <- ggplot(
  df_posterior,
  aes(x = Proporcion, y = Densidad)
) +

  geom_area(
    fill = "#A8DADC",
    alpha = 0.75
  ) +

  geom_line(
    color = "#1D3557",
    linewidth = 1.2
  ) +

  # Línea roja: H0
  geom_vline(
    xintercept = p0,
    color = "#E63946",
    linetype = "dashed",
    linewidth = 1
  ) +

  # Línea verde: media posterior
  geom_vline(
    xintercept = media_post,
    color = "#2A9D8F",
    linewidth = 1
  ) +

  # Eje X en porcentajes
  scale_x_continuous(
    labels = function(x) paste0(
      round(x * 100), "%"
    ),
    breaks = seq(0, 1, 0.1),
    limits = c(0, 1)
  ) +

  labs(
    title = "Análisis bayesiano: proporción de respuesta",
    subtitle = paste0(
      "BF10 = ", round(BF10, 3),
      " | Media posterior = ",
      round(media_post * 100, 2), "%",
      " | n = ", n
    ),
    x = "Proporción de pacientes que respondieron",
    y = "Densidad posterior",
    caption = paste0(
      "Línea roja: H0 = 50%",
      " | Línea verde: media posterior"
    )
  ) +

  theme_classic(base_size = 14) +

  theme(
    plot.title = element_text(size = 16),
    plot.subtitle = element_text(size = 12),
    plot.caption = element_text(size = 10)
  )

print(grafica_bayes_prop)


# =====================================================
# 8. GUARDAR GRÁFICA Y RESULTADOS
# =====================================================

dir.create(
  "Graficas_BTX101",
  showWarnings = FALSE
)

ggsave(
  "Graficas_BTX101/Bayesiana_Proporcion_BTX101.png",
  plot = grafica_bayes_prop,
  width = 9,
  height = 6,
  dpi = 300,
  bg = "white"
)

tabla_bayes_prop <- data.frame(
  n = n,
  Exitos = x,
  Proporcion_observada = p_hat,
  Media_posterior = media_post,
  IC95_inferior = IC_post[1],
  IC95_superior = IC_post[2],
  Prob_p_mayor_050 = prob_superior,
  BF10 = BF10,
  BF01 = BF01
)

print(tabla_bayes_prop)

write.csv(
  tabla_bayes_prop,
  "Graficas_BTX101/Resultados_Bayes_Proporcion.csv",
  row.names = FALSE
)





# =====================================================
# ANÁLISIS BINOMIAL BAYESIANO CON MCMC
# BTX-101 | PROPORCIÓN DE RESPUESTA
# 4 CADENAS + CONVERGENCIA
# =====================================================

library(ggplot2)
library(posterior)

set.seed(2026)

# =====================================================
# 1. SELECCIONAR VARIABLE
# =====================================================

respuesta <- datos$respuesta
respuesta <- respuesta[!is.na(respuesta)]

stopifnot(
  length(respuesta) > 0,
  all(respuesta %in% c(0, 1))
)

n <- length(respuesta)
x <- sum(respuesta == 1)

p_hat <- x / n
p0 <- 0.50

cat("\n--- DATOS OBSERVADOS ---\n")
cat("Pacientes =", n, "\n")
cat("Respondieron =", x, "\n")
cat("No respondieron =", n - x, "\n")
cat("Proporción observada =", p_hat, "\n")


# =====================================================
# 2. DISTRIBUCIÓN PREVIA Y POSTERIOR
# =====================================================

# H0: p = 0.50
# H1: p ~ Beta(1,1)

a <- 1
b <- 1

a_post <- a + x
b_post <- b + n - x

media_exacta <- a_post / (a_post + b_post)

IC_exacto <- qbeta(
  c(0.025, 0.975),
  a_post,
  b_post
)

# Factor de Bayes analítico
log_BF10 <- lbeta(a_post, b_post) -
  lbeta(a, b) -
  x * log(p0) -
  (n - x) * log1p(-p0)

BF10 <- exp(log_BF10)
BF01 <- 1 / BF10

cat("\n--- RESULTADOS EXACTOS ---\n")
cat("Posterior: Beta(", a_post, ",",
    b_post, ")\n")
cat("Media posterior =", round(media_exacta, 5), "\n")
cat("IC creíble 95% =", round(IC_exacto, 5), "\n")
cat("BF10 =", round(BF10, 4), "\n")
cat("BF01 =", round(BF01, 4), "\n")


# =====================================================
# 3. CONFIGURACIÓN MCMC
# =====================================================

n_cadenas <- 4
n_iter <- 6000
warmup <- 2000

# Puntos iniciales distintos
inicios <- c(0.15, 0.35, 0.65, 0.85)

# Tamaño de propuesta en escala logit
propuesta_sd <- 0.55


# =====================================================
# 4. FUNCIÓN LOG-POSTERIOR
# =====================================================

# Trabajaremos en escala logit para garantizar
# que todas las proporciones estén entre 0 y 1.

log_posterior_theta <- function(theta) {

  p <- plogis(theta)

  # Beta posterior + Jacobiano
  dbeta(
    p,
    shape1 = a_post,
    shape2 = b_post,
    log = TRUE
  ) +
    log(p) +
    log1p(-p)
}


# =====================================================
# 5. ALGORITMO METROPOLIS-HASTINGS
# =====================================================

metropolis_prop <- function(
    inicio,
    iteraciones,
    propuesta_sd
) {

  theta <- numeric(iteraciones)

  theta[1] <- qlogis(inicio)

  aceptadas <- 0

  for (i in 2:iteraciones) {

    actual <- theta[i - 1]

    propuesta <- rnorm(
      1,
      mean = actual,
      sd = propuesta_sd
    )

    log_alpha <- log_posterior_theta(propuesta) -
      log_posterior_theta(actual)

    if (log(runif(1)) < log_alpha) {

      theta[i] <- propuesta
      aceptadas <- aceptadas + 1

    } else {

      theta[i] <- actual
    }
  }

  list(
    p = plogis(theta),
    aceptacion = aceptadas / (iteraciones - 1)
  )
}


# =====================================================
# 6. EJECUTAR LAS CUATRO CADENAS
# =====================================================

resultados <- lapply(
  seq_len(n_cadenas),
  function(i) {

    metropolis_prop(
      inicio = inicios[i],
      iteraciones = n_iter,
      propuesta_sd = propuesta_sd
    )
  }
)


# =====================================================
# 7. ELIMINAR CALENTAMIENTO
# =====================================================

indices_post <- (warmup + 1):n_iter

cadenas_post <- sapply(
  resultados,
  function(r) r$p[indices_post]
)

colnames(cadenas_post) <- paste0(
  "Cadena_", 1:n_cadenas
)

cat("\nMuestras posteriores conservadas =",
    length(cadenas_post), "\n")


# =====================================================
# 8. CONVERTIR AL FORMATO DE posterior
# =====================================================

array_post <- array(
  cadenas_post,
  dim = c(
    nrow(cadenas_post),
    n_cadenas,
    1
  ),
  dimnames = list(
    NULL,
    paste0("Cadena_", 1:n_cadenas),
    "p"
  )
)

draws <- posterior::as_draws_array(array_post)


# =====================================================
# 9. DIAGNÓSTICOS DE CONVERGENCIA
# =====================================================

diagnosticos <- posterior::summarise_draws(
  draws,
  mean,
  sd,
  rhat,
  ess_bulk,
  ess_tail
)

cat("\n--- DIAGNÓSTICOS MCMC ---\n")

print(diagnosticos)

# Tasa de aceptación de cada cadena
tasas_aceptacion <- sapply(
  resultados,
  function(r) r$aceptacion
)

cat("\n--- TASAS DE ACEPTACIÓN ---\n")

print(round(tasas_aceptacion, 3))


# =====================================================
# 10. RESUMEN POSTERIOR MCMC
# =====================================================

muestras <- as.vector(cadenas_post)

media_MCMC <- mean(muestras)

IC_MCMC <- quantile(
  muestras,
  c(0.025, 0.975)
)

prob_superior <- mean(muestras > p0)

cat("\n--- RESULTADOS POSTERIORES MCMC ---\n")

cat("Media posterior MCMC =",
    round(media_MCMC, 5), "\n")

cat("IC creíble MCMC 95% =",
    round(IC_MCMC, 5), "\n")

cat("P(p > 0.50 | datos, H1) =",
    round(prob_superior, 5), "\n")

cat("\n--- COMPARACIÓN CON SOLUCIÓN EXACTA ---\n")

cat("Media exacta =", round(media_exacta, 5), "\n")
cat("Media MCMC =", round(media_MCMC, 5), "\n")


# =====================================================
# 11. PREPARAR DATOS PARA GRÁFICAS
# =====================================================

df_mcmc <- data.frame(
  Iteracion = rep(
    seq_len(nrow(cadenas_post)),
    times = n_cadenas
  ),

  Proporcion = as.vector(cadenas_post),

  Cadena = factor(
    rep(
      1:n_cadenas,
      each = nrow(cadenas_post)
    )
  )
)


# =====================================================
# 12. GRÁFICA 1: TRACE PLOTS
# =====================================================

grafica_trazas_prop <- ggplot(
  df_mcmc,
  aes(x = Iteracion, y = Proporcion)
) +

  geom_line(
    color = "#457B9D",
    linewidth = 0.25,
    alpha = 0.8
  ) +

  facet_wrap(
    ~ Cadena,
    ncol = 2
  ) +

  geom_hline(
    yintercept = media_exacta,
    color = "#E63946",
    linewidth = 0.8
  ) +

  scale_y_continuous(
    labels = function(z) paste0(
      round(z * 100), "%"
    )
  ) +

  labs(
    title = "Convergencia MCMC: respuesta a BTX-101",
    subtitle = "Cuatro cadenas de Metropolis-Hastings",
    x = "Iteración posterior",
    y = "Proporción poblacional de respuesta"
  ) +

  theme_classic(base_size = 13)

print(grafica_trazas_prop)


# =====================================================
# 13. GRÁFICA 2: DENSIDADES POR CADENA
# =====================================================

grafica_densidades_prop <- ggplot(
  df_mcmc,
  aes(
    x = Proporcion,
    color = Cadena
  )
) +

  geom_density(
    linewidth = 1
  ) +

  # H0: p = 0.50
  geom_vline(
    xintercept = p0,
    color = "#E63946",
    linetype = "dashed",
    linewidth = 1
  ) +

  scale_x_continuous(
    labels = function(z) paste0(
      round(z * 100), "%"
    )
  ) +

  labs(
    title = "Distribuciones posteriores por cadena",
    subtitle = paste0(
      "BF10 = ", round(BF10, 3),
      " | Media posterior = ",
      round(media_exacta * 100, 2),
      "%"
    ),
    x = "Proporción poblacional de respuesta",
    y = "Densidad posterior",
    color = "Cadena"
  ) +

  theme_classic(base_size = 13)

print(grafica_densidades_prop)


# =====================================================
# 14. GRÁFICA 3: POSTERIOR EXACTA VS MCMC
# =====================================================

df_exacta <- data.frame(
  Proporcion = seq(0.3, 0.85, length.out = 1000)
)

df_exacta$Densidad <- dbeta(
  df_exacta$Proporcion,
  a_post,
  b_post
)

grafica_comparacion <- ggplot() +

  geom_histogram(
    data = data.frame(Proporcion = muestras),
    aes(
      x = Proporcion,
      y = after_stat(density)
    ),
    bins = 40,
    fill = "#A8DADC",
    color = "white",
    alpha = 0.7
  ) +

  geom_line(
    data = df_exacta,
    aes(
      x = Proporcion,
      y = Densidad
    ),
    color = "#1D3557",
    linewidth = 1.2
  ) +

  geom_vline(
    xintercept = p0,
    color = "#E63946",
    linetype = "dashed",
    linewidth = 1
  ) +

  scale_x_continuous(
    labels = function(z) paste0(
      round(z * 100), "%"
    )
  ) +

  labs(
    title = "Posterior exacta vs simulación MCMC",
    subtitle = "Modelo Beta-Binomial | BTX-101",
    x = "Proporción poblacional de respuesta",
    y = "Densidad posterior",
    caption = paste0(
      "Curva azul: posterior exacta Beta(",
      a_post, ", ", b_post, ")"
    )
  ) +

  theme_classic(base_size = 13)

print(grafica_comparacion)


# =====================================================
# 15. GUARDAR GRÁFICAS Y DIAGNÓSTICOS
# =====================================================

dir.create(
  "Graficas_BTX101",
  showWarnings = FALSE
)

ggsave(
  "Graficas_BTX101/MCMC_Proporcion_Traceplot.png",
  plot = grafica_trazas_prop,
  width = 10,
  height = 7,
  dpi = 300,
  bg = "white"
)

ggsave(
  "Graficas_BTX101/MCMC_Proporcion_Densidades.png",
  plot = grafica_densidades_prop,
  width = 9,
  height = 6,
  dpi = 300,
  bg = "white"
)

ggsave(
  "Graficas_BTX101/MCMC_Proporcion_Comparacion.png",
  plot = grafica_comparacion,
  width = 9,
  height = 6,
  dpi = 300,
  bg = "white"
)

write.csv(
  as.data.frame(diagnosticos),
  "Graficas_BTX101/MCMC_Proporcion_Diagnosticos.csv",
  row.names = FALSE
)

cat("\nAnálisis MCMC finalizado.\n")







