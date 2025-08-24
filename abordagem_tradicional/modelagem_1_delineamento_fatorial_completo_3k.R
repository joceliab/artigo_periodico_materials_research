# =========================================
# Pacotes
# =========================================
req <- c("rsm","ggplot2","dplyr","patchwork","viridis","plotly","broom")
for(p in req){
  if(!suppressWarnings(require(p, character.only=TRUE))){
    install.packages(p, dependencies=TRUE)
    library(p, character.only=TRUE)
  }
}

# =========================================
# 1) Dados e DOE fatorial 3^4
# =========================================
DOE <- expand.grid(
  C = c(8, 10, 12),
  F = c(0.2, 0.5, 0.8),
  D = c(15, 17.5, 20),
  V = c(15, 20, 25)
)
DOE <- DOE[order(DOE$C, DOE$F, DOE$D, DOE$V), ]

DOE$diametro <- c(
  3.067, 2.980, 3.755, 2.664, 2.622, 3.039, 3.217, 2.971, 1.712,
  2.989, 3.293, 3.015, 3.363, 2.037, 2.899, 3.552, 2.829, 3.135,
  2.530, 2.618, 3.002, 2.799, 2.626, 3.051, 2.887, 3.213, 2.561,
  2.807, 1.885, 2.939, 3.667, 2.698, 2.804, 3.190, 3.120, 3.107,
  2.928, 2.661, 2.993, 2.661, 2.108, 2.565, 3.121, 2.699, 2.486,
  2.353, 1.957, 3.512, 2.425, 3.196, 2.933, 2.612, 1.871, 3.020,
  3.333, 3.165, 2.240, 4.026, 2.154, 2.787, 3.855, 3.061, 2.925,
  4.342, 3.101, 3.333, 3.300, 2.660, 4.007, 3.275, 3.108, 2.754,
  2.884, 2.406, 3.098, 3.438, 3.078, 4.119, 3.440, 3.185, 3.405
)

# =========================================
# 2) Codificação centrada e escalada
# =========================================
DOE_cod <- coded.data(
  DOE,
  C_cod ~ (C - 10)/2,
  F_cod ~ (F - 0.5)/0.3,
  D_cod ~ (D - 17.5)/2.5,
  V_cod ~ (V - 20)/5
)

# =========================================
# 3) Modelo completo (2ª ordem) em termos codificados
# =========================================
modelo_full <- rsm(diametro ~ SO(C_cod, F_cod, D_cod, V_cod), data = DOE_cod)

# =========================================
# 4) Seleção hierárquica automática de termos significativos
# =========================================
tab <- as.data.frame(coef(summary(modelo_full)))
tab$term <- rownames(tab)

tab_sig <- subset(tab, term != "(Intercept)" & `Pr(>|t|)` < 0.05)

lin_vars  <- character(0)
int_terms <- list()
quad_vars <- character(0)

for(nm in tab_sig$term){
  if(grepl(":", nm)){
    pr <- strsplit(nm, ":", fixed=TRUE)[[1]]
    int_terms <- append(int_terms, list(pr))
    lin_vars <- unique(c(lin_vars, pr))
  } else if(grepl("\\^2", nm)){
    var <- sub("^I\\((.*)\\^2\\)$", "\\1", nm)
    var <- sub("\\^2$", "", var)
    quad_vars <- unique(c(quad_vars, var))
    lin_vars  <- unique(c(lin_vars, var))
  } else {
    lin_vars <- unique(c(lin_vars, nm))
  }
}

lin_part  <- if(length(lin_vars)>0) paste(lin_vars, collapse=" + ") else "1"
int_part  <- if(length(int_terms)>0) paste(sapply(int_terms, function(p) paste(p, collapse=":")), collapse=" + ") else ""
quad_part <- if(length(quad_vars)>0) paste(paste0("I(", quad_vars, "^2)"), collapse=" + ") else ""

rhs <- paste(c(lin_part, int_part, quad_part), collapse=" + ")
rhs <- gsub("\\+\\s*$","",rhs)
rhs <- gsub("^\\s*\\+\\s*","",rhs)

f_red <- as.formula(paste("diametro ~", rhs))

# Ajuste do modelo reduzido com os termos selecionados
modelo_red <- lm(f_red, data = DOE_cod)

# =========================================
# 5) Métricas de ajuste para ambos modelos
# =========================================
metricas <- function(obj, yobs){
  pr <- as.numeric(predict(obj))
  res <- yobs - pr
  mse  <- mean(res^2)
  rmse <- sqrt(mse)
  mape <- mean(abs(res / yobs)) * 100
  data.frame(MSE=mse, RMSE=rmse, MAPE=mape, R2aj = summary(obj)$adj.r.squared,
             AIC=AIC(obj), BIC=BIC(obj), row.names = NULL)
}

met_full <- metricas(modelo_full, DOE_cod$diametro)
met_red  <- metricas(modelo_red,  DOE_cod$diametro)

comp <- rbind(
  cbind(Modelo="Completo_RSM", met_full),
  cbind(Modelo="Reduzido_LM",  met_red)
)
print(comp)

# =========================================
# 6) Escolha do melhor modelo (critério principal: BIC)
#     Empate no BIC -> menor AIC -> maior R²aj.
# =========================================
idx_best <- with(comp, order(BIC, AIC, -R2aj))[1]
modelo_escolhido <- comp$Modelo[idx_best]
cat("\n>>> Melhor modelo (critério: BIC):", modelo_escolhido, "\n")

# Converter o melhor para objeto 'rsm' (necessário p/ análise canônica e RSM)
if(modelo_escolhido == "Completo_RSM"){
  best_rsm <- modelo_full
} else {
  # reconstituir como rsm com a mesma fórmula reduzida
  best_rsm <- rsm(formula(f_red), data = DOE_cod)
}

# =========================================
# 7) Diagnósticos de resíduos (apenas melhor modelo)
# =========================================
pred_best  <- as.numeric(predict(best_rsm))
resid_best <- as.numeric(residuals(best_rsm))
aux <- data.frame(pred = pred_best, resid = resid_best, obs = DOE_cod$diametro)

g_res1 <- ggplot(aux, aes(x = pred, y = resid)) +
  geom_point() + geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title="Resíduos vs Ajustado", x="Ajustado", y="Resíduo") + theme_minimal()

g_res2 <- ggplot(aux, aes(sample = resid)) +
  stat_qq() + stat_qq_line() + labs(title="QQ-Plot dos Resíduos") + theme_minimal()

g_res3 <- ggplot(aux, aes(x = resid)) +
  geom_histogram(bins = 15, color="black") +
  labs(title="Histograma dos Resíduos", x="Resíduo", y="Frequência") + theme_minimal()

g_res4 <- ggplot(aux, aes(x = pred, y = obs)) +
  geom_point() + geom_abline(slope=1, intercept=0, linetype="dashed") +
  labs(title="Observado vs Ajustado", x="Ajustado", y="Observado") + theme_minimal()

print((g_res1 + g_res2) / (g_res3 + g_res4))

# =========================================
# 8) Predições em novos pontos (escala real), com IP
# =========================================
novos <- data.frame(
  C = c(9, 10, 11),
  F = c(0.3, 0.5, 0.7),
  D = rep(17.5, 3),
  V = c(18, 20, 22)
)
novos_cod <- within(novos, {
  C_cod <- (C - 10)/2
  F_cod <- (F - 0.5)/0.3
  D_cod <- (D - 17.5)/2.5
  V_cod <- (V - 20)/5
})

pi_best <- as.data.frame(predict(best_rsm, newdata = novos_cod, interval = "prediction"))
novos_out <- cbind(novos, pi_best)
cat("\n=== Predições (melhor modelo) ===\n"); print(novos_out)

# =========================================
# 9) Mapas de contorno 2D (apenas melhor modelo)
# =========================================
contorno_2D <- function(mod_rsm, xname, yname, fix_vals, nx=120, ny=120, titulo=""){
  ranges <- list(C=c(8,12), F=c(0.2,0.8), D=c(15,20), V=c(15,25))
  seqx <- seq(ranges[[xname]][1], ranges[[xname]][2], length.out=nx)
  seqy <- seq(ranges[[yname]][1], ranges[[yname]][2], length.out=ny)
  grid <- expand.grid(x = seqx, y = seqy)
  names(grid) <- c(xname, yname)
  for (nm in c("C","F","D","V")){
    if (!nm %in% c(xname,yname)) grid[[nm]] <- fix_vals[[nm]]
  }
  grid$C_cod <- (grid$C - 10)/2
  grid$F_cod <- (grid$F - 0.5)/0.3
  grid$D_cod <- (grid$D - 17.5)/2.5
  grid$V_cod <- (grid$V - 20)/5

  grid$pred <- as.numeric(predict(mod_rsm, newdata = grid))
  p_cont <- ggplot(grid, aes_string(x = xname, y = yname)) +
    geom_tile(aes(fill = pred)) +
    geom_contour(aes(z = pred), color="black", alpha=0.5) +
    scale_fill_viridis_c(option = "D") +
    labs(title = titulo, fill = "Diâmetro (µm)") + theme_minimal()

  p_disc <- ggplot(grid, aes_string(x = xname, y = yname, z = "pred")) +
    geom_contour_filled(bins = 20) +
    scale_fill_viridis_d(option = "D") +
    labs(title = paste0(titulo, " (faixas)"), fill = "Diâmetro (µm)") + theme_minimal()

  p_cont + p_disc
}

g_CV <- contorno_2D(best_rsm, "C","V", fix_vals = list(F=0.5, D=17.5), titulo="C × V")
g_FV <- contorno_2D(best_rsm, "F","V", fix_vals = list(C=10, D=17.5),  titulo="F × V")
g_DV <- contorno_2D(best_rsm, "D","V", fix_vals = list(C=10, F=0.5),  titulo="D × V")
g_CF <- contorno_2D(best_rsm, "C","F", fix_vals = list(D=17.5, V=20), titulo="C × F")
print((g_CV / g_FV) | (g_DV / g_CF))

# =========================================
# 10) Superfícies 3D (apenas melhor modelo)
# =========================================
superficie_3D <- function(mod_rsm, xname, yname, fix_vals, nx=60, ny=60, titulo=""){
  ranges <- list(C=c(8,12), F=c(0.2,0.8), D=c(15,20), V=c(15,25))
  seqx <- seq(ranges[[xname]][1], ranges[[xname]][2], length.out=nx)
  seqy <- seq(ranges[[yname]][1], ranges[[yname]][2], length.out=ny)
  grid <- expand.grid(x = seqx, y = seqy)
  names(grid) <- c(xname, yname)
  for (nm in c("C","F","D","V")){
    if (!nm %in% c(xname,yname)) grid[[nm]] <- fix_vals[[nm]]
  }
  grid$C_cod <- (grid$C - 10)/2
  grid$F_cod <- (grid$F - 0.5)/0.3
  grid$D_cod <- (grid$D - 17.5)/2.5
  grid$V_cod <- (grid$V - 20)/5

  grid$pred <- as.numeric(predict(mod_rsm, newdata = grid))
  # expand.grid gera blocos de y para cada x; organize z como [nx x ny]
  z <- matrix(grid$pred, nrow = nx, ncol = ny, byrow = TRUE)

  plot_ly(
    x = seqx, y = seqy, z = z
  ) %>%
    add_surface(contours = list(z = list(show=TRUE, usecolormap=TRUE, project=list(z=TRUE)))) %>%
    layout(
      title = titulo,
      scene = list(
        xaxis = list(title = xname),
        yaxis = list(title = yname),
        zaxis = list(title = "Diâmetro (µm)")
      )
    )
}

# exemplos 3D (apenas dois pares ilustrativos)
superficie_3D(best_rsm, "C","V", list(F=0.5, D=17.5), titulo="3D: C × V (melhor modelo)")
superficie_3D(best_rsm, "F","V", list(C=10, D=17.5),  titulo="3D: F × V (melhor modelo)")

# =========================================
# 11) Análise canônica do melhor modelo (objeto rsm)
#     Diagnostica a natureza da superfície de resposta (máximo, mínimo, sela)
# =========================================
cat("\n=== Análise canônica (melhor modelo) ===\n")
print(summary(best_rsm)$canonical)

# =========================================
# 12) Saída consolidada de métricas e equações
# =========================================
cat("\n=== Métricas de comparação ===\n"); print(comp)

cat("\n=== Equações em termos codificados ===\n")
eq_print <- function(obj){
  cf <- coef(obj)
  termos <- names(cf)
  out <- sprintf("ŷ = %.6f", cf[1])
  if(length(cf)>1){
    for(i in 2:length(cf)){
      s <- ifelse(cf[i]>=0, " + ", " - ")
      out <- paste0(out, s, sprintf("%.6f", abs(cf[i])), "·", termos[i])
    }
  }
  out
}
cat("\nModelo completo:\n", eq_print(modelo_full), "\n", sep="")
cat("\nModelo reduzido:\n", eq_print(modelo_red),  "\n", sep="")
cat("\n>>> Melhor modelo:", modelo_escolhido, "\n")


