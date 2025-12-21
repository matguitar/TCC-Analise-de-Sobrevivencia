options(scipen = 999)  # evita notação científica
set.seed(1)

library(survival)
library(nlme)   # para fdHess
library(MASS)   # para ginv

##-----------------------------------------------------------------
# Dados Infarto Agudo do Miocárdio
##-----------------------------------------------------------------
dados <- read.table("dados-infarto.txt")

# Tempo e status
y      <- dados$fim - dados$ini
status <- dados$status

# Covariáveis: idade (<65 / ≥65) e sexo
x1 	   <- factor(ifelse(dados$idade < 65, "(<65)", "(>=65)"))
x2 	   <- factor(dados$sexo)
dados      <- data.frame(dados, x1, x2)

#-----------------------------------------------------------------
# Matrizes de covariáveis
#-----------------------------------------------------------------
Xe_sigma <- model.matrix(~ 1, data = dados)         # para σ
Xe_mu    <- model.matrix(~ 1 + x1 + x2, data = dados)  # para μ
Xe_p0    <- model.matrix(~ 1 + x1 + x2, data = dados)  # para p₀
Xe_phi   <- model.matrix(~ 1, data = dados)         # para φ

#-----------------------------------------------------------------
# Função de verossimilhança: Lognormal Discreta + fração de cura
#-----------------------------------------------------------------
fvero <- function(theta) {
  n_sigma <- ncol(Xe_sigma)
  n_mu    <- ncol(Xe_mu)
  n_p0    <- ncol(Xe_p0)
  n_phi   <- ncol(Xe_phi)

  # índices dos parâmetros
  beta_sigma <- theta[1:n_sigma]
  beta_mu    <- theta[(n_sigma+1):(n_sigma+n_mu)]
  beta_p0    <- theta[(n_sigma+n_mu+1):(n_sigma+n_mu+n_p0)]
  beta_phi   <- theta[(n_sigma+n_mu+n_p0+1):(n_sigma+n_mu+n_p0+n_phi)]

  # transformações
  sigma <- exp(Xe_sigma %*% beta_sigma)
  mu    <- Xe_mu %*% beta_mu
  p0    <- 1/(1 + exp(-Xe_p0 %*% beta_p0))
  phi   <- exp(Xe_phi %*% beta_phi)
  eta   <- phi^(-1) * (p0^(-phi) - 1)

  # Função de sobrevivência Lognormal Discreta
  S_LND <- function(y, sigma, mu) {
    1 - pnorm((log(y) - mu) / sigma)
  }

  vF1   <- 1 - S_LND(y, sigma, mu)
  Spop1 <- (1 + phi * eta * vF1)^(-1/phi)

  vF2   <- 1 - S_LND(y + 1, sigma, mu)
  Spop2 <- (1 + phi * eta * vF2)^(-1/phi)

  fpop  <- Spop1 - Spop2

  # log-verossimilhança
  loglik <- sum(status * log(fpop) + (1 - status) * log(Spop1))
  return(loglik)
}

#-----------------------------------------------------------------
# Estimação
#-----------------------------------------------------------------
n_sigma <- ncol(Xe_sigma)
n_mu    <- ncol(Xe_mu)
n_p0    <- ncol(Xe_p0)
n_phi   <- ncol(Xe_phi)

theta0 <- rep(0, n_sigma + n_mu + n_p0 + n_phi)

mgg1 <- optim(theta0, fvero, method = "SANN",
              control = list(fnscale = -1, maxit = 10000))
theta0 <- mgg1$par

mgg <- optim(theta0, fvero, method = "BFGS",
             control = list(fnscale = -1, maxit = 1000))

#-----------------------------------------------------------------
# Parâmetros estimados
#-----------------------------------------------------------------
beta_sigma <- mgg$par[1:n_sigma]
beta_mu    <- mgg$par[(n_sigma+1):(n_sigma+n_mu)]
beta_p0    <- mgg$par[(n_sigma+n_mu+1):(n_sigma+n_mu+n_p0)]
beta_phi   <- mgg$par[(n_sigma+n_mu+n_p0+1):(n_sigma+n_mu+n_p0+n_phi)]

sigma <- exp(Xe_sigma %*% beta_sigma)
mu    <- Xe_mu %*% beta_mu
p0    <- 1/(1 + exp(-Xe_p0 %*% beta_p0))
phi   <- exp(Xe_phi %*% beta_phi)
eta   <- phi^(-1) * (p0^(-phi) - 1)

S_LND <- function(y, sigma, mu) {
  1 - pnorm((log(y) - mu) / sigma)
}

vF   <- 1 - S_LND(y + 1, sigma, mu)
Spop <- (1 + phi * eta * vF)^(-1/phi)
p0   <- (1 + phi * eta)^(-1/phi)

#-----------------------------------------------------------------
# Plot função de sobrevivência empírica vs modelo
#-----------------------------------------------------------------
#mKM <- survfit(Surv(y, status) ~ x1 + x2, se.fit = FALSE)
#plot(mKM, xlab = "Tempo (meses)", ylab = "Função de sobrevivência",
#     cex.axis = 1.5, cex.lab = 1.5)
#points(y, Spop, pch = 20, cex = 0.8)

#-----------------------------------------------------------------
# Erros padrão (pseudo-inversa da Hessiana)
#-----------------------------------------------------------------

Estimativa <- mgg$par
obsinf <- -fdHess(Estimativa, fvero)$Hessian
covmat <- ginv(obsinf)
setheta <- sqrt(pmax(0, diag(covmat)))
tvals <- Estimativa / setheta
pvals <- 2 * (1 - pnorm(abs(tvals)))

# nomes dos parâmetros
nomes <- c(paste0("b_sigma_", colnames(Xe_sigma)),
           paste0("b_mu_",    colnames(Xe_mu)),
           paste0("b_p0_",    colnames(Xe_p0)),
           paste0("b_phi_",   colnames(Xe_phi)))

mest <- cbind(
  Estimate = Estimativa,
  `s.e.` = setheta,
  `|t|` = abs(tvals),
  `p-value` = pvals
)
rownames(mest) <- nomes

print(mest, 6)

cormat <- cov2cor(covmat)
print(cormat, 3)

#-----------------------------------------------------------------
# AIC
#-----------------------------------------------------------------
2 * length(Estimativa) - 2 * fvero(Estimativa)

# Teste Log-Rank
teste_log_rank <- survdiff(Surv(y, status) ~ x1 + x2)
print(teste_log_rank)


library(hnp)
library(ggplot2)
#-----------------------------------------------------------------
# CÁLCULO DOS RESÍDUOS QUANTÍLICOS
#-----------------------------------------------------------------

#vF   <- 1 - S_LND(y + 1, sigma, mu)
#Spop <- (1 + phi * eta * vF)^(-1/phi)

# FUNÇÕES DO MODELO (Reutilizadas do seu código)
S_LND <- function(y, sigma, mu) {
  1 - pnorm((log(y) - mu) / sigma)
}


# Função de Distribuição Acumulada Populacional (F_pop = 1 - S_pop)
F_pop <- function(y, sigma, mu, phi, eta) {
  vF <- 1 - S_LND(y+1, sigma, mu) 
  S_pop_t <- (1 + phi * eta * vF)^(-1/phi)
  return(1-S_pop_t)
}

# Valores das Funções de Distribuição Acumulada
# F_pop no tempo de observação (y_i)
F_pop_y <- F_pop(y, sigma, mu, phi, eta)

# F_pop no tempo imediatamente anterior (y_i - 1). Substituímos 0 por um valor pequeno
# para evitar problemas com log(0) no cálculo de S_LND quando y=1.
#y_minus_1 <- pmax(y - 1, 1e-10)
y_minus_1 <- y - 1
F_pop_y_minus_1 <- F_pop(y_minus_1, sigma, mu, phi, eta)

epsilon <- 1e-6
zero=0

# Geração de números aleatórios U(0, 1) para a aleatorização
#set.seed(42) # Usando uma nova seed para replicabilidade dos resíduos quantílicos
rand_u <- runif(length(y),zero,1)
rand_u_prime <- runif(length(y),zero,1)

u_hat <- numeric(length(y))

for (i in 1:length(y)) {
  # 1. Indivíduos com Evento (status == 1):
  if (status[i] == 1) {
    # Amostramos U(F_pop(y-1), F_pop(y))
    lower_bound <- F_pop_y_minus_1[i]
    upper_bound <- F_pop_y[i]
    u_hat[i] <- runif(1,lower_bound,upper_bound)
  }
  # 2. Indivíduos Censurados (status == 0):
  else {
    # Amostramos U(F_pop(y), 1)
    lower_bound <- F_pop_y[i]
    upper_bound <- 1
    u_hat[i] <- runif(1,lower_bound,upper_bound)
  }
}

# 3. Resíduo Quantílico (r_Q_i)
# Função quantil da Normal Padrão: qnorm()
R_Quantile <- qnorm(u_hat)

#
Graph16=hnp(R_Quantile , halfnormal = F, xlab = 'Percentil da N(0,1)', ylab = 'Resíduos',main = 'Gráfico Normal de Probabilidades')

df <- data.frame(
  x = Graph16$x,
  lower = Graph16$lower,
  median = Graph16$median,
  upper = Graph16$upper,
  residuals = Graph16$residuals
)

env=ggplot(df, aes(x)) + geom_ribbon(aes(ymin = lower, ymax = upper),alpha = 0.4) +geom_point(aes(y = residuals),size=0.5) + geom_line(aes(y = median)) + ylab("Resíduos Quantílicos") + xlab("Quantis Teóricos")

env
