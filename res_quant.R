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

env=ggplot(df, aes(x)) + geom_ribbon(aes(ymin = lower, ymax = upper),alpha = 0.4) +geom_point(aes(y = residuals),size=0.5) + geom_line(aes(y = median)) + ylab("Resíduos Quantílicos") + xlab("Quantis teóricos")

env
