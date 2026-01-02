library(survival)
library(nlme)     
library(MASS)     
library(stats)    
library(extraDistr) # Para rdunif e rdunif

set.seed(42)

# =================================================================
# SECTION 0: PARÂMETROS E CONFIGURAÇÃO
# =================================================================

PARAM_FIXOS_ADJ <- list(
  beta_sigma_int = 0.577998, 
  beta_mu_int    = 2.753863,
  beta_mu_x1     = 1.937027,
  beta_mu_x2     = -0.887691,
  beta_p0_int    = 1.688387,
  beta_p0_x1     = -1.206960,
  beta_p0_x2     = 0.420279,
  beta_phi_int   = 2.297517 
)

THETA_TRUE <- unlist(PARAM_FIXOS_ADJ)
N_PARAM <- length(THETA_TRUE)
SAMPLE_SIZES <- c(500) # Teste focado em N=500
N_MONTE_CARLO <- 1000
P_SUSC_TARGET <- 0.05

# =================================================================
# SECTION 1: FUNÇÕES DE SOBREVIVÊNCIA E BUSCA DO TAU
# =================================================================

# S_latente conforme seu código que funcionou bem
S_latente <- function(w, mu, sigma) {
  1 - pnorm((log(pmax(w-1, 1e-7)) - mu) / sigma)
}

prob_censura_estavel <- function(tau, mu, sigma, n_size) {
  tau_int <- floor(tau)
  if (tau_int < 1) return(1)
  y_vals <- 1:tau_int
  # P(Y > C) = 1/tau * sum( S_pop(y)^n )
  # Aqui usamos a S_latente para o cálculo do tau
  sy <- S_latente(y_vals, mu, sigma)
  p_cens <- (1 / tau_int) * sum(exp(n_size * log(pmax(sy, 1e-100))))
  return(p_cens)
}

find_tau_robust <- function(mu, sigma, p_susc, n_size) {
  res <- tryCatch({
    uniroot(function(t) prob_censura_estavel(t, mu, sigma, n_size) - p_susc,
            interval = c(1, 15000), extendInt = "yes", tol = 1e-5)$root
  }, error = function(e) 2000)
  return(res)
}

# =================================================================
# SECTION 2: GERAÇÃO DE DADOS (CONFORME CÓDIGO DE SUCESSO)
# =================================================================

gen_data_lnd_nb <- function(N, params, p_susc) {
  x1 <- sample(c(0, 1), size = N, replace = TRUE)
  x2 <- sample(c(0, 1), size = N, replace = TRUE)
  Xe <- cbind(1, x1, x2)
  
  sigma <- exp(params[1])
  mu    <- as.numeric(Xe %*% params[2:4])
  p0    <- 1 / (1 + exp(-Xe %*% params[5:7]))
  phi   <- exp(params[8])
  
  eta <- (pmax(p0, 1e-7)^(-phi) - 1) / phi
  M   <- rnbinom(N, size = 1/phi, prob = 1/(1 + phi*eta))
  
  W <- rep(Inf, N)
  idx_susc <- which(M > 0)
  for (i in idx_susc) {
    U <- runif(M[i])
    # Lógica rlnd_disc: ceiling(exp(...)) - 1
    W[i] <- ceiling(exp(mu[i] + sigma * qnorm(1 - U))) - 1
  }
  W <- pmax(0, W)
  
  tau_i <- find_tau_robust(mean(mu), sigma, p_susc, N)
  # Censura Uniforme Discreta [1, Tau]
  C <- rdunif(N, 1, max(1, round(tau_i)))
  
  Y_obs <- pmin(W, C)
  Delta <- as.numeric(W <= C & W < Inf)
  Y_obs[Y_obs == Inf] <- C[Y_obs == Inf]
  
  return(data.frame(y = Y_obs, status = Delta, x1 = x1, x2 = x2))
}

# =================================================================
# SECTION 3: VEROSSIMILHANÇA (ESTRUTURA QUE FUNCIONOU BEM)
# =================================================================

fvero <- function(theta, data) {
  Xe <- cbind(1, data$x1, data$x2)
  
  sigma <- exp(theta[1])
  mu    <- as.numeric(Xe %*% theta[2:4])
  p0    <- 1 / (1 + exp(-as.numeric(Xe %*% theta[5:7])))
  p0    <- pmin(pmax(p0, 1e-6), 0.9999)
  phi   <- exp(theta[8])
  
  eta <- (p0^(-phi) - 1) / phi
  
  # Spop1 = S_pop(y) | Spop2 = S_pop(y+1)
  vF1 <- 1 - S_latente(data$y, mu, sigma)
  Spop1 <- (1 + phi * eta * vF1)^(-1 / phi)
  
  vF2 <- 1 - S_latente(data$y + 1, mu, sigma)
  Spop2 <- (1 + phi * eta * vF2)^(-1 / phi)
  
  # P(Y=y)
  fpop <- pmax(Spop1 - Spop2, 1e-16)
  
  # Likelihood Contribution: status=1 -> fpop | status=0 -> Spop1
  # EXATAMENTE COMO NO SEU CÓDIGO ANTERIOR
  loglik <- sum(data$status * log(fpop) + (1 - data$status) * log(pmax(Spop1, 1e-16)))
  
  if(!is.finite(loglik)) return(-1e20)
  return(loglik)
}

ml_optim <- function(theta_init, data) {
  # Uso do BFGS direto como no seu código de sucesso
  mgg <- optim(theta_init, fvero, data = data, method = "BFGS",
               control = list(fnscale = -1, maxit = 1000, reltol = 1e-10))
  
  if (mgg$convergence != 0) return(list(converged = FALSE))
  
  # Hessiana e ginv (como no código que funcionou)
  h <- tryCatch(fdHess(mgg$par, fvero, data = data)$Hessian, error = function(e) NA)
  if(any(is.na(h))) return(list(converged = FALSE))
  
  covmat <- tryCatch(ginv(-h), error = function(e) matrix(NA, N_PARAM, N_PARAM))
  setheta <- sqrt(pmax(0, diag(covmat)))
  
  if(any(is.na(setheta)) || any(setheta > 15)) return(list(converged = FALSE))
  
  return(list(par = mgg$par, se = setheta, converged = TRUE))
}

# =================================================================
# SECTION 4: LOOP DE EXECUÇÃO E MÉTRICAS
# =================================================================

# O loop gerará a tabela df_resultados com Bias, RMSE e CP.
# Execute para N=500 e observe a queda drástica do Bias em beta_mu_int.
# =================================================================
# SECTION 4: LOOP DE SIMULAÇÃO
# =================================================================

resultados <- list()

for (N in SAMPLE_SIZES) {
  cat(sprintf("\n--- Executando N = %d (%d repetições) ---\n", N, N_MONTE_CARLO))
  
  estimativas_rep <- matrix(NA, nrow = N_MONTE_CARLO, ncol = N_PARAM)
  se_rep <- matrix(NA, nrow = N_MONTE_CARLO, ncol = N_PARAM)
  censura_obs <- numeric(N_MONTE_CARLO)
  valid_runs <- 0
  
  for (run in 1:N_MONTE_CARLO) {
    set.seed(N + run)
    sim_data <- gen_data_lnd_nb(N, THETA_TRUE, P_SUSC)
    censura_obs[run] <- 1 - mean(sim_data$status)
    
    ml_result <- ml_optim(THETA_TRUE, sim_data)
    
    if (ml_result$converged) {
      valid_runs <- valid_runs + 1
      estimativas_rep[valid_runs, ] <- ml_result$par
      se_rep[valid_runs, ] <- ml_result$se
    }
    if (run %% 200 == 0) cat(sprintf("  Rep %d. Válidas: %d. Censura Média: %.2f%%\n", 
                                     run, valid_runs, mean(censura_obs[1:run])*100))
  }
  
  est_clean <- estimativas_rep[1:valid_runs, ]
  se_clean <- se_rep[1:valid_runs, ]
  
  # Métricas
  bias   <- colMeans(est_clean) - THETA_TRUE
  rmse   <- sqrt(colMeans(sweep(est_clean, 2, THETA_TRUE, "-")^2))
  avg_se <- colMeans(se_clean)
  aw     <- 2 * qnorm(0.975) * avg_se
  
  lower_ci <- est_clean - qnorm(0.975) * se_clean
  upper_ci <- est_clean + qnorm(0.975) * se_clean
  cp <- colMeans(sweep(lower_ci, 2, THETA_TRUE, "<") & sweep(upper_ci, 2, THETA_TRUE, ">"))
  
  resultados[[as.character(N)]] <- data.frame(
    N = N,
    Parameter = names(THETA_TRUE),
    True_Value = THETA_TRUE,
    Bias = bias,
    RMSE = rmse,
    Avg_SE = avg_se,
    AW = aw,
    CP = cp,
    Valid_Runs = valid_runs,
    Avg_Cens_Obs = mean(censura_obs)
  )
}

# Finalização
df_resultados <- do.call(rbind, resultados)
print(df_resultados)
