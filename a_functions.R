#####
get_St <- function(t, S_vt, q_SIA, mu_vt, k_vt, Bstar_vt, I_vt){
  
  St <- round(S_vt[t-1] + Bstar_vt[t] - I_vt[t] - q_SIA*mu_vt[t]*S_vt[k_vt[t]])
  return(St)
}

#####
get_S_vt <- function(S_1, q_SIA, mu_vt, k_vt, Bstar_vt, I_vt){
  Ts <- length(I_vt)
  S_vt <- rep(NA, Ts)
  S_vt[1] <- S_1
  for (t in 2:Ts){
    S_vt[t] <- get_St(t = t, 
                      S_vt = S_vt, 
                      q_SIA = q_SIA, 
                      mu_vt = mu_vt, 
                      k_vt = k_vt, 
                      Bstar_vt = Bstar_vt, 
                      I_vt = I_vt)
  }
  return(S_vt)
}

#####
get_lambdat <- function(t, I_vt, S_vt, N_vt, betaEN, gamma_vt, sins_vt, coss_vt, 
                        alpha = 0.975){
  
  lambdat <- exp(gamma_vt[1] + gamma_vt[2]*t + gamma_vt[3]*sins_vt[t] + gamma_vt[4]*coss_vt[t]) * I_vt[t-1]^alpha * S_vt[t-1] / N_vt[t-1] + exp(betaEN) * N_vt[t-1]
  
  if(lambdat < 0){
    print(c(t))
  }
  return(lambdat)
}

#####
get_log_dens_vt <- function(phi, gamma_vt, sins_vt, coss_vt, betaEN, I_vt, S_vt, N_vt, 
                            alpha = 0.975){
  Ts <- length(N_vt)
  lambda_vt <- log_dens_vt <- rep(NA, Ts)
  for (t in 2:Ts){
    lambda_vt[t] <- get_lambdat(t = t, 
                                I_vt = I_vt, 
                                S_vt = S_vt, 
                                N_vt = N_vt, 
                                betaEN = betaEN, 
                                gamma_vt = gamma_vt, 
                                sins_vt = sins_vt, 
                                coss_vt = coss_vt)
    log_dens_vt[t] <- dnbinom(x = I_vt[t], mu = lambda_vt[t], size = phi, log = T)
  }
  return(log_dens_vt)
}

#####
update_phi <- function(log_dens_vt_curr, phi_curr, 
                       gamma_vt, sins_vt, coss_vt, betaEN, I_vt, S_vt, N_vt, 
                       sig_phi, alpha = 0.975){
  
  Ts <- length(N_vt)
  
  log_post_dens_phi_curr <- sum(log_dens_vt_curr[2:Ts])
  
  phi_prop <- exp(rnorm(1, log(phi_curr), sig_phi))
  log_dens_vt_prop <- get_log_dens_vt(phi = phi_prop, 
                                      gamma_vt = gamma_vt, 
                                      sins_vt = sins_vt, 
                                      coss_vt = coss_vt, 
                                      betaEN = betaEN, 
                                      I_vt = I_vt, 
                                      S_vt = S_vt, 
                                      N_vt = N_vt)
  log_post_dens_phi_prop <- sum(log_dens_vt_prop[2:Ts])
  
  p <- min(exp(log_post_dens_phi_prop - log_post_dens_phi_curr), 1)
  a <- runif(1)
  if (a < p){
    phi_new <- phi_prop
    log_dens_vt_new <- log_dens_vt_prop
  }else{
    phi_new <- phi_curr
    log_dens_vt_new <- log_dens_vt_curr
  }
  return(list(phi = phi_new,
              log_dens_vt = log_dens_vt_new))
}

#####
update_gamma_j <- function(log_dens_vt_curr, j, gamma_vt_curr, 
                           phi, sins_vt, coss_vt, betaEN, I_vt, S_vt, N_vt, 
                           sig_gamma_vt, alpha = 0.975){
  
  Ts <- length(N_vt)
  
  log_post_dens_gamma_j_curr <- sum(log_dens_vt_curr[2:Ts]) + dnorm(x = gamma_vt_curr[j], mean = 0, sd = 20, log = T)
  
  gamma_j_curr <- gamma_vt_curr[j]
  gamma_j_prop <- rnorm(1, gamma_j_curr, sig_gamma_vt[j])
  gamma_vt_prop <- gamma_vt_curr
  gamma_vt_prop[j] <- gamma_j_prop
  
  log_dens_vt_prop <- get_log_dens_vt(phi = phi, 
                                      gamma_vt = gamma_vt_prop, 
                                      sins_vt = sins_vt, 
                                      coss_vt = coss_vt, 
                                      betaEN = betaEN, 
                                      I_vt = I_vt, 
                                      S_vt = S_vt, 
                                      N_vt = N_vt)
  log_post_dens_gamma_j_prop <- sum(log_dens_vt_prop[2:Ts]) + dnorm(x = gamma_vt_prop[j], mean = 0, sd = 20, log = T)
  
  p <- min(exp(log_post_dens_gamma_j_prop - log_post_dens_gamma_j_curr), 1)
  a <- runif(1)
  if (a < p){
    gamma_j_new <- gamma_j_prop
    log_dens_vt_new <- log_dens_vt_prop
  }else{
    gamma_j_new <- gamma_j_curr
    log_dens_vt_new <- log_dens_vt_curr
  }
  return(list(gamma_j = gamma_j_new,
              log_dens_vt = log_dens_vt_new))
}

#####
update_betaEN <- function(log_dens_vt_curr, betaEN_curr, 
                          phi, gamma_vt, sins_vt, coss_vt, I_vt, S_vt, N_vt, 
                          sig_betaEN, alpha = 0.975){
  
  Ts <- length(N_vt)
  
  log_post_dens_betaEN_curr <- sum(log_dens_vt_curr[2:Ts]) + dnorm(x = betaEN_curr, mean = 0, sd = 20, log = T)
  
  betaEN_prop <- rnorm(1, betaEN_curr, sig_betaEN)
  log_dens_vt_prop <- get_log_dens_vt(phi = phi, 
                                      gamma_vt = gamma_vt, 
                                      sins_vt = sins_vt, 
                                      coss_vt = coss_vt, 
                                      betaEN = betaEN_prop, 
                                      I_vt = I_vt, 
                                      S_vt = S_vt, 
                                      N_vt = N_vt)
  log_post_dens_betaEN_prop <- sum(log_dens_vt_prop[2:Ts]) + dnorm(x = betaEN_prop, mean = 0, sd = 20, log = T)
  
  p <- min(exp(log_post_dens_betaEN_prop - log_post_dens_betaEN_curr), 1)
  a <- runif(1)
  if (a < p){
    betaEN_new <- betaEN_prop
    log_dens_vt_new <- log_dens_vt_prop
  }else{
    betaEN_new <- betaEN_curr
    log_dens_vt_new <- log_dens_vt_curr
  }
  return(list(betaEN = betaEN_new,
              log_dens_vt = log_dens_vt_new))
}

#####
update_theta_gamma_12 <- function(log_dens_vt_curr, theta_curr, gamma_1_curr, gamma_2_curr,
                                 phi, gamma_vt, sins_vt, coss_vt, betaEN, I_vt, N_vt,
                                 q_SIA, mu_vt, k_vt, Bstar_vt, 
                                 cov_logit_theta_gamma12_mt, alpha = 0.975){
  
  Ts <- length(N_vt)
  
  log_post_dens_theta_gamma_12_curr <- sum(log_dens_vt_curr[2:Ts]) + 
    dbeta(x = theta_curr, shape1 = 2, shape2 = 18, log = T) + 
    dnorm(x = gamma_1_curr, mean = 0, sd = 20, log = T) + 
    dnorm(x = gamma_2_curr, mean = 0, sd = 20, log = T)
  
  Sigma_prop <- cov_logit_theta_gamma12_mt
  
  logit_theta_gamma12_prop <- mvrnorm(mu = c(logit(theta_curr), gamma_1_curr, gamma_2_curr),
                                      Sigma = Sigma_prop)
  
  theta_prop <- expit(logit_theta_gamma12_prop[1])
  gamma_1_prop <- logit_theta_gamma12_prop[2]
  gamma_2_prop <- logit_theta_gamma12_prop[3]
  
  S_1_prop <-  N_vt[1] * theta_prop
  S_vt_prop <- get_S_vt(S_1 = S_1_prop,
                        q_SIA = q_SIA , 
                        mu_vt = mu_vt, 
                        k_vt = k_vt, 
                        Bstar_vt = Bstar_vt, 
                        I_vt = I_vt)
  
  if(sum(c(min(S_vt_prop) <= 0, S_1_prop > N_vt[1], S_vt_prop[1:(Ts-1)] < I_vt[2:Ts], S_vt_prop > N_vt)) > 0){
    theta_new <- theta_curr
    gamma_1_new <- gamma_1_curr
    gamma_2_new <- gamma_2_curr
    log_dens_vt_new <- log_dens_vt_curr
  }else{
    gamma_vt[1] <- gamma_1_prop
    gamma_vt[2] <- gamma_2_prop
    log_dens_vt_prop <- get_log_dens_vt(phi = phi, 
                                        gamma_vt = gamma_vt, 
                                        sins_vt = sins_vt, 
                                        coss_vt = coss_vt, 
                                        betaEN = betaEN, 
                                        I_vt = I_vt, 
                                        S_vt = S_vt_prop, 
                                        N_vt = N_vt)
    log_post_dens_theta_gamma_12_prop <- sum(log_dens_vt_prop[2:Ts]) + 
      dbeta(x = theta_prop, shape1 = 2, shape2 = 18, log = T) + 
      dnorm(x = gamma_1_prop, mean = 0, sd = 20, log = T) + 
      dnorm(x = gamma_2_prop, mean = 0, sd = 20, log = T)
    
    p <- min(exp(log_post_dens_theta_gamma_12_prop - log_post_dens_theta_gamma_12_curr), 1)
    a <- runif(1)
    if (a < p){
      theta_new <- theta_prop
      gamma_1_new <- gamma_1_prop
      gamma_2_new <- gamma_2_prop
      log_dens_vt_new <- log_dens_vt_prop
    }else{
      theta_new <- theta_curr
      gamma_1_new <- gamma_1_curr
      gamma_2_new <- gamma_2_curr
      log_dens_vt_new <- log_dens_vt_curr
    }
  }
  
  return(list(theta = theta_new,
              gamma_1 = gamma_1_new,
              gamma_2 = gamma_2_new,
              log_dens_vt = log_dens_vt_new))
}

#####
update_I_1 <- function(log_dens_vt_curr, I_vt_curr, S_vt_curr, 
                       rho, q_SIA, mu_vt, k_vt, gamma_vt, sins_vt, coss_vt,
                       phi, betaEN, Bstar_vt, N_vt, C_vt, 
                       sig_I_1, alpha = 0.975){
  Ts <- length(N_vt)
  
  I_1_curr <- I_vt_curr[1]
  log_post_dens_I_1_curr <- dbinom(x = C_vt[1], size = I_1_curr + I_vt_curr[2], prob = rho, log = T) + sum(log_dens_vt_curr[2:Ts])
  
  delta <- sample(size = 1, x = c(1:sig_I_1))
  I_1_prop <- sample(size = 1, x = c(I_1_curr-delta, I_1_curr+delta))
  I_vt_prop <- I_vt_curr
  I_vt_prop[1] <- I_1_prop
  S_vt_prop <- get_S_vt(S_1 = S_vt_curr[1],
                        q_SIA = q_SIA, 
                        mu_vt = mu_vt, 
                        k_vt = k_vt, 
                        Bstar_vt = Bstar_vt, 
                        I_vt = I_vt_prop)
  
  I_vt_prop_agg <- I_vt_prop[1:(length(I_vt_prop)-1)] + I_vt_prop[2:(length(I_vt_prop))]
  I_vt_prop_agg <- I_vt_prop_agg[seq(1, length(I_vt_prop), 2)]
  
  if(sum(c(I_vt_prop < 0, S_vt_prop[1:(Ts-1)] < I_vt_prop[2:Ts], S_vt_prop > N_vt, min(S_vt_prop) <= 0, I_vt_prop_agg < C_vt)) > 0){
    I_1_new <- I_1_curr
    log_dens_vt_new <- log_dens_vt_curr
    S_vt_new <- S_vt_curr
  }else{
    log_dens_vt_prop <- get_log_dens_vt(phi = phi, 
                                        gamma_vt = gamma_vt, 
                                        sins_vt = sins_vt, 
                                        coss_vt = coss_vt, 
                                        betaEN = betaEN, 
                                        I_vt = I_vt_prop, 
                                        S_vt = S_vt_prop, 
                                        N_vt = N_vt)
    log_post_dens_I_1_prop <- dbinom(x = C_vt[1], size = I_1_prop + I_vt_prop[2], prob = rho, log = T) + sum(log_dens_vt_prop[2:Ts])
    
    p <- min(exp(log_post_dens_I_1_prop - log_post_dens_I_1_curr), 1)
    a <- runif(1)
    if (a < p){
      I_1_new <- I_1_prop
      log_dens_vt_new <- log_dens_vt_prop
      S_vt_new <- S_vt_prop
    }else{
      I_1_new <- I_1_curr
      log_dens_vt_new <- log_dens_vt_curr
      S_vt_new <- S_vt_curr
    }
  }
  
  return(list(I_1 = I_1_new,
              log_dens_vt = log_dens_vt_new,
              S_vt = S_vt_new))
}

#####
update_I_t <- function(log_dens_vt_curr, I_vt_curr, S_vt_curr, t, 
                       rho, q_SIA, mu_vt, k_vt, gamma_vt, sins_vt, coss_vt,
                       phi, betaEN, Bstar_vt, N_vt, C_vt, 
                       sig_I_t, alpha = 0.975){
  
  Ts <- length(N_vt)
  
  I_t_curr <- I_vt_curr[t]
  
  if (t%%2 == 0){
    log_post_dens_I_t_curr <- dbinom(x = C_vt[t/2], size = I_t_curr + I_vt_curr[t-1], prob = rho, log = T) + sum(log_dens_vt_curr[t:Ts])
  }else{
    log_post_dens_I_t_curr <- dbinom(x = C_vt[(t+1)/2], size = I_t_curr + I_vt_curr[t+1], prob = rho, log = T) + sum(log_dens_vt_curr[t:Ts])
  }
  
  delta <- sample(size = 1, x = c(1:sig_I_t))
  I_t_prop <- sample(size = 1, x = c(I_t_curr-delta, I_t_curr+delta))
  I_vt_prop <- I_vt_curr
  I_vt_prop[t] <- I_t_prop
  S_vt_prop <- S_vt_curr
  for (j in t:Ts){
    S_vt_prop[j] <- get_St(t = j, 
                           S_vt = S_vt_prop, 
                           q_SIA = q_SIA, 
                           mu_vt = mu_vt, 
                           k_vt = k_vt, 
                           Bstar_vt = Bstar_vt, 
                           I_vt = I_vt_prop)
  }
  
  I_vt_prop_agg <- I_vt_prop[1:(length(I_vt_prop)-1)] + I_vt_prop[2:(length(I_vt_prop))]
  I_vt_prop_agg <- I_vt_prop_agg[seq(1, length(I_vt_prop), 2)]
  
  if(sum(c(I_vt_prop < 0, I_vt_prop_agg < C_vt, S_vt_prop[1:(Ts-1)] < I_vt_prop[2:Ts], S_vt_prop > N_vt, min(S_vt_prop) <= 0)) > 0){
    I_t_new <- I_t_curr
    log_dens_vt_new <- log_dens_vt_curr
    S_vt_new <- S_vt_curr
  }else{
    log_dens_vt_prop <- log_dens_vt_curr
    lambda_vt_prop <-  rep(NA, Ts)
    for (j in t:Ts){
      lambda_vt_prop[j] <- get_lambdat(t = j, 
                                       I_vt = I_vt_prop, 
                                       S_vt = S_vt_prop, 
                                       N_vt = N_vt, 
                                       betaEN = betaEN, 
                                       gamma_vt = gamma_vt, 
                                       sins_vt = sins_vt, 
                                       coss_vt = coss_vt)
      log_dens_vt_prop[j] <- dnbinom(x = I_vt_prop[j], mu = lambda_vt_prop[j], size = phi, log = T)
    }
    if (t%%2 == 0){
      log_post_dens_I_t_prop <- dbinom(x = C_vt[t/2], size = I_t_prop + I_vt_prop[t-1], prob = rho, log = T) + sum(log_dens_vt_prop[t:Ts])
    }else{
      log_post_dens_I_t_prop <- dbinom(x = C_vt[(t+1)/2], size = I_t_prop + I_vt_prop[t+1], prob = rho, log = T) + sum(log_dens_vt_prop[t:Ts])
    }
    
    p <- min(exp(log_post_dens_I_t_prop - log_post_dens_I_t_curr), 1)
    a <- runif(1)
    if (a < p){
      I_t_new <- I_t_prop
      log_dens_vt_new <- log_dens_vt_prop
      S_vt_new <- S_vt_prop
    }else{
      I_t_new <- I_t_curr
      log_dens_vt_new <- log_dens_vt_curr
      S_vt_new <- S_vt_curr
    }
  }
  
  return(list(I_t = I_t_new,
              log_dens_vt = log_dens_vt_new,
              S_vt = S_vt_new))
}

#####
update_q_SIA <- function(log_dens_vt_curr, q_SIA_curr,
                         phi, gamma_vt, sins_vt, coss_vt, betaEN, 
                         mu_vt, k_vt, Bstar_vt, I_vt, N_vt, theta,
                         sig_q_SIA, alpha = 0.975){
  
  Ts <- length(N_vt)
  
  log_post_dens_q_SIA_curr <- sum(log_dens_vt_curr[2:Ts]) + dbeta(x = q_SIA_curr, shape1 = 2, shape2 = 2, log = T)
  
  q_SIA_prop <- rnorm(1, q_SIA_curr, sig_q_SIA)
  S_vt_prop <- get_S_vt(S_1 = theta*N_vt[1], 
                        q_SIA = q_SIA_prop, 
                        mu_vt = mu_vt, 
                        k_vt = k_vt, 
                        Bstar_vt = Bstar_vt, 
                        I_vt = I_vt)
  
  if(sum(c(q_SIA_prop < 0, q_SIA_prop > 1, S_vt_prop[1:(Ts-1)] < I_vt[2:Ts], S_vt_prop > N_vt, min(S_vt_prop) <= 0)) > 0){
    q_SIA_new <- q_SIA_curr
    log_dens_vt_new <- log_dens_vt_curr
  }else{
    log_dens_vt_prop <- get_log_dens_vt(phi = phi, 
                                        gamma_vt = gamma_vt, 
                                        sins_vt = sins_vt, 
                                        coss_vt = coss_vt, 
                                        betaEN = betaEN, 
                                        I_vt = I_vt, 
                                        S_vt = S_vt_prop, 
                                        N_vt = N_vt)
    log_post_dens_q_SIA_prop <- sum(log_dens_vt_prop[2:Ts]) + dbeta(x = q_SIA_prop, shape1 = 2, shape2 = 2, log = T)
    
    p <- min(exp(log_post_dens_q_SIA_prop - log_post_dens_q_SIA_curr), 1)
    a <- runif(1)
    if (a < p){
      q_SIA_new <- q_SIA_prop
      log_dens_vt_new <- log_dens_vt_prop
    }else{
      q_SIA_new <- q_SIA_curr
      log_dens_vt_new <- log_dens_vt_curr
    }
  }
  return(list(q_SIA = q_SIA_new,
              log_dens_vt = log_dens_vt_new))
}

#####
save.image(file = "functions.RData")

