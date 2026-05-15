library(ggplot2)
library(MASS)
library(sandwich)

logit <- function(x){log(x/(1-x))}
expit <- function(x){exp(x)/(1+exp(x))}

cov_logit_theta_gamma12_mt <- readRDS("cov_logit_theta_gamma12_mt.rds")
data_tr_lt <- readRDS("data_tr_lt.rds")

load("functions.RData")

rho_samp=1

n_iter <- 20000

############################################################
# stage 2: estimate other parameter via MCMC
############################################################

data_lt <- list(C_vt = data_tr_lt$C_vt,
                N_vt = data_tr_lt$N_vt,
                Bstar_vt = data_tr_lt$Bstar_vt,
                mu_vt = data_tr_lt$mu_vt,
                k_vt = data_tr_lt$k_vt,
                t_vt = data_tr_lt$t_vt,
                sins_vt = data_tr_lt$sins_vt, 
                coss_vt = data_tr_lt$coss_vt,
                kappa_est = data_tr_lt$kappa_est,
                kappa_se = data_tr_lt$kappa_se,
                rho_est = data_tr_lt$rho_est)

init_lt <- list(I_vt = round(rep(data_lt$C_vt/2, each = 2)/data_tr_lt$rho_samp_vt[rho_samp]*0.95),
                theta = 0.1, 
                gamma_vt = c(0, 0, 0, 0), 
                betaEN = 0,
                phi = 1,
                q_SIA = 0,
                kappa = data_tr_lt$kappa_samp_vt[rho_samp])

I_vt_init_agg <- init_lt$I_vt[1:(length(init_lt$I_vt)-1)] + init_lt$I_vt[2:(length(init_lt$I_vt))]
I_vt_init_agg <- I_vt_init_agg[seq(1, length(init_lt$I_vt), 2)]
sum(I_vt_init_agg < data_lt$C_vt)

sig_par_lt <- list(sig_I_1 = 250, sig_I_t = 250,
                   sig_phi = 0.65,
                   sig_gamma_vt = c(0.3, 0.005, 0.4, 0.4),
                   sig_betaEN = 0.65,
                   sig_q_SIA = 0.3,
                   sig_theta = 0.04)

############################################################

pmt=proc.time()

# get mcmc samples
seed <- 12345
set.seed(seed)

Ts <- length(data_lt$N_vt)

# create matrix to store posterior samples
phi_sample_vt <- rep(NA, n_iter)
gamma_sample_mt <- matrix(NA, nrow = 4, ncol = n_iter)
betaEN_sample_vt <- rep(NA, n_iter)
theta_sample_vt <- rep(NA, n_iter)
I_sample_mt <- matrix(NA, nrow = Ts, ncol = n_iter)
S_sample_mt <- matrix(NA, nrow = Ts, ncol = n_iter)
q_SIA_sample_vt <- rep(NA, n_iter)
kappa_sample_vt <- rep(NA, n_iter)
rho_sample_vt <- rep(NA, n_iter)

# fill in initial values
phi_sample_vt[1] <- init_lt$phi
gamma_sample_mt[, 1] <- init_lt$gamma_vt
betaEN_sample_vt[1] <- init_lt$betaEN
theta_sample_vt[1] <-  init_lt$theta
I_sample_mt[, 1] <- init_lt$I_vt
S_sample_mt[, 1] <- get_S_vt(S_1 = init_lt$theta*data_lt$N_vt[1], 
                             q_SIA = init_lt$q_SIA, 
                             mu_vt = data_lt$mu_vt, 
                             k_vt = data_lt$k_vt, 
                             Bstar_vt = data_lt$Bstar_vt, 
                             I_vt = init_lt$I_vt)
q_SIA_sample_vt[1] <- init_lt$q_SIA
kappa_sample_vt[1] <- init_lt$kappa
rho_sample_vt[1] <- 1/kappa_sample_vt[1]

# acceptance rate tracker
accpt_phi <- 0
accpt_gamma_vt <- rep(0, 4)
accpt_betaEN <- 0
accpt_theta <- 0
accpt_I_vt <- rep(0, Ts)
accpt_q_SIA <- 0

for(iter in 2:n_iter){
  # iter <- 2
  
  # get curr log_dens_vt
  log_dens_vt_curr <- get_log_dens_vt(phi = phi_sample_vt[iter-1], 
                                      gamma_vt = gamma_sample_mt[, iter-1], 
                                      sins_vt = data_lt$sins_vt, 
                                      coss_vt = data_lt$coss_vt, 
                                      betaEN = betaEN_sample_vt[iter-1], 
                                      I_vt = I_sample_mt[, iter-1], 
                                      S_vt = S_sample_mt[, iter-1], 
                                      N_vt = data_lt$N_vt)
  
  # update phi
  update_phi_lt <- update_phi(log_dens_vt_curr = log_dens_vt_curr, 
                              phi_curr = phi_sample_vt[iter-1], 
                              gamma_vt = gamma_sample_mt[, iter-1], 
                              sins_vt = data_lt$sins_vt, 
                              coss_vt = data_lt$coss_vt, 
                              betaEN = betaEN_sample_vt[iter-1], 
                              I_vt = I_sample_mt[, iter-1],
                              S_vt = S_sample_mt[, iter-1],
                              N_vt = data_lt$N_vt,
                              sig_phi = sig_par_lt$sig_phi)
  
  phi_sample_vt[iter] <- update_phi_lt$phi
  accpt_phi <- accpt_phi + as.integer(phi_sample_vt[iter] != phi_sample_vt[iter-1])
  log_dens_vt_curr <- update_phi_lt$log_dens_vt
  
  # update gamma_vt
  gamma_sample_mt[, iter] <- gamma_sample_mt[, iter-1]
  for (j in 3:4){
    # j <- 1
    update_gamma_j_lt <- update_gamma_j(log_dens_vt_curr = log_dens_vt_curr, 
                                        j = j, 
                                        gamma_vt_curr = gamma_sample_mt[, iter], 
                                        phi = phi_sample_vt[iter], 
                                        sins_vt = data_lt$sins_vt, 
                                        coss_vt = data_lt$coss_vt, 
                                        betaEN = betaEN_sample_vt[iter-1], 
                                        I_vt = I_sample_mt[, iter-1], 
                                        S_vt = S_sample_mt[, iter-1], 
                                        N_vt = data_lt$N_vt, 
                                        sig_gamma_vt = sig_par_lt$sig_gamma_vt)
    
    gamma_sample_mt[j, iter] <- update_gamma_j_lt$gamma_j
    log_dens_vt_curr <- update_gamma_j_lt$log_dens_vt
  }
  accpt_gamma_vt[3:4] <- accpt_gamma_vt[3:4] + as.integer(gamma_sample_mt[3:4, iter] != gamma_sample_mt[3:4, iter-1])
  
  # update betaEN
  update_betaEN_lt <- update_betaEN(log_dens_vt_curr = log_dens_vt_curr, 
                                    betaEN_curr = betaEN_sample_vt[iter-1], 
                                    phi = phi_sample_vt[iter], 
                                    gamma_vt = gamma_sample_mt[, iter], 
                                    sins_vt = data_lt$sins_vt, 
                                    coss_vt = data_lt$coss_vt, 
                                    I_vt = I_sample_mt[, iter-1], 
                                    S_vt = S_sample_mt[, iter-1], 
                                    N_vt = data_lt$N_vt, 
                                    sig_betaEN = sig_par_lt$sig_betaEN)
  
  betaEN_sample_vt[iter] <- update_betaEN_lt$betaEN
  accpt_betaEN <- accpt_betaEN + as.integer(betaEN_sample_vt[iter] != betaEN_sample_vt[iter-1])
  log_dens_vt_curr <- update_betaEN_lt$log_dens_vt
  
  # update q_SIA
  update_q_SIA_lt <- update_q_SIA(log_dens_vt_curr = log_dens_vt_curr, 
                                  q_SIA_curr = q_SIA_sample_vt[iter-1],
                                  phi = phi_sample_vt[iter], 
                                  gamma_vt = gamma_sample_mt[, iter], 
                                  sins_vt = data_lt$sins_vt, 
                                  coss_vt = data_lt$coss_vt,
                                  betaEN = betaEN_sample_vt[iter], 
                                  mu_vt = data_lt$mu_vt, 
                                  k_vt = data_lt$k_vt, 
                                  Bstar_vt = data_lt$Bstar_vt, 
                                  I_vt = I_sample_mt[, iter-1],
                                  N_vt = data_lt$N_vt, 
                                  theta = theta_sample_vt[iter-1],
                                  sig_q_SIA = sig_par_lt$sig_q_SIA)
  
  q_SIA_sample_vt[iter] <- update_q_SIA_lt$q_SIA
  accpt_q_SIA <- accpt_q_SIA + as.integer(q_SIA_sample_vt[iter] != q_SIA_sample_vt[iter-1])
  log_dens_vt_curr <- update_q_SIA_lt$log_dens_vt
  
  # update theta, gamma_1, gamma_2
  update_theta_gamma_12_lt <- update_theta_gamma_12(log_dens_vt_curr = log_dens_vt_curr, 
                                  theta_curr = theta_sample_vt[iter-1], 
                                  gamma_1_curr = gamma_sample_mt[1, iter], 
                                  gamma_2_curr = gamma_sample_mt[2, iter], 
                                  phi = phi_sample_vt[iter], 
                                  gamma_vt = gamma_sample_mt[, iter], 
                                  sins_vt = data_lt$sins_vt, 
                                  coss_vt = data_lt$coss_vt,
                                  betaEN = betaEN_sample_vt[iter], 
                                  I_vt = I_sample_mt[, iter-1],
                                  N_vt = data_lt$N_vt, 
                                  q_SIA = q_SIA_sample_vt[iter],
                                  mu_vt = data_lt$mu_vt, 
                                  k_vt = data_lt$k_vt, 
                                  Bstar_vt = data_lt$Bstar_vt, 
                                  # sig_theta = sig_par_lt$sig_theta,
                                  # sig_gamma_1 = sig_par_lt$sig_gamma_vt[1],
                                  # sig_gamma_2 = sig_par_lt$sig_gamma_vt[2],
                                  cov_logit_theta_gamma12_mt = 0.5*cov_logit_theta_gamma12_mt)
  
  theta_sample_vt[iter] <- update_theta_gamma_12_lt$theta
  gamma_sample_mt[1, iter] <- update_theta_gamma_12_lt$gamma_1
  gamma_sample_mt[2, iter] <- update_theta_gamma_12_lt$gamma_2
  log_dens_vt_curr <- update_theta_gamma_12_lt$log_dens_vt
  
  accpt_theta <- accpt_theta + as.integer(theta_sample_vt[iter] != theta_sample_vt[iter-1])
  accpt_gamma_vt[1] <- accpt_gamma_vt[1] + as.integer(gamma_sample_mt[1, iter] != gamma_sample_mt[1, iter-1])
  accpt_gamma_vt[2] <- accpt_gamma_vt[2] + as.integer(gamma_sample_mt[2, iter] != gamma_sample_mt[2, iter-1])
  
  # update S_vt
  S_sample_mt[, iter] <- get_S_vt(S_1 = theta_sample_vt[iter]*data_lt$N_vt[1], 
                                  q_SIA = q_SIA_sample_vt[iter], 
                                  mu_vt = data_lt$mu_vt, 
                                  k_vt = data_lt$k_vt, 
                                  Bstar_vt = data_lt$Bstar_vt, 
                                  I_vt = I_sample_mt[, iter-1])
  
  # update kappa_vt & rho_vt
  kappa_sample_vt[iter] <- data_tr_lt$kappa_samp_vt[rho_samp]
  rho_sample_vt[iter] <- min(1/kappa_sample_vt[iter], 0.99)
  
  # update I_vt
  I_sample_mt[, iter] <- I_sample_mt[, iter-1]
  
  # update I_1
  update_I_1_lt <- update_I_1(log_dens_vt_curr = log_dens_vt_curr, 
                              I_vt_curr = I_sample_mt[, iter], 
                              S_vt_curr = S_sample_mt[, iter], 
                              rho = rho_sample_vt[iter], 
                              q_SIA = q_SIA_sample_vt[iter], 
                              mu_vt = data_lt$mu_vt, 
                              k_vt = data_lt$k_vt, 
                              gamma_vt = gamma_sample_mt[, iter], 
                              sins_vt = data_lt$sins_vt, 
                              coss_vt = data_lt$coss_vt,  
                              phi = phi_sample_vt[iter], 
                              betaEN = betaEN_sample_vt[iter], 
                              Bstar_vt = data_lt$Bstar_vt, 
                              N_vt = data_lt$N_vt, 
                              C_vt = data_lt$C_vt, 
                              sig_I_1 = sig_par_lt$sig_I_1)
  
  I_sample_mt[1, iter] <- update_I_1_lt$I_1
  S_sample_mt[, iter] <- update_I_1_lt$S_vt
  log_dens_vt_curr <- update_I_1_lt$log_dens_vt
  
  # update I_2, ..., I_T
  for (t in 2:Ts){
    # t <- 2
    update_I_t_lt <- update_I_t(log_dens_vt_curr = log_dens_vt_curr, 
                                I_vt_curr = I_sample_mt[, iter], 
                                S_vt_curr = S_sample_mt[, iter], 
                                t = t, 
                                rho = rho_sample_vt[iter], 
                                q_SIA = q_SIA_sample_vt[iter], 
                                mu_vt = data_lt$mu_vt, 
                                k_vt = data_lt$k_vt, 
                                gamma_vt = gamma_sample_mt[, iter], 
                                sins_vt = data_lt$sins_vt, 
                                coss_vt = data_lt$coss_vt, 
                                phi = phi_sample_vt[iter], 
                                betaEN = betaEN_sample_vt[iter], 
                                Bstar_vt = data_lt$Bstar_vt, 
                                N_vt = data_lt$N_vt, 
                                C_vt = data_lt$C_vt, 
                                sig_I_t = sig_par_lt$sig_I_t)
    
    I_sample_mt[t, iter] <- update_I_t_lt$I_t
    S_sample_mt[, iter] <- update_I_t_lt$S_vt
    log_dens_vt_curr <- update_I_t_lt$log_dens_vt
  }
  accpt_I_vt <- accpt_I_vt + as.integer(I_sample_mt[, iter] != I_sample_mt[, iter-1])
  
  if(iter%%100 == 0){
    print(paste("iter = ", iter, sep = ""))
  }
}

proc.time()-pmt

# save workspace
save.image(file = paste0("results_rhosamp", rho_samp, "_image.RData"))

