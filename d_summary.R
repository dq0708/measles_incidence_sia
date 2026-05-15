rm(list=ls())

library(data.table)
library(xtable)
library(ggplot2)
library(scoringRules)

set.seed(12345)

logit <- function(x){log(x/(1-x))}
expit <- function(x){exp(x)/(1+exp(x))}

########################################
# save all post samples
########################################

n_rho_samp <- 100

comp_dt <- post_I_sample_mt <- post_S_sample_mt<- NULL

for (rho_samp in 1:n_rho_samp){
  # rho_samp <- 1
  
  load(paste0("results_rhosamp", rho_samp, "_image.RData"))
  
  n_burn <- n_iter/2 + 1
  
  comp_dt <- rbind(comp_dt,
                       data.table(rho_samp = rho_samp,
                                  rho = as.numeric(rho_sample_vt[seq(n_burn, n_iter, 5)]),
                                  theta = as.numeric(theta_sample_vt[seq(n_burn, n_iter, 5)]),
                                  gamma_1 = as.numeric(gamma_sample_mt[1, seq(n_burn, n_iter, 5)]),
                                  gamma_2 = as.numeric(gamma_sample_mt[2, seq(n_burn, n_iter, 5)]),
                                  gamma_3 = as.numeric(gamma_sample_mt[3, seq(n_burn, n_iter, 5)]),
                                  gamma_4 = as.numeric(gamma_sample_mt[4, seq(n_burn, n_iter, 5)]),
                                  betaEN = as.numeric(betaEN_sample_vt[seq(n_burn, n_iter, 5)]),
                                  phi = as.numeric(phi_sample_vt[seq(n_burn, n_iter, 5)]),
                                  q_SIA = as.numeric(q_SIA_sample_vt[seq(n_burn, n_iter, 5)])))
  
  post_I_sample_mt <- rbind(post_I_sample_mt, 
                       t(as.matrix(I_sample_mt[, seq(n_burn, n_iter, 5)])))
  
  post_S_sample_mt <- rbind(post_S_sample_mt, 
                       t(as.matrix(S_sample_mt[, seq(n_burn, n_iter, 5)])))
  
  print(rho_samp)
  
}

########################################
# posterior table
########################################

mean_col <- c(mean(as.numeric(comp_dt[, rho])),
              mean(as.numeric(comp_dt[, theta])),
              mean(as.numeric(comp_dt[, gamma_1])),
              mean(as.numeric(comp_dt[, gamma_2])),
              mean(as.numeric(comp_dt[, gamma_3])),
              mean(as.numeric(comp_dt[, gamma_4])),
              mean(as.numeric(comp_dt[, betaEN])),
              mean(as.numeric(comp_dt[, q_SIA])), 
              mean(as.numeric(comp_dt[, phi])))

median_col <- c(median(as.numeric(comp_dt[, rho])),
                median(as.numeric(comp_dt[, theta])),
                median(as.numeric(comp_dt[, gamma_1])),
                median(as.numeric(comp_dt[, gamma_2])),
                median(as.numeric(comp_dt[, gamma_3])),
                median(as.numeric(comp_dt[, gamma_4])),
                median(as.numeric(comp_dt[, betaEN])),
                median(as.numeric(comp_dt[, q_SIA])), 
                median(as.numeric(comp_dt[, phi])))

lower_col <- c(quantile(as.numeric(comp_dt[, rho]), 0.025),
               quantile(as.numeric(comp_dt[, theta]), 0.025),
               quantile(as.numeric(comp_dt[, gamma_1]), 0.025),
               quantile(as.numeric(comp_dt[, gamma_2]), 0.025),
               quantile(as.numeric(comp_dt[, gamma_3]), 0.025),
               quantile(as.numeric(comp_dt[, gamma_4]), 0.025),
               quantile(as.numeric(comp_dt[, betaEN]), 0.025),
               quantile(as.numeric(comp_dt[, q_SIA]), 0.025), 
               quantile(as.numeric(comp_dt[, phi]), 0.025))

upper_col <- c(quantile(as.numeric(comp_dt[, rho]), 0.975),
               quantile(as.numeric(comp_dt[, theta]), 0.975),
               quantile(as.numeric(comp_dt[, gamma_1]), 0.975),
               quantile(as.numeric(comp_dt[, gamma_2]), 0.975),
               quantile(as.numeric(comp_dt[, gamma_3]), 0.975),
               quantile(as.numeric(comp_dt[, gamma_4]), 0.975),
               quantile(as.numeric(comp_dt[, betaEN]), 0.975),
               quantile(as.numeric(comp_dt[, q_SIA]), 0.975), 
               quantile(as.numeric(comp_dt[, phi]), 0.975))

param_col <- c("rho", "theta", "gamma_1", "gamma_2", "gamma_3", "gamma_4", 
               "beta_EN", "q_SIA", "phi")

summary_table <- data.frame(parameter = param_col,
                            # post_mean = mean_col,
                            post_med = median_col,
                            lower = lower_col,
                            upper = upper_col)


