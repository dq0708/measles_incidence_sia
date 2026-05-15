library(ggplot2)
library(MASS)
library(sandwich)

logit <- function(x){log(x/(1-x))}
expit <- function(x){exp(x)/(1+exp(x))}

load("example_data.RData")

############################################################

Ts <- length(c(data_tr_lt$N_vt, data_test_lt$N_vt))
Ts_tr <- length(data_tr_lt$N_vt)
C_vt <- c(data_tr_lt$C_vt, data_test_lt$C_vt)
t_SIA <- 71
N_vt <- c(data_tr_lt$N_vt, data_test_lt$N_vt)
Bstar_vt <- c(data_tr_lt$Bstar_vt, data_test_lt$Bstar_vt)

############################################################
# stage 1: estimate rho via wls
############################################################

T_wls <- (t_SIA-1)/2
Bstar_vt_tr <- data_tr_lt$Bstar_vt
C_vt_tr <- data_tr_lt$C_vt

Bstar_cum <- rep(0, T_wls)
C_cum <- rep(0, T_wls)

for (m in 1:T_wls){
  Bstar_cum[m] <- sum(Bstar_vt_tr[1:(2*m)])
  C_cum[m] <- sum(C_vt_tr[1:m])
}

m_wls <- c(7:35)
Bstar_cum_final <- Bstar_cum[m_wls]- Bstar_cum[m_wls[1]-1]
C_cum_final <- C_cum[m_wls]- C_cum[m_wls[1]-1]

###

ols_fit <- lm(Bstar_cum_final ~ 1 + C_cum_final)
sandwich_se <- diag(vcovHC(ols_fit , type = "HC"))^0.5

kappa_est <- as.numeric(coef(ols_fit)[2])
kappa_se <- as.numeric(sandwich_se[2])
kappa_low <- as.numeric(kappa_est - 1.96*kappa_se)
kappa_up <- as.numeric(kappa_est + 1.96*kappa_se)

rho_est <- as.numeric(1/kappa_est)
rho_low <- as.numeric(1/kappa_up)
rho_up <- as.numeric(1/kappa_low)

# print(c(rho_est, rho_low, rho_up))

###

data_tr_lt$rho_est <- rho_est
data_tr_lt$rho_low <- rho_low
data_tr_lt$rho_up <- rho_up

data_tr_lt$kappa_est <- kappa_est
data_tr_lt$kappa_low <- kappa_low
data_tr_lt$kappa_up <- kappa_up
data_tr_lt$kappa_se <- kappa_se

##############################
# draw random samples from estimated rho distribution
##############################

seed <- 12345
set.seed(seed)

n_rho_samp <- 100

kappa_samp_vt <- rnorm(n_rho_samp, kappa_est, sd = kappa_se)
rho_samp_vt <- 1/kappa_samp_vt
rho_samp_vt[rho_samp_vt >= 0.99] <- 0.99

data_tr_lt$kappa_samp_vt <- kappa_samp_vt
data_tr_lt$rho_samp_vt <- rho_samp_vt

##############################

saveRDS(data_tr_lt, "data_tr_lt.rds")
saveRDS(data_test_lt, "data_test_lt.rds")

