



fcfun_LR_EUR_USD<-function (series, horizon, nsim) {
   X <- series
   X[, c("US_CPI", "EZ_CPI")] <- log(X[, c("US_CPI", "EZ_CPI")])
   UIP_ST = X[, "EZ_ST"] - X[, "US_ST"]
   PPP_CPI = X[, "EUR_USD"] + X[, "US_CPI"] - X[, "EZ_CPI"]
   mods <- list(UIP_ST = lm(UIP_ST ~ 1), PPP_CPI = lm(PPP_CPI ~ 
                                                         1))
   ect <- do.call(cbind, lapply(mods, function(m) {
      ts(residuals(m), start = start(X), frequency = frequency(X))
   }))
   ce_vars <- lapply(mods, function(m) c(1, -coef(m)))
   ec_lag <- window(lag(ect, k = -1), end = end(ect))
   Beta <- matrix(0, NCOL(ect), NCOL(X), dimnames = list(colnames(ect), 
                                                         colnames(X)))
   Beta["UIP_ST", c("US_ST", "EZ_ST")] <- c(-1, 1)
   Beta["PPP_CPI", c("EZ_CPI", "US_CPI", "EUR_USD")] <- c(-1, 
                                                          1, 1)
   p <- 3
   constraint_LR <- macroForecasts:::makeVarxConstraintMatrix(X, p, ec_lag, 
                                                              include.mean = TRUE)
   constraint_LR[paste(rep(c("EZ_CPI_", "EZ_ST_"), p), 1:p, 
                       sep = ""), c("F1", "US_CPI", "US_ST")] <- FALSE
   constraint_LR[c("ect_PPP_CPI", "ect_UIP_ST"), "F1"] <- FALSE
   constraint_LR["drift", c("F1", "US_ST", "EUR_USD")] <- FALSE
   w <- 300
   utils::capture.output(VEC_test <- MTS::VARX(utils::tail(diff(X), 
                                                           w), utils::tail(ec_lag, w), p = p, m = 0, include.mean = TRUE, 
                                               fixed = constraint_LR))
   LT <- VEC_test$Ph0
   Alphas <- VEC_test$beta
   Pi <- Alphas %*% Beta
   Phi_comp <- macroForecasts:::makeVarxPhiCompanionMatrix(Pi, VEC_test$Phi)
   constants_ce <- sapply(ce_vars, function(x) unname(x["(Intercept)"]))
   trend_ce <- c(0, 0)
   mu_constant <- Alphas %*% constants_ce
   delta <- Alphas %*% trend_ce
   u_t_VEC <- ts(residuals(VEC_test), end = end(X), frequency = frequency(X))
   spec_F1 <- rugarch::ugarchspec(variance.model = list(garchOrder = c(1, 
                                                                       1), model = "sGARCH"), mean.model = list(armaOrder = c(0, 
                                                                                                                              0), include.mean = FALSE), distribution.model = "sged")
   spec_US_ST <- rugarch::ugarchspec(variance.model = list(garchOrder = c(1, 
                                                                          1), model = "sGARCH"), mean.model = list(armaOrder = c(0, 
                                                                                                                                 0), include.mean = FALSE), distribution.model = "sged")
   spec_EZ_ST <- rugarch::ugarchspec(variance.model = list(garchOrder = c(1, 
                                                                          1), model = "sGARCH"), mean.model = list(armaOrder = c(0, 
                                                                                                                                 0), include.mean = FALSE), distribution.model = "sstd")
   mod_F1 <- rugarch::ugarchfit(spec_F1, u_t_VEC[, "F1"], solver = "hybrid")
   mod_US_ST <- rugarch::ugarchfit(spec_US_ST, u_t_VEC[, "US_ST"], 
                                   solver = "hybrid")
   mod_EZ_ST <- rugarch::ugarchfit(spec_EZ_ST, u_t_VEC[, "EZ_ST"], 
                                   solver = "hybrid")
   z_t_VEC <- scale(u_t_VEC, center = TRUE, scale = TRUE)
   z_F1 <- ts(mod_F1@fit$z, start = start(u_t_VEC), frequency = frequency(u_t_VEC))
   z_US_ST <- ts(mod_US_ST@fit$z, start = start(u_t_VEC), frequency = frequency(u_t_VEC))
   z_EZ_ST <- ts(mod_EZ_ST@fit$z, start = start(u_t_VEC), frequency = frequency(u_t_VEC))
   z_t_VEC[, c("F1", "US_ST", "EZ_ST")] <- window(cbind(z_F1, 
                                                        z_US_ST, z_EZ_ST), start = start(u_t_VEC))
   z_t_VEC <- scale(z_t_VEC, center = TRUE, scale = TRUE)
   eps_z <- z_t_VEC[sample.int(NROW(z_t_VEC), size = nsim * 
                                  horizon, replace = TRUE), ]
   eps_z <- aperm(array(eps_z, dim = c(nsim, horizon, NCOL(z_t_VEC)), 
                        dimnames = list(NULL, NULL, colnames(z_t_VEC))), c(1, 
                                                                           3, 2))
   innov_VEC <- array(0, dim = dim(eps_z), dimnames = dimnames(eps_z))
   innov_VEC[, "F1", ] <- t(rugarch::ugarchsim(mod_F1, n.sim = horizon, 
                                               m.sim = nsim, custom.dist = list(name = "sample", distfit = t(eps_z[, 
                                                                                                                   "F1", ])), startMethod = "sample")@simulation$seriesSim)
   innov_VEC[, "US_ST", ] <- t(rugarch::ugarchsim(mod_US_ST, 
                                                  n.sim = horizon, m.sim = nsim, custom.dist = list(name = "sample", 
                                                                                                    distfit = t(eps_z[, "US_ST", ])), startMethod = "sample")@simulation$seriesSim)
   innov_VEC[, "EZ_ST", ] <- t(rugarch::ugarchsim(mod_EZ_ST, 
                                                  n.sim = horizon, m.sim = nsim, custom.dist = list(name = "sample", 
                                                                                                    distfit = t(eps_z[, "EZ_ST", ])), startMethod = "sample")@simulation$seriesSim)
   for (z in colnames(innov_VEC)) {
      if (!(z %in% c("F1", "US_ST", "EZ_ST"))) {
         innov_VEC[, z, ] <- eps_z[, z, ] * sd(u_t_VEC[, 
                                                       z])
      }
   }
   k_VEC <- NCOL(X)
   store_VEC_full <- matrix(0, nsim, k_VEC * (horizon + p + 
                                                 1))
   store_VEC_full[, 1:((p + 1) * k_VEC)] <- rep(t(utils::tail(X, 
                                                              n = p + 1)), each = nsim)
   sim_store <- array(0, dim = c(nsim, NCOL(X), horizon), dimnames = list(NULL, 
                                                                          colnames(X), NULL))
   ELB <- -0.0075
   for (j in 1:horizon) {
      store_VEC_full[, ((j + p) * k_VEC + 1):((j + p) * k_VEC + 
                                                 k_VEC)] <- rep(LT + mu_constant + delta * (NROW(X) + 
                                                                                               j), each = nsim) + store_VEC_full[, ((j - 1) * k_VEC + 
                                                                                                                                       1):((j + p - 1) * k_VEC + k_VEC)] %*% t(Phi_comp) + 
         innov_VEC[, , j]
      store_VEC_full[, (j + p) * k_VEC + 3] <- pmax(store_VEC_full[, 
                                                                   (j + p) * k_VEC + 3], ELB)
      store_VEC_full[, (j + p) * k_VEC + 5] <- pmax(store_VEC_full[, 
                                                                   (j + p) * k_VEC + 5], ELB)
      sim_store[, , j] <- store_VEC_full[, ((j + p) * k_VEC + 
                                               1):((j + p) * k_VEC + k_VEC)]
      sim_store[, c("US_CPI", "EZ_CPI"), j] <- exp(sim_store[, 
                                                             c("US_CPI", "EZ_CPI"), j])
   }
   sim_store
}

getEUR<-function (horizon, nsim, historicalFrom = NULL,end_month=NULL) {
   
   
   # a) Variables:
   
   US_ST<-window(macroForecasts::makeUniformMonthly(inqr::getMacrobondTimeseries("usrate0204")),start=1980)/ 100
   EZ_ST_d<-window(macroForecasts::makeUniformMonthly(inqr::getMacrobondTimeseries("de3mgov")),start=1980)/ 100
   EZ_ST<-window(macroForecasts::makeUniformMonthly(inqr::getMacrobondTimeseries("de3mgov_m")),start=1980)/ 100
   EUR_USD<-window(macroForecasts::makeUniformMonthly(inqr::getMacrobondTimeseries("eur")),start=1980)
   EZ_CPI_NSA<-window(macroForecasts::makeUniformMonthly(inqr::getMacrobondTimeseries("eupric0001")),start=1980)
   US_CPI<-window(macroForecasts::makeUniformMonthly(inqr::getMacrobondTimeseries("uspric2156")),start=1980)
   
   
   EZ_CPI <- seasonal::final(seasonal::seas(EZ_CPI_NSA,transform.function="log",arima.model=c(9,1,0,0,1,1),regression.variables=c("td","easter[5]"),outlier=NULL,x11=""))
   US_CPI <- macroForecasts:::forecastMissingUS_CPI(na.omit(US_CPI),until=end(EUR_USD))
   
   # Re-index to 2010 = 100
   
   EZ_CPI_2010 <- EZ_CPI / mean(window(EZ_CPI,start=2010,end=c(2010,12))) * 100
   EZ_CPI_2010<-macroForecasts:::forecastMissingUS_CPI(EZ_CPI_2010,until=end(EUR_USD))
   US_CPI_2010 <- US_CPI / mean(window(US_CPI,start=2010,end=c(2010,12))) * 100
   US_CPI_2010<-macroForecasts:::forecastMissingUS_CPI(na.omit(US_CPI_2010),until=end(EUR_USD))
   
   F1 <- macroForecasts::F1global_alt()
   
   
   
   Endo <- cbind(F1,US_CPI_2010,US_ST,EZ_CPI_2010,EZ_ST,EUR_USD)
   Endo[nrow(Endo),"EZ_ST"]<-EZ_ST_d[length(EZ_ST_d)]
   colnames(Endo) <- c("F1","US_CPI","US_ST","EZ_CPI","EZ_ST","EUR_USD")
   
   
   Endo <- window(Endo, start = 1980,end(EUR_USD))
   
   if (is.null(end_month)==FALSE){
      Endo<-window(Endo,end=end_month)
      
   }
   
   
   
   
   
   Endo<-imputeTS::na_kalman(Endo,"auto.arima")
   
   
   
   if (is.null(historicalFrom)) {
      prob <- macroForecasts:::fcfun_LR_EUR_USD(Endo, horizon, nsim)[, "EUR_USD", 
                                                                     , drop = FALSE]
      return(tsims::tsims(Endo[, "EUR_USD", drop = FALSE], 
                          prob))
   }
   else {
      if (is.na(historicalFrom)) {
         historicalFrom <- as.Date("2002-01-01")
      }
      first_fc <- which(zoo::as.Date(Endo) == historicalFrom)
      fcfun_LR_EUR_USD.part <- roostr::simwrap.partialsave(macroForecasts:::fcfun_LR_EUR_USD, 
                                                           Endo[, colnames(Endo) != "EUR_USD"], colnames(Endo))
      perf <- roostr::computeErrorDistributions(Endo[, "EUR_USD", 
                                                     drop = FALSE], fcfun_LR_EUR_USD.part, first_fc, 
                                                horizon, nsim)
      return(perf)
   }
}

EUR<-getEUR(60,20000,historicalFrom = NULL,end_month = NULL)

plot(EUR)
