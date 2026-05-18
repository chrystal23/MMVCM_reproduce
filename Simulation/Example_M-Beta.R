
# rm(list = ls())
# gc()
# 

func_path = '/path/to/your/code/' # './MMVCM_functions'

source(paste0(func_path, 'simulation_function.R'))
source(paste0(func_path, 'compare.R'))

library('Rcpp')
sourceCpp(paste0(func_path, 'CPPlwls2d_s1.cpp'))
sourceCpp(paste0(func_path, 'CPPlwls2d_s2.cpp'))
sourceCpp(paste0(func_path, 'CPPlwls2d_s2_LR.cpp'))

closeAllConnections()




library(mvtnorm)

# Generate n samples from a 2D Beta distribution with Gaussian copula dependence
r2beta_gausscop <- function(n,
                            shape1 = c(0.8, 0.8),  # s1: c(0.5, 0.5) s2: c(0.8, 0.8)  s3: c(2, 2)
                            shape2 = c(1, 1.2),  # s1: c(1, 1.5)  s2: c(1, 1.2)  s3: c(3, 4)
                            rho = 0.5) {
  # Check correlation
  if (abs(rho) >= 1) stop("rho must be in (-1, 1)")
  
  # Correlation matrix for Gaussian copula
  Sigma <- matrix(c(1, rho,
                    rho, 1), nrow = 2)
  
  # Step 1: latent Gaussian samples
  z <- rmvnorm(n, mean = c(0, 0), sigma = Sigma)
  
  # Step 2: convert to uniforms
  u <- pnorm(z)
  
  # Step 3: convert to Beta margins
  x1 <- qbeta(u[, 1], shape1 = shape1[1], shape2 = shape2[1])
  x2 <- qbeta(u[, 2], shape1 = shape1[2], shape2 = shape2[2])
  
  out <- cbind(x1, x2)
  colnames(out) <- c("x1", "x2")
  return(out)
}





### Example M-Beta

load('exM_model.RData')

runtime = 200



###### n = 500


# specify bandwidth candidates
bw.seq1 = seq(0.01, 0.2, by = 0.01)
bw.seq2 = seq(0.15, 0.30, by = 0.03)
closeAllConnections()

cat('\n\n#############################################
#############################################
Example 3  n = 500', '\n\n')

Res_ex3_n500 = NULL

for (rid in 1:runtime) {
  
  cat('\n\n############################################# Ex3 n500  NO.', rid, '\n')
  ### make data
  
  dat <- makedata_fun(n = 500, alpfun = model$alpfun, beta = model$beta,
                      rdist = r2beta_gausscop, t_corr = T,
                      covrate = 0.2, noisd = sqrt(0.5))
  
  tin <- dat$tin
  yin <- dat$yin
  xin <- dat$xin
  zin <- dat$zin
  
  
  ### estimation
  
  resCP = try( VCMCP(tin = tin, yin = yin, xin = xin, zin = zin, wi = NULL,
                        zeta = zeta.seq, cutoff = max, alpha = 0.05, kernel='epan',
                        bw.seq1 = bw.seq1, bw.seq2 = bw.seq2,
                        NbGrid = 51, nRegGrid = 101, kFolds = 5,
                        npoly = 1, nder = 0, refined = TRUE,
                        hkappa = 1.2, method_opt='average',
                        h2_sepcv = F) )
  
  Res_ex3_n500[[rid]] = resCP
  
  
}

save(Res_ex3_n500, file = 'Res_Ex3_tBeta.RData')





###### n = 1000

bw.seq1 = seq(0.01, 0.2, by = 0.01)
bw.seq2 = seq(0.15, 0.3, by = 0.03)
zeta.seq = NULL 

cat('\n\n#############################################
#############################################
Example 3  n = 1000', '\n\n')

Res_ex3_n1000 = NULL

for (rid in 1:runtime) {

  cat('\n\n############################################# Ex3 n1000  NO.', rid, '\n')
  ### make data

  dat <- makedata_fun(n = 1000, alpfun = model$alpfun, beta = model$beta,
                      rdist = r2beta_gausscop, t_corr = T,
                         covrate = 0.2, noisd = sqrt(0.5))

  tin <- dat$tin
  yin <- dat$yin
  xin <- dat$xin
  zin <- dat$zin


  ### estimation

  resCP = try( VCMCP(tin = tin, yin = yin, xin = xin, zin = zin, wi = NULL,
                   zeta = zeta.seq, cutoff = max, alpha = 0.05, kernel='epan',
                   bw.seq1 = bw.seq1, bw.seq2 = bw.seq2,
                   NbGrid = 51, nRegGrid = 101, kFolds = 5,
                   npoly = 1, nder = 0, refined = TRUE,
                   hkappa = 1.2, method_opt='average',
                   h2_sepcv = F) )

  Res_ex3_n1000[[rid]] = resCP


}


save(Res_ex3_n500, Res_ex3_n1000, file = 'Res_Ex3_tBeta.RData')





###### n = 2000

bw.seq1 = seq(0.01, 0.2, by = 0.01)
bw.seq2 = seq(0.10, 0.25, by = 0.03)

cat('\n\n#############################################
#############################################
Example 3  n = 2000', '\n\n')

Res_ex3_n2000 = NULL

for (rid in 1:runtime) {

  cat('\n\n############################################# Ex3 n2000  NO.', rid, '\n')
  ### make data

  dat <- makedata_fun(n = 2000, alpfun = model$alpfun, beta = model$beta,
                      rdist = r2beta_gausscop, t_corr = T,
                         covrate = 0.2, noisd = sqrt(0.5))

  tin <- dat$tin
  yin <- dat$yin
  xin <- dat$xin
  zin <- dat$zin


  ### estimation

  resCP = try( VCMCP(tin = tin, yin = yin, xin = xin, zin = zin, wi = NULL,
                        zeta = zeta.seq, cutoff = max, alpha = 0.05, kernel='epan',
                        bw.seq1 = bw.seq1, bw.seq2 = bw.seq2,
                        NbGrid = 51, nRegGrid = 101, kFolds = 5,
                        npoly = 1, nder = 0, refined = TRUE,
                        hkappa = 1.6, method_opt='average',
                        h2_sepcv = F) )

  Res_ex3_n2000[[rid]] = resCP


}


save(Res_ex3_n500, Res_ex3_n1000, Res_ex3_n2000, file = 'Res_Ex3_tBeta.RData')








############### Add Testing Results



Res = Res_ex3_n2000 

for (i in 1:length(Res)) {
  
  resCP = Res[[i]]
  
  dat.te <- makedata_fun_2d(n = dim(resCP$obsGrid)[1]/5, m = model$m, alpfun = model$alpfun, beta = model$beta, 
                            covrate = 0.2, noisd = sqrt(0.5))
  
  tin.te <- dat.te$tin
  yin.te <- dat.te$yin
  xin.te <- dat.te$xin
  zin.te <- dat.te$zin
  
  ordt = order(tin.te[,1])
  tin.te = tin.te[ordt, , drop=F]
  yin.te = yin.te[ordt]
  xin.te = xin.te[ordt, , drop=F]
  zin.te = zin.te[ordt, , drop=F]
  
  resCP$tin.te = tin.te
  resCP$yin.te = yin.te
  resCP$xin.te = xin.te
  resCP$zin.te = zin.te
  
  library('gstat')
  library('sp')
  cord_idw = resCP$obsGrid
  alp_idw = sapply(resCP$alp_est, function(x) x$alp.hat)
  
  # grd_idw = expand.grid(t1 = seq(min(cord_idw[,1]), max(cord_idw[,1]), length.out=101),
  #                       t2 = seq(min(cord_idw[,2]), max(cord_idw[,2]), length.out=101))
  grd_idw = as.data.frame(tin.te)
  colnames(grd_idw) = c('t1', 't2')
  coordinates(grd_idw) = c('t1', 't2')
  
  alp_plot = c()  # alp.te
  
  for (j in 1:dim(alp_idw)[2]) {
    dat_idw = cbind(cord_idw, alp_idw[,j])
    dat_idw = dat_idw[!is.na(dat_idw[,3]),]
    dat_idw = as.data.frame(dat_idw)
    colnames(dat_idw) = c('t1', 't2', 'alpha')
    
    # lattice::cloud(alpha~t1*t2, data = dat_idw, zlab=list(rot=90), scales = list(arrows=F))
    
    coordinates(dat_idw) = c('t1', 't2')
    new_alp = idw(alpha ~ 1, dat_idw, grd_idw, nmax=dim(resCP$obsGrid)[1]/100)
    alp_plot = cbind(alp_plot, new_alp$var1.pred)
    
    resCP$alp_est[[j]]$alp_te.hat = new_alp$var1.pred
    
    # pdat = as.data.frame(cbind(tin.te, new_alp$var1.pred))
    # colnames(pdat) = c('t1', 't2', 'alpha')
    # lattice::cloud(alpha~t1*t2, data = pdat, zlab=list(rot=90), scales = list(arrows=F))
  }
  
  colnames(alp_plot) = paste0('alpha', 1:dim(alp_idw)[2]-1)
  resCP$alp_te.hat = alp_plot
  
  alpha_te.hat <- alp_plot
  yhat.te <- c(zin.te %*% resCP$beta_hat + rowSums(alpha_te.hat * cbind(1, xin.te)))
  
  yid = which(!is.na(yhat.te))
  
  mse.te <- mean((yhat.te[yid]-yin.te[yid])^2)
  resCP$mse.te = mse.te
  
  Res[[i]] = resCP
  
  print(mse.te)
  print(i)
  
}

Res_ex3_n2000 = Res


