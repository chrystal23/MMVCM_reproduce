
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



#################################################################################

## Example 2

load(file = 'exS_model.RData')

runtime = 200



###### n = 500

bw.seq1 = seq(0.01, 0.2, by = 0.01)
bw.seq2 = seq(0.1, 0.2, length.out = 6)
zeta.seq = seq(0.4, 0.6, length.out=6)

cat('\n\n#############################################
#############################################
Example 2  n = 500', '\n\n')

# Res_ex2_n500 = NULL

for (rid in 1:runtime) {
  
  cat('\n\n############################################# Ex2 n500  NO.', rid, '\n')
  
  ### make data
  dat <- makedata_fun(n = 500, m = model$m, alpfun = model$alpfun, beta = model$beta, 
                         covrate = 0.2, noisd = sqrt(0.5))
  
  tin <- dat$tin
  yin <- dat$yin
  xin <- dat$xin
  zin <- dat$zin
  
  # h2_sepcv = ifelse(rid<=9, F, T)
  
  resCP = VCMCP(tin = tin, yin = yin, xin = xin, zin = zin, wi = NULL, 
                   zeta = zeta.seq, cutoff = max, alpha = 0.05, kernel='epan',
                   bw.seq1 = bw.seq1, bw.seq2 = bw.seq2, 
                   NbGrid = 201, nRegGrid = 201, kFolds = 5, 
                   npoly = 1, nder = 0, refined = TRUE,
                   hkappa = hkappa, method_opt='average',
                   h2_sepcv = h2_sepcv)
  
  Res_ex2_n500[[rid]] = resCP
  
}

save(Res_ex2_n500, file = 'Res_Ex2.RData')






### Run competing methods

cat('\n\n############################################# 
############################################# 
Example 2 compare  n = 500', '\n\n')

n = 500

Res_SVCoef_ex2_n500 = NULL
Res_VCoef_ex2_n500 = NULL
Res_lm_ex2_n500 = NULL
Res_seqMS_ex2_n500 = NULL

for (rid in 1:runtime) {
  
  cat('\n\n############################################# NO.', rid, '\n')
  ### make data
  
  dat <- makedata_fun_2d(n = n, m = model$m, alpfun = model$alpfun, beta = model$beta,
                         rdist = rtrunc_gaussian_2d, t_corr = T,
                         covrate = 0.2, noisd = sqrt(0.5))
  
  tin <- dat$tin
  yin <- dat$yin
  xin <- dat$xin
  zin <- dat$zin
  
  
  ### estimation
  
  res_SVCoef <- try( SVCoef_func_2d(tin = tin, yin = yin, xin = xin, zin = zin, wi = NULL,
                                    bw.seq1 = bw.seq1, bw.seq2 = bw.seq2, 
                                    nRegGrid = 101, kernel='epan',
                                    npoly = 1, nder = 0, kFolds = 5, 
                                    method_opt='average') )
  res_VCoef <- try( VCoef_func_2d(tin = tin, yin = yin, xin = cbind(xin, zin), wi = NULL,
                                  bw.seq = bw.seq2, 
                                  nRegGrid = 101, kernel='epan',
                                  npoly = 1, nder = 0, kFolds = 5) )
  lmdat <- data.frame(cbind(yin, xin, zin))
  names(lmdat)[1] <- 'Y'
  res_lm <- try( lm(Y~., data = lmdat) )
  res_lm$mse <- mean(res_lm$residuals^2)
  
  res_seqMS <- try( seqMS(dat = cbind(yin, tin[,1], 1, zin, xin)) )
  
  Res_SVCoef_ex2_n500[[rid]] = res_SVCoef
  Res_VCoef_ex2_n500[[rid]] = res_VCoef
  Res_lm_ex2_n500[[rid]] = res_lm
  Res_seqMS_ex2_n500[[rid]] = res_seqMS
  
  
  ######
  
  dat.te <- makedata_fun_2d(n = n/5, m = model$m, alpfun = model$alpfun, beta = model$beta, 
                            rdist = rtrunc_gaussian_2d, t_corr = T,
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
  
  
  ###
  
  resCP = res_SVCoef
  
  resCP$tin.te = tin.te
  resCP$yin.te = yin.te
  resCP$xin.te = xin.te
  resCP$zin.te = zin.te
  
  library('gstat')
  library('sp')
  cord_idw = resCP$obsGrid
  alp_idw = resCP$alpha_hat
  
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
    
    # pdat = as.data.frame(cbind(tin.te, new_alp$var1.pred))
    # colnames(pdat) = c('t1', 't2', 'alpha')
    # lattice::cloud(alpha~t1*t2, data = pdat, zlab=list(rot=90), scales = list(arrows=F))
  }
  
  colnames(alp_plot) = paste0('alpha', 1:dim(alp_idw)[2]-1)
  resCP$alpha_hat.te = alp_plot
  
  alpha_te.hat <- alp_plot
  yhat.te <- c(zin.te %*% resCP$beta_hat + rowSums(alpha_te.hat * cbind(1, xin.te)))
  
  yid = which(!is.na(yhat.te))
  
  mse.te <- mean((yhat.te[yid]-yin.te[yid])^2)
  resCP$mse.te = mse.te
  
  Res_SVCoef_ex2_n500[[rid]] = resCP
  
  ###
  
  resCP = res_VCoef
  
  resCP$tin.te = tin.te
  resCP$yin.te = yin.te
  resCP$xin.te = cbind(xin.te, zin.te)
  
  library('gstat')
  library('sp')
  cord_idw = resCP$obsGrid
  alp_idw = resCP$alpha_hat
  
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
    
    # pdat = as.data.frame(cbind(tin.te, new_alp$var1.pred))
    # colnames(pdat) = c('t1', 't2', 'alpha')
    # lattice::cloud(alpha~t1*t2, data = pdat, zlab=list(rot=90), scales = list(arrows=F))
  }
  
  colnames(alp_plot) = paste0('alpha', 1:dim(alp_idw)[2]-1)
  resCP$alpha_hat.te = alp_plot
  
  alpha_te.hat <- alp_plot
  yhat.te <- rowSums(alpha_te.hat * cbind(1, resCP$xin.te))
  
  yid = which(!is.na(yhat.te))
  
  mse.te <- mean((yhat.te[yid]-yin.te[yid])^2)
  resCP$mse.te = mse.te
  
  Res_VCoef_ex2_n500[[rid]] = resCP
  
  ###
  
  res_seqMS.te <- try( seqMS.predict(dat = cbind(yin.te, tin.te[,1], 1, zin.te, xin.te),
                                     thresholds = res_seqMS$thresholds,
                                     beta = res_seqMS$beta) )
  
  res_seqMS = c(res_seqMS, res_seqMS.te)
  res_seqMS$obsGrid = tin
  
  Res_seqMS_ex2_n500[[rid]] = res_seqMS
}


save(Res_SVCoef_ex2_n500, Res_VCoef_ex2_n500, Res_lm_ex2_n500, Res_seqMS_ex2_n500,
     file = 'Res_Ex2_compare.RData')








###### n = 1000

bw.seq1 = seq(0.01, 0.2, by = 0.01)
bw.seq2 = seq(0.1, 0.15, length.out = 6)
zeta.seq = seq(0.4, 0.6, length.out=6)

cat('\n\n#############################################
#############################################
Example 2  n = 1000', '\n\n')

# Res_ex2_n1000 = NULL

for (rid in 1:runtime) {
  
  cat('\n\n############################################# Ex2 n1000  NO.', rid, '\n')
  
  ### make data
  dat <- makedata_fun(n = 1000, m = model$m, alpfun = model$alpfun, beta = model$beta,
                      covrate = 0.2, noisd = sqrt(0.5))
  
  tin <- dat$tin
  yin <- dat$yin
  xin <- dat$xin
  zin <- dat$zin
  
  # h2_sepcv = ifelse(rid<=9, F, T)
  
  resCP = VCMCP(tin = tin, yin = yin, xin = xin, zin = zin, wi = NULL,
                zeta = zeta.seq, cutoff = max, alpha = 0.05, kernel='epan',
                bw.seq1 = bw.seq1, bw.seq2 = bw.seq2,
                NbGrid = 201, nRegGrid = 201, kFolds = 5,
                npoly = 1, nder = 0, refined = TRUE,
                hkappa = hkappa, method_opt='average',
                h2_sepcv = h2_sepcv)
  
  Res_ex2_n1000[[rid]] = resCP
  
}

save(Res_ex2_n500, Res_ex2_n1000, file = 'Res_Ex2.RData')






### Run competing methods

cat('\n\n############################################# 
############################################# 
Example 2 compare  n = 1000', '\n\n')

n = 1000

Res_SVCoef_ex2_n1000 = NULL
Res_VCoef_ex2_n1000 = NULL
Res_lm_ex2_n1000 = NULL
Res_seqMS_ex2_n1000 = NULL

for (rid in 1:runtime) {
  
  cat('\n\n############################################# NO.', rid, '\n')
  ### make data
  
  dat <- makedata_fun_2d(n = n, m = model$m, alpfun = model$alpfun, beta = model$beta,
                         rdist = rtrunc_gaussian_2d, t_corr = T,
                         covrate = 0.2, noisd = sqrt(0.5))
  
  tin <- dat$tin
  yin <- dat$yin
  xin <- dat$xin
  zin <- dat$zin
  
  
  ### estimation
  
  res_SVCoef <- try( SVCoef_func_2d(tin = tin, yin = yin, xin = xin, zin = zin, wi = NULL,
                                    bw.seq1 = bw.seq1, bw.seq2 = bw.seq2, 
                                    nRegGrid = 101, kernel='epan',
                                    npoly = 1, nder = 0, kFolds = 5, 
                                    method_opt='average') )
  res_VCoef <- try( VCoef_func_2d(tin = tin, yin = yin, xin = cbind(xin, zin), wi = NULL,
                                  bw.seq = bw.seq2, 
                                  nRegGrid = 101, kernel='epan',
                                  npoly = 1, nder = 0, kFolds = 5) )
  lmdat <- data.frame(cbind(yin, xin, zin))
  names(lmdat)[1] <- 'Y'
  res_lm <- try( lm(Y~., data = lmdat) )
  res_lm$mse <- mean(res_lm$residuals^2)
  
  res_seqMS <- try( seqMS(dat = cbind(yin, tin[,1], 1, zin, xin)) )
  
  Res_SVCoef_ex2_n1000[[rid]] = res_SVCoef
  Res_VCoef_ex2_n1000[[rid]] = res_VCoef
  Res_lm_ex2_n1000[[rid]] = res_lm
  Res_seqMS_ex2_n1000[[rid]] = res_seqMS
  
  
  ######
  
  dat.te <- makedata_fun_2d(n = n/5, m = model$m, alpfun = model$alpfun, beta = model$beta, 
                            rdist = rtrunc_gaussian_2d, t_corr = T,
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
  
  
  ###
  
  resCP = res_SVCoef
  
  resCP$tin.te = tin.te
  resCP$yin.te = yin.te
  resCP$xin.te = xin.te
  resCP$zin.te = zin.te
  
  library('gstat')
  library('sp')
  cord_idw = resCP$obsGrid
  alp_idw = resCP$alpha_hat
  
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
    
    # pdat = as.data.frame(cbind(tin.te, new_alp$var1.pred))
    # colnames(pdat) = c('t1', 't2', 'alpha')
    # lattice::cloud(alpha~t1*t2, data = pdat, zlab=list(rot=90), scales = list(arrows=F))
  }
  
  colnames(alp_plot) = paste0('alpha', 1:dim(alp_idw)[2]-1)
  resCP$alpha_hat.te = alp_plot
  
  alpha_te.hat <- alp_plot
  yhat.te <- c(zin.te %*% resCP$beta_hat + rowSums(alpha_te.hat * cbind(1, xin.te)))
  
  yid = which(!is.na(yhat.te))
  
  mse.te <- mean((yhat.te[yid]-yin.te[yid])^2)
  resCP$mse.te = mse.te
  
  Res_SVCoef_ex2_n1000[[rid]] = resCP
  
  ###
  
  resCP = res_VCoef
  
  resCP$tin.te = tin.te
  resCP$yin.te = yin.te
  resCP$xin.te = cbind(xin.te, zin.te)
  
  library('gstat')
  library('sp')
  cord_idw = resCP$obsGrid
  alp_idw = resCP$alpha_hat
  
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
    
    # pdat = as.data.frame(cbind(tin.te, new_alp$var1.pred))
    # colnames(pdat) = c('t1', 't2', 'alpha')
    # lattice::cloud(alpha~t1*t2, data = pdat, zlab=list(rot=90), scales = list(arrows=F))
  }
  
  colnames(alp_plot) = paste0('alpha', 1:dim(alp_idw)[2]-1)
  resCP$alpha_hat.te = alp_plot
  
  alpha_te.hat <- alp_plot
  yhat.te <- rowSums(alpha_te.hat * cbind(1, resCP$xin.te))
  
  yid = which(!is.na(yhat.te))
  
  mse.te <- mean((yhat.te[yid]-yin.te[yid])^2)
  resCP$mse.te = mse.te
  
  Res_VCoef_ex2_n1000[[rid]] = resCP
  
  ###
  
  res_seqMS.te <- try( seqMS.predict(dat = cbind(yin.te, tin.te[,1], 1, zin.te, xin.te),
                                     thresholds = res_seqMS$thresholds,
                                     beta = res_seqMS$beta) )
  
  res_seqMS = c(res_seqMS, res_seqMS.te)
  res_seqMS$obsGrid = tin
  
  Res_seqMS_ex2_n1000[[rid]] = res_seqMS
}


save(Res_SVCoef_ex2_n500, Res_VCoef_ex2_n500, Res_lm_ex2_n500, Res_seqMS_ex2_n500,
     Res_SVCoef_ex2_n1000, Res_VCoef_ex2_n1000, Res_lm_ex2_n1000, Res_seqMS_ex2_n1000,
     file = 'Res_Ex2_compare.RData')













###### n = 2000

bw.seq1 = seq(0.01, 0.2, by = 0.01)
bw.seq2 = seq(0.05, 0.15, length.out = 6)
zeta.seq = seq(0.3, 0.5, length.out=6)

cat('\n\n#############################################
#############################################
Example 2  n = 2000', '\n\n')

Res_ex2_n2000 = NULL

for (rid in 1:runtime) {
  
  cat('\n\n############################################# Ex2 n2000  NO.', rid, '\n')
  
  ### make data
  dat <- makedata_fun(n = 2000, m = model$m, alpfun = model$alpfun, beta = model$beta, 
                         covrate = 0.2, noisd = sqrt(0.5))
  
  tin <- dat$tin
  yin <- dat$yin
  xin <- dat$xin
  zin <- dat$zin
  
  # h2_sepcv = ifelse(rid<=9, F, T)
  
  resCP = VCMCP(tin = tin, yin = yin, xin = xin, zin = zin, wi = NULL, 
                   zeta = zeta.seq, cutoff = max, alpha = 0.05, kernel='epan',
                   bw.seq1 = bw.seq1, bw.seq2 = bw.seq2, 
                   NbGrid = 201, nRegGrid = 201, kFolds = 5, 
                   npoly = 1, nder = 0, refined = TRUE,
                   hkappa = hkappa, method_opt='average',
                   h2_sepcv = h2_sepcv)
  
  Res_ex2_n2000[[rid]] = resCP
  
}

save(Res_ex2_n500, Res_ex2_n1000, Res_ex2_n2000, file = 'Res_Ex2.RData')







### Run competing methods


cat('\n\n############################################# 
############################################# 
Example 2 compare  n = 2000', '\n\n')

n = 2000

Res_SVCoef_ex2_n2000 = NULL
Res_VCoef_ex2_n2000 = NULL
Res_lm_ex2_n2000 = NULL
Res_seqMS_ex2_n2000 = NULL

for (rid in 1:runtime) {
  
  cat('\n\n############################################# NO.', rid, '\n')
  ### make data
  
  dat <- makedata_fun_2d(n = n, m = model$m, alpfun = model$alpfun, beta = model$beta,
                         rdist = rtrunc_gaussian_2d, t_corr = T,
                         covrate = 0.2, noisd = sqrt(0.5))
  
  tin <- dat$tin
  yin <- dat$yin
  xin <- dat$xin
  zin <- dat$zin
  
  
  ### estimation
  
  res_SVCoef <- try( SVCoef_func_2d(tin = tin, yin = yin, xin = xin, zin = zin, wi = NULL,
                                    bw.seq1 = bw.seq1, bw.seq2 = bw.seq2, 
                                    nRegGrid = 101, kernel='epan',
                                    npoly = 1, nder = 0, kFolds = 5, 
                                    method_opt='average') )
  res_VCoef <- try( VCoef_func_2d(tin = tin, yin = yin, xin = cbind(xin, zin), wi = NULL,
                                  bw.seq = bw.seq2, 
                                  nRegGrid = 101, kernel='epan',
                                  npoly = 1, nder = 0, kFolds = 5) )
  lmdat <- data.frame(cbind(yin, xin, zin))
  names(lmdat)[1] <- 'Y'
  res_lm <- try( lm(Y~., data = lmdat) )
  res_lm$mse <- mean(res_lm$residuals^2)
  
  res_seqMS <- try( seqMS(dat = cbind(yin, tin[,1], 1, zin, xin)) )
  
  Res_SVCoef_ex2_n2000[[rid]] = res_SVCoef
  Res_VCoef_ex2_n2000[[rid]] = res_VCoef
  Res_lm_ex2_n2000[[rid]] = res_lm
  Res_seqMS_ex2_n2000[[rid]] = res_seqMS
  
  
  ######
  
  dat.te <- makedata_fun_2d(n = n/5, m = model$m, alpfun = model$alpfun, beta = model$beta, 
                            rdist = rtrunc_gaussian_2d, t_corr = T,
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
  
  
  ###
  
  resCP = res_SVCoef
  
  resCP$tin.te = tin.te
  resCP$yin.te = yin.te
  resCP$xin.te = xin.te
  resCP$zin.te = zin.te
  
  library('gstat')
  library('sp')
  cord_idw = resCP$obsGrid
  alp_idw = resCP$alpha_hat
  
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
    
    # pdat = as.data.frame(cbind(tin.te, new_alp$var1.pred))
    # colnames(pdat) = c('t1', 't2', 'alpha')
    # lattice::cloud(alpha~t1*t2, data = pdat, zlab=list(rot=90), scales = list(arrows=F))
  }
  
  colnames(alp_plot) = paste0('alpha', 1:dim(alp_idw)[2]-1)
  resCP$alpha_hat.te = alp_plot
  
  alpha_te.hat <- alp_plot
  yhat.te <- c(zin.te %*% resCP$beta_hat + rowSums(alpha_te.hat * cbind(1, xin.te)))
  
  yid = which(!is.na(yhat.te))
  
  mse.te <- mean((yhat.te[yid]-yin.te[yid])^2)
  resCP$mse.te = mse.te
  
  Res_SVCoef_ex2_n2000[[rid]] = resCP
  
  ###
  
  resCP = res_VCoef
  
  resCP$tin.te = tin.te
  resCP$yin.te = yin.te
  resCP$xin.te = cbind(xin.te, zin.te)
  
  library('gstat')
  library('sp')
  cord_idw = resCP$obsGrid
  alp_idw = resCP$alpha_hat
  
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
    
    # pdat = as.data.frame(cbind(tin.te, new_alp$var1.pred))
    # colnames(pdat) = c('t1', 't2', 'alpha')
    # lattice::cloud(alpha~t1*t2, data = pdat, zlab=list(rot=90), scales = list(arrows=F))
  }
  
  colnames(alp_plot) = paste0('alpha', 1:dim(alp_idw)[2]-1)
  resCP$alpha_hat.te = alp_plot
  
  alpha_te.hat <- alp_plot
  yhat.te <- rowSums(alpha_te.hat * cbind(1, resCP$xin.te))
  
  yid = which(!is.na(yhat.te))
  
  mse.te <- mean((yhat.te[yid]-yin.te[yid])^2)
  resCP$mse.te = mse.te
  
  Res_VCoef_ex2_n2000[[rid]] = resCP
  
  ###
  
  res_seqMS.te <- try( seqMS.predict(dat = cbind(yin.te, tin.te[,1], 1, zin.te, xin.te),
                                     thresholds = res_seqMS$thresholds,
                                     beta = res_seqMS$beta) )
  
  res_seqMS = c(res_seqMS, res_seqMS.te)
  res_seqMS$obsGrid = tin
  
  Res_seqMS_ex2_n2000[[rid]] = res_seqMS
}


save(Res_SVCoef_ex2_n500, Res_VCoef_ex2_n500, Res_lm_ex2_n500, Res_seqMS_ex2_n500,
     Res_SVCoef_ex2_n1000, Res_VCoef_ex2_n1000, Res_lm_ex2_n1000, Res_seqMS_ex2_n1000,
     Res_SVCoef_ex2_n2000, Res_VCoef_ex2_n2000, Res_lm_ex2_n2000, Res_seqMS_ex2_n2000,
     file = 'Res_Ex2_compare.RData')














############### Add Testing Results



Res = Res_ex2_n2000 

for (i in 1:length(Res)) {
  
  resCP = Res[[i]]
  
  dat.te <- makedata_fun(n = dim(resCP$obsGrid)[1]/5, m = model$m, alpfun = model$alpfun, beta = model$beta, 
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

Res_ex2_n2000 = Res




