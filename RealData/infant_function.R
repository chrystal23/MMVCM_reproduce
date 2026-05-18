

func_path = '/path/to/your/code/' # './RealData'

source(paste0(func_path, 'simulation_function.R'))
source(paste0(func_path, 'compare.R'))

Rcpp::sourceCpp(paste0(func_path_int, 'CPPlwls2d_s1_real.cpp'))
Rcpp::sourceCpp(paste0(func_path_int, 'CPPlwls2d_s2_real.cpp'))
Rcpp::sourceCpp(paste0(func_path_int, 'CPPlwls2d_s2_LR_real.cpp'))



## main function


VCMCP_2d_real <- function(tin, yin, xin, zin, wi = NULL,
                          zeta.seq = NULL, cutoff = max, alpha = 0.05, 
                          bw.seq1 = NULL, bw.seq2 = NULL, bw.seq22 = NULL,
                          nRegGrid = 101, NbGrid = 101, kernel='epan',
                          npoly = 1, nder = 0, kFolds = 5, refined = FALSE,
                          hkappa = 2, method_opt=c('average','profile'),
                          h2_sepcv = F){
  
  # firsttsbreaks <- Sys.time()
  # Check the data validity for further analysis
  # CheckData(Ly,Lt)
  
  # inputData <- HandleNumericsAndNAN(Ly,Lt);
  # Ly <-  inputData$Ly;
  # Lt <-  inputData$Lt;
  
  n = dim(xin)[1]
  p = dim(xin)[2]
  q = dim(zin)[2]
  m = dim(tin)[2]
  
  method_opt <- match.arg(method_opt)
  
  # weight
  if (is.null(wi)) {
    win = rep(1, n)
    win = win/sum(win)
  } else {
    if (length(wi) != n){
      stop("The length of provided weight wi does not match data dimension.")
    }
    if (sum(wi>=0) < n){
      stop("Weight wi needs to be a vector of nonnegative numbers.")
    }
    wi = wi[order(tin[,1])]
    win = wi/sum(wi)
  }
  yin = yin[order(tin[,1])]
  xin = xin[order(tin[,1]), , drop=F]
  zin = zin[order(tin[,1]), , drop=F]
  tin = tin[order(tin[,1]), , drop=F]
  
  # Set the options structure members that are still NULL
  # optns <- list(methodBwMu = 'CV',
  #               kFoldMuCov = kFolds,
  #               nRegGrid = nRegGrid,
  #               error = TRUE,
  #               kernel= kernel,
  #               verbose=FALSE)
  # optns <- SetOptions(Ly, Lt, optns);
  # optns$useBW1SE = FALSE;
  # # optns$dataType = 'Sparse'
  # 
  # ### Check the options validity for the PCA function.
  # numOfCurves = length(Ly);
  # CheckOptions(Lt, optns, numOfCurves)
  
  
  
  ##### Stage 1 Estimation for the linear part
  
  ## k-fold cross validation to choose h_1
  
  cv1 = array(Inf, dim = c(length(bw.seq1), length(bw.seq1), kFolds));
  
  theFolds = SimpleFolds(1:n, kFolds)
  
  for (j in 1:length(bw.seq1)) {
    
    for (k in 1:length(bw.seq1)) {
      
      for (i in 1:kFolds) {
        
        ttest <- tin[theFolds[[i]], , drop=F]
        ytest <- yin[theFolds[[i]]]
        xtest <- xin[theFolds[[i]], , drop=F]
        ztest <- zin[theFolds[[i]], , drop=F]
        ttrain <- tin[-theFolds[[i]], , drop=F]
        ytrain <- yin[-theFolds[[i]]]
        xtrain <- xin[-theFolds[[i]], , drop=F]
        ztrain <- zin[-theFolds[[i]], , drop=F]
        wtrain <- win[-theFolds[[i]]]
        
        ## method 1
        if (method_opt == 'average'){
          coef_cv1 = tryCatch(
            CPPlwls2d_s1_real(bw = c(bw.seq1[j], bw.seq1[k]), kernel_type = kernel, win = wtrain,
                              tin = ttrain, yin = ytrain, xin = xtrain, zin = ztrain, tout = ttest, npoly = npoly),
            error=function(err) {
              warning('Invalid bandwidth during stage 1 CV. Try enlarging the window size.')
              return(Inf)
            })
          nan_rate <- length(which(is.nan(coef_cv1[,1]))) / length(ytest)
          # if (nan_rate > 0.2) cat('Warning: Bandwidth: h1 =', bw.seq1[j], bw.seq1[k], 'too small for stage 1 CV.')
          tid <- which(!is.nan(coef_cv1[,1]))
          beta_cv1 = colMeans(coef_cv1[tid, -c(1:(p+1)), drop=F])
          alp_cv1 = coef_cv1[tid, 1:(p+1), drop=F]
        }
        
        # ## method 2
        # if (method_opt == 'profile'){
        #   beta_cv1 = tryCatch(
        #     CPPlwls1d_s1_v2(bw = bw.seq1[j], kernel_type = kernel, win = wtrain,
        #                     tin = ttrain, yin = ytrain, xin = xtrain, zin = ztrain, npoly = npoly),
        #     error=function(err) {
        #       warning('Invalid bandwidth during stage 1 CV. Try enlarging the window size.')
        #       return(Inf)
        #     })
        #   alp_cv1 = CPPlwls1d_s2(bw = bw.seq1[j], kernel_type = kernel, win = wtrain,
        #                          tin = ttrain, yin = c(ytrain-ztrain %*% beta_cv1),
        #                          xin = xtrain, tout = ttest, npoly = npoly)
        # }
        
        yhat_cv1 = rowMeans(cbind(1, xtest[tid, , drop=F]) * alp_cv1) + ztest[tid, , drop=F] %*% beta_cv1
        
        cv1[j,k,i] = sum((ytest[tid] - c(yhat_cv1))^2)
        if (nan_rate >= 0.5) cv1[j,k,i] = Inf
        # cv[j,k,i] = trapzRcpp(xout, (obs - muout)^2)
        # print(cv)
        if(is.na(cv1[j,k,i]) || is.nan(cv1[j,k,i])){
          cv1[j,k,i] = Inf
        }
      }
    }
  }
  
  
  if(min(cv1) == Inf){
    stop("All bandwidths resulted in infinite CV costs. (Stage 1)")
  }
  
  cvMean = apply(cv1, c(1,2), mean)
  cvMeanid = which(cvMean == min(cvMean), arr.ind=TRUE)
  if (length(dim(cvMeanid))>1 && dim(cvMeanid)[1]>1) cvMeanid = cvMeanid[ceiling(dim(cvMeanid)[1]/2),]
  bopt1 = bw.seq1[cvMeanid];
  names(bopt1)[1:2] = c('h_1')
  
  
  ## use the chosen bandwidth h_1 to estimate beta
  
  h_1= unname(bopt1)
  
  ## method 1
  if (method_opt == 'average'){
    beta_est = CPPlwls2d_s1_real(bw = h_1, kernel_type = kernel, win = win,
                                 tin = tin, yin = yin, xin = xin, zin = zin,
                                 tout = tin, npoly = npoly)[,-c(1:(p+1)), drop=F]
    tid <- which(!is.nan(beta_est[,1]))
    beta_hat = colMeans(beta_est[tid, , drop=F])
  }
  
  # ## method 2
  # if (method_opt == 'profile'){
  #   beta_hat = CPPlwls1d_s1_v2(bw = h_1, kernel_type = kernel, win = win,
  #                              tin = tin, yin = yin, xin = xin, zin = zin,
  #                              npoly = npoly)
  # }
  
  yin_s2 = c(yin - zin %*% beta_hat)
  
  cat('\n', 'Stage 1 completed! \n')
  cat('Select h_1 = ', h_1, ';', 'beta_hat:', beta_hat, '\n')
  
  
  
  ##### Stage 2 Estimation for the nonparametric part
  
  # Initialize: k-fold cross validation to choose h_2
  h2.seq = bw.seq22
  cv2 = array(Inf, dim = c(length(h2.seq), length(h2.seq), kFolds));
  theFolds = SimpleFolds(1:n, kFolds)
  
  for (j in 1:length(h2.seq)) {
    
    for (k in 1:length(h2.seq)) {
      
      for (i in 1:kFolds) {
        
        ttest <- tin[theFolds[[i]], , drop=F]
        ytest <- yin_s2[theFolds[[i]]]
        xtest <- xin[theFolds[[i]], , drop=F]
        ttrain <- tin[-theFolds[[i]], , drop=F]
        ytrain <- yin_s2[-theFolds[[i]]]
        xtrain <- xin[-theFolds[[i]], , drop=F]
        wtrain <- win[-theFolds[[i]]]
        
        coef_cv2 = tryCatch(
          CPPlwls2d_s2_real(bw = c(h2.seq[j], h2.seq[k]), 
                            kernel_type = kernel, win = wtrain,
                            tin = ttrain, yin = ytrain, xin = xtrain, 
                            tout = ttest, npoly = npoly),
          error=function(err) {
            warning('Invalid bandwidth during stage 2 CV for h_2. Try enlarging the window size.')
            return(Inf)
          })
        nan_rate <- length(which(is.nan(coef_cv2[,1]))) / length(ytest)
        if (nan_rate > 0.2) cat('Warning: Bandwidth: h2 =', h2.seq[j], h2.seq[k], 'too small for stage 2 CV.')
        tid <- which(!is.nan(coef_cv2[,1]))
        alp_cv2 = coef_cv2[tid, 1:(p+1), drop=F]
        
        yhat_cv2 = rowMeans(cbind(1, xtest[tid, , drop=F]) * alp_cv2)
        
        cv2[j,k,i] = sum((ytest[tid] - c(yhat_cv2))^2)
        # print(cv2)
        if(is.na(cv2[j,k,i]) || is.nan(cv2[j,k,i]) || nan_rate>0.05){
          cv2[j,k,i] = Inf;
        }
      }
    }
  }
  
  if(min(cv2) == Inf){
    stop("All bandwidths resulted in infinite CV costs. (Stage 2, h_2)")
  }
  
  cvMean = apply(cv2, c(1,2), mean)
  cvMeanid = which(cvMean == min(cvMean), arr.ind=TRUE)
  if (length(dim(cvMeanid))>1 && dim(cvMeanid)[1]>1) cvMeanid = cvMeanid[floor(dim(cvMeanid)[1]/2),]
  bopt2 = h2.seq[cvMeanid];
  names(bopt2)[1:2] = c('h_2')
  
  h2_init = unname(bopt2)
  
  cat('Initial h_2 = ', h2_init, '\n')
  
  ###
  # firsttsCVmu <- Sys.time()
  ## mean function
  if(length(bw.seq2) == 1) {
    h_tau = rep(bw.seq2, m)
    h_d = rep(2 * bw.seq2, m)
    h_2 = h_tau+(h_d-h_tau)/3
    zeta = zetaFun(bw = bw.seq2, alpha = alpha, Lt = Lt, Ly = Ly, optns = optns,
                   cutoff = max)
    zeta = rep(zeta, m)
  } else if(!h2_sepcv) {
    
    ## select by CV procedure
    cat('Tuning parameter selection (CV: h_tau, zeta, h_d, h_2) for Stage 2 ......\n')
    sink(file = 'output.txt', append = T, split = F)
    h_tau = h_d = h_2 = zeta = c()
    for (ell in 1:m) {
      tunings = CVbandwidth_2d_real(bw.seq2 = bw.seq2, bw.seq22 = bw.seq22,
                                    zeta.seq = zeta.seq, win = win,
                                    tin = tin, yin = yin_s2, xin = xin,
                                    single_ell = ell, h2_init = h2_init,
                                    npoly = npoly, nder = nder, kernel = kernel,
                                    NbGrid = NbGrid, nRegGrid = nRegGrid,
                                    refined = refined, kFolds = kFolds,
                                    alpha = alpha, cutoff = cutoff,
                                    hkappa = hkappa)
      h_tau = c(h_tau, unname(tunings$bopt['h_tau']))
      h_d = c(h_d, unname(tunings$bopt['h_d']))
      zeta = c(zeta, tunings$zeta)
      h_2 = c(h_2, unname(tunings$bopt['h_2']))
    }
    
    closeAllConnections()
    sink(file = 'output.txt', append = T, split = T)
    
  } else {
    
    ## select by CV procedure
    cat('Tuning parameter selection (CV: h_tau, zeta, h_d) for Stage 2 ......\n')
    sink(file = 'output.txt', append = T, split = F)
    h_tau = h_d = h_2 = zeta = c()
    for (ell in 1:m) {
      tunings = CVbandwidth_2d_noh2_real(bw.seq2 = bw.seq2, bw.seq22 = bw.seq22, 
                                         zeta.seq = zeta.seq, win = win,
                                         tin = tin, yin = yin_s2, xin = xin,
                                         single_ell = ell, h2_init = h2_init,
                                         npoly = npoly, nder = nder, kernel = kernel,
                                         NbGrid = NbGrid, nRegGrid = nRegGrid,
                                         refined = refined, kFolds = kFolds,
                                         alpha = alpha, cutoff = cutoff,
                                         hkappa = hkappa)
      h_tau = c(h_tau, unname(tunings$bopt['h_tau']))
      h_d = c(h_d, unname(tunings$bopt['h_d']))
      zeta = c(zeta, tunings$zeta)
    }
    
    
    # k-fold cross validation to choose h_2
    h2.seq = sapply(1:m, function(ell) seq(from = h_tau[ell]+(h_d[ell]-h_tau[ell])/10, 
                                           to = h_d[ell]-(h_d[ell]-h_tau[ell])/3*2,
                                           length.out = length(bw.seq22)))  
    cv2 = array(Inf, dim = c(length(bw.seq22), length(bw.seq22), kFolds));
    theFolds = SimpleFolds(1:n, kFolds)
    
    for (j in 1:length(bw.seq22)){
      
      for (k in 1:length(bw.seq22)){
        
        for (i in 1:kFolds) {
          
          ttest <- tin[theFolds[[i]], , drop=F]
          ytest <- yin_s2[theFolds[[i]]]
          xtest <- xin[theFolds[[i]], , drop=F]
          ttrain <- tin[-theFolds[[i]], , drop=F]
          ytrain <- yin_s2[-theFolds[[i]]]
          xtrain <- xin[-theFolds[[i]], , drop=F]
          wtrain <- win[-theFolds[[i]]]
          
          muout = tryCatch(
            CoefJump_2d_real(tin = ttrain, yin = ytrain, xin = xtrain, win = wtrain,
                             tout = ttest, xout = xtest,
                             single_ell = NULL,
                             h_tau = h_tau, h_d = h_d, zeta = zeta, 
                             h_2 = c(h2.seq[j,1], h2.seq[k,2]),
                             refined = refined, npoly=npoly, nder= nder,
                             kernel = kernel,  NbGrid = NbGrid,
                             hkappa = hkappa, silent = T)$muout,
            error=function(err) {
              warning('Invalid bandwidth during stage 2 CV for h_2. Try enlarging the window size. h_2=', bw.seq22[j], '\n')
              return(Inf)
            })
          nan_rate <- length(which(is.nan(muout))) / length(ytest)
          if (nan_rate > 0.2) cat('Warning: Bandwidth: h_2=', c(h2.seq[j,1], h2.seq[k,2]), 
                                  'too small for varying coef CV.')
          tid <- which(!is.nan(muout))
          
          cv2[j,k,i] = sum((ytest[tid] - muout[tid])^2)
          # print(cv2)
          if(is.na(cv2[j,k,i]) || is.nan(cv2[j,k,i]) || nan_rate>0.2){
            cv2[j,k,i] = Inf;
          }
          
        }
      }
    }
    
    if(min(cv2) == Inf){
      stop("All bandwidths resulted in infinite CV costs. (Stage 2, h_2)")
    }
    
    cvMean = apply(cv2, c(1,2), mean)
    cvMeanid = which(cvMean == min(cvMean), arr.ind=TRUE)
    if (length(dim(cvMeanid))>1 && dim(cvMeanid)[1]>1) cvMeanid = cvMeanid[floor(dim(cvMeanid)[1]/2),]
    bopt2 = h2.seq[cbind(c(cvMeanid), 1:m)];
    names(bopt2)[1:2] = c('h_2')
    
    h_2 = unname(bopt2)
    
    
    closeAllConnections()
    sink(file = 'output.txt', append = T, split = T)
    
  }
  
  
  
  ##
  cat(' ', '\n')
  cat('Final selected: ', 'h_tau = ', h_tau, ';', 'h_d = ', h_d, 'h_2 = ', h_2, '\n')
  cat('and selected: ', 'zeta = ', zeta, '\n')
  
  # h_tau = 0.05
  # h_2 = 0.07
  # h_d = 0.17
  # zeta = 0.5
  
  ###
  obsGrid = tin
  ## cut in the region [h, 1-h]
  regGrid = sapply(1:m, function(ell){
    seq( max(min(obsGrid[,ell]), h_tau[ell]), 
         min(max(obsGrid[,ell]), 1 - h_tau[ell]), 
         length.out = nRegGrid)
  })
  regGrid = as.matrix(expand.grid(as.data.frame(regGrid)))
  
  
  # Get the mean function using the bandwidth estimated above:
  smcObj = CoefJump_2d_real(tin = tin, yin = yin_s2, xin = xin, win = win, 
                            tout = regGrid, xout = NULL,
                            single_ell = NULL,
                            h_tau = h_tau, h_d = h_d, zeta = zeta, h_2 = h_2,
                            refined = refined, npoly = npoly, nder = nder, 
                            kernel = kernel, NbGrid = NbGrid,
                            hkappa = hkappa)
  
  alpha_hat <- sapply(smcObj$alp_est, function(z) z$alp.hat)
  yhat <- c(zin %*% beta_hat + rowSums(alpha_hat * cbind(1, xin)))
  
  mse <- mean((yhat[!is.na(yhat)]-yin[!is.na(yhat)])^2)
  
  cat('MSE: ', mse, '\n')
  
  # mu = smcObj$mu
  # muWork = smcObj$muout
  # mu_jumpsize = smcObj$mu_jumpsize
  # rho = smcObj$rho_d
  # lasttsCVmu <- Sys.time()
  
  
  
  # firsttsCov <- Sys.time() #First time-stamp for calculation of the covariance
  # ## Covariance function and sigma2
  # scsObj = GetSmoothedCovarSurface(y = Ly, t = Lt, mu = mu, obsGrid = obsGrid,
  #                                  regGrid = regGrid, optns = optns)
  # sigma2 <- scsObj[['sigma2']]
  # 
  # # Get the results for the eigen-analysis
  # eigObj = GetEigenAnalysisResults(smoothCov = scsObj$smoothCov, regGrid, optns, muWork = muWork)
  # fittedCov = eigObj$fittedCov
  # 
  # lasttsCov <- Sys.time()
  
  
  #
  # timestamps = c(lasttsCVmu, lasttsCov, firsttsbreaks, firsttsCVmu, firsttsCov)
  
  timings = NULL
  
  # if(is.null(timestamps)) {
  #   timings = NULL;
  # } else {
  #   timestamps = c(Sys.time(), timestamps)
  #   timings = round(digits=3, timestamps[1:3]-timestamps[4:6]);
  #   names(timings) <- c('total', 'mu', 'cov')
  # }
  
  
  ##
  res <- list(beta_hat = beta_hat, alp_est = smcObj$alp_est, yhat = yhat, mse = mse,
              h_1 = h_1, h_tau = h_tau, h_d = h_d, h_2 = h_2,
              zeta = zeta, rho_d = unname(smcObj$rho_d),
              tin = tin, yin = yin, xin = xin, zin = zin, win = win,
              bw.seq1 = bw.seq1, bw.seq2 = bw.seq2, bw.seq22 = bw.seq22,
              obsGrid = obsGrid,
              workGrid = regGrid,
              timings = timings)
  
  return(res)
  
  
}






VCMCP_2d_real_bw <- function(tin, yin, xin, zin, 
                             h_1, h_tau, h_d, h_2, zeta,
                             wi = NULL,
                             cutoff = max, alpha = 0.05, 
                             nRegGrid = 101, NbGrid = 101, kernel='epan',
                             npoly = 1, nder = 0, kFolds = 5, refined = FALSE,
                             hkappa = 2, method_opt=c('average','profile'),
                             h2_sepcv = F){
  
  # firsttsbreaks <- Sys.time()
  # Check the data validity for further analysis
  # CheckData(Ly,Lt)
  
  # inputData <- HandleNumericsAndNAN(Ly,Lt);
  # Ly <-  inputData$Ly;
  # Lt <-  inputData$Lt;
  
  n = dim(xin)[1]
  p = dim(xin)[2]
  q = dim(zin)[2]
  m = dim(tin)[2]
  
  method_opt <- match.arg(method_opt)
  
  # weight
  if (is.null(wi)) {
    win = rep(1, n)
    win = win/sum(win)
  } else {
    if (length(wi) != n){
      stop("The length of provided weight wi does not match data dimension.")
    }
    if (sum(wi>=0) < n){
      stop("Weight wi needs to be a vector of nonnegative numbers.")
    }
    wi = wi[order(tin[,1])]
    win = wi/sum(wi)
  }
  yin = yin[order(tin[,1])]
  xin = xin[order(tin[,1]), , drop=F]
  zin = zin[order(tin[,1]), , drop=F]
  tin = tin[order(tin[,1]), , drop=F]
  
  
  h_1= h_1
  
  ## method 1
  if (method_opt == 'average'){
    beta_est = CPPlwls2d_s1_real(bw = h_1, kernel_type = kernel, win = win,
                                 tin = tin, yin = yin, xin = xin, zin = zin,
                                 tout = tin, npoly = npoly)[,-c(1:(p+1)), drop=F]
    tid <- which(!is.nan(beta_est[,1]))
    beta_hat = colMeans(beta_est[tid, , drop=F])
  }
  
  # ## method 2
  # if (method_opt == 'profile'){
  #   beta_hat = CPPlwls1d_s1_v2(bw = h_1, kernel_type = kernel, win = win,
  #                              tin = tin, yin = yin, xin = xin, zin = zin,
  #                              npoly = npoly)
  # }
  
  yin_s2 = c(yin - zin %*% beta_hat)
  
  cat('\n', 'Stage 1 completed! \n')
  cat('Select h_1 = ', h_1, ';', 'beta_hat:', beta_hat, '\n')
  
  
  
  ##### Stage 2 Estimation for the nonparametric part
  
  ##
  cat(' ', '\n')
  cat('Known: ', 'h_tau = ', h_tau, ';', 'h_d = ', h_d, 'h_2 = ', h_2, '\n')
  cat('Known: ', 'zeta = ', zeta, '\n')
  
  # h_tau = 0.05
  # h_2 = 0.07
  # h_d = 0.17
  # zeta = 0.5
  
  ###
  obsGrid = tin
  ## cut in the region [h, 1-h]
  regGrid = sapply(1:m, function(ell){
    seq( max(min(obsGrid[,ell]), h_tau[ell]), 
         min(max(obsGrid[,ell]), 1 - h_tau[ell]), 
         length.out = nRegGrid)
  })
  regGrid = as.matrix(expand.grid(as.data.frame(regGrid)))
  
  
  # Get the mean function using the bandwidth estimated above:
  smcObj = CoefJump_2d_real(tin = tin, yin = yin_s2, xin = xin, win = win, 
                            tout = regGrid, xout = NULL,
                            single_ell = NULL,
                            h_tau = h_tau, h_d = h_d, zeta = zeta, h_2 = h_2,
                            refined = refined, npoly = npoly, nder = nder, 
                            kernel = kernel, NbGrid = NbGrid,
                            hkappa = hkappa)
  
  alpha_hat <- sapply(smcObj$alp_est, function(z) z$alp.hat)
  yhat <- c(zin %*% beta_hat + rowSums(alpha_hat * cbind(1, xin)))
  
  mse <- mean((yhat[!is.na(yhat)]-yin[!is.na(yhat)])^2)
  
  cat('MSE: ', mse, '\n')
  
  # mu = smcObj$mu
  # muWork = smcObj$muout
  # mu_jumpsize = smcObj$mu_jumpsize
  # rho = smcObj$rho_d
  # lasttsCVmu <- Sys.time()
  
  
  
  # firsttsCov <- Sys.time() #First time-stamp for calculation of the covariance
  # ## Covariance function and sigma2
  # scsObj = GetSmoothedCovarSurface(y = Ly, t = Lt, mu = mu, obsGrid = obsGrid,
  #                                  regGrid = regGrid, optns = optns)
  # sigma2 <- scsObj[['sigma2']]
  # 
  # # Get the results for the eigen-analysis
  # eigObj = GetEigenAnalysisResults(smoothCov = scsObj$smoothCov, regGrid, optns, muWork = muWork)
  # fittedCov = eigObj$fittedCov
  # 
  # lasttsCov <- Sys.time()
  
  
  #
  # timestamps = c(lasttsCVmu, lasttsCov, firsttsbreaks, firsttsCVmu, firsttsCov)
  
  timings = NULL
  
  # if(is.null(timestamps)) {
  #   timings = NULL;
  # } else {
  #   timestamps = c(Sys.time(), timestamps)
  #   timings = round(digits=3, timestamps[1:3]-timestamps[4:6]);
  #   names(timings) <- c('total', 'mu', 'cov')
  # }
  
  
  ##
  res <- list(beta_hat = beta_hat, alp_est = smcObj$alp_est, yhat = yhat, mse = mse,
              h_1 = h_1, h_tau = h_tau, h_d = h_d, h_2 = h_2,
              zeta = zeta, rho_d = unname(smcObj$rho_d),
              tin = tin, yin = yin, xin = xin, zin = zin, win = win,
              obsGrid = obsGrid,
              workGrid = regGrid,
              timings = timings)
  
  return(res)
  
  
}




######################## Create k folds

SimpleFolds <- function(yy, k=10) {
  if (length(yy) > 1)
    allSamp <- sample(yy)
  else
    allSamp <- yy
  
  n <- length(yy)
  nEach <- n %/% k
  samp <- list()
  length(samp) <- k
  for (i in seq_along(samp)) {
    if (nEach > 0)
      samp[[i]] <- allSamp[1:nEach + (i - 1) * nEach]
    else
      samp[[i]] <- numeric(0)
  }
  restSamp <- allSamp[seq(nEach * k + 1, length(allSamp), length.out=length(allSamp) - nEach * k)]
  restInd <- sample(k, length(restSamp))
  for (i in seq_along(restInd)) {
    sampInd <- restInd[i]
    samp[[sampInd]] <- c(samp[[sampInd]], restSamp[i])
  }
  
  return(samp)
}










zetaFun_2d <- function(bw, alpha = 0.05, tin, yin, xin, win, 
                       npoly = 1, nder = 0, kernel = 'epan',
                       nRegGrid = 101, cutoff = max){
  
  n <- dim(xin)[1]
  
  ### epanechnikov kernel
  K_Epa = function(u) 0.75 * (1 - u^2) * (abs(u) <= 1);
  K_Epa_r = function(u, r) 0.75 * (1 - u^2) *u^r * (abs(u) <= 1);
  nu0 = integrate(K_Epa_r, 0, 1, r=0)$value
  nu1 = integrate(K_Epa_r, 0, 1, r=1)$value
  nu2 = integrate(K_Epa_r, 0, 1, r=2)$value
  K_Epa2 = function(u) (0.75 * (1 - u^2))^2 * (abs(u) <= 1);
  K_Epa2_r = function(u, r) (0.75 * (1 - u^2))^2 *u^r * (abs(u) <= 1);
  mu0 = integrate(K_Epa2_r, 0, 1, r=0)$value
  mu1 = integrate(K_Epa2_r, 0, 1, r=1)$value
  mu2 = integrate(K_Epa2_r, 0, 1, r=2)$value
  
  obsGrid = tin
  regGrid = seq( max(min(obsGrid[,1]), bw), min(max(obsGrid[,1]), 1-bw), length.out = nRegGrid)
  
  coef_hat <- CPPlwls2d_s2(bw = bw, 
                           tin = tin, yin = yin, xin = xin, win = win, 
                           tout = obsGrid, kernel_type = kernel, 
                           npoly = npoly, nder = nder)
  yhat <- c(rowSums(coef_hat * cbind(1, xin)))
  sigma2 <- mean((yhat-yin)^2)
  
  df = density(tin, kernel = 'rectangular')
  f.T <- cutoff( approx(df$x, df$y, xout=regGrid)$y )
  
  Gamma <- t(cbind(1,xin)) %*% cbind(1,xin) / n
  Lambdaj <- 2*sigma2/f.T * 
    (nu2^2*mu0 - 2*nu1*nu2*mu1 + nu1^2*mu2)/(nu0*nu2-nu1^2)^2 *
    cutoff(1/diag(Gamma))
  zeta <- stats::qnorm(1-alpha / 2) * sqrt(abs(Lambdaj)/(n*bw))  
  
  return(zeta)
  
}








CVbandwidth_2d_real <- function(bw.seq2 = NULL, bw.seq22 = NULL,
                                zeta.seq = NULL, win,
                                tin, yin, xin,
                                single_ell = NULL, h2_init,
                                npoly = 1, nder = 0, kernel = 'epan',
                                NbGrid, nRegGrid, 
                                refined = T, kFolds = 5, 
                                alpha = 0.05, cutoff = max,
                                hkappa = 2){
  
  m = dim(tin)[2]
  n = dim(xin)[1]
  p = dim(xin)[2]
  
  zlen = ifelse(is.null(zeta.seq), 6, length(zeta.seq))
  
  cv2 = array(Inf, dim = c(length(bw.seq2), zlen, length(bw.seq2), length(bw.seq22), kFolds));
  zeta.pool = matrix(0, nrow = length(bw.seq2), ncol = zlen)
  
  theFolds = SimpleFolds(1:n, kFolds)
  
  cat('\n\n#############',
      'CV procedure: \n')
  
  for (j in 1:length(bw.seq2)) {
    
    cat('h_tau = ', bw.seq2[j], '\n')
    
    if (is.null(zeta.seq)) {
      ## get cut-off value from the hypothesis test
      zeta.pool[j,] = zetaFun_2d(bw = bw.seq2[j], alpha = alpha,
                                 tin = tin, yin = yin, xin = xin, win = win,
                                 npoly = npoly, nder = nder, kernel = kernel,
                                 nRegGrid = nRegGrid, cutoff = cutoff)
      if (m == 1){
        ztt <- seq(-zeta.pool[j,1]/2, zeta.pool[j,1]/3, length.out=zlen)
      } else {
        ztt <- seq(zeta.pool[j,1]/3, zeta.pool[j,1], length.out=zlen)
      }
      zeta.pool[j,] = zeta.pool[j,] + ztt
      
    } else {
      zeta.pool[j,] = zeta.seq
    }
    
    for (zj in 1:zlen) {
      
      cat('--zeta = ', zeta.pool[j,zj], '\n')
      
      # for (k in 1:(nbw-1))
      for (k in 1:length(bw.seq2)) {
        
        if (2*bw.seq2[k] <= bw.seq2[j]) next
        cat('----h_d = ', 2*bw.seq2[k], '\n')
        
        for (k2 in 1:length(bw.seq22)) {
          
          if (2*bw.seq2[k] <= bw.seq22[k2]) next
          cat('------h_2 = ', bw.seq22[k2], ':  ')
          
          for (i in 1:kFolds) {
            
            ttest <- tin[theFolds[[i]], , drop=F]
            ytest <- yin[theFolds[[i]]]
            xtest <- xin[theFolds[[i]], , drop=F]
            ttrain <- tin[-theFolds[[i]], , drop=F]
            ytrain <- yin[-theFolds[[i]]]
            xtrain <- xin[-theFolds[[i]], , drop=F]
            wtrain <- win[-theFolds[[i]]]
            
            if (is.null(single_ell)) {
              muout = tryCatch(
                CoefJump_2d_real(tin = ttrain, yin = ytrain, xin = xtrain, win = wtrain, 
                                 tout = ttest, xout = xtest,
                                 single_ell = NULL,
                                 h_tau = rep(bw.seq2[j],m), h_d = rep(2*bw.seq2[k],m), 
                                 zeta = rep(zeta.pool[j,zj],m), h_2 = rep(bw.seq22[k2],m),
                                 refined = refined, npoly=npoly, nder= nder, 
                                 kernel = kernel,  NbGrid = NbGrid, 
                                 hkappa = hkappa, silent = F)$muout,
                error=function(err) {
                  warning('Invalid bandwidth during stage 2 CV. Try enlarging the window size. h_tau=', bw.seq2[j], 
                          ' h_d=', 2*bw.seq2[k], 'the', i, '-th fold \n')
                  return(Inf)
                })
            } else {
              ell = single_ell
              h_tau_t = h_d_t = h_2_t = h2_init
              zeta_t = rep(1, m)
              h_tau_t[ell] = bw.seq2[j]
              h_d_t[ell] = 2*bw.seq2[k]
              zeta_t[ell] = zeta.pool[j,zj]
              h_2_t[ell] = bw.seq22[k2]
              muout = tryCatch(
                CoefJump_2d_real(tin = ttrain, yin = ytrain, xin = xtrain, win = wtrain, 
                                 tout = ttest, xout = xtest,
                                 single_ell = single_ell,
                                 h_tau = h_tau_t, h_d = h_d_t, 
                                 zeta = zeta_t, h_2 = h_2_t,
                                 refined = refined, npoly=npoly, nder= nder, 
                                 kernel = kernel,  NbGrid = NbGrid, 
                                 hkappa = hkappa, silent = F)$muout,
                error=function(err) {
                  warning('Invalid bandwidth during stage 2 CV. Try enlarging the window size. h_tau=', bw.seq2[j], 
                          ' h_d=', 2*bw.seq2[k], 'the', i, '-th fold \n')
                  return(Inf)
                })
            }
            nan_rate <- length(which(is.nan(muout))) / length(ytest)
            if (nan_rate > 0.2) cat('Warning: Bandwidth: h_tau=', bw.seq2[j], 
                                    ' h_d=', 2*bw.seq2[k], 'too small for varying coef CV.')
            tid <- which(!is.nan(muout))
            
            cv2[j,zj,k,k2,i] = sum((ytest[tid] - muout[tid])^2)
            
            if(is.na(cv2[j,zj,k,k2,i]) || nan_rate > 0.05){
              cv2[j,zj,k,k2,i] = Inf;
            }
            
            cat('MSE: ', mean((ytest[tid] - muout[tid])^2), '\n')
            cat('==')
            
          }
          cat('\n')
        }
      }
    }
  }
  
  
  if(min(cv2) == Inf){
    stop("All bandwidths resulted in infinite CV costs. (Stage 2)")
  } else cat("CV completed!!! (Stage 2)")
  
  cvMean = apply(cv2, c(1,2,3,4), mean)
  cvMeanid = which(cvMean == min(cvMean), arr.ind=TRUE)
  if (length(dim(cvMeanid))>1 && dim(cvMeanid)[1]>1) cvMeanid = cvMeanid[ceiling(dim(cvMeanid)[1]/2),] # cvMeanid[dim(cvMeanid)[1],]
  bopt2 = bw.seq2[cvMeanid[c(1,3,4)]];
  names(bopt2) = c('h_tau', 'h_d', 'h_2')
  bopt2['h_d'] = 2 * bopt2['h_d']
  bopt2['h_2'] = bw.seq22[cvMeanid[4]]
  #zeta = zeta.pool[cvMeanid[1]]
  zeta = zeta.pool[cvMeanid[1],cvMeanid[2]]
  
  boptList <- list('bopt' = bopt2, 'zeta' = zeta, 'cvMean' = cvMean)
  
  return(boptList)
  
}



CVbandwidth_2d_noh2_real <- function(bw.seq2 = NULL, 
                                     zeta.seq = NULL, win,
                                     tin, yin, xin,
                                     single_ell = NULL, h2_init,
                                     npoly = 1, nder = 0, kernel = 'epan',
                                     NbGrid, nRegGrid, 
                                     refined = T, kFolds = 5, 
                                     alpha = 0.05, cutoff = max,
                                     hkappa = 2){
  
  m = dim(tin)[2]
  n = dim(xin)[1]
  p = dim(xin)[2]
  
  zlen = ifelse(is.null(zeta.seq), 6, length(zeta.seq))
  
  cv2 = array(Inf, dim = c(length(bw.seq2), zlen, length(bw.seq2), kFolds));
  zeta.pool = matrix(0, nrow = length(bw.seq2), ncol = zlen)
  
  theFolds = SimpleFolds(1:n, kFolds)
  
  cat('\n\n#############',
      'CV procedure: \n')
  
  for (j in 1:length(bw.seq2)) {
    
    cat('h_tau = ', bw.seq2[j], '\n')
    
    if (is.null(zeta.seq)) {
      ## get cut-off value from the hypothesis test
      zeta.pool[j,] = zetaFun_2d(bw = bw.seq2[j], alpha = alpha,
                                 tin = tin, yin = yin, xin = xin, win = win,
                                 npoly = npoly, nder = nder, kernel = kernel,
                                 nRegGrid = nRegGrid, cutoff = cutoff)
      if (m == 1){
        ztt <- seq(-zeta.pool[j,1]/2, zeta.pool[j,1]/3, length.out=zlen)
      } else {
        ztt <- seq(zeta.pool[j,1]/3, zeta.pool[j,1], length.out=zlen)
      }
      zeta.pool[j,] = zeta.pool[j,] + ztt
      
    } else {
      zeta.pool[j,] = zeta.seq
    }
    
    for (zj in 1:zlen) {
      
      cat('--zeta = ', zeta.pool[j,zj], '\n')
      
      # for (k in 1:(nbw-1))
      for (k in 1:length(bw.seq2)) {
        
        if (2*bw.seq2[k] <= bw.seq2[j]) next
        cat('----h_d = ', 2*bw.seq2[k], ': ')
        
        for (i in 1:kFolds) {
          
          ttest <- tin[theFolds[[i]], , drop=F]
          ytest <- yin[theFolds[[i]]]
          xtest <- xin[theFolds[[i]], , drop=F]
          ttrain <- tin[-theFolds[[i]], , drop=F]
          ytrain <- yin[-theFolds[[i]]]
          xtrain <- xin[-theFolds[[i]], , drop=F]
          wtrain <- win[-theFolds[[i]]]
          
          if (is.null(single_ell)){
            muout = tryCatch(
              CoefJump_2d_real(tin = ttrain, yin = ytrain, xin = xtrain, win = wtrain, 
                               tout = ttest, xout = xtest,
                               h_tau = bw.seq2[j], h_d = 2*bw.seq2[k], zeta = zeta.pool[j,zj], 
                               h_2 = (bw.seq2[j] + 2*bw.seq2[k])/2,
                               refined = refined, npoly=npoly, nder= nder, 
                               kernel = kernel,  NbGrid = NbGrid, 
                               hkappa = hkappa, silent = F)$muout,
              error=function(err) {
                warning('Invalid bandwidth during stage 2 CV. Try enlarging the window size. h_tau=', bw.seq2[j], 
                        ' h_d=', bw.seq2[k], ' the ', i, '-th fold \n')
                return(Inf)
              })
          } else {
            ell = single_ell
            h_tau_t = h_d_t = h_2_t = h2_init
            zeta_t = rep(1, m)
            h_tau_t[ell] = bw.seq2[j]
            h_d_t[ell] = 2*bw.seq2[k]
            zeta_t[ell] = zeta.pool[j,zj]
            # h_2_t[ell] = (bw.seq2[j] + 2*bw.seq2[k])/2
            muout = tryCatch(
              CoefJump_2d_real(tin = ttrain, yin = ytrain, xin = xtrain, win = wtrain, 
                               tout = ttest, xout = xtest,
                               single_ell = single_ell,
                               h_tau = h_tau_t, h_d = h_d_t, 
                               zeta = zeta_t, h_2 = h_2_t,
                               refined = refined, npoly=npoly, nder= nder, 
                               kernel = kernel,  NbGrid = NbGrid, 
                               hkappa = hkappa, silent = F)$muout,
              error=function(err) {
                warning('Invalid bandwidth during stage 2 CV. Try enlarging the window size. h_tau=', bw.seq2[j], 
                        ' h_d=', 2*bw.seq2[k], 'the', i, '-th fold \n')
                return(Inf)
              })
          }
          nan_rate <- length(which(is.nan(muout))) / length(ytest)
          if (nan_rate > 0.2) cat('Warning: Bandwidth: h_tau=', bw.seq2[j], 
                                  ' h_d=', 2*bw.seq2[k], 'too small for varying coef CV.')
          tid <- which(!is.nan(muout))
          
          cv2[j,zj,k,i] = sum((ytest[tid] - muout[tid])^2)
          # cv2[j,k,i] = trapzRcpp(xout, (obs - muout)^2)
          # print(cv2)
          if(is.na(cv2[j,zj,k,i]) || nan_rate > 0.05){
            cv2[j,zj,k,i] = Inf;
          }
          
          cat('MSE: ', mean((ytest[tid] - muout[tid])^2), '\n')
          cat('==')
          
        }
        cat('\n')
      }
    }
  }
  
  
  if(min(cv2) == Inf){
    stop("All bandwidths resulted in infinite CV costs. (Stage 2)")
  } else cat("CV completed!!! (Stage 2)")
  
  cvMean = apply(cv2, c(1,2,3), mean)
  cvMeanid = which(cvMean == min(cvMean), arr.ind=TRUE)
  if (length(dim(cvMeanid))>1 && dim(cvMeanid)[1]>1) cvMeanid = cvMeanid[floor(dim(cvMeanid)[1]/2),] # cvMeanid[dim(cvMeanid)[1],]
  bopt2 = bw.seq2[cvMeanid[c(1,3)]];
  names(bopt2) = c('h_tau', 'h_d')
  bopt2['h_d'] = 2 * bopt2['h_d']
  #zeta = zeta.pool[cvMeanid[1]]
  zeta = zeta.pool[cvMeanid[1],cvMeanid[2]]
  
  boptList <- list('bopt' = bopt2, 'zeta' = zeta, 'cvMean' = cvMean)
  
  return(boptList)
  
}





