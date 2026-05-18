

func_path_int = '/path/to/your/code/'

Rcpp::sourceCpp(paste0(func_path_int, 'CPPlwls2d_s1.cpp'))
Rcpp::sourceCpp(paste0(func_path_int, 'CPPlwls2d_s2.cpp'))
Rcpp::sourceCpp(paste0(func_path_int, 'CPPlwls2d_s2_LR.cpp'))





## generate simulation data



makedata_fun <- function(n=5000, m=2, rdist=runif, alpfun, beta, covrate=0.5,
                            t_corr = F,
                            noidist=rnorm, noisd=1, downsample = F,
                            draw = T) {
  
  if(n < 2){
    stop("Samples of size 1 are irrelevant.")
  }
  if(!is.function(rdist)){
    stop("'rdist' needs to be a function.")
  }
  if(sum(!sapply(alpfun, is.function))>0){
    stop("'alpfun' needs to be a list of functions.")
  }
  if(length(alpfun)<2){
    stop("Need to provide at least 2 functional coefficients (one intercept, one slope).")
  }
  if (!is.numeric(beta)) {
    stop("Constant coefficients 'beta' need to be numeric.")
  }
  if (!is.numeric(noisd) || noisd < 0) {
    stop("'sigma' needs to be a nonnegative number.")
  }
  
  
  if (downsample){
    if (m > 1) {
      NN = max(1000, ceiling(n^{1/m}))
      tin = matrix(rdist(NN*m), nrow = NN, ncol = m)
      tin = apply(tin, 2, sort)
      tin = expand.grid(as.data.frame(tin))
      tin = tin[sort(sample(NN^m, n, replace = F)), ]
      tin = tin[, c(m,1)]
      tin = as.matrix(tin)
    } else tin <- matrix(sort(rdist(n)), nrow = n, ncol = 1)
  } else {
    if (m > 1) {
      if  (!t_corr){
        tin <- matrix(NA, nrow = n, ncol = m)
        tin[,1] = sort(rdist(n))
        for (ell in 2:m) {
          tin[,ell] = rdist(n)
        }
      } else {
        tin = rdist(n)
      } 
    } else tin <- matrix(sort(rdist(n)), nrow = n, ncol = 1)
  }
  
  p <- length(alpfun) - 1
  q <- length(beta)
  
  xcov <- covrate^abs(row(diag(p))-col(diag(p)))
  xin <- mvtnorm::rmvnorm(n, sigma = xcov)
  
  zcov <- covrate^abs(row(diag(q))-col(diag(q)))
  zin <- mvtnorm::rmvnorm(n, sigma = zcov)
  
  varycoef <- sapply(alpfun, function(f){
    apply(tin, 1, f)
  })
  
  if (draw){
    if (m == 2) {
      pdat <- data.frame(t1=tin[,1], t2=tin[,2], 
                         alpha0=varycoef[,1], alpha1=varycoef[,2], alpha2=varycoef[,3])
      plot0 <- lattice::cloud(alpha0~t1*t2, data = pdat, zlab=list(rot=90), scales = list(arrows=F))
      plot1 <- lattice::cloud(alpha1~t1*t2, data = pdat, zlab=list(rot=90), scales = list(arrows=F))
      plot2 <- lattice::cloud(alpha2~t1*t2, data = pdat, zlab=list(rot=90), scales = list(arrows=F))
      a <- cowplot::plot_grid(plot0, plot1, plot2, nrow = 2)
      print(a)
      # library(ggplot2)
      # ggplot(pdat, aes(x=t1, y=t2)) + geom_point(aes(color=alpha)) +
      #   scale_color_gradient(low = 'lightyellow', high = 'darkblue')
      #scale_color_viridis_c(option='E',direction = 1)
    } else if (m == 1) {
      par(mfrow=c(1, dim(varycoef)[2]))
      for (j in 1:dim(varycoef)[2]){
        plot(c(tin), varycoef[,j])
      }
    }
  }
  
  y_m <- rowSums(varycoef * cbind(1, xin)) + 
    zin %*% beta
  y_sigma <- noidist(n, sd=noisd)
  yin <- y_m + y_sigma
  
  cat('signal var:', var(y_m), '\n')
  cat('noise var:', var(y_sigma), '\n')
  
  
  return(list(tin=tin, yin=yin, xin=xin, zin=zin, 
              varycoef=varycoef, beta=beta))
  
}








## main function

VCMCP <- function(tin, yin, xin, zin, wi = NULL,
                     zeta = NULL, cutoff = max, alpha = 0.05, 
                     bw.seq1 = NULL, bw.seq2 = NULL,
                     nRegGrid = 101, NbGrid = 101, kernel='epan',
                     RefiningStage = FALSE, max_iter = 10, tol_hausdorff = 1e-4, bw_amp = 1/3,
                     npoly = 1, nder = 0, kFolds = 5, 
                     refined = FALSE,
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
  
  cv1 = array(Inf, dim = c(length(bw.seq1), kFolds));
  
  theFolds = SimpleFolds(1:n, kFolds)
  
  for (j in 1:length(bw.seq1)) {
    
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
          CPPlwls2d_s1(bw = bw.seq1[j], kernel_type = kernel, win = wtrain,
                       tin = ttrain, yin = ytrain, xin = xtrain, zin = ztrain, tout = ttest, npoly = npoly),
          error=function(err) {
            warning('Invalid bandwidth during stage 1 CV. Try enlarging the window size.')
            return(Inf)
          })
        nan_rate <- length(which(is.nan(coef_cv1[,1]))) / length(ytest)
        if (nan_rate > 0.2) cat('Warning: Bandwidth: h1 =', bw.seq1[j], 'too small for stage 1 CV.')
        tid <- which(!is.nan(coef_cv1[,1]))
        beta_cv1 = colMeans(coef_cv1[tid, -c(1:(p+1)), drop=F])
        alp_cv1 = coef_cv1[tid, 1:(p+1), drop=F]
      }
      
      ## method 2
      if (method_opt == 'profile'){
        beta_cv1 = tryCatch(
          CPPlwls1d_s1_v2(bw = bw.seq1[j], kernel_type = kernel, win = wtrain,
                          tin = ttrain, yin = ytrain, xin = xtrain, zin = ztrain, npoly = npoly),
          error=function(err) {
            warning('Invalid bandwidth during stage 1 CV. Try enlarging the window size.')
            return(Inf)
          })
        alp_cv1 = CPPlwls1d_s2(bw = bw.seq1[j], kernel_type = kernel, win = wtrain,
                               tin = ttrain, yin = c(ytrain-ztrain %*% beta_cv1),
                               xin = xtrain, tout = ttest, npoly = npoly)
      }
      
      yhat_cv1 = rowMeans(cbind(1, xtest[tid, , drop=F]) * alp_cv1) + ztest[tid, , drop=F] %*% beta_cv1
      
      cv1[j,i] = sum((ytest[tid] - c(yhat_cv1))^2)
      if (nan_rate >= 0.5) cv1[j,i] = Inf  # only 0.95 for Perturbing h_1, 0.5 for other cases
      # cv[j,k,i] = trapzRcpp(xout, (obs - muout)^2) 
      # print(cv)
      if(is.na(cv1[j,i]) || is.nan(cv1[j,i])){
        cv1[j,i] = Inf
      }
      
    }
  }
  
  
  if(min(cv1) == Inf){
    stop("All bandwidths resulted in infinite CV costs. (Stage 1)")
  }
  
  cvMean = apply(cv1, c(1), mean)
  cvMeanid = which(cvMean == min(cvMean), arr.ind=TRUE)
  bopt1 = bw.seq1[cvMeanid];
  names(bopt1) = c('h_1')
  
  
  ## use the chosen bandwidth h_1 to estimate beta
  
  h_1= unname(bopt1)
  
  ## method 1
  if (method_opt == 'average'){
    beta_est = CPPlwls2d_s1(bw = h_1, kernel_type = kernel, win = win,
                            tin = tin, yin = yin, xin = xin, zin = zin,
                            tout = tin, npoly = npoly)[,-c(1:(p+1)), drop=F]
    tid <- which(!is.nan(beta_est[,1]))
    beta_hat = colMeans(beta_est[tid, , drop=F])
  }
  
  ## method 2
  if (method_opt == 'profile'){
    beta_hat = CPPlwls1d_s1_v2(bw = h_1, kernel_type = kernel, win = win,
                               tin = tin, yin = yin, xin = xin, zin = zin,
                               npoly = npoly)
  }
  
  yin_s2 = c(yin - zin %*% beta_hat)
  
  cat('\n', 'Stage 1 completed! \n')
  cat('Select h_1 = ', h_1, ';', 'beta_hat:', beta_hat, '\n')
  
  
  ##### Stage 2 Estimation for the nonparametric part
  
  ###
  # firsttsCVmu <- Sys.time()
  ## mean function
  if(length(bw.seq2) == 1) {
    h_tau = bw.seq2
    h_d = 2 * bw.seq2
    h_2 = h_tau
    zeta = ifelse(is.null(zeta.seq),  
                  zetaFun(bw = bw.seq2, alpha = alpha,
                             tin = tin, yin = yin, xin = xin, win = win,
                             npoly = npoly, nder = nder, kernel = kernel,
                             nRegGrid = nRegGrid, cutoff = cutoff), 
                  zeta.seq[1])
  } else if(!h2_sepcv) {
    
    ## select by CV procedure
    cat('Tuning parameter selection (CV: h_tau, zeta, h_d, h_2) for Stage 2 ......\n')
    # sink(file = 'output.txt', append = T, split = F)
    tunings = CVbandwidth(bw.seq2 = bw.seq2, zeta = zeta, win = win,
                             tin = tin, yin = yin_s2, xin = xin,
                             npoly = npoly, nder = nder, kernel = kernel,
                             NbGrid = NbGrid, nRegGrid = nRegGrid,
                             refined = refined, kFolds = kFolds,
                             alpha = alpha, cutoff = cutoff,
                             hkappa = hkappa)
    # closeAllConnections()
    # sink(file = 'output.txt', append = T, split = T)
    h_tau = unname(tunings$bopt['h_tau'])
    h_d = unname(tunings$bopt['h_d'])
    zeta = tunings$zeta
    h_2 = unname(tunings$bopt['h_2'])
    
  } else {
    
    ## select by CV procedure
    cat('Tuning parameter selection (CV: h_tau, zeta, h_d) for Stage 2 ......\n')
    sink(file = 'output.txt', append = T, split = F)
    tunings = CVbandwidth_noh2(bw.seq2 = bw.seq2, zeta = zeta, win = win,
                                  tin = tin, yin = yin_s2, xin = xin,
                                  npoly = npoly, nder = nder, kernel = kernel,
                                  NbGrid = NbGrid, nRegGrid = nRegGrid,
                                  refined = refined, kFolds = kFolds,
                                  alpha = alpha, cutoff = cutoff,
                                  hkappa = hkappa)
    h_tau = unname(tunings$bopt['h_tau'])
    h_d = unname(tunings$bopt['h_d'])
    zeta = tunings$zeta
    cat('h_tau = ', h_tau, ';', 'h_d = ', h_d, '\n')
    
    # k-fold cross validation to choose h_2
    h2.seq = seq(from = h_tau+(h_d-h_tau)/6, to = h_d-(h_d-h_tau)/2,
                 length.out = length(bw.seq2))
    cv2 = array(Inf, dim = c(length(h2.seq), kFolds));
    theFolds = SimpleFolds(1:n, kFolds)
    
    for (j in 1:length(h2.seq)) {
      
      for (i in 1:kFolds) {
        
        ttest <- tin[theFolds[[i]], , drop=F]
        ytest <- yin_s2[theFolds[[i]]]
        xtest <- xin[theFolds[[i]], , drop=F]
        ttrain <- tin[-theFolds[[i]], , drop=F]
        ytrain <- yin_s2[-theFolds[[i]]]
        xtrain <- xin[-theFolds[[i]], , drop=F]
        wtrain <- win[-theFolds[[i]]]
        
        muout = tryCatch(
          CoefJump(tin = ttrain, yin = ytrain, xin = xtrain, win = wtrain,
                      tout = ttest, xout = xtest,
                      h_tau = h_tau, h_d = h_d, zeta = zeta, h_2 = h2.seq[j],
                      refined = refined, npoly=npoly, nder= nder,
                      kernel = kernel,  NbGrid = NbGrid,
                      hkappa = hkappa, silent = T)$muout,
          error=function(err) {
            warning('Invalid bandwidth during stage 2 CV for h_2. Try enlarging the window size. h_2=', bw.seq2[j], '\n')
            return(Inf)
          })
        nan_rate <- length(which(is.nan(muout))) / length(ytest)
        if (nan_rate > 0.2) cat('Warning: Bandwidth: h_2=', h2.seq[j], 'too small for varying coef CV.')
        tid <- which(!is.nan(muout))
        
        cv2[j,i] = sum((ytest[tid] - muout[tid])^2)
        # print(cv2)
        if(is.na(cv2[j,i]) || is.nan(cv2[j,i]) || nan_rate>0.2){
          cv2[j,i] = Inf;
        }
        
      }
    }
    
    if(min(cv2) == Inf){
      stop("All bandwidths resulted in infinite CV costs. (Stage 2, h_2)")
    }
    
    cvMean = apply(cv2, c(1), mean)
    cvMeanid = which(cvMean == min(cvMean), arr.ind=TRUE)
    if (length(dim(cvMeanid))>1 && dim(cvMeanid)[1]>1) cvMeanid = cvMeanid[floor(dim(cvMeanid)[1]/2),]
    bopt2 = h2.seq[cvMeanid];
    names(bopt2) = c('h_2')
    
    h_2 = unname(bopt2)
    
    # closeAllConnections()
    # sink(file = 'output.txt', append = T, split = T)
    
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
  regGrid = apply(obsGrid, 2, function(x){
    seq( max(min(x), h_tau), min(max(x), 1- h_tau), length.out = nRegGrid)
  })
  regGrid = as.matrix(expand.grid(as.data.frame(regGrid)))
  
  
  # Get the mean function using the bandwidth estimated above:
  smcObj = CoefJump(tin = tin, yin = yin_s2, xin = xin, win = win, 
                       tout = regGrid, xout = NULL,
                       h_tau = h_tau, h_d = h_d, zeta = zeta, h_2 = h_2,
                       RefiningStage = RefiningStage, max_iter = max_iter, tol_hausdorff = tol_hausdorff, bw_amp = bw_amp,
                       npoly = npoly, nder = nder, 
                       refined = refined, 
                       kernel = kernel, NbGrid = NbGrid,
                       hkappa = hkappa)
  
  alpha_hat <- sapply(smcObj$alp_est, function(z) z$alp.hat)
  yhat <- c(zin %*% beta_hat + rowSums(alpha_hat * cbind(1, xin)))
  
  mse <- mean((yhat[!is.na(yhat)]-yin[!is.na(yhat)])^2)
  
  cat('MSE: ', mse, '\n')
  
  
  timings = NULL
  
  
  ##
  res <- list(beta_hat = beta_hat, alp_est = smcObj$alp_est, yhat = yhat, mse = mse,
              h_1 = h_1, h_tau = h_tau, h_d = h_d, h_2 = h_2,
              zeta = zeta, rho_d = unname(smcObj$rho_d),
              tin = tin, yin = yin, xin = xin, zin = zin, win = win,
              bw.seq1 = bw.seq1, bw.seq2 = bw.seq2,
              obsGrid = obsGrid,
              workGrid = regGrid,
              iter = smcObj$iter,
              timings = timings)
  
  return(res)
  
  
}







##### User-specified bandwidths

VCMCP_bw <- function(tin, yin, xin, zin, 
                        h_1, h_tau, h_d, h_2, zeta,
                        wi = NULL, 
                        cutoff = max, alpha = 0.05, 
                        nRegGrid = 101, NbGrid = 101, kernel='epan',
                        RefiningStage = FALSE, max_iter = 10, tol_hausdorff = 1e-4,
                        npoly = 1, nder = 0, kFolds = 5, 
                        refined = FALSE,
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
  
  
  
  ##### Stage 1 Estimation for the linear part
  
  h_1= h_1
  
  ## method 1
  if (method_opt == 'average'){
    beta_est = CPPlwls2d_s1(bw = h_1, kernel_type = kernel, win = win,
                            tin = tin, yin = yin, xin = xin, zin = zin,
                            tout = tin, npoly = npoly)[,-c(1:(p+1)), drop=F]
    tid <- which(!is.nan(beta_est[,1]))
    beta_hat = colMeans(beta_est[tid, , drop=F])
  }
  
  ## method 2
  if (method_opt == 'profile'){
    beta_hat = CPPlwls1d_s1_v2(bw = h_1, kernel_type = kernel, win = win,
                               tin = tin, yin = yin, xin = xin, zin = zin,
                               npoly = npoly)
  }
  
  yin_s2 = c(yin - zin %*% beta_hat)
  
  cat('\n', 'Stage 1 completed! \n')
  cat('Known h_1 = ', h_1, ';', 'beta_hat:', beta_hat, '\n')
  
  
  ##### Stage 2 Estimation for the nonparametric part
  
  ##
  cat(' ', '\n')
  cat('Known: ', 'h_tau = ', h_tau, ';', 'h_d = ', h_d, 'h_2 = ', h_2, '\n')
  cat('Known: ', 'zeta = ', zeta, '\n')
  
  
  ###
  obsGrid = tin
  ## cut in the region [h, 1-h]
  regGrid = apply(obsGrid, 2, function(x){
    seq( max(min(x), h_tau), min(max(x), 1- h_tau), length.out = nRegGrid)
  })
  regGrid = as.matrix(expand.grid(as.data.frame(regGrid)))
  
  
  # Get the mean function using the bandwidth estimated above:
  smcObj = CoefJump(tin = tin, yin = yin_s2, xin = xin, win = win, 
                       tout = regGrid, xout = NULL,
                       h_tau = h_tau, h_d = h_d, zeta = zeta, h_2 = h_2,
                       RefiningStage = RefiningStage, max_iter = max_iter, tol_hausdorff = tol_hausdorff, bw_amp = bw_amp,
                       npoly = npoly, nder = nder, 
                       refined = refined, 
                       kernel = kernel, NbGrid = NbGrid,
                       hkappa = hkappa)
  
  alpha_hat <- sapply(smcObj$alp_est, function(z) z$alp.hat)
  yhat <- c(zin %*% beta_hat + rowSums(alpha_hat * cbind(1, xin)))
  
  mse <- mean((yhat[!is.na(yhat)]-yin[!is.na(yhat)])^2)
  
  cat('MSE: ', mse, '\n')
  
 
  
  timings = NULL

  
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










zetaFun <- function(bw, alpha = 0.05, tin, yin, xin, win, 
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
  yid = which(!is.nan(yhat))
  sigma2 <- mean((yhat[yid]-yin[yid])^2)
  
  df = density(tin, kernel = 'rectangular')
  f.T <- cutoff( approx(df$x, df$y, xout=regGrid)$y )
  
  Gamma <- t(cbind(1,xin)) %*% cbind(1,xin) / n
  Lambdaj <- 2*sigma2/f.T * 
    (nu2^2*mu0 - 2*nu1*nu2*mu1 + nu1^2*mu2)/(nu0*nu2-nu1^2)^2 *
    cutoff(1/diag(Gamma))
  zeta <- stats::qnorm(1-alpha / 2) * sqrt(abs(Lambdaj)/(n*bw))  
  
  return(zeta)
  
}








CVbandwidth <- function(bw.seq2 = NULL, zeta = NULL, win = win,
                           tin = tin, yin = yin, xin = xin,
                           npoly = 1, nder = 0, kernel = 'epan',
                           NbGrid, nRegGrid, 
                           refined = T, kFolds = 5, 
                           alpha = 0.05, cutoff = max,
                           hkappa = 2){
  
  m = dim(tin)[2]
  n = dim(xin)[1]
  p = dim(xin)[2]
  
  zlen = ifelse(is.null(zeta), 6, length(zeta))
  
  cv2 = array(Inf, dim = c(length(bw.seq2), zlen, length(bw.seq2), length(bw.seq2), kFolds));
  zeta.seq = matrix(0, nrow = length(bw.seq2), ncol = zlen)
  
  theFolds = SimpleFolds(1:n, kFolds)
  
  # cat('\n\n#############',
  #     'CV procedure: \n')
  
  for (j in 1:length(bw.seq2)) {
    
    # cat('h_tau = ', bw.seq2[j], '\n')
    
    if (is.null(zeta)) {
      ## get cut-off value from the hypothesis test
      zeta.seq[j,] = zetaFun(bw = bw.seq2[j], alpha = alpha,
                                tin = tin, yin = yin, xin = xin, win = win,
                                npoly = npoly, nder = nder, kernel = kernel,
                                nRegGrid = nRegGrid, cutoff = cutoff)
      # if (is.nan(zeta.seq[j,1])) zeta.seq[j,1] = 0.4
      if (m == 1){
        ztt <- seq(-zeta.seq[j,1]/2, zeta.seq[j,1]/3, length.out=zlen)
      } else {
        ztt <- seq(zeta.seq[j,1]/3, zeta.seq[j,1], length.out=zlen)
      }
      zeta.seq[j,] = zeta.seq[j,] + ztt
      
    } else {
      zeta.seq[j,] = zeta
    }
    
    for (zj in 1:zlen) {
      
      # cat('--zeta = ', zeta.seq[j,zj], '\n')
      
      # for (k in 1:(nbw-1))
      for (k in 1:length(bw.seq2)) {
        
        if (2*bw.seq2[k] <= bw.seq2[j]) next
        # cat('----h_d = ', 2*bw.seq2[k], '\n')
        
        for (k2 in j:length(bw.seq2)) {
          
          if (2*bw.seq2[k] <= bw.seq2[k2]) next
          # cat('------h_2 = ', bw.seq2[k2], ':  ')
          
          for (i in 1:kFolds) {
            
            ttest <- tin[theFolds[[i]], , drop=F]
            ytest <- yin[theFolds[[i]]]
            xtest <- xin[theFolds[[i]], , drop=F]
            ttrain <- tin[-theFolds[[i]], , drop=F]
            ytrain <- yin[-theFolds[[i]]]
            xtrain <- xin[-theFolds[[i]], , drop=F]
            wtrain <- win[-theFolds[[i]]]
            
            muout = tryCatch(
              CoefJump(tin = ttrain, yin = ytrain, xin = xtrain, win = wtrain, 
                          tout = ttest, xout = xtest,
                          h_tau = bw.seq2[j], h_d = 2*bw.seq2[k], zeta = zeta.seq[j,zj], h_2 = bw.seq2[k2],
                          refined = refined, npoly=npoly, nder= nder, 
                          kernel = kernel,  NbGrid = NbGrid, 
                          hkappa = hkappa, silent = T)$muout,
              error=function(err) {
                warning('Invalid bandwidth during stage 2 CV. Try enlarging the window size. h_tau=', bw.seq2[j], 
                        ' h_d=', 2*bw.seq2[k], 'the', i, '-th fold \n')
                return(Inf)
              })
            nan_rate <- length(which(is.nan(muout))) / length(ytest)
            # if (nan_rate > 0.2) cat('Warning: Bandwidth: h_tau=', bw.seq2[j], 
            #                         ' h_d=', 2*bw.seq2[k], 'too small for varying coef CV.')
            tid <- which(!is.nan(muout))
            
            cv2[j,zj,k,k2,i] = sum((ytest[tid] - muout[tid])^2)
            # cv2[j,k,i] = trapzRcpp(xout, (obs - muout)^2)
            # print(cv2)
            if(is.na(cv2[j,zj,k,k2,i]) || nan_rate > 0.5){
              cv2[j,zj,k,k2,i] = Inf;
            }
            
            # cat('MSE: ', mean(ytest[tid] - muout[tid])^2, '\n')
            # cat('==')
            
          }
          # cat('\n')
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
  #zeta = zeta.seq[cvMeanid[1]]
  zeta = zeta.seq[cvMeanid[1],cvMeanid[2]]
  
  boptList <- list('bopt' = bopt2, 'zeta' = zeta, 'cvMean' = cvMean)
  
  return(boptList)
  
}



CVbandwidth_noh2 <- function(bw.seq2 = NULL, zeta = NULL, win = win,
                                tin = tin, yin = yin, xin = xin,
                                npoly = 1, nder = 0, kernel = 'epan',
                                NbGrid, nRegGrid, 
                                refined = T, kFolds = 5, 
                                alpha = 0.05, cutoff = max,
                                hkappa = 2){
  
  m = dim(tin)[2]
  n = dim(xin)[1]
  p = dim(xin)[2]
  
  zlen = ifelse(is.null(zeta), 6, length(zeta))
  
  cv2 = array(Inf, dim = c(length(bw.seq2), zlen, length(bw.seq2), kFolds));
  zeta.seq = matrix(0, nrow = length(bw.seq2), ncol = zlen)
  
  theFolds = SimpleFolds(1:n, kFolds)
  
  # cat('\n\n#############',
  #     'CV procedure: \n')
  
  for (j in 1:length(bw.seq2)) {
    
    # cat('h_tau = ', bw.seq2[j], '\n')
    
    if (is.null(zeta)) {
      ## get cut-off value from the hypothesis test
      zeta.seq[j,] = zetaFun(bw = bw.seq2[j], alpha = alpha,
                                tin = tin, yin = yin, xin = xin, win = win,
                                npoly = npoly, nder = nder, kernel = kernel,
                                nRegGrid = nRegGrid, cutoff = cutoff)
      if (m == 1){
        ztt <- seq(-zeta.seq[j,1]/2, zeta.seq[j,1]/3, length.out=zlen)
      } else {
        ztt <- seq(zeta.seq[j,1]/3, zeta.seq[j,1], length.out=zlen)
      }
      zeta.seq[j,] = zeta.seq[j,] + ztt
      
    } else {
      zeta.seq[j,] = zeta
    }
    
    for (zj in 1:zlen) {
      
      # cat('--zeta = ', zeta.seq[j,zj], '\n')
      
      # for (k in 1:(nbw-1))
      for (k in 1:length(bw.seq2)) {
        
        if (2*bw.seq2[k] <= bw.seq2[j]) next
        # cat('----h_d = ', 2*bw.seq2[k], ': ')
        
        for (i in 1:kFolds) {
          
          ttest <- tin[theFolds[[i]], , drop=F]
          ytest <- yin[theFolds[[i]]]
          xtest <- xin[theFolds[[i]], , drop=F]
          ttrain <- tin[-theFolds[[i]], , drop=F]
          ytrain <- yin[-theFolds[[i]]]
          xtrain <- xin[-theFolds[[i]], , drop=F]
          wtrain <- win[-theFolds[[i]]]
          
          muout = tryCatch(
            CoefJump(tin = ttrain, yin = ytrain, xin = xtrain, win = wtrain, 
                        tout = ttest, xout = xtest,
                        h_tau = bw.seq2[j], h_d = 2*bw.seq2[k], zeta = zeta.seq[j,zj], 
                        h_2 = (bw.seq2[j] + 2*bw.seq2[k])/2,
                        refined = refined, npoly=npoly, nder= nder, 
                        kernel = kernel,  NbGrid = NbGrid, 
                        hkappa = hkappa, silent = T)$muout,
            error=function(err) {
              warning('Invalid bandwidth during stage 2 CV. Try enlarging the window size. h_tau=', bw.seq2[j], 
                      ' h_d=', bw.seq2[k], ' the ', i, '-th fold \n')
              return(Inf)
            })
          nan_rate <- length(which(is.nan(muout))) / length(ytest)
          # if (nan_rate > 0.2) cat('Warning: Bandwidth: h_tau=', bw.seq2[j], 
          #                         ' h_d=', 2*bw.seq2[k], 'too small for varying coef CV.')
          tid <- which(!is.nan(muout))
          
          cv2[j,zj,k,i] = sum((ytest[tid] - muout[tid])^2)
          # cv2[j,k,i] = trapzRcpp(xout, (obs - muout)^2)
          # print(cv2)
          if(is.na(cv2[j,zj,k,i]) || nan_rate > 0.5){
            cv2[j,zj,k,i] = Inf;
          }
          
          # cat('MSE: ', mean(ytest[tid] - muout[tid])^2, '\n')
          # cat('==')
          
        }
        # cat('\n')
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
  #zeta = zeta.seq[cvMeanid[1]]
  zeta = zeta.seq[cvMeanid[1],cvMeanid[2]]
  
  boptList <- list('bopt' = bopt2, 'zeta' = zeta, 'cvMean' = cvMean)
  
  return(boptList)
  
}






### mean estimation with multiple breaks
CoefJump <- function(tin, yin, xin, win, 
                        tout, xout = NULL,
                        h_tau, h_d, zeta, h_2,
                        RefiningStage = F, max_iter = 10, tol_hausdorff = 1e-4, bw_amp = 1/3,
                        npoly = 1, nder = 0, 
                        kernel = 'epan', NbGrid = 101,
                        refined = FALSE,
                        hkappa = 2, silent = FALSE) {
  
  timings = NULL
  
  n = dim(xin)[1]
  p = dim(xin)[2]
  m = dim(tin)[2]
  
  yyin = yin
  
  res = vector(mode = 'list', length = p+1)
  names(res) = paste0('alpha', c(0:p))
  
  rho_d = NULL
  
  # ---------------------------------------------------------
  # PHASE 1: SEARCHING PHASE & CACHING
  # ---------------------------------------------------------
  cached_grids = list() # Store grids and differences for the refining phase
  
  if (!silent) cat('-----------------------------------------\n', 
                   'Searching Phase\n',
                   '---------------\n')
  
  for (ell in 1:m) {
    tin_dm = tin
    tin_dm[, c(1, ell)] = tin[, c(ell, 1)]
    tord = order(tin_dm[, 1])
    tin_dm = tin_dm[tord, ,drop=F]
    xin_dm = xin[tord, ,drop=F]
    yin_dm = yin[tord]
    win_dm = win[tord]
    
    # Generate basic grids for jump detect:
    obsGrid = tin_dm;
    if(is.null(NbGrid)){
      jumpGrid = obsGrid
    } else {
      jumpGrid = apply(obsGrid, 2, function(x){
        seq( max(min(x), h_tau), min(max(x), 1- h_tau), length.out = NbGrid)
      })
      D_h = jumpGrid[, 1]
      jumpGrid[, c(1, m)] = jumpGrid[, c(m, 1)]
      jumpGrid = as.matrix(expand.grid(as.data.frame(jumpGrid)))
      jumpGrid[, c(1, m)] = jumpGrid[, c(m, 1)] 
    }
    
    ## Local linear estimate based on one-sided kernel
    alp_est_lr = CPPlwls2d_s2_LR(bw = h_tau, kernel_type = kernel, win = win_dm,
                                 tin = tin_dm, yin = yin_dm, xin = xin_dm, tout = jumpGrid, npoly = npoly)
    alp_est_left = alp_est_lr[,1:(p+1), drop=F]
    alp_est_right = alp_est_lr[,-c(1:(p+1)), drop=F]
    
    # Cache the grid and the absolute differences for Phase 2
    cached_grids[[ell]] = list(
      jumpGrid = jumpGrid,
      D_h = D_h,
      diff_abs = abs(alp_est_right - alp_est_left),
      tin_dm = tin_dm, xin_dm = xin_dm, yin_dm = yin_dm, win_dm = win_dm
    )
    
    for (j in 1:(p+1)) {
      res[[j]][[ell]] = vector(mode = 'list', length = 3)
      names(res[[j]][[ell]]) = c('jumptime','jumpsize','jumpsize_h_tau')
      
      alpj_diff_abs = cached_grids[[ell]]$diff_abs[, j]
      alpj_diff_abs = sapply(1:NbGrid, function(ii) mean( alpj_diff_abs[jumpGrid[,1] == D_h[ii]] ))
      
      alpj_diff_time = D_h[order(-alpj_diff_abs)]
      alpj_diff_size = alpj_diff_abs[order(-alpj_diff_abs)]
      
      timeindex = which(abs(alpj_diff_size) > zeta)
      ll=length(timeindex)
      
      alpj_jumptime = NULL
      alpj_jumpsize_h_tau = NULL
      if (ll > 0) {
        alpj_diff_time=alpj_diff_time[timeindex]
        alpj_diff_size=alpj_diff_size[timeindex]
        
        alpj_jumptime=alpj_diff_time[1]
        alpj_jumpsize=alpj_diff_size[1]
        for (i in 2:ll){
          index=which(abs(alpj_diff_time-alpj_jumptime[i-1]) > hkappa * h_tau)
          if (length(index) > 0){ 
            alpj_diff_time = alpj_diff_time[index]
            alpj_diff_size = alpj_diff_size[index]
            alpj_jumptime=append(alpj_jumptime, alpj_diff_time[1])
            alpj_jumpsize=append(alpj_jumpsize, alpj_diff_size[1])
          } else {
            break
          }
        }
        alpj_jumpsize_h_tau = alpj_jumpsize[order(alpj_jumptime)]
        alpj_jumptime = sort(alpj_jumptime)
        
        # Compute jump sizes using h_d
        rho_d = 0.25*(h_tau^2 + sqrt(log(n)/(n*h_tau))*h_tau)
        names(rho_d) = c("rho")
        jumpset_l = c(alpj_jumptime - rho_d)
        jumpset_r = c(alpj_jumptime + rho_d)
        
        othergrid = apply(tin[, -c(1), drop=F], 2, function(x){
          seq( max(min(x), h_tau), min(max(x), 1- h_tau), length.out = NbGrid)
        }, simplify = F)
        
        jumpset_l_grid = as.matrix(expand.grid(c(othergrid, list(jumpset_l))))
        jumpset_r_grid = as.matrix(expand.grid(c(othergrid, list(jumpset_r))))
        if (m>1) {
          jumpset_l_grid = jumpset_l_grid[,c(m, 1:(m-1))]
          jumpset_r_grid = jumpset_r_grid[,c(m, 1:(m-1))]
        }
        
        alpj_est_lr_l = CPPlwls2d_s2_LR(bw = h_d, kernel_type = kernel, win = win,
                                        tin = cached_grids[[ell]]$tin_dm, yin = cached_grids[[ell]]$yin_dm, 
                                        xin = cached_grids[[ell]]$xin_dm, tout = jumpset_l_grid, npoly = npoly)
        alpj_est_lr_r = CPPlwls2d_s2_LR(bw = h_d, kernel_type = kernel, win = win,
                                        tin = cached_grids[[ell]]$tin_dm, yin = cached_grids[[ell]]$yin_dm, 
                                        xin = cached_grids[[ell]]$xin_dm, tout = jumpset_r_grid, npoly = npoly)
        
        alpj_jumpsize_h_d = alpj_est_lr_r[,p+1+j] - alpj_est_lr_l[,j]
        alpj_jumpsize_h_d = sapply(1:length(jumpset_l), 
                                   function(ii) mean( alpj_jumpsize_h_d[cached_grids[[ell]]$jumpGrid[,1] == cached_grids[[ell]]$D_h[ii]] ))
        
        ## refine the jump location
        if(refined){
          alpj_jumptime = alpj_jumptime[abs(alpj_jumpsize_h_d) > zeta]
          alpj_jumpsize_h_d = alpj_jumpsize_h_d[abs(alpj_jumpsize_h_d) > zeta]
        }
        
      }
      
      if(length(alpj_jumptime)==0){
        res[[j]][[ell]]$jumptime = NULL
        res[[j]][[ell]]$jumpsize = NULL
        res[[j]][[ell]]$jumpsize_h_tau = alpj_jumpsize_h_tau
        if (!silent) cat('  alpha_', j-1, 'no.', ell, 'dimension --- no change point after refine.', '\n')
      } else {
        res[[j]][[ell]]$jumptime = alpj_jumptime
        res[[j]][[ell]]$jumpsize = alpj_jumpsize_h_d
        res[[j]][[ell]]$jumpsize_h_tau = alpj_jumpsize_h_tau
        if (!silent) cat('  alpha_', j-1, 'no.', ell, 'dimension --- change point at:', 
                         res[[j]][[ell]]$jumptime, '  jump size:', res[[j]][[ell]]$jumpsize, '\n')
      }
      
    }
  }
  
  iter = NULL
  
  # ---------------------------------------------------------
  # PHASE 2: REFINING PHASE
  # ---------------------------------------------------------
  if (RefiningStage && m >= 2) {
    
    if (!silent) cat('-----------------------------------------\n', 
                     'Refining Phase\n',
                     '---------------\n')
    
    for (iter in 1:max_iter) {
      
      if (!silent) cat('Iteration No.', iter, '\n')
      
      max_dist = 0 
      
      # Store previous iteration's jumps for ALL j and ell to calculate Hausdorff distance later
      old_jumps <- lapply(1:(p+1), function(j) {
        lapply(1:m, function(ell) res[[j]][[ell]]$jumptime)
      })
      
      for (ell in 1:m) {
        
        for (j in 1:(p+1)) {
          
          current_grid = cached_grids[[ell]]$jumpGrid 
          current_diff_abs = cached_grids[[ell]]$diff_abs[, j] 
          D_h = cached_grids[[ell]]$D_h
          
          # Construct integration domain: exclude h_tau neighborhoods of jumps in OTHER dimensions
          valid_idx = rep(TRUE, nrow(current_grid))
          for (k in 1:m) {
            if (k == ell) next
            
            jumps_k = res[[j]][[k]]$jumptime
            if (!is.null(jumps_k)) {
              # Map original dimension k to the swapped columns in jumpGrid
              col_idx = ifelse(k == 1, ell, ifelse(k == ell, 1, k))
              
              for (tau in jumps_k) {
                valid_idx[abs(current_grid[, col_idx] - tau) <= h_tau * bw_amp] = FALSE
              }
            }
          }
          
          # Integrate over restricted domain
          refined_diff_abs = sapply(1:length(D_h), function(ii) {
            idx = which(current_grid[, 1] == D_h[ii] & valid_idx)
            if (length(idx) > 0) mean(current_diff_abs[idx]) else 0
          })
          
          # Find peaks again
          alpj_diff_time = D_h[order(-refined_diff_abs)]
          alpj_diff_size = refined_diff_abs[order(-refined_diff_abs)]
          
          timeindex = which(abs(alpj_diff_size) > zeta)
          ll_idx = length(timeindex)
          
          alpj_jumptime = NULL
          alpj_jumpsize_h_tau = NULL
          if (ll_idx > 0) {
            alpj_diff_time = alpj_diff_time[timeindex]
            alpj_diff_size = alpj_diff_size[timeindex]
            
            alpj_jumptime = alpj_diff_time[1]
            for (i in 2:ll_idx) {
              index = which(abs(alpj_diff_time - alpj_jumptime[i-1]) > hkappa * h_tau)
              if (length(index) > 0) {
                alpj_diff_time = alpj_diff_time[index]
                alpj_diff_size = alpj_diff_size[index]
                alpj_jumptime = append(alpj_jumptime, alpj_diff_time[1])
                alpj_jumpsize = append(alpj_jumpsize, alpj_diff_size[1])
              } else { 
                break 
              }
            }
            alpj_jumpsize_h_tau = alpj_jumpsize[order(alpj_jumptime)]
            alpj_jumptime = sort(alpj_jumptime)
            
            # Compute jump sizes using h_d
            rho_d = 0.25*(h_tau^2 + sqrt(log(n)/(n*h_tau))*h_tau)
            names(rho_d) = c("rho")
            jumpset_l = c(alpj_jumptime - rho_d)
            jumpset_r = c(alpj_jumptime + rho_d)
            
            othergrid = apply(tin[, -c(1), drop=F], 2, function(x){
              seq( max(min(x), h_tau), min(max(x), 1- h_tau), length.out = NbGrid)
            }, simplify = F)
            
            jumpset_l_grid = as.matrix(expand.grid(c(othergrid, list(jumpset_l))))
            jumpset_r_grid = as.matrix(expand.grid(c(othergrid, list(jumpset_r))))
            if (m>1) {
              jumpset_l_grid = jumpset_l_grid[,c(m, 1:(m-1))]
              jumpset_r_grid = jumpset_r_grid[,c(m, 1:(m-1))]
            }
            
            alpj_est_lr_l = CPPlwls2d_s2_LR(bw = h_d, kernel_type = kernel, win = win,
                                            tin = cached_grids[[ell]]$tin_dm, yin = cached_grids[[ell]]$yin_dm, 
                                            xin = cached_grids[[ell]]$xin_dm, tout = jumpset_l_grid, npoly = npoly)
            alpj_est_lr_r = CPPlwls2d_s2_LR(bw = h_d, kernel_type = kernel, win = win,
                                            tin = cached_grids[[ell]]$tin_dm, yin = cached_grids[[ell]]$yin_dm, 
                                            xin = cached_grids[[ell]]$xin_dm, tout = jumpset_r_grid, npoly = npoly)
            
            alpj_jumpsize_h_d = alpj_est_lr_r[,p+1+j] - alpj_est_lr_l[,j]
            alpj_jumpsize_h_d = sapply(1:length(jumpset_l), 
                                       function(ii) mean( alpj_jumpsize_h_d[cached_grids[[ell]]$jumpGrid[,1] == cached_grids[[ell]]$D_h[ii]] ))
            
            ## refine the jump location
            if(refined){
              alpj_jumptime = alpj_jumptime[abs(alpj_jumpsize_h_d) > zeta]
              alpj_jumpsize_h_d = alpj_jumpsize_h_d[abs(alpj_jumpsize_h_d) > zeta]
            }
          }
          
          if(length(alpj_jumptime)==0){
            res[[j]][[ell]]$jumptime = NULL
            res[[j]][[ell]]$jumpsize = NULL
            res[[j]][[ell]]$jumpsize_h_tau = alpj_jumpsize_h_tau
            if (!silent) cat('  alpha_', j-1, 'no.', ell, 'dimension --- no change point after refine.', '\n')
          } else {
            res[[j]][[ell]]$jumptime = alpj_jumptime
            res[[j]][[ell]]$jumpsize = alpj_jumpsize_h_d
            res[[j]][[ell]]$jumpsize_h_tau = alpj_jumpsize_h_tau
            if (!silent) cat('  alpha_', j-1, 'no.', ell, 'dimension --- change point at:', 
                             res[[j]][[ell]]$jumptime, '  jump size:', res[[j]][[ell]]$jumpsize, '\n')
          }
        }
        
      }
      
      # Calculate Hausdorff distance for convergence
      for (j in 1:(p+1)) {
        for (ell in 1:m) {
          t_old = old_jumps[[j]][[ell]]
          t_new = res[[j]][[ell]]$jumptime
          
          if (is.null(t_old) && is.null(t_new)) { dist = 0 }
          else if (is.null(t_old) || is.null(t_new)) { dist = Inf }
          else {
            dist1 = max(sapply(t_old, function(x) min(abs(x - t_new))))
            dist2 = max(sapply(t_new, function(x) min(abs(x - t_old))))
            dist = max(dist1, dist2)
          }
          max_dist = max(max_dist, dist)
        }
      }
      
      if (!silent) cat('  Hausdorff distance =', max_dist, '\n')
      
      if (max_dist < tol_hausdorff) {
        if (!silent) cat('Refining phase converged in', iter, 'iterations.\n')
        break
      }
    }
  }
  
  # ---------------------------------------------------------
  # PHASE 3: JUMP SIZE ESTIMATION & UPDATING yyin
  # ---------------------------------------------------------
  
  if (!silent) cat('-----------------------------------------\n', 
                   'Final Estimation\n',
                   '---------------\n')
  
  for (ell in 1:m) {
    for (j in 1:(p+1)) {
      alpj_jumptime = res[[j]][[ell]]$jumptime
      
      if (is.null(alpj_jumptime)) {
        if (!silent) cat(' alpha_', j-1, 'no.', ell, 'dimension --- no change point.\n')
      } else {
        rho_d = 0.25*(h_tau^2 + sqrt(log(n)/(n*h_tau))*h_tau)
        names(rho_d) = c("rho")
        jumpset_l = c(alpj_jumptime - rho_d)
        jumpset_r = c(alpj_jumptime + rho_d)
        
        othergrid = apply(tin[, -c(1), drop=F], 2, function(x){
          seq( max(min(x), h_tau), min(max(x), 1- h_tau), length.out = NbGrid)
        }, simplify = F)
        
        jumpset_l_grid = as.matrix(expand.grid(c(othergrid, list(jumpset_l))))
        jumpset_r_grid = as.matrix(expand.grid(c(othergrid, list(jumpset_r))))
        if (m>1) {
          jumpset_l_grid = jumpset_l_grid[,c(m, 1:(m-1))]
          jumpset_r_grid = jumpset_r_grid[,c(m, 1:(m-1))]
        }
        
        alpj_est_lr_l = CPPlwls2d_s2_LR(bw = h_d, kernel_type = kernel, win = win,
                                        tin = cached_grids[[ell]]$tin_dm, yin = cached_grids[[ell]]$yin_dm, 
                                        xin = cached_grids[[ell]]$xin_dm, tout = jumpset_l_grid, npoly = npoly)
        alpj_est_lr_r = CPPlwls2d_s2_LR(bw = h_d, kernel_type = kernel, win = win,
                                        tin = cached_grids[[ell]]$tin_dm, yin = cached_grids[[ell]]$yin_dm, 
                                        xin = cached_grids[[ell]]$xin_dm, tout = jumpset_r_grid, npoly = npoly)
        
        alpj_jumpsize_h_d = alpj_est_lr_r[,p+1+j] - alpj_est_lr_l[,j]
        alpj_jumpsize_h_d = sapply(1:length(jumpset_l), 
                                   function(ii) mean( alpj_jumpsize_h_d[cached_grids[[ell]]$jumpGrid[,1] == cached_grids[[ell]]$D_h[ii]] ))
        
        ## refine the jump location
        if(refined){
          alpj_jumptime = alpj_jumptime[abs(alpj_jumpsize_h_d) > zeta]
          alpj_jumpsize_h_d = alpj_jumpsize_h_d[abs(alpj_jumpsize_h_d) > zeta]
        }
        
        res[[j]][[ell]]$jumptime = alpj_jumptime
        res[[j]][[ell]]$jumpsize = alpj_jumpsize_h_d
        
        if(length(alpj_jumptime)==0){
          if (!silent) cat(' alpha_', j-1, 'no.', ell, 'dimension --- no change point after refine.', '\n')
        } else {
          if (!silent) cat(' alpha_', j-1, 'no.', ell, 'dimension --- change point at:', 
                           alpj_jumptime, '  jump size:', alpj_jumpsize_h_d, '\n')
        }
        
        
        # Finally subtract the refined jumps from yyin
        if (is.null(alpj_jumptime)) {
          yyin = yyin
        } else{
          ysubstr = sapply(tin[, ell], 
                           function(z) sum(alpj_jumpsize_h_d*(z >= alpj_jumptime)))
          yyin = yyin - ysubstr * cbind(1,xin)[,j]
        }
      }
    }
  }
  
  # ---------------------------------------------------------
  # PHASE 4: FINAL SMOOTHING
  # ---------------------------------------------------------
  gamma = CPPlwls2d_s2(bw = h_2, kernel_type = kernel, win = win,
                       tin = tin, yin = yyin, xin = xin, tout = tin, npoly = npoly)
  
  gamma_out = CPPlwls2d_s2(bw = h_2, kernel_type = kernel, win = win,
                           tin = tin, yin = yyin, xin = xin, tout = tout, npoly = npoly)
  
  for (j in 1:(p+1)) {
    
    res[[j]]$gamma.hat = gamma[,j]
    res[[j]]$gamma.hat_out = gamma_out[,j]
    
    jump = rep(0, dim(tin)[1])
    jump_out = rep(0, dim(tout)[1])
    
    for (ell in 1:m) {
      if (!is.null(res[[j]][[ell]]$jumptime))
        jump = jump + sapply(tin[, ell], function(z) sum(res[[j]][[ell]]$jumpsize * (z >= res[[j]][[ell]]$jumptime)))
      if (!is.null(res[[j]][[ell]]$jumptime))
        jump_out = jump_out + sapply(tout[, ell], function(z) sum(res[[j]][[ell]]$jumpsize * (z >= res[[j]][[ell]]$jumptime)))
    }
    
    res[[j]]$alp.hat = gamma[,j] + jump
    res[[j]]$alp.hat_out = gamma_out[,j] + jump_out
    
  }
  
  if (is.null(xout)) {
    muout = NULL
  } else {
    alp_out = sapply(res, function(z) z$alp.hat_out)
    muout = rowSums(alp_out * cbind(1, xout)) 
  }
  
  return(list(alp_est = res, muout = muout, xout = xout, tout = tout,
              obsGrid = tin, win=win, h_tau = h_tau, h_d = h_d, zeta = zeta, 
              rho_d = rho_d, h_2 = h_2, 
              iter = iter,
              timings = timings))
}



