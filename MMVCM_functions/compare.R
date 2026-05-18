

func_path_int = '/path/to/your/code/'

Rcpp::sourceCpp(paste0(func_path_int, 'CPPlwls2d_s1.cpp'))
Rcpp::sourceCpp(paste0(func_path_int, 'CPPlwls2d_s2.cpp'))

## main functions

SVCoef_func_2d <- function(tin, yin, xin, zin, wi = NULL,
                           bw.seq1 = NULL, bw.seq2 = NULL,
                           nRegGrid = 101, kernel='epan',
                           npoly = 1, nder = 0, kFolds = 5, 
                           method_opt=c('average','profile')){
  
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
        # if (nan_rate > 0.2) cat('Warning: Bandwidth: h1 =', bw.seq1[j], 'too small for stage 1 CV.')
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
      if (nan_rate > 0.05) cv1[j,i] = Inf
      
      if(is.na(cv1[j,i]) || is.nan(cv1[j,i])){
        cv1[j,i] = Inf;
      }
      
    }
  }
  
  
  if(min(cv1) == Inf){
    stop("All bandwidths resulted in infinite CV costs. (Stage 1)")
  }
  
  cvMean = apply(cv1, c(1), mean)
  cvMeanid = which(cvMean == min(cvMean), arr.ind=TRUE)
  if (length(cvMeanid)>1) cvMeanid = cvMeanid[ceiling(length(cvMeanid)/2)]
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
  
  
  
  ##### Stage 2 Estimation for the nonparametric part
  
  ###
  # firsttsCVmu <- Sys.time()
  ## mean function
  if(length(bw.seq2) == 1) {
    h_2 = bw.seq2
  } else {
    ## select by CV procedure
    tunings = CVbandwidth_varycoef_2d(bw.seq = bw.seq2, win = win,
                                      tin = tin, yin = yin_s2, xin = xin,
                                      npoly = npoly, nder = nder, 
                                      kFolds = kFolds, kernel = kernel)
    h_2 = unname(tunings$bopt[1])
  }
  
  ##
  cat(' ', '\n')
  cat('Final selected: ', 'h_1 = ', h_1, ';',  'h_2 = ', h_2, '\n')
  cat('beta_hat:', beta_hat, '\n')
  
  ###
  obsGrid = tin
  ## cut in the region [h, 1-h]
  regGrid = apply(obsGrid, 2, function(x){
    seq( max(min(x), h_2), min(max(x), 1- h_2), length.out = nRegGrid)
  })
  regGrid = as.matrix(expand.grid(as.data.frame(regGrid)))
  
  # Get the mean function using the bandwith estimated above:
  alpha_hat = CPPlwls2d_s2(bw = h_2, 
                           tin = tin, yin = yin_s2, xin = xin, win = win, 
                           tout = obsGrid, kernel_type = kernel, 
                           npoly = npoly)
  
  alpha_hat_out = CPPlwls2d_s2(bw = h_2, 
                               tin = tin, yin = yin_s2, xin = xin, win = win, 
                               tout = regGrid, kernel_type = kernel, 
                               npoly = npoly)
  
  yhat <- c(zin %*% beta_hat + rowSums(alpha_hat * cbind(1, xin)))
  
  mse <- mean((yhat[!is.na(yhat)]-yin[!is.na(yhat)])^2)
  
  cat('MSE: ', mse, '\n')
  
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
  res <- list(beta_hat = beta_hat, yhat = yhat, mse = mse,
              alpha_hat = alpha_hat, alpha_hat_out = alpha_hat_out,
              h_1 = h_1, h_2 = h_2,
              tin = tin, yin = yin, xin = xin, zin = zin, win = win,
              bw.seq1 = bw.seq1, bw.seq2 = bw.seq2,
              obsGrid = obsGrid,
              workGrid = regGrid,
              timings = timings)
  
  return(res)
  
  
}






VCoef_func_2d <- function(tin, yin, xin, wi = NULL,
                          bw.seq = NULL,
                          nRegGrid = 101, kernel='epan',
                          npoly = 1, nder = 0, kFolds = 5){
  
  # firsttsbreaks <- Sys.time()
  # Check the data validity for further analysis
  # CheckData(Ly,Lt)
  
  # inputData <- HandleNumericsAndNAN(Ly,Lt);
  # Ly <-  inputData$Ly;
  # Lt <-  inputData$Lt;
  
  n = dim(xin)[1]
  p = dim(xin)[2]
  m = dim(tin)[2]
  
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
  tin = tin[order(tin[,1]), , drop=F]
  
  
  ###
  # firsttsCVmu <- Sys.time()
  ## mean function
  if(length(bw.seq) == 1) {
    h = bw.seq
  } else {
    ## select by CV procedure
    tunings = CVbandwidth_varycoef_2d(bw.seq = bw.seq, win = win,
                                      tin = tin, yin = yin, xin = xin,
                                      npoly = npoly, nder = nder, 
                                      kFolds = kFolds, kernel = kernel)
    h = unname(tunings$bopt[1])
  }
  
  ##
  cat(' ', '\n')
  cat('Final selected: ', 'h = ', h, '\n')
  
  ###
  obsGrid = tin
  ## cut in the region [h, 1-h]
  regGrid = apply(obsGrid, 2, function(x){
    seq( max(min(x), h), min(max(x), 1- h), length.out = nRegGrid)
  })
  regGrid = as.matrix(expand.grid(as.data.frame(regGrid)))
  
  # Get the mean function using the bandwith estimated above:
  alpha_hat = CPPlwls2d_s2(bw = h, 
                           tin = tin, yin = yin, xin = xin, win = win, 
                           tout = obsGrid, kernel_type = kernel, 
                           npoly = npoly)
  
  alpha_hat_out = CPPlwls2d_s2(bw = h, 
                               tin = tin, yin = yin, xin = xin, win = win, 
                               tout = regGrid, kernel_type = kernel, 
                               npoly = npoly)
  
  yhat <- c(rowSums(alpha_hat * cbind(1, xin)))
  
  mse <- mean((yhat[!is.na(yhat)]-yin[!is.na(yhat)])^2)
  
  cat('MSE: ', mse, '\n')
  
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
  res <- list(yhat = yhat, mse = mse,
              alpha_hat = alpha_hat, alpha_hat_out = alpha_hat_out,
              h = h,
              tin = tin, yin = yin, xin = xin, win = win,
              bw.seq = bw.seq,
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




######################## k-fold cross validation to select h_tau and h_d
CVbandwidth_varycoef_2d <- function(bw.seq = NULL, win = win,
                                    tin = tin, yin = yin, xin = xin,
                                    npoly=1, nder=0, kFolds=5, kernel='epan'){
  
  
  m = dim(tin)[2]
  n = dim(xin)[1]
  p = dim(xin)[2]
  
  cv2 = array(Inf, dim = c(length(bw.seq), kFolds));
  
  theFolds = SimpleFolds(1:n, kFolds)
  
  for (j in 1:length(bw.seq)) {
    
    # cat('CV procedure: h_2 = ', bw.seq[j], '\n')
    
    for (i in 1:kFolds) {
      
      ttest <- tin[theFolds[[i]], , drop=F]
      ytest <- yin[theFolds[[i]]]
      xtest <- xin[theFolds[[i]], , drop=F]
      ttrain <- tin[-theFolds[[i]], , drop=F]
      ytrain <- yin[-theFolds[[i]]]
      xtrain <- xin[-theFolds[[i]], , drop=F]
      wtrain <- win[-theFolds[[i]]]
      
      coef_cv2 = tryCatch(
        CPPlwls2d_s2(bw = bw.seq[j], 
                     tin = ttrain, yin = ytrain, xin = xtrain, win = wtrain, 
                     tout = ttest, kernel_type = kernel, 
                     npoly = npoly, nder = nder),
        error=function(err) {
          warning('Invalid bandwidth during stage 2 CV. Try enlarging the window size.')
          return(Inf)
        })
      nan_rate <- length(which(is.nan(coef_cv2[,1]))) / length(ytest)
      # if (nan_rate > 0.2) cat('Warning: Bandwidth: h2 =', bw.seq[j], 'too small for varying coef CV.')
      tid <- which(!is.nan(coef_cv2[,1]))
      alp_cv2 = coef_cv2[tid, , drop=F]
      
      yhat_cv2 = rowMeans(cbind(1, xtest[tid, , drop=F]) * alp_cv2) 
      
      cv2[j,i] = sum((ytest[tid] - c(yhat_cv2))^2)
      if (nan_rate > 0.01) cv2[j,i] = Inf
      # cv2[j,k,i] = trapzRcpp(xout, (obs - muout)^2)
      # print(cv2)
      if(is.na(cv2[j,i])){
        cv2[j,i] = Inf;
      }
      
    }
    
  }
  
  
  if(min(cv2) == Inf){
    stop("All bandwidths resulted in infinite CV costs. (Stage 2)")
  }
  
  cvMean = apply(cv2, c(1), mean)
  cvMeanid = which(cvMean == min(cvMean))
  if (length(cvMeanid)>1) cvMeanid = cvMeanid[length(cvMeanid)]
  bopt2 = bw.seq[cvMeanid];
  names(bopt2) = c('h_2')
  
  boptList <- list('bopt' = bopt2)
  
  return(boptList)
  
}






################# seqMS

IC1<-function(dat,r)
{
  
  y=as.matrix(dat[,1])
  z=dat[,2]
  X=dat[,-(1:2)]
  
  n=dim(X)[1]
  p=dim(X)[2]
  
  ind=which(z<=r)
  
  if (length(ind)<=p+1 | length(ind)>=n-p-1) {return(Inf)}
  
  X1=X
  X2=X
  X1[-ind,]=0
  X2[ind,]=0
  
  S=t(y)%*%y-t(y)%*%X1%*%solve(t(X1)%*%X1+diag((0.1)^10,p,p))%*%t(X1)%*%y-t(y)%*%X2%*%solve(t(X2)%*%X2+diag((0.1)^10,p,p))%*%t(X2)%*%y
  ic=log(S)+log(n)/n*2*p
  
  if (is.nan(ic)) {ic=Inf}
  
  return(ic)
  
}


IC0<-function(dat)
{
  
  y=as.matrix(dat[,1])
  z=dat[,2]
  X=dat[,-(1:2)]
  
  n=dim(X)[1]
  p=dim(X)[2]
  
  
  S=t(y)%*%y-t(y)%*%X%*%solve(t(X)%*%X+diag((0.1)^100,p,p))%*%t(X)%*%y
  ic=log(S)+log(n)/n*p
  
  return(ic)
  
}


seqMS<-function(dat,bound=c(-2,2))
{
  
  y=as.matrix(dat[,1])
  z=dat[,2]
  X=dat[,-(1:2), drop=F]
  
  n=dim(X)[1]
  p=dim(X)[2]
  
  thre=c()
  segs=list(dat)
  bounds=list(bound)
  
  while(length(segs)>0){
    
    cur=segs[[1]]
    cur_bound=bounds[[1]]
    lwb=cur_bound[1]
    upb=cur_bound[2]
    segs=segs[-1]
    bounds=bounds[-1]
    
    rs=seq(lwb,upb,length.out = 100)
    ic1s=rs
    
    for (i in 1:100){
      ic1s[i]=IC1(cur,rs[i])
    }
    
    ind=which(ic1s==min(ic1s))
    ic1=ic1s[ind][1]
    r=rs[ind][1]
    
    ic0=IC0(cur)
    
    if (ic0<=ic1) {next}
    
    thre=append(thre,r)
    
    lseg=cur[which(cur[,2]<=r),]
    rseg=cur[which(cur[,2]>r),]
    segs[[length(segs)+1]]=lseg
    segs[[length(segs)+1]]=rseg
    bounds[[length(bounds)+1]]=c(lwb,r)
    bounds[[length(bounds)+1]]=c(r,upb)
    
  }
  
  
  m=length(thre)
  
  thre=sort(thre)
  thresholds=thre
  thre=c(-Inf,thre,Inf)
  
  beta=matrix(0,nrow=m+1,ncol=p)
  alpha = matrix(0,nrow=n,ncol=p)
  res=rep(0,n)
  yhat = rep(0,n)
  
  for (i in 1:(m+1)){
    
    ind=which(z<=thre[i+1] & z>thre[i])
    Xi=X[ind,]
    yi=as.matrix(y[ind,])
    
    betai=solve(t(Xi)%*%Xi+diag((0.1)^10,p,p))%*%t(Xi)%*%yi
    beta[i,]=betai
    
    for (ii in 1:length(ind)) {
      alpha[ind[ii],] = betai
    }
    
    yhat[ind] = Xi%*%betai
    res[ind] = yi-yhat[ind]
    
  }
  
  mse = mean(res^2)
  
  return(list(num=m,thresholds=thresholds,beta=beta,alpha=alpha,
              y=y, t=z, X=X, yhat = yhat, 
              res=res, mse = mse))
  
}


seqMS.predict<-function(dat,thresholds,beta)
{
  
  y=as.matrix(dat[,1])
  z=dat[,2]
  X=dat[,-(1:2)]
  
  n=dim(X)[1]
  p=dim(X)[2]
  
  m=length(thresholds)
  thre=c(-Inf,thresholds,Inf)
  
  res=rep(0,n)
  yhat = rep(0,n)
  
  for (i in 1:(m+1)){
    
    ind=which(z<=thre[i+1] & z>thre[i])
    
    Xi=X[ind,]
    yi=as.matrix(y[ind,])
    
    yhat[ind] = Xi%*%beta[i,]
    res[ind] = yi-yhat[ind]
    
  }
  
  mse.te = mean(res^2)
  
  return(list(y.te=y, t.te=z, X.te=X, yhat.te = yhat, 
              res.te = res, mse.te = mse.te))
  
}

