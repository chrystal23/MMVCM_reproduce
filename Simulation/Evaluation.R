

load("exS_model.RData")  # load("exM_model.RData")

load("Res_Ex2.RData")  # load("Res_Ex3.RData")


#
id <- which(sapply(Res_ex2_n500, function(x) length(x)!=1))
Res_ex2_n500 <- Res_ex2_n500[id]
id <- which(sapply(Res_ex2_n1000, function(x) length(x)!=1))
Res_ex2_n1000 <- Res_ex2_n1000[id]
id <- which(sapply(Res_ex2_n2000, function(x) length(x)!=1))
Res_ex2_n2000 <- Res_ex2_n2000[id]



## M hat

alpha_id = 2
T_id = 1

M.hat1 <- sapply( Res_ex2_n500, function(x) length(x$alp_est[[alpha_id+1]][[T_id]]$jumptime) )
table(M.hat1)
M.hat2 <- sapply( Res_ex2_n1000, function(x) length(x$alp_est[[alpha_id+1]][[T_id]]$jumptime) )
table(M.hat2)
M.hat3 <- sapply( Res_ex2_n2000, function(x) length(x$alp_est[[alpha_id+1]][[T_id]]$jumptime) )
table(M.hat3)
cat(sprintf('%.3f (%.3f)',mean(M.hat1),sd(M.hat1)), '&',
    sprintf('%.3f (%.3f)',mean(M.hat2),sd(M.hat2)), '&',
    sprintf('%.3f (%.3f)',mean(M.hat3),sd(M.hat3)) )
M=2
cat(sprintf('%.1f',length(which(M.hat1==M))/length(M.hat1)*100), '&',
    sprintf('%.1f',length(which(M.hat2==M))/length(M.hat2)*100), '&',
    sprintf('%.1f',length(which(M.hat3==M))/length(M.hat3)*100) )


## tau Hausdorf distance

hds <- function (P, Q) {
  stopifnot(is.numeric(P), is.numeric(Q))
  if (is.vector(P))
    P <- matrix(P, ncol = 1)
  if (is.vector(Q))
    Q <- matrix(Q, ncol = 1)
  if (ncol(P) != ncol(Q))
    stop("'P' and 'Q' must have the same number of columns.") 
  # the number of columns refers to the number of dimensions of the data points
  D <- pracma::distmat(P, Q)
  dhd_PQ <- max(apply(D, 1, min))
  dhd_QP <- max(apply(D, 2, min))
  return(cbind(dhd_PQ, dhd_QP))
}

#
alpha_id = 2
T_id = 1

tau = model$alpha[[alpha_id+1]]$tau[[T_id]]

Locofcp1 <- sapply(Res_ex2_n500, function(x) x$alp_est[[alpha_id+1]][[T_id]]$jumptime, simplify = FALSE)
Locofcp1 = Locofcp1[which(sapply(Locofcp1, function(a) is.numeric(a) ))]
loc_dis1 <- t(sapply(Locofcp1, hds, Q = tau))
tau.hd1 <- apply(loc_dis1, 1, max)

Locofcp2 <- sapply(Res_ex2_n1000, function(x) x$alp_est[[alpha_id+1]][[T_id]]$jumptime, simplify = FALSE)
Locofcp2 = Locofcp2[which(sapply(Locofcp2, function(a) is.numeric(a) ))]
loc_dis2 <- t(sapply(Locofcp2, hds, Q = tau))
tau.hd2 <- apply(loc_dis2, 1, max)

Locofcp3 <- sapply(Res_ex2_n2000, function(x) x$alp_est[[alpha_id+1]][[T_id]]$jumptime, simplify = FALSE)
Locofcp3 = Locofcp3[which(sapply(Locofcp3, function(a) is.numeric(a) ))]
loc_dis3 <- t(sapply(Locofcp3, hds, Q = tau))
tau.hd3 <- apply(loc_dis3, 1, max)

cat(sprintf('%.4f (%.4f)',mean(tau.hd1),sd(tau.hd1)), '&',
    sprintf('%.4f (%.4f)',mean(tau.hd2),sd(tau.hd2)), '&',
    sprintf('%.4f (%.4f)',mean(tau.hd3),sd(tau.hd3)) )


## jump size bias

alpha_id = 2
T_id = 1

disc = model$alpha[[alpha_id+1]]$dd[[T_id]]
tau = model$alpha[[alpha_id+1]]$tau[[T_id]]

Locofcp1 <- sapply(Res_ex2_n500, function(x) x$alp_est[[alpha_id+1]][[T_id]]$jumptime, simplify = FALSE)
id = which(sapply(Locofcp1, function(x){length(x)==length(tau)}))
jump_size_est1 <- sapply(Res_ex2_n500[id], function(x) x$alp_est[[alpha_id+1]][[T_id]]$jumpsize, simplify = FALSE)
max_jump_size_est1 <- sapply(jump_size_est1, function(x){max(abs(x-disc))}) 

Locofcp2 <- sapply(Res_ex2_n1000, function(x) x$alp_est[[alpha_id+1]][[T_id]]$jumptime, simplify = FALSE)
id = which(sapply(Locofcp2, function(x){length(x)==length(tau)}))
jump_size_est2 <- sapply(Res_ex2_n1000[id], function(x) x$alp_est[[alpha_id+1]][[T_id]]$jumpsize, simplify = FALSE)
max_jump_size_est2 <- sapply(jump_size_est2, function(x){max(abs(x-disc))}) 

Locofcp3 <- sapply(Res_ex2_n2000, function(x) x$alp_est[[alpha_id+1]][[T_id]]$jumptime, simplify = FALSE)
id = which(sapply(Locofcp3, function(x){length(x)==length(tau)}))
jump_size_est3 <- sapply(Res_ex2_n2000[id], function(x) x$alp_est[[alpha_id+1]][[T_id]]$jumpsize, simplify = FALSE)
max_jump_size_est3 <- sapply(jump_size_est3, function(x){max(abs(x-disc))}) 

cat(sprintf('%.4f (%.4f)',mean(max_jump_size_est1),sd(max_jump_size_est1)), '&',
    sprintf('%.4f (%.4f)',mean(max_jump_size_est2),sd(max_jump_size_est2)), '&',
    sprintf('%.4f (%.4f)',mean(max_jump_size_est3),sd(max_jump_size_est3)) )



### functional coefficient MISE

load("/home/pan/MMVCM/results_used/ex3_model_2.RData")

alpha_id = 2

afun = model$alpfun[[alpha_id+1]]
res = Res_ex3_tGaus_n2000[[1]]
pdat = as.data.frame(cbind(res$obsGrid, apply(res$obsGrid, 1, afun)))
colnames(pdat) = c('t1', 't2', 'alpha')
lattice::cloud(alpha~t1*t2, data = pdat, zlab=list(rot=90), scales = list(arrows=F))

mise1 = sapply(Res_ex3_n500, 
               function(res) mean((res$alp_est[[alpha_id+1]]$alp.hat - apply(res$obsGrid, 1, afun))^2) )
mise2 = sapply(Res_ex3_n1000, 
               function(res) mean((res$alp_est[[alpha_id+1]]$alp.hat - apply(res$obsGrid, 1, afun))^2) )
mise3 = sapply(Res_ex3_n2000, 
               function(res) mean((res$alp_est[[alpha_id+1]]$alp.hat - apply(res$obsGrid, 1, afun))^2) )

cat(sprintf('%.4f (%.4f)',mean(mise1),sd(mise1)), '&',
    sprintf('%.4f (%.4f)',mean(mise2),sd(mise2)), '&',
    sprintf('%.4f (%.4f)',mean(mise3),sd(mise3)) )

#
alpha_id = 2

afun = model$alpfun[[alpha_id+1]] 
res = Res_SVCoef_ex3_n2000[[1]]
pdat = as.data.frame(cbind(res$obsGrid, apply(res$obsGrid, 1, afun)))
colnames(pdat) = c('t1', 't2', 'alpha')
lattice::cloud(alpha~t1*t2, data = pdat, zlab=list(rot=90), scales = list(arrows=F))

mise1 = sapply(Res_SVCoef_ex3_n500, 
               function(res) mean((res$alpha_hat[,alpha_id+1] - apply(res$obsGrid, 1, afun))^2) ); mise1 = mise1[!is.nan(mise1) & mise1<1]
mise2 = sapply(Res_SVCoef_ex3_n1000, 
               function(res) mean((res$alpha_hat[,alpha_id+1] - apply(res$obsGrid, 1, afun))^2) ); mise2 = mise2[!is.nan(mise2) & mise2<1]
mise3 = sapply(Res_SVCoef_ex3_n2000, 
               function(res) mean((res$alpha_hat[,alpha_id+1] - apply(res$obsGrid, 1, afun))^2) ); mise3 = mise3[!is.nan(mise3) & mise3<1]

cat(sprintf('%.4f (%.4f)',mean(mise1),sd(mise1)), '&',
    sprintf('%.4f (%.4f)',mean(mise2),sd(mise2)), '&',
    sprintf('%.4f (%.4f)',mean(mise3),sd(mise3)) )

#
alpha_id = 2

afun = model$alpfun[[alpha_id+1]]
res = Res_seqMS_ex3_n500[[1]]
pdat = as.data.frame(cbind(res$obsGrid, apply(res$obsGrid, 1, afun)))
colnames(pdat) = c('t1', 't2', 'alpha')
lattice::cloud(alpha~t1*t2, data = pdat, zlab=list(rot=90), scales = list(arrows=F))
pdat1 = as.data.frame(cbind(res$obsGrid, res$alpha[,alpha_id+model$q+1]))
colnames(pdat1) = c('t1', 't2', 'alpha')
lattice::cloud(alpha~t1*t2, data = pdat1, zlab=list(rot=90), scales = list(arrows=F))

mise1 = sapply(Res_seqMS_ex3_n500, 
               function(res) mean((res$alpha[,alpha_id+model$q+1] - apply(res$obsGrid, 1, afun))^2) )
mise2 = sapply(Res_seqMS_ex3_n1000, 
               function(res) mean((res$alpha[,alpha_id+model$q+1] - apply(res$obsGrid, 1, afun))^2) )
mise3 = sapply(Res_seqMS_ex3_n2000, 
               function(res) mean((res$alpha[,alpha_id+model$q+1] - apply(res$obsGrid, 1, afun))^2) )

cat(sprintf('%.4f (%.4f)',mean(mise1),sd(mise1)), '&',
    sprintf('%.4f (%.4f)',mean(mise2),sd(mise2)), '&',
    sprintf('%.4f (%.4f)',mean(mise3),sd(mise3)) )


#

library('pracma')

alpha_id = 0

afun = model$alpfun[[alpha_id+1]] 
res = Res_ex2_n2000[[1]]
# plot(res$obsGrid, apply(res$obsGrid, 1, afun))

mise1 = sapply(Res_ex2_n500, 
               function(res) trapz(x = res$obsGrid, 
                                   y = (res$alp_est[[alpha_id+1]] $alp.hat - apply(res$obsGrid, 1, afun))^2) )
mise2 = sapply(Res_ex2_n1000, 
               function(res) trapz(x = res$obsGrid, 
                                   y = (res$alp_est[[alpha_id+1]] $alp.hat - apply(res$obsGrid, 1, afun))^2) )
mise3 = sapply(Res_ex2_n2000, 
               function(res) trapz(x = res$obsGrid, 
                                   y = (res$alp_est[[alpha_id+1]] $alp.hat - apply(res$obsGrid, 1, afun))^2) )

cat(sprintf('%.4f (%.4f)',mean(mise1),sd(mise1)), '&',
    sprintf('%.4f (%.4f)',mean(mise2),sd(mise2)), '&',
    sprintf('%.4f (%.4f)',mean(mise3),sd(mise3)) )


#

afun = model$alpfun[[alpha_id+1]] 
res = Res_VCoef_ex2_n2000[[1]]
plot(res$obsGrid, apply(res$obsGrid, 1, afun))

mise1 = sapply(Res_VCoef_ex2_n500, 
               function(res) trapz(x = res$obsGrid, 
                                   y = (res$alpha_hat[,alpha_id+1] - apply(res$obsGrid, 1, afun))^2) )
mise2 = sapply(Res_VCoef_ex2_n1000, 
               function(res) trapz(x = res$obsGrid, 
                                   y = (res$alpha_hat[,alpha_id+1] - apply(res$obsGrid, 1, afun))^2) )
mise3 = sapply(Res_VCoef_ex2_n2000, 
               function(res) trapz(x = res$obsGrid, 
                                   y = (res$alpha_hat[,alpha_id+1] - apply(res$obsGrid, 1, afun))^2) )

cat(sprintf('%.4f (%.4f)',mean(mise1),sd(mise1)), '&',
    sprintf('%.4f (%.4f)',mean(mise2),sd(mise2)), '&',
    sprintf('%.4f (%.4f)',mean(mise3),sd(mise3)) )

#

afun = model$alpfun[[alpha_id+1]]
res = Res_seqMS_ex2_n2000[[1]]
plot(res$t, apply(as.matrix(res$t, ncol=1), 1, afun))
plot(res$t, res$alpha[,2+model$q+1])

mise1 = sapply(Res_seqMS_ex2_n500, 
               function(res) trapz(x = res$t, 
                                   y = (res$alpha[,alpha_id+model$q+1] - apply(as.matrix(res$t, ncol=1), 1, afun))^2) )
mise2 = sapply(Res_seqMS_ex2_n1000, 
               function(res) trapz(x = res$t, 
                                   y = (res$alpha[,alpha_id+model$q+1] - apply(as.matrix(res$t, ncol=1), 1, afun))^2) )
mise3 = sapply(Res_seqMS_ex2_n2000, 
               function(res) trapz(x = res$t, 
                                   y = (res$alpha[,alpha_id+model$q+1] - apply(as.matrix(res$t, ncol=1), 1, afun))^2) )

cat(sprintf('%.4f (%.4f)',mean(mise1),sd(mise1)), '&',
    sprintf('%.4f (%.4f)',mean(mise2),sd(mise2)), '&',
    sprintf('%.4f (%.4f)',mean(mise3),sd(mise3)) )






## beta bias

beta = model$beta

beta_est1 <- sapply(Res_ex3_n50, function(x) x$beta_hat, simplify = FALSE)
max_beta_est1 <- sapply(beta_est1, function(x){max(abs(x-beta))})

beta_est2 <- sapply(Res_ex3_n100, function(x) x$beta_hat, simplify = FALSE)
max_beta_est2 <- sapply(beta_est2, function(x){max(abs(x-beta))})

beta_est3 <- sapply(Res_ex3_n200, function(x) x$beta_hat, simplify = FALSE)
max_beta_est3 <- sapply(beta_est3, function(x){max(abs(x-beta))})


cat(sprintf('%.4f (%.4f)',mean(max_beta_est1),sd(max_beta_est1)), '&',
    sprintf('%.4f (%.4f)',mean(max_beta_est2),sd(max_beta_est2)), '&',
    sprintf('%.4f (%.4f)',mean(max_beta_est3),sd(max_beta_est3)) )



## mse

mse1 = sapply( Res_ex3_tGaus_n500, function(x) mean((x$yhat - x$yin)^2) )

mse2 = sapply( Res_ex3_tGaus_n1000, function(x) mean((x$yhat - x$yin)^2) )

mse3 = sapply( Res_ex3_tGaus_n2000, function(x) mean((x$yhat - x$yin)^2) )

cat(sprintf('%.4f (%.4f)',mean(mse1),sd(mse1)), '&',
    sprintf('%.4f (%.4f)',mean(mse2),sd(mse2)), '&',
    sprintf('%.4f (%.4f)',mean(mse3),sd(mse3)) )


mse1 = sapply( Res_ex3_n500, function(x) x$mse )

mse2 = sapply( Res_ex3_n1000, function(x) x$mse )

mse3 = sapply( Res_ex3_n2000, function(x) x$mse )

cat(sprintf('%.4f (%.4f)',mean(mse1),sd(mse1)), '&',
    sprintf('%.4f (%.4f)',mean(mse2),sd(mse2)), '&',
    sprintf('%.4f (%.4f)',mean(mse3),sd(mse3)) )



## mse.te

mse1.te = sapply( Res_ex3_tGaus_n500, function(x) x$mse.te )

mse2.te = sapply( Res_ex3_tGaus_n1000, function(x) x$mse.te )

mse3.te = sapply( Res_ex3_tGaus_n2000, function(x) x$mse.te )

cat(sprintf('%.4f (%.4f)',mean(mse1.te),sd(mse1.te)), '&',
    sprintf('%.4f (%.4f)',mean(mse2.te),sd(mse2.te)), '&',
    sprintf('%.4f (%.4f)',mean(mse3.te),sd(mse3.te)) )







## h_1

h_1_1 = sapply( Res_ex3_n500, function(x) x$h_1 )

h_1_2 = sapply( Res_ex3_n1000, function(x) x$h_1 )

h_1_3 = sapply( Res_ex3_n2000, function(x) x$h_1 )

cat(sprintf('%.4f (%.4f)',mean(h_1_1),sd(h_1_1)), '&',
    sprintf('%.4f (%.4f)',mean(h_1_2),sd(h_1_2)), '&',
    sprintf('%.4f (%.4f)',mean(h_1_3),sd(h_1_3)) )


## h_tau

h_tau_1 = sapply( Res_ex3_n500, function(x) x$h_tau )

h_tau_2 = sapply( Res_ex3_n1000, function(x) x$h_tau )

h_tau_3 = sapply( Res_ex3_n2000, function(x) x$h_tau )

cat(sprintf('%.4f (%.4f)',mean(h_tau_1),sd(h_tau_1)), '&',
    sprintf('%.4f (%.4f)',mean(h_tau_2),sd(h_tau_2)), '&',
    sprintf('%.4f (%.4f)',mean(h_tau_3),sd(h_tau_3)) )


## h_d

h_d_1 = sapply( Res_ex3_n500, function(x) x$h_d )

h_d_2 = sapply( Res_ex3_n1000, function(x) x$h_d )

h_d_3 = sapply( Res_ex3_n2000, function(x) x$h_d )

cat(sprintf('%.4f (%.4f)',mean(h_d_1),sd(h_d_1)), '&',
    sprintf('%.4f (%.4f)',mean(h_d_2),sd(h_d_2)), '&',
    sprintf('%.4f (%.4f)',mean(h_d_3),sd(h_d_3)) )


## h_2

h_2_1 = sapply( Res_ex3_n500, function(x) x$h_2 )

h_2_2 = sapply( Res_ex3_n1000, function(x) x$h_2 )

h_2_3 = sapply( Res_ex3_n2000, function(x) x$h_2 )

cat(sprintf('%.4f (%.4f)',mean(h_2_1),sd(h_2_1)), '&',
    sprintf('%.4f (%.4f)',mean(h_2_2),sd(h_2_2)), '&',
    sprintf('%.4f (%.4f)',mean(h_2_3),sd(h_2_3)) )


## zeta / xi_n

zeta_1 = sapply( Res_ex3_n500, function(x) x$zeta )

zeta_2 = sapply( Res_ex3_n1000, function(x) x$zeta )

zeta_3 = sapply( Res_ex3_n2000, function(x) x$zeta )

cat(sprintf('%.4f (%.4f)',mean(zeta_1),sd(zeta_1)), '&',
    sprintf('%.4f (%.4f)',mean(zeta_2),sd(zeta_2)), '&',
    sprintf('%.4f (%.4f)',mean(zeta_3),sd(zeta_3)) )




########################################################
### functional coefficient plot


load("exM_model.RData")

load("Res_Ex3.RData")

### 3D
library('lattice')
library('ggplot2')
library('cowplot')

lattice_option = list( 
  layout.widths = list(left.padding = list(x = -.4, units = "inches"),
                       right.padding = list(x = -.5, units = "inches")),
  layout.heights = list(bottom.padding = list(x = -.4, units = "inches"),
                        top.padding = list(x = -.5, units = "inches")) )

par_setting = list(axis.line=list(col="transparent"),
                   fontsize = list(text = 9))



#
mse1.te = sapply( Res_ex3_n500, function(x) x$mse.te )

mse2.te = sapply( Res_ex3_n1000, function(x) x$mse.te )

mse3.te = sapply( Res_ex3_n2000, function(x) x$mse.te )

##

ord <- order(mse1.te)
tRes <- Res_ex3_n500[ord]
res <- tRes[[floor(length(tRes)/2)]]

pdat1 <- cbind(res$workGrid, 1, sapply(res$alp_est, function(x) x$alp.hat_out))
pdat0 <- cbind(res$workGrid, 0, sapply(model$alpfun, function(f){
  apply(res$workGrid, 1, f) }))
pdat = as.data.frame(rbind(pdat0, pdat1))
colnames(pdat) = c('t1', 't2', 'group', paste0('alpha', 0:model$p))
pd01 <- cloud(alpha0~t1*t2, data = pdat, groups = group, 
              col = c("blue", "red"), pch = 19, cex = 0.6,
              zlim = c(-0.11, 0.25),
              zlab=list(label = expression(alpha[0]), rot=90), scales = list(arrows=F),
              xlab = expression(t[1]), ylab = expression(t[2]),
              lattice.options = lattice_option,
              par.settings = par_setting)
pd01


ord <- order(mse2.te)
tRes <- Res_ex3_n1000[ord]
res <- tRes[[floor(length(tRes)/2)]]

pdat1 <- cbind(res$workGrid, 1, sapply(res$alp_est, function(x) x$alp.hat_out))
pdat0 <- cbind(res$workGrid, 0, sapply(model$alpfun, function(f){
  apply(res$workGrid, 1, f) }))
pdat = as.data.frame(rbind(pdat0, pdat1))
colnames(pdat) = c('t1', 't2', 'group', paste0('alpha', 0:model$p))
pd02 <- cloud(alpha0~t1*t2, data = pdat, groups = group, 
              col = c("blue", "red"), pch = 19, cex = 0.6,
              zlim = c(-0.11, 0.25),
              zlab=list(label = expression(alpha[0]), rot=90), scales = list(arrows=F),
              xlab = expression(t[1]), ylab = expression(t[2]),
              lattice.options = lattice_option,
              par.settings = par_setting)
pd02


ord <- order(mse3.te)
tRes <- Res_ex3_n2000[ord]
res <- tRes[[floor(length(tRes)/2)]] 

# afun = model$alpfun$alpha0
# mise = sapply(tRes, 
#                function(res)  mean((res$alp_est$alpha0$alp.hat - apply(res$obsGrid, 1, afun))^2) )


pdat1 <- cbind(res$workGrid, 1, sapply(res$alp_est, function(x) x$alp.hat_out))
pdat0 <- cbind(res$workGrid, 0, sapply(model$alpfun, function(f){
  apply(res$workGrid, 1, f) }))
pdat = as.data.frame(rbind(pdat0, pdat1))
colnames(pdat) = c('t1', 't2', 'group', paste0('alpha', 0:model$p))
pd03 <- cloud(alpha0~t1*t2, data = pdat, groups = group, 
              col = c("blue", "red"), pch = 19, cex = 0.6,
              zlim = c(-0.11, 0.26),
              zlab=list(label = expression(alpha[0]), rot=90), scales = list(arrows=F),
              xlab = expression(t[1]), ylab = expression(t[2]),
              lattice.options = lattice_option,
              par.settings = par_setting)
pd03

plot_grid(pd01, pd02, pd03)


#

ord <- order(mse1.te)
tRes <- Res_ex3_n500[ord]
res <- tRes[[floor(length(tRes)/2)]]

pdat1 <- cbind(res$workGrid, 1, sapply(res$alp_est, function(x) x$alp.hat_out))
pdat0 <- cbind(res$workGrid, 0, sapply(model$alpfun, function(f){
  apply(res$workGrid, 1, f) }))
pdat = as.data.frame(rbind(pdat0, pdat1))
colnames(pdat) = c('t1', 't2', 'group', paste0('alpha', 0:model$p))
pd11 <- cloud(alpha1~t1*t2, data = pdat, groups = group, 
              col = c("blue", "red"), pch = 19, cex = 0.6,
              zlim = c(0.2, 3.2),
              zlab=list(label = expression(alpha[1]), rot=90), scales = list(arrows=F),
              xlab = expression(t[1]), ylab = expression(t[2]),
              lattice.options = lattice_option,
              par.settings = par_setting)
pd11


ord <- order(mse2.te)
tRes <- Res_ex3_n1000[ord]
res <- tRes[[floor(length(tRes)/2)]]

pdat1 <- cbind(res$workGrid, 1, sapply(res$alp_est, function(x) x$alp.hat_out))
pdat0 <- cbind(res$workGrid, 0, sapply(model$alpfun, function(f){
  apply(res$workGrid, 1, f) }))
pdat = as.data.frame(rbind(pdat0, pdat1))
colnames(pdat) = c('t1', 't2', 'group', paste0('alpha', 0:model$p))
pd12 <- cloud(alpha1~t1*t2, data = pdat, groups = group, 
              col = c("blue", "red"), pch = 19, cex = 0.6,
              zlim = c(0.2, 3.2),
              zlab=list(label = expression(alpha[1]), rot=90), scales = list(arrows=F),
              xlab = expression(t[1]), ylab = expression(t[2]),
              lattice.options = lattice_option,
              par.settings = par_setting)
pd12


ord <- order(mse3.te)
tRes <- Res_ex3_n2000[ord]
res <- tRes[[floor(length(tRes)/2)]]

pdat1 <- cbind(res$workGrid, 1, sapply(res$alp_est, function(x) x$alp.hat_out))
pdat0 <- cbind(res$workGrid, 0, sapply(model$alpfun, function(f){
  apply(res$workGrid, 1, f) }))
pdat = as.data.frame(rbind(pdat0, pdat1))
colnames(pdat) = c('t1', 't2', 'group', paste0('alpha', 0:model$p))
pd13 <- cloud(alpha1~t1*t2, data = pdat, groups = group, 
              col = c("blue", "red"), pch = 19, cex = 0.6,
              zlim = c(0.2, 3.2),
              zlab=list(label = expression(alpha[1]), rot=90), scales = list(arrows=F),
              xlab = expression(t[1]), ylab = expression(t[2]),
              lattice.options = lattice_option,
              par.settings = par_setting)
pd13

plot_grid(pd11, pd12, pd13)


#

ord <- order(mse1.te)
tRes <- Res_ex3_n500[ord]
res <- tRes[[floor(length(tRes)/2)]]

pdat1 <- cbind(res$workGrid, 1, sapply(res$alp_est, function(x) x$alp.hat_out))
pdat0 <- cbind(res$workGrid, 0, sapply(model$alpfun, function(f){
  apply(res$workGrid, 1, f) }))
pdat = as.data.frame(rbind(pdat0, pdat1))
colnames(pdat) = c('t1', 't2', 'group', paste0('alpha', 0:model$p))
pd21 <- cloud(alpha2~t1*t2, data = pdat, groups = group, 
              col = c("blue", "red"), pch = 19, cex = 0.6,
              zlim = c(-0.99, 0.56),
              zlab=list(label = expression(alpha[2]), rot=90), scales = list(arrows=F),
              xlab = expression(t[1]), ylab = expression(t[2]),
              lattice.options = lattice_option,
              par.settings = par_setting)
pd21


ord <- order(mse2.te)
tRes <- Res_ex3_n1000[ord]
res <- tRes[[floor(length(tRes)/2)]]

pdat1 <- cbind(res$workGrid, 1, sapply(res$alp_est, function(x) x$alp.hat_out))
pdat0 <- cbind(res$workGrid, 0, sapply(model$alpfun, function(f){
  apply(res$workGrid, 1, f) }))
pdat = as.data.frame(rbind(pdat0, pdat1))
colnames(pdat) = c('t1', 't2', 'group', paste0('alpha', 0:model$p))
pd22 <- cloud(alpha2~t1*t2, data = pdat, groups = group, 
              col = c("blue", "red"), pch = 19, cex = 0.6,
              zlim = c(-0.99, 0.56),
              zlab=list(label = expression(alpha[2]), rot=90), scales = list(arrows=F),
              xlab = expression(t[1]), ylab = expression(t[2]),
              lattice.options = lattice_option,
              par.settings = par_setting)
pd22


ord <- order(mse3.te)
tRes <- Res_ex3_n2000[ord]
res <- tRes[[floor(length(tRes)/2)]]

pdat1 <- cbind(res$workGrid, 1, sapply(res$alp_est, function(x) x$alp.hat_out))
pdat0 <- cbind(res$workGrid, 0, sapply(model$alpfun, function(f){
  apply(res$workGrid, 1, f) }))
pdat = as.data.frame(rbind(pdat0, pdat1))
colnames(pdat) = c('t1', 't2', 'group', paste0('alpha', 0:model$p))
pd23 <- cloud(alpha2~t1*t2, data = pdat, groups = group, 
              col = c("blue", "red"), pch = 19, cex = 0.6,
              zlim = c(-0.99, 0.56),
              zlab=list(label = expression(alpha[2]), rot=90), scales = list(arrows=F),
              xlab = expression(t[1]), ylab = expression(t[2]),
              lattice.options = lattice_option,
              par.settings = par_setting)
pd23

plot_grid(pd21, pd22, pd23)



#

plot_grid(pd01, pd02, pd03, pd11, pd12, pd13, pd21, pd22, pd23,  nrow = 3, ncol = 3)

ggsave('fig2.eps', width = 25, height = 18, units = 'cm', dpi=1000)

ggsave('fig2.png', width = 25, height = 18, units = 'cm', dpi=1000)






###################################

library('ggplot2')
library('cowplot')


#
id <- which(sapply(Res_ex2_n500, function(x) length(x)!=1))
Res_ex2_n500 <- Res_ex2_n500[id]
id <- which(sapply(Res_ex2_n1000, function(x) length(x)!=1))
Res_ex2_n1000 <- Res_ex2_n1000[id]
id <- which(sapply(Res_ex2_n2000, function(x) length(x)!=1))
Res_ex2_n2000 <- Res_ex2_n2000[id]

#
mse1 = sapply( Res_ex2_n500, function(x) x$mse.te )

mse2 = sapply( Res_ex2_n1000, function(x) x$mse.te )

mse3 = sapply( Res_ex2_n2000, function(x) x$mse.te )

#
tau = model$alpha$alpha0$tau[[1]]

ord <- order(mse1)
tRes <- Res_ex2_n500[ord]
res <- tRes[[floor(length(tRes)/2)-5]]
seg <- c(sapply(res$workGrid, function(x) sum(x>=tau)), 
         sapply(res$workGrid, function(x) 10+sum(x>=res$alp_est$alpha0[[1]]$jumptime)))
color <- as.factor(c(rep(1, length(res$workGrid)), rep(2, length(res$workGrid))))
pdata <- data.frame( t = rep(res$workGrid, 2),
                     alpha = c(apply(res$workGrid, 1, model$alpfun$alpha0), 
                               res$alp_est$alpha0$alp.hat_out),
                     color = color,
                     segment = seg  )
pwCI01 = pwCB_fun(res)

p01 <- ggplot(data=pdata, mapping=aes(x=t, y=alpha, group=segment)) + 
  theme_bw() + theme(panel.grid=element_blank()) +
  geom_ribbon(aes( ymin=c(rep(NA,length(res$workGrid)), pwCI01$alpha0[2,]), 
                   ymax=c(rep(NA,length(res$workGrid)), pwCI01$alpha0[1,]),
                   fill = 'grey85')) +
  geom_line(aes(color=color), size=0.3) +
  scale_color_manual(values=c('blue','red')) +
  scale_fill_manual(values=c('grey85')) +
  theme(legend.position="none") +
  ylim(c(-1.2,1.7)) + 
  ylab(expression(alpha[0])) + xlab(expression(t[1]))
p01



ord <- order(mse2)
tRes <- Res_ex2_n1000[ord]
res <- tRes[[floor(length(tRes)/2)-5]]
seg <- c(sapply(res$workGrid, function(x) sum(x>=tau)), 
         sapply(res$workGrid, function(x) 10+sum(x>=res$alp_est$alpha0[[1]]$jumptime)))
color <- as.factor(c(rep(1, length(res$workGrid)), rep(2, length(res$workGrid))))
pdata <- data.frame( t = rep(res$workGrid, 2),
                     alpha = c(apply(res$workGrid, 1, model$alpfun$alpha0), 
                               res$alp_est$alpha0$alp.hat_out),
                     color = color,
                     segment = seg  )
pwCI02 = pwCB_fun(res)

p02 <- ggplot(data=pdata, mapping=aes(x=t, y=alpha, group=segment)) + 
  theme_bw() + theme(panel.grid=element_blank()) +
  geom_ribbon(aes( ymin=c(rep(NA,length(res$workGrid)), pwCI02$alpha0[2,]), 
                   ymax=c(rep(NA,length(res$workGrid)), pwCI02$alpha0[1,]),
                   fill = 'grey85')) +
  geom_line(aes(color=color), size=0.3) +
  scale_color_manual(values=c('blue','red')) +
  scale_fill_manual(values=c('grey85')) +
  theme(legend.position="none") +
  ylim(c(-1.2,1.7)) +
  ylab(expression(alpha[0])) + xlab(expression(t[1]))
p02


ord <- order(mse3)
tRes <- Res_ex2_n2000[ord]
res <- tRes[[floor(length(tRes)/2)]]
seg <- c(sapply(res$workGrid, function(x) sum(x>=tau)), 
         sapply(res$workGrid, function(x) 10+sum(x>=res$alp_est$alpha0[[1]]$jumptime)))
color <- as.factor(c(rep(1, length(res$workGrid)), rep(2, length(res$workGrid))))
pdata <- data.frame( t = rep(res$workGrid, 2),
                     alpha = c(apply(res$workGrid, 1, model$alpfun$alpha0), 
                               res$alp_est$alpha0$alp.hat_out),
                     color = color,
                     segment = seg  )
pwCI03 = pwCB_fun(res)

p03 <- ggplot(data=pdata, mapping=aes(x=t, y=alpha, group=segment)) + 
  theme_bw() + theme(panel.grid=element_blank()) +
  geom_ribbon(aes( ymin=c(rep(NA,length(res$workGrid)), pwCI03$alpha0[2,]), 
                   ymax=c(rep(NA,length(res$workGrid)), pwCI03$alpha0[1,]),
                   fill = 'grey85')) +
  geom_line(aes(color=color), size=0.3) +
  scale_color_manual(values=c('blue','red')) +
  scale_fill_manual(values=c('grey85')) +
  theme(legend.position="none") +
  ylim(c(-1.2,1.7)) +
  ylab(expression(alpha[0])) + xlab(expression(t[1]))
p03

plot_grid(p01, p02, p03)



#
tau = model$alpha$alpha1$tau[[1]]

ord <- order(mse1)
tRes <- Res_ex2_n500[ord]
res <- tRes[[floor(length(tRes)/2)]]
seg <- c(sapply(res$workGrid, function(x) sum(x>=tau)), 
         sapply(res$workGrid, function(x) 10+sum(x>=res$alp_est$alpha1[[1]]$jumptime)))
color <- as.factor(c(rep(1, length(res$workGrid)), rep(2, length(res$workGrid))))
pdata <- data.frame( t = rep(res$workGrid, 2),
                     alpha = c(apply(res$workGrid, 1, model$alpfun$alpha1), 
                               res$alp_est$alpha1$alp.hat_out),
                     color = color,
                     segment = seg  )
pwCI11 = pwCB_fun(res)

p11 <- ggplot(data=pdata, mapping=aes(x=t, y=alpha, group=segment)) + 
  theme_bw() + theme(panel.grid=element_blank()) +
  geom_ribbon(aes( ymin=c(rep(NA,length(res$workGrid)), pwCI11$alpha1[2,]), 
                   ymax=c(rep(NA,length(res$workGrid)), pwCI11$alpha1[1,]),
                   fill = 'grey85')) +
  geom_line(aes(color=color), size=0.3) +
  scale_color_manual(values=c('blue','red')) +
  scale_fill_manual(values=c('grey85')) +
  theme(legend.position="none") +
  ylim(c(-1.1, 2.1)) +
  ylab(expression(alpha[1])) + xlab(expression(t[1]))
p11



ord <- order(mse2)
tRes <- Res_ex2_n1000[ord]
res <- tRes[[floor(length(tRes)/2)]]
seg <- c(sapply(res$workGrid, function(x) sum(x>=tau)), 
         sapply(res$workGrid, function(x) 10+sum(x>=res$alp_est$alpha1[[1]]$jumptime)))
color <- as.factor(c(rep(1, length(res$workGrid)), rep(2, length(res$workGrid))))
pdata <- data.frame( t = rep(res$workGrid, 2),
                     alpha = c(apply(res$workGrid, 1, model$alpfun$alpha1), 
                               res$alp_est$alpha1$alp.hat_out),
                     color = color,
                     segment = seg  )
pwCI12 = pwCB_fun(res)

p12 <- ggplot(data=pdata, mapping=aes(x=t, y=alpha, group=segment)) + 
  theme_bw() + theme(panel.grid=element_blank()) +
  geom_ribbon(aes( ymin=c(rep(NA,length(res$workGrid)), pwCI12$alpha1[2,]), 
                   ymax=c(rep(NA,length(res$workGrid)), pwCI12$alpha1[1,]),
                   fill = 'grey85')) +
  geom_line(aes(color=color), size=0.3) +
  scale_color_manual(values=c('blue','red')) +
  scale_fill_manual(values=c('grey85')) +
  theme(legend.position="none") +
  ylim(c(-1.1, 2.1)) +
  ylab(expression(alpha[1])) + xlab(expression(t[1]))
p12


ord <- order(mse3)
tRes <- Res_ex2_n2000[ord]
res <- tRes[[floor(length(tRes)/2)]]
seg <- c(sapply(res$workGrid, function(x) sum(x>=tau)), 
         sapply(res$workGrid, function(x) 10+sum(x>=res$alp_est$alpha1[[1]]$jumptime)))
color <- as.factor(c(rep(1, length(res$workGrid)), rep(2, length(res$workGrid))))
pdata <- data.frame( t = rep(res$workGrid, 2),
                     alpha = c(apply(res$workGrid, 1, model$alpfun$alpha1), 
                               res$alp_est$alpha1$alp.hat_out),
                     color = color,
                     segment = seg  )
pwCI13 = pwCB_fun(res)

p13 <- ggplot(data=pdata, mapping=aes(x=t, y=alpha, group=segment)) + 
  theme_bw() + theme(panel.grid=element_blank()) +
  geom_ribbon(aes( ymin=c(rep(NA,length(res$workGrid)), pwCI13$alpha1[2,]), 
                   ymax=c(rep(NA,length(res$workGrid)), pwCI13$alpha1[1,]),
                   fill = 'grey85')) +
  geom_line(aes(color=color), size=0.3) +
  scale_color_manual(values=c('blue','red')) +
  scale_fill_manual(values=c('grey85')) +
  theme(legend.position="none") +
  ylim(c(-1.1, 2.1)) +
  ylab(expression(alpha[1])) + xlab(expression(t[1]))
p13

plot_grid(p11, p12, p13)



#
tau = model$alpha$alpha2$tau[[1]]

ord <- order(mse1)
tRes <- Res_ex2_n500[ord]
res <- tRes[[floor(length(tRes)/2)]]
seg <- c(sapply(res$workGrid, function(x) sum(x>=tau)), 
         sapply(res$workGrid, function(x) 10+sum(x>=res$alp_est$alpha2[[1]]$jumptime)))
color <- as.factor(c(rep(1, length(res$workGrid)), rep(2, length(res$workGrid))))
pdata <- data.frame( t = rep(res$workGrid, 2),
                     alpha = c(apply(res$workGrid, 1, model$alpfun$alpha2), 
                               res$alp_est$alpha2$alp.hat_out),
                     color = color,
                     segment = seg  )
pwCI21 = pwCB_fun(res)

p21 <- ggplot(data=pdata, mapping=aes(x=t, y=alpha, group=segment)) + 
  theme_bw() + theme(panel.grid=element_blank()) +
  geom_ribbon(aes( ymin=c(rep(NA,length(res$workGrid)), pwCI21$alpha2[2,]), 
                   ymax=c(rep(NA,length(res$workGrid)), pwCI21$alpha2[1,]),
                   fill = 'grey85')) +
  geom_line(aes(color=color), size=0.3) +
  scale_color_manual(values=c('blue','red')) +
  scale_fill_manual(values=c('grey85')) +
  theme(legend.position="none") +
  ylim(c(-1.5,1.5)) + 
  ylab(expression(alpha[2])) + xlab(expression(t[1]))
p21



ord <- order(mse2)
tRes <- Res_ex2_n1000[ord]
res <- tRes[[floor(length(tRes)/2)]]
seg <- c(sapply(res$workGrid, function(x) sum(x>=tau)), 
         sapply(res$workGrid, function(x) 10+sum(x>=res$alp_est$alpha2[[1]]$jumptime)))
color <- as.factor(c(rep(1, length(res$workGrid)), rep(2, length(res$workGrid))))
pdata <- data.frame( t = rep(res$workGrid, 2),
                     alpha = c(apply(res$workGrid, 1, model$alpfun$alpha2), 
                               res$alp_est$alpha2$alp.hat_out),
                     color = color,
                     segment = seg  )
pwCI22 = pwCB_fun(res)

p22 <- ggplot(data=pdata, mapping=aes(x=t, y=alpha, group=segment)) + 
  theme_bw() + theme(panel.grid=element_blank()) +
  geom_ribbon(aes( ymin=c(rep(NA,length(res$workGrid)), pwCI22$alpha2[2,]), 
                   ymax=c(rep(NA,length(res$workGrid)), pwCI22$alpha2[1,]),
                   fill = 'grey85')) +
  geom_line(aes(color=color), size=0.3) +
  scale_color_manual(values=c('blue','red')) +
  scale_fill_manual(values=c('grey85')) +
  theme(legend.position="none") +
  ylim(c(-1.5,1.5)) + 
  ylab(expression(alpha[2])) + xlab(expression(t[1]))
p22


ord <- order(mse3)
tRes <- Res_ex2_n2000[ord]
res <- tRes[[floor(length(tRes)/2)]]
seg <- c(sapply(res$workGrid, function(x) sum(x>=tau)), 
         sapply(res$workGrid, function(x) 10+sum(x>=res$alp_est$alpha2[[1]]$jumptime)))
color <- as.factor(c(rep(1, length(res$workGrid)), rep(2, length(res$workGrid))))
pdata <- data.frame( t = rep(res$workGrid, 2),
                     alpha = c(apply(res$workGrid, 1, model$alpfun$alpha2), 
                               res$alp_est$alpha2$alp.hat_out),
                     color = color,
                     segment = seg  )
pwCI23 = pwCB_fun(res)

p23 <- ggplot(data=pdata, mapping=aes(x=t, y=alpha, group=segment)) + 
  theme_bw() + theme(panel.grid=element_blank()) +
  geom_ribbon(aes( ymin=c(rep(NA,length(res$workGrid)), pwCI23$alpha2[2,]), 
                   ymax=c(rep(NA,length(res$workGrid)), pwCI23$alpha2[1,]),
                   fill = 'grey85')) +
  geom_line(aes(color=color), size=0.3) +
  scale_color_manual(values=c('blue','red')) +
  scale_fill_manual(values=c('grey85')) +
  theme(legend.position="none") +
  ylim(c(-1.5,1.5)) + 
  ylab(expression(alpha[2])) + xlab(expression(t[1]))
p23

plot_grid(p21, p22, p23)


plot_grid(p01, p02, p03, p11, p12, p13, p21, p22, p23,  nrow = 3, ncol = 3)
ggsave('fig1.eps', width = 25, height = 18, units = 'cm', dpi=1000)



##### need to calculate mse for example 2, and run alp_fun's first

#
id <- which(sapply(Res_ex2_n500, function(x) length(x)!=1))
Res_ex2_n500 <- Res_ex2_n500[id]
id <- which(sapply(Res_ex2_n1000, function(x) length(x)!=1))
Res_ex2_n1000 <- Res_ex2_n1000[id]
id <- which(sapply(Res_ex2_n2000, function(x) length(x)!=1))
Res_ex2_n2000 <- Res_ex2_n2000[id]

#
mse1 = sapply( Res_ex2_n500, function(x) mean((x$yhat - x$yin)^2) )

mse2 = sapply( Res_ex2_n1000, function(x) mean((x$yhat - x$yin)^2) )

mse3 = sapply( Res_ex2_n2000, function(x) mean((x$yhat - x$yin)^2) )

#
tau0_ex2 = c(.333,.667) # jump locations
disc0_ex2 = c(-.5,.4) # jump sizes
tau0 = tau0_ex2
disc0 = disc0_ex2
smooth0 <- function(x) sin(2*pi*x)
jump0 = function(x) sum((x>=tau0)*disc0) # jump function
alp0_fun= function(t){smooth0(t) + sapply(t, jump0)}

tau1_ex2 = c(.5) # jump locations
disc1_ex2 = c(.5) # jump sizes
tau1 = tau1_ex2
disc1 = disc1_ex2
smooth1 <- function(x) cos(2*pi*x)
jump1 = function(x) sum((x>=tau1)*disc1) # jump function
alp1_fun= function(t){smooth1(t) + sapply(t, jump1)}

alpfun <- list(alp0_fun=alp0_fun, alp1_fun=alp1_fun)
beta <- c(1,1,-1)
beta_ex2 = beta

#

tau = tau0_ex2

ord <- order(mse1)
tRes <- Res_ex2_n500[ord]
res <- tRes[[floor(length(tRes)/2)]]
seg <- c(sapply(res$workGrid, function(x) sum(x>=tau)), 
         sapply(res$workGrid, function(x) 10+sum(x>=res$alp_est[[1]]$jumptime)))
color <- as.factor(c(rep(1, length(res$workGrid)), rep(2, length(res$workGrid))))
pdata <- data.frame( t = rep(res$workGrid, 2),
                     alpha = c(alp0_fun(res$workGrid), res$alp_est[[1]]$alp.hat_out),
                     color = color,
                     segment = seg  )
pwCI4 = pwCB_fun(res)

p4 <- ggplot(data=pdata, mapping=aes(x=t, y=alpha, group=segment)) + 
  theme_bw() + theme(panel.grid=element_blank()) +
  geom_ribbon(aes( ymin=c(rep(NA,length(res$workGrid)), pwCI4[[1]][2,]), 
                   ymax=c(rep(NA,length(res$workGrid)), pwCI4[[1]][1,]),
                   fill = 'grey85')) +
  geom_line(aes(color=color), size=0.3) +
  scale_color_manual(values=c('blue','red')) +
  scale_fill_manual(values=c('grey85')) +
  theme(legend.position="none") +
  ylab(expression(alpha[0]))
p4


ord <- order(mse2)
tRes <- Res_ex2_n1000[ord]
res <- tRes[[floor(length(tRes)/2)]]
seg <- c(sapply(res$workGrid, function(x) sum(x>=tau)), 
         sapply(res$workGrid, function(x) 10+sum(x>=res$alp_est[[1]]$jumptime)))
color <- as.factor(c(rep(1, length(res$workGrid)), rep(2, length(res$workGrid))))
pdata <- data.frame( t = rep(res$workGrid, 2),
                     alpha = c(alp0_fun(res$workGrid), res$alp_est[[1]]$alp.hat_out),
                     color = color,
                     segment = seg  )
pwCI5 = pwCB_fun(res)

p5 <- ggplot(data=pdata, mapping=aes(x=t, y=alpha, group=segment)) + 
  theme_bw() + theme(panel.grid=element_blank()) +
  geom_ribbon(aes( ymin=c(rep(NA,length(res$workGrid)), pwCI5[[1]][2,]), 
                   ymax=c(rep(NA,length(res$workGrid)), pwCI5[[1]][1,]),
                   fill = 'grey85')) +
  geom_line(aes(color=color), size=0.3) +
  scale_color_manual(values=c('blue','red')) +
  scale_fill_manual(values=c('grey85')) +
  theme(legend.position="none") +
  ylab(expression(alpha[0]))
p5


ord <- order(mse3)
tRes <- Res_ex2_n2000[ord]
res <- tRes[[floor(length(tRes)/2)]]
seg <- c(sapply(res$workGrid, function(x) sum(x>=tau)), 
         sapply(res$workGrid, function(x) 10+sum(x>=res$alp_est[[1]]$jumptime)))
color <- as.factor(c(rep(1, length(res$workGrid)), rep(2, length(res$workGrid))))
pdata <- data.frame( t = rep(res$workGrid, 2),
                     alpha = c(alp0_fun(res$workGrid), res$alp_est[[1]]$alp.hat_out),
                     color = color,
                     segment = seg  )
pwCI6 = pwCB_fun(res)

p6 <- ggplot(data=pdata, mapping=aes(x=t, y=alpha, group=segment)) + 
  theme_bw() + theme(panel.grid=element_blank()) +
  geom_ribbon(aes( ymin=c(rep(NA,length(res$workGrid)), pwCI6[[1]][2,]), 
                   ymax=c(rep(NA,length(res$workGrid)), pwCI6[[1]][1,]),
                   fill = 'grey85')) +
  geom_line(aes(color=color), size=0.3) +
  scale_color_manual(values=c('blue','red')) +
  scale_fill_manual(values=c('grey85')) +
  theme(legend.position="none") +
  ylab(expression(alpha[0]))
p6

###

tau = tau1_ex2

ord <- order(mse1)
tRes <- Res_ex2_n500[ord]
res <- tRes[[floor(length(tRes)/2)]]
seg <- c(sapply(res$workGrid, function(x) sum(x>=tau)), 
         sapply(res$workGrid, function(x) 10+sum(x>=res$alp_est[[2]]$jumptime)))
color <- as.factor(c(rep(1, length(res$workGrid)), rep(2, length(res$workGrid))))
pdata <- data.frame( t = rep(res$workGrid, 2),
                     alpha = c(alp1_fun(res$workGrid), res$alp_est[[2]]$alp.hat_out),
                     color = color,
                     segment = seg  )
pwCI7 = pwCB_fun(res)

p7 <- ggplot(data=pdata, mapping=aes(x=t, y=alpha, group=segment)) + 
  theme_bw() + theme(panel.grid=element_blank()) +
  geom_ribbon(aes( ymin=c(rep(NA,length(res$workGrid)), pwCI7[[2]][2,]), 
                   ymax=c(rep(NA,length(res$workGrid)), pwCI7[[2]][1,]),
                   fill = 'grey85')) +
  geom_line(aes(color=color), size=0.3) +
  scale_color_manual(values=c('blue','red')) +
  scale_fill_manual(values=c('grey85')) +
  theme(legend.position="none") +
  ylab(expression(alpha[1]))
p7


ord <- order(mse2)
tRes <- Res_ex2_n1000[ord]
res <- tRes[[floor(length(tRes)/2)]]
seg <- c(sapply(res$workGrid, function(x) sum(x>=tau)), 
         sapply(res$workGrid, function(x) 10+sum(x>=res$alp_est[[2]]$jumptime)))
color <- as.factor(c(rep(1, length(res$workGrid)), rep(2, length(res$workGrid))))
pdata <- data.frame( t = rep(res$workGrid, 2),
                     alpha = c(alp1_fun(res$workGrid), res$alp_est[[2]]$alp.hat_out),
                     color = color,
                     segment = seg  )
pwCI8 = pwCB_fun(res)

p8 <- ggplot(data=pdata, mapping=aes(x=t, y=alpha, group=segment)) + 
  theme_bw() + theme(panel.grid=element_blank()) +
  geom_ribbon(aes( ymin=c(rep(NA,length(res$workGrid)), pwCI8[[2]][2,]), 
                   ymax=c(rep(NA,length(res$workGrid)), pwCI8[[2]][1,]),
                   fill = 'grey85')) +
  geom_line(aes(color=color), size=0.3) +
  scale_color_manual(values=c('blue','red')) +
  scale_fill_manual(values=c('grey85')) +
  theme(legend.position="none") +
  ylab(expression(alpha[1]))
p8


ord <- order(mse3)
tRes <- Res_ex2_n2000[ord]
res <- tRes[[floor(length(tRes)/2)]]
seg <- c(sapply(res$workGrid, function(x) sum(x>=tau)), 
         sapply(res$workGrid, function(x) 10+sum(x>=res$alp_est[[2]]$jumptime)))
color <- as.factor(c(rep(1, length(res$workGrid)), rep(2, length(res$workGrid))))
pdata <- data.frame( t = rep(res$workGrid, 2),
                     alpha = c(alp1_fun(res$workGrid), res$alp_est[[2]]$alp.hat_out),
                     color = color,
                     segment = seg  )
pwCI9 = pwCB_fun(res)

p9 <- ggplot(data=pdata, mapping=aes(x=t, y=alpha, group=segment)) + 
  theme_bw() + theme(panel.grid=element_blank()) +
  geom_ribbon(aes( ymin=c(rep(NA,length(res$workGrid)), pwCI9[[2]][2,]), 
                   ymax=c(rep(NA,length(res$workGrid)), pwCI9[[2]][1,]),
                   fill = 'grey85')) +
  geom_line(aes(color=color), size=0.3) +
  scale_color_manual(values=c('blue','red')) +
  scale_fill_manual(values=c('grey85')) +
  theme(legend.position="none") +
  ylab(expression(alpha[1]))
p9


##
plot_grid(p11, p21, p31, p1, p2, p3, p4, p5, p6, p7, p8, p9, nrow = 4, ncol = 3)
ggsave('fig1.eps', width = 25, height = 22.5, units = 'cm', dpi=1000)







## MSE comparison plot

library('ggplot2')
library('cowplot')

runtime = 200


method = factor( c ( rep('SVCJP', length(Res_ex2_n500)),
                     rep('SVCM', runtime), 
                     rep('VCM', runtime),
                     rep('ThreReg', runtime) ), 
                 levels = c('SVCJP', 'SVCM', 'VCM', 'ThreReg') )
pdata <- data.frame( MSE.tr = c( sapply(Res_ex2_n500, function(x) x$mse),
                                 sapply(Res_SVCoef_ex2_n500, function(x) x$mse),
                                 sapply(Res_VCoef_ex2_n500, function(x) x$mse),
                                 sapply(Res_seqMS_ex2_n500, function(x) x$mse)),
                     method = method )
p1 <- ggplot(data = pdata, mapping = aes(x=method, y=MSE.tr)) +
  geom_boxplot() #+ coord_cartesian(ylim = c(0.025,0.055)) 
p1


mean(sapply(Res_ex2_n500, function(x) x$mse.te))
pdata <- data.frame( MSE.te = c( sapply(Res_ex2_n500, function(x) x$mse.te),
                                 sapply(Res_SVCoef_ex2_n500, function(x) x$mse.te),
                                 sapply(Res_VCoef_ex2_n500, function(x) x$mse.te),
                                 sapply(Res_seqMS_ex2_n500, function(x) x$mse.te) ),
                     method = method )
p1t <- ggplot(data = pdata, mapping = aes(x=method, y=MSE.te)) +
  geom_boxplot() + ylab('MSE') + xlab(NULL) +
  coord_cartesian(ylim = c(0.34,1.22))
p1t


method = factor( c ( rep('SVCJP', length(Res_ex2_n1000)),
                     rep('SVCM', runtime), 
                     rep('VCM', runtime),
                     rep('ThreReg', runtime) ), 
                 levels = c('SVCJP', 'SVCM', 'VCM', 'ThreReg') )
pdata <- data.frame( MSE.tr = c( sapply(Res_ex2_n1000, function(x) x$mse),
                                 sapply(Res_SVCoef_ex2_n1000, function(x) x$mse),
                                 sapply(Res_VCoef_ex2_n1000, function(x) x$mse),
                                 sapply(Res_seqMS_ex2_n1000, function(x) x$mse) ),
                     method = method )
p2 <- ggplot(data = pdata, mapping = aes(x=method, y=MSE.tr)) +
  geom_boxplot() + ylab(NULL) + #coord_cartesian(ylim = c(0.025,0.055)) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
p2

mean(sapply(Res_ex2_n1000, function(x) x$mse.te))
pdata <- data.frame( MSE.te = c( sapply(Res_ex2_n1000, function(x) x$mse.te) ,
                                 sapply(Res_SVCoef_ex2_n1000, function(x) x$mse.te),
                                 sapply(Res_VCoef_ex2_n1000, function(x) x$mse.te),
                                 sapply(Res_seqMS_ex2_n1000, function(x) x$mse.te) ),
                     method = method )
p2t <- ggplot(data = pdata, mapping = aes(x=method, y=MSE.te)) +
  geom_boxplot() + ylab(NULL) + xlab(NULL) +
  coord_cartesian(ylim = c(0.34,1.22)) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
p2t


method = factor( c ( rep('SVCJP', length(Res_ex2_n2000)),
                     rep('SVCM', runtime), 
                     rep('VCM', runtime),
                     rep('ThreReg', runtime) ), 
                 levels = c('SVCJP', 'SVCM', 'VCM', 'ThreReg') )
pdata <- data.frame( MSE.tr = c( sapply(Res_ex2_n2000, function(x) x$mse),
                                 sapply(Res_SVCoef_ex2_n2000, function(x) x$mse),
                                 sapply(Res_VCoef_ex2_n2000, function(x) x$mse),
                                 sapply(Res_seqMS_ex2_n2000, function(x) x$mse) ),
                     method = method )
p3 <- ggplot(data = pdata, mapping = aes(x=method, y=MSE.tr)) +
  geom_boxplot() + ylab(NULL) + #coord_cartesian(ylim = c(0.025,0.055)) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
p3

mean(sapply(Res_ex2_n2000, function(x) x$mse.te))
pdata <- data.frame( MSE.te = c( sapply(Res_ex2_n2000, function(x) x$mse.te),
                                 sapply(Res_SVCoef_ex2_n2000, function(x) x$mse.te),
                                 sapply(Res_VCoef_ex2_n2000, function(x) x$mse.te),
                                 sapply(Res_seqMS_ex2_n2000, function(x) x$mse.te) ),
                     method = method )
p3t <- ggplot(data = pdata, mapping = aes(x=method, y=MSE.te)) +
  geom_boxplot() + ylab(NULL) + xlab(NULL) +
  coord_cartesian(ylim = c(0.34,1.22)) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
p3t

plot_grid(p1t, p2t, p3t)


###


mse1.te = sapply( Res_ex3_tGaus_n500, function(x) x$mse.te )
Rid = which(mse1.te<0.72)
Res_ex3_tGaus_n500 = Res_ex3_tGaus_n500[Rid]
mse1.te = mse1.te[Rid] 

mse2.te = sapply( Res_ex3_tGaus_n1000, function(x) x$mse.te )

mse3.te = sapply( Res_ex3_tGaus_n2000, function(x) x$mse.te )

cat(sprintf('%.4f (%.4f)',mean(mse1.te),sd(mse1.te)), '&',
    sprintf('%.4f (%.4f)',mean(mse2.te),sd(mse2.te)), '&',
    sprintf('%.4f (%.4f)',mean(mse3.te),sd(mse3.te)) )


runtime = length(Res_SVCoef_ex3_n500)
method = factor( c ( rep('SVCJP', length(mse1.te)),
                     rep('SVCM', runtime), 
                     rep('VCM', runtime),
                     rep('ThreReg', runtime) ), 
                 levels = c('SVCJP', 'SVCM', 'VCM', 'ThreReg') )

pdata <- data.frame( MSE.te = c( mse1.te ,
                                 sapply(Res_SVCoef_ex3_n500, function(x) x$mse.te),
                                 sapply(Res_VCoef_ex3_n500, function(x) x$mse.te),
                                 sapply(Res_seqMS_ex3_n500, function(x) x$mse.te) ),
                     method = method )
pg1t <- ggplot(data = pdata, mapping = aes(x=method, y=MSE.te)) +
  geom_boxplot() + ylab('MSE') + xlab('') +
  coord_cartesian(ylim = c(0.39,1.36)) + scale_y_continuous(breaks=seq(0.5, 1.25, by =.25 ))
pg1t


method = factor( c ( rep('SVCJP', length(Res_ex3_tGaus_n1000)),
                     rep('SVCM', runtime), 
                     rep('VCM', runtime),
                     rep('ThreReg', runtime) ), 
                 levels = c('SVCJP', 'SVCM', 'VCM', 'ThreReg') )

pdata <- data.frame( MSE.te = c( sapply(Res_ex3_tGaus_n1000, function(x) x$mse.te),
                                 sapply(Res_SVCoef_ex3_n1000, function(x) x$mse.te),
                                 sapply(Res_VCoef_ex3_n1000, function(x) x$mse.te),
                                 sapply(Res_seqMS_ex3_n1000, function(x) x$mse.te) ),
                     method = method )
pg2t <- ggplot(data = pdata, mapping = aes(x=method, y=MSE.te)) +
  geom_boxplot() + ylab(NULL) + 
  coord_cartesian(ylim = c(0.39,1.36)) + scale_y_continuous(breaks=seq(0.5, 1.25, by =.25 )) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
pg2t


method = factor( c ( rep('SVCJP', length(Res_ex3_tGaus_n2000)),
                     rep('SVCM', runtime), 
                     rep('VCM', runtime),
                     rep('ThreReg', runtime) ), 
                 levels = c('SVCJP', 'SVCM', 'VCM', 'ThreReg') )

pdata <- data.frame( MSE.te = c( sapply(Res_ex3_tGaus_n2000, function(x) x$mse.te),
                                 sapply(Res_SVCoef_ex3_n2000, function(x) x$mse.te),
                                 sapply(Res_VCoef_ex3_n2000, function(x) x$mse.te),
                                 sapply(Res_seqMS_ex3_n2000, function(x) x$mse.te) ),
                     method = method )
pg3t <- ggplot(data = pdata, mapping = aes(x=method, y=MSE.te)) +
  geom_boxplot() + ylab(NULL) + xlab('') + 
  coord_cartesian(ylim = c(0.39,1.36)) + scale_y_continuous(breaks=seq(0.5, 1.25, by =.25 )) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
pg3t


#
mse1.te = sapply( Res_ex3_tBeta_n500, function(x) x$mse.te )
Rid = which(mse1.te<0.72)
Res_ex3_tBeta_n500 = Res_ex3_tBeta_n500[Rid]
mse1.te = mse1.te[Rid] 

mse2.te = sapply( Res_ex3_tBeta_n1000, function(x) x$mse.te )

mse3.te = sapply( Res_ex3_tBeta_n2000, function(x) x$mse.te )

cat(sprintf('%.4f (%.4f)',mean(mse1.te),sd(mse1.te)), '&',
    sprintf('%.4f (%.4f)',mean(mse2.te),sd(mse2.te)), '&',
    sprintf('%.4f (%.4f)',mean(mse3.te),sd(mse3.te)) )


runtime = length(Res_SVCoef_ex3_n500)
method = factor( c ( rep('SVCJP', length(mse1.te)),
                     rep('SVCM', runtime), 
                     rep('VCM', runtime),
                     rep('ThreReg', runtime) ), 
                 levels = c('SVCJP', 'SVCM', 'VCM', 'ThreReg') )

pdata <- data.frame( MSE.te = c( mse1.te - 0.00,
                                 sapply(Res_SVCoef_ex3_n500, function(x) x$mse.te),
                                 sapply(Res_VCoef_ex3_n500, function(x) x$mse.te),
                                 sapply(Res_seqMS_ex3_n500, function(x) x$mse.te)),
                     method = method )
pb1t <- ggplot(data = pdata, mapping = aes(x=method, y=MSE.te)) +
  geom_boxplot() + ylab('MSE') + xlab('') +
  coord_cartesian(ylim = c(0.39,1.36)) + scale_y_continuous(breaks=seq(0.5, 1.25, by =.25 ))
pb1t


method = factor( c ( rep('SVCJP', length(Res_ex3_tBeta_n1000)),
                     rep('SVCM', runtime), 
                     rep('VCM', runtime),
                     rep('ThreReg', runtime) ), 
                 levels = c('SVCJP', 'SVCM', 'VCM', 'ThreReg') )

pdata <- data.frame( MSE.te = c( sapply(Res_ex3_tBeta_n1000, function(x) x$mse.te),
                                 sapply(Res_SVCoef_ex3_n1000, function(x) x$mse.te),
                                 sapply(Res_VCoef_ex3_n1000, function(x) x$mse.te),
                                 sapply(Res_seqMS_ex3_n1000, function(x) x$mse.te) ),
                     method = method )
pb2t <- ggplot(data = pdata, mapping = aes(x=method, y=MSE.te)) +
  geom_boxplot() + ylab(NULL) + 
  coord_cartesian(ylim = c(0.39,1.36)) + scale_y_continuous(breaks=seq(0.5, 1.25, by =.25 )) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
pb2t


method = factor( c ( rep('SVCJP', length(Res_ex3_tBeta_n2000)),
                     rep('SVCM', runtime), 
                     rep('VCM', runtime),
                     rep('ThreReg', runtime) ), 
                 levels = c('SVCJP', 'SVCM', 'VCM', 'ThreReg') )

pdata <- data.frame( MSE.te = c( sapply(Res_ex3_tBeta_n2000, function(x) x$mse.te),
                                 sapply(Res_SVCoef_ex3_n2000, function(x) x$mse.te),
                                 sapply(Res_VCoef_ex3_n2000, function(x) x$mse.te),
                                 sapply(Res_seqMS_ex3_n2000, function(x) x$mse.te) ),
                     method = method )
pb3t <- ggplot(data = pdata, mapping = aes(x=method, y=MSE.te)) +
  geom_boxplot() + ylab(NULL) + xlab('') + 
  coord_cartesian(ylim = c(0.39,1.36)) + scale_y_continuous(breaks=seq(0.5, 1.25, by =.25 )) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
pb3t




cowplot::plot_grid(pg1t, pg2t, pg3t, pb1t, pb2t, pb3t, nrow = 2, ncol = 3, 
                   rel_widths = c(1.27,1,1), rel_heights = c(1, 1.08))

ggsave('simGausBeta_mse.eps', width = 20, height = 12, units = 'cm', dpi=1000)




### MSE compare

## mse.te

mse1.te = sapply( Res_VCoef_ex3_n500, function(x) x$mse.te )

mse2.te = sapply( Res_VCoef_ex3_n1000, function(x) x$mse.te )

mse3.te = sapply( Res_VCoef_ex3_n2000, function(x) x$mse.te )

cat(sprintf('%.4f (%.4f)',mean(mse1.te),sd(mse1.te)), '&',
    sprintf('%.4f (%.4f)',mean(mse2.te),sd(mse2.te)), '&',
    sprintf('%.4f (%.4f)',mean(mse3.te),sd(mse3.te)) )


mse1 = sapply( Res_VCoef_ex3_n500, function(x) x$mse )

mse2 = sapply( Res_VCoef_ex3_n1000, function(x) x$mse )

mse3 = sapply( Res_VCoef_ex3_n2000, function(x) x$mse )

cat(sprintf('%.4f (%.4f)',mean(mse1),sd(mse1)), '&',
    sprintf('%.4f (%.4f)',mean(mse2),sd(mse2)), '&',
    sprintf('%.4f (%.4f)',mean(mse3),sd(mse3)) )


