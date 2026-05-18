
rm(list = ls())
gc()

func_path = '/path/to/your/code/' # './RealData'

source(paste0(func_path, 'infant_function.R'))
source(paste0(func_path, 'compare.R'))


rawdata <- read.csv(file = './original/Dataset_maternal_mental_health_infant_sleep.csv',
                    sep = ';', nrows = 410)

id_missing = which(rawdata$Sleep_night_duration_bb1 == '99:99')
rawdata <- rawdata[-id_missing,]

X <- rawdata[, sapply(names(rawdata), function(x) {
  grepl('Age', x) || grepl('Education', x) || grepl('age', x) || 
    grepl('CBTS', x) || grepl('EPDS', x) || grepl('HADS', x) || grepl('IBQ', x)
})]

ind_CBTS <- which( sapply(names(X), function(x){
  grepl('CBTS', x)
}) )
ind_EPDS <- which( sapply(names(X), function(x){
  grepl('EPDS', x)
}) )
ind_HADS <- which( sapply(names(X), function(x){
  grepl('HADS', x)
}) )
ind_IBQ <- which( sapply(names(X), function(x){
  grepl('IBQ', x)
}) )

X[,ind_IBQ][is.na(X[,ind_IBQ])] = 0

X$EPDS = rowSums(X[,ind_EPDS])
X$HADS = rowSums(X[,ind_HADS])
X$CBTS = rowSums(X[,ind_CBTS])
X$IBQ = rowSums(X[,ind_IBQ])

X = X[,-c(ind_EPDS,ind_HADS,ind_CBTS,ind_IBQ)]

names(X)[grepl('Gestation', names(X))] = 'Gestation'
names(X)[names(X)=='Age'] = 'Age_M'
names(X)[names(X)=='Age_bb'] = 'Age_B'
X$Age_B = X$Age_B * 3

summary(X)

for (j in 1:dim(X)[2]) cat(names(X)[j], round(mean(X[[j]]),3), '&', round(sd(X[[j]]),3), '\n')

X$EPDS = sqrt(X$EPDS)
X$HADS = sqrt(X$HADS)
X$CBTS = sqrt(X$CBTS)

par(mfrow=c(3,3))
for (j in 1:dim(X)[2]) hist(X[,j], main=names(X)[j], breaks=20)

Tmean = mean(X$Age_M)
Tsd = sd(X$Age_M)

X = as.data.frame(scale(X))

Y <- rawdata$Sleep_night_duration_bb1
Y = sapply(Y, function(x) {
  t = unlist(strsplit(x,':'))
  #as.integer(t[1])*60 + as.integer(t[2])
  as.integer(t[1]) + as.integer(t[2])/60
})

# Y <- rawdata$night_awakening_number_bb1
# Y = sqrt(Y)

Y = unname(Y)
Y = c(scale(Y))

hist(Y, breaks = 20)



ZZ = as.data.frame(cbind(Y,X))
lm1 <- lm(Y~., data = ZZ)
summary(lm1)
mse_lm = mean(lm1$residuals^2)
R2_lm = mean((lm1$fitted.values-mean(Y))^2) / var(Y)

for (j in 1:dim(X)[2]) cat(colnames(X)[j],cor(X[,j],Y), '\n')





######################

X = as.matrix(X)

Tid = which(colnames(X) == 'Age_M') 
Xin = X[,-Tid]
Tin = unname(X[,Tid])
Yin = Y

Tmin = min(Tin)
Tmax = max(Tin)
Tin = (Tin - min(Tin))/(max(Tin) - min(Tin))
hist(Tin)


bw.seq1 = seq(0.11, 0.20, by = 0.01)
bw.seq2 = seq(0.25, 0.35, length.out=6)
zeta.seq = NULL 




### Varying coefficient model

res_VCoef <- VCoef_func_2d(tin = Tin, yin = Yin, xin = Xin, wi = NULL,
                           bw.seq = bw.seq2, 
                           nRegGrid = 101, kernel='epan',
                           npoly = 1, nder = 0, kFolds = 5)
par( mfrow=c( 3, ceiling((dim(Xin)[2]+1)/3) ) )
for (j in 1:(dim(Xin)[2]+1)){
  plot(x=res_VCoef$workGrid, y=res_VCoef$alpha_hat_out[,j], type = 'l',
       main = c('intercept', colnames(Xin))[j], ylab = 'alpha', xlab = 'T')
}
mse_VCoef = res_VCoef$mse
R2_VCoef = mean((res_VCoef$yhat-mean(Yin))^2) / var(Yin)

for (j in 1:(dim(Xin)[2]+1)){
  vt = res_VCoef$alpha_hat_out[,j]
  vt1 = vt[-1]
  vt2 = vt[-length(vt)]
  total_var = sum(abs(vt1-vt2))
  
  cat(c('intercept', colnames(Xin))[j], total_var, var(vt), '\n')
}




### Semi varying coefficient model

Zid <- which(colnames(Xin) %in% c('Education') )  # , 'Age_B'

Xin_CP <- Xin[,-Zid, drop=F]
Zin_CP <- Xin[,Zid, drop=F]

res_SVCoef <- SVCoef_func_2d(tin = Tin, yin = Yin, xin = Xin_CP, zin = Zin_CP, wi = NULL,
                             bw.seq1 = bw.seq1, bw.seq2 = bw.seq2, 
                             nRegGrid = 101, kernel='epan',
                             npoly = 1, nder = 0, kFolds = 5, 
                             method_opt='average')
par( mfrow=c( 3, ceiling((dim(Xin_CP)[2]+1)/3) ) )
for (j in 1:(dim(Xin_CP)[2]+1)){
  plot(x=res_SVCoef$workGrid, y=res_SVCoef$alpha_hat_out[,j], type = 'l',
       main = c('intercept', colnames(Xin_CP))[j], ylab = 'alpha', xlab = 'T')
}
# par( mfrow=c( 3, ceiling((dim(Xin_CP)[2]+1)/3) ) )
# for (j in 1:(dim(Xin_CP)[2]+1)){
#   plot(x=res_SVCoef$obsGrid, y=res_SVCoef$alpha_hat[,j], type = 'l',
#        main = c('intercept', colnames(Xin_CP))[j], ylab = 'alpha', xlab = 'T')
# }
mse_SVCoef = res_SVCoef$mse
R2_SVCoef = mean((res_SVCoef$yhat-mean(Yin))^2) / var(Yin)



### Our method

rangeT = c(quantile(Tin,0.2), quantile(Tin,0.8))
hist(Tin, breaks = 20)
  

resCP = VCMCP_2d(tin = Tin, yin = Yin, xin = Xin_CP, zin = Zin_CP, wi = NULL, 
                 zeta = zeta.seq, cutoff = max, alpha = 0.05, kernel='epan',
                 bw.seq1 = bw.seq1, bw.seq2 = bw.seq2, 
                 NbGrid = 101, nRegGrid = 101, kFolds = 5, 
                 npoly = 1, nder = 0, refined = TRUE,
                 hkappa = 1.2, h2_sepcv = F,
                 method_opt='average')


save(resCP, file = 'resCP_plot.RData')










##################### Boostrap

source(paste0('/home/pan/MMVCM/infant/', 'VCMCP_2d_real.R'))

# load("/home/pan/MMVCM/infant/resCP_plot.RData")
load("/home/pan/MMVCM/infant/res_try.RData")


{
  nm_order1 = 15
  nm_order2 = 16
  
  bt_time = 1000
  ifreplace = F
  
  
  resCP = res_try  # resCP_plot_final
  
  tin <- resCP$tin
  yin <- resCP$yin
  xin <- resCP$xin
  zin <- resCP$zin
  
  n = dim(xin)[1]
  p = dim(xin)[2]
  
  h_1 = resCP$h_1
  h_tau = resCP$h_tau
  h_d = resCP$h_d
  h_2 = resCP$h_2
  zeta = resCP$zeta
  
  sapply(resCP$alp_est, function(x) x[[1]]$jumptime)
  sapply(resCP$alp_est, function(x) x[[2]]$jumptime)
  
  ### reproduce the full sample estimation
  
  res = VCMCP_2d_real_bw(tin = tin, yin = yin, xin = xin, zin = zin, 
                         h_1 = h_1, h_tau = h_tau, h_d = h_d, h_2 = h_2, zeta = zeta,
                         wi = NULL, 
                         cutoff = max, alpha = 0.05, kernel='epan',
                         NbGrid = 101, nRegGrid = 101, kFolds = 5, 
                         npoly = 1, nder = 0, refined = TRUE,
                         hkappa = 1.2, h2_sepcv = F,
                         method_opt='average')
  
  
  #
  nm = round(n ^ (nm_order1/nm_order2))
  
  h_1_BT = h_1 * (log(nm) / nm) ^ (1/4) / (log(n) / n) ^ (1/4)   # nm ^ (-0.175) / n ^ (-0.175)
  h_tau_BT = h_tau * (log(nm) / nm) ^ (1/5) / (log(n) / n) ^ (1/5)   # nm ^ (-0.225) / n ^ (-0.225)
  h_d_BT = h_d * (log(nm) / nm) ^ (1/5) / (log(n) / n) ^ (1/5)   # nm ^ (-0.175) / n ^ (-0.175)
  h_2_BT = h_tau_BT  # h_2 * (log(nm) / nm) ^ (1/3) / (log(n) / n) ^ (1/3)
  zeta_BT = zeta * ((h_tau_BT + sqrt(log(nm)/(nm*h_tau_BT)) + h_1_BT + sqrt(log(nm)/(nm*h_1_BT))) / 
                      (h_tau + sqrt(log(n)/(n*h_tau)) + h_1 + sqrt(log(n)/(n*h_1)))) ^ (1/2)
  
  bt_res = NULL
  
  bt_round = 1
  while (bt_round <= bt_time){
    
    cat('\n')
    cat('########### The NO.', bt_round, ' Boostrap begins...\n')
    
    bt_id = sort(sample(1:n, nm, replace = ifreplace))
    
    Xm = xin[bt_id,,drop = F]
    Ym = yin[bt_id]
    Zm = zin[bt_id,,drop = F]
    Tm = tin[bt_id,,drop = F]
    
    ###
    # bw.seq1 = seq(0.10, 0.15, by = 0.01)
    # bw.seq2 = seq(0.20, 0.35, by = 0.03)
    # 
    # resBT = VCMCP_2d(tin = Tm, yin = Ym, xin = Xm, zin = Zm, wi = NULL, 
    #                  zeta = NULL, cutoff = max, alpha = 0.05, kernel='epan',
    #                  bw.seq1 = bw.seq1, bw.seq2 = bw.seq2, 
    #                  NbGrid = 51, nRegGrid = 101, kFolds = 5, 
    #                  npoly = 1, nder = 0, refined = TRUE,
    #                  hkappa = hkappa, method_opt='average',
    #                  h2_sepcv = F)
    ###
    
    bt_res[[bt_round]] = try(VCMCP_2d_real_bw(tin = Tm, yin = Ym, xin = Xm, zin = Zm, 
                                              h_1 = h_1_BT, h_tau = h_tau_BT, h_d = h_d_BT, h_2 = h_2_BT, zeta = zeta_BT,
                                              wi = NULL,
                                              cutoff = max, alpha = 0.05, kernel='epan',
                                              NbGrid = 101, nRegGrid = 101, kFolds = 5, 
                                              npoly = 1, nder = 0, refined = TRUE,
                                              hkappa = 1.2, h2_sepcv = F,
                                              method_opt='average'))
    if (!is.character(bt_res[[bt_round]])) bt_round = bt_round + 1
    
  }
  
  
  resCP$BT_res = bt_res
  
  
  save(resCP, file = paste0('BTres_infant_mord', nm_order1, nm_order2, '.RData'))
  
}












############### Visualization

load("resCP_plot.RData")

##
resCP = resCP_plot_final
f_agem = ecdf(X_org$Age_M)
f_gest = ecdf(X_org$Gestation)
resCP$workGrid_org = cbind(approx(x=f_agem(X_org$Age_M), y=X_org$Age_M, xout = resCP$workGrid[,1])$y,
                           approx(x=f_gest(X_org$Gestation), y=X_org$Gestation, xout = resCP$workGrid[,2])$y)

approx(x=f_agem(X_org$Age_M), y=X_org$Age_M, xout = resCP$alp_est$alpha2[[1]]$jumptime)$y
approx(x=f_gest(X_org$Gestation), y=X_org$Gestation, xout = resCP$alp_est$alpha0[[2]]$jumptime)$y

library('gstat')
library('sp')
cord_idw = resCP$workGrid_org
alp_idw = sapply(resCP$alp_est, function(x) x$alp.hat_out)
grd_idw = expand.grid(t1 = seq(min(cord_idw[,1]), max(cord_idw[,1]), length.out=201),
                      t2 = seq(min(cord_idw[,2]), max(cord_idw[,2]), length.out=201))
resCP$grid_plot = as.matrix(grd_idw)
coordinates(grd_idw) = c('t1', 't2')
alp_plot = c()
for (j in 1:dim(alp_idw)[2]) {
  dat_idw = cbind(cord_idw, alp_idw[,j])
  dat_idw = as.data.frame(dat_idw)
  colnames(dat_idw) = c('t1', 't2', 'alpha')
  coordinates(dat_idw) = c('t1', 't2')
  new_alp = idw(alpha ~ 1, dat_idw, grd_idw, nmax=5, idp=20)
  alp_plot = cbind(alp_plot, new_alp$var1.pred)
}
colnames(alp_plot) = paste0('alpha', 1:dim(alp_idw)[2]-1)
resCP$alp_plot = alp_plot


# pdat <- as.data.frame(cbind(resCP$grid_plot, resCP$alp_plot))
pdat <- data.frame(t1=resCP$workGrid_org[,1], t2=resCP$workGrid_org[,2], 
                   alpha0=resCP$alp_est$alpha0$alp.hat_out, 
                   alpha1=resCP$alp_est$alpha1$alp.hat_out, 
                   alpha2=resCP$alp_est$alpha2$alp.hat_out,
                   alpha3=resCP$alp_est$alpha3$alp.hat_out, 
                   alpha4=resCP$alp_est$alpha4$alp.hat_out,
                   alpha5=resCP$alp_est$alpha5$alp.hat_out)
plot0 <- lattice::cloud(alpha0~t1*t2, data = pdat, zlab=list(NULL, rot=90), 
                        par.settings = list(
                          axis.line = list(col = "transparent"),
                          fontsize = list(text = 12, points = 10)
                        ),
                        scales = list(arrows=F), main = 'intercept',
                        zlim = c(-1,1), screen = list(z = 30, x = -65, y = 0)) 
plot0
plot1 <- lattice::cloud(alpha1~t1*t2, data = pdat, zlab=list(NULL, rot=90), 
                        par.settings = list(
                          axis.line = list(col = "transparent"),
                          fontsize = list(text = 12, points = 10)
                        ),
                        scales = list(arrows=F), main = colnames(Xin_CP)[1],
                        zlim = c(-1,1), screen = list(z = 30, x = -65, y = 0)) 
plot1
plot2 <- lattice::cloud(alpha2~t1*t2, data = pdat, zlab=list(NULL, rot=90), 
                        par.settings = list(
                          axis.line = list(col = "transparent"),
                          fontsize = list(text = 12, points = 10)
                        ),
                        scales = list(arrows=F), main = colnames(Xin_CP)[2],
                        zlim = c(-1,1), screen = list(z = 30, x = -65, y = 0)) 
plot2
plot3 <- lattice::cloud(alpha3~t1*t2, data = pdat, zlab=list(NULL, rot=90), 
                        par.settings = list(
                          axis.line = list(col = "transparent"),
                          fontsize = list(text = 12, points = 10)
                        ),
                        scales = list(arrows=F), main = colnames(Xin_CP)[3],
                        zlim = c(-1,1), screen = list(z = 30, x = -65, y = 0)) 
plot3
plot4 <- lattice::cloud(alpha4~t1*t2, data = pdat, zlab=list(NULL, rot=90), 
                        par.settings = list(
                          axis.line = list(col = "transparent"),
                          fontsize = list(text = 12, points = 10)
                        ),
                        scales = list(arrows=F), main = colnames(Xin_CP)[4],
                        zlim = c(-1,1), screen = list(z = 30, x = -65, y = 0)) 
plot4
plot5 <- lattice::cloud(alpha5~t1*t2, data = pdat, zlab=list(NULL, rot=90), 
                        par.settings = list(
                          axis.line = list(col = "transparent"),
                          fontsize = list(text = 12, points = 10)
                        ),
                        scales = list(arrows=F), main = colnames(Xin_CP)[5],
                        zlim = c(-1,1), screen = list(z = 30, x = -65, y = 0)) 
plot5

cowplot::plot_grid(plot2, plot3, plot4, nrow = 1) 
ggsave('1.eps', last_plot(), dpi = 1000,
       units = 'cm', width = 40, height = 15)


library(ggplot2)
library(cowplot)

pdat <- as.data.frame(cbind(resCP$grid_plot, resCP$alp_plot))
# pdat <- data.frame(t1=resCP$workGrid_org[,1], t2=resCP$workGrid_org[,2], 
#                    alpha0=resCP$alp_est$alpha0$alp.hat_out, 
#                    alpha1=resCP$alp_est$alpha1$alp.hat_out, 
#                    alpha2=resCP$alp_est$alpha2$alp.hat_out,
#                    alpha3=resCP$alp_est$alpha3$alp.hat_out, 
#                    alpha4=resCP$alp_est$alpha4$alp.hat_out,
#                    alpha5=resCP$alp_est$alpha5$alp.hat_out)
colmin = min(pdat[,-c(1,2)])
colmax = max(pdat[,-c(1,2)])
pid = which(!is.na(pdat$alpha0) & abs(pdat$alpha3)<10)
pdat = pdat[pid,]

orig_data <- data.frame(t1 = X_org$Age_M, t2 = X_org$Gestation)

my_breaks <- seq(min(pdat$t1, na.rm = TRUE), max(pdat$t1, na.rm = TRUE), length.out = 7)
hist_bottom <- ggplot(orig_data, aes(x=t1)) +
  # geom_histogram(fill = "grey60", color = "white", bins = 30) +
  geom_histogram(breaks = my_breaks, fill = "grey60", color = "white") +
  coord_cartesian(xlim = c(26.85, 32)) + # MUST match heatmap xlim
  theme_void() + # Removes all axes, background, and text from the histogram
  theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0))
hist_bottom

my_breaks2 <- seq(min(pdat$t2, na.rm = TRUE), max(pdat$t2, na.rm = TRUE), length.out = 6)
hist_left <- ggplot(orig_data, aes(y=t2)) +
  # geom_histogram(fill = "grey60", color = "white", bins = 30) +
  geom_histogram(breaks = my_breaks2, fill = "grey60", color = "white") +
  coord_cartesian(ylim = c(37.96, 40.3)) + # MUST match heatmap xlim
  theme_void() + # Removes all axes, background, and text from the histogram
  theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0))
hist_left

heat0 <- ggplot(pdat, aes(x=t1, y=t2)) + geom_point(aes(color=alpha0), size=1) +
  #  scale_color_gradientn(colors=c('red', 'yellow', 'green','blue')) +
  scale_color_viridis_c(name=expression(alpha), direction = 1, option = 'H', limits = c(colmin, colmax)) +
  # geom_rug(data = orig_data, aes(x=t1), inherit.aes = FALSE, sides = "b", color = "grey60", alpha = 0.3, linewidth = 40) +
  xlab(expression(t[1])) + ylab(expression(t[2])) +
  annotate('text', x=26.65, y=39.97, label = '39.97', col='red', size=2.3) + 
  # annotate('text', x=27.9, y=37.8, label = '27.9', col='red', size=2.3) + 
  coord_cartesian(ylim=c(37.96, 40.3), xlim = c(26.85, 32), clip='off') +
  ggtitle('intercept') + theme_classic() + theme(plot.title = element_text(hjust = 0.5)) 
heat0

heat2_with_hist <- plot_grid(
  heat0,          # Your existing heatmap
  hist_bottom,    # The new clean histogram
  ncol = 1,       # Stack them vertically
  align = "v",    # Vertically align the plot panels
  axis = "lr",    # Lock the left and right axes together
  rel_heights = c(1, 0.15) # The heatmap takes up 85% of the height, histogram 15%
)

layout <- "
AB
#C
"
final_heat0 <- hist_left + heat0 + hist_bottom + 
  plot_layout(
    design = layout, 
    widths = c(1.5, 10),  # hist_left gets 1.5 parts width, heatmap gets 10
    heights = c(10, 1.5)  # heatmap gets 10 parts height, hist_bottom gets 1.5
  )

final_heat0


heat1 <- ggplot(pdat, aes(x=t1, y=t2)) + geom_point(aes(color=alpha1), size=1) +
  scale_color_viridis_c(name=expression(alpha), direction = 1, option = 'H', limits = c(colmin, colmax)) +
  xlab(expression(t[1])) + ylab(expression(t[2])) +
  # annotate('text', x=26.65, y=39.97, label = '39.97', col='red', size=2.3) + 
  # annotate('text', x=27.9, y=37.8, label = '27.9', col='red', size=2.3) + 
  coord_cartesian(ylim=c(37.96, 40.3), xlim = c(26.85, 32), clip='off') +
  ggtitle(colnames(Xin_CP)[1]) + theme_classic() + theme(plot.title = element_text(hjust = 0.5)) 
heat1
heat2 <- ggplot(pdat, aes(x=t1, y=t2)) + geom_point(aes(color=alpha2), size=1) +
  scale_color_viridis_c(name=expression(alpha), direction = 1, option = 'H', limits = c(colmin, colmax)) +
  xlab(expression(t[1])) + ylab(expression(t[2])) +
  # annotate('text', x=26.65, y=39.97, label = '39.97', col='red', size=2.3) + 
  annotate('text', x=27.98, y=37.9, label = '27.98', col='red', size=2.3) +
  annotate('text', x=29.76, y=37.9, label = '29.76', col='red', size=2.3) +
  coord_cartesian(ylim=c(37.96, 40.3), xlim = c(26.85, 32), clip='off') +
  ggtitle(colnames(Xin_CP)[2]) + theme_classic() + theme(plot.title = element_text(hjust = 0.5)) 
heat2
heat3 <- ggplot(pdat, aes(x=t1, y=t2)) + geom_point(aes(color=alpha3), size=1) +
  scale_color_viridis_c(name=expression(alpha), direction = 1, option = 'H', limits = c(colmin, colmax)) +
  xlab(expression(t[1])) + ylab(expression(t[2])) +
  annotate('text', x=26.65, y=39.97, label = '39.97', col='red', size=2.3) + 
  annotate('text', x=26.65, y=39.37, label = '39.37', col='red', size=2.3) + 
  annotate('text', x=27.77, y=37.9, label = '27.77', col='red', size=2.3) +
  coord_cartesian(ylim=c(37.96, 40.3), xlim = c(26.85, 32), clip='off') +
  ggtitle(colnames(Xin_CP)[3]) + theme_classic() + theme(plot.title = element_text(hjust = 0.5)) 
heat3
heat4 <- ggplot(pdat, aes(x=t1, y=t2)) + geom_point(aes(color=alpha4), size=1) +
  scale_color_viridis_c(name=expression(alpha), direction = 1, option = 'H', limits = c(colmin, colmax)) +
  xlab(expression(t[1])) + ylab(expression(t[2])) +
  annotate('text', x=26.65, y=39.37, label = '39.37', col='red', size=2.3) + 
  annotate('text', x=27.98, y=37.9, label = '27.98', col='red', size=2.3) +
  annotate('text', x=29.76, y=37.9, label = '29.76', col='red', size=2.3) +
  coord_cartesian(ylim=c(37.96, 40.3), xlim = c(26.85, 32), clip='off') +
  ggtitle(colnames(Xin_CP)[4]) + theme_classic() + theme(plot.title = element_text(hjust = 0.5)) 
heat4
heat5 <- ggplot(pdat, aes(x=t1, y=t2)) + geom_point(aes(color=alpha5), size=1) +
  scale_color_viridis_c(name=expression(alpha), direction = 1, option = 'H', limits = c(colmin, colmax)) +
  ggtitle(colnames(Xin_CP)[5]) + theme_classic() + theme(plot.title = element_text(hjust = 0.5)) 
heat5

library('cowplot')
comb = plot_grid(heat0+ theme(legend.position = 'none'),
                 heat1+ theme(legend.position = 'none'), 
                 heat2+ theme(legend.position = 'none'),
                 heat3+ theme(legend.position = 'none'), 
                 heat4+ theme(legend.position = 'none'), 
                 heat5+ theme(legend.position = 'none'), nrow = 2) 
plot_grid(comb, get_legend(heat0), rel_widths = c(12,1)) 
ggsave('2.eps', last_plot(), dpi = 1000,
       units = 'cm', width = 42, height = 26)


my_design <- "
ABCD
EFGH
#IJK
"
final_grid <- hist_left + heat0 + heat1 + heat2 + 
  hist_left + heat3 + heat4 + heat5 + 
  hist_bottom + hist_bottom + hist_bottom +
  
  plot_layout(
    design = my_design,
    widths = c(1, 10, 10, 10), 
    heights = c(10, 10, 1),
    
    guides = "collect" 
  )

final_grid

ggsave('2.pdf', last_plot(), dpi = 1000,
       units = 'cm', width = 42, height = 26)


###

library('ggplot2')
library('cowplot')

grid11=which(resCP$obsGrid[,1]<=0.2960)
grid12=which(resCP$obsGrid[,1]>0.2960 & resCP$obsGrid[,1]<=0.5068)
grid13=which(resCP$obsGrid[,1]>0.5068)

xx = cbind(resCP$zin, resCP$xin)

drawdata <- data.frame(covariate=rep(colnames(xx),3), 
                       group=c(rep('below 28 years old (112)',dim(xx)[2]),
                               rep('28~29 years old (66)',dim(xx)[2]),
                               rep('above 29 years old (231)',dim(xx)[2])), 
                       mean=c(apply(xx[grid11,],2,mean),
                              apply(xx[grid12,],2,mean),
                              apply(xx[grid13,],2,mean) ))
drawdata$covariate <- factor(drawdata$covariate, levels = colnames(xx))
p=ggplot(data=drawdata, aes(x=covariate,y=mean,fill=group))
bar1 = p+geom_col(position = 'dodge')+
  theme(axis.text.x = element_text(angle = 60,vjust=1,hjust=1))+
  scale_fill_manual(name = 'group (by maternal age)',
                    values = c('skyblue2','green','yellow')) +
  theme(legend.title = element_text(size=9),
        legend.key.width = unit(0.4, "cm"),
        legend.key.height = unit(0.4, "cm"),
        legend.text = element_text(size = 8))
bar1





grid21=which(resCP$obsGrid[,2]<=0.5612)
grid22=which(resCP$obsGrid[,2]>0.5612 & resCP$obsGrid[,2]<=0.7652)
grid23=which(resCP$obsGrid[,2]>0.7652)


xx = cbind(resCP$zin, resCP$xin)

drawdata <- data.frame(covariate=rep(colnames(xx),3), 
                       group=c(rep('below 39.4 weeks (225)',dim(xx)[2]),
                               rep('39.4~40 weeks (77)',dim(xx)[2]),
                               rep('above 40 weeks (107)',dim(xx)[2])), 
                       mean=c(apply(xx[grid21,],2,mean),
                              apply(xx[grid22,],2,mean),
                              apply(xx[grid23,],2,mean) ))
drawdata$covariate <- factor(drawdata$covariate, levels = colnames(xx))
p=ggplot(data=drawdata, aes(x=covariate,y=mean,fill=group))
bar2 = p+geom_col(position = 'dodge')+
  theme(axis.text.x = element_text(angle = 60,vjust=1,hjust=1))+
  scale_fill_manual(name = 'group (by gestation length)',
                    values = c('pink2','orange','purple1')) +
  theme(legend.title = element_text(size=9),
        legend.key.width = unit(0.4, "cm"),
        legend.key.height = unit(0.4, "cm"),
        legend.text = element_text(size = 8))
bar2

library('patchwork')
bar1+bar2 + plot_layout(guides = 'collect')

# plot_grid(bar1+theme(legend.position = 'none'),
#           bar2+theme(legend.position = 'none'),
#           plot_grid(get_legend(bar1),get_legend(bar2), ncol=1), 
#           nrow = 1, rel_widths = c(1,1,0.5))


ggsave('3.eps', width = 26, height = 12, units = 'cm', dpi=1000)



###############


f_agem = ecdf(X_org$Age_M)
f_gest = ecdf(X_org$Gestation)
res_SVCoef$workGrid_org = cbind(approx(x=f_agem(X_org$Age_M), y=X_org$Age_M, xout = res_SVCoef$workGrid[,1])$y,
                                approx(x=f_gest(X_org$Gestation), y=X_org$Gestation, xout = res_SVCoef$workGrid[,2])$y)
# a=approx(x=f_agem(X_org$Age_M), y=X_org$Age_M, xout = res_SVCoef$tin[,1])$y
# aa = cbind(X_org$Age_M[order(X_org$Age_M)], a)


library('gstat')
library('sp')
cord_idw = res_SVCoef$workGrid_org
alp_idw = res_SVCoef$alpha_hat_out
grd_idw = as.data.frame(resCP$grid_plot) 
res_SVCoef$grid_plot = resCP$grid_plot
coordinates(grd_idw) = c('t1', 't2')
alp_plot = c()
for (j in 1:dim(alp_idw)[2]) {
  dat_idw = cbind(cord_idw, alp_idw[,j])
  dat_idw = as.data.frame(dat_idw)
  colnames(dat_idw) = c('t1', 't2', 'alpha')
  coordinates(dat_idw) = c('t1', 't2')
  new_alp = idw(alpha ~ 1, dat_idw, grd_idw, nmax=5, idp=20)
  alp_plot = cbind(alp_plot, new_alp$var1.pred)
}
colnames(alp_plot) = paste0('alpha', 1:dim(alp_idw)[2]-1)
res_SVCoef$alp_plot = alp_plot


pdat <- as.data.frame(cbind(res_SVCoef$grid_plot, res_SVCoef$alp_plot))

colmin = min(resCP$alp_plot)
colmax = max(resCP$alp_plot)
pid = which(!is.na(pdat$alpha0) & abs(pdat$alpha3)<10)
pdat = pdat[pid,]
pdat$alpha2[pdat$alpha2<colmin] = colmin
heat0 <- ggplot(pdat, aes(x=t1, y=t2)) + geom_point(aes(color=alpha0), size=1) +
  #  scale_color_gradientn(colors=c('red', 'yellow', 'green','blue')) +
  scale_color_viridis_c(name=expression(alpha), direction = 1, option = 'H', limits = c(colmin, colmax)) +
  xlab(expression(t[1])) + ylab(expression(t[2])) +
  coord_cartesian(ylim=c(37.96, 40.3), xlim = c(26.85, 32), clip='off') +
  ggtitle('intercept') + theme_classic() + theme(plot.title = element_text(hjust = 0.5)) 
heat0
heat1 <- ggplot(pdat, aes(x=t1, y=t2)) + geom_point(aes(color=alpha1), size=1) +
  scale_color_viridis_c(direction = 1, option = 'H', limits = c(colmin, colmax)) +
  xlab(expression(t[1])) + ylab(expression(t[2])) +
  coord_cartesian(ylim=c(37.96, 40.3), xlim = c(26.85, 32), clip='off') +
  ggtitle(colnames(Xin_CP)[1]) + theme_classic() + theme(plot.title = element_text(hjust = 0.5)) 
heat1
heat2 <- ggplot(pdat, aes(x=t1, y=t2)) + geom_point(aes(color=alpha2), size=1) +
  scale_color_viridis_c(direction = 1, option = 'H', limits = c(colmin, colmax)) +
  xlab(expression(t[1])) + ylab(expression(t[2])) +
  coord_cartesian(ylim=c(37.96, 40.3), xlim = c(26.85, 32), clip='off') +
  ggtitle(colnames(Xin_CP)[2]) + theme_classic() + theme(plot.title = element_text(hjust = 0.5)) 
heat2
heat3 <- ggplot(pdat, aes(x=t1, y=t2)) + geom_point(aes(color=alpha3), size=1) +
  scale_color_viridis_c(direction = 1, option = 'H', limits = c(colmin, colmax)) +
  xlab(expression(t[1])) + ylab(expression(t[2])) +
  coord_cartesian(ylim=c(37.96, 40.3), xlim = c(26.85, 32), clip='off') +
  ggtitle(colnames(Xin_CP)[3]) + theme_classic() + theme(plot.title = element_text(hjust = 0.5)) 
heat3
heat4 <- ggplot(pdat, aes(x=t1, y=t2)) + geom_point(aes(color=alpha4), size=1) +
  scale_color_viridis_c(direction = 1, option = 'H', limits = c(colmin, colmax)) +
  xlab(expression(t[1])) + ylab(expression(t[2])) +
  coord_cartesian(ylim=c(37.96, 40.3), xlim = c(26.85, 32), clip='off') +
  ggtitle(colnames(Xin_CP)[4]) + theme_classic() + theme(plot.title = element_text(hjust = 0.5)) 
heat4
heat5 <- ggplot(pdat, aes(x=t1, y=t2)) + geom_point(aes(color=alpha5), size=1) +
  scale_color_viridis_c(direction = 1, option = 'H', limits = c(colmin, colmax)) +
  ggtitle(colnames(Xin_CP)[5]) + theme_classic() + theme(plot.title = element_text(hjust = 0.5)) 
heat5

library('cowplot')
comb = plot_grid(heat0+ theme(legend.position = 'none'),
                 heat1+ theme(legend.position = 'none'), 
                 heat2+ theme(legend.position = 'none'),
                 heat3+ theme(legend.position = 'none'), 
                 heat4+ theme(legend.position = 'none'), 
                 heat5+ theme(legend.position = 'none'), nrow = 2) 
plot_grid(comb, get_legend(heat0), rel_widths = c(12,1)) 
ggsave('2_svcoef.eps', last_plot(), dpi = 1000,
       units = 'cm', width = 42, height = 26)



