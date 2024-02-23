library(rstan)
rstan_options(auto_write = TRUE)
library(doParallel)
library(GGally)
library(ggpubr)
library(ggplot2)
library(ggjoy)
library(bayesplot)
library(dplyr)
library(tidyr)

#prob eff to simulate
pe = c(0.38,0.55,0.71,0.83)
#prob tox to simulate 
pt = c(0.05,0.12,0.27,0.5)

psi=1

#1=N,2=T,3=E,4=B
#fgm
fgm = function(p_e, p_t, psi){
  vec = (1-p_e) * (1-p_t) + psi * (1-p_e) * (1-p_t) * p_e * p_t
  vec[2]=  (1-p_e) - vec[1] 
  vec[3]=  (1-p_t) - vec[1] 
  vec[4] = 1 - (1-p_e) - (1-p_t) + vec[1] 
  return(vec)
}  

doses = c(rep(1:2,each=6),rep(2:4,each=9))
nsim=5000
npats=30
num_doses = 4


dat = matrix(0, ncol=npats,nrow=nsim)
set.seed=123;
for (i in 1:nsim) {
dat[i,1:6] = sample(1:4, 6, replace=TRUE, prob=fgm(pe[1],pt[1],psi))
dat[i,7:12] = sample(1:4, 6, replace=TRUE, prob=fgm(pe[2],pt[2],psi))
dat[i,13:21] =  sample(1:4, 9, replace=TRUE, prob=fgm(pe[3],pt[3],psi))
dat[i,22:30] =  sample(1:4, 9, replace=TRUE, prob=fgm(pe[4],pt[4],psi))
}
eff= matrix(1L* dat %in% c(3,4), ncol=npats)
tox= matrix(1L* dat %in% c(2,4), ncol=npats)
n=c(6,6,9,9)
x_e = x_t =  matrix(0,nrow=nsim, ncol=num_doses)

x_e[,1] = rowSums(eff[,1:6])
x_e[,2] = rowSums(eff[,7:12])
x_e[,3] = rowSums(eff[,13:21])
x_e[,4] = rowSums(eff[,22:30])
x_t[,1] = rowSums(tox[,1:6])
x_t[,2] = rowSums(tox[,7:12])
x_t[,3] = rowSums(tox[,13:21])
x_t[,4] = rowSums(tox[,22:30])

#fit FGM cop
setwd("C:/Users/medahala/OneDrive - University of Leeds/PhD/R Code/copula_sim")
# 
# FGM_stan = stan_model(file= 'EffTox.stan')
# ind_stan = stan_model(file= 'EffTox_independent_vector.stan')
# 
# 
# one_model_fit_and_save_IND = function(i, dat_xe, dat_xt) {
#   saveRDS( 
#     suppressWarnings(rstan::extract(rstan::sampling(ind_stan, 
#                                                     data=list(num_doses = 4,
#                                                               X=c(0,1,2,3), X2=c(0,1,4,9), 
#                                                               n = c(6,6,9,9),
#                                                               alpha_mean = -3, alpha_sd = 3,
#                                                               beta_mean = 1, beta_sd = 2,
#                                                               gamma_mean = -1, gamma_sd = 3,
#                                                               zeta_mean = 1, zeta_sd =2,
#                                                               eta_mean = 0, eta_sd = 0.25, 
#                                                               x_e = dat_xe[i,], x_t=dat_xt[i,]), 
#                                                     iter= 3000,
#                                                     warmup=1000,
#                                                     chains=3L,
#                                                     cores = 1L,
#                                                     seed = "12758694",
#                                                     refresh=0), pars=c("alpha","beta","gamma","zeta","eta")))
#     ,file=paste0("./stan_fits/IND_Log", i,".rds"), compress=FALSE
#   )
# }
# 
# ncores=parallel::detectCores()
# #fit all stan models for efficacy
# cl = makeCluster(ncores)
# registerDoParallel(cl)
# garb = foreach(i = 1:nsim, .inorder=FALSE) %dopar% {
#   one_model_fit_and_save_IND(i, x_e, x_t)
# }
# stopCluster(cl)

##############################################################
########fgm
##############################################################
# 
# one_model_fit_and_save_FGM = function(i, eff_dat, tox_dat) {
#   saveRDS( 
#     suppressWarnings(rstan::extract(rstan::sampling(FGM_stan, 
#                                                     data=list(num_doses = 4, num_patients=30,
#                                                               coded_doses=c(0,1,2,3), coded_doses_squ=c(0,1,4,9), 
#                                                               alpha_mean = -3, alpha_sd = 3,
#                                                               beta_mean = 1, beta_sd = 2,
#                                                               gamma_mean = -1, gamma_sd = 3,
#                                                               zeta_mean = 1, zeta_sd =2,
#                                                               eta_mean = 0, eta_sd = 0.25,
#                                                               psi_l=-1, psi_u=1,
#                                                               doses = c(rep(1,6),rep(2,6),rep(3,9),rep(4,9)),
#                                                               eff = eff_dat[i,], tox=tox_dat[i,]), 
#                                                     iter= 3000,
#                                                     warmup=1000,
#                                                     chains=3L,
#                                                     cores = 1L,
#                                                     seed = "12758694",
#                                                     refresh=0),pars=c("alpha","beta","gamma","zeta","eta","psi")))
#     ,file=paste0("./stan_fits/FGM_Log", i,".rds"), compress=FALSE
#   )
# }
# 
# ncores=parallel::detectCores()
# #fit all stan models for efficacy
# cl = makeCluster(ncores)
# registerDoParallel(cl)
# garb = foreach(i = 1:nsim, .inorder=FALSE) %dopar% {
#   one_model_fit_and_save_FGM(i, eff, tox)
# }
# stopCluster(cl)
# 
# 

m_FGM_alpha=m_FGM_beta=m_FGM_gamma=m_FGM_zeta=
sd_FGM_alpha=sd_FGM_beta=sd_FGM_gamma=sd_FGM_zeta=
m_IND_alpha=m_IND_beta=m_IND_gamma=m_IND_zeta=
sd_IND_alpha=sd_IND_beta=sd_IND_gamma=sd_IND_zeta=rep(NA,nsim)


for (i in 1:nsim){
  dat=readRDS(file=paste0("./stan_fits/FGM_Log", i,".rds"))

  m_FGM_alpha[i]= mean(dat$alpha)
  m_FGM_beta[i]= mean(dat$beta)
  m_FGM_gamma[i]= mean(dat$gamma)
  m_FGM_zeta[i]= mean(dat$zeta)
  sd_FGM_alpha[i]= sd(dat$alpha)
  sd_FGM_beta[i]= sd(dat$beta)
  sd_FGM_gamma[i]= sd(dat$gamma)
  sd_FGM_zeta[i]= sd(dat$zeta)
  
  dat=readRDS(file=paste0("./stan_fits/IND_Log", i,".rds"))
  
  m_IND_alpha[i]= mean(dat$alpha)
  m_IND_beta[i]= mean(dat$beta)
  m_IND_gamma[i]= mean(dat$gamma)
  m_IND_zeta[i]= mean(dat$zeta)
  sd_IND_alpha[i]= sd(dat$alpha)
  sd_IND_beta[i]= sd(dat$beta)
  sd_IND_gamma[i]= sd(dat$gamma)
  sd_IND_zeta[i]= sd(dat$zeta)
}


BIG =        tibble(data.frame(sm_alpha =   m_FGM_alpha -  m_IND_alpha,
                               sm_beta  =   m_FGM_beta  -  m_IND_beta,
                               sm_gamma =   m_FGM_gamma -  m_IND_gamma,
                               sm_zeta =   m_FGM_zeta -  m_IND_zeta,
                               ssd_alpha =   sd_FGM_alpha/sd_IND_alpha,
                               ssd_beta =   sd_FGM_beta/sd_IND_beta,
                               ssd_gamma =   sd_FGM_gamma/sd_IND_gamma,
                               ssd_zeta =   sd_FGM_zeta/sd_IND_zeta))
                               
                               
#alpha
p1<-ggplot(BIG, aes(x=sm_alpha)) +
  geom_density() +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank() 
  ) + xlab("Mean Difference") + ylab("Density")
p2<-ggplot(BIG, aes(x=ssd_alpha)) +
  geom_density() +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank() 
  ) + xlab("Ratio of Standard deviation") + ylab("Density")


#beta
p3<-ggplot(BIG, aes(x=sm_beta)) +
  geom_density() +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank() 
  ) + xlab("Mean Difference") + ylab("Density")
p4<-ggplot(BIG, aes(x=ssd_beta)) +
  geom_density() +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank() 
  ) + xlab("Ratio of Standard deviation") + ylab("Density")

a=ggarrange(p1,p2,p3,p4)

output_dir="C:/Users/medahala/OneDrive - University of Leeds/Dropbox/Apps/Overleaf/Thesis/Figures_Tables/Chapter_2/" 

ggsave (file = paste0(output_dir,"fig_multdose_tox.pdf"), plot=a, width = 5, height = 4.5,   units = "in" , scale=1.75)



p1_<-ggplot(BIG, aes(x=sm_gamma)) +
  geom_density() +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank() 
  ) + xlab("Mean Difference") + ylab("Density")
p2_<-ggplot(BIG, aes(x=ssd_gamma)) +
  geom_density() +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank() 
  ) + xlab("Ratio of Standard deviation") + ylab("Density")

#zeta
p3_<-ggplot(BIG, aes(x=sm_zeta)) +
  geom_density() +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank() 
  ) + xlab("Mean Difference") + ylab("Density")
p4_<-ggplot(BIG, aes(x=ssd_zeta)) +
  geom_density() +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank() 
  ) + xlab("Ratio of Standard deviation") + ylab("Density")

b=ggarrange(p1_,p2_,p3_,p4_)

output_dir="C:/Users/medahala/OneDrive - University of Leeds/Dropbox/Apps/Overleaf/Thesis/Figures_Tables/Chapter_2/" 

ggsave (file = paste0(output_dir,"fig_multdose_eff.pdf"), plot=b, width = 5, height = 4.5,   units = "in" , scale=1.75)



