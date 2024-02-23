#This program runs the single dose example

library(rstan)
rstan_options(auto_write = TRUE)
library(doParallel)
library(GGally)
library(ggpubr)
library(gtools)
library(dplyr)
library(tidyr)
library(kableExtra)

#all possible combinations of data
res1 <- combinations(n= 4, r = 20, repeats.allowed=TRUE)
#note this is a combination as order isn't important

eff = matrix( 1 * res1 %in% c(4,2), ncol=20)
tox = matrix( 1 * res1 %in% c(4,3), ncol=20)
setwd("C:/Users/medahala/OneDrive - University of Leeds/PhD/R Code/copula_sim")

# uncomment to run code 
# FGM_stan = stan_model(file= 'FGM.stan')
# gaus_stan= stan_model(file= 'Gaus.stan')
# ind_stan = stan_model(file= 'ind.stan')
# 
# #function that will open one external data file, fit a stan model and save outputs 
# one_model_fit_and_save = function(stan_model_object, lab, i, dat_e, dat_t) {
#   saveRDS( 
#     suppressWarnings(rstan::extract(rstan::sampling(stan_model_object, 
#                                                     data=   list(
#                                                       num_patients=20,  
#                                                       peshape1=1,
#                                                       peshape2=1,
#                                                       ptshape1=1,
#                                                       ptshape2=1,
#                                                       psi_l=-1,
#                                                       psi_u=1,
#                                                       eff =  dat_e[i,],
#                                                       tox =  dat_t[i,]
#                                                     ),
#                                                     iter= 10000,
#                                                     warmup=1000,
#                                                     chains=5L,
#                                                      seed="23475852",
#                                                     cores = 1L,
#                                                     refresh=0),pars=c("pe","pt","tau")))
#     ,file=paste0("./stan_fits/",lab,"_", i,".rds"), compress=FALSE
#   )
#   invisible(gc())
# }
# 
# ncores=parallel::detectCores()
# #fit all stan models for efficacy
# cl = makeCluster(ncores)
# registerDoParallel(cl)
# garb = foreach(i = 1:nrow(eff), .inorder=FALSE) %dopar% {
#   one_model_fit_and_save(FGM_stan,"FGM",i, eff, tox)
# }
# stopCluster(cl)
# 
# #fit all stan models for efficacy
# cl = makeCluster(ncores)
# registerDoParallel(cl)
# garb = foreach(i = 1:1:nrow(eff), .inorder=FALSE) %dopar% {
#   one_model_fit_and_save(gaus_stan,"Gaus",i, eff, tox)
# }
# stopCluster(cl)


N=nrow(eff)
m_FGM_pe=m_FGM_pt=m_FGM_tau=sd_FGM_pe=sd_FGM_pt=sd_FGM_tau=FGM_tau_l=FGM_tau_u=rep(NA,N) 
m_G_pe=m_G_pt=m_G_tau=sd_G_pe=sd_G_pt=sd_G_tau=G_tau_l=G_tau_u=rep(NA,N)

#calcualte summary statistics
for (i in 1:N){
dat=readRDS(file=paste0("./stan_fits/FGM_", i,".rds"))
m_FGM_pe[i]= mean(dat$pe)
m_FGM_pt[i]= mean(dat$pt)
m_FGM_tau[i]= mean(dat$tau)
FGM_tau_l[i] = quantile(dat$tau,0.05)
FGM_tau_u[i] = quantile(dat$tau,0.95)
sd_FGM_pe[i]= sd(dat$pe)
sd_FGM_pt[i]= sd(dat$pt)
sd_FGM_tau[i]= sd(dat$tau)
dat=readRDS(file=paste0("./stan_fits/Gaus_", i,".rds"))
m_G_pe[i]= mean(dat$pe)
m_G_pt[i]= mean(dat$pt)
m_G_tau[i]= mean(dat$tau)
G_tau_l[i] = quantile(dat$tau,0.05)
G_tau_u[i] = quantile(dat$tau,0.95)
sd_G_pe[i]= sd(dat$pe)
sd_G_pt[i]= sd(dat$pt)
sd_G_tau[i]= sd(dat$tau)
}

alpha_e = 1 + rowSums(eff)
alpha_t = 1 + rowSums(tox)  
beta_e = 1 + 20 - rowSums(eff)
beta_t = 1 + 20 - rowSums(tox)
m_I_pe = alpha_e / (alpha_e+beta_e)
m_I_pt = alpha_t / (alpha_t+beta_t)
sd_I_pe = sqrt( alpha_e * beta_e / ((alpha_e+beta_e +1)*(alpha_e + beta_e)^2 ) ) 
sd_I_pt = sqrt( alpha_t * beta_t / ((alpha_t + beta_t +1)*(alpha_t + beta_t)^2 ) ) 

cor=rep(NA,N)
for (i in 1:N){
  cor[i] = cor(x=eff[i,], y=tox[i,], method="kendall")
}


#comboine all summary stats
BIG = tibble(rbind( data.frame(m_pe = m_FGM_pe, m_pt = m_FGM_pt, m_tau= m_FGM_tau, sd_pe = sd_FGM_pe, sd_pt=sd_FGM_pt, sd_tau=sd_FGM_tau,
             sm_pe = m_FGM_pe - m_I_pe,
             sm_pt=  m_FGM_pt - m_I_pt, 
             ssd_pe = sd_FGM_pe/sd_I_pe, 
             ssd_pt = sd_FGM_pt/sd_I_pt,
             tau_l = FGM_tau_l,  tau_u=FGM_tau_u,
             mod = rep(1,N), cor=cor) ,
             data.frame(m_pe = m_G_pe, m_pt = m_G_pt, m_tau= m_G_tau, sd_pe = sd_G_pe, sd_pt=sd_G_pt, sd_tau=sd_G_tau,
             sm_pe = m_G_pe - m_I_pe,
             sm_pt=  m_G_pt - m_I_pt, 
             ssd_pe = sd_G_pe/sd_I_pe, 
             ssd_pt = sd_G_pt/sd_I_pt,
             tau_l = G_tau_l,  tau_u=G_tau_u,
             mod = rep(2,N), cor=cor)
))


BIG$m = factor(BIG$mod, labels = c("FGM","Gaussian"))

BIG$sssd_pe =log(BIG$ssd_pe)

ggpairs(BIG[BIG$mod==1,],          # Data frame
        columns = c("sm_pt", "ssd_pt")) # Columns


ggpairs(BIG[BIG$mod==2,],          # Data frame
        columns = c("sm_pe", "sssd_pe", "m_tau")) # Columns

#effiacy
p1<-ggplot(BIG, aes(x=sm_pe)) +
  geom_density(adjust=1.5)+facet_grid(m ~ .) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank() 
  ) + xlab("Mean Difference") + ylab("Density")
  

p2<-ggplot(BIG, aes(x=ssd_pe)) +
  geom_density(adjust=1.5)+facet_grid(m ~ .) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank() 
  ) + xlab("Ratio of Standard deviation") + ylab("Density")


a= ggarrange(p1,p2)

output_dir="C:/Users/medahala/OneDrive - University of Leeds/Dropbox/Apps/Overleaf/Thesis/Figures_Tables/Chapter_2/" 

ggsave (file = paste0(output_dir,"fig_single_d_cop_pe.pdf"), plot=a, width = 5, height = 4.25,   units = "in" , scale=1.75)


#tox
p1<-ggplot(BIG, aes(x=sm_pt)) +
  geom_density(adjust=1.5)+facet_grid(m ~ .) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank() 
  ) + xlab("Mean Difference") + ylab("Density")


p2<-ggplot(BIG, aes(x=ssd_pt)) +
  geom_density(adjust=1.5)+facet_grid(m ~ .) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank() 
  ) + xlab("Ratio of Standard deviation") + ylab("Density")

b= ggarrange(p1,p2)

output_dir="C:/Users/medahala/OneDrive - University of Leeds/Dropbox/Apps/Overleaf/Thesis/Figures_Tables/Chapter_2/" 

ggsave (file = paste0(output_dir,"fig_single_d_cop_pt.pdf"), plot=b, width = 5, height = 4.25,   units = "in" , scale=1.75)



#investigate=BIG[abs(BIG$sm_pe)>0.01,]


N = rowSums(res1==1)
E = rowSums(res1==2)
T = rowSums(res1==3)
B = rowSums(res1==4)
X = matrix(c(N,E,T,B),ncol = 4 )

#p = c((1-e)*(1-t),e*(1-t),(1-e)*t,e*t)



#fgm
fgm = function(p_e, p_t, psi){
  vec = (1-p_e) * (1-p_t) + psi * (1-p_e) * (1-p_t) * p_e * p_t
  vec[2]=  (1-p_e) - vec[1] 
  vec[3]=  (1-p_t) - vec[1] 
  vec[4] = 1 - (1-p_e) - (1-p_t) + vec[1] 
  return(vec)
}  

normal = function(p_e, p_t, psi){
  vec = pbivnorm::pbivnorm(qnorm(1-p_e),qnorm(1-p_t),rho=psi)
  vec[2]=  (1-p_e) - vec[1] 
  vec[3]=  (1-p_t) - vec[1] 
  vec[4] = 1 - (1-p_e) - (1-p_t) + vec[1] 
  return(vec)
}  

p = fgm(0.7,0.3,0.8)

chance = function(X,p)
  {
 exp( X%*%log(p) - log(factorial(X))%*%rep(1,4) + log(factorial(rowSums(X))))
}

weight=chance(X,p)

BIG$weight=c(weight,weight)


#effiacy
p1<-ggplot(BIG, aes(x=sm_pe, weight=weight)) +
  geom_density(adjust=1.5)+facet_grid(m ~ .) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank() 
  ) + xlab("Mean Difference") + ylab("Density") + xlim(c(-0.015,0.015 )) 


p2<-ggplot(BIG, aes(x=ssd_pe, weight=weight)) +
  geom_density(adjust=1.5)+facet_grid(m ~ .) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank() 
  ) + xlab("Ratio of Standard deviation") + ylab("Density") + xlim(c(0.975,1.025 )) 

c=ggarrange(p1,p2)

ggsave (file = paste0(output_dir,"fig_single_d_cop_pe_ex.pdf"), plot=c, width = 5, height = 4.25,   units = "in" , scale=1.75)



#tox
p1<-ggplot(BIG, aes(x=sm_pt, weight=weight)) +
  geom_density(adjust=1.5)+facet_grid(m ~ .) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank() 
  ) + xlab("Mean Difference") + ylab("Density") + xlim(c(-0.015,0.015 )) 


p2<-ggplot(BIG, aes(x=ssd_pt, weight=weight)) +
  geom_density(adjust=1.5)+facet_grid(m ~ .) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank() 
  ) + xlab("Ratio of Standard deviation") + ylab("Density") + xlim(c(0.975,1.025 )) 

d=ggarrange(p1,p2)

ggsave (file = paste0(output_dir,"fig_single_d_cop_pt_ex.pdf"), plot=d, width = 5, height = 4.25,   units = "in" , scale=1.75)


p = fgm(0.7,0.3,0)

chance = function(X,p)
{
  exp( X%*%log(p) - log(factorial(X))%*%rep(1,4) + log(factorial(rowSums(X))))
}

weight=chance(X,p)

BIG$weight2=c(weight,weight)


#effiacy
p1<-ggplot(BIG, aes(x=sm_pe, weight=weight2)) +
  geom_density(adjust=1.5)+facet_grid(m ~ .) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank() 
  ) + xlab("Mean Difference") + ylab("Density") + xlim(c(-0.015,0.015 )) 


p2<-ggplot(BIG, aes(x=ssd_pe, weight=weight2)) +
  geom_density(adjust=1.5)+facet_grid(m ~ .) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank() 
  ) + xlab("Ratio of Standard deviation") + ylab("Density") + xlim(c(0.975,1.025 )) 

e=ggarrange(p1,p2)

ggsave (file = paste0(output_dir,"fig_single_d_cop_pe_ex1.pdf"), plot=e, width = 5, height = 4.25,   units = "in" , scale=1.75)


#effiacy
p1<-ggplot(BIG, aes(x=sm_pt, weight=weight)) +
  geom_density(adjust=1.5)+facet_grid(m ~ .) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank() 
  ) + xlab("Mean Difference") + ylab("Density") + xlim(c(-0.015,0.015 )) 


p2<-ggplot(BIG, aes(x=ssd_pt, weight=weight)) +
  geom_density(adjust=1.5)+facet_grid(m ~ .) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank() 
  ) + xlab("Ratio of Standard deviation") + ylab("Density") + xlim(c(0.975,1.025 )) 

f=ggarrange(p1,p2)

ggsave (file = paste0(output_dir,"fig_single_d_cop_pt_ex1.pdf"), plot=f, width = 5, height = 4.25,   units = "in" , scale=1.75)


ggplot(BIG, aes(x=weight)) +
  geom_density() +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank() 
  ) + xlab("Mean Difference") + ylab("Density")


#correlation plot
big_tau = tibble(T_FGM= BIG$m_tau[BIG$mod==1], cor= cor, T_G= BIG$m_tau[BIG$mod==2]) 


n_eff=rowSums(eff)
n_tox=rowSums(tox)
N = rowSums(matrix(res1 %in% c(1), ncol=20 ))
big_tau$ex2 = factor(floor(0.5 * N))

#ex_eff = 2 * (abs(n_eff - 10) > 8)
#ex = ex_tox + ex_eff

#big_tau$ex = factor(ex, labels=c("non_ex","ex_tox","ex_eff", "ex_both"))

#look at how correlation works
ggpairs(big_tau,          # Data frame
        columns = c("T_FGM", "cor")) # Columns


#plot of big

cor_plot = ggplot(BIG, aes(x=cor, y=m_tau)) +
  geom_point(size = 0.5) + facet_grid(. ~ m) +
 xlab(expression(Data ~tau[b])) + ylab(expression(Model~mean ~tau[b]))

ggsave (file = paste0(output_dir,"cop_correlation.pdf"), plot=cor_plot, width = 5, height = 2.5,   units = "in" , scale=1.75)


BIG$tau_ls=sign(BIG$tau_l) 
BIG$tau_us=sign(BIG$tau_u) 

BIG$tau_ex =  sign( BIG$tau_ls + BIG$tau_us )^2


tau_ul=ggplot(BIG, aes(x=cor, y=m_tau, ymin=tau_l, ymax=tau_u)) +
  geom_pointrange(size=0.2, aes(color=tau_ex)) + facet_grid(. ~ m) +
  xlab(expression(Data ~tau[b])) + ylab(expression(Model~CI~tau[b])) +
    theme(legend.position = "none")
   
 ggsave (file = paste0(output_dir,"cop_correlation_ul.pdf"), plot=tau_ul, width = 5, height = 2.5,   units = "in" , scale=1.75)
 
 


#chance of exclusion

sum(BIG$weight[BIG$tau_ex==1 & BIG$mod==1]) 

sum(BIG$weight[BIG$tau_ex==1 & BIG$mod==2]) 

sum(BIG$weight2[BIG$tau_ex==1 & BIG$mod==1]) 

sum(BIG$weight2[BIG$tau_ex==1 & BIG$mod==2]) 

p = normal(0.7,0.3,0.69)
weight=chance(X,p)
BIG$weight3=c(weight,weight)
sum(BIG$weight3[BIG$tau_ex==1 & BIG$mod==2]) 


p=1-0.7
q=1-0.3
gaus_r=pbivnorm::pbivnorm(qnorm(1-p),qnorm(1-q),rho=0.69)
 (gaus_r-(1-p)*(1-q))/(sqrt(p*(1-p)*q*(1-q)))





#prior plots for kendals tau

pe = runif(10000)
pt = runif(100000)
theta = runif(100000, min=-1,max=1)

p=1-pe
q=1-pt

gaus_r=pbivnorm::pbivnorm(qnorm(1-p),qnorm(1-q),rho=theta)
gaus_tau= (gaus_r-(1-p)*(1-q))/(sqrt(p*(1-p)*q*(1-q)))




fgm_r=(1-p)*(1-q)*(1+theta*p*q)
fgm_tau=((fgm_r-(1-p)*(1-q))/(sqrt(p*(1-p)*q*(1-q))))

plot2d=rbind(data.frame(tau=fgm_tau, m="FGM"), data.frame(tau=gaus_tau, m="Gaussian"))


theta = ggplot(plot2d, aes(x=tau)) +
  geom_density(adjust = 2) + facet_grid(. ~ m) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank() 
  ) + xlab(expression(~tau[b])) + ylab("Density")


ggsave (file = paste0(output_dir,"tau_prior.pdf"), plot=theta, width = 5, height = 2.5,   units = "in" , scale=1.75)





