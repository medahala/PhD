#packages used in program
library(rstansim)
library(rstan)
library(ggplot2)
library(ggjoy)
library(bayesplot)
library(gtools)
library(dplyr)
library(tidyr)
library(kableExtra)

#directory to save simulations
setwd("C:/Users/medahala/OneDrive - University of Leeds/PhD/R Code/copula_sim")


#generate every possible combination of outcomes for 20 patients
res1 <- combinations(n= 4, r = 20, repeats.allowed=TRUE)
#note this is a combination as order isn't important

#for this simulation only want those with 14 efficacy responses
# and 6 toxicty responses 
N = rowSums(res1==1)
E = rowSums(res1==2)
T = rowSums(res1==3)
B = rowSums(res1==4)

#Single example of 12 effiacy responses 6 toxicity respons
mat=res1[(B+E==14) & (B+T==6),]

#number of both effiacy and toxicty events
matB=rowSums(mat==4)

#stan models compile
FGM_stan = stan_model(file= 'FGM.stan')
gaus_stan= stan_model(file= 'Gaus.stan')
ind_stan = stan_model(file= 'ind.stan')


#fit the models for gaussian and FGM copulas
 for (i in 0:6){
in_dat <- list(
  num_patients=20,  
  peshape1=1,
  peshape2=1,
  ptshape1=1,
  ptshape2=1,
  psi_l=-1,
  psi_u=1,
  eff =  1 * mat[which(matB==i),] %in% c(4,2),
  tox = 1 * mat[which(matB==i),] %in% c(4,3)
)
assign(paste0("FGM_fit_",i), as.data.frame( sampling(FGM_stan, data = in_dat, chains=5 , iter=10000, warmup=1000, thin=1, seed="1263465893")))
assign(paste0("Gaus_fit_",i), as.data.frame( sampling(gaus_stan, data = in_dat, chains=5 , iter=10000, warmup=1000, thin=1, seed="2376592304")))
 
}
 
#fit the indpendent model
ind_fit = as.data.frame(as.data.frame(sampling(ind_stan, data = in_dat, chains=5 , iter=10000, warmup=1000, thin=1, seed="485875934")))
#note this doesn't need fitting with a stan model as known paramteric form for posterior
#this is done as easier to combine and calculate summary stats the same as other models 


### get the kendals tau estimates the data
cor=rep(0,7)
  for (i in 0:6){ 
  eff =  1 * mat[which(matB==i), ] %in% c(4,2)
  tox = 1 * mat[which(matB==i), ] %in% c(4,3)
  cor[(i+1)] = cor(x=eff, y=tox, method="kendall")
}
  
ind_fit$lp__=NULL
ind_fit$psi=NA
ind_fit$tau=NA
  

#combine the data        
dat =  rbind(ind_fit, FGM_fit_0[,-c(5)], FGM_fit_1[,-c(5)], FGM_fit_2[,-c(5)], FGM_fit_3[,-c(5)], FGM_fit_4[,-c(5)], 
                      FGM_fit_5[,-c(5)], FGM_fit_6[,-c(5)],
             Gaus_fit_0[,-c(5)], Gaus_fit_1[,-c(5)], Gaus_fit_2[,-c(5)], Gaus_fit_3[,-c(5)], Gaus_fit_4[,-c(5)], 
             Gaus_fit_5[,-c(5)], Gaus_fit_6[,-c(5)])


dat$fit = rep(1:15, each=45000)
dat$g = factor(dat$fit,  labels = c("Independent","FGM Copula, Ties=0 ","FGM Copula, Ties=1","FGM Copula, Ties=2","FGM Copula, Ties=3","FGM Copula, Ties=4","FGM Copula, Ties=5","FGM Copula, Ties=6",
                                                  "Gaussian Copula, Ties=0 ","Gaussian Copula, Ties=1","Gaussian Copula, Ties=2","Gaussian Copula, Ties=3","Gaussian Copula, Ties=4","Gaussian Copula, Ties=5","Gaussian Copula, Ties=6") )



#Table for Output - cacluate sumnmary stattistics
a = tibble(dat) %>%
    group_by(g) %>%
    summarise(pe_10=quantile(pe,probs=0.05), pe_25=quantile(pe,probs=0.25), pe_50=quantile(pe,probs=0.5), pe_75=quantile(pe,probs=0.75), pe_90=quantile(pe,probs=0.95), 
              pt_10=quantile(pt,probs=0.05), pt_25=quantile(pt,probs=0.25), pt_50=quantile(pt,probs=0.5), pt_75=quantile(pt,probs=0.75), pt_90=quantile(pt,probs=0.95),
              t_10=quantile(tau,probs=0.05, na.rm = TRUE), t_25=quantile(tau,probs=0.25, na.rm = TRUE), t_50=quantile(tau,probs=0.5, na.rm = TRUE), 
              t_75=quantile(tau,probs=0.75, na.rm = TRUE), t_90=quantile(tau,probs=0.95, na.rm = TRUE))

a$samp_taub= c(NA,cor,cor)

#make pretty
tab_out = a %>%
  mutate(p_e = paste0(sprintf("%.2f", pe_50),' (',sprintf("%.2f", pe_10),', ',sprintf("%.2f", pe_90),')'),
         p_t = paste0(sprintf("%.2f", pt_50),' (',sprintf("%.2f", pt_10),', ',sprintf("%.2f", pt_90),')'),
         t =   ifelse(is.na(samp_taub),NA, paste0(sprintf("%.2f", t_50),' (',sprintf("%.2f", t_10),', ',sprintf("%.2f", t_90),')')),
         st = ifelse(is.na(samp_taub),NA,sprintf("%.2f", samp_taub)) ) %>%
         select(g, p_e,p_t,t,st)

#colnames for output
colnames(tab_out) = c("Model fit and data","Median $\\pi_E$ (90\\% CI)","Median $\\pi_T$ (90\\% CI)", "Median $\\tau_b$ (90\\% CI)", "Data $\\tau_b$" )

#create table for latex
p = kable(tab_out,  format = "latex", escape = FALSE, booktabs=TRUE,
          caption = "Needs editing",
          label = "Small_example") %>% kable_styling(position = "center")


output_dir="C:/Users/medahala/OneDrive - University of Leeds/Dropbox/Apps/Overleaf/Thesis/Figures_Tables/Chapter_2/" 

#save
sink(paste0(output_dir,'small_example.tex'), append=FALSE)
print(p, append=T,table.placement = "h",
        caption.placement="bottom", hline.after=seq(from=-1,to=nrow(p),by=1, na.print = ""))
sink()

#calculate probabilities for example in text with known sampling distributions
fgm = function(p_e, p_t, psi){
  vec = (1-p_e) * (1-p_t) + psi * (1-p_e) * (1-p_t) * p_e * p_t
  vec[2]=  (1-p_e) - vec[1] 
  vec[3]=  (1-p_t) - vec[1] 
  vec[4] = 1 - (1-p_e) - (1-p_t) + vec[1] 
  return(vec)
} 

#probabilities for the 4 outcomes
p=fgm(0.7,0.3,0.8)
p2=fgm(0.7,0.3,0)

N_ = rowSums(mat==1)
E_ = rowSums(mat==2)
T_ = rowSums(mat==3)
B_ = rowSums(mat==4)

X = matrix(c(N_,E_,T_,B_),ncol = 4 )

chance = function(X,p)
{
  exp( X%*%log(p) - log(factorial(X))%*%rep(1,4) + log(factorial(rowSums(X))))
}

weight=chance(X,p)

out = weight / sum(weight) 

round(100*out)/100

weight2=chance(X,p2)
out2 = weight2 / sum(weight2) 
round(100*out2)/100

#calculate tau_b fro the example
p=0.3; q=0.7; theta=0.8
fgm_r=(1-p)*(1-q)*(1+theta*p*q)
((fgm_r-(1-p)*(1-q))/(sqrt(p*(1-p)*q*(1-q))))




