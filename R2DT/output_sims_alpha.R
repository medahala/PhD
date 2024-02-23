#remove everything with command below and restart R with(Ctrl+Shift+F10)
#rm(list = ls())
setwd("XXXXXX")


###################################################################
library(rstan)
rstan_options(auto_write = TRUE)
library(shinystan)
library(trialr)
library(doParallel)
library(Rcpp)


#open the simulation function
source("F_simulations_cache_sep.R")


###############################################################################
#Data generation process for simulated model
###############################################################################
# open previously generated scenarios
scens_eff = readRDS("./data/scenerios.rds")[[1]]
scens_tox = readRDS("./data/scenerios.rds")[[2]]

###################################################################################################################################
#Model constants
########################################################################

########################################################################
#load the Stan models
 mod_e = stan_model(file ='../Eff_logisitc.stan') 
 mod_t = stan_model(file ='../Tox_logisitc.stan')  
# 
# # obtain priors for model 
# 
# #expected sample size
# ess_T = 1
# ess_E = 1
# 
# #Mean vectors
# prior_tox=round(colMeans(scens_tox[1:6,]),2)
# prior_eff=round(colMeans(scens_eff[1:6,]),2)
# 
# #actual doses
# actual_doses=c(20,30,40,50)
# #coded doses
# cd=log(actual_doses)- mean(log(actual_doses))
# 
# #priors defined from ellicited values
# priors <- trialr::get_efftox_priors(
#   doses = actual_doses,
#   pi_T = prior_tox, ess_T = ess_T,
#   pi_E = prior_eff, ess_E = ess_E
# )
# 
# #no correlation in model
# priors$psi_mean=priors$psi_sd=NULL
# 
#  
# 
# #add number of doses and coded doses
# 
# mod_constants = append(priors, list(num_doses=4,
#                                     coded_doses = cd))
# 
# 
# saveRDS(mod_constants,"./data/mod_constants.rds")

#split constants to be efficacy and toxicity
mod_constants_e= readRDS("./data/mod_constants.rds")[-c(1:4)]
mod_constants_t= readRDS("./data/mod_constants.rds")[-c(5:10)]   



#############################################################################################
#Parameters for decision models
###########################################################################################
n_decision_models=7

#vector giving information for which function to apply 
d_func = rep("",0)

d_func[c(1,2,7)] = "doseDet::prospect_pow_addC" 
d_func[c(3,5)] = "doseDet::prospect_pow_add_u_C" 
d_func[c(4)] = "doseDet::prospect_pow_u_C" 
d_func[6] = "doseDet::efftox_add" 


#-----------------------------------------------
#specify all constants for decision models
#---------------------------------------------
#matrix has 1 row for each decision model and columns for each parameter within that decision model
X = matrix(NA, ncol=19, nrow = n_decision_models)

#name all columns so can easily find constant som
colnames(X)= c( "ref_p_e", "lam_e", "pow_g_e", "pow_l_e","ref_p_t", "lam_t",
                "pow_g_t", "pow_l_t", "k1", "k2", "pt_thres", "pe_thres", "pt_cut", "pe_cut","pi_tox",
                "pi_eff","r", "ustop", "pstop")

#thresholds and cut values are the same for all decisions models (when used)
X[c(1,2,6,7) ,"pt_thres"]=0.4
X[c(1,2,6,7),"pe_thres"]=0.5
X[c(1,2,6,7),"pt_cut" ]=0.05
X[c(1,2,6,7),"pe_cut"]=0.05


#reference points are constant for all scenerios except 6 efftox
X[c(-6),"ref_p_e"]=0.5 
X[c(-6),"ref_p_t"]=0.35 

#k1 and k2
X[c(-6),"k1"]= 0.25
X[c(7),"k1"]= 0.5


X[c(-6),"k2"]= 0.15
X[c(7),"k2"]= 0.3

#loss aversion parameters for efficacy and toxicity
X[c(-6),c("lam_e", "lam_t")]= 2
X[c(2,5,7),c("lam_e", "lam_t")]= 1


#power coefficents for r2dt are 0.7 ordinary efftox is 1 
X[c(-6),c("pow_g_e", "pow_l_e", "pow_g_t", "pow_l_t")] = 0.7  
X[c(2,5,7),c("pow_g_e", "pow_l_e", "pow_g_t", "pow_l_t")] = 1   


x=c(3)
ustopi=c(0)
args_f=X
for (i in 1:1)
ustopi[i]=
    u=R.utils::doCall(doseDet::UdetC, ep=0.5, tp=0.35, h_treated=5, 
                      er=args_f[x[i],"ref_p_e"],   el=args_f[x[i],"lam_e"],      eag=args_f[x[i],"pow_g_e"],   
                      eal=args_f[x[i],"pow_l_e"],  tr=args_f[x[i],"ref_p_t"],    tl=args_f[x[i],"lam_t"],   
                      tag=args_f[x[i],"pow_g_t"],  tal=args_f[x[i],"pow_l_t"],   k1=args_f[x[i],"k1"],        
                      k2=args_f[x[i],"k2"])
    

X[c(3,4),"ustop"]  =ustopi
X[c(3,4),"pstop"] = 0.05

#5 is different
x=c(5)
ustopi=c(0)
args_f=X
for (i in 1:1)
  ustopi[i]=
  u=R.utils::doCall(doseDet::UdetC, ep=0.5, tp=0.35, h_treated=5, 
                    er=args_f[x[i],"ref_p_e"],   el=args_f[x[i],"lam_e"],      eag=args_f[x[i],"pow_g_e"],   
                    eal=args_f[x[i],"pow_l_e"],  tr=args_f[x[i],"ref_p_t"],    tl=args_f[x[i],"lam_t"],   
                    tag=args_f[x[i],"pow_g_t"],  tal=args_f[x[i],"pow_l_t"],   k1=args_f[x[i],"k1"],        
                    k2=args_f[x[i],"k2"])


X[5,"ustop"]  =ustopi
X[5,"pstop"] = 0.05



#######################
#efftox parameters
#######################
#want contour to go through contour from utlity based design at references for efficacy and toxicity
#need three points to solve through#
#(eff_star, tox_star), (eff0,0) and (1,tox1)
eff_star = X[2,"ref_p_e"]
tox_star = X[2,"ref_p_t"]
#other points need solving
#utility function has form ax+by+cxy=D
#find D at (eff_star,tox_star) then solve above equation by rearranging i.e 
#x = (D - by)/(a+cy)  
#y = (D-  ax)/(b+cy) 
#set the contour to go through 0.5  and 0.35 i.e. referecne for efficacy and toxicity  
k1_efftox = X[2,"k1"]
k2_efftox = X[2,"k2"] 
D =  k1_efftox *eff_star  + k2_efftox * (1 - tox_star)   +((1- k1_efftox - k2_efftox) * eff_star  * (1-tox_star))
eff0 = (D - k2_efftox) / (k1_efftox + (1 - k1_efftox - k2_efftox) )  
tox1 = 1 - ((D - k1_efftox) / (k2_efftox + (1 - k1_efftox - k2_efftox) ))

#add these points to paramters list
X[6,"pi_tox"] = tox1
X[6,"pi_eff"] = eff0
X[6,"r"]= efftox_solve_p(eff0 = eff0, tox1 = tox1, eff_star = eff_star, tox_star = tox_star)



d_func=rep(d_func,3)

X_75 = X
X_75[,"pstop"] = X_75[,"pstop"] + 0.025
X_75[,"pt_cut"] = X_75[,"pt_cut"] + 0.025
X_75[,"pe_cut"] = X_75[,"pe_cut"] + 0.025

X_1 = X
X_1[,"pstop"] = X_1[,"pstop"] + 0.05
X_1[,"pt_cut"] = X_1[,"pt_cut"] + 0.05
X_1[,"pe_cut"] = X_1[,"pe_cut"] + 0.05


D_args=rbind(X,X_75,X_1)

saveRDS(D_args, "./data/decision_constants.rds")

rm(X)

#####################################################################################################################
#simulation
# #########################################################################################################################
# 
# #clear up everthing other than what is needed for next step
 rm(list=setdiff(ls(), list("mod_e","mod_t", "mod_constants_e","mod_constants_t","D_args","d_func", "run_sim")))
# # 
# # 
# # # 
sim_1 =run_sim(
  #list of arguments
  ndose=4L, # number of doses
  csize=3L, # cohort size
  maxcohort=20, # maximum number of cohorts
  nsims=2000,   # number of simulations for each decision model
  startdose=1, #starting dose for the trial
  simdat_="./data/simulated_data.rds",  # simulated data array
  stan_eff_model=mod_e,
  stan_tox_model=mod_t,
  model_args_e=mod_constants_e, # constants for proability model
  model_args_t=mod_constants_t,
  decision_args=D_args,  #matrix of decision function and associated decision parameters
  decision_fnames=d_func,
  ncores=detectCores(),
  path="./data/sim_allcohorts_2.rds" #this just collects the dataset after each cohort
)

