
setwd("C:/Users/medahala/OneDrive - University of Leeds/PhD/R code/ill_example/four_dose/paper_sims")

#generate scenerios function
scenereo_pat_sim = function(nsims, max_pats, pe, pt, psi){
out=array(NA, dim=c(max_pats, length(pe), nsims, 2))
for (d in 1:length(pe)){
p_t = pt[d]
p_e = pe[d]
p_T = (1-p_e)*p_t    -  p_t*p_e*(1-p_t)*(1-p_e)*psi #tox only
p_E = p_e*(1-p_t)    -  p_t*p_e*(1-p_t)*(1-p_e)*psi #eff only
p_B = p_e*p_t        +  p_t*p_e*(1-p_t)*(1-p_e)*psi #eff and tox
p_N = (1-p_e)*(1-p_t)+  p_t*p_e*(1-p_t)*(1-p_e)*psi #neither
x = sample(c("T","E","B","N"), max_pats * nsims, replace=TRUE, prob=c(p_T,p_E,p_B,p_N))
out[,d,,1] = x %in% c("E","B")
out[,d,,2] = x %in% c("T","B")
}
return(out)
} 



###############################################################################
#Data generation process for simulated model
###############################################################################

#trial doses
actual_doses=c(20,30,40,50)
cd=log(actual_doses)- mean(log(actual_doses))

#function to fit linear odds model to given probabilities at given doses -  a vector of length 2
#returns vector of probabilities
#will fit a straight line and  
p_expand = function(p,d,c_d=cd){
  X=cbind(rep(1,length(c_d)), c_d, c_d^2)
  a <- X[d,1:length(p)]
  b = log(p/(1-p))
  betas=solve(a, b)
  return(1/(1+exp(-X[,1:length(p)] %*% betas)))
}
#scenerio 0
t_0 =round( p_expand(p=c(0.05,0.15), d=c(1,4)) , 2)
e_0 =round( p_expand(p=c(0.3,0.85), d=c(1,4))  , 2)

scen_0=list(p_t=t_0,p_e=e_0)

#scenerio 0a
t_0a =round( p_expand(p=c(0.05,0.15), d=c(1,4)) , 2)
e_0a =round( p_expand(p=c(0.45,0.55), d=c(2,4))  , 2)
scen_0a=list(p_t=t_0a,p_e=e_0a)


#scenerio1
t_1 =round( p_expand(p=c(0.05,0.35), d=c(1,4)) , 2)
e_1 =round( p_expand(p=c(0.3,0.85), d=c(1,4))  , 2)

scen_1=list(p_t=t_1,p_e=e_1)

#scenerio2
#same toxicty
e_2 =round( p_expand(p=c(0.45,0.55), d=c(2,4))  , 2)
scen_2=list(p_t=t_1,p_e=e_2)


#scenerio4 - now the last one (7)
e_4 =round( p_expand(p=c(0.6,0.7,0.7), d=c(2,3,4))  , 2)

scen_4=list(p_t=t_1,p_e=e_4)

fit <- lm(e_4 ~ poly(cd, 2))   ## polynomial of degree 3
#plot(cd, e_4)  ## scatter plot (colour: black)

x0 <- seq(min(cd), max(cd), length = 20)  ## prediction grid
y0 <- predict.lm(fit, newdata = list(cd = x0))  ## predicted values
#lines(x0, y0, col = 2)  ## add regression curve (colour: re


#scenerio5
t_5 = round( p_expand(p=c(0.35,0.42), d=c(1,2)) , 2)
e_5 =round( p_expand(p=c(0.55,0.9), d=c(1,4))  , 2)

scen_5=list(p_t=t_5,p_e=e_5)

#scenerio6
t_6 = round( p_expand(p=c(0.35,0.42), d=c(2,3)) , 2)
e_6 = round( p_expand(p=c(0.6,0.62), d=c(1,2)) , 2)

scen_6=list(p_t=t_6, p_e=e_6)


#scenario 8 - dose level 2 is optimum
t_8 = round( p_expand(p=c(0.35,0.5), d=c(2,3)) , 2)
e_8 = round( p_expand(p=c(0.6,0.7,0.7), d=c(2,3,4)) , 2)
scen_8=list(p_t=t_8, p_e=e_8)

#scenario 9 - all too toxic
t_9 = round( p_expand(p=c(0.45,0.7), d=c(1,4)) , 2)
e_9 =round( p_expand(p=c(0.55,0.9), d=c(1,4))  , 2)
scen_9=list(p_t=t_9, p_e=e_9)
scen_9

#scenerio 10 - no e
t_10 =round( p_expand(p=c(0.05,0.15), d=c(1,4)) , 2)
e_10 =round( p_expand(p=c(0.3,0.45), d=c(2,4))  , 2)

scen_10 = list(p_t=t_10, p_e=e_10)

scens_tox= matrix(c(scen_0$p_t, scen_0a$p_t, scen_1$p_t, scen_2$p_t,scen_5$p_t,scen_6$p_t, scen_4$p_t, scen_8$p_t,
                    scen_9$p_t, scen_10$p_t), 
                  nrow=10, byrow=TRUE)

scens_eff= matrix(c(scen_0$p_e, scen_0a$p_e, scen_1$p_e,scen_2$p_e,scen_5$p_e,scen_6$p_e, scen_4$p_e, scen_8$p_e,
                    scen_9$p_e, scen_10$p_e ), 
                  nrow=10, byrow=TRUE)

# 
# par(mfrow=c(3,2))
# 
# for (i in 1:10){
#   
#   plot(scens_tox[i,], ylim=c(0,1), ann=FALSE, type="o", col="red", lwd=3, xaxt="n")
#   lines(scens_eff[i,], type="o", pch=22, col="green", lwd=3)
#   text(scens_tox[i,], labels = scens_tox[i,], cex = 0.6, pos = 3)
#   text(scens_eff[i,], labels = scens_eff[i,], cex = 0.6, pos = 3)
#   axis(side=1,at=c(1,2,3,4))
#   title(xlab="Dose",ylab="Probability", main = paste("Scenerio",i))
#   abline(h=0.5, col="green", lty=2, lwd=1)
#   abline(h=0.4, col="red", lty=2, lwd=1)
#   
# }
# # 

saveRDS(list(scens_eff,scens_tox), "./data/scenerios.rds")


###################################################################################################################################
#Simulate data
########################################################################
#diminesions for array of form c(number of patients, number of doses, number of sims, 2 for effan tox, number of scenerios)

 sim_dat=array(0L, dim=c(60, 4, 3000, 2, 10))
 set.seed(123456789)
 for (i in 1:10){
   sim_dat[,,,,i]= scenereo_pat_sim(nsims=3000, max_pats = 60, pe=scens_eff[i,], pt = scens_tox[i,], psi=0)
 }
  saveRDS(sim_dat, "./data/simulated_data.rds")


