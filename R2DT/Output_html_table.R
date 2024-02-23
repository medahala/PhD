
library(dplyr)
library(tidyr)
library(kableExtra)
library(DT)
library(ggplot2)

setwd("XXXXXXXXX")

###################################################################


###############################################################################
#Data generation process for simulated model
###############################################################################
# open previously generated scenarios
scens_eff = readRDS("./data/scenerios.rds")[[1]]
scens_tox = readRDS("./data/scenerios.rds")[[2]]


###################################################################################################################################
#Model constants
########################################################################

mod_constants= readRDS("./data/mod_constants.rds")
mod_constants_e= readRDS("./data/mod_constants.rds")[-c(1:4)]
mod_constants_t= readRDS("./data/mod_constants.rds")[-c(5:10)]   


#############################################################################################
#Parameters for decision models
###########################################################################################
D_args=readRDS("./data/decision_constants.rds")


#####################################################################################################################
#Output
#########################################################################################################################

#These are plots of each scenerio
par(mfrow=c(2,2))
for (i in 7:10){
  
  plot(scens_tox[i,], ylim=c(0,1), ann=FALSE, type="o", col="red", lwd=3, xaxt="n")
  lines(scens_eff[i,], type="o", pch=22, col="green", lwd=3)
  text(scens_tox[i,], labels = scens_tox[i,], cex = 0.9, pos = 3)
  text(scens_eff[i,], labels = scens_eff[i,], cex = 0.9, pos = 3)
  axis(side=1,at=c(1,2,3,4))
  title(xlab="Dose",ylab="Probability", main = paste("Scenerio",i))
  abline(h=0.5, col="green", lty=2, lwd=1)
  abline(h=0.4, col="red", lty=2, lwd=1)

}




# 

#write a simple latex table for model outputs





#####



#this is needed in the chunk below 
cohortsize=3

#read in, unlist and combine data for each cohort, add variable for total sample size
results=tibble(do.call(rbind.data.frame,readRDS("./Data/sim_allcohorts_2.rds"))) %>%
  mutate(ssize = rowSums(select(., starts_with("n_")))) %>%
#ssize is not quite right since trials that stop will have a duplicate ssize- to solve create a duplication number and add  
   group_by(ssize, scenerio, dec_func, sim) %>%
  # add row number which works per group due to prior grouping
   mutate(duplicateID = row_number()) %>%
  # ungroup to prevent unexpected behaviour down stream
    ungroup() %>%
    mutate(ssize = ssize + (duplicateID-1)*cohortsize)



#selection proportions
selection =   results %>%
  group_by(ssize, scenerio, dec_func, .drop = FALSE) %>%  
  count(dose) %>%        # now required with changes to dplyr::count()
  mutate(prop = round(100*prop.table(n),1)) %>%
  mutate(name= paste0('sel_',dose)) %>%
  select(-c(n, dose))%>%
  pivot_wider(names_from = name, values_from = prop, values_fill = 0)
# the values will be the result of the fight


#selection proportions
n_treated = results %>%
  pivot_longer(cols = starts_with("n_"), 
               names_to = "numbers_dose",
               names_prefix = "n_" , 
               values_to = "count") %>%
  mutate(name= paste0('npat_',numbers_dose)) %>%
  group_by(ssize, scenerio, dec_func, name) %>%
  summarise(mpat=mean(count)) %>%
  mutate(mpat=round(mpat,1) ) %>%
  pivot_wider(names_from = name, values_from = mpat,  values_fill = 0) 


a= full_join(selection, n_treated, by = c("ssize","scenerio","dec_func") )

#find the utility values for each scenerio and decision function


n_decfunc_u=length(unique(a$dec_func))
n_scen_u= length(unique(a$scenerio))
#calculate u

numdoses= length(scens_eff[1,]) 

efftox_ucalc= function(ep, tp, pie, pit, r){
  min = 1 - (( (1) / (1 - pie))^ r + ( 1 / pit)^ r)^ (1/r)
  u = 1 - (( (1 - ep) / (1 - pie))^ r + ( tp / pit)^ r)^ (1/r)
  u2 = (u - min)/(1-min)
return(u2)
}

#vector giving information for which function to apply 
d_func = rep("doseDet::UdetC",n_decfunc_u)
d_func[c(6,13,20)] = "efftox_ucalc" 
args_f = D_args

for (j in 1:n_scen_u){
  for (i in 1:n_decfunc_u){
    
    u = R.utils::doCall(  eval(str2expression(d_func[i])), ep=as.matrix(scens_eff[j,]), tp=as.matrix(scens_tox[j,]), h_treated=5, 
                        er=args_f[i,"ref_p_e"],   el=args_f[i,"lam_e"],      eag=args_f[i,"pow_g_e"],   
                        eal=args_f[i,"pow_l_e"],  tr=args_f[i,"ref_p_t"],    tl=args_f[i,"lam_t"],   
                        tag=args_f[i,"pow_g_t"],  tal=args_f[i,"pow_l_t"],   k1=args_f[i,"k1"],        
                        k2=args_f[i,"k2"],        ptt=args_f[i,"pt_thres"],  
                        pet=args_f[i,"pe_thres"], ptc=args_f[i,"pt_cut"],    pec=args_f[i,"pe_cut"],
                        r=args_f[i,"r"],  pie=args_f[i,"pi_eff"], pit=args_f[i,"pi_tox"], ustop=args_f[i,"ustop"],
                        pstop = args_f[i,"pstop"], .ignoreUnusedArgs=TRUE) 
    
    u2=round(u,3)
    dfi = data.frame(scenerio=j, dec_func=i, u_1= u2[1], u_2= u2[2], u_3= u2[3], u_4= u2[4])  
    
    if (j==1&i==1) {
      df=dfi
    } else{
      df=rbind(df,dfi)
    } 
    
  }}

b = tibble(df )





res_= left_join(a, b, by = c("scenerio","dec_func") )  %>%
  select (ssize,scenerio, dec_func, 
          u_1, u_2, u_3, u_4,
          sel_0, sel_1, sel_2, sel_3, sel_4,
          npat_1, npat_2, npat_3, npat_4)



res = res_ %>%
  mutate(ssize=as.factor(ssize), dec_func = as.factor(dec_func), scenerio = as.factor(scenerio))

sketch = htmltools::withTags(table(
  class = 'display',
  thead(
    tr(
      th(rowspan = 2, 'Sample size'),
      th(rowspan = 2, 'Scenerio'),
      th(rowspan = 2, 'Decision Function'),
      th(colspan = 4, 'Utility of Dose'),
      th(colspan = 5, 'Selection %'),
      th(colspan = 4, 'Number of Patients')
      
    ),
    tr(
      lapply(c(c(1:4),c(0,1:4),c(1:4)) , th)
    )
  )
))



datatable(res, rownames = FALSE, filter = 'top', container=sketch, options=list(pageength = 12)) %>%
  formatStyle(c(7,12), `border-right` = "solid 2px")
colnames= c("Sample Size","Scenerio","Decision Function", c(1:4),c(0,1:4),c(1:4)) 



