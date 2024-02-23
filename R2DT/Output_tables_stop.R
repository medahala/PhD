
library(dplyr)
library(tidyr)
library(kableExtra)
library(DT)
library(ggplot2)

setwd("C:/Users/medahala/OneDrive - University of Leeds/PhD/R code/ill_example/four_dose/paper_sims")

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
#D_args=readRDS("./data/decision_constants.rds")

D_args=readRDS("./data/decision_constants_stop.rds")


#####################################################################################################################
#Output
#########################################################################################################################


#this is needed in the chunk below 
cohortsize=3

#read in, unlist and combine data for each cohort, add variable for total sample size
results_m=tibble(do.call(rbind.data.frame,readRDS("./Data/sim_allcohorts_2.rds"))) %>%
  filter(dec_func %in% c(1)) %>%
  mutate(ssize = rowSums(select(., starts_with("n_")))) %>%
#ssize is not quite right since trials that stop will have a duplicate ssize- to solve create a duplication number and add  
   group_by(ssize, scenerio, dec_func, sim) %>%
  # add row number which works per group due to prior grouping
   mutate(duplicateID = row_number()) %>%
  # ungroup to prevent unexpected behaviour down stream
    ungroup() %>%
    mutate(ssize = ssize + (duplicateID-1)*cohortsize)

results_s=tibble(do.call(rbind.data.frame,readRDS("./Data/sim_allcohorts_stop.rds"))) %>%
  filter(dec_func %in% c(19:24)) %>%
  mutate(ssize = rowSums(select(., starts_with("n_")))) %>%
  #ssize is not quite right since trials that stop will have a duplicate ssize- to solve create a duplication number and add  
  group_by(ssize, scenerio, dec_func, sim) %>%
  # add row number which works per group due to prior grouping
  mutate(duplicateID = row_number(), dec_func = dec_func-17) %>%
  # ungroup to prevent unexpected behaviour down stream
  ungroup() %>%
  mutate(ssize = ssize + (duplicateID-1)*cohortsize)


results=rbind(results_m, results_s) %>%
  filter( ssize %in% c(45))



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
#d_func[c(6,13,20)] = "efftox_ucalc" 
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

##############################################################################################################
#outputs 
##############################################################################################################
#############################

# #expected sample size
startdose = 1
ess_T = 1
ess_E = 1
# 
# #Mean vectors
prior_tox=round(colMeans(scens_tox[1:6,]),2)
prior_eff=round(colMeans(scens_eff[1:6,]),2)
# 
# #actual doses
actual_doses=c(20,30,40,50)
# #coded doses
cd=log(actual_doses)- mean(log(actual_doses))

#############################

#install.packages("prettyR")
#library(prettyR)

#the decsion functions are repliucates based upon a threshold prbailities
#i.e 1:7 is 5% 8:15 0.075 and 16:21 0.1

dec_func_labels = c(" ",
                    "R2DT (1)",
                    "R2DT (3i)",
                    "R2DT (3ii)",
                    "R2DT (3iii)",
                    "R2DT (4i)",
                    "R2DT (4ii)",
                    "R2DT (4iii)"
                    )
  


res_$max = apply(res_[, c("u_1","u_2","u_3","u_4")], 1, which.max)


res1 = res_ %>%
         mutate(c_1 = paste0('[',ifelse(max==1,cell_spec(sprintf("%.2f", round(u_1,2)),"latex",bold=T),sprintf("%.2f", round(u_1,2))),"] ",
                      round(sel_1,1),
                      ' (',round(npat_1,1),')'),
                c_2 = paste0('[',ifelse(max==2,cell_spec(sprintf("%.2f", round(u_2,2)),"latex",bold=T),sprintf("%.2f", round(u_2,2))),"] ",
                             round(sel_2,1),
                             ' (',round(npat_2,1),')'),
                c_3 = paste0('[',ifelse(max==3,cell_spec(sprintf("%.2f", round(u_3,2)),"latex",bold=T),sprintf("%.2f", round(u_3,2))),"] ",
                             round(sel_3,1),
                             ' (',round(npat_3,1),')'),
                c_4 = paste0('[',ifelse(max==4,cell_spec(sprintf("%.2f", round(u_4,2)),"latex",bold=T),sprintf("%.2f", round(u_4,2))),"] ",
                             round(sel_4,1),
                             ' (',round(npat_4,1),')'),
                c_0 = paste(round(sel_0,1))
                              )

#merge scenerio data
scens = matrix("",nrow=10,ncol=4)
for (i in 1:10) {
  for (j in 1:4){
    scens [i,j] = paste0('(',scens_eff[i,j],', ',scens_tox[i,j],')')
  }
}
scens_t=as.data.frame(scens)
colnames(scens_t) = c("c_1","c_2","c_3","c_4")
scens_t$scenerio = 1:nrow(scens_t)
scens_t$method = factor(0, levels = " ")
scens_t$dec_func = 0
scens_t$c_0 = " "

scens_t$ssize=0


res2 = bind_rows(res1, tibble(scens_t))

res2 = res2 %>% arrange(scenerio, ssize, dec_func)
    

###########################################################
##########Stopping Rule
##########################################################

stop_tab = res2 %>%
  mutate(method = factor(dec_func, levels = 0:7,
                         labels =  dec_func_labels)) %>%
  group_by() %>%
  select(scenerio, method, c_1, c_2, c_3, c_4, c_0)


nscen=length(unique(stop_tab$scenerio))
nmethod=length(unique(stop_tab$method)) 
grouping_index_1 = rep(nmethod,nscen)
grouping_index_2 = rep(nmethod,nscen)
names(grouping_index_1)=  paste0("\\textbf{Scenario ",unique(stop_tab$scenerio),"}", " $(\\\\pi_E, \\\\pi_T)$") 
scenerio_rows=which(stop_tab$method==" ")

stop_tab_out = stop_tab %>% select(-scenerio)

colnames(stop_tab_out)=c("Method", "20","30","40","50", "NDS")


stop_tab_out_1=stop_tab_out[1:40,]
stop_tab_out_2=stop_tab_out[41:80,]

p_1 = kable(stop_tab_out_1,  booktabs=T, escape = FALSE, format = "latex" , 
          caption = "Stopping rules of R2DT: data of form: [utility at scenerio probability $(\\pi_E, \\pi_T)$] percentage selection (average number of patients treated). Percentage of trials with no dose selected abbreviated to NDS",
          label = "R2DT_stop"
) %>% kable_styling(position = "center") %>%
  add_header_above(c(" " = 1, "Dose (mg/kg)" = 4, " " = 1)) %>%
  pack_rows(index = grouping_index_1[1:5], latex_align = "c", bold=FALSE, escape = FALSE) %>%
  row_spec(scenerio_rows[1:5], align = "c",  hline_after = TRUE) %>%
  row_spec(0, align=c('l', rep('c', 5)))

p_2 = kable(stop_tab_out_2,  booktabs=T, escape = FALSE, format = "latex" , 
            caption = "Stopping rules of R2DT: data of form: [utility at scenerio probability $(\\pi_E, \\pi_T)$] percentage selection (average number of patients treated). Percentage of trials with no dose selected abbreviated to NDS",
            label = "R2DT_stop"
) %>% kable_styling(position = "center") %>%
  add_header_above(c(" " = 1, "Dose (mg/kg)" = 4, " " = 1)) %>%
  pack_rows(index = grouping_index_1[6:10], latex_align = "c", bold=FALSE, escape = FALSE) %>%
  row_spec(scenerio_rows[1:5], align = "c",  hline_after = TRUE) %>%
  row_spec(0, align=c('l', rep('c', 5)))


output_dir="C:/Users/medahala/OneDrive - University of Leeds/Dropbox/Apps/Overleaf/Thesis/Figures_Tables/Chapter_4/" 


sink(paste0(output_dir,'R2DT_stop_1.tex'), append=FALSE)
print(p_1, append=T,table.placement = "h",
      caption.placement="bottom", hline.after=seq(from=-1,to=nrow(p),by=1))
sink()

sink(paste0(output_dir,'R2DT_stop_2.tex'), append=FALSE)
print(p_2, append=T,table.placement = "h",
      caption.placement="bottom", hline.after=seq(from=-1,to=nrow(p),by=1))
sink()




