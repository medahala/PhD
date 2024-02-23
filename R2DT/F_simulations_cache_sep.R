#######################################################################################
#Simulating function for one probability model and but potentially multiple simulated datasets and decision models
#######################################################################################


run_sim = function(
  #list of arguments
  ndose=4L, # number of doses
  csize=3, # cohort size
  maxcohort=10, # maximu, number of cohorts
  nsims=2000,   # number of simulations for each decsion model
  startdose=1, #starting dose for the trial 
  simdat_="../data/simulated_data.rds",  # list of simulated data datasets 
  
  stan_eff_model=mod_e,
  stan_tox_model=mod_t,
  # Stan Proability model
  model_args_e=mod_constants_e, # constants for proability model
  model_args_t=mod_constants_t,
  decision_args=D_args,  #matrix of decision function and asociated decision parameters
  decision_fnames=d_func,
  ncores=detectCores(),
  path="./data/sim_combined_test.rds"
){
  
  ##########################################################################
  #  Setup dataframe structure and extract constants used later
  ########################################################################  
#create an empty external list to save to ongoing
  saveRDS(list(NULL),path)  
  dir.create("./temp")    
 
   sim_dat=readRDS(simdat_)
    
  ptm_= proc.time()
  #dimensions for the simulation array of the form c(maxpatinets,doses,simulations,eff/tox(2),scenereo) 
  dim_sim=dim(sim_dat)
  
  #number of decision models being applied for particular probability model
  n_dec_funcs=length(decision_fnames)
  #number of scenereos applied
  n_scenerios=dim_sim[5]
  
  #set up dataframe structure
  df = data.frame(matrix(0L,ncol = (3*ndose + 5), nrow = n_dec_funcs * nsims * n_scenerios))
  df_colnames = c(paste0("xe_",1:ndose), paste0("xt_",1:ndose),paste0("n_",1:ndose),"scenerio",  "dec_func", "sim", "dose","stop")
  colnames(df) = df_colnames
  
  #add idnetification numbers to dataset
  df$sim=rep(1:nsims,n_dec_funcs*n_scenerios)
  df$dec_func=rep(rep(c(1:n_dec_funcs), each=nsims), n_scenerios)
  df$scenerio=rep(1:n_scenerios, each=n_dec_funcs*nsims)
  
  #set the start dose for all simulations
  df$dose = startdose
  
  rm(startdose, n_dec_funcs, n_scenerios, sim_dat)
  
  #loop for each cohort    
  for (c in 1:maxcohort) {
    #for update on time
    print(paste0("Cohort ",c," of ",maxcohort))
    #--------------------------------------------------    
    #iteratively add data previously simulated to dataframe
    #--------------------------------------------------     
    

    
    #read in simulation datasets
   sim_dat=readRDS(simdat_)

    patnums =  (c-1)*csize+1:csize #patient numbers for cohort  
    
    # this calculates the total number events for given patient numbers, dose, simulation number and scenerio
    eff = as.integer(colSums(sim_dat[patnums,,,1,], dims=1)[matrix(c(ifelse(df$dose==0, 1, df$dose),df$sim, df$scenerio), ncol = 3)])
    tox = as.integer(colSums(sim_dat[patnums,,,2,], dims=1)[matrix(c(ifelse(df$dose==0, 1, df$dose),df$sim, df$scenerio), ncol = 3)]) 
    #if the trial has stopped, dose level 1 is returned this is just so R can use vectorised functions and information isn't used at next stage 
    
  #update data frame i.e if current dose is 1 then  xe_1+eff events, X_t + tox events and n_1 + number in cohort.    
      for (d in 1:ndose){
      df[, paste0('xe_',d)] =  df[, paste0('xe_',d)] + (eff * (df$dose==d) )
      df[, paste0('xt_',d)] =  df[, paste0('xt_',d)] + (tox * (df$dose==d) )
      df[, paste0('n_',d)] =  df[, paste0('n_',d)] + (as.integer(csize) * (df$dose==d) )
    }

    rm(sim_dat, eff, tox, patnums)
    invisible(gc())
    #----------------------------------------------------------------
    #Probability model fitting
    #----------------------------------------------------------------    
    #create a list of data  for each model fit in stan and save to external file
  dir.create("./temp/stan_args")    
   
#function to write a single model fit data argument to file      
    stan_rds_writer = function (j, z, prefix, ncols=rep(c(FALSE,TRUE),each=ndose), model_args){
      id = z[ncols] > 0
      x = unlist(model_args["coded_doses"])[id]
      dat = append( model_args,
        list(
          num_patients = sum(z[ncols]),
          num_used_doses = sum(id),
          counts = as.array(z[!ncols][id]),
          dose_counts = as.array(z[ncols][id]),
          X = as.array(x)
        )
      )
      if(prefix=="E") {dat$X2 = as.array(x ^ 2)}  
      saveRDS(dat,
              file=paste0("./temp/stan_args/stan_", prefix, j,".rds") )
    }

    #reduce dataframe down to unique rows for each model fit
    #each model fit (inclusing data) is then written to single RDS file which stan will fit later 

#efficacy  
    x_n_c_e = unique(df[df$stop==0, grepl("xe_|n_", colnames(df))])
    nrow_c_e = nrow(x_n_c_e) 
    invisible(lapply(c(1:nrow_c_e), function(i) stan_rds_writer(j=i, z=x_n_c_e[i, ], prefix="E", model_args =  model_args_e)))
    x_n_c_e$eff_mod_id=1:nrow_c_e
#toxicity    
    x_n_c_t = unique(df[df$stop==0, grepl("xt_|n_", colnames(df))])
    nrow_c_t = nrow(x_n_c_t) 
    invisible(lapply(c(1:nrow_c_t), function(i) stan_rds_writer(j=i, z=x_n_c_t[i, ], prefix="T", model_args =  model_args_t)))
    x_n_c_t$tox_mod_id=1:nrow_c_t
#merge model numbers back to df       
 df = merge(merge(df,x_n_c_e, all.x=TRUE), x_n_c_t, all.x = TRUE)
  
rm(x_n_c_e, x_n_c_t)   
    


ptm= proc.time()        
cat(paste0(nrow_c_e + nrow_c_t," models"))  

#fit model and recomendations for all combinations
#parrallel computing faster but only for later cohorts as setting up (approx 1sec) takes longer than actually running code
#if number of times the model needs fitting is less than 100 then just run in serial


suppressWarnings(dir.create("./temp/stan_fits"))

#function that will open one external data file, fit a stan model and save outputs 
one_model_fit_and_save = function(stan_model_object, lid, param_out, j){
  saveRDS( 
    suppressWarnings(rstan::extract(rstan::sampling(stan_model_object, 
                                                    data=readRDS(file=paste0("./temp/stan_args/stan_",lid,j,".rds")),
                                                    iter= 3000,
                                                    warmup=1000,
                                                    chains=1L,
                                                    cores = 1L,
                                                    refresh=0),pars=c(param_out)))
    ,file=paste0("./temp/stan_fits/stan_",lid,j,".rds"), compress=FALSE
  )
  invisible(gc())
}

#fit all stan models for efficacy
cl = makeCluster(ncores)
registerDoParallel(cl)
garb = foreach(i = 1:nrow_c_e, .inorder=FALSE) %dopar% {
 one_model_fit_and_save(stan_eff_model,"E",'p_e',i)
}
stopCluster(cl)
rm(garb); invisible(gc())
#fit all stan models for toxicity
cl = makeCluster(ncores)
registerDoParallel(cl)
garb =  foreach(i = 1:nrow_c_t, .inorder=FALSE) %dopar% {
  one_model_fit_and_save(stan_tox_model,"T",'p_t',i)
}
stopCluster(cl)
rm(garb); invisible(gc())

time=as.numeric((proc.time() - ptm)[3]); print(paste0(floor(time/60)," minutes and ",round(time %% 60), " seconds"))       

    #----------------------------------------------------------------
    #Decision model fitting
    #----------------------------------------------------------------
    #calculate details of highest dose treated at previously

# calculate highest dose treated or every row in df
df$h_treated = apply(df[,which(grepl("n_", colnames(df)))], 1, function(x) max(which(x>0)))

#create a unique matrix of decision functions to be fit for unique efficacy and toxicity models
dec =  data.matrix(unique(df[df$stop==0,colnames(df) %in% c("dec_func","eff_mod_id","tox_mod_id","h_treated")]  ))


#function specification to apply decision model given     
    call_dfunc = function(x, f_name=decision_fnames,  args_f = decision_args ){
             
      mfite = readRDS(file=paste0("./temp/stan_fits/stan_E",x[2],".rds"))$p_e
      mfitt = readRDS(file=paste0("./temp/stan_fits/stan_T",x[3],".rds"))$p_t
      
     c(x,R.utils::doCall(eval(str2expression(f_name[x[1]])), ep=mfite, tp=mfitt, h_treated=x[4], 
                      er=args_f[x[1],"ref_p_e"],   el=args_f[x[1],"lam_e"],      eag=args_f[x[1],"pow_g_e"],   
                      eal=args_f[x[1],"pow_l_e"],  tr=args_f[x[1],"ref_p_t"],    tl=args_f[x[1],"lam_t"],   
                      tag=args_f[x[1],"pow_g_t"],  tal=args_f[x[1],"pow_l_t"],   k1=args_f[x[1],"k1"],        
                      k2=args_f[x[1],"k2"],        ptt=args_f[x[1],"pt_thres"],  
                      pet=args_f[x[1],"pe_thres"], ptc=args_f[x[1],"pt_cut"],    pec=args_f[x[1],"pe_cut"],
                      r=args_f[x[1],"r"],  pie=args_f[x[1],"pi_eff"], pit=args_f[x[1],"pi_tox"], ustop=args_f[x[1],"ustop"],
                      pstop = args_f[x[1],"pstop"], .ignoreUnusedArgs=TRUE) )

    }   
    

    
    ###############
    #main loop to apply decsion model to given efficicacy and tox model 
    #fit this in parrallel 
    ptm= proc.time()  
    cat(paste0(nrow(dec)," decisions"))      
    cl = makeCluster(ncores)
    registerDoParallel(cl)
    
    d_outs = foreach(j = 1:nrow(dec) ,  .combine = rbind,  .inorder=FALSE)  %dopar% call_dfunc(x=dec[j,])
    
    
    stopCluster(cl) 
    time=as.numeric((proc.time() - ptm)[3]); print(paste0(floor(time/60)," minutes and ",round(time %% 60), " seconds"))      
    
    
    #-------------------------------------------------------------
    #merge which dose to treat next cohort to unduplicated dataset
    #----------------------------------------------------------------------------------
    colnames(d_outs)[5]="newdose"
    df = merge(df,as.data.frame(d_outs), all.x=TRUE)
    df$dose = ifelse(df$stop==0,df$newdose, df$dose)                
    df$stop = ifelse(df$dose==0,1,0)                
    df = df[order(df$scenerio, df$dec_func, df$sim),df_colnames] # reorder to decsion number and keep only  
    row.names(df)=NULL  #drop any new rownames so rownames correspond with current ordering 
 
    rm(d_outs, dec )
    
#save cohorts you go to external file   
  temp=readRDS(path)  
  temp[[c]]=df
  saveRDS(temp,path)
  rm(temp)
#delete model data and fits  
  unlink("./temp/stan_fits", recursive = TRUE)
  unlink("./temp/stan_args", recursive = TRUE)
  }
  
  time=as.numeric((proc.time() - ptm_)[3])
  print(paste0("Total Time: ", floor(time/60)," minutes and ",round(time %% 60), " seconds"))
  unlink("./temp", recursive = TRUE)
  return(df)
}








