library(dplyr)
library(ggplot2)
library(tidyverse)

source("./initialize_ABM_v.R")
source("./helper_functions_tradeoff.R")

stochastic_NH <- function(parms, Ns, delta_t, t, num_staff){

  beta = parms[["beta"]]
  beta.s = parms[["beta.s"]]
  beta.rm = parms[["beta.rm"]]
  sigma1 = parms[["sigma1"]]
  sigma2 = parms[["sigma2"]]
  gamma = parms[["gamma"]]
  I.C_time=parms[["I.C_time"]]
  I.C = ifelse(t>I.C_time,parms[["I.C"]],0) # make I.C dependent on I.C_time
  k.HH = parms[["k.HH"]]
  k.HR=parms[["k.HR"]]
  k.RH=parms[["k.RH"]]
  k.RR=parms[["k.RR"]]
  alpha.r= parms[["alpha.r"]]
  alpha.hcw= parms[["alpha.hcw"]]
  ID= parms[["id.I"]]
  int.r = parms[["int.r"]]
  int.hcw = parms[["int.hcw"]]
  int.rooms = parms[["int.rooms"]]
  mu.C = parms[["mu.C"]]
  mu.NC = parms[["mu.NC"]]
  ppe = parms[["ppe"]]
  prop_rhcwR = parms[["prop_rhcwR"]]
  VL.threshold.r = parms[["VL.threshold.r"]]
  VL.threshold.hcw = parms[["VL.threshold.hcw"]]
  int_time = parms[["int_time"]]
  VE_s1_r1 = parms[["VE_s1_r1"]]
  VE_p1_r1 = parms[["VE_p1_r1"]]
  VE_i1_r1 = parms[["VE_i1_r1"]]
  VE_s1_hcw1 = parms[["VE_s1_hcw1"]]
  VE_p1_hcw1 = parms[["VE_p1_hcw1"]]
  VE_i1_hcw1 = parms[["VE_i1_hcw1"]]
  VE_s2_r1 = parms[["VE_s2_r1"]]
  VE_p2_r1 = parms[["VE_p2_r1"]]
  VE_i2_r1 = parms[["VE_i2_r1"]]
  VE_s2_hcw1 = parms[["VE_s2_hcw1"]]
  VE_p2_hcw1 = parms[["VE_p2_hcw1"]]
  VE_i2_hcw1 = parms[["VE_i2_hcw1"]]
  VE_s1_r2 = parms[["VE_s1_r2"]]
  VE_p1_r2 = parms[["VE_p1_r2"]]
  VE_i1_r2 = parms[["VE_i1_r2"]]
  VE_s1_hcw2 = parms[["VE_s1_hcw2"]]
  VE_p1_hcw2 = parms[["VE_p1_hcw2"]]
  VE_i1_hcw2 = parms[["VE_i1_hcw2"]]
  VE_s2_r2 = parms[["VE_s2_r2"]]
  VE_p2_r2 = parms[["VE_p2_r2"]]
  VE_i2_r2 = parms[["VE_i2_r2"]]
  VE_s2_hcw2 = parms[["VE_s2_hcw2"]]
  VE_p2_hcw2 = parms[["VE_p2_hcw2"]]
  VE_i2_hcw2 = parms[["VE_i2_hcw2"]]
  t_first_dose = parms[["t_first_dose"]]
  #pr_each_day = 1/length(t_first_dose),
  second_dose_delay = parms[["second_dose_delay"]]
  immune_delay = parms[["immune_delay"]]
  vax_coverage_r = parms[["vax_coverage_r"]]
  vax_coverage_hcw = parms[["vax_coverage_hcw"]]
  vax_available=parms[["vax_available"]]
  refusal_hcw= parms[["refusal_hcw"]]
  refusal_r = parms[["refusal_r"]]
  comm.vax = parms[["comm.vax"]]
  visitors = parms[["visitors"]]
  int.vax = parms[["int.vax"]]
  
  VL.threshold.r <- parms[["VL.threshold.r"]]    # threshold for detectable VL (put arbitrary # in for now)
  VL.threshold.hcw <- parms[["VL.threshold.hcw"]]  
  test_freq_hcw <- parms[["test_freq_hcw"]] # staff testing frequency (days)
  test_freq_r <- parms[["test_freq_r"]]  # resident testing frequency (days)
  test_delay_hcw <- parms[["test_delay_hcw"]] # staff test delay (days)
  test_delay_r <- parms[["test_delay_r"]]  # resident test delay (days)
  
  wane_day_r = parms[["wane_day_r"]] # day waning starts for residents
  wane_rate_r = parms[["wane_rate_r"]] # daily waning rate for residents 
  wane_day_hcw = parms[["wane_day_hcw"]] # day waning starts for residents
  wane_rate_hcw = parms[["wane_rate_hcw"]] # daily waning rate for residents 

  
  S.rNC=Ns[["S.rNC"]]
  E.rNC=Ns[["E.rNC"]]
  A.rNC=Ns[["A.rNC"]]
  I.rNC=Ns[["I.rNC"]]
  R.rNC=Ns[["R.rNC"]]
  I.rC=Ns[["I.rC"]]
  R.rC=Ns[["R.rC"]]
  S.hcwNC=Ns[["S.hcwNC"]]
  E.hcwNC=Ns[["E.hcwNC"]]
  A.hcwNC=Ns[["A.hcwNC"]]
  I.hcwNC=Ns[["I.hcwNC"]]
  R.hcwNC=Ns[["R.hcwNC"]]
  S.hcwC=Ns[["S.hcwC"]]
  E.hcwC=Ns[["E.hcwC"]]
  A.hcwC=Ns[["A.hcwC"]]
  I.hcwC=Ns[["I.hcwC"]]
  R.hcwC=Ns[["R.hcwC"]]
  I.hcwH=Ns[["I.hcwH"]]
  cum_inc_r=Ns[["cum_inc_r"]]
  cum_inc_hcw=Ns[["cum_inc_hcw"]]
  cum_inc_vax1=Ns[["cum_inc_vax1"]]
  cum_inc_vax2=Ns[["cum_inc_vax2"]]
  cum_inc_community=Ns[["cum_inc_community"]]
  total=Ns[["total"]]
  rooms=Ns[["rooms"]]
  vax_count1=Ns[["vax_count1"]]
  vax_count2=Ns[["vax_count2"]]
  t_last_case=Ns[["t_last_case"]]
  

  ## update testing strategy (check if we are using outbreak testing)
  if(!is.na(t_last_case)){
    if(t < (t_last_case + 14)){
      
      VL.threshold.r <- 5    # threshold for detectable VL (put arbitrary # in for now)
      VL.threshold.hcw <- 5  
      test_freq_hcw <- 3 # staff testing frequency (days)
      test_freq_r <- 3  # resident testing frequency (days)
      test_delay_hcw <- 0 # staff test delay (days)
      test_delay_r <-  0  # resident test delay (days)
      
      testing <- "outbreak"
      
    }else{ 
      
      VL.threshold.r <- parms[["VL.threshold.r"]]    # threshold for detectable VL (put arbitrary # in for now)
      VL.threshold.hcw <- parms[["VL.threshold.hcw"]]  
      test_freq_hcw <- parms[["test_freq_hcw"]] # staff testing frequency (days)
      test_freq_r <- parms[["test_freq_r"]]  # resident testing frequency (days)
      test_delay_hcw <- parms[["test_delay_hcw"]] # staff test delay (days)
      test_delay_r <- parms[["test_delay_r"]]  # resident test delay (days)
      
      testing <- "non-outbreak"
      
    }
  } else{
    testing <- "non-outbreak"
  }
  
  testing_parameters <- list("VL.threshold.r" = VL.threshold.r,
                             "VL.threshold.hcw" = VL.threshold.hcw,
                             "test_freq_hcw" = test_freq_hcw,
                             "test_freq_r" = test_freq_r,
                             "test_delay_hcw" = test_delay_hcw,
                             "test_delay_r" = test_delay_r,
                             "testing" = testing)
  
  #cat(t, t_last_case, test_freq_hcw, "\n")

  N.rNC <- nrow(S.rNC) + nrow(E.rNC) + nrow(A.rNC) + nrow(I.rNC) + nrow(R.rNC)
  N.rC <- ifelse(nrow(I.rC) + nrow(R.rC) > 0, nrow(I.rC) + nrow(R.rC), 1)
  N.rNC + N.rC

  N.hcwNC <-  nrow(S.hcwNC) + nrow(E.hcwNC) + nrow(A.hcwNC) + nrow(I.hcwNC) + nrow(R.hcwNC)
  N.hcwC <- nrow(S.hcwC) + nrow(E.hcwC) + nrow(A.hcwC) + nrow(I.hcwC) + nrow(R.hcwC)
  N.hcwNC + N.hcwC

  N.rC <- nrow(I.rC) + nrow(R.rC)
  #cat("pre-vax",N.rNC + N.rC,"\n")
  ## vaccinate residents
  vaccinate.S.rNC <- vaccinate(S.rNC, parms, "residents", vax_count1, vax_count2, t)
  S.rNC <- vaccinate.S.rNC[[1]]
  
  #cat("post-vax",N.rNC + N.rC,"\n")

  ## vaccinate staff
  vaccinate.S.hcwNC <- vaccinate(S.hcwNC, parms, "staff", vax_count1, vax_count2, t)
  S.hcwNC <- vaccinate.S.hcwNC[[1]]
  vax_count1 <- vaccinate.S.hcwNC[[2]]
  vax_count2 <- vaccinate.S.hcwNC[[3]]

  vaccinate.E.hcwNC <- vaccinate(E.hcwNC, parms, "staff", vax_count1, vax_count2, t)
  E.hcwNC <- vaccinate.E.hcwNC[[1]]
  vax_count1 <- vaccinate.E.hcwNC[[2]]
  vax_count2 <- vaccinate.E.hcwNC[[3]]

  vaccinate.A.hcwNC <- vaccinate(A.hcwNC, parms, "staff", vax_count1, vax_count2, t)
  A.hcwNC <- vaccinate.A.hcwNC[[1]]
  vax_count1 <- vaccinate.A.hcwNC[[2]]
  vax_count2 <- vaccinate.A.hcwNC[[3]]

  vaccinate.I.hcwNC <- vaccinate(I.hcwNC, parms, "staff", vax_count1, vax_count2, t)
  I.hcwNC <- vaccinate.I.hcwNC[[1]]
  vax_count1 <- vaccinate.I.hcwNC[[2]]
  vax_count2 <- vaccinate.I.hcwNC[[3]]

  vaccinate.R.hcwNC <- vaccinate(R.hcwNC, parms, "staff", vax_count1, vax_count2, t)
  R.hcwNC <- vaccinate.R.hcwNC[[1]]
  vax_count1 <- vaccinate.R.hcwNC[[2]]
  vax_count2 <- vaccinate.R.hcwNC[[3]]

  vaccinate.S.hcwC <- vaccinate(S.hcwC, parms, "staff", vax_count1, vax_count2, t)
  S.hcwC <- vaccinate.S.hcwC[[1]]
  vax_count1 <- vaccinate.S.hcwC[[2]]
  vax_count2 <- vaccinate.S.hcwC[[3]]

  vaccinate.E.hcwC <- vaccinate(E.hcwC, parms, "staff", vax_count1, vax_count2, t)
  E.hcwC <- vaccinate.E.hcwC[[1]]
  vax_count1 <- vaccinate.E.hcwC[[2]]
  vax_count2 <- vaccinate.E.hcwC[[3]]

  vaccinate.A.hcwC <- vaccinate(A.hcwC, parms, "staff", vax_count1, vax_count2, t)
  A.hcwC <- vaccinate.A.hcwC[[1]]
  vax_count1 <- vaccinate.A.hcwC[[2]]
  vax_count2 <- vaccinate.A.hcwC[[3]]

  vaccinate.I.hcwC <- vaccinate(I.hcwC, parms, "staff", vax_count1, vax_count2, t)
  I.hcwC <- vaccinate.I.hcwC[[1]]
  vax_count1 <- vaccinate.I.hcwC[[2]]
  vax_count2 <- vaccinate.I.hcwC[[3]]

  vaccinate.R.hcwC <- vaccinate(R.hcwC, parms, "staff", vax_count1, vax_count2, t)
  R.hcwC <- vaccinate.R.hcwC[[1]]
  vax_count1 <- vaccinate.R.hcwC[[2]]
  vax_count2 <- vaccinate.R.hcwC[[3]]
  
  res_coverage1 = sum(S.rNC %>% subset(V_doses == 1) %>% nrow(),
                      E.rNC %>% subset(V_doses == 1) %>% nrow(),
                      A.rNC %>% subset(V_doses == 1) %>% nrow(),
                      I.rNC %>% subset(V_doses == 1) %>% nrow(),
                      R.rNC %>% subset(V_doses == 1) %>% nrow(),
                      I.rC %>% subset(V_doses == 1) %>% nrow(),
                      R.rC %>% subset(V_doses == 1) %>% nrow())/100
  
  res_coverage2 = sum(S.rNC %>% subset(V_doses == 2) %>% nrow(),
                      E.rNC %>% subset(V_doses == 2) %>% nrow(),
                      A.rNC %>% subset(V_doses == 2) %>% nrow(),
                      I.rNC %>% subset(V_doses == 2) %>% nrow(),
                      R.rNC %>% subset(V_doses == 2) %>% nrow(),
                      I.rC %>% subset(V_doses == 2) %>% nrow(),
                      R.rC %>% subset(V_doses == 2) %>% nrow())/100

  
   
  # move residents from A/I to R and NC to C
  recover.A.rNC <- recover_or_test_r(A.rNC,parms, "A", testing_parameters)
  R.rNC <- rbind(R.rNC,recover.A.rNC[[1]])    # move A.rNC to R.rNC if recovered
  I.rC <- rbind(I.rC,recover.A.rNC[[2]])      # move A.rNC to I.rC if identified as positive
  A.rNC <- recover.A.rNC[[3]]                 # keep rest of A.rNC in A.rNC

  rooms %>%
    mutate(Res1 = case_when(Res1 %in% recover.A.rNC[[2]]$ID ~ NA_real_, TRUE~Res1),
           Res2 = case_when(Res2 %in% recover.A.rNC[[2]]$ID ~ NA_real_, TRUE~Res2)) -> rooms


  recover.I.rNC <- recover_or_test_r(I.rNC,parms, "I", testing_parameters)
  R.rNC <- rbind(R.rNC,recover.I.rNC[[1]])     # move I.rNC to R.rNC if recovered
  I.rC <- rbind(I.rC,recover.I.rNC[[2]])      # move I.rNC to I.rC if identified as posiitive
  I.rNC <- recover.I.rNC[[3]]                 # keep rest of I.rNC in I.rNC

  rooms %>%
    mutate(Res1 = case_when(Res1 %in% recover.I.rNC[[2]]$ID ~ NA_real_, TRUE~Res1),
           Res2 = case_when(Res2 %in% recover.I.rNC[[2]]$ID ~ NA_real_, TRUE~Res2)) -> rooms

  recover.I.rC <- recover_I.rC(I.rC, parms)

  if (int.r==1 & t>int_time){
    R.rNC <- rbind(R.rNC,recover.I.rC[[1]]) # recover I.rC
    R.room <- recover.I.rC[[1]]$ID          # need room assignment
  } else {
    R.rC <- rbind(R.rC,recover.I.rC[[1]])   # recover I.rC
    R.room <- NULL
  }
  I.rC <- recover.I.rC[[2]] # keep rest in I.rC

  #cat("recover",(nrow(S.rNC) + nrow(E.rNC) + nrow(A.rNC) + nrow(I.rNC) + nrow(R.rNC) + nrow(I.rC) + nrow(R.rC)),t,"\n")
  
  # ratio of staff to residents
  ratio <- (N.hcwNC + N.hcwC)/(N.rNC + N.rC)

  # move hcw from A/I to R and NC or C to home
  recover.A.hcwNC <- recover_or_test_hcw(A.hcwNC,parms,total, "A", ratio, testing_parameters)
  R.hcwNC <- rbind(R.hcwNC,recover.A.hcwNC[[1]]) # move A.hcwNC to R.hcwNC if recover
  I.hcwH <- rbind(I.hcwH,recover.A.hcwNC[[2]]) # move A.hcwNC to I.hcwH if test positive
  A.hcwNC <- recover.A.hcwNC[[3]] # keep rest of A.hcwNC in A.hcwNC
  E.hcwNC <- rbind(E.hcwNC,recover.A.hcwNC[[4]])
  S.hcwNC <- rbind(S.hcwNC,recover.A.hcwNC[[5]])
  R.hcwNC <- rbind(R.hcwNC,recover.A.hcwNC[[6]])
  total <- total + nrow(recover.A.hcwNC[[4]]) + nrow(recover.A.hcwNC[[5]]) + nrow(recover.A.hcwNC[[6]])

  recover.I.hcwNC <- recover_or_test_hcw(I.hcwNC,parms,total, "I", ratio, testing_parameters)
  R.hcwNC <- rbind(R.hcwNC,recover.I.hcwNC[[1]]) # move I.hcwNC to R.hcwNC if recover
  I.hcwH <- rbind(I.hcwH,recover.I.hcwNC[[2]]) # move I.hcwNC to I.hcwH if test positive
  I.hcwNC <- recover.I.hcwNC[[3]] # keep rest of I.hcwNC in I.hcwNC
  E.hcwNC <- rbind(E.hcwNC,recover.I.hcwNC[[4]])
  S.hcwNC <- rbind(S.hcwNC,recover.I.hcwNC[[5]])
  R.hcwNC <- rbind(R.hcwNC,recover.I.hcwNC[[6]])
  total <- total + nrow(recover.I.hcwNC[[4]]) + nrow(recover.I.hcwNC[[5]]) + nrow(recover.I.hcwNC[[6]])

  recover.A.hcwC <- recover_or_test_hcw(A.hcwC,parms,total, "A", ratio, testing_parameters)
  R.hcwC <- rbind(R.hcwC,recover.A.hcwC[[1]]) # move A.hcwC to R.hcwC if recover
  I.hcwH <- rbind(I.hcwH,recover.A.hcwC[[2]]) # move A.hcwC to I.hcwH if test positive
  A.hcwC <- recover.A.hcwC[[3]] # keep rest of A.hcwC in A.hcwC
  E.hcwC <- rbind(E.hcwC,recover.A.hcwC[[4]])
  S.hcwC <- rbind(S.hcwC,recover.A.hcwC[[5]])
  R.hcwC <- rbind(R.hcwC,recover.A.hcwC[[6]])
  total <- total + nrow(recover.A.hcwC[[4]]) + nrow(recover.A.hcwC[[5]]) + nrow(recover.A.hcwC[[6]])

  recover.I.hcwC <- recover_or_test_hcw(I.hcwC,parms,total, "I", ratio, testing_parameters)
  R.hcwC <- rbind(R.hcwC,recover.I.hcwC[[1]]) # move I.hcwC to R.hcwC if recover
  I.hcwH <- rbind(I.hcwH,recover.I.hcwC[[2]]) # move I.hcwC to I.hcwH if test positive
  I.hcwC <- recover.I.hcwC[[3]] # keep rest of I.hcwC in I.hcwC
  E.hcwC <- rbind(E.hcwC,recover.I.hcwC[[4]])
  S.hcwC <- rbind(S.hcwC,recover.I.hcwC[[5]])
  R.hcwC <- rbind(R.hcwC,recover.I.hcwC[[6]])
  total <- total + nrow(recover.I.hcwC[[4]]) + nrow(recover.I.hcwC[[5]]) + nrow(recover.I.hcwC[[6]])

  
  ## check if there were any new infections identified today
  
  if(nrow(recover.A.rNC[[2]])>0 | nrow(recover.I.rNC[[2]])>0 | nrow(recover.A.hcwNC[[2]])>0 | nrow(recover.A.hcwC[[2]])>0 | nrow(recover.I.hcwNC[[2]])>0 | nrow(recover.I.hcwC[[2]])>0){
    t_last_case <- t
  }
  
  # move infected hcw that are home back to one of Rs depending on intervention
  recover.home <- recover_home_hcw(I.hcwH)
  if (int.hcw==1 & t>int_time){
    R.hcwNC <- rbind(R.hcwNC,recover.home[[1]])
    # remove temporary hcw
    if (nrow(recover.home[[1]])>0){
      num_remove <- nrow(recover.home[[1]])
      temp <- unlist(c(S.hcwNC %>% subset(ID>10000) %>% dplyr::select(ID),
                       S.hcwC %>% subset(ID>10000) %>% dplyr::select(ID),
                       E.hcwNC %>% subset(ID>10000) %>% dplyr::select(ID),
                       E.hcwC %>% subset(ID>10000) %>% dplyr::select(ID),
                       A.hcwNC %>% subset(ID>10000) %>% dplyr::select(ID),
                       A.hcwC %>% subset(ID>10000) %>% dplyr::select(ID),
                       I.hcwNC %>% subset(ID>10000) %>% dplyr::select(ID),
                       I.hcwC %>% subset(ID>10000) %>% dplyr::select(ID),
                       R.hcwNC %>% subset(ID>10000) %>% dplyr::select(ID),
                       R.hcwC %>% subset(ID>10000) %>% dplyr::select(ID)), use.names=FALSE)
      num_remove <- min(num_remove,length(temp))
      if (length(temp)>0){
        if (length(temp)>1){
          temp_remove <- cbind("ID"=sample(temp,num_remove,replace=FALSE))
        } else{
          temp_remove <- temp
        }
        S.hcwNC %>%
          subset(!(ID %in% temp_remove)) -> S.hcwNC
        S.hcwC %>%
          subset(!(ID %in% temp_remove)) -> S.hcwC
        E.hcwNC %>%
          subset(!(ID %in% temp_remove)) -> E.hcwNC
        E.hcwC %>%
          subset(!(ID %in% temp_remove)) -> E.hcwC
        A.hcwNC %>%
          subset(!(ID %in% temp_remove)) -> A.hcwNC
        A.hcwC %>%
          subset(!(ID %in% temp_remove)) -> A.hcwC
        I.hcwNC %>%
          subset(!(ID %in% temp_remove)) -> I.hcwNC
        I.hcwC %>%
          subset(!(ID %in% temp_remove)) -> I.hcwC
        R.hcwNC %>%
          subset(!(ID %in% temp_remove)) -> R.hcwNC
        R.hcwC %>%
          subset(!(ID %in% temp_remove)) -> R.hcwC
      }
    }
  } else{
    R.hcwC <- rbind(R.hcwC,recover.home[[1]])
    # remove temporary hcw
    if (nrow(recover.home[[1]])>0){
      num_remove <- nrow(recover.home[[1]])
      temp <- unlist(c(S.hcwNC %>% subset(ID>10000) %>% dplyr::select(ID),
                       S.hcwC %>% subset(ID>10000) %>% dplyr::select(ID),
                       E.hcwNC %>% subset(ID>10000) %>% dplyr::select(ID),
                       E.hcwC %>% subset(ID>10000) %>% dplyr::select(ID),
                       A.hcwNC %>% subset(ID>10000) %>% dplyr::select(ID),
                       A.hcwC %>% subset(ID>10000) %>% dplyr::select(ID),
                       I.hcwNC %>% subset(ID>10000) %>% dplyr::select(ID),
                       I.hcwC %>% subset(ID>10000) %>% dplyr::select(ID),
                       R.hcwNC %>% subset(ID>10000) %>% dplyr::select(ID),
                       R.hcwC %>% subset(ID>10000) %>% dplyr::select(ID)), use.names=FALSE)
      num_remove <- min(num_remove,length(temp))
      if (length(temp)>0){
        if (length(temp)>1){
          temp_remove <- cbind("ID"=sample(temp,num_remove,replace=FALSE))
        } else{
          temp_remove <- temp
        }
        S.hcwNC %>%
          subset(!(ID %in% temp_remove)) -> S.hcwNC
        S.hcwC %>%
          subset(!(ID %in% temp_remove)) -> S.hcwC
        E.hcwNC %>%
          subset(!(ID %in% temp_remove)) -> E.hcwNC
        E.hcwC %>%
          subset(!(ID %in% temp_remove)) -> E.hcwC
        A.hcwNC %>%
          subset(!(ID %in% temp_remove)) -> A.hcwNC
        A.hcwC %>%
          subset(!(ID %in% temp_remove)) -> A.hcwC
        I.hcwNC %>%
          subset(!(ID %in% temp_remove)) -> I.hcwNC
        I.hcwC %>%
          subset(!(ID %in% temp_remove)) -> I.hcwC
        R.hcwNC %>%
          subset(!(ID %in% temp_remove)) -> R.hcwNC
        R.hcwC %>%
          subset(!(ID %in% temp_remove)) -> R.hcwC
      }
    }
  }

  #cat(nrow(recover.home[[1]]),temp,nrow(S.hcwNC),nrow(E.hcwNC),nrow(A.hcwNC),nrow(I.hcwNC),nrow(R.hcwNC),"\n")
  I.hcwH <- recover.home[[2]]

  # moving hcw
  prop.rNC <- (nrow(S.rNC) + nrow(E.rNC) + nrow(A.rNC) + nrow(I.rNC) + nrow(R.rNC))/
    (nrow(S.rNC) + nrow(E.rNC) + nrow(A.rNC) + nrow(I.rNC) + nrow(R.rNC) + nrow(I.rC) + nrow(R.rC))
  prop.hcwNC <- (nrow(S.hcwNC) + nrow(E.hcwNC) + nrow(A.hcwNC) +  nrow(I.hcwNC) + nrow(R.hcwNC))/
    (nrow(S.hcwNC) + nrow(E.hcwNC) + nrow(A.hcwNC) +  nrow(I.hcwNC) + nrow(R.hcwNC) +
       nrow(S.hcwC) + nrow(E.hcwC) + nrow(A.hcwC) + nrow(I.hcwC) + nrow(R.hcwC))

  moved_hcw <- move_hcw(S.hcwNC,
                        E.hcwNC,
                        A.hcwNC,
                        I.hcwNC,
                        R.hcwNC,
                        S.hcwC,
                        E.hcwC,
                        A.hcwC,
                        I.hcwC,
                        R.hcwC,
                        prop.rNC,
                        prop.hcwNC)
  S.hcwNC <- moved_hcw[[1]]
  E.hcwNC <- moved_hcw[[2]]
  A.hcwNC <- moved_hcw[[3]]
  I.hcwNC <- moved_hcw[[4]]
  R.hcwNC <- moved_hcw[[5]]
  S.hcwC <- moved_hcw[[6]]
  E.hcwC <- moved_hcw[[7]]
  A.hcwC <- moved_hcw[[8]]
  I.hcwC <- moved_hcw[[9]]
  R.hcwC <- moved_hcw[[10]]

  # E -> A/I
  infect_E.rNC <- E_to_I_r(E.rNC,parms,testing_parameters)
  A.rNC <- bind_rows(A.rNC,infect_E.rNC[[1]])
  I.rNC <- bind_rows(I.rNC,infect_E.rNC[[2]])
  E.rNC <- infect_E.rNC[[3]]
  
  asymptomatic_inc_r_0 <- nrow(infect_E.rNC[[1]] %>% subset(V_doses == 0 | is.na(V_doses)))
  asymptomatic_inc_r_1 <- nrow(infect_E.rNC[[1]] %>% subset(V_doses == 1))
  asymptomatic_inc_r_2 <- nrow(infect_E.rNC[[1]] %>% subset(V_doses == 2))
  
  symptomatic_inc_r_0 <- nrow(infect_E.rNC[[2]] %>% subset(V_doses == 0 | is.na(V_doses)))
  symptomatic_inc_r_1 <- nrow(infect_E.rNC[[2]] %>% subset(V_doses == 1))
  symptomatic_inc_r_2 <- nrow(infect_E.rNC[[2]] %>% subset(V_doses == 2))

  #cat("infect",(nrow(S.rNC) + nrow(E.rNC) + nrow(A.rNC) + nrow(I.rNC) + nrow(R.rNC) + nrow(I.rC) + nrow(R.rC)),t,"\n")

  infect_E.hcwNC <- E_to_I_hcw(E.hcwNC,parms,testing_parameters)
  A.hcwNC <- rbind(A.hcwNC,infect_E.hcwNC[[1]])
  I.hcwNC <- rbind(I.hcwNC,infect_E.hcwNC[[2]])
  E.hcwNC <- infect_E.hcwNC[[3]]

  infect_E.hcwC <- E_to_I_hcw(E.hcwC,parms,testing_parameters)
  A.hcwC <- rbind(A.hcwC,infect_E.hcwC[[1]])
  I.hcwC <- rbind(I.hcwC,infect_E.hcwC[[2]])
  E.hcwC <- infect_E.hcwC[[3]]
  
  asymptomatic_inc_hcw_0 <- nrow(infect_E.hcwNC[[1]] %>% subset(V_doses == 0 | is.na(V_doses)))
  asymptomatic_inc_hcw_1 <- nrow(infect_E.hcwNC[[1]] %>% subset(V_doses == 1))
  asymptomatic_inc_hcw_2 <- nrow(infect_E.hcwNC[[1]] %>% subset(V_doses == 2))
  
  symptomatic_inc_hcw_0 <- nrow(infect_E.hcwNC[[2]] %>% subset(V_doses == 0 | is.na(V_doses)))
  symptomatic_inc_hcw_1 <- nrow(infect_E.hcwNC[[2]] %>% subset(V_doses == 1))
  symptomatic_inc_hcw_2 <- nrow(infect_E.hcwNC[[2]] %>% subset(V_doses == 2))
  
  # update contact rates based on ratio
  k.HR <- k.HR*num_staff/(N.hcwNC + N.hcwC)
  #print(k.HR)

  # S -> E
  #cat("expose",S.rNC.exposed,nrow(S.rNC),nrow(E.rNC),(nrow(S.rNC) + nrow(E.rNC) + nrow(A.rNC) + nrow(I.rNC) + nrow(R.rNC) + nrow(I.rC) + nrow(R.rC)),t,"\n")

  S.rNC2 <- S.rNC

  inc_vax1 <- 0
  inc_vax2 <- 0

  daily_visitors = 0
  
  if(nrow(S.rNC)>=1){
    for (i in 1:nrow(S.rNC)){
      roommate.exp <- expose_roommate(rooms,S.rNC$ID[i],A.rNC,I.rNC)
      
      if(int.vax == 1){
        if(S.rNC$VE_s[i]==VE_s2_r){
          visitor.exp <- rbinom(1,1,visitors)
        }else{visitor.exp <- 0}
      }else{visitor.exp <- rbinom(1,1,visitors)}
      
      daily_visitors = (daily_visitors + visitor.exp)
      
      visitor.infected <- rbinom(1,visitor.exp,I.C*10)
      
      # adjust betas based on vaccination status
      beta.rm.v <- beta.rm*(1-S.rNC$VE_s[i])
      beta.s.v <- beta.s*(1-S.rNC$VE_s[i])
      beta.v <- parms[["beta"]]*(1-S.rNC$VE_s[i]) # reduce susceptibility by VEs (assuming leaky here)

      exp.prob <- beta.rm.v*roommate.exp + 
        beta.rm.v*visitor.infected + 
        beta.v*(ifelse(N.rNC>0,k.RR*(sum(I.rNC$Infectiousness*(1-I.rNC$VE_i))+sum(A.rNC$Infectiousness*(1-A.rNC$VE_i)))/N.rNC,0)) +
        beta.s.v*(ifelse(N.hcwNC>0,k.RH*(sum(I.hcwNC$Infectiousness*(1-I.hcwNC$VE_i))+sum(A.hcwNC$Infectiousness*(1-A.hcwNC$VE_i)))/N.hcwNC,0))
      
      if (exp.prob >= 1){
        exp.prob <- 1
      } else if (exp.prob <0){
        exp.prob <- 0
      }

      S.rNC.exposed <- rbinom(1,1,exp.prob)
      #cat(i,S.rNC.exposed,"\n")
      if (S.rNC.exposed==1){
        E.rNC <- rbind(E.rNC,cbind(ID=S.rNC$ID[i],VL=0,arrival=S.rNC$arrival[i],departure=S.rNC$departure[i],
                                       V_doses = S.rNC$V_doses[i],V_type = S.rNC$V_type[i],V_refused = S.rNC$V_refused[i],V1_t = S.rNC$V1_t[i],V2_t = S.rNC$V2_t[i],VE_i = S.rNC$VE_i[i],VE_s = S.rNC$VE_s[i],VE_p = S.rNC$VE_p[i],
                                       Inc.pd=round(runif(1, parms[["sigma1"]], parms[["sigma2"]])), Days=0))

        ## update incidence counter
        if(!is.na(S.rNC$V1_t[i])){
          if(t>=(S.rNC$V1_t[i] + parms[["immune_delay"]]) & t < (S.rNC$V2_t[i] + parms[["immune_delay"]])){
            inc_vax1 = inc_vax1 + 1
          }else if (t> (S.rNC$V1_t[i] + parms[["immune_delay"]])){
            inc_vax2 = inc_vax2 + 1
          }
        }

        ## remove from Susceptibles
        S.rNC2 %>%
          subset(ID != S.rNC$ID[i]) -> S.rNC2

      }
    }
  }

  expose.S.rNC <- nrow(S.rNC) - nrow(S.rNC2)
  S.rNC <- S.rNC2


  #cat("expose",S.rNC.exposed,nrow(S.rNC),nrow(E.rNC),(nrow(S.rNC) + nrow(E.rNC) + nrow(A.rNC) + nrow(I.rNC) + nrow(R.rNC) + nrow(I.rC) + nrow(R.rC)),t,"\n")

  ## community introductions
  # S.hcwNC.exposed <- rbinom(1,nrow(S.hcwNC), I.C)
  # expose.S.hcwNC <- infect(S.hcwNC,parms,S.hcwNC.exposed)
  # E.hcwNC <- rbind(E.hcwNC,expose.S.hcwNC[[1]])
  # S.hcwNC <- expose.S.hcwNC[[2]]
  #
  # S.hcwC.exposed <- rbinom(1,nrow(S.hcwC), I.C)
  # expose.S.hcwC <- infect(S.hcwC,parms,S.hcwC.exposed)
  # E.hcwC <- rbind(E.hcwC,expose.S.hcwC[[1]])
  # S.hcwC <- expose.S.hcwC[[2]]

  ## tally total community introductions
  #inc_community <- nrow(expose.S.hcwNC[[1]]) + nrow(expose.S.hcwC[[1]])


  ## staff transmission (community and nursing homes)

  inc_community <- 0


  S.hcwNC2 <- S.hcwNC

  if(nrow(S.hcwNC2)>=1){
    for (i in 1:nrow(S.hcwNC)){

      ## first see if infected in community introductions
      S.hcwNC.exposed <- rbinom(1,1, I.C*(1-S.hcwNC$VE_s[i])) # could also combine this with exp.prob below?

      if (S.hcwNC.exposed != 1){ # if not infected in community, check if infected in nursing home

        # adjust betas based on vaccination status
        beta.s.v <- beta.s*(1-S.hcwNC$VE_s[i])
        beta.v <- beta*(1-S.hcwNC$VE_s[i]) # reduce susceptibility by VEs (assuming leaky here)

        exp.prob <- beta.v*(ifelse(N.hcwNC>0,k.HH*(sum(I.hcwNC$Infectiousness*(1-I.hcwNC$VE_i))+sum(A.hcwNC$Infectiousness*(1-A.hcwNC$VE_i)))/N.hcwNC,0)) +
          beta.s.v*(ifelse(N.rNC>0,k.HR*(sum(I.rNC$Infectiousness*(1-I.rNC$VE_i))+sum(A.rNC$Infectiousness*(1-A.rNC$VE_i)))/N.rNC,0))

        if (exp.prob >= 1){
          exp.prob <- 1
        } else if (exp.prob <0){
          exp.prob <- 0
        }

        S.hcwNC.exposed <- rbinom(1,1,exp.prob)

      } else{
        inc_community <- inc_community + 1 # add to comm intro tracker
      }

      if (S.hcwNC.exposed==1){
        E.hcwNC <- rbind(E.hcwNC,cbind(ID=S.hcwNC$ID[i],VL=0,
                                   V_doses = S.hcwNC$V_doses[i],V_type = S.hcwNC$V_type[i],V_refused = S.hcwNC$V_refused[i],V1_t = S.hcwNC$V1_t[i],V2_t = S.hcwNC$V2_t[i],VE_i = S.hcwNC$VE_i[i],VE_s = S.hcwNC$VE_s[i],VE_p = S.hcwNC$VE_p[i],
                                   Inc.pd=round(runif(1, parms[["sigma1"]], parms[["sigma2"]])), Days=0))
        S.hcwNC2 %>%
          subset(ID != S.hcwNC$ID[i]) -> S.hcwNC2
      }
    }
  }

  expose.S.hcwNC <- nrow(S.hcwNC) - nrow(S.hcwNC2)
  S.hcwNC <- S.hcwNC2



  S.hcwC2 <- S.hcwC

  if(nrow(S.hcwC2)>=1){
    for (i in 1:nrow(S.hcwC)){

      ## first see if infected in community introductions
      S.hcwC.exposed <- rbinom(1,1, I.C*(1-S.hcwC$VE_s[i]))

      if (S.hcwC.exposed != 1){ # if not infected in community, check if infected in nursing home

        # adjust betas based on vaccination status
        beta.s.v <- beta.s*(1-S.hcwC$VE_s[i])
        beta.v <- beta*(1-S.hcwC$VE_s[i]) # reduce susceptibility by VEs (assuming leaky here)

        exp.prob <- beta.v*ppe*(ifelse(N.hcwC>0,k.HH*(sum(I.hcwC$Infectiousness*(1-I.hcwC$VE_i))+sum(A.hcwC$Infectiousness*(1-A.hcwC$VE_i)))/N.hcwC,0)) +
          beta.s.v*(ifelse(N.rC>0,k.HR*(sum(I.rC$Infectiousness*(1-I.rC$VE_i)))/N.rNC,0))

        if (exp.prob >= 1){
          exp.prob <- 1
        } else if (exp.prob <0){
          exp.prob <- 0
        }

        S.hcwC.exposed <- rbinom(1,1,exp.prob)

      } else{
        inc_community <- inc_community + 1 # add to comm intro tracker
      }

      if (S.hcwC.exposed==1){
        E.hcwC <- rbind(E.hcwC,cbind(ID=S.hcwC$ID[i],VL=0,
                                       V_doses = S.hcwC$V_doses[i],V_type = S.hcwC$V_type[i],V_refused = S.hcwC$V_refused[i],V1_t = S.hcwC$V1_t[i],V2_t = S.hcwC$V2_t[i],VE_i = S.hcwC$VE_i[i],VE_s = S.hcwC$VE_s[i],VE_p = S.hcwC$VE_p[i],
                                       Inc.pd=round(runif(1, parms[["sigma1"]], parms[["sigma2"]])), Days=0))
        S.hcwC2 %>%
          subset(ID != S.hcwC$ID[i]) -> S.hcwC2
      }
    }
  }

  expose.S.hcwC <- nrow(S.hcwC) - nrow(S.hcwC2)
  S.hcwC <- S.hcwC2



  N.rNC <- nrow(S.rNC) + nrow(E.rNC) + nrow(A.rNC) + nrow(I.rNC) + nrow(R.rNC)
  N.rC <- nrow(I.rC) + nrow(R.rC)
  N.rNC + N.rC
  #cat("pre-death",N.rNC + N.rC,"\n")

  # death (and discharges)
  death.S.rNC <- death(S.rNC,mu.NC,parms,total,t,dt)
  S.rNC <- death.S.rNC[[1]]
  S.rNC <- bind_rows(S.rNC,death.S.rNC[[3]]) # new entry
  total <- total + nrow(death.S.rNC[[3]])
  
  #cat("S.rNC",nrow(death.S.rNC[[3]]),nrow(death.S.rNC[[4]]),"\n")

  death.E.rNC <- death(E.rNC,mu.NC,parms,total,t,dt)
  E.rNC <- death.E.rNC[[1]]
  S.rNC <- bind_rows(S.rNC,death.E.rNC[[3]]) # new entry
  total <- total + nrow(death.E.rNC[[3]])
  
  #cat("E.rNC",nrow(death.E.rNC[[3]]),nrow(death.E.rNC[[4]]),"\n")
  

  death.A.rNC <- death(A.rNC,mu.NC,parms,total,t,dt)
  A.rNC <- death.A.rNC[[1]]
  S.rNC <- bind_rows(S.rNC,death.A.rNC[[3]]) # new entry
  total <- total + nrow(death.A.rNC[[3]])
  
  #cat("A.rNC",nrow(death.A.rNC[[3]]),nrow(death.A.rNC[[4]]),"\n")
  

  death.I.rNC <- death(I.rNC,mu.NC,parms,total,t,dt)
  I.rNC <- death.I.rNC[[1]]
  S.rNC <- bind_rows(S.rNC,death.I.rNC[[3]]) # new entry
  total <- total + nrow(death.I.rNC[[3]])

  #cat("I.rNC",nrow(death.I.rNC[[3]]),nrow(death.I.rNC[[4]]),"\n")
  

  death.R.rNC <- death(R.rNC,mu.NC,parms,total,t,dt)
  R.rNC <- death.R.rNC[[1]]
  S.rNC <- bind_rows(S.rNC,death.R.rNC[[3]]) # new entry
  total <- total + nrow(death.R.rNC[[3]])
  
  #cat("R.rNC",nrow(death.R.rNC[[3]]),nrow(death.R.rNC[[4]]),"\n")
  

  death.I.rC <- covid_death(I.rC,mu.C,parms,total,t,dt)
  I.rC <- death.I.rC[[1]]
  S.rNC <- bind_rows(S.rNC,death.I.rC[[3]]) # new entry
  total <- total + nrow(death.I.rC[[3]])
  
  #cat("I.rC",nrow(death.I.rC[[3]]),nrow(death.I.rC[[4]]),"\n")
  

  death.R.rC <- death(R.rC,mu.NC,parms,total,t,dt)
  R.rC <- death.R.rC[[1]]
  S.rNC <- bind_rows(S.rNC,death.R.rC[[3]]) # new entry
  total <- total + nrow(death.R.rC[[3]])
  
  #cat("R.rC",nrow(death.R.rC[[3]]),nrow(death.R.rC[[4]]),"\n")
  
  
  N.rNC <- nrow(S.rNC) + nrow(E.rNC) + nrow(A.rNC) + nrow(I.rNC) + nrow(R.rNC)
  N.rC <- nrow(I.rC) + nrow(R.rC)
  N.rNC + N.rC
  #cat("post-death",N.rNC + N.rC,"\n")
  

  # update rooms
  new_entry <- c(death.S.rNC[[3]]$ID,death.E.rNC[[3]]$ID,death.A.rNC[[3]]$ID,death.I.rNC[[3]]$ID,death.R.rNC[[3]]$ID,
                 death.I.rC[[3]]$ID,death.R.rC[[3]]$ID)

  new_dead <- c(death.S.rNC[[4]]$ID,death.E.rNC[[4]]$ID,death.A.rNC[[4]]$ID,death.I.rNC[[4]]$ID,death.R.rNC[[4]]$ID,
                death.I.rC[[4]]$ID,death.R.rC[[4]]$ID) # also includes discharges

  rooms %>%
    mutate(Res1 = case_when(Res1 %in% new_dead ~ NA_real_, TRUE~Res1),
           Res2 = case_when(Res2 %in% new_dead ~ NA_real_, TRUE~Res2)) -> rooms

  if (length(c(new_entry,R.room))>0){
    rooms <- assign_rooms(rooms,new_entry,R.room,int.rooms,t, parms)
  }

  #cat(new_entry,new_entry %in% rooms$Res1 | new_entry %in% rooms$Res2,"\n")

  #cat("death",(nrow(S.rNC) + nrow(E.rNC) + nrow(A.rNC) + nrow(I.rNC) + nrow(R.rNC) + nrow(I.rC) + nrow(R.rC)),t,"\n")

  # change VL for recovered
  R.rNC <- recover_VL(R.rNC)
  R.hcwNC <- recover_VL(R.hcwNC)
  R.rC <- recover_VL(R.rC)
  R.hcwC <- recover_VL(R.hcwC)

  mortality <- death.S.rNC[[2]] + death.E.rNC[[2]] + death.A.rNC[[2]] + death.I.rNC[[2]] + death.R.rNC[[2]] +
    death.I.rC[[2]] + death.R.rC[[2]]

  discharges <- death.S.rNC[[5]] + death.E.rNC[[5]] + death.A.rNC[[5]] + death.I.rNC[[5]] + death.R.rNC[[5]] +
    death.I.rC[[5]] + death.R.rC[[5]]

  inc_r <- expose.S.rNC
  inc_hcw <- expose.S.hcwNC + expose.S.hcwC

  cum_inc_r = (cum_inc_r + inc_r)
  cum_inc_hcw = (cum_inc_hcw + inc_hcw)
  cum_inc_vax1 = (cum_inc_vax1 + inc_vax1)
  cum_inc_vax2 = (cum_inc_vax2 + inc_vax2)
  cum_inc_community = (cum_inc_community + inc_community)
  
  #cat(t,res_coverage1,res_coverage2,"\n")

  final <- list("S.rNC"=S.rNC,
                "E.rNC"=E.rNC,
                "A.rNC"=A.rNC,
                "I.rNC"=I.rNC,
                "R.rNC"=R.rNC,
                "I.rC"=I.rC,
                "R.rC"=R.rC,
                "S.hcwNC"=S.hcwNC,
                "E.hcwNC"=E.hcwNC,
                "A.hcwNC"=A.hcwNC,
                "I.hcwNC"=I.hcwNC,
                "R.hcwNC"=R.hcwNC,
                "S.hcwC"=S.hcwC,
                "E.hcwC"=E.hcwC,
                "A.hcwC"=A.hcwC,
                "I.hcwC"=I.hcwC,
                "R.hcwC"=R.hcwC,
                "I.hcwH"=I.hcwH,
                "inc_r"=inc_r,
                "inc_hcw"=inc_hcw,
                "cum_inc_r"=cum_inc_r,
                "cum_inc_hcw"=cum_inc_hcw,
                "cum_inc_vax1"=cum_inc_vax1,
                "cum_inc_vax2"=cum_inc_vax2,
                "cum_inc_community"=cum_inc_community,
                "mortality"=mortality,
                "discharges"=discharges,
                "total"=total,
                "rooms"=rooms,
                "vax_count1"=vax_count1,
                "vax_count2"=vax_count2, 
                "asymptomatic_inc_r_0"=asymptomatic_inc_r_0,
                "asymptomatic_inc_r_1"=asymptomatic_inc_r_1,
                "asymptomatic_inc_r_2"=asymptomatic_inc_r_2,
                "symptomatic_inc_r_0"=symptomatic_inc_r_0,
                "symptomatic_inc_r_1"=symptomatic_inc_r_1,
                "symptomatic_inc_r_2"=symptomatic_inc_r_2,
                "asymptomatic_inc_hcw_0"=asymptomatic_inc_hcw_0,
                "asymptomatic_inc_hcw_1"=asymptomatic_inc_hcw_1,
                "asymptomatic_inc_hcw_2"=asymptomatic_inc_hcw_2,
                "symptomatic_inc_hcw_0"=symptomatic_inc_hcw_0,
                "symptomatic_inc_hcw_1"=symptomatic_inc_hcw_1,
                "symptomatic_inc_hcw_2"=symptomatic_inc_hcw_2,
                "daily_visitors"=daily_visitors,
                "res_coverage1"=res_coverage1,
                "res_coverage2"=res_coverage2,
                "t_last_case"=t_last_case
  )
  
  return(final)

}

inits=c(S.rNC.init = 100,
        E.rNC.init = 0,
        A.rNC.init = 0,
        I.rNC.init = 0,
        R.rNC.init = 0,
        I.rC.init  = 0,
        R.rC.init = 0,
        S.hcwNC.init = 99,
        E.hcwNC.init = 0,
        A.hcwNC.init = 0,
        I.hcwNC.init = 0,
        R.hcwNC.init = 0,
        S.hcwC.init = 1,
        E.hcwC.init = 0,
        A.hcwC.init = 0,
        I.hcwC.init = 0,
        R.hcwC.init = 0,
        I.hcwH.init = 0 #,
        # inc_r=0,
        # inc_hcw=0,
        # cum_inc_r=0,
        # cum_inc_hcw=0,
        # cum_inc_community=0,
        # mortality=0,
        # discharges = 0,
        # total=0,
        # vax_count1=0,
        # vax_count2=0
        )

num_staff <- sum(inits[8:18])
ratio <- sum(inits[8:18])/sum(inits[1:7])

args=(commandArgs(TRUE))

for (i in 1:length(args)) {
  eval (parse (text = args[[i]] ))
}

nsim <- 1
j <- j

### run simulations ----
parms <- list(beta=0.0372,    ## beta
           beta.s=(1-(1-0.0372)**3), ## beta for interactions between staff and residents
           beta.rm=(1-(1-0.0372)**10),  ## beta for roommates
           sigma1=1.5,      ## incubation period shortest
           sigma2=4.49,      ## incubation period longest
           gamma=14,     ## infectious period
           I.C=0.03,    ## probability of infection from community -- community prevalence?
           k.HH=2,       ## n contacts between staff
           k.RH=6,       ## n staff contacted by each resident, per day
           k.HR=6,       ## n residents contacted by each staff member, per day
           k.RR=6,       ## n daily contacts between residents besides roommates
           alpha.hcw = 0.4, # proportion hcw asymptomatic
           alpha.r = 0.2,   # proportion residents asymptomatic
           id.I= 2,         # duration of pre-symptomatic transmission
           ppe=0.05,        # reduction in beta with ppe
           prop_rhcwR=0.2,    # probability replacement hcw recovered
           mu.C=0.02,       # COVID mortality
           mu.NC=0,     # regular mortality
           int.r=1,         # 1 if go to new_recovered.rNC, 0 if go to new_recovered.rC
           int.rooms=1,     # if 1 put susceptible and recovered together in NC
           int.hcw=1,       # 1 if recovered work with NC and 0 if work with C
           VL.threshold.r=0,  # threshold for detectable VL (put arbitrary # in for now)
           VL.threshold.hcw=0,
           test_freq_hcw=7, # staff testing frequency (days)
           test_freq_r=7,   # resident testing frequency (days)
           test_delay_hcw=2, # staff test delay (days)
           test_delay_r=2,   # resident test delay (days)
           int_time=1,
           I.C_time=16,  ## 37 if waiting till after second dose
           VE_s1_r1 = 0.35, # mRNA
           VE_p1_r1 = 0.16,
           VE_i1_r1 = 0.9,
           VE_s1_hcw1 = 0.35,
           VE_p1_hcw1 = 0.16,
           VE_i1_hcw1 = 0.9,
           VE_s2_r1 = 0.7,
           VE_p2_r1 = 2/3,
           VE_i2_r1 = 0.9,
           VE_s2_hcw1 = 0.7,
           VE_p2_hcw1 = 2/3,
           VE_i2_hcw1 = 0.9,
           VE_s1_r2 = 0.72, # J&J
           VE_p1_r2 = 0,
           VE_i1_r2 = 0.9,
           VE_s1_hcw2 = 0.72,
           VE_p1_hcw2 = 0,
           VE_i1_hcw2 = 0.9,
           VE_s2_r2 = 0.72,
           VE_p2_r2 = 0,
           VE_i2_r2 = 0.9,
           VE_s2_hcw2 = 0.72,
           VE_p2_hcw2 = 0,
           VE_i2_hcw2 = 0.9,
           vax1_prob = 1,
           vax1_prob_admit = 1,
           vax_gap = 0,
           t_first_dose = c(1),
           #pr_each_day = 1/length(t_first_dose),
           second_dose_delay = 21,
           immune_delay = 14,
           vax_coverage_r = 0.8,
           vax_coverage_hcw = 0.5,
           vax_available=1000000,
           refusal_hcw= 0,
           refusal_r = 0,
           shortage = TRUE,
           shortage_threshold = 0.50*ratio,
           comm.vax = 0.8, # kept it same as coverage_r for consistency 
           visitors = 0.1,  ## c(0, 0.3, 0.6)
           int.vax = 0,     ## if 1, intervention relaxation only applies to vaccinated
           wane_day_r = 1000, # day waning starts for residents
           wane_rate_r = 0.005, # daily waning rate for residents 
           wane_day_hcw = 1000, # day waning starts for residents
           wane_rate_hcw = 0.001 # daily waning rate for residents 
  )

t_step <- 1

dt <- seq(1,76,t_step)  ## 2 months after vaccination
#dt <- seq(0,180,t_step)



## loop over testing strategies

testing_scenarios <- as.data.frame(cbind(c("PCR", "PCR", "Antigen","PCR","PCR","Antigen","Antigen","PCR","PCR","None",
                                           "PCR","Antigen","Antigen","None","None","None","None"),
                                         c("PCR", "Antigen", "Antigen","PCR","PCR","Antigen","Antigen","PCR","PCR", "None",
                                           "PCR","Antigen", "Antigen","Antigen", "Antigen","Antigen", "Antigen")))
rownames(testing_scenarios) <- c("weekly_PCR", "daily_staff_antigen", "daily_antigen","weekly_PCR_slow","twice_weekly_PCR",
                                 "weekly_antigen","antigen_highLOD","daily_PCR", "fast_PCR","none",
                                 "2 day both PCR","2 day both","3 day both",
                                 "1 day staff", "2 day staff","3 day staff", "7 day staff")
colnames(testing_scenarios) <- c("Res", "HCW")
testing_scenarios$LOD.r <- c(3,3,5,3,3,5,7,3,3,1000,3,5,5,1000,1000,1000,1000)
testing_scenarios$freq.r <- c(7,7,1,7,3,7,1,1,7,1000,2,2,3,1000,1000,1000,1000)
testing_scenarios$delay.r <- c(2,2,0,7,2,0,0,2,1,1000,2,0,0,1000,1000,1000,1000)
testing_scenarios$LOD.hcw <- c(3,5,5,3,3,5,7,3,3,1000,3,5,5,5,5,5,5)
testing_scenarios$freq.hcw <- c(7,1,1,7,3,7,1,1,7,1000,2,2,3,1,2,3,7)
testing_scenarios$delay.hcw <- c(2,0,0,7,2,0,0,2,1,1000,2,0,0,0,0,0,0)

testing_scenarios <- testing_scenarios[c(16),] # subset to testing scenarios
#testing_scenarios <- testing_scenarios[c(10),] # subset to testing scenarios




#for (j in j_list){
  
res_master <- NULL
VL_master <- NULL

#IC_list <- c(0.0001,0.001,0.003)   # did 3 since not sure we'll see much difference with low or medium

IC_list <- c(0.001, 0.003)    


for(VE_s2_r1 in c(0.5)){ #0.7
  
  if (VE_s2_r1 == 0.5){
    parms[["VE_s1_r1"]] = VE_s2_r1
    parms[["VE_s1_hcw1"]] = 0.7
    parms[["VE_s2_hcw1"]] = 0.7
    
    parms[["VE_p1_r1"]] = 0.6
    parms[["VE_p2_r1"]] = 0.6
    parms[["VE_p1_hcw1"]] = 2/3
    parms[["VE_p2_hcw1"]] = 2/3
    
    parms[["VE_s1_r2"]] = VE_s2_r1
    parms[["VE_s1_hcw2"]] = 0.7
    parms[["VE_s2_hcw2"]] = 0.7
    
    parms[["VE_p1_r2"]] = 0.6
    parms[["VE_p2_r2"]] = 0.6
    parms[["VE_p1_hcw2"]] = 2/3
    parms[["VE_p2_hcw2"]] = 2/3
    
    boost_list <- c(0.5, 0.6, 0.7, 0.8, 0.9)
    
    
  } else{
    parms[["VE_s1_r1"]] = VE_s2_r1
    parms[["VE_s1_hcw1"]] = VE_s2_r1
    parms[["VE_s2_hcw1"]] = VE_s2_r1
    
    parms[["VE_p1_r1"]] = 2/3
    parms[["VE_p2_r1"]] = 2/3
    parms[["VE_p1_hcw1"]] = 2/3
    parms[["VE_p2_hcw1"]] = 2/3
    
    parms[["VE_s1_r2"]] = VE_s2_r1
    parms[["VE_s1_hcw2"]] = VE_s2_r1
    parms[["VE_s2_hcw2"]] = VE_s2_r1
    
    parms[["VE_p1_r2"]] = 2/3
    parms[["VE_p2_r2"]] = 2/3
    parms[["VE_p1_hcw2"]] = 2/3
    parms[["VE_p2_hcw2"]] = 2/3
    
    boost_list <- c(0.7, 0.8, 0.9)
      
  }
  
  for(boost in boost_list){
    
    parms[["VE_s2_r1"]] =  boost
    parms[["VE_s2_r2"]] =  boost

    for(ic in IC_list){
      
      parms[["I.C"]] = ic
      

        for(testing in 1:nrow(testing_scenarios)){
          
          test_scenario <- rownames(testing_scenarios)[testing]
          parms[["VL.threshold.r"]] <-  testing_scenarios[test_scenario,"LOD.r"]
          parms[["test_freq_r"]] <-  testing_scenarios[test_scenario,"freq.r"]
          parms[["test_delay_r"]] <-  testing_scenarios[test_scenario,"delay.r"]
          
          parms[["VL.threshold.hcw"]] <- testing_scenarios[test_scenario,"LOD.hcw"]
          parms[["test_freq_hcw"]] <- testing_scenarios[test_scenario,"freq.hcw"]
          parms[["test_delay_hcw"]] <- testing_scenarios[test_scenario,"delay.hcw"]
       
            
            for(vax_coverage_hcw in c(0.4,0.5,0.6,0.7,0.8,0.9,1)){ #0.40,0.5,0.6,0.7,0.8
              
              
              parms[["vax_coverage_hcw"]] <- vax_coverage_hcw
              
                vacc_scenario <- NA
                VE <- NA
                
                parms[["int.r"]] <- 1
                parms[["int.rooms"]] <- 1
                parms[["int.hcw"]] <- 0
                Intervention <- "Resident Cohorting"
          
                for (sim in 1:nsim){
                  cat(j,VE_s2_r1,ic,vax_coverage_hcw,"\n")
                  
                  Ns <- initialize(inits,parms,dt)
                  res <- as.data.frame(matrix(nrow=length(dt),ncol=length(Ns)))
                  VLs <- NULL
                  bug_staff <- NULL
                  bug_resident <- NULL
                  bug_r <-0
                  bug_s <- 0
                  for(i in 1:length(dt)){
                    #cat(i,"\n")
                    debug <- Ns
                    final <- stochastic_NH(parms,Ns,t_step,(i)*t_step,num_staff)
                    Ns <- final
                    res[i,] <- c(i, nrow(Ns[["S.rNC"]]), nrow(Ns[["E.rNC"]]), nrow(Ns[["A.rNC"]]), nrow(Ns[["I.rNC"]]), nrow(Ns[["R.rNC"]]),
                                 nrow(Ns[["I.rC"]]), nrow(Ns[["R.rC"]]),
                                 nrow(Ns[["S.hcwNC"]]), nrow(Ns[["E.hcwNC"]]), nrow(Ns[["A.hcwNC"]]), nrow(Ns[["I.hcwNC"]]), nrow(Ns[["R.hcwNC"]]),
                                 nrow(Ns[["S.hcwC"]]), nrow(Ns[["E.hcwC"]]), nrow(Ns[["A.hcwC"]]), nrow(Ns[["I.hcwC"]]), nrow(Ns[["R.hcwC"]]),
                                 nrow(Ns[["I.hcwH"]]), Ns[["inc_r"]],  Ns[["inc_hcw"]], 
                                 Ns[["cum_inc_r"]], Ns[["cum_inc_hcw"]], Ns[["cum_inc_vax1"]], Ns[["cum_inc_vax2"]],Ns[["cum_inc_community"]],
                                 Ns[["mortality"]],Ns[["discharges"]],Ns[["total"]],Ns[["vax_count1"]],Ns[["vax_count2"]], 
                                 Ns[["asymptomatic_inc_r_0"]],Ns[["asymptomatic_inc_r_1"]],Ns[["asymptomatic_inc_r_2"]], 
                                 Ns[["symptomatic_inc_r_0"]],Ns[["symptomatic_inc_r_1"]],Ns[["symptomatic_inc_r_2"]],  
                                 Ns[["asymptomatic_inc_hcw_0"]],Ns[["asymptomatic_inc_hcw_1"]],Ns[["asymptomatic_inc_hcw_2"]], 
                                 Ns[["symptomatic_inc_hcw_0"]],Ns[["symptomatic_inc_hcw_1"]],Ns[["symptomatic_inc_hcw_2"]], 
                                 Ns[["daily_visitors"]],Ns[["res_coverage1"]],Ns[["res_coverage2"]],Ns[["t_last_case"]])
          
                    staff <- sum(res[i,c(9:18)])
                    residents <- sum(res[i,2:8])
                    if (staff!=100 & bug_s==0){
                      bug_staff <- debug
                      #cat("staff",i,staff,"\n")
                      bug_s <- 1
                    }
                    if (residents!=100 & bug_r==0){
                      bug_resident <- debug
                      cat("residents",i,residents,"\n")
                      bug_r <- 1
                    }
                    VLs <- rbind(VLs, get_VL(final, i))
                  }
          
                  
                  res %>%
                    add_column("Sim"=j) %>%
                    add_column("beta"=parms[["beta"]]) %>%
                    add_column("Testing"=test_scenario) %>%
                    add_column("VE2"=VE_s2_r1) %>%
                    add_column("VEstaff"=parms[["VE_s2_hcw1"]]) %>%
                    add_column("vax_coverage_hcw"=parms[["vax_coverage_hcw"]]) %>%
                    add_column("booster"=boost) %>%
                    add_column("VE_i2_r1"=parms[["VE_i2_r1"]]) %>%
                    add_column("I.C" =  parms[["I.C"]]) %>%
                    add_column("visitors" = parms[["visitors"]]) %>%
                    add_column("comm.vax" = parms[["comm.vax"]]) %>%
                    add_column("sigma" =parms[["sigma1"]]) %>%
                    bind_rows(res_master) -> res_master
                  
                  
                }
            }
          }
        }
  }
}


colnames(res_master) <- c("time", "S.rNC", "E.rNC", "A.rNC", "I.rNC", "R.rNC",
                          "I.rC", "R.rC",
                          "S.hcwNC", "E.hcwNC", "A.hcwNC", "I.hcwNC", "R.hcwNC",
                          "S.hcwC", "E.hcwC", "A.hcwC", "I.hcwC", "R.hcwC", "I.hcwH",
                          "inc_r", "inc_hcw", "cum_inc_r", "cum_inc_hcw", "cum_inc_vax1", "cum_inc_vax2",
                          "cum_inc_community", "mortality","discharges","total",
                          "vax_count1", "vax_count2", 
                          "asymptomatic_inc_r_0","asymptomatic_inc_r_1","asymptomatic_inc_r_2",  
                          "symptomatic_inc_r_0","symptomatic_inc_r_1","symptomatic_inc_r_2",
                          "asymptomatic_inc_hcw_0","asymptomatic_inc_hcw_1","asymptomatic_inc_hcw_2", 
                          "symptomatic_inc_hcw_0","symptomatic_inc_hcw_1","symptomatic_inc_hcw_2", 
                          "daily_visitors","res_coverage1","res_coverage2", "t_last_case",
                          "Sim","beta", "Testing","VE_s2_r1","vax_coverage_hcw","booster","VE2","VEstaff", "I.C","visitors","comm.vax", "sigma")


write.csv(res_master, paste0(j,"_tradeoff.csv"))



