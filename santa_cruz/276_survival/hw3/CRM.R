rm(list=ls(all=TRUE))
set.seed(5554)
source("fn_CRM.R")

#########################################################################################################
###  I implement CRM with power model  p(d)= (p^0_d)^a where p^0_d: skeleton (initial value of p(d)) and a: parameter
###  I put p(a)=Gamma(a_gam, b_gam)
#########################################################################################################

nsims = 100
out = rep(list(matrix(0, nsims, 7)), 3)

###  True probability --- p(d)
# TR_p <- c(0.05, 0.15, 0.30, 0.45, 0.60)  ### scenario 1
# TR_p <- c(0.05, 0.10, 0.20, 0.30, 0.50)  ### scenario 2
# TR_p <- c(0.15, 0.30, 0.45, 0.60, 0.85)  ### scenario 3

TR_vec = list(c(0.05, 0.15, 0.30, 0.45, 0.60),
    c(0.05, 0.10, 0.20, 0.30, 0.50),
    c(0.15, 0.30, 0.45, 0.60, 0.85))

for (j in 1:length(TR_vec)){
    for (i in 1:nsims){
        TR_p = TR_vec[[j]]

        ### design parameters
        design_para <- NULL
        design_para$m <- 5  ### # of dose levels
        design_para$chrt_size <- 3  ## how many subjects for one group (cohort)
        design_para$N_p <- 20    ### max # of subjects in a trial
        design_para$N_chrt <- design_para$N_p/design_para$chrt_size   ### how many cohorts
        design_para$TTL <- 0.3  ## target toxicity level
        
        ## hyperparameters
        hyper <- NULL
        hyper$d_0 <- c(0.05, 0.15, 0.30, 0.45, 0.60)  ## skeleton
        hyper$a_gam <- 1
        hyper$b_gam <- 1
        
        ## MCMC parameter
        n_burn <- 3000
        n_sam <- 3000
        
        #####################################################
        cur_dat <- NULL
        #####################################################
        ###  FIRST cohort
        cur_dat$a <- 1  ### set a=1 for chrt 1
        cur_dat$d_star <- which.min(abs(hyper$d_0 - design_para$TTL))  ## current MTD
        cur_dat$dose <- rep(cur_dat$d_star, design_para$chrt_size)  ## record dose level for chrt 1 (current MTD)
        cur_dat$y <- fn.sam.dat(TR_p[cur_dat$d_star], design_para$chrt_size)  ## simulate y using the true prob  --- 1: DLT, 0: NO DLT
        
        
        #####################################################
        ###  SECOND cohort  --- the LAST cohort
        
        for(i_chrt in 2:design_para$N_chrt)
        {
            ## use the current data and find the posterior dist of a
            ## then find p(d) based on the updated a and find the current MTD
            tmp <- fn.post.a(n_sam, n_burn, cur_dat$a, cur_dat$y, cur_dat$dose, hyper$d_0, hyper$a_gam, hyper$b_gam, design_para$TTL)
            cur_dat$a <- tmp$a  ## we use the post mean as the initial value of a for the next MCMC running
            cur_dat$d_star <- tmp$d_star
            
            ### assign the next cohort to the current MTD and simulate y for the next chrt
            cur_dat$dose <- c(cur_dat$dose, rep(cur_dat$d_star, design_para$chrt_size))  ## record dose level for chrt 1 (current MTD)
            cur_dat$y <- c(cur_dat$y, fn.sam.dat(TR_p[cur_dat$d_star], design_para$chrt_size))  ## simulate y using the true prob  --- 1: DLT, 0: NO DLT
        }
        
        
        ## use ALL data and find the posterior dist of a
        ## then find p(d) based on the updated a and find the *FINAL* MTD
        tmp <- fn.post.a(n_sam, n_burn, cur_dat$a, cur_dat$y, cur_dat$dose, hyper$d_0, hyper$a_gam, hyper$b_gam, design_para$TTL)
        cur_dat$a <- tmp$a  ## we use the post mean as the initial value of a for the next MCMC running
        cur_dat$d_star <- tmp$d_star
        
        out[[j]][i, unique(cur_dat$dose)] = as.numeric(table(cur_dat$dose))
        out[[j]][i, 6] = mean(cur_dat$y)
        out[[j]][i, 7] = cur_dat$d_star
        }
    }

# percentage of patients treated at the five dose levels
cri1 = matrix(0, 3, 5)
cri1[1,] = round(apply(out[[1]][,1:5], 2, sum) / (20 * nsims), 2)
cri1[2,] = round(apply(out[[2]][,1:5], 2, sum) / (20 * nsims), 2)
cri1[3,] = round(apply(out[[3]][,1:5], 2, sum) / (20 * nsims), 2)

# percentage of trials recommending the true mtd as the mtd
cri2 = matrix(0, 3, 5)
cri2[1, sort(unique(out[[1]][,7]))] = table(out[[1]][,7]) / nsims
cri2[2, sort(unique(out[[2]][,7]))] = table(out[[2]][,7]) / nsims
cri2[3, sort(unique(out[[3]][,7]))] = table(out[[3]][,7]) / nsims

# overall percent of dlt
cri3 = double(3)
cri3[1] = round(mean(out[[1]][,6]), 2)
cri3[2] = round(mean(out[[2]][,6]), 2)
cri3[3] = round(mean(out[[3]][,6]), 2)

