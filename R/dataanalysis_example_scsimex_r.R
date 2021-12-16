##############################################
# Functions for SC-SIMEX                     #
##############################################

get_W <- function(VectorID_to_coar, dist.subj.efs.truecoord.list, dist.subj.efs.zipcentcoord.list){
  new.list <- dist.subj.efs.truecoord.list
  new.list[VectorID_to_coar] <- dist.subj.efs.zipcentcoord.list[VectorID_to_coar]
  w.mat <- matrix(NA, nrow(my.subj.data),Q)
  for(j in 1:Q){
    a <- unlist(new.list[which(ServiceType==j)])
    w.mat[,j] <- tabulate(a, nbins = nrow(my.subj.data))
    rm(a)
  }
  return(w.mat)
}

#Function that samples business id to be coarsened in addition to the observed data 
sample_id_add_coarsening<- function(add.coarsened_num, efs, Q, sample_ID_purchased, J){
  sample_ID=NA
  if (Q>1) {
    start <- 1
    end <- efs[1]
    
    for (j in 1:Q) {
      sample_ID <- c(sample_ID,sample(setdiff(start:end,sample_ID_purchased),add.coarsened_num[j]))
      start <- start+efs[j]
      end <- end+efs[j+1]
    }
  } else {
    sample_ID <- sample(setdiff(1:J,sample_ID_purchased),add.coarsened_num)
  }
  return(sample_ID[-1])
}


get_beta <- function(add.coarsened_num, boot_subj_data){
  add_coar_efs_id <- NULL
  add_coar_efs_id <- sample_id_add_coarsening(add.coarsened_num, efs, Q, unlist(init.coarsen_id, use.names=FALSE), J)
  
  #replace true coords to zip centroid
  coar_efs_id <- c(add_coar_efs_id, init.coarsen_id, use.names=FALSE)
  W_b <- get_W(coar_efs_id, dist.subj.efs.truecoord.list, dist.subj.efs.zipcentcoord.list)[boot_subj_data$id,]
  
  
  #generate betas
  stratum1 <- boot_subj_data$numcent==0
  stratum2 <- boot_subj_data$numcent>0
  
  formula1 <- as.formula("boot_subj_data[stratum1,]$y ~ boot_subj_data[stratum1,]$z + W_b[stratum1,]")
  formula2 <- as.formula("boot_subj_data[stratum2,]$y ~ boot_subj_data[stratum2,]$z + W_b[stratum2,]")
  regression1 <- lm( formula1 ) 
  regression2 <- lm( formula2 ) 
  out1 <- regression1$coefficients
  out2 <- regression2$coefficients
  out <- list(stratum_without = out1,
              stratum_with = out2)
  return(out)
}


sc_simex <- function(my.subj.data, K, S, extrapolation=c("wcubic", "wquad")){
simex.with.k <- matrix(NA, nrow=K, ncol=length(PredictorName))
simex.without.k <- matrix(NA, nrow=K, ncol=length(PredictorName))

pb = txtProgressBar(style=3)
for(k in 1:K){
  #resampling subjects
  boot_id <- sample(my.subj.data$id, replace=T)
  boot_subj_data <- my.subj.data[boot_id,]
  
  avgbeta_sim_with <- NULL; avgbeta_sim_without <- NULL
  varbeta_sim_with <- NULL; varbeta_sim_without <- NULL
  for (s in 1:S)
  {
    init.coarsen_num <- ceiling(efs*lambda.naive)
    add.coarsened_num <- ceiling(efs*(lambda[s]-lambda.naive)) #number of businesses to be coarsened additionally
    
    
    #for each lambda_i, get coefficient estimates from W_b 
    for (q in 1:Q) { if( (init.coarsened_num[q]+add.coarsened_num[q])>efs[q] ){add.coarsened_num[q]=efs[q]-init.coarsened_num[q]} }
    
    betamat <- replicate(B, get_beta(add.coarsened_num, boot_subj_data)) 
    avg_without <- apply(do.call(rbind,betamat[rownames(betamat)=="stratum_without",]), 2, function(x){mean(x, na.rm=T)})
    var_without <- apply(do.call(rbind,betamat[rownames(betamat)=="stratum_without",]), 2, function(x){var(x, na.rm=T)} )
    avg_with <- apply(do.call(rbind,betamat[rownames(betamat)=="stratum_with",]), 2, function(x){mean(x, na.rm=T)})
    var_with <- apply(do.call(rbind,betamat[rownames(betamat)=="stratum_with",]), 2, function(x){var(x, na.rm=T)} )
    avgbeta_sim_without <- rbind(avgbeta_sim_without, c(lambda[s], avg_without))
    varbeta_sim_without <- rbind(varbeta_sim_without, c(lambda[s], var_without))
    avgbeta_sim_with <- rbind(avgbeta_sim_with, c(lambda[s], avg_with))
    varbeta_sim_with <- rbind(varbeta_sim_with, c(lambda[s], var_with))
  }

  
   if(extrapolation == 'wquad'){
     for(i in 1:length(PredictorName)){
       fit <- lm(avgbeta_sim_with[,(i+1)]~lambda + I(lambda^2), 
                 weights = 1/(varbeta_sim_with[,(i+1)]+1e-3), na.action = na.exclude)
       simex.with.k[k,i] <- predict(fit, data.frame(lambda = 0))
     }
     
     for(i in 1:length(PredictorName)){
       fit <- lm(avgbeta_sim_without[,(i+1)]~lambda + I(lambda^2), 
                 weights = 1/(varbeta_sim_without[,(i+1)]+1e-3), na.action = na.exclude)
       simex.without.k[k,i] <- predict(fit, data.frame(lambda = 0))
     }
   }else if(extrapolation == 'wcubic'){
    for(i in 1:length(PredictorName)){
      fit <- lm(avgbeta_sim_with[,(i+1)]~lambda + I(lambda^2)+ I(lambda^3), 
                weights = 1/(varbeta_sim_with[,(i+1)]+1e-3), na.action = na.exclude)
      simex.with.k[k,i] <- predict(fit, data.frame(lambda = 0))
    }

    for(i in 1:length(PredictorName)){
      fit <- lm(avgbeta_sim_without[,(i+1)]~lambda + I(lambda^2)+ I(lambda^3), 
                weights = 1/(varbeta_sim_without[,(i+1)]+1e-3), na.action = na.exclude)
      simex.without.k[k,i] <- predict(fit, data.frame(lambda = 0))
    }
   }
  setTxtProgressBar(pb,k/K)
}
  close(pb)
  

  empvar_with <- apply(simex.with.k, 2, function(x)var(x,na.rm=T))
  empvar_without <- apply(simex.without.k, 2, function(x)var(x,na.rm=T))
  
  w_without <- (1/empvar_without) / (1/empvar_with + 1/empvar_without)
  w_with <- (1/empvar_with) / (1/empvar_with + 1/empvar_without)
  b_without <- colMeans(simex.without.k, na.rm = TRUE)
  b_with <- colMeans(simex.with.k, na.rm = TRUE)
  
  
  comb_simex <- b_without*w_without + b_with*w_with
  v_comb_simex <- w_without^2*empvar_without + w_with^2*empvar_with 
  
  alpha=0.05
  ub <- comb_simex + qnorm(1-alpha/2,lower.tail = T)*sqrt(v_comb_simex)
  lb <- comb_simex - qnorm(1-alpha/2,lower.tail = T)*sqrt(v_comb_simex)
  resultscsimex <- data.frame(comb_simex, sqrt(v_comb_simex), lb, ub)
  return(resultscsimex)
}



