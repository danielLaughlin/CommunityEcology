
traitspace <- function(data, multi.mod, env_p, PT.Sk, N = 100, avg.out = TRUE, norm = TRUE, parallel = FALSE, mc.cores = NULL){
  
  # Arguments
  # multi.mod : trait-environment regression model
  # env_p : the environmental data for prediction
  # PT.Sk: the trait distribution of species
  # avg.out: make the last MCMC integration step or not.
  # parallel: parallelization option to speed up the calculation
  # N: number of draws
  # norm: transform into relative abundances or not
  
  require(mclust)
  require(mvtnorm)
  # Setting up parallelisation if the option was selected.
  if (parallel){
    require(parallel)
    if (is.null(mc.cores)){
      mc.cores = detectCores() - 1
    }}
  #Extraction des données du modèle trait-environnement
  pred_multi <- predict(multi.mod,newdata=env_p,interval="prediction",level=0.95)
  row.names(pred_multi) <- row.names(env_p)
  cov_mat <- cov(as.matrix(multi.mod$residuals))
  
  # Construction of the simulated data table
  ## STEP 2a: Drawing samples from P(T/E)
  env_prd <- env_p[gl(nrow(env_p), N),]
  tr_mean_pred <- pred_multi[gl(nrow(env_p), N), ]
  if (all(dim(cov_mat) == c(1,1))){
    tr_mean_pred <- tr_mean_pred[,1,drop = F]
    trait_sample <- as.data.frame(apply(tr_mean_pred, 1, function(x) rmvnorm(1, x, cov_mat)))
  }else{
    trait_sample <- as.data.frame(t(apply(tr_mean_pred, 1, function(x) rmvnorm(1, x, cov_mat))))
    }
 
  ## computing(P(T/E))
  if (parallel){
    P_T_E <-do.call(rbind, mclapply(1:nrow(trait_sample), function(j) dmvnorm(trait_sample[j,],tr_mean_pred[j,], cov_mat), mc.cores = mc.cores))
  }else{
    P_T_E <-do.call(rbind, lapply(1:nrow(trait_sample), function(j) dmvnorm(trait_sample[j,],tr_mean_pred[j,], cov_mat)))
  }
  row.names(P_T_E) <- row.names(env_prd)
  
  ## Step 2b: Computing the likelihood P(T/Sk) using Mclust done earlier
  if (all(dim(cov_mat) == c(1,1))){
    if (parallel){
      P_T_S <- mclapply(PT.Sk, function(pdf){
        mclust::dens(modelName=pdf$modelName,data=trait_sample[,1],parameters=pdf$parameters)
      }, mc.cores = mc.cores)}
    else{
      P_T_S <- lapply(PT.Sk, function(pdf){
        mclust::dens(modelName=pdf$modelName,data=trait_sample[,1],parameters=pdf$parameters)
      })
    }
  }else{
    if (parallel){
      P_T_S <- mclapply(PT.Sk, function(pdf){
        mclust::dens(modelName=pdf$modelName,data=trait_sample,parameters=pdf$parameters)
      }, mc.cores = mc.cores)}
    else{
      P_T_S <- lapply(PT.Sk, function(pdf){
        mclust::dens(modelName=pdf$modelName,data=trait_sample,parameters=pdf$parameters)
      })
    }
  }

  
  P_T_S <- do.call(cbind, P_T_S)
  colnames(P_T_S) <- names(PT.Sk)
  ## Step 2c: Computing posterior P(Sk/T,E)using Bayes theorem
  # P_T_S_pr = P_T_S/ncol(P_T_S)  #No need to multiply likelihood by flat prior - numerator in Bayes thm
  # P_S_T_E <- exp(sweep(log(P_T_S_pr), 1, log(rowSums(P_T_S_pr)), "-"))
  P_S_T_E <- exp(sweep(log(P_T_S), 1, log(rowSums(P_T_S)), "-"))
  
  ## Step 2d: Posterior P(Sk/T) by integrating out T's (with log)
  P_S_E_all <- exp(sweep(log(P_S_T_E), 1, log(P_T_E), "+"))

  if (!avg.out){
    colnames(P_S_E_all) <- names(PT.Sk)
    if (is.null(row.names(env_p))) row.names(env_p) <- as.character(1:nrow(env_p))
    row.names(P_S_E_all) <- paste(row.names(env_p)[gl(nrow(env_p), N)],rep(1:N, nrow(env_p)) , sep ="_")
    return(P_S_E_all )
  }else{
    ### Step 2d.1: Original Traitspace predictions
    #Monte Carlo integration across trait samples
    P_S_E_all.sp <- split(as.data.frame(P_S_E_all), as.factor(gl(nrow(env_p), N)))
    if (parallel){
      P_S_E_traitspace_unnorm <- do.call(rbind, mclapply(P_S_E_all.sp, function(x) apply(x, 2, mean), mc.cores = mc.cores))
    }else{
      P_S_E_traitspace_unnorm <- do.call(rbind, lapply(P_S_E_all.sp, function(x) apply(x, 2, mean)))
    }
    row.names(P_S_E_traitspace_unnorm) <- row.names(env_p)
    if (norm){
    P_S_E_traitspace <- sweep(P_S_E_traitspace_unnorm, 1, rowSums(P_S_E_traitspace_unnorm) ,"/" )
    }else{
      P_S_E_traitspace <- P_S_E_traitspace_unnorm
    }
    return(P_S_E_traitspace)
  }
}


nullModel <- function(obs, result, permutations){
  require(topicmodels)
  dis.obs <-  mean(sapply(1:nrow(obs), function(j) distHellinger(t(as.matrix(obs[j,])), t(as.matrix(result[j,])))))
  
  NM.distri <- sapply(1:permutations, function(i){
    res.rand <- result[shuffle(nrow(obs)),]
    mean(sapply(1:nrow(obs), function(j) distHellinger(t(as.matrix(obs[j,])), t(as.matrix(res.rand[j,])))))
  })
  SES <- (dis.obs - mean(NM.distri))/sd(NM.distri)
  pval <- sum(dis.obs > NM.distri)/permutations
  
  list(SES = SES, pval = pval)
}

