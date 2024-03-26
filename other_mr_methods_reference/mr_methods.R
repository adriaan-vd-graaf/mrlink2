## This file contains reference functions for PCA-MR and MR-IVW-LD
## These functions accept separate dataframes for an exposure and an outcome, ordered for the same snps, and an LD matrix,
## also ordered with the same snps. They will output causal estimates in terms of alpha estimate, se(alpha) and a p value
## The mr_ivw_ld function also requires you to provide a list of instruments (that indices of the dataframes and ld matrix)


mr_pca <- function(exp_df, out_df, ld_mat){
  betaXG <- t(exp_df$beta)
  betaYG <- t(out_df$beta)
  sebetaXG <- t(as.vector(exp_df$ses))
  sebetaYG <- t(as.vector(out_df$ses))
  rho <- ld_mat

  start_time <- Sys.time()
  ### PCA MR uses all SNPs, so doing it first before selecting for instruments.

  Phi <- as.vector(betaXG/sebetaYG)%o%as.vector(betaXG/sebetaYG)*rho
  pca_obj <- prcomp(Phi, scale=FALSE) ## this takes 20 secs per run.. very long...


  K <- which(cumsum(pca_obj$sdev^2/sum((pca_obj$sdev^2)))>0.99)[1]
        # K is number of principal components to include in analysis
        # this code includes principal components to explain 99% of variance in the risk factor
  betaXG0 <- as.numeric(betaXG%*%pca_obj$rotation[,1:K])
  betaYG0 <- as.numeric(betaYG%*%pca_obj$rotation[,1:K])
  Omega <- as.vector(sebetaYG)%o%as.vector(sebetaYG) * rho
  pcOmega <- t(pca_obj$rotation[,1:K])%*%Omega%*%pca_obj$rotation[,1:K]

  inv_pcOmega <- solve(pcOmega)
  intermediate_bx_inv_omega_bx <- solve(t(betaXG0)%*%inv_pcOmega%*%betaXG0)
  beta_IVWcorrel.pc <- intermediate_bx_inv_omega_bx*t(betaXG0)%*%inv_pcOmega%*%betaYG0
  se_IVWcorrel.fixed.pc <- sqrt(intermediate_bx_inv_omega_bx)

  return(c(beta_IVWcorrel.pc, # pca_mr estimate
           se_IVWcorrel.fixed.pc,  ## standard error of the PCA MR estimate
           2*pnorm(-1*abs(beta_IVWcorrel.pc/se_IVWcorrel.fixed.pc)) ## p value of the PCA MR estimate
      ))
}

mr_ivw_ld <- function(exp_df, out_df, instruments, ld_mat){
    instruments <- this_sim[(this_sim$exp_or_outcome == 'exposure')]$selected_as_instrument
    rho <- ld_mat[instruments,instruments]
    betaXG <- as.vector(exp_df$beta)
    betaYG <- as.vector(out_df$beta)
    sebetaXG <- as.vector(exp_df$ses)
    sebetaYG <- as.vector(out_df$ses)

    Omega <- sebetaYG%o%sebetaYG*rho
    inv_omega <- solve(Omega)
    bx_invomega_bx <- solve(t(betaXG)%*%inv_omega%*%betaXG)

    beta_IVWcorrel <- bx_invomega_bx*t(betaXG)%*%inv_omega%*%betaYG
    se_IVWcorrel.fixed<- sqrt(bx_invomega_bx)
    resid <- betaYG-beta_IVWcorrel*betaXG


     return(c(beta_IVWcorrel, # mr_ivw_ld estimate
           se_IVWcorrel.fixed,  ## standard error of themr_ivw_ld estimate
           2*pnorm(-1*abs(beta_IVWcorrel.pc/se_IVWcorrel.fixed)) ## p value of the PCA MR estimate
      ))
}

