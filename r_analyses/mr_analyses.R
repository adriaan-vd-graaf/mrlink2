require('data.table')
# require('TwoSampleMR')

mr_ivw <- function(b_exp, b_out, se_exp, se_out, parameters=default_parameters())
{
	if(sum(!is.na(b_exp) & !is.na(b_out) & !is.na(se_exp) & !is.na(se_out)) < 2)
	return(list(b=NA, se=NA, pval=NA, nsnp=NA))

	ivw.res <- summary(lm(b_out ~ -1 + b_exp, weights = 1/se_out^2))
	b <- ivw.res$coef["b_exp","Estimate"]
	se <- ivw.res$coef["b_exp","Std. Error"]/min(1,ivw.res$sigma) #sigma is the residual standard error
	pval <- 2 * pnorm(abs(b/se), lower.tail=FALSE)
	Q_df <- length(b_exp) - 1
	Q <- ivw.res$sigma^2 * Q_df
	Q_pval <- pchisq(Q, Q_df, lower.tail=FALSE)
	# from formula phi =  Q/DF rearranged to to Q = phi*DF, where phi is sigma^2
	# Q.ivw<-sum((1/(se_out/b_exp)^2)*(b_out/b_exp-ivw.reg.beta)^2)
	return(list(b = b, se = se, pval = pval, nsnp=length(b_exp), Q = Q, Q_df = Q_df, Q_pval = Q_pval))
}


mr_wald_ratio <- function(b_exp, b_out, se_exp, se_out, parameters) {
  if (length(b_exp) > 1) {
    return(list(b = NA, se = NA, pval = NA, nsnp = NA))
  }
  b <- b_out / b_exp
  se <- se_out / abs(b_exp)
  pval <- pnorm(abs(b) / se, lower.tail = FALSE) * 2
  return(list(b = b, se = se, pval = pval, nsnp = 1))
}


mr_egger_regression <- function(b_exp, b_out, se_exp, se_out, parameters) {
  stopifnot(length(b_exp) == length(b_out))
  stopifnot(length(se_exp) == length(se_out))
  stopifnot(length(b_exp) == length(se_out))
  nulllist <- list(b = NA, se = NA, pval = NA, nsnp = NA, b_i = NA, se_i = NA, pval_i = NA, Q = NA, Q_df = NA, Q_pval = NA, mod = NA, smod = NA, dat = NA)
  if (sum(!is.na(b_exp) & !is.na(b_out) & !is.na(se_exp) & !is.na(se_out)) < 3) {
    return(nulllist)
  }
  sign0 <- function(x) {
    x[x == 0] <- 1
    return(sign(x))
  }
  to_flip <- sign0(b_exp) == -1
  b_out = b_out * sign0(b_exp)
  b_exp = abs(b_exp)
  dat <- data.frame(b_out = b_out, b_exp = b_exp, se_exp = se_exp, se_out = se_out, flipped = to_flip)
  mod <- lm(`~`(b_out, b_exp), weights = 1 / se_out ^ 2)
  smod <- summary(mod)
  if (nrow(coefficients(smod)) > 1) {
    b <- coefficients(smod)[2, 1]
    se <- coefficients(smod)[2, 2] / min(1, smod$sigma)
    pval <- 2 * pt(abs(b / se), length(b_exp) - 2, lower.tail = FALSE)
    b_i <- coefficients(smod)[1, 1]
    se_i <- coefficients(smod)[1, 2] / min(1, smod$sigma)
    pval_i <- 2 * pt(abs(b_i / se_i), length(b_exp) - 2, lower.tail = FALSE)
    Q <- smod$sigma ^ 2 * (length(b_exp) - 2)
    Q_df <- length(b_exp) - 2
    Q_pval <- pchisq(Q, Q_df, lower.tail = FALSE)
  } else {
    warning("Collinearities in MR Egger, try LD pruning the exposure variables.")
    return(nulllist)
  }
  return(list(b = b, se = se, pval = pval, nsnp = length(b_exp), b_i = b_i, se_i = se_i, pval_i = pval_i, Q = Q, Q_df = Q_df, Q_pval = Q_pval, mod = smod, dat = dat))
}



args <- commandArgs(trailingOnly=TRUE)

input_file <- args[1]
output_file <- args[2]
ld_mat_file <- args[3]

full_df <- fread(input_file)
ld_mat <- readBin(ld_mat_file, what='double', n=1e9)
dim(ld_mat) <- c(sqrt(length(ld_mat)), sqrt(length(ld_mat))) ## gotta love R for allowing this...
snp_names <- paste0('snp', 1:nrow(ld_mat))
colnames(ld_mat) <- snp_names
rownames(ld_mat) <- snp_names

simulations <- unique(full_df$simulation)

results <- matrix(NA, nrow=length(simulations), ncol=11)


results_df <-as.data.frame(matrix(nrow=length(simulations),ncol=14))
colnames(results_df) <- c('simulation_idx', 'n_instruments',
                          'beta_ivw', 'se_ivw', 'p_ivw',
                          'beta_egger', 'se_egger', 'p_egger',
                          'beta_ivw_r', 'se_ivw_r', 'p_ivw_r',
                          'beta_pca', 'se_pca', 'p_pca'
                          )

for(simulation_idx in simulations){
  this_sim <- full_df[full_df$simulation == simulation_idx,]

  exp_df <- this_sim[(this_sim$exp_or_outcome == 'exposure'),]
  out_df  <- this_sim[(this_sim$exp_or_outcome == 'outcome'),]

  results_df[simulation_idx+1, 'simulation_idx'] <- simulation_idx

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

  results_df[simulation_idx+1,'beta_pca'] <- beta_IVWcorrel.pc
  results_df[simulation_idx+1,'se_pca'] <- se_IVWcorrel.fixed.pc
  results_df[simulation_idx+1,'p_pca'] <- 2*pnorm(-1*abs(beta_IVWcorrel.pc/se_IVWcorrel.fixed.pc))

  ### the other methods
  exp_df <- this_sim[(this_sim$exp_or_outcome == 'exposure') & (this_sim$selected_as_instrument),]
  out_df  <- this_sim[(this_sim$exp_or_outcome == 'outcome') & (this_sim$selected_as_instrument),]

  BetaXG <- exp_df$beta
  BetaYG <- out_df$beta
  sebetaXG <- as.vector(exp_df$ses)
  sebetaYG <- as.vector(out_df$ses)
  rho <- ld_mat


  results_df[simulation_idx+1, 'n_instruments'] <- dim(exp_df)[1]
  if(dim(exp_df)[1] == 1){
    ## IVW
    ivw_results <- mr_wald_ratio(exp_df$beta, out_df$beta, se_exp=exp_df$ses, se_out=out_df$ses)
    results_df[simulation_idx+1,'beta_ivw'] <- ivw_results$b
    results_df[simulation_idx+1,'se_ivw'] <- ivw_results$se
    results_df[simulation_idx+1,'p_ivw'] <- ivw_results$pval
  }

  if(dim(exp_df)[1] >1){
    ## IVW
    ivw_results <- mr_ivw(exp_df$beta, out_df$beta, se_exp=exp_df$ses, se_out=out_df$ses)

    results_df[simulation_idx+1,'beta_ivw'] <- ivw_results$b
    results_df[simulation_idx+1,'se_ivw'] <- ivw_results$se
    results_df[simulation_idx+1,'p_ivw'] <- ivw_results$pval
  }

  if(dim(exp_df)[1] >= 3){
  ## Egger
    egger_results <- mr_egger_regression(exp_df$beta, out_df$beta, se_exp=exp_df$ses, se_out=out_df$ses)
    results_df[simulation_idx+1,'beta_egger'] <- egger_results$b
    results_df[simulation_idx+1,'se_egger'] <- egger_results$se
    results_df[simulation_idx+1,'p_egger'] <- egger_results$pval
  }

  if(dim(exp_df)[1] >= 1){
    ## IVW correlated
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
    se_IVWcorrel.random <- sqrt(bx_invomega_bx)*max(sqrt(t(resid)%*%inv_omega%*%resid/(length(betaXG)-1)),1)

    results_df[simulation_idx+1,'beta_ivw_r'] <- beta_IVWcorrel
    results_df[simulation_idx+1,'se_ivw_r'] <- se_IVWcorrel.fixed
    results_df[simulation_idx+1,'p_ivw_r'] <- 2*pnorm(-1*abs(beta_IVWcorrel/se_IVWcorrel.fixed))
  }


}

data.table::fwrite(results_df, file=output_file, sep='\t')