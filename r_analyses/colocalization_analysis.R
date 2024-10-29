require('coloc')
require('data.table')

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


for(simulation_idx in simulations){
  this_sim <- full_df[full_df$simulation == simulation_idx,]

  exp_df <- this_sim[this_sim$exp_or_outcome == 'exposure',]
  out_df  <- this_sim[this_sim$exp_or_outcome == 'outcome',]

  exp_varbeta <- (abs(exp_df$beta)/abs(qnorm(exp_df$p_values/2)))^2
  out_varbeta <- (abs(out_df$beta)/abs(qnorm(out_df$p_values/2)))^2
  filter <- (exp_varbeta != 0) & (out_varbeta != 0)

  exp_list <- list(beta=exp_df$beta[filter], varbeta=exp_varbeta[filter],
                   MAF=exp_df$MAF[filter], pvalues=exp_df$p_values[filter], N=exp_df$N[1], type="quant", LD=ld_mat,
                   snp=snp_names[filter])
  out_list <- list(beta=out_df$beta[filter], varbeta=out_varbeta[filter],
                   MAF=out_df$MAF[filter], pvalues=out_df$p_values[filter], N=out_df$N[1], type="quant", LD=ld_mat,
                   snp=snp_names[filter])



  coloc_results <- coloc.abf(exp_list, out_list)
  susie_results <- coloc.susie(runsusie(exp_list, n=10000, maxit=100), runsusie(out_list, n=300000, maxit=100))
  susie_results <- susie_results$summary[which.max(susie_results$summary$PP.H4.abf),]

  if(!is.null(susie_results)){
    susie_names <- names(susie_results)
  }

  results[simulation_idx+1,] <- c(as.numeric(coloc_results$summary),
                                  as.numeric(susie_results)[4:8])
}


if(!is.null(susie_results)){
  cat('debug1\n')
  colnames(results) <- c(names(coloc_results$summary), paste0('susie_max_', susie_names[4:8]))
}else{
  cat('debug2\n')
  cat(results)
  colnames(results) <- c(names(coloc_results$summary), paste0('susie_max_', names(coloc_results$summary)[1:5]))
}
data.table::fwrite(results, file=output_file, sep='\t')
