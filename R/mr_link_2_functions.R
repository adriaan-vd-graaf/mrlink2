mr_link2_loglik_reference_v2 <- function(th, lam, c_x, c_y, n_x, n_y) {
  # Convert n_x and n_y to float
  n_x <- as.numeric(n_x)
  n_y <- as.numeric(n_y)
  a <- th[1]
  tX <- abs(th[2])
  tY <- abs(th[3])

  Dyy <- 1 / (n_y * lam + tY)

  if (a != 0) {
    Dxx <- 1 / (exp(log(a^2 * n_y + n_x) + log(lam)) + tX -
                exp(log(a^2 * n_y^2 * lam^2) - log(n_y * lam + tY)))
    Dxy <- -Dxx * a * exp(log(n_y * lam) - log(n_y * lam + tY))
    Dyy <- Dyy + exp(log(Dxx * (a^2 * n_y^2 * lam^2)) - (2 * log(n_y * lam + tY)))
    asq_ny_sq_lam_sq_div_ny_lam_ty <- exp(log(a^2 * n_y^2 * lam^2) - log(n_y * lam + tY))
  } else {
    Dxx <- 1 / (exp(log(n_x) + log(lam)) + tX)
    Dxy <- -Dxx * a * exp(log(n_y * lam) - log(n_y * lam + tY))
    asq_ny_sq_lam_sq_div_ny_lam_ty <- 0 * lam
  }

  dX <- n_x * c_x + a * n_y * c_y
  dY <- n_y * c_y
  m <- length(c_x)

  loglik <- -m * log(2 * pi) -
            (1 / 2) * sum(log((a^2 * n_y + n_x) * lam + tX - asq_ny_sq_lam_sq_div_ny_lam_ty)) -
            (1 / 2) * sum(log(n_y * lam + tY)) +
            (1 / 2) * (sum(dX^2 * Dxx) + 2 * sum(dX * dY * Dxy) + sum(dY^2 * Dyy)) -
            (n_x / 2) * sum((c_x^2) / lam) -
            (n_y / 2) * sum((c_y^2) / lam) +
            (m / 2) * (log(n_x) + log(n_y)) - sum(log(lam)) + (m / 2) * (log(tX) + log(tY))

  return(-loglik)
}

mr_link2_loglik_alpha_h0 <- function(th, lam, cX, cY, nX, nY) {
  return(mr_link2_loglik_reference_v2(c(0, th[1], th[2]), lam, cX, cY, nX, nY))
}

mr_link2_loglik_sigma_y_h0 <- function(th, lam, c_x, c_y, n_x, n_y) {
  n_x <- as.numeric(n_x)
  n_y <- as.numeric(n_y)
  a <- th[1]
  tX <- abs(th[2])

  Dyy <- rep(0, length(lam))

  if (a != 0) {
    Dxx <- 1 / (exp(log(a^2 * n_y + n_x) + log(lam)) + tX)
    Dxy <- rep(0, length(lam))
    asq_ny_sq_lam_sq_div_ny_lam_ty <- rep(0, length(lam))
  } else {
    Dxx <- 1 / (exp(log(n_x) + log(lam)) + tX)
    Dxy <- rep(0, length(lam))
    asq_ny_sq_lam_sq_div_ny_lam_ty <- rep(0, length(lam))
  }

  dX <- n_x * c_x + a * n_y * c_y
  dY <- n_y * c_y
  m <- length(c_x)

  loglik <- -m * log(2 * pi) -
            (1 / 2) * sum(log((a^2 * n_y + n_x) * lam + tX - asq_ny_sq_lam_sq_div_ny_lam_ty)) +
            (1 / 2) * (sum(dX^2 * Dxx) + 2 * sum(dX * dY * Dxy) + sum(dY^2 * Dyy)) -
            (n_x / 2) * sum((c_x^2) / lam) -
            (n_y / 2) * sum((c_y^2) / lam) +
            (m / 2) * (log(n_x) + log(n_y)) - sum(log(lam)) + (m / 2) * log(tX)

  return(-loglik)
}



mr_link2 <- function(selected_eigenvalues, selected_eigenvectors,
                     exposure_betas, outcome_betas,
                     n_exp, n_out, sigma_exp_guess, sigma_out_guess) {

  start_time <- Sys.time()

  # Define optimization options
  method <- "Nelder-Mead"
  control <- list(maxit = 100000)

  # Calculate c_x and c_y
  c_x <- t(selected_eigenvectors) %*% exposure_betas
  c_y <- t(selected_eigenvectors) %*% outcome_betas

  max_sigma <- sqrt(.Machine$double.xmax)

  # Alpha h0 estimation
  alpha_h0_guesses <- list(c(sigma_exp_guess, sigma_out_guess),
                           c(max_sigma, max_sigma),
                           c(1, max_sigma),
                           c(1e3, 1e3))

  alpha_h0_results <- optim(par = alpha_h0_guesses[[1]],
                            fn = mr_link2_loglik_alpha_h0,
                            gr = NULL,
                            selected_eigenvalues, c_x, c_y, n_exp, n_out,
                            method = method, control = control)

  # Loop over other guesses
  for (alpha_h0_guess in alpha_h0_guesses[-1]) {
    if (alpha_h0_results$convergence == 0) break

    new_alpha_h0_results <- optim(par = alpha_h0_guess,
                                  fn = mr_link2_loglik_alpha_h0,
                                  gr = NULL,
                                  selected_eigenvalues, c_x, c_y, n_exp, n_out,
                                  method = method, control = control)

    if (alpha_h0_results$value >= new_alpha_h0_results$value) {
      alpha_h0_results <- new_alpha_h0_results
    }
  }

  # Sigma_y estimation
  sigma_y_guesses <- list(c(0.0, sigma_exp_guess),
                          c(1.0, sigma_exp_guess),
                          c(0.0, alpha_h0_results$par[1]),
                          c(1.0, alpha_h0_results$par[1]),
                          c(0.0, max_sigma),
                          c(1e-10, max_sigma))

  sigma_y_h0_results <- optim(par = sigma_y_guesses[[1]],
                              fn = mr_link2_loglik_sigma_y_h0,
                              gr = NULL,
                              selected_eigenvalues, c_x, c_y, n_exp, n_out,
                              method = method, control = control)

  for (sigma_y_guess in sigma_y_guesses[-1]) {
    if (sigma_y_h0_results$convergence == 0) break

    new_sigma_y_h0_results <- optim(par = sigma_y_guess,
                                    fn = mr_link2_loglik_sigma_y_h0,
                                    gr = NULL,
                                    selected_eigenvalues, c_x, c_y, n_exp, n_out,
                                    method = method, control = control)

    if (new_sigma_y_h0_results$value < sigma_y_h0_results$value) {
      sigma_y_h0_results <- new_sigma_y_h0_results
    }
  }

  # Ha estimation
  ha_guesses <- list(c(0.0, alpha_h0_results$par[1], alpha_h0_results$par[2]),
                     c(sigma_y_h0_results$par[1], sigma_y_h0_results$par[2], sqrt(.Machine$double.xmax)),
                     c(1.0, alpha_h0_results$par[1], alpha_h0_results$par[2]),
                     c(1e-10, max_sigma, max_sigma))

  ha_results <- optim(par = ha_guesses[[1]],
                      fn = mr_link2_loglik_reference_v2,
                      gr = NULL,
                      selected_eigenvalues, c_x, c_y, n_exp, n_out,
                      method = method, control = control)

  for (ha_guess in ha_guesses[-1]) {
    if (ha_results$convergence == 0) break

    new_ha_result <- optim(par = ha_guess,
                           fn = mr_link2_loglik_reference_v2,
                           gr = NULL,
                           selected_eigenvalues, c_x, c_y, n_exp, n_out,
                           method = method, control = control)

    if (new_ha_result$value < ha_results$value) {
      ha_results <- new_ha_result
    }
  }

  # Likelihood Ratio Test and Estimation
  alpha <- ha_results$par[1]
  alpha_chi_sq <- 2 * (alpha_h0_results$value - ha_results$value)
  alpha_p_val <- pchisq(alpha_chi_sq, df = 1, lower.tail = FALSE)
  z_alpha <- ifelse(alpha_chi_sq <= 0, 0.0, sign(alpha) * sqrt(alpha_chi_sq))
  se_alpha <- ifelse(z_alpha != 0, alpha / z_alpha, NA)

  sigma_y <- 1 / abs(ha_results$par[3])
  sigma_y_chi_sq <- 2 * (sigma_y_h0_results$value - ha_results$value)
  sigma_y_p_val <- pchisq(sigma_y_chi_sq, df = 1, lower.tail = FALSE)
  z_sigma_y <- ifelse(sigma_y_chi_sq <= 0, 0.0, sqrt(sigma_y_chi_sq))
  se_sigma_y <- ifelse(z_sigma_y != 0, sigma_y / z_sigma_y, NA)

  # Return results as a list (equivalent to Python's dictionary)
  list(
    alpha = alpha,
    `se(alpha)` = se_alpha,
    `p(alpha)` = alpha_p_val,
    sigma_y = sigma_y,
    `se(sigma_y)` = se_sigma_y,
    `p(sigma_y)` = sigma_y_p_val,
    sigma_x = 1 / abs(ha_results$par[2]),
    alpha_h0_sigma_x = 1 / abs(alpha_h0_results$par[1]),
    alpha_h0_sigma_y = 1 / abs(alpha_h0_results$par[2]),
    alpha_h0_loglik = alpha_h0_results$value,
    sigma_y_h0_alpha = sigma_y_h0_results$par[1],
    sigma_y_h0_sigma_x = 1 / abs(sigma_y_h0_results$par[2]),
    sigma_y_h0_loglik = sigma_y_h0_results$value,
    ha_loglik = ha_results$value,
    optim_alpha_h0_success = alpha_h0_results$convergence == 0,
    optim_alpha_h0_nit = alpha_h0_results$counts,
    optim_sigma_y_h0_success = sigma_y_h0_results$convergence == 0,
    optim_sigma_y_h0_nit = sigma_y_h0_results$counts,
    optim_ha_success = ha_results$convergence == 0,
    optim_ha_nit = ha_results$counts,
    function_time = as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  )
}