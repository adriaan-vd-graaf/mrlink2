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