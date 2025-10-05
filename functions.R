
#Core functions for bandwidth selection project


# =======================================================
# 1. Regression function m(x)
# =======================================================
m_true <- function(x) {
  sin(1 / (x/3 + 0.1))
}

# =======================================================
# 2. Data generation
# Generate covariates X ~ Beta(alpha, beta)
# Generate Y = m(X) + epsilon, epsilon ~ N(0, sigma^2)
# =======================================================
simulate_data <- function(n, alpha, beta, sigma2) {
  set.seed(123) #for reproducibility
  X <- rbeta(n, alpha, beta)
  m <- m_true(X)
  Y <- m + rnorm(n, mean = 0, sd = sqrt(sigma2))
  data.frame(X = X, Y = Y)
}

# =======================================================
# 3. Blockwise polynomial fitting (degree 4)
# Divide support into N equal-width blocks
# Fit polynomial of degree 4 within each block:
#   Y = b0 + b1*X + b2*X^2 + b3*X^3 + b4*X^4 + error
# Compute blockwise residuals and second derivatives
# Return estimates of theta22(N), sigma2(N), and RSS(N)
# =======================================================
blockwise_fit <- function(X, Y, N) {
  n <- length(Y)
  block_id <- cut(X, breaks = seq(0, 1, length.out = N+1), labels = FALSE, 
                  include.lowest = TRUE)
  
  # Storage
  second_deriv_sq <- numeric(n)
  residuals_sq <- numeric(n)
  RSS_total <- 0
  
  for (j in 1:N) {
    idx <- which(block_id == j)
    if (length(idx) < 6) next  # need at least 6 obs for non null residual
    
    # Fit polynomial of degree 4
    df <- data.frame(x = X[idx], y = Y[idx])
    fit <- lm(y ~ x + I(x^2) + I(x^3) + I(x^4), data = df)
    
    # Second derivative: m''(x) = 2 b2 + 6 b3 x + 12 b4 x^2
    b <- coef(fit)
    sec_deriv <- 2*b["I(x^2)"] + 6*b["I(x^3)"]*X[idx] + 12*b["I(x^4)"]*X[idx]^2
    second_deriv_sq[idx] <- sec_deriv^2
    
    
    # Residuals
    res <- residuals(fit)
    residuals_sq[idx] <- res^2
    RSS_total <- RSS_total + sum(res^2)
  }
  
  # Estimates
  theta22_hat <- mean(second_deriv_sq, na.rm = TRUE)
  df_resid <- n - 5*N
  sigma2_hat <- if (df_resid > 0) sum(residuals_sq, na.rm = TRUE)/df_resid else NA
  
  list(theta22 = theta22_hat, sigma2 = sigma2_hat, RSS = RSS_total, df = df_resid)
}

# =======================================================
# 4. Bandwidth formula
# =======================================================
h_amise <- function(n, sigma2_hat, theta22_hat, support_length = 1) {
  # AMISE optimal bandwidth:
  # h_AMISE = n^(-1/5) * ( (35 * sigma^2 * |supp(X)|) / theta22 )^(1/5)
  if (is.na(sigma2_hat) || is.na(theta22_hat) || theta22_hat <= 0) return(NA)
  n^(-1/5) * ((35 * sigma2_hat * support_length) / theta22_hat)^(1/5)
}

# =======================================================
# 5. Cp(N) criterion and optimal N
# =======================================================
cp_table <- function(X, Y) {
  # Compute Cp(N) for N = 1,...,Nmax
  # Nmax = max(min(floor(n/20),5),1)
  n <- length(Y)
  Nmax <- max(min(floor(n/20), 5), 1)
  
  results <- data.frame(N = 1:Nmax, Cp = NA, RSS = NA, df = NA)
  
  # First compute RSS(N) for each N
  for (N in 1:Nmax) {
    est <- blockwise_fit(X, Y, N)
    results$RSS[N] <- est$RSS
    results$df[N]  <- est$df
  }
  
  # Scaling denominator: RSS(Nmax)/(n - 5 Nmax)
  denom <- results$RSS[results$N == Nmax] / (n - 5*Nmax)

  if (is.na(denom) || denom <= 0) {
  warning("Invalid denominator in Cp calculation. Returning NA for Cp values.")
  results$Cp <- NA
  return(results)
  }

  # Cp(N) = RSS(N)/denom - (n - 10N)
  results$Cp <- results$RSS/denom - (n - 10*results$N)
  
  results
}

optimal_N <- function(X, Y) {
  # Select N_opt = argmin Cp(N)
  tab <- cp_table(X, Y)
  best_row <- which.min(tab$Cp)
  list(N_opt = tab$N[best_row], Cp_opt = tab$Cp[best_row], table = tab)
}
