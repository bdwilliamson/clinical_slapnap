# utility functions for simulation 1: assess relative efficiency of using extra Env sequences

# compute the covariance based on R-squared
# @param r2 the r-squared
# @param var_x the variance of x (log10 PAR score)
# @param var_y the variance of y (log10 IC-80)
# @return the value of Sigma_{12} (i.e., the covariance between x and y)
get_sigma <- function(r2 = 0.345, var_y = 1, var_x = 1) {
    return(1 - r2 / (var_y * var_x))
}

# generate continuous data
# @param n the sample size
# @param epsilon the fraction to augment with auxiliary Env sequences
# @param r2 the r-squared
# @param var_y the variance of y (log10 IC-80)
# @param mu the mean of y
# @return the mean of Y using (X, Y) and (W, X, Y)
est_continuous <- function(n = 100, epsilon = 0.2, r2 = 0.345, var_y = 1, mu = 0) {
    # generate x, y
    cov <- get_sigma(r2 = r2, var_y = var_y, var_x = var_y)
    Sigma <- matrix(c(var_y, cov, cov, var_y), nrow = 2, byrow = TRUE)
    xy <- MASS::mvrnorm(n = n, mu = rep(mu, 2), Sigma = Sigma)
    w <- MASS::mvrnorm(n = n * (1 + epsilon), mu = rep(mu, 2), Sigma = Sigma)[, 2]
    mn_1 <- mean(c(xy[, 1], xy[, 2]))
    mn_2 <- mean(c(xy[, 1], xy[, 2], w))
    return(list(xy = mn_1, wxy = mn_2))
}

# compute the mean of X in the Y = 1 group based on AUC (Y = IC-80 < 1; X = log10 PAR score)
# @param auc the AUC
# @param mu0 the mean in the Y = 0 group
# @param sigma0 the variance in the Y = 0 group
# @param sigma1 the variance in the Y = 1 group
# @return the mean in the Y = 1 group
get_mu <- function(auc = 0.744, mu0 = -0.32, sigma0 = 0.005, sigma1 = 0.005) {
    return(sqrt(sigma0 + sigma1) * (-1) * qnorm(1 - auc) + mu0)
}

# generate a binary dataset
# @param n the sample size
# @param epsilon the fraction to augment with auxiliary Env sequences
# @param auc the AUC
# @param mu0 the mean in the Y = 0 group
# @param sigma0 the variance in the Y = 0 group
# @param sigma1 the variance in the Y = 1 group
# @return estimated prob of Y = 1 using (X, Y) and (W, X, Y)
est_binary <- function(n = 100, epsilon = 0.2, auc = 0.744, mu0 = -0.32, sigma0 = 0.005, sigma1 = 0.005, p_y = 0.5) {
    mu1 <- get_mu(auc = auc, mu0 = mu0, sigma0 = sigma0, sigma1 = sigma1)
    n2 <- n * (1 + epsilon)
    # generate Y
    y <- rbinom(n = n, size = 1, prob = p_y)
    y2 <- rbinom(n = n2, size = 1, prob = p_y)
    # generate X, W | Y
    x <- vector("numeric", length = n)
    w <- vector("numeric", length = n * (1 + epsilon))
    x[y == 0] <- rnorm(n = n - sum(y), mean = mu0, sd = sqrt(sigma0))
    x[y == 1] <- rnorm(n = sum(y), mean = mu1, sd = sqrt(sigma1))
    w[y2 == 0] <- rnorm(n = n2 - sum(y2), mean = mu0, sd = sqrt(sigma0))
    w[y2 == 1] <- rnorm(n = sum(y2), mean = mu1, sd = sqrt(sigma1))
    mn_1 <- mean(c(10 ^ x > 0.5, y))
    mn_2 <- mean(c(10^x > 0.5, 10^w > 0.5, y))
    return(list(xy = mn_1, wxy = mn_2))
}
