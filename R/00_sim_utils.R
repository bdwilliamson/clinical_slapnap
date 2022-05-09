# utility functions for simulation 1: assess relative efficiency of using extra Env sequences

# for simulation 1 -------------------------------------------------------------

# compute the covariance based on R-squared
# @param r2 the r-squared
# @param var_x the variance of x (log10 PAR score)
# @param var_y the variance of y (log10 IC-80)
# @return the value of Sigma_{12} (i.e., the covariance between x and y)
get_sigma <- function(r2 = 0.345, var_y = 1, var_x = 1) {
    return((1 - r2) * (var_y * var_x))
}

# generate continuous data
# @param n the sample size
# @param epsilon the fraction to augment with auxiliary Env sequences
# @param r2 the r-squared
# @param var_y the variance of y (log10 IC-80)
# @param mu the mean of y
# @return a data.frame (W, R, RY)
gen_data_continuous <- function(n = 100, epsilon = 0.2, r2 = 0.345, var_y = 1, mu = 0) {
    # generate w, y for everyone
    cov <- get_sigma(r2 = r2, var_y = var_y, var_x = var_y)
    if (cov > var_y) {
        cov <- var_y - 0.2
    }
    Sigma <- matrix(c(var_y, cov, cov, var_y), nrow = 2, byrow = TRUE)
    # note that this is (Y, W)
    n2 <- n * (1 + epsilon)
    dat <- MASS::mvrnorm(n = n2, mu = rep(mu, 2), Sigma = Sigma)
    w <- dat[, 2]
    y <- dat[, 1]
    # generate R
    r1 <- sample(1:n2, n, replace = FALSE)
    r <- ifelse(1:n2 %in% r1, 1, 0)
    newdat <- data.frame(w = w, r = r, y = ifelse(r == 1, y, 0), ystar = y)
    return(newdat)
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
expit <- function(x) exp(x) / (1 + exp(x))

# generate a binary dataset
# @param n the sample size
# @param epsilon the fraction to augment with auxiliary Env sequences
# @param auc the AUC
# @param mu0 the mean in the Y = 0 group
# @param sigma0 the variance in the Y = 0 group
# @param sigma1 the variance in the Y = 1 group
# @return a data.frame (W, R, RY)
gen_data_binary <- function(n = 100, epsilon = 0.2, auc = 0.744, mu0 = -0.32, sigma0 = 0.005, sigma1 = 0.005, p_y = 0.5) {
    mu1 <- get_mu(auc = auc, mu0 = mu0, sigma0 = sigma0, sigma1 = sigma1)
    n2 <- n * (1 + epsilon)
    # generate Y
    y <- rbinom(n = n2, size = 1, prob = p_y)
    # generate W
    w <- vector("numeric", length = n2)
    w[y == 0] <- rnorm(n = n2 - sum(y), mean = mu0, sd = sqrt(sigma0))
    w[y == 1] <- rnorm(n = sum(y), mean = mu1, sd = sqrt(sigma1))
    # generate R
    r1 <- sample(1:n2, size = n, replace = FALSE)
    r <- ifelse(1:n2 %in% r1, 1, 0)
    newdat <- data.frame(w = expit(w), r = r, y = ifelse(r == 1, y, 0), ystar = y)
    return(newdat)
}

# estimate E(Y | W)
est_g <- function(dat, type = "binary") {
    if (grepl("binary", type)) {
        fam <- binomial()
    } else {
        fam <- gaussian()
    }
    g <- glm(y ~ w, data = dat, subset = (dat$r == 1), family = fam)
    preds <- predict(g, newdata = dat, type = "response")
    return(preds)
}
est_theta <- function(dat, preds, lambda, augmented = TRUE) {
    if (!augmented) {
        preds <- preds[dat$r == 1]
        dat <- dat[dat$r == 1, ]
    }
    theta_n <- (1 / lambda) * mean(dat$r * (dat$y - preds), na.rm = TRUE) + mean(preds)
    return(theta_n)
}

# for simulation 1b
est_var <- function(dat, theta, preds, lambda, augmented = TRUE) {
    if (!augmented) {
        y <- dat$y[dat$r == 1]
        r <- rep(1, length(y))
        preds <- preds[dat$r == 1]
    } else {
        y <- dat$y
        r <- dat$r
    }
    eif <- (r / lambda) * (y - theta) - ((r - lambda) / lambda) * (preds - theta)
    return(mean(eif ^ 2, na.rm = TRUE))
}

sim1b_boot_stat <- function(data, indices, augmented, outcome_type, lambda) {
    dat <- data[indices, ]
    if (!augmented) {
        dat <- subset(dat, r == 1)
    }
    y <- dat$y
    r <- dat$r
    w <- dat$w
    if (grepl("binary", outcome_type)) {
        fam <- binomial()
    } else {
        fam <- gaussian()
    }
    g_n <- glm(y ~ w, data = dat, subset = (dat$r == 1), family = fam)
    preds <- predict(g_n, newdata = dat, type = "response")
    theta_n <- (1 / lambda) * mean(dat$r * (dat$y - preds), na.rm = TRUE) + mean(preds)
    return(theta_n)
}

# get one set of estimates
get_ests <- function(mc_id = 1, n = 100, epsilon = 0.2, point_est = 0.2, datatype = "binary", params = list( mu0 = -0.32, sigma0 = 0.005, sigma1 = 0.005, p_y = 0.5), augmented = TRUE) {
    # generate data; a data.frame (W, R, RY)
    # estimate E(Y | W)
    if (datatype == "binary") {
        dat <- gen_data_binary(n = n, epsilon = epsilon, auc = point_est, mu0 = params$mu0, sigma0 = params$sigma0, sigma1 = params$sigma1, p_y = params$p_y)
    } else {
        dat <- gen_data_continuous(n = n, epsilon = epsilon, r2 = point_est, var_y = params$var_y, mu = params$mu)
    }
    ystar <- dat$ystar
    dat <- dat %>% select(-ystar)
    # estimate lambda
    lambda_n <- n / (n * (1 + epsilon))
    # estimate g
    g_n <- est_g(dat = dat, type = datatype)
    # estimate theta
    theta <- est_theta(dat = dat, preds = g_n, lambda = lambda_n, augmented = FALSE)
    theta_aug <- est_theta(dat = dat, preds = g_n, lambda = lambda_n, augmented = TRUE)
    true_mn_y <- mean(ystar)
    # make a little tibble
    ret <- tibble::tibble(mc_id = mc_id, n = n, epsilon = epsilon,
                          point_est = point_est, augmented = c(FALSE, TRUE),
                          datatype = datatype, est = c(theta, theta_aug),
                          truth = true_mn_y)
    return(ret)
}

# for simulation 2 -------------------------------------------------------------
prob_v <- function(v, alpha, gamma, delta) {
    gamma * exp(alpha * (1 - v) - delta * v)
}
objective_function <- function(par, gamma, delta) {
    abs(sum(unlist(lapply(as.list(0:1), function(v) prob_v(v, alpha = par, gamma = gamma, delta = delta)))) - 1)
}
# generate a dataset
# @param n the sample size
# @param J the number of AA positions
# @param gamma the vector of gamma values for the important sites
# @param delta the treatment effect (log hazard ratio, in isolation from genotype)
# @param lambda the non-genotype-specific event rate
# @param eos end of study time
# @param positions the truly important positions
# @return a dataset with outcomes and gp120 AA sites
gen_data_sim2 <- function(n = 1000, J = 100, gamma, delta = log(1 - .29), lambda = -log(.85)/3,
                          eos = 3, positions = c(1:26)) {
    # generate treatment data
    a <- rbinom(n, 1, 0.5)
    n1 <- sum(a == 1)
    n0 <- sum(a == 0)
    # generate outcome data according to the model
    t <- vector("numeric", length = n)
    t[a == 1] <- rexp(n1, rate = lambda * exp(delta))
    t[a == 0] <- rexp(n0, rate = lambda)
    y <- as.numeric(t <= eos)
    t <- pmin(t, eos)
    # generate marks
    v <- matrix(NA, nrow = n, ncol = J)
    for (j in seq_len(J)) {
        v[a == 0, j] <- rbinom(n0, 1, gamma)
        if (j %in% positions) {
            # alpha <- optim(par = 0, fn = objective_function, method = "L-BFGS-B",
            #                gamma = gamma[j], delta = delta)$par
            v[a == 1, j] <- rbinom(n1, 1, gamma[j] * exp((-1)*delta))
        } else {
            v[a == 1, j] <- rbinom(n1, 1, gamma[j])
        }
    }
    # combine
    dat <- data.frame(y = y, t = t, a = a, v)
    return(dat)
}

# run the procedure once
# @param n the sample size
# @param J the number of AA positions
# @param gamma the vector of gamma values for the important sites
# @param delta the treatment effect (log hazard ratio, in isolation from genotype)
# @param lambda the non-genotype-specific event rate
# @param eos end of study time
# @param site_scanning should we look at all sites (TRUE) or the sites specified in 'positions' (FALSE)?
# @param positions a vector of AA positions to look at for sieve analysis
##         (only used in the sieve analysis if site_scanning = FALSE)
# @return whether or not each important site was truly detected
run_sim2_once <- function(mc_id = 1, n = 1000, J = 100, gamma, delta = log(1 - .29),
                          lambda = -log(.85)/3, eos = 3, site_scanning = TRUE,
                          positions = c(1:26)) {
    # generate data
    dat <- gen_data_sim2(n = n, J = J, gamma = gamma, delta = delta, lambda = lambda,
                         eos = eos, positions = positions)
    # do sieve analysis
    if (site_scanning) {
        # run at each AA position
        all_pvals <- vector("numeric", length = J)
        for (j in seq_len(J)) {
            this_mark <- dat %>% pull(!!paste0("X", j))
            this_mark[dat$y == 0] <- NA
            this_result <- sievePH::sievePH(eventTime = dat$t, eventInd = dat$y,
                                            mark = this_mark, tx = dat$a)
            this_summ <- summary(this_result, markGrid = 0)
            all_pvals[j] <- this_summ$pWald.HRconstant.2sided
        }
    } else {
        # run at only the prespecified positions
        all_pvals <- vector("numeric", length = length(positions))
        for (j in seq_len(length(positions))) {
            this_mark <- dat %>% pull(!!paste0("X", positions[j]))
            this_mark[dat$y == 0] <- NA
            this_result <- sievePH::sievePH(eventTime = dat$t, eventInd = dat$y,
                                            mark = this_mark, tx = dat$a)
            this_summ <- summary(this_result, markGrid = 0)
            all_pvals[j] <- this_summ$pWald.HRconstant.2sided
        }
    }
    adj_p <- p.adjust(all_pvals, method = "holm")
    if (site_scanning) {
        ret_p <- adj_p[positions]
    } else {
        ret_p <- adj_p
    }
    ret <- tibble::tibble(mc_id = mc_id, n = n, analysis = ifelse(site_scanning, "Site-scanning", "Priority"),
                          position = positions, p_val = ret_p, reject = ret_p < 0.05)
    return(ret)
}
