# utility functions for simulation 1: assess relative efficiency of using extra Env sequences

# for simulation 0 -------------------------------------------------------------
# get a datset of size n
# @param data the dataset (VRC01 catnap data)
# @param n the sample size
# @return a dataset of size n sampled from the VRC01 catnap data
subsample_dataset <- function(data, n = 10) {
  samp <- sample(1:nrow(data), size = n, replace = FALSE)
  return(data[samp, ])
}

# get CV-MSE of a lasso predicting the outcome
# @param data the dataset
# @param K the number of folds
# @return the CV-MSE
get_lasso_cvmse <- function(data, K = 5) {
  mses <- vector("numeric", length = K)
  y <- data$ic80
  x <- as.matrix(data %>% select(-ic80))
  cc_y <- y[complete.cases(y)]
  cc_x <- x[complete.cases(y), ]
  folds <- sample(rep(1:K, length = nrow(cc_x)))
  for (k in 1:K) {
    this_fit <- cv.glmnet(x = cc_x[folds != k, ], y = cc_y[folds != k], nfolds = 5)
    these_preds <- predict(this_fit, newx = cc_x[folds == k, ], s = this_fit$lambda.min)
    mses[k] <- mean((cc_y[folds == k] - these_preds) ^ 2)
  }
  return(mean(mses))
}

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
# @param var1 the variance of y (log10 IC-80)
# @param var2 the variance of the PAR scores
# @param mu1 the mean of y
# @param mu2 the mean of the PAR scores
# @return a data.frame (W, R, RY)
gen_data_continuous <- function(n = 100, epsilon = 0.2, r2 = 0.345, var1 = 1, mu1 = 0,
                                var2 = 1, mu2 = 0) {
    # generate w, y for everyone
    cov <- get_sigma(r2 = r2, var_y = var1, var_x = var2)
    if (cov > var1) {
        cov <- var1 - 0.2
    }
    Sigma <- matrix(c(var1, cov, cov, var2), nrow = 2, byrow = TRUE)
    # note that this is (Y, W)
    n2 <- n * (1 + epsilon)
    dat <- MASS::mvrnorm(n = n2, mu = c(mu1, mu2), Sigma = Sigma)
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
estimating_equation <- function(dat, preds, theta, lambda) {
  term <- ifelse(is.na(dat$y), (-1) * (dat$r - lambda) / lambda * (preds - theta),
                 dat$r / lambda * (dat$y - theta))
  sum(term)
}

# sim1b_boot_stat <- function(data, indices, augmented, outcome_type, lambda) {
sim1b_boot_stat <- function(data, indices, augmented, outcome_type, g_n) {
  dat <- data[indices, ]
  lambda <- sum(dat$r == 1) / nrow(dat)
  if (!augmented) {
    dat <- subset(dat, r == 1)
  } 
  y <- dat$y
  y2 <- ifelse(is.na(y), 0, y)
  r <- dat$r
  w <- dat$w
  if (grepl("binary", outcome_type)) {
    fam <- binomial()
  } else {
    fam <- gaussian()
  }
  # g_n <- glm(y ~ w, data = dat, subset = (dat$r == 1), family = fam)
  preds <- predict(g_n, newdata = dat, type = "response")
  theta_n <- (1 / lambda) * mean(r * (y2 - preds), na.rm = TRUE) + mean(preds)
  return(theta_n)
}

# get one set of estimates
get_ests <- function(mc_id = 1, n = 100, epsilon = 0.2, point_est = 0.2, 
                     datatype = "binary", params = list( mu0 = -0.32, sigma0 = 0.005, sigma1 = 0.005, p_y = 0.5), 
                     augmented = TRUE) {
    # generate data; a data.frame (W, R, RY)
    # estimate E(Y | W)
    if (datatype == "binary") {
        dat <- gen_data_binary(n = n, epsilon = epsilon, auc = point_est, mu0 = params$mu0, sigma0 = params$sigma0, sigma1 = params$sigma1, p_y = params$p_y)
        dat2 <- gen_data_binary(n = 1e6, epsilon = epsilon, auc = point_est, mu0 = params$mu0, sigma0 = params$sigma0, sigma1 = params$sigma1, p_y = params$p_y)
    } else {
        dat <- gen_data_continuous(n = n, epsilon = epsilon, r2 = point_est, var1 = params$var1, mu1 = params$mu1,
                                   var2 = params$var2, mu2 = params$mu2)
        dat2 <- gen_data_continuous(n = 1e6, epsilon = epsilon, r2 = point_est, var1 = params$var1, mu1 = params$mu1,
                                   var2 = params$var2, mu2 = params$mu2)
    }
    true_mn_y <- mean(dat2$ystar)
    ystar <- dat$ystar
    dat <- dat %>% select(-ystar)
    # estimate lambda
    lambda_n <- n / (n * (1 + epsilon))
    # estimate g
    g_n <- est_g(dat = dat, type = datatype)
    # estimate theta
    theta <- est_theta(dat = dat, preds = g_n, lambda = lambda_n, augmented = FALSE)
    theta_aug <- est_theta(dat = dat, preds = g_n, lambda = lambda_n, augmented = TRUE)
    theta_oracle <- mean(ystar)
    # make a little tibble
    ret <- tibble::tibble(mc_id = mc_id, n = n, epsilon = epsilon,
                          point_est = point_est, estimator = c("non-augmented", "augmented", "oracle"),
                          datatype = datatype, est = c(theta, theta_aug, theta_oracle),
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
get_pe_1 <- function(pe_overall, pe_0, gamma) {
  1 - exp( (log(1 - pe_overall) - (1 - gamma) * log(1 - pe_0)) / gamma)
}
# generate a dataset
# @param n the sample size
# @param all_positions the positions on Env
# @param gamma_0 proportions of resistant genotypes at each site in the placebo arm
# @param gamma_1 proportions of resistant genotypes at each site in the VRC01 arm
# @param beta the treatment effect (log hazard ratio, in isolation from genotype)
# @param p0 the mean non-genotype-specific event rate
# @param q the probability of censoring
# @param eos end of study time
# @param positions the truly important positions
# @return a dataset with outcomes and gp120 AA sites
# gen_data_sim2_v1 <- function(n = 1000, all_positions = 1:100, gamma_0, gamma_1, p0 = 67 / 1535,
#                           beta = log(1 - 0.18), q = 769 / 4611,
#                           eos = 80, positions = c(1:26)) {
#     # generate treatment data
#     a <- rbinom(n, 1, 0.5)
#     n1 <- sum(a == 1)
#     n0 <- sum(a == 0)
#     # generate censoring
#     C <- runif(n, min = 0, max = eos / q)
#     # generate outcome data according to the model
#     t <- vector("numeric", length = n)
#     rate_placebo <- (-1) * log(1 - p0) / eos
#     t[a == 1] <- rexp(n1, rate = rate_placebo * exp(beta))
#     t[a == 0] <- rexp(n0, rate = rate_placebo)
#     # determine observed endpoints
#     y <- as.numeric(t <= eos & t <= C)
#     n01 <- sum(a == 0 & y == 1)
#     n11 <- sum(a == 1 & y == 1)
#     t <- pmin(t, eos)
#     # generate marks; note that if we don't have an endpoint (and thus a sequence),
#     # the marks are NA
#     J <- length(all_positions)
#     v <- matrix(NA, nrow = n, ncol = J)
#     for (j in seq_len(J)) {
#         v[a == 0 & y == 1, j] <- rbinom(n01, 1, gamma_0[j])
#         v[a == 1 & y == 1, j] <- rbinom(n11, 1, gamma_1[j])
#     }
#     colnames(v) <- paste0("X", all_positions)
#     v_df <- data.frame(v)
#     # combine
#     dat <- data.frame(y = y, t = t, a = a, v)
#     return(dat)
# }
gen_data_sim2_simple <- function(n = 1000, lambda_0 = 0.018, 
                                 alpha = log(1 - .7),
                                 q = .1, eos = 365 * 2) {
    # generate treatment data
    a <- rbinom(n, 1, 0.5)
    # generate censoring
    C <- runif(n, min = 0, max = eos / q)
    t <- rexp(n, rate = lambda_0 * exp(a * alpha)) * 365 
    # determine observed endpoints
    y <- as.numeric(t <= eos & t <= C)
    obstime <- pmin(t, C, eos)
    delta <- (obstime == C)
    dat <- data.frame(y = y, t = obstime, a = a, delta = delta)
    return(dat)
}
# generate a dataset
# @param n the sample size
# @param all_positions the positions on Env
# @param gamma_0 proportions of resistant genotypes at each site in the placebo arm
# @param pe_max the maximum PE of VRC01 for a sensitive virus
# @param beta the treatment effect (log hazard ratio, in isolation from genotype)
# @param alpha a vector of genotype-specific effects (an interaction between treatment and
#              presence of the sensitive genotype at a given residue)
# @param lambda_0 the incidence rate (per person-year) in the AMP placebo arm
# @param q the probability of censoring
# @param eos the final study date (in days)
# @param positions the truly important positions
# @param position the position(s) of interest that actually have sieve effects
# @return a dataset with outcomes and gp120 AA sites
gen_data_sim2 <- function(n = 1000, all_positions = 1:100, gamma_0, lambda_0 = 3.04,
                          pe_max = 0.95, beta = log(1 - 0.18), pe_0 = 0, pe_1 = .89,
                          q = 769 / 4611, eos = 365 * 2,
                          positions = c(1:26),
                          position = 1) {
    # generate treatment data
    a <- rbinom(n, 1, 0.5)
    # generate censoring
    C <- runif(n, min = 0, max = eos / q)
    # generate latent cause-specific infection times
    t0 <- rexp(n, rate = lambda_0 * (1 - gamma_0[position]) * exp(a * log(1 - pe_0))) * 365 # hazard among those with S230 = 0
    t1 <- rexp(n, rate = lambda_0 * gamma_0[position] * exp(a * log(1 - pe_1))) * 365 # hazard among those with S230 = 1
    t <- pmin(t0, t1)
    # generate latent mark variable denoting sensitive vs other at each position
    S <- t(replicate(n, rbinom(length(all_positions), 1, gamma_0)))
    S[, position] <- as.numeric(t == t1)
    # determine observed endpoints
    y <- as.numeric(t <= eos & t <= C)
    obstime <- pmin(t, C, eos)
    delta <- (obstime == C)
    # generate marks; note that if we don't have an endpoint (and thus a sequence),
    # the marks are NA
    J <- length(all_positions)
    v <- matrix(NA, nrow = n, ncol = J)
    v[y == 1, ] <- S[y == 1, ]
    colnames(v) <- paste0("X", all_positions)
    v_df <- data.frame(v)
    # combine
    dat <- data.frame(y = y, t = obstime, a = a, delta = delta, v)
    return(dat)
}

# compute the number of events at each time point,
# return the closest time point that yields the desired number of events
# @param t a vector of survival times
# @param n_events the number of events we're stopping at
get_time_point <- function(t = rep(1, 100), C = rep(1, 100), n_events = 88) {
    all_y <- do.call(cbind, lapply(1:200 * 7, function(x) as.numeric(t <= x & t <= C)))
    time_point <- which.min(abs(colSums(all_y) - n_events))
    return(time_point)
}

# get p-value from Lunn & McNeil test
# @param dat the dataset
# @param mark the mark value
lunn_mcneil_pval <- function(dat, mark, package = "sievePH") {
    if (package == "sievePH") {
        this_sievePH <- sievePH::sievePH(eventTime = dat$t, eventInd = dat$y,
                                         mark = mark, tx = dat$a)
        pval <- summary(this_sievePH, markGrid = 0)$pWald.HRconstant.2sided
    } else {
        this_result <- lunnMcneilTest(flrtime = dat$t, flrstatus = dat$y,
                                      flrtype = mark + 1, Vx = dat$a)
        pval <- this_result$coefficients[3, 5]
    }
    return(pval)
}

run_sim2_simple_once <- function(mc_id = 1, n = 1000, lambda_0 = 0.018,
                                 alpha = log(1 - .7), q = .1, eos = 365 * 2) {
    dat <- gen_data_sim2_simple(n = n, lambda_0 = lambda_0, alpha = alpha, q = q, eos = eos)
    surv_fit <- coxph(Surv(time = t, event = y) ~ a, data = dat)
    ret <- tibble::tibble(mc_id = mc_id, n = n, 
                          hr_p = summary(surv_fit)$waldtest[3],
                          n_events = sum(dat$y),
                          n_0 = sum(dat$y[dat$a == 0]),
                          n_1 = sum(dat$y[dat$a == 1]))
}
# run the procedure once
# @param n the sample size
# @param all_positions the unique Env positions
# @param gammas a tibble with the positions, proportions of resistant genotypes at each site in the placebo arm, and
#               proportions of resistant genotypes at each site in the VRC01 arm
# @param beta the treatment effect (log hazard ratio, in isolation from genotype)
# @param alpha a vector of genotype-specific effects (an interaction between treatment and
#              presence of the sensitive genotype at a given residue)
# @param lambda_0 the incidence rate (per person-year) in the AMP placebo arm
# @param q the probability of censoring
# @param eos the final study date (in days)
# @param positions the truly important positions (only used in the sieve analysis if site_scanning = FALSE)
# @param position the position(s) of interest that actually have sieve effects
# @param site_scanning should we look at all sites (TRUE) or the sites specified in 'positions' (FALSE)?
# @param pe_0 prevention efficacy when S230 = 0
# @param pe_1 prevention efficacy when S230 = 1
# @param minvar_screen how many observed marks do we need for stability
# @param debug return a simple summary, for debugging
# @param package which code to use for testing?
# @return whether or not each important site was truly detected
run_sim2_once <- function(mc_id = 1, n = 1000, all_positions = 1:100, gammas,
                          lambda_0 = 3.04, beta = log(1 - 0.18),
                          q = 769 / 4611, eos = 365 * 2, site_scanning = TRUE, 
                          pe_0 = 0, pe_1 = .89,
                          positions = c(1:26), position = 1,
                          minvar_screen = 10, debug = FALSE, package = "none") {
    # generate data
    dat <- gen_data_sim2(n = n, all_positions = all_positions, lambda_0 = lambda_0,
                         gamma_0 = gammas$prop_placebo,
                         beta = beta, pe_0 = pe_0, pe_1 = pe_1, q = q,
                         eos = eos, positions = positions, position = position)
    cc_dat <- dat[complete.cases(dat), ]
    y_summ <- data.frame(mc_id = mc_id, n = sum(dat$y), n00 = sum(cc_dat$y[cc_dat$a == 0 & cc_dat$X230 == 0]),
                         n01 = sum(cc_dat$y[cc_dat$a == 0 & cc_dat$X230 == 1]),
                         n10 = sum(cc_dat$y[cc_dat$a == 1 & cc_dat$X230 == 0]),
                         n11 = sum(cc_dat$y[cc_dat$a == 1 & cc_dat$X230 == 1]))
    if (debug) {
        # get number of events
        return(y_summ)
    }
    # survival analysis for primary analysis
    surv_fit <- survival::coxph(Surv(time = t, event = y) ~ a, data = dat)
    overall_p <- summary(surv_fit)$waldtest[3]
    # apply minimum variability filter
    n_cases_with_each_mark <- apply(dat[dat$y == 1, 4:ncol(dat)], 2, function(x) sum(dat$y[x == 1], na.rm = TRUE))
    dat2 <- dat[, c(1:3, which(n_cases_with_each_mark >= minvar_screen) + 3)]
    analysis_positions <- names(dat2)[4:ncol(dat2)]
    minvar_positions <- positions[positions %in% (which(n_cases_with_each_mark >= minvar_screen) + 3)]
    # do sieve analysis
    if (site_scanning) {
        # run at each AA position
        all_pvals <- vector("numeric", length = length(analysis_positions))
        for (j in seq_len(length(analysis_positions))) {
            this_mark <- dat2 %>% pull(!!analysis_positions[j])
            this_mark[dat2$y == 0] <- NA
            if (sum(this_mark[!is.na(this_mark)] == 0) <= minvar_screen |
                sum(this_mark[!is.na(this_mark)] == 1) <= minvar_screen) {
                all_pvals[j] <- 1
            } else {
                all_pvals[j] <- lunn_mcneil_pval(dat = dat2, mark = this_mark, package = package)
            }
        }
    } else {
        # run at only the prespecified positions
        all_pvals <- vector("numeric", length = length(minvar_positions))
        for (j in seq_len(length(minvar_positions))) {
            this_mark <- dat2 %>% pull(!!paste0("X", minvar_positions[j]))
            this_mark[dat2$y == 0] <- NA
            if (sum(this_mark[!is.na(this_mark)] == 0) <= minvar_screen |
                sum(this_mark[!is.na(this_mark)] == 1) <= minvar_screen) {
                all_pvals[j] <- 1
            } else {
                all_pvals[j] <- lunn_mcneil_pval(dat = dat2, mark = this_mark, package = package)
            }
        }
    }
    p_order <- order(all_pvals, decreasing = FALSE)
    ordered_p <- all_pvals[p_order]
    adj_p_init <- p.adjust(ordered_p, method = "holm")
    adj_p <- adj_p_init[order(p_order)]
    if (site_scanning) {
        # ret_p <- adj_p[positions]
        ret_p <- adj_p[analysis_positions == paste0("X", 230)]
        unadjusted_p <- all_pvals[analysis_positions == paste0("X", 230)]
    } else {
        # ret_p <- adj_p
        ret_p <- adj_p[positions == 230]
        unadjusted_p = all_pvals[positions == 230]
    }
    ret <- tibble::tibble(mc_id = mc_id, n = n, analysis = ifelse(site_scanning, "Site-scanning", "Priority"),
                          pe_0 = pe_0,
                          position = 230, p_val = ret_p, reject = ret_p < 0.05, 
                          unadjusted_p_val = unadjusted_p,
                          hr_p = overall_p,
                          n_events = y_summ$n, n00 = y_summ$n00, n01 = y_summ$n01,
                          n10 = y_summ$n10, n11 = y_summ$n11)
    return(ret)
}
