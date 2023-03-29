# boosted algorithms

# gradient boosting function using h2o package
SL.h2oboost <- function(Y, X, newX, family, obsWeights = rep(1, length(Y)), ...)

{
    SuperLearner:::.SL.require("h2o")

    # Set GBM parameters to test
    hyper_parameters <- list(ntrees = list(1000),
                              max_depth = list(2,4,5,6),
                              learn_rate = list(0.05, 0.1, 0.2),
                              col_sample_rate = list(0.1, 0.2, 0.3))

    # Bind vector of outcome and covariate matrix together;
    # if Y is binary, make it a factor first
    Y <- switch((family$family == "binomial") + 1, Y, as.factor(Y))
    dat <- cbind(Y, X)

    # Convert dat to an h2o object
    dat.hex <- as.h2o(dat)

    # set up GBM hyperparameters
    if (family$family == "binomial") {
        h2o_dist <- "bernoulli"
        h2o_metric <- "AUC"
    } else if (family$family == "gaussian") {
        h2o_dist <- "gaussian"
        h2o_metric <- "MSE"
    } else {
        stop("The entered family isn't currently supported. Please enter one of 'binomial' or 'gaussian'.")
    }

    # search over the grid
    gbm.model <- h2o::h2o.grid("gbm",
                               hyper_params = hyper_parameters,
                               training_frame = dat.hex,
                               y = "Y",
                               distribution = h2o_dist,
                               nfolds = 5,
                               balance_classes = (h2o_dist == "bernoulli"),
                               max_after_balance_size = 5,
                               fold_assignment = ifelse(h2o_dist == "bernoulli",
                                                        "Stratified", "AUTO"),
                               stopping_metric = toupper(h2o_metric),
                               stopping_rounds = 3,
                               stopping_tolerance = 0.001,
                               max_runtime_secs = 60,
                               parallelism = 0)

    # get the models from the grid and sort by metric
    grid <- h2o::h2o.getGrid(gbm.model@grid_id,
                             sort_by = tolower(h2o_metric),
                             decreasing = (h2o_dist == "bernoulli"))

    # Save best parameters
    best.max_depth <- as.numeric(grid@summary_table[1, ]$max_depth)
    best.learn_rate <- as.numeric(grid@summary_table[1, ]$learn_rate)
    best.col_sample_rate <- as.numeric(grid@summary_table[1, ]$col_sample_rate)

    # Remove all models in grid to save memory
    h2o.removeAll(retained_elements = c(dat.hex))
    rm(gbm.model, grid)

    # Call garbage collection
    h2o:::.h2o.garbageCollect()
    h2o:::.h2o.garbageCollect()
    h2o:::.h2o.garbageCollect()

    # Train the model with best hyperparameters
    gbm.final.model <- h2o::h2o.gbm(training_frame = dat.hex,
                                    y = "Y",
                                    distribution = h2o_dist,
                                    stopping_metric = toupper(h2o_metric),
                                    stopping_rounds = 3,
                                    stopping_tolerance = 0.001,
                                    ntrees = 1000,
                                    max_depth = best.max_depth,
                                    learn_rate = best.learn_rate,
                                    col_sample_rate = best.col_sample_rate)

    # Convert newdata to h2o object
    newX.hex <- as.h2o(newX)
    # Get predictions
    pred.raw <- h2o::h2o.predict(object = gbm.final.model,
                              newdata = newX.hex)

    # Extract predicted probabilities
    if(family$family == "gaussian"){
      pred <- as.numeric(as.vector(pred.raw))
    }
    else if(family$family == "binomial"){
      pred <- as.numeric(as.vector(pred.raw$p1))
    }
    # Make fit an object with the filepath we need to reload the h2o object
    fit <- list(object = gbm.final.model)
    class(fit) <- c("SL.h2oboost")
    out = list(pred = pred, fit = fit)

    return(out)
}

# Predict method for h2oboost
predict.SL.h2oboost <- function(object, newdata, ...)
{
    SuperLearner:::.SL.require("h2o")
    L <- list(...)
    # convert data to h2o object
    newdata.hex <- h2o::as.h2o(newdata)
    # Get predictions
    pred.raw <- h2o::h2o.predict(object = object$object,
                              newdata = newdata.hex)
    # Extract predicted probabilites
    if (L$family$family == "gaussian"){
      pred <- as.numeric(as.vector(pred.raw))
    }
    else if (L$family$family == "binomial"){
      pred <- as.numeric(as.vector(pred.raw$p1))
    }
    pred
}
descr_SL.h2oboost <- "boosted regression trees with (maximum depth, learning rate, column sampling rate) selected by 5-fold CV over the grid $(2, 4, 5, 6)\\times(.05, .1, .2)\\times(.1, .2, .3)$"

SL.xgboost.corrected <- function (Y, X, newX, family, obsWeights = rep(1, length(Y)), id, ntrees = 1000,
    max_depth = 4, shrinkage = 0.1, minobspernode = 10, params = list(),
    nthread = 1, verbose = 0, save_period = NULL, ...)
{
    SuperLearner:::.SL.require("xgboost")
    if (packageVersion("xgboost") < 0.6)
        stop("SL.xgboost requires xgboost version >= 0.6, try help('SL.xgboost') for details")
    if (!is.matrix(X)) {
        X = model.matrix(~. - 1, X)
    }
    xgmat = xgboost::xgb.DMatrix(data = X, label = Y, weight = obsWeights)
    if (family$family == "gaussian") {
        model = xgboost::xgboost(data = xgmat, objective = "reg:squarederror",
            nrounds = ntrees, max_depth = max_depth, min_child_weight = minobspernode,
            eta = shrinkage, verbose = verbose, nthread = nthread,
            params = params, save_period = save_period)
    }
    if (family$family == "binomial") {
        model = xgboost::xgboost(data = xgmat, objective = "binary:logistic",
            nrounds = ntrees, max_depth = max_depth, min_child_weight = minobspernode, eval_metric = "logloss",
            eta = shrinkage, verbose = verbose, nthread = nthread,
            params = params, save_period = save_period)
    }
    if (family$family == "multinomial") {
        model = xgboost::xgboost(data = xgmat, objective = "multi:softmax",
            nrounds = ntrees, max_depth = max_depth, min_child_weight = minobspernode,
            eta = shrinkage, verbose = verbose, num_class = length(unique(Y)),
            nthread = nthread, params = params, save_period = save_period)
    }
    if (!is.matrix(newX)) {
        newX = model.matrix(~. - 1, newX)
    }
    pred = predict(model, newdata = newX)
    fit = list(object = model)
    class(fit) = c("SL.xgboost")
    out = list(pred = pred, fit = fit)
    return(out)
}

# SL.xgboost.2 <- function(..., max_depth = 2){
# 	SL.xgboost.corrected(..., max_depth = max_depth)
# }
SL.xgboost.4 <- function(..., max_depth = 4){
	SL.xgboost.corrected(..., max_depth = max_depth)
}
# SL.xgboost.6 <- function(..., max_depth = 6){
# 	SL.xgboost.corrected(..., max_depth = max_depth)
# }
SL.xgboost.8 <- function(..., max_depth = 8){
	SL.xgboost.corrected(..., max_depth = max_depth)
}
SL.xgboost.12 <- function(..., max_depth = 12){
    SL.xgboost.corrected(..., max_depth = max_depth)
}
descr_SL.xgboost <- "boosted regression trees with maximum depth of "
descr_SL.xgboost.2 <- paste0(descr_SL.xgboost, 2)
descr_SL.xgboost.4 <- paste0(descr_SL.xgboost, 4)
descr_SL.xgboost.6 <- paste0(descr_SL.xgboost, 6)
descr_SL.xgboost.8 <- paste0(descr_SL.xgboost, 8)
descr_SL.xgboost.12 <- paste0(descr_SL.xgboost, 12)

# random forests
SL.ranger.imp <- function (Y, X, newX, family, obsWeights = rep(1, length(Y)), num.trees = 500, mtry = floor(sqrt(ncol(X))),
    write.forest = TRUE, probability = family$family == "binomial",
    min.node.size = ifelse(family$family == "gaussian", 5, 1),
    replace = TRUE, sample.fraction = ifelse(replace, 1, 0.632),
    num.threads = 1, verbose = TRUE, ...) {
    SuperLearner:::.SL.require("ranger")
    if (family$family == "binomial") {
        Y = as.factor(Y)
    }
    if (is.matrix(X)) {
        X = data.frame(X)
    }
    fit <- ranger::ranger(`_Y` ~ ., data = cbind(`_Y` = Y, X),
        num.trees = num.trees, mtry = mtry, min.node.size = min.node.size,
        replace = replace, sample.fraction = sample.fraction,
        case.weights = obsWeights, write.forest = write.forest,
        probability = probability, num.threads = num.threads,
        verbose = verbose, importance = "impurity")
    pred <- predict(fit, data = newX)$predictions
    if (family$family == "binomial") {
        pred = pred[, "1"]
    }
    fit <- list(object = fit, verbose = verbose)
    class(fit) <- c("SL.ranger")
    out <- list(pred = pred, fit = fit)
    return(out)
}
SL.ranger.reg <- function(..., X, mtry = floor(sqrt(ncol(X)))){
	SL.ranger.imp(..., X = X, mtry = mtry)
}

SL.ranger.small <- function(..., X, mtry = floor(sqrt(ncol(X)) * 1/2)){
	SL.ranger.imp(..., X  = X, mtry = mtry)
}

SL.ranger.large <- function(..., X, mtry = floor(sqrt(ncol(X)) * 2)){
	SL.ranger.imp(..., X = X, mtry = mtry)
}
descr_SL.ranger.imp <- "random forest with `mtry` equal to "
descr_SL.ranger.reg <- paste0(descr_SL.ranger.imp, "square root of number of predictors")
descr_SL.ranger.small <- paste0(descr_SL.ranger.imp, "one-half times square root of number of predictors")
descr_SL.ranger.large <- paste0(descr_SL.ranger.imp, "two times square root of number of predictors")

# function used to do smarter CV for glmnet
get_fold_id <- function(Y){
  fold_id <- rep(0, length(Y))
  wiY0 <- which(Y == 0)
  wiY1 <- which(Y == 1)
  #if <4 cases, no cv
  if(length(wiY1) == 4){
    #if exactly 4 cases, 4-fold cv
    #1 case per fold
    fold <- 1:4
    fold_id[sample(wiY1)] <- fold
    fold_id[sample(wiY0)] <- rep(fold, length = length(wiY0))
  }else{
    #if >=5 cases, 5 fold cv
    #cases split as evenly as possible
    fold <- 1:5
    fold_id[sample(wiY1)] <- rep(fold, length = length(wiY1))
    fold_id[sample(wiY0)] <- rep(fold, length = length(wiY0))
  }
  return(fold_id)
}


# function to have more robust behavior in SL.glmnet
SL.glmnet.0 <- function (Y, X, newX, family, obsWeights = rep(1, length(Y)), id, alpha = 1, nfolds = 5,
    nlambda = 100, useMin = TRUE, loss = "deviance", ...) {
    SuperLearner:::.SL.require("glmnet")
    if (!is.matrix(X)) {
        X <- model.matrix(~-1 + ., X)
        newX <- model.matrix(~-1 + ., newX)
    }

    if(family$family == "binomial"){
        fold_id <- get_fold_id(Y)
        nfolds <- max(fold_id)
        if(nfolds != 0){
            fitCV <- glmnet::cv.glmnet(x = X, y = Y, weights = obsWeights,
                lambda = NULL, type.measure = loss, nfolds = nfolds,
                foldid = fold_id, family = family$family, alpha = alpha, nlambda = nlambda,
                ...)
            pred <- predict(fitCV, newx = newX, type = "response", s = ifelse(useMin,
                "lambda.min", "lambda.1se"))
            fit <- list(object = fitCV, useMin = useMin)
            class(fit) <- "SL.glmnet"
        }else{
            # if fewer than 3 cases, just use mean
            meanY <- weighted.mean(Y, w = obsWeights)
            pred <- rep.int(meanY, times = nrow(newX))
            fit <- list(object = meanY)
            out <- list(pred = pred, fit = fit)
            class(fit) <- c("SL.mean")
        }
    }else{
        fitCV <- glmnet::cv.glmnet(x = X, y = Y, weights = obsWeights,
                lambda = NULL, type.measure = loss, nfolds = nfolds,
                family = family$family, alpha = alpha, nlambda = nlambda,
                ...)
            pred <- predict(fitCV, newx = newX, type = "response", s = ifelse(useMin,
                "lambda.min", "lambda.1se"))
            fit <- list(object = fitCV, useMin = useMin)
            class(fit) <- "SL.glmnet"
    }
    out <- list(pred = pred, fit = fit)
    return(out)
}

# lasso
SL.glmnet.50 <- function(..., alpha = 0.5){
	SL.glmnet.0(..., alpha = alpha)
}
SL.glmnet.25 <- function(..., alpha = 0.25){
	SL.glmnet.0(..., alpha = alpha)
}
SL.glmnet.75 <- function(..., alpha = 0.75){
	SL.glmnet.0(..., alpha = alpha)
}

descr_SL.glmnet <- "elastic net with $\\lambda$ selected by 5-fold CV and $\\alpha$ equal to "
descr_SL.glmnet.50 <- paste0(descr_SL.glmnet, "0.5")
descr_SL.glmnet.25 <- paste0(descr_SL.glmnet, "0.25")
descr_SL.glmnet.75 <- paste0(descr_SL.glmnet, "0.75")
descr_SL.glmnet.0 <- paste0(descr_SL.glmnet, "1 (i.e., the lasso)")

descr_SL.mean <- "intercept only regression"
descr_SL.glm <- "main terms generalized linear model"

#' Temporary fix for convex combination method mean squared error
#' Relative to existing implementation, we reduce the tolerance at which
#' we declare predictions from a given algorithm the same as another
tmp_method.CC_LS <- function ()
{
    computeCoef = function(Z, Y, libraryNames, verbose, obsWeights,
        errorsInLibrary = NULL, ...) {
        cvRisk <- apply(Z, 2, function(x) mean(obsWeights * (x -
            Y)^2))
        names(cvRisk) <- libraryNames
        compute <- function(x, y, wt = rep(1, length(y))) {
            wX <- sqrt(wt) * x
            wY <- sqrt(wt) * y
            D <- crossprod(wX)
            d <- crossprod(wX, wY)
            A <- cbind(rep(1, ncol(wX)), diag(ncol(wX)))
            bvec <- c(1, rep(0, ncol(wX)))
            fit <- tryCatch({quadprog::solve.QP(Dmat = D, dvec = d, Amat = A,
                bvec = bvec, meq = 1)
          }, error = function(e){
            out <- list()
            class(out) <- "error"
            out
          })
            invisible(fit)
        }
        modZ <- Z
        naCols <- which(apply(Z, 2, function(z) {
            all(z == 0)
        }))
        anyNACols <- length(naCols) > 0
        if (anyNACols) {
            warning(paste0(paste0(libraryNames[naCols], collapse = ", "),
                " have NAs.", "Removing from super learner."))
        }
        tol <- 4
        dupCols <- which(duplicated(round(Z, tol), MARGIN = 2))
        anyDupCols <- length(dupCols) > 0
        if (anyDupCols) {
            warning(paste0(paste0(libraryNames[dupCols], collapse = ", "),
                " are duplicates of previous learners.", " Removing from super learner."))
        }
        if (anyDupCols | anyNACols) {
            rmCols <- unique(c(naCols, dupCols))
            modZ <- Z[, -rmCols, drop = FALSE]
        }
        fit <- compute(x = modZ, y = Y, wt = obsWeights)
        if(class(fit) != "error"){
          coef <- fit$solution
        }else{
          coef <- rep(0, ncol(Z))
          coef[which.min(cvRisk)] <- 1
        }
        if (anyNA(coef)) {
            warning("Some algorithms have weights of NA, setting to 0.")
            coef[is.na(coef)] = 0
        }
        if(class(fit) != "error"){
          if (anyDupCols | anyNACols) {
              ind <- c(seq_along(coef), rmCols - 0.5)
              coef <- c(coef, rep(0, length(rmCols)))
              coef <- coef[order(ind)]
          }
          coef[coef < 1e-04] <- 0
          coef <- coef/sum(coef)
        }
        if (!sum(coef) > 0)
            warning("All algorithms have zero weight", call. = FALSE)
        list(cvRisk = cvRisk, coef = coef, optimizer = fit)
    }
    computePred = function(predY, coef, ...) {
        predY %*% matrix(coef)
    }
    out <- list(require = "quadprog", computeCoef = computeCoef,
        computePred = computePred)
    invisible(out)
}


#' Temporary fix for convex combination method negative log-likelihood loss
#' Relative to existing implementation, we reduce the tolerance at which
#' we declare predictions from a given algorithm the same as another.
#' Note that because of the way \code{SuperLearner} is structure, one needs to
#' install the optimization software separately.
tmp_method.CC_nloglik <- function ()
{
    computePred = function(predY, coef, control, ...) {
        if (sum(coef != 0) == 0) {
            stop("All metalearner coefficients are zero, cannot compute prediction.")
        }
        stats::plogis(trimLogit(predY[, coef != 0], trim = control$trimLogit) %*%
            matrix(coef[coef != 0]))
    }
    computeCoef = function(Z, Y, libraryNames, obsWeights, control,
        verbose, ...) {
        tol <- 4
        dupCols <- which(duplicated(round(Z, tol), MARGIN = 2))
        anyDupCols <- length(dupCols) > 0
        modZ <- Z
        if (anyDupCols) {
            warning(paste0(paste0(libraryNames[dupCols], collapse = ", "),
                " are duplicates of previous learners.", " Removing from super learner."))
            modZ <- modZ[, -dupCols, drop = FALSE]
        }
        modlogitZ <- trimLogit(modZ, control$trimLogit)
        logitZ <- trimLogit(Z, control$trimLogit)
        cvRisk <- apply(logitZ, 2, function(x) -sum(2 * obsWeights *
            ifelse(Y, stats::plogis(x, log.p = TRUE), stats::plogis(x, log.p = TRUE,
                lower.tail = FALSE))))
        names(cvRisk) <- libraryNames
        obj_and_grad <- function(y, x, w = NULL) {
            y <- y
            x <- x
            function(beta) {
                xB <- x %*% cbind(beta)
                loglik <- y * stats::plogis(xB, log.p = TRUE) + (1 -
                  y) * stats::plogis(xB, log.p = TRUE, lower.tail = FALSE)
                if (!is.null(w))
                  loglik <- loglik * w
                obj <- -2 * sum(loglik)
                p <- stats::plogis(xB)
                grad <- if (is.null(w))
                  2 * crossprod(x, cbind(p - y))
                else 2 * crossprod(x, w * cbind(p - y))
                list(objective = obj, gradient = grad)
            }
        }
        lower_bounds = rep(0, ncol(modZ))
        upper_bounds = rep(1, ncol(modZ))
        if (anyNA(cvRisk)) {
            upper_bounds[is.na(cvRisk)] = 0
        }
        r <- tryCatch({nloptr::nloptr(x0 = rep(1/ncol(modZ), ncol(modZ)),
            eval_f = obj_and_grad(Y, modlogitZ), lb = lower_bounds,
            ub = upper_bounds, eval_g_eq = function(beta) (sum(beta) -
                1), eval_jac_g_eq = function(beta) rep(1, length(beta)),
            opts = list(algorithm = "NLOPT_LD_SLSQP", xtol_abs = 1e-08))
        }, error = function(e){
          out <- list()
          class(out) <- "error"
          out
        })
        if (r$status < 1 || r$status > 4) {
            warning(r$message)
        }
        if(class(r) != "error"){
          coef <- r$solution
        }else{
          coef <- rep(0, ncol(Z))
          coef[which.min(cvRisk)] <- 1
        }
        if (anyNA(coef)) {
            warning("Some algorithms have weights of NA, setting to 0.")
            coef[is.na(coef)] <- 0
        }
        if (anyDupCols) {
            ind <- c(seq_along(coef), dupCols - 0.5)
            coef <- c(coef, rep(0, length(dupCols)))
            coef <- coef[order(ind)]
        }
        coef[coef < 1e-04] <- 0
        coef <- coef/sum(coef)
        out <- list(cvRisk = cvRisk, coef = coef, optimizer = r)
        return(out)
    }
    list(require = "nloptr", computeCoef = computeCoef, computePred = computePred)
}

# based on user-inputted options, make a vector describing the super learner library
make_sl_library_vector <- function(opts){
    default_library <- NULL
    # check if h2o boost is requested
    if ("h2oboost" %in% opts$learners){
        default_library <- c(default_library, "SL.h2oboost")
    }
    # check if rf is requested
    if ("rf" %in% opts$learners) {
      if (opts$cvtune) {
        default_library <- c(default_library, "SL.ranger.small", "SL.ranger.reg", "SL.ranger.large")
      } else {
        if (length(opts$learners) == 1) {
            default_library <- "SL.ranger.reg"
        } else {
            default_library <- c(default_library, "SL.ranger.reg")
        }
      }
    }
    # check if xgboost is requested
    if("xgboost" %in% opts$learners){
      if (opts$cvtune) {
        # default_library <- c(default_library, "SL.xgboost.2", "SL.xgboost.4", "SL.xgboost.6", "SL.xgboost.8")
        default_library <- c(default_library, "SL.xgboost.4", "SL.xgboost.8", "SL.xgboost.12")
      } else {
        if (length(opts$learners) == 1) {
            default_library <- "SL.xgboost.4"
        } else {
            default_library <- c(default_library, "SL.xgboost.4")
        }
      }
    }
    # check if elastic net is requested
    if("lasso" %in% opts$learners){
      if (opts$cvtune) {
        default_library <- c(default_library, "SL.glmnet.0", "SL.glmnet.25", "SL.glmnet.50", "SL.glmnet.75")
      } else {
        if (length(opts$learners) == 1) {
            default_library <- "SL.glmnet.0"
        } else {
            default_library <- c(default_library, "SL.glmnet.0")
        }
      }
    }

    # check if screening is required
    if(!all(opts$var_thresh == 0)){
        # read into memory
        for(i in opts$var_thresh){
            make_screen_wrapper(var_thresh = i)
        }
        # include in list
        which_xgboost <- c(grep("xgboost", default_library), grep("h2oboost", default_library))
        if(length(which_xgboost) > 0){
            default_library_for_grid <- default_library[-which_xgboost]
        }else{
            default_library_for_grid <- default_library
        }
        learn_grid <- expand.grid(learner=default_library_for_grid, screen=paste0("var_thresh_", opts$var_thresh))
        if(length(which_xgboost) > 0){
            learn_grid <- rbind(learn_grid, data.frame(learner=default_library[which_xgboost], screen="All"))
        }
        default_library <- as.list(as.data.frame(t(learn_grid), stringsAsFactors = FALSE))
    }

    # if fitting a super learner, throw in SL.mean
    if(all(opts$var_thresh == 0)){
        if(length(opts$learners) > 1){
          default_library <- c(default_library, "SL.mean")
        }
    }else{
        if(length(opts$learners) > 1){
          default_library <- c(default_library, list(c("SL.mean", "All")))
        }
    }
    return(default_library)
}

var_thresh_general <- function(Y, X, family, obsWeights, var_thresh, ...){
    include <- rep(TRUE, ncol(X))
    n <- length(X[,1])
    for(i in 1:ncol(X)){
        if(all(X[,i]) %in% c(0,1)){
            sum_Xi <- sum(X[,i])
            if(sum_Xi < var_thresh | sum_Xi > n - var_thresh){
                include[i] <- FALSE
            }
        }
    }
    return(include)
}

#' read a screen wrapper into global environment
make_screen_wrapper <- function(var_thresh){
    eval(parse(text=paste0(
        "var_thresh_", var_thresh, "<<-",
          "function(..., var_thresh = ", as.numeric(var_thresh), "){",
             "var_thresh_general(..., var_thresh = var_thresh)",
          "}"
    )))
}

#' function to run super learner and cv super learner on a single outcome
#' @param complete_dat the full dataset
#' @param outcome_name String name of outcome
#' @param pred_names Vector of string names of predictor variables
#' @param opts List of options outputted by get_global_options()
#' @param save_dir name of directory to save results to
#' @param fit_name name of fits (defaults to fit_<outcome_name>.rds)
#' @param cv_fit_name name of CV fits (defaults to cvfit_<outcome_name>.rds)
#' @param save_full_object Flag for whether or not to save the full fitted object, or just the fitted values
#' @param full_fit is it the full fit (TRUE) or a reduced fit (FALSE)?
#' @param SL.library the library to pass to SuperLearner
#' @param ... additional arguments to pass to individual algorithms or the SuperLearner
sl_one_outcome <- function(complete_dat, outcome_name,
                           pred_names,
                           opts,
                           h2o_here,
                           save_dir = "/home/slfits/",
                           fit_name = paste0("fit_", outcome_name, ".rds"),
                           cv_fit_name = paste0("cvfit_", outcome_name, ".rds"),
                           save_full_object = TRUE,
                           full_fit = TRUE,
                           ss_folds = rep(1, opts$nfolds),
                           SL.library = "SL.mean",
                           call_out = FALSE,
                           ...){
  # three cases to worry about in terms of how to deal with missing data
  # 1. same_subset requested, but only studying ic80 or only studying ic50-derived outcomes,
  #    in which case, we will ignore your same_subset request because you're only studying
  #    endpoints based entirely on ic50 or entirely on ic80
  # 2. same_subset requested and studying mix of ic50(-derived) and ic80 endpoints, in which
  #    case we will subset down to sequences with both ic50 and ic80
  # 3. no same_subset requested, in which case we use all data available on each outcome
  if(!opts$same_subset | !(("ic80" %in% opts$outcomes | "iip" %in% opts$outcomes) & length(opts$outcomes) > 1)){
    complete_cases_idx <- complete.cases(complete_dat[,c(outcome_name,pred_names)])
  } else {
    complete_cases_idx <- complete.cases(complete_dat)
  }
  # subset data to only complete outcome and pred_names
  dat <- complete_dat[complete_cases_idx, ]
  pred <- dat[ , pred_names]
  # grab args
  L <- list(...)
  non_cv_L <- L
  non_cv_L$cvControl$validRows <- NULL
  # make names for the full learner, cv learner, fitted values, cv fitted values, etc.
  learner_name <- gsub("fit_", "learner_", fit_name)
  fitted_name <- gsub("fit_", "fitted_", fit_name)
  cv_fitted_name <- gsub("cvfit_", "cvfitted_", cv_fit_name)
  cv_preds_name <- gsub("cvfitted_", "cvpreds_", gsub("cvfit_", "cvpreds_", cv_fit_name))
  cv_folds_name <- gsub("cvfitted_", "cvfolds_", gsub("cvfit_", "cvfolds_", cv_fit_name))
  se_ss <- rep(ifelse(full_fit, 1, 2), L$cvControl$V)
  cv_learner_name <- gsub("cvfit_", "cvlearner_", cv_fit_name)
  # unless we do not want to tune using CV nor evaluate performance using CV,
  # we will make a call to super learner
  if (!(opts$cvtune == FALSE & opts$cvperf == FALSE)) {
    if(call_out){
        print(paste0("Fitting super learner with", SL.library, collapse = ","))
    }else{
        fit <- SuperLearner(Y = dat[ , outcome_name], X = pred, SL.library = SL.library, family = non_cv_L$family, method = non_cv_L$method, cvControl = non_cv_L$cvControl)
        # since cv super learner saves this as output and we need some parsimony later...
        fit$Y <- dat[ , outcome_name]
        if (save_full_object) {
            saveRDS(fit, file = paste0(save_dir, fit_name))
            saveRDS(fit, file = paste0(save_dir, learner_name))
            if (h2o_here){
            # save h2o models
            h2o.saveModel(object = fit$fitLibrary$SL.h2oboost_All$object,
                              path = paste0(save_dir),
                              force = TRUE)
            }
        }
        if (length(opts$learners) > 1) {
          # save super learner predictions
          saveRDS(fit$SL.predict, file = paste0(save_dir, fitted_name))
          # save super learner weights
          saveRDS(fit$coef, file = paste0(save_dir, "slweights_", fit_name))
        } else if (length(opts$learners) == 1) {
            # save learner predictions (these are cv-selected if cvtune = TRUE)
            saveRDS(fit$library.predict[, which.min(fit$cvRisk)], file = paste0(save_dir, fitted_name))
            # save learner
            if (save_full_object) {
                saveRDS(fit$fitLibrary[[which.min(fit$cvRisk)]]$object, file = paste0(save_dir, learner_name))
            }
        } else {
            # don't save anything
        }
    }
  } else {
    # if more than one learner is specified or more than one var_thresh is specified,
    # then we're fitting a super learner
    if (length(opts$learners) > 1 | length(opts$var_thresh) > 1) {
        if(call_out){
            print(paste0("Fitting SuperLearner with ", SL.library, collapse = ","))
        }else{
            fit <- SuperLearner(Y = dat[ , outcome_name], X = pred, SL.library = SL.library, family = non_cv_L$family, method = non_cv_L$method, cvControl = non_cv_L$cvControl)
            fit$Y <- dat[ , outcome_name]
            if (save_full_object) {
                saveRDS(fit, file = paste0(save_dir, fit_name))
                saveRDS(fit, file = paste0(save_dir, learner_name))
                if (h2o_here){
                    # save h2o models
                    h2o.saveModel(object = fit$fitLibrary$SL.h2oboost_All$object,
                                  path = paste0(save_dir),
                                  force = TRUE)
                }
            }
            # save super learner predictions
            saveRDS(fit$SL.predict, file = paste0(save_dir, fitted_name))
            # save super learner weights
            saveRDS(fit$coef, file = paste0(save_dir, "slweights_", fit_name))
        }
    } else {
        if(all(opts$var_thresh == 0)){
            include <- rep(TRUE, ncol(X))
            this_learner <- SL.library
        }else{
            include <- do.call(paste0("var_thresh_", opts$var_thresh), args = list(X = pred))
            this_learner <- SL.library[[1]][1]
        }
        if(call_out){
            print(paste0("Fitting single algorithm ", ifelse(all(include), "with no screening", "with screening")))
        }else{
            fit <- do.call(SL.library, args = c(L[!grepl("cvControl", names(L)) & !grepl("method", names(L))], list(Y = dat[ , outcome_name], X = pred[,include, drop = FALSE], newX = pred[,include, drop = FALSE])))
            saveRDS(fit$pred, file = paste0(save_dir, fitted_name))
            # this will be an object with class native to what the individual learner is
            # i.e., if rf is desired, it'll be ranger object
            if (save_full_object) {
                saveRDS(fit$fit$object, file = paste0(save_dir, fit_name))
                saveRDS(fit$fit$object, file = paste0(save_dir, learner_name))
                if (h2o_here){
                    # save h2o models
                    h2o.saveModel(object = fit$fit$object,
                                  path = paste0(save_dir),
                                  force = TRUE)
                }
            }
        }
    }
  }
  # now copy cvControl to innerCvControl to pass to CV.SL
  inner_cv_control <- L$cvControl
  inner_cv_control$validRows <- NULL
  L$innerCvControl <- list(inner_cv_control)
  new_arg_list <- c(list(Y = dat[ , outcome_name], X = pred, SL.library = SL.library), L)
  if (length(opts$learners) == 1 & opts$cvtune & opts$cvperf) {
    # in this case, we want to report back only results for discrete super learner
    # since we're trying to pick out e.g., the best single random forest fit
    # note that if either opts$cvtune AND/OR opts$cvperf is FALSE then everything we need
    # will already be in the SuperLearner fit object.
    if(call_out){
        print(paste0("Fitting CV.SuperLearner with ", SL.library, collapse = ","))
    }else{
        cv_fit <- do.call(CV.SuperLearner, new_arg_list)
        if (save_full_object) {
            saveRDS(cv_fit, file = paste0(save_dir, cv_fit_name))
            saveRDS(cv_fit, file = paste0(save_dir, cv_learner_name))
        }
        saveRDS(cv_fit$discreteSL.predict, file = paste0(save_dir, cv_fitted_name))
        saveRDS(cv_fit$folds, file = paste0(save_dir, cv_folds_name))
        vimp_cf_folds <- vimp::get_cv_sl_folds(cv_fit$folds)
        cv_preds <- vimp::extract_sampled_split_predictions(cvsl_obj = cv_fit, sample_splitting = TRUE, sample_splitting_folds = ss_folds, full = full_fit)
        saveRDS(cv_preds, file = paste0(save_dir, cv_preds_name))
    }
  } else if ((length(opts$learners) > 1 | length(opts$var_thresh) > 1) & opts$cvperf) {
    if(call_out){
        print(paste0("Fitting CV.SuperLearner with ", SL.library, collapse = ","))
    }else{
        # if multiple learners or multiple var_thresh then we need a CV superlearner
        # unless cvperf = FALSE
        cv_fit <- do.call(CV.SuperLearner, new_arg_list)
        if (save_full_object) {
            saveRDS(cv_fit, file = paste0(save_dir, cv_fit_name))
            saveRDS(cv_fit, file = paste0(save_dir, cv_learner_name))
        }
        saveRDS(cv_fit$SL.predict, file = paste0(save_dir, cv_fitted_name))
        saveRDS(cv_fit$folds, file = paste0(save_dir, cv_folds_name))
        vimp_cf_folds <- vimp::get_cv_sl_folds(cv_fit$folds)
        cv_preds <- vimp::extract_sampled_split_predictions(cvsl_obj = cv_fit, sample_splitting = TRUE, sample_splitting_folds = ss_folds, full = full_fit)
        saveRDS(cv_preds, file = paste0(save_dir, cv_preds_name))
    }
} else {
    # do nothing
}
  return(invisible(NULL))
}
