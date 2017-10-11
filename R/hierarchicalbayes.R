#' @importFrom rstan rstan_options stan extract sampling
#' @importFrom flipChoice ReduceStanFitSize
hierarchicalBayesMaxDiff <- function(dat, n.iterations = 500, n.chains = 8, max.tree.depth = 10,
                                     adapt.delta = 0.8, is.tricked = TRUE, seed = 123,
                                     keep.samples = FALSE, n.classes = 1)
{
    # We want to replace this call with a proper integration of rstan into this package
    require(rstan)

    stan.dat <- createStanData(dat, n.classes, is.tricked)

    # allows Stan chains to run in parallel on multiprocessor machines
    options(mc.cores = parallel::detectCores())

    if (.Platform$OS.type == "unix")
    {
        # Loads a precompiled stan model called mod from sysdata.rda to avoid recompiling.
        # The R code used to generate mod on a linux machine is:
        # mod <- rstan::stan_model(model_code = model.code)
        # devtools::use_data(mod, internal = TRUE, overwrite = TRUE)
        # where model.code is the stan code as a string.
        # Ideally we would want to recompile when the package is built (similar to Rcpp)
        m <- if (n.classes > 1) mod.mix else mod
        stan.fit <- sampling(m, data = stan.dat, chains = n.chains, iter = n.iterations, seed = seed,
                             control = list(max_treedepth = max.tree.depth, adapt_delta = adapt.delta))
    }
    else # windows
    {
        rstan_options(auto_write = TRUE) # saves a compiled Stan object to avoid recompiling next time
        stan.file <- if (n.classes > 1) "exec/mixtureofnormals.stan" else "exec/hb.stan"
        stan.fit <- stan(file = stan.file, data = stan.dat, iter = n.iterations,
                         chains = n.chains, seed = seed,
                         control = list(max_treedepth = max.tree.depth, adapt_delta = adapt.delta))
    }

    result <- list()
    result$respondent.parameters <- computeRespParsStan(stan.fit, dat)
    result$stan.fit <- if (keep.samples) stan.fit else ReduceStanFitSize(stan.fit)
    class(result) <- "FitMaxDiff"
    result
}

createStanData <- function(dat, n.classes, is.tricked)
{
    n.choices <- ncol(dat$X.in)
    n.alternatives <- dat$n.alternatives
    n.respondents <- dat$n
    n.questions <- dat$n.questions
    n.questions.left.in <- dat$n.questions.left.in

    X <- array(dim = c(n.respondents, n.questions.left.in, n.choices, n.alternatives))
    Y.best <- array(1, dim = c(n.respondents, n.questions.left.in))
    Y.worst <- array(n.choices, dim = c(n.respondents, n.questions.left.in))

    for (n in 1:n.respondents)
    {
        for (q in 1:n.questions.left.in)
        {
            X[n, q, , ] <- matrix(0, nrow = n.choices, ncol = n.alternatives)
            for (i in 1:n.choices)
                X[n, q, i, dat$X.in[(n - 1) * n.questions.left.in + q, i]] <- 1
        }
    }

    stan.dat <- list(C = n.choices,
                     K = n.alternatives,
                     R = n.respondents,
                     S = n.questions.left.in,
                     YB = Y.best,
                     YW = Y.worst,
                     X = X,
                     logit_type = ifelse(is.tricked, 1, 2))

    if (n.classes > 1)
        stan.dat$P <- n.classes

    stan.dat
}

computeRespParsStan <- function(stan.fit, dat)
{
    n.respondents <- dat$n
    n.alternatives <- dat$n.alternatives
    alternative.names <- dat$alternative.names
    subset <- dat$subset

    beta <- extract(stan.fit, pars=c("beta"))$beta
    resp.pars.subset <- matrix(NA, n.respondents, n.alternatives)
    if (length(dim(beta)) == 4) # n.classes > 1
    {
        pp <- exp(extract(stan.fit, pars=c("posterior_prob"))$posterior_prob)
        for (i in 1:n.respondents)
            pp[, i, ] <- pp[, i, ] / rowSums(pp[, i, ])

        for (j in 1:n.alternatives)
            for (i in 1:n.respondents)
                resp.pars.subset[i, j] <- mean(beta[, i, , j] * pp[, i, ])
        resp.pars.subset
    }
    else
        resp.pars.subset <- colMeans(beta, dims = 1)

    result <- matrix(NA, nrow = length(subset), ncol = ncol(resp.pars.subset))
    result[subset, ] <- resp.pars.subset
    colnames(result) <- alternative.names
    result
}
