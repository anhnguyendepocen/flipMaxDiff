#' @importFrom rstan rstan_options stan extract sampling
#' @importFrom flipChoice ReduceStanFitSize ComputeRespPars
hierarchicalBayesMaxDiff <- function(dat, n.iterations = 500, n.chains = 8, max.tree.depth = 10,
                                     adapt.delta = 0.8, is.tricked = TRUE, seed = 123,
                                     keep.samples = FALSE, n.classes = 1, include.stanfit = TRUE,
                                     normal.covariance = "Full", is.unit.test = FALSE)
{
    # We want to replace this call with a proper integration of rstan into this package
    require(rstan)

    stan.dat <- createStanData(dat, n.classes, is.tricked, normal.covariance)

    # allows Stan chains to run in parallel on multiprocessor machines
    options(mc.cores = parallel::detectCores())

    if (is.unit.test)
    {
        stan.fit <- stan(model_code = modelCode(), data = stan.dat, iter = n.iterations,
                         chains = n.chains, seed = seed,
                         control = list(max_treedepth = max.tree.depth, adapt_delta = adapt.delta))
    }
    else if (.Platform$OS.type == "unix")
    {
        # Loads a precompiled stan model called mod from sysdata.rda to avoid recompiling.
        # The R code used to generate mod on a linux machine is:
        # mod <- rstan::stan_model(model_code = model.code)
        # devtools::use_data(mod, internal = TRUE, overwrite = TRUE)
        # where model.code is the stan code as a string.
        # Ideally we would want to recompile when the package is built (similar to Rcpp)
        m <- stanModel(n.classes, normal.covariance)
        stan.fit <- sampling(m, data = stan.dat, chains = n.chains, iter = n.iterations, seed = seed,
                             control = list(max_treedepth = max.tree.depth, adapt_delta = adapt.delta))
    }
    else # windows
    {
        rstan_options(auto_write = TRUE) # saves a compiled Stan object to avoid recompiling next time
        stan.file <- stanFileName(n.classes, normal.covariance)
        stan.fit <- stan(file = stan.file, data = stan.dat, iter = n.iterations,
                         chains = n.chains, seed = seed,
                         control = list(max_treedepth = max.tree.depth, adapt_delta = adapt.delta))
    }

    result <- list()
    result$respondent.parameters <- ComputeRespPars(stan.fit, dat$alternative.names, dat$subset)
    if (include.stanfit)
    {
        result$stan.fit <- if (keep.samples) stan.fit else ReduceStanFitSize(stan.fit)
        result$beta.draws <- extract(stan.fit, pars=c("beta"))$beta
    }
    class(result) <- "FitMaxDiff"
    result
}

createStanData <- function(dat, n.classes, is.tricked, normal.covariance)
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

    if (normal.covariance == "Diagonal")
        stan.dat$U <- n.alternatives
    else if (normal.covariance == "Spherical")
        stan.dat$U <- 1

    stan.dat
}

stanFileName <- function(n.classes, normal.covariance)
{
    if (n.classes == 1)
    {
        if (normal.covariance == "Full")
            "exec/hb.stan"
        else
            "exec/diagonal.stan"
    }
    else
    {
        if (normal.covariance == "Full")
            "exec/mixtureofnormals.stan"
        else
            "exec/diagonalmixture.stan"
    }
}

stanModel <- function(n.classes, normal.covariance)
{
    if (n.classes == 1)
    {
        if (normal.covariance == "Full")
            mod
        else
            mod.diag
    }
    else
    {
        if (normal.covariance == "Full")
            mod.mix
        else
            mod.mix.diag
    }
}

modelCode <- function()
{
    "
functions {
    real tricked_logit_lpmf(int[] y, vector xb)
    {
    real result = 0;
    result = result + categorical_logit_lpmf(y[1] | xb);
    result = result + categorical_logit_lpmf(y[2] | -xb);
    return result;
    }

    real rank_ordered_logit_lpmf(int[] y, vector xb, int[,] combinations)
    {
    int n_combinations = dims(combinations)[1];
    int n_choices = rows(xb);
    real sum_combinations = 0;
    vector[n_choices - 2] others;
    int c = 1;

    for (k in 1:n_choices)
    {
    if (k != y[1] && k != y[2])
    {
    others[c] = xb[k];
    c = c + 1;
    }
    }

    for (i in 1:n_combinations)
    {
    real log_phi_dot_omega;
    int m = combinations[i, 1];
    log_phi_dot_omega = log_sum_exp(others[combinations[i, 2:(m + 1)]]);
    sum_combinations = sum_combinations + (-1) ^ m / (1 + exp(log_phi_dot_omega - xb[y[2]]));
    }
    if (sum_combinations <= 0)
    sum_combinations = 1e-16;

    return categorical_logit_lpmf(y[1] | xb) + log(sum_combinations);
    }

    // Stan doesn't have a power function for positive powers that returns int
    int integer_power(int a, int b)
    {
    int result = 1;
    for (i in 1:b)
    result = result * a;
    return result;
    }

    int[] next_combination(int[] comb)
    {
    int n_items = size(comb);
    int next_comb[n_items] = comb;
    for (i in 1:n_items) {
    next_comb[i] = -comb[i];
    if (next_comb[i] == 1)
    break;
    }
    return next_comb;
    }

    int[,] generate_combinations(int n_items)
    {
    int n_combinations = integer_power(2, n_items);
    int combinations[n_combinations, n_items + 1]; // first column contains counts
    int comb[n_items];
    int c = 1;
    int j;
    for (i in 1:n_items)
    comb[i] = -1;

    while (c <= n_combinations) {
    j = 2;
    for (i in 1:n_items)
    {
    if (comb[i] == 1)
    {
    combinations[c, j] = i;
    j = j + 1;
    }
    }
    combinations[c, 1] = j - 2;
    comb = next_combination(comb);
    c = c + 1;
    }
    return combinations;
    }
}

data {
int<lower=2> C; // Number of alternatives (choices) in each question
int<lower=1> K; // Number of alternatives
int<lower=1> R; // Number of respondents
int<lower=1> S; // Number of questions per respondent
int<lower=1,upper=C> YB[R, S]; // best choices
int<lower=1,upper=C> YW[R, S]; // worst choices
matrix[C, K] X[R, S]; // matrix of attributes for each obs
int logit_type; // 1: tricked logit, 2: rank-ordered logit
}

transformed data {
int combinations[integer_power(2, C - 2), C - 1] = generate_combinations(C - 2);
}

parameters {
vector[K - 1] theta_raw;
cholesky_factor_corr[K] L_omega;
vector<lower=0, upper=pi()/2>[K] sigma_unif;
vector[K] standard_normal[R];
}

transformed parameters {
vector<lower=0>[K] sigma;
matrix[K, K] L_sigma;
vector[C] XB[R, S];
vector[K] theta; // sums to zero
vector[K] beta[R];

sigma = 2.5 * tan(sigma_unif);
L_sigma = diag_pre_multiply(sigma, L_omega);

theta[1] = -sum(theta_raw);
for (k in 1:(K - 1))
theta[k + 1] = theta_raw[k];

for (r in 1:R)
{
    beta[r] = theta + L_sigma * standard_normal[r];
    for (s in 1:S)
    XB[r,s] = X[r,s] * beta[r];
}
}

model {
//priors
theta_raw ~ normal(0, 10);
L_omega ~ lkj_corr_cholesky(4);

for (r in 1:R)
standard_normal[r] ~ normal(0, 1);

//likelihood
for (r in 1:R) {
for (s in 1:S) {
int Y[2];
Y[1] = YB[r, s];
Y[2] = YW[r, s];
if (logit_type == 1)
Y ~ tricked_logit(XB[r, s]);
else
Y ~ rank_ordered_logit(XB[r, s], combinations);
}
}
}
"
}
