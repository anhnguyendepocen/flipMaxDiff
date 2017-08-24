data {
    int<lower=2> C; // Number of alternatives (choices) in each scenario
    int<lower=1> K; // Number of alternatives
    int<lower=1> R; // Number of respondents
    int<lower=1> S; // Number of scenarios per respondent
    int<lower=1,upper=C> YB[R, S]; // best choices
    int<lower=1,upper=C> YW[R, S]; // worst choices
    matrix[C, K] X[R, S]; // matrix of attributes for each obs
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
    vector[K] theta;
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
            YB[r,s] ~ categorical_logit(XB[r,s]);
            YW[r,s] ~ categorical_logit(-XB[r,s]);
        }
    }
}
