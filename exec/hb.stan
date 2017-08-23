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
    vector[K - 1] Theta_raw;
    cholesky_factor_corr[K] L_Omega;
    vector<lower=0, upper=pi()/2>[K] L_sigma_unif;
    vector[K] standard_normal[R];
}

transformed parameters {
    vector<lower=0>[K] L_sigma;
    matrix[K, K] L_Sigma;
    vector[C] XB[R, S];
    vector[K] Theta;
    vector[K] Beta[R];

    L_sigma = 2.5 * tan(L_sigma_unif);
    L_Sigma = diag_pre_multiply(L_sigma, L_Omega);

    Theta[1] = -sum(Theta_raw);
    for (k in 1:(K - 1))
        Theta[k + 1] = Theta_raw[k];

    for (r in 1:R)
    {
        Beta[r] = Theta + L_Sigma * standard_normal[r];
        for (s in 1:S)
            XB[r,s] = X[r,s] * Beta[r];
    }
}

model {
    //priors
    Theta_raw ~ normal(0, 10);
    L_Omega ~ lkj_corr_cholesky(4);

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
