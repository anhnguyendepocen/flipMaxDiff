data {
    int<lower=2> C; // Number of alternatives (choices) in each scenario
    int<lower=1> K; // Number of alternatives
    int<lower=1> R; // Number of respondents
    int<lower=1> S; // Number of scenarios per respondent
    int<lower=1> P; // Number of classes
    int<lower=1,upper=C> YB[R, S]; // best choices
    int<lower=1,upper=C> YW[R, S]; // worst choices
    matrix[C, K] X[R, S]; // matrix of attributes for each obs
}

parameters {
    vector[K - 1] theta_raw[P];
    cholesky_factor_corr[K] L_omega[P];
    vector<lower=0, upper=pi()/2>[K] sigma_unif[P];
    vector[K] standard_normal[R, P];
    simplex[P] class_weights;
}

transformed parameters {
    vector<lower=0>[K] sigma[P];
    matrix[K, K] L_sigma[P];
    vector[K] theta[P]; // sums to zero
    vector[K] beta[R, P];
    vector[P] posterior_prob[R];

    for (p in 1:P)
    {
        sigma[p] = 2.5 * tan(sigma_unif[p]);
        L_sigma[p] = diag_pre_multiply(sigma[p], L_omega[p]);

        theta[p, 1] = -sum(theta_raw[p]);
        for (k in 1:(K - 1))
            theta[p, k + 1] = theta_raw[p, k];

        for (r in 1:R)
            beta[r, p] = theta[p] + L_sigma[p] * standard_normal[r, p];
    }

    for (r in 1:R)
    {
        for (p in 1:P)
        {
            real prob = log(class_weights[p]);
            for (s in 1:S)
            {
                prob = prob + categorical_logit_lpmf(YB[r,s] | X[r, s] * beta[r, p]);
                prob = prob + categorical_logit_lpmf(YW[r,s] | -X[r, s] * beta[r, p]);
            }
            posterior_prob[r, p] = prob;
        }
    }
}

model {
    for (p in 1:P)
    {
        theta_raw[p] ~ normal(0, 10);
        L_omega[p] ~ lkj_corr_cholesky(4);
        for (r in 1:R)
            standard_normal[r, p] ~ normal(0, 1);
    }

    //likelihood
    for (r in 1:R)
        target += log_sum_exp(posterior_prob[r]);
}
