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
    int<lower=1> U; // Number of standard deviation parameters
    vector[K - 1] prior_sd; // Prior sd for theta_raw
}

transformed data {
    int combinations[integer_power(2, C - 2), C - 1] = generate_combinations(C - 2);
}

parameters {
    vector[K - 1] theta_raw;
    vector<lower=0, upper=pi()/2>[U] sigma_unif;
    vector[K] standard_normal[R];
}

transformed parameters {
    vector<lower=0>[K] sigma;
    vector[C] XB[R, S];
    vector[K] theta; // sums to zero
    vector[K] beta[R];

    if (U == 1)
    {
        real sigma_value = 2.5 * tan(sigma_unif[1]);
        for (k in 1:K)
            sigma[k] = sigma_value;
    }
    else
        sigma = 2.5 * tan(sigma_unif);

    theta[1] = -sum(theta_raw);
    for (k in 1:(K - 1))
        theta[k + 1] = theta_raw[k];

    for (r in 1:R)
    {
        beta[r] = theta + sigma .* standard_normal[r];
        for (s in 1:S)
            XB[r,s] = X[r,s] * beta[r];
    }
}

model {
    //priors
    for (k in 1:(K - 1))
        theta_raw ~ normal(0, prior_sd[k]);

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
