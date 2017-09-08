#' @importFrom bayesm rhierMnlRwMixture
hierarchicalBayesGibbsMaxDiff <- function(dat, n.draws, seed = 123, n.classes = 1)
{
    set.seed(seed)

    n.choices <- ncol(dat$X.in)
    n.alternatives <- dat$n.alternatives
    n.respondents <- dat$n
    n.questions.left.in <- dat$n.questions.in
    n.qc <- n.questions.left.in * n.choices

    data <- list()
    for (i in 1:n.respondents)
    {
        X <- matrix(0, nrow = 2 * n.qc, ncol = n.alternatives - 1)
        for (j in 1:n.questions.left.in)
        {
            X.row <- dat$X.in[(i - 1) * n.questions.left.in + j, ]
            for (k in 1:n.choices)
            {
                ind <- (j - 1) * n.choices + ((k + i - 2) %% n.choices) + 1
                if (X.row[k] > 1)
                {
                    X[ind, X.row[k] - 1] <- 1
                    X[ind + n.qc, X.row[k] - 1] <- -1
                }
            }
        }
        best.i <- ((i - 1) %% n.choices) + 1
        worst.i <- ((i + n.choices - 2) %% n.choices) + 1
        y <- c(rep(best.i, n.questions.left.in), rep(worst.i, n.questions.left.in))
        data[[i]] <- list(y = y, X = X)
    }

    out <- invisible(rhierMnlRwMixture(Data = list(lgtdata = data, p = n.choices),
                             Prior = list(ncomp = n.classes),
                             Mcmc = list(R = n.draws, nprint = 0)))

    resp.pars <- matrix(NA, nrow = n.respondents, ncol = n.alternatives)
    draws.to.use <- ceiling(0.5 * n.draws):n.draws
    for (i in 1:n.respondents)
        resp.pars[i, ] <- c(0, rowMeans(out$betadraw[i, , draws.to.use]))
    respondent.parameters <- matrix(NA, nrow = length(dat$subset), ncol = ncol(resp.pars))
    respondent.parameters[dat$subset, ] <- resp.pars
    colnames(respondent.parameters) <- dat$alternative.names

    result <- list(respondent.parameters = respondent.parameters, bayesm.output = out)
    class(result) <- "FitMaxDiff"
    result
}
