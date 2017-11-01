context("HB Gibbs")

tech.data <- suppressWarnings(foreign::read.spss("http://wiki.q-researchsoftware.com/images/f/f1/Technology_2017.sav", to.data.frame = TRUE))
tech.design <- read.csv("http://wiki.q-researchsoftware.com/images/7/78/Technology_MaxDiff_Design.csv")
best <- tech.data[, c("Q5a_left", "Q5b_left", "Q5c_left", "Q5d_left", "Q5e_left", "Q5f_left")]
worst <- tech.data[, c("Q5a_right", "Q5b_right", "Q5c_right", "Q5d_right", "Q5e_right", "Q5f_right")]
names <- c("Apple", "Microsoft", "IBM", "Google", "Intel", "Samsung", "Sony", "Dell", "Yahoo", "Nokia")

test_that("Gibbs", {
    result <- FitMaxDiff(design = tech.design, version = rep(1, nrow(best)), best = best, worst = worst,
                         alternative.names = names, algorithm = "HB-bayesm", hb.iterations = 1000)
    expect_error(print(result), NA)
})

test_that("Gibbs mixture of normals", {
    result <- FitMaxDiff(design = tech.design, version = rep(1, nrow(best)), best = best, worst = worst,
                         alternative.names = names, algorithm = "HB-bayesm", hb.iterations = 1000, n.classes = 2)
    expect_error(print(result), NA)
})

test_that("Gibbs filter", {
    sub <- unclass(tech.data$Q2) <= 3
    result <- FitMaxDiff(design = tech.design, version = rep(1, nrow(best)), best = best, worst = worst,
                         alternative.names = names, algorithm = "HB-bayesm", hb.iterations = 1000,
                         subset = sub)
    expect_error(print(result), NA)
})

test_that("Gibbs weights", {
    wgt <- tech.data$RESPNUM
    expect_error(FitMaxDiff(design = tech.design, version = rep(1, nrow(best)), best = best, worst = worst,
                         alternative.names = names, algorithm = "HB-bayesm", hb.iterations = 1000,
                         weights = wgt), "Weights are not able to be applied for Hierarchical Bayes.")
})
