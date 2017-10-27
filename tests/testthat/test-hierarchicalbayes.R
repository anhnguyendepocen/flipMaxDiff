context("Hierarchical Bayes")

# We currently leave out unit tests due to issues with getting it to run successfully.
# We would like to get unit tests working when rstan is properly integrated into this package.
# Hierarchical Bayes is currently tested via our internal regression testing system.

tech.data <- suppressWarnings(foreign::read.spss("http://wiki.q-researchsoftware.com/images/f/f1/Technology_2017.sav", to.data.frame = TRUE))
tech.design <- read.csv("http://wiki.q-researchsoftware.com/images/7/78/Technology_MaxDiff_Design.csv")
best <- tech.data[, c("Q5a_left", "Q5b_left", "Q5c_left", "Q5d_left", "Q5e_left", "Q5f_left")]
worst <- tech.data[, c("Q5a_right", "Q5b_right", "Q5c_right", "Q5d_right", "Q5e_right", "Q5f_right")]
names <- c("Apple", "Microsoft", "IBM", "Google", "Intel", "Samsung", "Sony", "Dell", "Yahoo", "Nokia")

test_that("HB", {
    result <- FitMaxDiff(design = tech.design, version = rep(1, nrow(best)), best = best, worst = worst,
                         alternative.names = names, algorithm = "HB", hb.iterations = 200, hb.chains = 2)
    expect_error(print(result), NA)
    expect_equal(result$in.sample.accuracy, 0.933222958057395)
})

test_that("HB cross validation", {
    result <- FitMaxDiff(design = tech.design, version = rep(1, nrow(best)), best = best, worst = worst,
                         alternative.names = names, algorithm = "HB", hb.iterations = 200, hb.chains = 2,
                         tasks.left.out = 2)
    expect_error(print(result), NA)
    expect_equal(result$out.sample.accuracy, 0.640728476821192)
})

test_that("HB filter", {
    sub <- unclass(tech.data$Q2) <= 3
    result <- FitMaxDiff(design = tech.design, version = rep(1, nrow(best)), best = best, worst = worst,
                         alternative.names = names, algorithm = "HB", hb.iterations = 200, hb.chains = 2,
                         subset = sub)
    expect_error(print(result), NA)
    expect_equal(result$in.sample.accuracy, 0.908602150537634)
})

test_that("HB weights", {
    wgt <- tech.data$RESPNUM
    expect_error(FitMaxDiff(design = tech.design, version = rep(1, nrow(best)), best = best, worst = worst,
                         alternative.names = names, algorithm = "HB", hb.iterations = 200, hb.chains = 2,
                         weights = wgt), "Weights are not able to be applied for Hierarchical Bayes.")
})

test_that("HB rank-ordered logit", {
    result <- FitMaxDiff(design = tech.design, version = rep(1, nrow(best)), best = best, worst = worst,
                         alternative.names = names, algorithm = "HB", hb.iterations = 200, hb.chains = 2,
                         is.tricked = FALSE)
    expect_error(print(result), NA)
    expect_equal(result$in.sample.accuracy, 0.922185430463576)
})

test_that("HB 2 classes", {
    result <- FitMaxDiff(design = tech.design, version = rep(1, nrow(best)), best = best, worst = worst,
                         alternative.names = names, algorithm = "HB", hb.iterations = 200, hb.chains = 2,
                         n.classes = 2)
    expect_error(print(result), NA)
    expect_equal(result$in.sample.accuracy, 0.931567328918322)
})

test_that("HB diagonal", {
    result <- FitMaxDiff(design = tech.design, version = rep(1, nrow(best)), best = best, worst = worst,
                         alternative.names = names, algorithm = "HB", hb.iterations = 200, hb.chains = 2,
                         normal.covariance = "Diagonal")
    expect_error(print(result), NA)
    expect_equal(result$in.sample.accuracy, 0.931567328918322)
})

test_that("HB diagonal 2 classes", {
    result <- FitMaxDiff(design = tech.design, version = rep(1, nrow(best)), best = best, worst = worst,
                         alternative.names = names, algorithm = "HB", hb.iterations = 200, hb.chains = 2,
                         normal.covariance = "Diagonal", n.classes = 2)
    expect_error(print(result), NA)
    expect_equal(result$in.sample.accuracy, 0.929911699779249)
})

test_that("HB spherical", {
    result <- FitMaxDiff(design = tech.design, version = rep(1, nrow(best)), best = best, worst = worst,
                         alternative.names = names, algorithm = "HB", hb.iterations = 200, hb.chains = 2,
                         normal.covariance = "Spherical")
    expect_error(print(result), NA)
    expect_equal(result$in.sample.accuracy, 0.934878587196468)
})

