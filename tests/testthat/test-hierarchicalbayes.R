context("Hierarchical Bayes")

tech.data <- suppressWarnings(foreign::read.spss("http://wiki.q-researchsoftware.com/images/f/f1/Technology_2017.sav", to.data.frame = TRUE))
tech.design <- read.csv("http://wiki.q-researchsoftware.com/images/7/78/Technology_MaxDiff_Design.csv")
best <- tech.data[, c("Q5a_left", "Q5b_left", "Q5c_left", "Q5d_left", "Q5e_left", "Q5f_left")]
worst <- tech.data[, c("Q5a_right", "Q5b_right", "Q5c_right", "Q5d_right", "Q5e_right", "Q5f_right")]
names <- c("Apple", "Microsoft", "IBM", "Google", "Intel", "Samsung", "Sony", "Dell", "Yahoo", "Nokia")

test_that("HB", {
    result <- FitMaxDiff(design = tech.design, version = rep(1, nrow(best)),
                         best = best, worst = worst,
                         alternative.names = names, algorithm = "HB",
                         hb.iterations = 10, hb.chains = 1,
                         hb.warnings = FALSE)
    expect_error(print(result), NA)
    # expect_equal(result$in.sample.accuracy, 0.933222958057395)
})

test_that("HB cross validation", {
    result <- FitMaxDiff(design = tech.design, version = rep(1, nrow(best)),
                         best = best, worst = worst, alternative.names = names,
                         algorithm = "HB", hb.iterations = 10, hb.chains = 1,
                         tasks.left.out = 2, hb.warnings = FALSE)
    expect_error(print(result), NA)
    # expect_equal(result$out.sample.accuracy, 0.640728476821192)
})

test_that("HB filter", {
    sub <- unclass(tech.data$Q2) <= 3
    result <- FitMaxDiff(design = tech.design, version = rep(1, nrow(best)),
                         best = best, worst = worst, alternative.names = names,
                         algorithm = "HB", hb.iterations = 10, hb.chains = 1,
                         subset = sub, hb.warnings = FALSE)
    expect_error(print(result), NA)
    # expect_equal(result$in.sample.accuracy, 0.908602150537634)
})

test_that("HB weights", {
    wgt <- tech.data$RESPNUM
    expect_error(FitMaxDiff(design = tech.design, version = rep(1, nrow(best)),
                            best = best, worst = worst,
                            alternative.names = names, algorithm = "HB",
                            hb.iterations = 10, hb.chains = 1, weights = wgt,
                            hb.warnings = FALSE),
                 "Weights are not able to be applied for Hierarchical Bayes.")
})

test_that("HB rank-ordered logit", {
    result <- FitMaxDiff(design = tech.design, version = rep(1, nrow(best)),
                         best = best, worst = worst,
                         alternative.names = names, algorithm = "HB",
                         hb.iterations = 10, hb.chains = 1,
                         is.tricked = FALSE, hb.warnings = FALSE)
    expect_error(print(result), NA)
    # expect_equal(result$in.sample.accuracy, 0.922185430463576)
})

test_that("HB 2 classes", {
    result <- FitMaxDiff(design = tech.design, version = rep(1, nrow(best)),
                         best = best, worst = worst,
                         alternative.names = names, algorithm = "HB",
                         hb.iterations = 10, hb.chains = 1,
                         n.classes = 2, hb.warnings = FALSE)
    expect_error(print(result), NA)
    # expect_equal(result$in.sample.accuracy, 0.931567328918322)
})

test_that("HB diagonal", {
    result <- FitMaxDiff(design = tech.design, version = rep(1, nrow(best)),
                         best = best, worst = worst,
                         alternative.names = names, algorithm = "HB",
                         hb.iterations = 10, hb.chains = 1,
                         normal.covariance = "Diagonal", hb.warnings = FALSE)
    expect_error(print(result), NA)
    # expect_equal(result$in.sample.accuracy, 0.931567328918322)
})

test_that("HB diagonal 2 classes", {
    result <- FitMaxDiff(design = tech.design, version = rep(1, nrow(best)),
                         best = best, worst = worst,
                         alternative.names = names, algorithm = "HB",
                         hb.iterations = 10, hb.chains = 1,
                         normal.covariance = "Diagonal", n.classes = 2,
                         hb.warnings = FALSE)
    expect_error(print(result), NA)
    # expect_equal(result$in.sample.accuracy, 0.929911699779249)
})

test_that("HB spherical", {
    result <- FitMaxDiff(design = tech.design, version = rep(1, nrow(best)),
                         best = best, worst = worst,
                         alternative.names = names, algorithm = "HB",
                         hb.iterations = 10, hb.chains = 1,
                         normal.covariance = "Spherical", hb.warnings = FALSE)
    expect_error(print(result), NA)
    # expect_equal(result$in.sample.accuracy, 0.934878587196468)
})
