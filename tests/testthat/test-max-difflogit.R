context("MaxDiff logit")

tech.data <- suppressWarnings(foreign::read.spss("http://wiki.q-researchsoftware.com/images/f/f1/Technology_2017.sav", to.data.frame = TRUE))
tech.design <- read.csv("http://wiki.q-researchsoftware.com/images/7/78/Technology_MaxDiff_Design.csv")
best <- tech.data[, c("Q5a_left", "Q5b_left", "Q5c_left", "Q5d_left",
                      "Q5e_left", "Q5f_left")]
worst <- tech.data[, c("Q5a_right", "Q5b_right", "Q5c_right", "Q5d_right",
                       "Q5e_right", "Q5f_right")]
names <- c("Apple", "Microsoft", "IBM", "Google", "Intel", "Samsung", "Sony",
           "Dell", "Yahoo", "Nokia")

test_that("Estimating logit parameters", {
    # Aggregate
    result <- FitMaxDiff(design = tech.design, version = rep(1, nrow(best)),
                         best = best, worst = worst, alternative.names = names,
                         is.tricked = FALSE)
    q.solution <- c(0.429913893438757, 0.339324277160547, -0.592482136137115,
                    0.8011237515597, -0.229733198090468, 0.459222489083277,
                    0.346698245594735, -0.440425502323207, -0.676080608209652,
                    -0.437561212076572)
    expect_equal(as.vector(result$coef), q.solution, tolerance = 0.0001)
    q.probs <- c(0.1358834315826, 0.1241144627641, 0.0488817685569,
                 0.1969619087625, 0.07025650843667, 0.1399245639353,
                 0.125033634106, 0.05690964591415, 0.04496150786246,
                 0.05707256807926)
    expect_equal(as.vector(result$class.preference.shares), q.probs,
                 tolerance = 0.0001)
    expect_error(print(result), NA)
    # Subset
    sub <- unclass(tech.data$Q2) <= 3
    result = FitMaxDiff(design = tech.design, version = rep(1, nrow(best)),
                        best = best, worst = worst, subset = sub,
                        alternative.names = names, is.tricked = FALSE)
    q.solution <- c(0.639695850685144, 0.548058173617228, -0.77639815378736,
                    1.18103309873161, -0.203097154503305, 0.0653366440742617,
                    0.157628024173681, -0.217544928475791, -0.890759970929614,
                    -0.503951583585859)
    expect_equal(as.vector(result$coef), q.solution, tolerance = 0.0001)
    # Subset and weight
    wgt <- tech.data$RESPNUM
    result <- FitMaxDiff(design = tech.design, version = rep(1, nrow(best)),
                         best = best, worst = worst, subset = sub,
                         weights = wgt, alternative.names = names,
                         is.tricked = FALSE)
    q.solution <- c(0.609489480351596, 0.473675078979805, -0.700733198610589,
                    1.03864583211041, -0.152305341107266, 0.076248739711427,
                    0.143219761355078, -0.202161480977952, -0.746942503957314,
                    -0.539136367855194)
    expect_equal(as.vector(result$coef), q.solution, tolerance = 0.00001)
    expect_equal(result$log.likelihood, -804.218, tolerance = 0.00001)

    expect_error(print(result), NA)

    # Cross validation
    result <- FitMaxDiff(design = tech.design, version = rep(1, nrow(best)),
                         best = best, worst = worst, alternative.names = names,
                         tasks.left.out = 2, is.tricked = TRUE)
    expect_error(print(result), NA)
    result <- FitMaxDiff(design = tech.design, version = rep(1, nrow(best)),
                         best = best, worst = worst, subset = sub,
                         weights = wgt, alternative.names = names,
                         tasks.left.out = 2, is.tricked = TRUE)
    expect_error(print(result), NA)
})

test_that("Latent class", {
    result <- FitMaxDiff(design = tech.design, version = rep(1, nrow(best)),
                         best = best, worst = worst, alternative.names = names,
                         n.classes = 2, is.tricked = FALSE)
    q.solution <- structure(c(2.55663493043737, 0.254412767640765,
                              -0.879102328586252, 0.840897783568564,
                              -0.476085243132221, 0.124556705874788,
                              0.135910368623237, -0.723184398493219,
                              -0.949129959810483, -0.884910626122545,
                              -1.01893136552894, 0.368312981702728,
                              -0.588559048171224, 0.831916954150721,
                              -0.0852423154677537, 0.957958780145827,
                              0.614930522695327, -0.315523724341152,
                              -0.689858256388732, -0.0750045287968016),
                            .Dim = c(10L, 2L))
    expect_equal(unname(result$coef), q.solution, tolerance = 0.01)
    expect_equal(result$log.likelihood, -4482.242, tolerance = 0.00001)
    q.probs <- structure(c(0.6118, 0.061, 0.0196, 0.1097, 0.0294, 0.0535,
                           0.0542, 0.0229, 0.0183, 0.0195, 0.0297, 0.1186,
                           0.0455, 0.1884, 0.0753,0.2138, 0.1516, 0.0598,
                           0.0412, 0.076), .Dim = c(10L, 2L))
    expect_equal(unname(result$class.preference.shares[, 1:2]), q.probs,
                 tolerance = 0.001)
    expect_error(print(result), NA)

    # Subset and weight
    sub <- unclass(tech.data$Q2) <= 3
    wgt <- tech.data$RESPNUM
    result <- FitMaxDiff(design = tech.design, version = rep(1, nrow(best)),
                         best = best, worst = worst, alternative.names = names,
                         n.classes = 2, weight = wgt, subset = sub,
                         is.tricked = FALSE)

    q.solution <- structure(c(2.2431358350092, 0.908145231913672,
                              -1.46834513563836, 1.72638477986375,
                              -0.482695445878375, 0.816295752009797,
                              -0.195367337130644, -0.925926041343033,
                              -1.52071774771516, -1.10090989109084,
                              -0.438796295404717, 0.0843811626415983,
                              -0.142267069468349, 0.609674028570444,
                              0.0032184699077446, -0.522916831318683,
                              0.474076349475965, 0.432650975133136,
                              -0.229604374528004, -0.270416415009135),
                            .Dim = c(10L, 2L))
    # This is quite lax right now. Fit will improve with the addition of a
    # numerical optimization stage.
    expect_equal(unname(result$coef), q.solution, tolerance = 0.01)
    expect_equal(result$log.likelihood, -730.724, tolerance = 0.01)

    expect_error(print(result), NA)

    expect_error(print(FitMaxDiff(design = tech.design,
                                  version = rep(1, nrow(best)),
                                  best = best,
                                  worst = worst,
                                  alternative.names = names,
                                  n.classes = 2,
                                  output = "Classes"),
                                  is.tricked = TRUE), NA)

    # Cross validation
    result <- FitMaxDiff(design = tech.design, version = rep(1, nrow(best)),
                         best = best, worst = worst, alternative.names = names,
                         n.classes = 2, tasks.left.out = 3, is.tricked = TRUE)
    expect_error(print(result), NA)
    result <- FitMaxDiff(design = tech.design, version = rep(1, nrow(best)),
                         best = best, worst = worst, alternative.names = names,
                         n.classes = 2, weight = wgt, subset = sub,
                         tasks.left.out = 3, is.tricked = TRUE)
    expect_error(print(result), NA)
})

test_that("Checking some of the inputs", {
    # No Version column in design
    expect_error(FitMaxDiff(design = tech.design[, -1],
                            version = rep(1, nrow(best)),
                            best = best, worst = worst,
                            alternative.names = names,
                            is.tricked = TRUE), NA)
    # No Task column in design
    expect_error(FitMaxDiff(design = tech.design[, -2],
                            version = rep(1, nrow(best)),
                            best = best, worst = worst,
                            alternative.names = names,
                            is.tricked = TRUE), NA)
    # Neither version nor Task column in design
    expect_error(FitMaxDiff(design = tech.design[, -1:-2],
                            version = rep(1, nrow(best)),
                            best = best, worst = worst,
                            alternative.names = names,
                            is.tricked = TRUE), NA)
    # No version
    expect_error(FitMaxDiff(design = tech.design[, -1:-2],
                            best = best, worst = worst,
                            alternative.names = names,
                            is.tricked = TRUE), NA)
    # inconsistent version information
    expect_error(FitMaxDiff(design = tech.design[, -1:-2],
                            version = 1:nrow(best),
                            best = best, worst = worst,
                            alternative.names = names,
                            is.tricked = TRUE))
    des <- tech.design
    des$Version <- 3
    expect_error(FitMaxDiff(design = des, best = best, worst = worst,
                            alternative.names = names, is.tricked = TRUE))

    # Incorrect names as a string
    alt.names <- "A,B,C,D,E,F,G,h,I,J"
    expect_error(suppressWarnings(FitMaxDiff(design = tech.design,
                                             best = best, worst = worst,
                                             alternative.names = alt.names,
                                             is.tricked = TRUE)))
    # Correct names as a string
    nms = paste(names, collapse = ", ")
    expect_error(FitMaxDiff(design = tech.design, best = best, worst = worst,
                            alternative.names = nms, is.tricked = TRUE), NA)
})

test_that("Varying coefficients", {
    # Aggregate
    ll.aggregate <- FitMaxDiff(design = tech.design,
                               version = rep(1, nrow(best)),
                               best = best, worst = worst,
                               alternative.names = names,
                               characteristics = NULL,
                               n.classes = 1,
                               is.tricked = TRUE)$log.likelihood
    # Gender varying coefficient
    expect_error(FitMaxDiff(design = tech.design,
                            version = rep(1, nrow(best)),
                            best = best, worst = worst,
                            alternative.names = names,
                            characteristics = data.frame(tech.data$Q1),
                            n.classes = 1, lc = FALSE, is.tricked = TRUE))
    # Apple varying coefficient
    ll.apple <- FitMaxDiff(design = tech.design, version = rep(1, nrow(best)),
                           best = best, worst = worst,
                           alternative.names = names,
                           characteristics = data.frame(tech.data$Q3_01),
                           n.classes = 1, lc = FALSE,
                           is.tricked = TRUE)$log.likelihood
    expect_true(ll.apple > ll.aggregate)
    # Apple varying coefficient with boosting
    d <- data.frame(tech.data$Q3_01)
    ll.apple.boosting.5.classes <- FitMaxDiff(design = tech.design,
                                              version = rep(1, nrow(best)),
                                              best = best, worst = worst,
                                              alternative.names = names,
                                              characteristics = d,
                                              n.classes = 5,
                                              is.tricked = TRUE)$log.likelihood
    ll.5.classes <- FitMaxDiff(design = tech.design,
                               version = rep(1, nrow(best)),
                               best = best, worst = worst,
                               alternative.names = names,
                               n.classes = 5,
                               is.tricked = TRUE)$log.likelihood
    expect_true(ll.apple.boosting.5.classes > ll.5.classes)

    d <- data.frame(tech.data$Q3_01)
    ll.apple.boosting.1.class <- FitMaxDiff(design = tech.design,
                                            version = rep(1, nrow(best)),
                                            best = best, worst = worst,
                                            alternative.names = names,
                                            characteristics = d,
                                            n.classes = 1,
                                            is.tricked = TRUE)$log.likelihood
    expect_true(ll.apple.boosting.1.class > ll.apple)
    expect_true(ll.apple.boosting.1.class > ll.aggregate)

    # Cross validation
    d <- data.frame(tech.data$Q3_01)
    apple.boosting.2.class <- FitMaxDiff(design = tech.design,
                                         version = rep(1, nrow(best)),
                                         best = best, worst = worst,
                                         alternative.names = names,
                                         characteristics = d,
                                         n.classes = 2, tasks.left.out = 1,
                                         is.tricked = TRUE)
    expect_error(print(apple.boosting.2.class), NA)

    #Subset
    sub <- unclass(tech.data$Q2) <= 3
    result <- FitMaxDiff(design = tech.design, version = rep(1, nrow(best)),
                         best = best, worst = worst, alternative.names = names,
                         characteristics = data.frame(tech.data$Q3_01),
                         n.classes = 2, tasks.left.out = 4, subset = sub,
                         is.tricked = TRUE)
    expect_error(print(result), NA)
})


test_that("Saving variables", {
    sub <- c(FALSE, FALSE, rep(TRUE, 100), rep(FALSE, 200))
    # Posterior probabilities.
    lc.3 <- FitMaxDiff(design = tech.design, version = rep(1, nrow(best)),
                       best = best, worst = worst, alternative.names = names,
                       n.classes = 3, subset = sub, is.tricked = TRUE)
    pp <- lc.3$posterior.probabilities[sub, ]
    expect_equal(ncol(pp), 3)
    expect_equal(sd(apply(pp, 1, sum)), 0)
    # Individual-level paramters
    pars <- RespondentParameters(lc.3)
    expect_equal(sum(!is.na(pars[, 1])), sum(sub))
    # Segment memberships
    m <- table(Memberships(lc.3))
    expect_equal(length(m), 3)
    expect_equal(sum(m), sum(sub))
})

test_that("Experimental designs with versions", {
    # President.
    dat <- suppressWarnings(foreign::read.spss("http://wiki.q-researchsoftware.com/images/6/61/President.sav", to.data.frame = TRUE))
    des <- read.csv("http://wiki.q-researchsoftware.com/images/9/9d/PresidentialDesign.csv")
    best <- dat[, c("MDmost_1", "MDmost_2", "MDmost_3", "MDmost_4", "MDmost_5",
                    "MDmost_6", "MDmost_7", "MDmost_8", "MDmost_9",
                    "MDmost_10")]
    worst <- dat[, c("MDleast_1", "MDleast_2", "MDleast_3" , "MDleast_4",
                     "MDleast_5", "MDleast_6",  "MDleast_7", "MDleast_8",
                     "MDleast_9", "MDleast_10")]
    names <- c("Decent/ethical", "Plain-speaking", "Healthy",
               "Successful in business", "Good in a crisis",
               "Experienced in government", "Concerned for minorities",
               "Understands economics", "Concerned about global warming",
               "Concerned about poverty", "Has served in the military",
               "Multilingual", "Entertaining", "Male",
               "From a traditional American background", "Christian")

    # We expect an error here because the alternative
    # "Concerned for minorities" should be
    # "Concerned for the welfare of minorities" instead.
    expect_error(suppressWarnings(FitMaxDiff(design = des, dat$MDversion,
                                             best = best, worst = worst,
                                             alternative.names = names,
                                             is.tricked = TRUE)))

    names <- c("Decent/ethical", "Plain-speaking", "Healthy",
               "Successful in business", "Good in a crisis",
               "Experienced in government",
               "Concerned for the welfare of minorities",
               "Understands economics", "Concerned about global warming",
               "Concerned about poverty", "Has served in the military",
               "Multilingual", "Entertaining", "Male",
               "From a traditional American background", "Christian")
    expect_error(FitMaxDiff(design = des, dat$MDversion, best = best,
                            worst = worst, alternative.names = names,
                            is.tricked = TRUE), NA)
    expect_error(FitMaxDiff(design = des, dat$MDversion, best = best,
                            worst = worst, alternative.names = names,
                            n.classes = 3, is.tricked = TRUE), NA)
    # Example from Q wiki
    # dat <- suppressWarnings(foreign::read.spss("http://wiki.q-researchsoftware.com/images/6/66/MaxDiffSetupExample.sav", to.data.frame = TRUE))
    # des <- read.csv("http://wiki.q-researchsoftware.com/images/2/24/ExampleMaxDiffDesign.csv")
    # best <- dat[, paste0("Q11task", 1:13,"most")]
    # worst <- dat[, paste0("Q11task", 1:13,"least")]
    # names <- paste("Brand", 1:13)
    # expect_error(FitMaxDiff(design = des, dat$MDversion, best = best, worst = worst, alternative.names = names))

})
