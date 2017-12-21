context("data")

tech.data = suppressWarnings(foreign::read.spss("http://wiki.q-researchsoftware.com/images/f/f1/Technology_2017.sav", to.data.frame = TRUE))
tech.design = read.csv("http://wiki.q-researchsoftware.com/images/7/78/Technology_MaxDiff_Design.csv")
best = tech.data[, c("Q5a_left", "Q5b_left", "Q5c_left", "Q5d_left", "Q5e_left", "Q5f_left")]
worst = tech.data[, c("Q5a_right", "Q5b_right", "Q5c_right", "Q5d_right", "Q5e_right", "Q5f_right")]
names <- c("Apple", "Microsoft", "IBM", "Google", "Intel", "Samsung", "Sony", "Dell", "Yahoo", "Nokia")

pres.data <- suppressWarnings(foreign::read.spss("http://wiki.q-researchsoftware.com/images/6/61/President.sav", to.data.frame = TRUE))
pres.design <- read.csv("http://wiki.q-researchsoftware.com/images/9/9d/PresidentialDesign.csv")
pres.best <- dat[, c("MDmost_1", "MDmost_2", "MDmost_3", "MDmost_4", "MDmost_5",
                "MDmost_6", "MDmost_7", "MDmost_8", "MDmost_9",
                "MDmost_10")]
pres.worst <- dat[, c("MDleast_1", "MDleast_2", "MDleast_3" , "MDleast_4",
                 "MDleast_5", "MDleast_6",  "MDleast_7", "MDleast_8",
                 "MDleast_9", "MDleast_10")]
pres.names <- c("Decent/ethical", "Plain-speaking", "Healthy",
           "Successful in business", "Good in a crisis",
           "Experienced in government",
           "Concerned for the welfare of minorities",
           "Understands economics", "Concerned about global warming",
           "Concerned about poverty", "Has served in the military",
           "Multilingual", "Entertaining", "Male",
           "From a traditional American background", "Christian")

test_that("Reading data works", {
    # expect_error(IntegrateDesignAndData(design = tech.design, version = rep(1, nrow(best)), best = best, worst = worst), NA)

    list.design = MaxDiffDesign(number.alternatives = 10, number.questions = 6, alternatives.per.question = 5, n.repeats = 1)
    expect_error(cleanAndCheckData(design = list.design, best = best, worst = worst, alternative.names = names))
    list.design$design <- tech.design
    expect_error(cleanAndCheckData(design = list.design, best = best, worst = worst, alternative.names = names), NA)

    binary.design = list.design$binary.design
    binary.design[binary.design == 1] <- 0
    for (r in 1:nrow(binary.design))
        for (c in 1:5)
            binary.design[r, c] = tech.design[r, c + 2]
    expect_error(cleanAndCheckData(design = binary.design, best = best, worst = worst, alternative.names = names), NA)
})

test_that("Extract alternative names", {
    extracted <- extractAlternativeNames(tech.design, best, worst)
    tech.names <- c("Apple", "Microsoft", "IBM", "Google", "Intel", "Samsung ",
               "Sony", "Dell", "Yahoo", "Nokia")
    expect_equal(extracted, tech.names)


    expect_equal(extractAlternativeNames(pres.design, pres.best, pres.worst), pres.names)

    reduced.levels <- levels(pres.best[[3]])[-13] # remove unused level
    pres.best[[3]] <- factor(as.numeric(pres.best[[3]]))
    levels(pres.best[[3]]) <- reduced.levels
    expect_error(extractAlternativeNames(pres.design, pres.best, pres.worst))
})

