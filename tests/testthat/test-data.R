context("data")

tech.data = suppressWarnings(foreign::read.spss("http://wiki.q-researchsoftware.com/images/f/f1/Technology_2017.sav", to.data.frame = TRUE))
tech.design = read.csv("http://wiki.q-researchsoftware.com/images/7/78/Technology_MaxDiff_Design.csv")
best = tech.data[, c("Q5a_left", "Q5b_left", "Q5c_left", "Q5d_left", "Q5e_left", "Q5f_left")]
worst = tech.data[, c("Q5a_right", "Q5b_right", "Q5c_right", "Q5d_right", "Q5e_right", "Q5f_right")]
names <- c("Apple", "Microsoft", "IBM", "Google", "Intel", "Samsung", "Sony", "Dell", "Yahoo", "Nokia")

pres.data <- suppressWarnings(foreign::read.spss("http://wiki.q-researchsoftware.com/images/6/61/President.sav", to.data.frame = TRUE))
pres.design <- read.csv("http://wiki.q-researchsoftware.com/images/9/9d/PresidentialDesign.csv")
pres.best <- pres.data[, c("MDmost_1", "MDmost_2", "MDmost_3", "MDmost_4", "MDmost_5",
                "MDmost_6", "MDmost_7", "MDmost_8", "MDmost_9",
                "MDmost_10")]
pres.worst <- pres.data[, c("MDleast_1", "MDleast_2", "MDleast_3" , "MDleast_4",
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


    expect_equal(extractAlternativeNames(pres.design, pres.best, pres.worst),
                 pres.names)

    reduced.levels <- levels(pres.best[[3]])[-13] # remove unused level
    pres.best[[3]] <- factor(as.numeric(pres.best[[3]]))
    levels(pres.best[[3]]) <- reduced.levels
    expect_error(extractAlternativeNames(pres.design, pres.best, pres.worst))
})

test_that("Construct design", {
    design.alternatives <- structure(list(
        alt1 = structure(c(1L, 1L, 2L, 3L, 2L, 1L),
                         class = "factor",
                         .Label = c("1", "2", "3", "4", "5", "6", "7", "8",
                                    "9", "10"),
                         questiontype = "PickOneMulti",
                         name = "Alt.1", label = "1", question = "Alt"),
        alt2 = structure(c(2L, 7L, 5L, 4L, 4L, 3L),
                         class = "factor",
                         .Label = c("1", "2", "3", "4", "5", "6", "7", "8",
                                    "9", "10"),
                         questiontype = "PickOneMulti",
                         name = "Alt.2", label = "2", question = "Alt"),
        alt3 = structure(c(3L, 8L, 6L, 5L, 6L, 5L),
                         class = "factor",
                         .Label = c("1", "2", "3", "4", "5", "6", "7", "8",
                                    "9", "10"),
                         questiontype = "PickOneMulti",
                         name = "Alt.3", label = "3", question = "Alt"),
        alt4 = structure(c(4L, 9L, 7L, 7L, 8L, 6L),
                         class = "factor",
                         .Label = c("1", "2", "3", "4", "5", "6", "7", "8",
                                    "9", "10"),
                         questiontype = "PickOneMulti",
                         name = "Alt.4", label = "4", question = "Alt"),
        alt5 = structure(c(10L, 10L, 10L, 8L, 9L, 9L),
                         class = "factor",
                         .Label = c("1", "2", "3", "4", "5", "6", "7", "8",
                                    "9", "10"),
                         questiontype = "PickOneMulti",
                         name = "Alt.5", label = "5", question = "Alt")),
        .Names = c("alt1", "alt2", "alt3", "alt4", "alt5" ),
        row.names = c(NA, -6L), class = "data.frame")

    design.version <- structure(c(1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L,
                                  2L),
                                class = "factor", .Label = c("1", "2"),
                                questiontype = "PickOne", name = "Version",
                                label = "Version", question = "Version")

    expect_equal(constructDesign(design.alternatives, NULL),
                 structure(list(Version = c(1, 1, 1, 1, 1, 1),
                                Task = 1:6,
                                Alt.1 = c(1, 1, 2, 3, 2, 1),
                                Alt.2 = c(2, 7, 5, 4, 4, 3),
                                Alt.3 = c(3, 8, 6, 5, 6, 5),
                                Alt.4 = c(4, 9, 7, 7, 8, 6),
                                Alt.5 = c(10, 10, 10, 8, 9, 9)),
                           .Names = c("Version", "Task", "Alt.1",
                                      "Alt.2", "Alt.3", "Alt.4", "Alt.5"),
                           row.names = c(NA, -6L), class = "data.frame"))

    expect_equal(constructDesign(rbind(design.alternatives,
                                       design.alternatives), design.version),
        structure(list(Version = c(1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2),
        Task = c(1L, 2L, 3L, 4L, 5L, 6L, 1L, 2L, 3L, 4L, 5L, 6L),
        Alt.1 = c(1, 1, 2, 3, 2, 1, 1, 1, 2, 3, 2, 1),
        Alt.2 = c(2, 7, 5, 4, 4, 3, 2, 7, 5, 4, 4, 3),
        Alt.3 = c(3, 8, 6, 5, 6, 5, 3, 8, 6, 5, 6, 5),
        Alt.4 = c(4, 9, 7, 7, 8, 6, 4, 9, 7, 7, 8, 6),
        Alt.5 = c(10, 10, 10, 8, 9, 9, 10, 10, 10, 8, 9, 9)),
        .Names = c("Version", "Task", "Alt.1", "Alt.2",
                   "Alt.3", "Alt.4", "Alt.5"),
        row.names = c(NA, -12L), class = "data.frame"))

    invalid.design.version <- structure(c(1L, 1L, 1L, 1L, 1L),
                                class = "factor", .Label = c("1", "2"),
                                questiontype = "PickOne", name = "Version",
                                label = "Version", question = "Version")

    expect_error(constructDesign(design.alternatives,
                                 invalid.design.version),
                 paste0("The design version is invalid as its length does ",
                        "not match the design alternatives."))
})
