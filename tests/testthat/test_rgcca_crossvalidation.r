set.seed(1)

data("Russett")
blocks <- list(
    agriculture = Russett[, seq(3)],
    industry = Russett[, 4:5],
    politic = Russett[, 6:11])

test_structure_cv <- function(res, scores, nrow = 47){
    expect_equal(length(res), 3)
    expect_is(res, "cv")
    expect_is(res$rgcca, "rgcca")
    pred <- res$preds
    expect_is(pred, "list")
    expect_is(pred[[1]], "matrix")
    expect_true(all(sapply(pred, NCOL) == 2))
    expect_true(all(sapply(pred, NROW) == nrow))
    expect_identical(round(res$scores, 4), scores)
}

test_that("rgcca_cv_default", {
        rgcca_out <- rgcca(blocks, response = 2)
        test_structure_cv(
            rgcca_crossvalidation(rgcca_out, n_cores = 1),
            0.0923)
        test_structure_cv(
            rgcca_crossvalidation(
                rgcca_out,
                blocks = blocks,
                n_cores = 1), 
            scores = 0.0923)
    }
)

test_that("rgcca_cv_with_args", {
    rgcca_out <- rgcca(blocks, response = 1)
    test_structure_cv(
        rgcca_crossvalidation(
            rgcca_out,
            validation = "kfold",
            k = 5,
            n_cores = 1),
        0.1122)
    # test_structure_cv(
    #     rgcca_crossvalidation(
    #         rgcca_out, 
    #         validation = "test", 
    #         n_cores = 1), 
    #     0.1057) # TODO : warnings
    }
)

test_that("rgcca_cv_withNA", {

    RussettWithNA <- Russett
    RussettWithNA[1:2,1:3] <- NA
    RussettWithNA[3,4:5] <- NA
    RussettWithNA[3,1] <- NA
    blocksNA <- list(
        agriculture = RussettWithNA[, seq(3)], 
        industry = RussettWithNA[, 4:5],
        politic = RussettWithNA[, 6:11])

    # cross validation
    rgcca_out <- rgcca(blocksNA, response = 3)
    # avec la method complete -> ne fonctionne pas
    test_structure_cv(
        rgcca_crossvalidation(rgcca_out, n_cores = 1),
        0.0938, 
        nrow = 44)
    }
)