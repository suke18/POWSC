# simulation of expression data
library(POWSC)
data("es_mef_sce")
sce = es_mef_sce[, colData(es_mef_sce)$cellTypes == "fibro"]
set.seed(12)
rix = sample(1:nrow(sce), 500)
sce = sce[rix, ]


# test the input data. Alternatively, it could a class of "matrix".
test_that('Check the input sce', {
    expect_that(sce, is_a('SingleCellExperiment'))
})



# test the estimated parameter
test_that("Check the output of Est2Phase", {
    est_Paras = Est2Phase(sce)
    # gene-wise
    expect_equal(length(est_Paras$mu), nrow(sce))
    expect_equal(length(est_Paras$sd), nrow(sce))
    # cell-wise
    expect_equal(length(est_Paras$sf), ncol(sce))
    expect_equal(length(est_Paras$lambda), ncol(sce))

})

# test the main function runPOWSC
test_that('Check the output of runPOWSC', {
    est_Paras = Est2Phase(sce)
    sim_size = c(100, 200) # A numeric vector
    pow_rslt = runPOWSC(sim_size = sim_size, est_Paras = est_Paras,per_DE=0.05, DE_Method = "MAST", Cell_Type = "PW")
    pow_rslt1 = pow_rslt[[1]]
    pow_rslt2 = pow_rslt[[2]]

    expect_true(all(names(pow_rslt1[[2]]) == names(pow_rslt1[[2]])))
    expect_true(all(names(pow_rslt2[[2]]) == names(pow_rslt2[[2]])))

    expect_equal(length(pow_rslt1[[1]][[1]]), length(pow_rslt2[[1]][[1]]))
    expect_equal(length(pow_rslt1[[2]][[1]]), length(pow_rslt2[[2]][[1]]))

})


