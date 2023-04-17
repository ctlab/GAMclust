data("metObjs")


test_that("preClustering works", {

  repeats <- colnames(metObjs$exprsMetab)

  exprsMetab.centers <- preClustering(gene.exprs = metObjs$exprsMetab,
                                      repeats = repeats,
                                      initial.number.of.clusters = 8,
                                      show.initial.clustering = TRUE)

  expect_true(is.matrix(exprsMetab.centers))
  expect_gt(sum(dim(exprsMetab.centers)), expected = 2)

})



test_that("gamClustering works", {

  cplex.dir <- "/opt/ibm/ILOG/CPLEX_Studio1271"
  solver <- mwcsr::virgo_solver(cplex_dir = cplex.dir)

  repeats <- colnames(metObjs$exprsMetab)

  exprsMetab.centers <- preClustering(gene.exprs = metObjs$exprsMetab,
                                      repeats = repeats,
                                      initial.number.of.clusters = 8,
                                      show.initial.clustering = TRUE)

  work.dir <- file.path(tempdir(), "testGAMclustering")
  dir.create(work.dir)

  results <- gamClustering(gene.exprs = metObjs$exprsMetab,
                           repeats = repeats,
                           cur.centers = exprsMetab.centers,
                           edge.table = metObjs$globalEdgeTable,
                           batch.solver = seq_batch_solver(solver),
                           show.intermediate.clustering = TRUE,
                           verbose = TRUE,
                           collect.stats = TRUE,
                           work.dir = work.dir)

   expect_length(results, 5)
  expect_true(all(
    names(results) %in% c("modules","nets",
                          "patterns.pos", "patterns.all",
                          "iter.stats")))
  expect_true(all(sapply(results$modules, class) == "igraph"))
  expect_false(any(is.na(results$patterns.pos)))

  unlink(work.dir, recursive = T)

})

