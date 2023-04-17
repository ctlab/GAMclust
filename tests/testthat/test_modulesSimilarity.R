data("metObjs")
data("metObjs_atoms")

network.annotation <- readRDS(url("http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/org.Mm.eg.gatom.anno.rds"))

test_that("modulesSimilarity works", {

  cplex.dir <- "/opt/ibm/ILOG/CPLEX_Studio1271"
  solver <- mwcsr::virgo_solver(cplex_dir = cplex.dir)

  # clustering 1
  work.dir1 <- file.path(tempdir(), "testSimilarity1")
  dir.create(work.dir1)
  repeats <- colnames(metObjs$exprsMetab)

  exprsMetab.centers <- preClustering(gene.exprs = metObjs$exprsMetab,
                                      repeats = repeats,
                                      initial.number.of.clusters = 8,
                                      show.initial.clustering = TRUE)

  results <- gamClustering(gene.exprs = metObjs$exprsMetab,
                           repeats = repeats,
                           cur.centers = exprsMetab.centers,
                           edge.table = metObjs$globalEdgeTable,
                           batch.solver = seq_batch_solver(solver),
                           show.intermediate.clustering = TRUE,
                           verbose = TRUE,
                           collect.stats = TRUE,
                           work.dir = work.dir1)

  getGeneTables(modules = results$modules,
                nets = results$nets,
                patterns = results$patterns.pos,
                gene.exprs = metObjs$exprsMetab,
                network.annotation = network.annotation,
                work.dir = work.dir1)

  # clustering 1
  work.dir2 <- file.path(tempdir(), "testSimilarity2")
  dir.create(work.dir2)
  repeats <- colnames(metObjs.atoms$exprsMetab)

  exprsMetab.centers <- preClustering(gene.exprs = metObjs.atoms$exprsMetab,
                                      repeats = repeats,
                                      initial.number.of.clusters = 8,
                                      show.initial.clustering = TRUE)

  results <- gamClustering(gene.exprs = metObjs.atoms$exprsMetab,
                           repeats = repeats,
                           cur.centers = exprsMetab.centers,
                           edge.table = metObjs.atoms$globalEdgeTable,
                           batch.solver = seq_batch_solver(solver),
                           show.intermediate.clustering = TRUE,
                           verbose = TRUE,
                           collect.stats = TRUE,
                           work.dir = work.dir2)

  getGeneTables(modules = results$modules,
                nets = results$nets,
                patterns = results$patterns.pos,
                gene.exprs = metObjs.atoms$exprsMetab,
                network.annotation = network.annotation,
                work.dir = work.dir2)


  r <- evaluate_promise(modulesSimilarity(dir1 = work.dir1,
                                          dir2 = work.dir2,
                                          name1 = "analysisOfDataset1",
                                          name2 = "analysisOfDataset1_changedParameters",
                                          same.data = TRUE,
                                          use.genes.with.pos.score = TRUE,
                                          work.dir = work.dir1,
                                          file.name = "compareTo_changedParameters.png"))

  expect_true(r$result == file.path(work.dir1, "compareTo_changedParameters.png"))
  expect_true(length(r$warnings) == 0)

  unlink(work.dir1, recursive = T)
  unlink(work.dir2, recursive = T)

})
