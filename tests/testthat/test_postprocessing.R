options(timeout = 600)
data("GAMresults_ex")
data("GAMresults_atoms_ex")
data("bulkRNAdata_ex_good_prep")

network.annotation.mm <- readRDS(url("http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/org.Mm.eg.gatom.anno.rds"))
metabolites.annotation.kegg <- readRDS(url("http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/met.kegg.db.rds"))

work.dir <- file.path(tempdir(), "testGAMclustering")
dir.create(work.dir, showWarnings = F)


test_that("getGraphs works", {

  r <- evaluate_promise(getGraphs(modules = results$modules,
                                  network.annotation = network.annotation.mm,
                                  metabolites.annotation = metabolites.annotation.kegg,
                                  seed.for.layout = 42,
                                  work.dir = work.dir))
  expect_true(is.null(r$result))
  expect_true(all(grepl("Graphs for module \\d are built", r$messages)))

})


test_that("getGraphs works with atom-based networks", {

  r <- evaluate_promise(getGraphs(modules = results.atom$modules,
                                  network.annotation = network.annotation.mm,
                                  metabolites.annotation = metabolites.annotation.kegg,
                                  seed.for.layout = 42,
                                  work.dir = work.dir))
  expect_true(is.null(r$result))
  expect_true(all(grepl("Graphs for module \\d are built", r$messages)))

})


test_that("getGeneTables works", {

  r <- evaluate_promise(getGeneTables(modules = results$modules,
                                      nets = results$nets,
                                      patterns = results$patterns.pos,
                                      gene.exprs = E.prep,
                                      network.annotation = network.annotation.mm,
                                      work.dir = work.dir))
  expect_length(r$result, length(results$modules))
  expect_length(r$warnings, 0)
  expect_true(all(grepl("Gene tables for module \\d are produced", r$messages)))

})


test_that("getAnnotationTables works", {

  r <- evaluate_promise(getAnnotationTables(network.annotation = network.annotation.mm,
                                            nets = results$nets,
                                            work.dir = work.dir))
  expect_true(is.null(r$result))
  expect_true(all(grepl("Pathway annotation for module \\d is produced", r$messages)))

})


test_that("getAnnotationHeatmap works", {

  r <- evaluate_promise(getAnnotationHeatmap(work.dir = work.dir))
  expect_true(class(r$result) == "pheatmap")
  expect_length(r$warnings, 0)
  expect_true(
    sum(grepl("Processing module \\d", r$messages)) == sum(grepl("Module size: \\d", r$messages))
    )

})


unlink(work.dir, recursive = T)
