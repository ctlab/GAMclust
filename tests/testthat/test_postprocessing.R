data("GAMresults_ex")
data("GAMresults_atoms_ex")
data("metObjs")

network.annotation.mm <- readRDS(url("http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/org.Mm.eg.gatom.anno.rds"))
metabolites.annotation.kegg <- readRDS(url("http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/met.kegg.db.rds"))

work.dir <- file.path(tempdir(), "gamTest")
dir.create(work.dir, showWarnings = F)



test_that("getGraphs works", {

  r <- evaluate_promise(getGraphs(modules = results$modules,
                                  network.annotation = network.annotation.mm,
                                  metabolites.annotation = metabolites.annotation.kegg,
                                  seed.for.layout = 22,
                                  work.dir = work.dir))
  expect_true(is.null(r$result))
  expect_true(all(grepl("Graphs for module \\d are built", r$messages)))

})



test_that("getGraphs works with renameAtoms", {

  r <- evaluate_promise(getGraphs(modules = results.atom$modules,
                                  network.annotation = network.annotation.mm,
                                  metabolites.annotation = metabolites.annotation.kegg,
                                  seed.for.layout = 22,
                                  work.dir = work.dir))
  expect_true(is.null(r$result))
  expect_true(all(grepl("Graphs for module \\d are built", r$messages)))

})



test_that("getHeatmaps works", {

  r <- evaluate_promise(getHeatmaps(modules = results$modules,
                                    patterns = results$patterns.pos,
                                    gene.exprs = metObjs$exprsMetab,
                                    order = 1:ncol(metObjs$exprsMetab),
                                    data.annotation = 1:ncol(results$patterns.pos),
                                    network.annotation = network.annotation.mm,
                                    width = 25, height = 10,
                                    work.dir = work.dir))
  expect_true(is.null(r$result))
  expect_length(r$warnings, 0)
  expect_true(grepl("Heatmap for patterns is produced", r$messages)[1])
  expect_true(all(grepl("Heatmap for module \\d is produced", r$messages)[-c(1)]))

})



test_that("getGeneTables works", {

  r <- evaluate_promise(getGeneTables(modules = results$modules,
                                      nets = results$nets,
                                      patterns = results$patterns.pos,
                                      gene.exprs = metObjs$exprsMetab,
                                      network.annotation = network.annotation.mm,
                                      work.dir = work.dir))
  expect_true(is.null(r$result))
  expect_length(r$warnings, 0)
  expect_true(all(grepl("Gene tables for module \\d are produced", r$messages)))

})



test_that("getModulesAnnotation works", {

  r <- evaluate_promise(getModulesAnnotation(network.annotation = network.annotation.mm,
                                             nets = results$nets,
                                             work.dir = work.dir))
  expect_true(is.null(r$result))
  expect_length(r$warnings, 0)
  expect_true(all(grepl("Pathway annotation for module \\d is produced", r$messages)))

})



test_that("heatmapModulesAnnotation works", {

  r <- evaluate_promise(heatmapModulesAnnotation(work.dir = work.dir,
                                                 threshold = 0,
                                                 file_name = "Modules_heatmap.png"))
  expect_true(class(r$result) == "pheatmap")
  expect_length(r$warnings, 0)
  expect_true(
    sum(grepl("Processing module \\d", r$messages)) == sum(grepl("Module size: \\d", r$messages))
    )

})



unlink(work.dir, recursive = T)
