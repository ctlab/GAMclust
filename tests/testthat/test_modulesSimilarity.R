options(bitmapType = "cairo")
options(timeout = 600)
data("GAMresults_ex")
data("GAMresults_atoms_ex")
data("bulkRNAdata_ex_good_prep")

network.annotation.mm <- readRDS(url("http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/org.Mm.eg.gatom.anno.rds"))


test_that("modulesSimilarity works", {
  
  work.dir <- file.path(tempdir(), "testGAMclustering")
  dir.create(work.dir)

  # clustering 1
  work.dir1 <- file.path(tempdir(), "testGAMclustering1")
  dir.create(work.dir1)
  getGeneTables(modules = results$modules,
                nets = results$nets,
                patterns = results$patterns.pos,
                gene.exprs = E.prep,
                network.annotation = network.annotation.mm,
                work.dir = work.dir1)
  dir.create(sprintf("%s/stats", work.dir1), showWarnings=FALSE, recursive=TRUE)
  write.tsv(results$patterns.pos, file = sprintf("%s/stats/patterns_pos.tsv", work.dir1))
  write.tsv(results$patterns.all, file = sprintf("%s/stats/patterns_all.tsv", work.dir1))

  # clustering 2
  work.dir2 <- file.path(tempdir(), "testGAMclustering2")
  dir.create(work.dir2)
  getGeneTables(modules = results.atom$modules,
                nets = results.atom$nets,
                patterns = results.atom$patterns.pos,
                gene.exprs = E.prep,
                network.annotation = network.annotation.mm,
                work.dir = work.dir2)
  dir.create(sprintf("%s/stats", work.dir2), showWarnings=FALSE, recursive=TRUE)
  write.tsv(results.atom$patterns.pos, file = sprintf("%s/stats/patterns_pos.tsv", work.dir2))
  write.tsv(results.atom$patterns.all, file = sprintf("%s/stats/patterns_all.tsv", work.dir2))


  r <- evaluate_promise(modulesSimilarity(dir1 = work.dir1,
                                          dir2 = work.dir2,
                                          name1 = "new",
                                          name2 = "old",
                                          same.data = TRUE,
                                          use.genes.with.pos.score = TRUE,
                                          work.dir = work.dir,
                                          file.name = "comparison.png"))

  expect_true(r$result == file.path(work.dir, "comparison.png"))
  expect_true(length(r$warnings) == 0)
  
  unlink(work.dir, recursive = T)
  unlink(work.dir1, recursive = T)
  unlink(work.dir2, recursive = T)

})
