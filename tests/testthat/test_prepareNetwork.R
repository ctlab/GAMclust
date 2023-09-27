options(timeout = 600)
data("bulkRNAdata_ex_good")
data("bulkRNAdata_ex_bad")
data("scRNAdata_small_ex")

met.to.filter <- data.table::fread(system.file("mets2mask.lst", package="GAMclust"))$ID
# KEGG:
network.kegg <- readRDS(url("http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/network.kegg.rds"))
# Rhea
network.rhea <- readRDS(url("http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/network.rhea.rds"))

network.annotation.mm <- readRDS(url("http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/org.Mm.eg.gatom.anno.rds"))
network.annotation.hs <- readRDS(url("http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/org.Hs.eg.gatom.anno.rds"))

E.prep <- prepareData(E = Biobase::exprs(data.good),
                      gene.id.type = "Symbol",
                      keep.top.genes = 12000,
                      use.PCA = TRUE,
                      use.PCA.n = ncol(Biobase::exprs(data.good)) - 1,
                      repeats = seq_len(ncol(Biobase::exprs(data.good))),
                      network.annotation = network.annotation.mm)

E.prep.bad <- prepareData(E = Biobase::exprs(data.bad),
                          gene.id.type = "Symbol",
                          keep.top.genes = 12000,
                          use.PCA = TRUE,
                          use.PCA.n = ncol(Biobase::exprs(data.bad)) - 1,
                          repeats = seq_len(ncol(Biobase::exprs(data.bad))),
                          network.annotation = network.annotation.mm)

E.prep.hs <- prepareData(E = Biobase::exprs(sc.data),
                         gene.id.type = "Entrez",
                         keep.top.genes = 12000,
                         use.PCA = TRUE,
                         use.PCA.n = ncol(Biobase::exprs(sc.data)) - 1,
                         repeats = seq_len(ncol(Biobase::exprs(sc.data))),
                         network.annotation = network.annotation.hs)


test_that("prepareNetwork works for Mus musculus (KEGG), metabolites", {
  
  r <- evaluate_promise(prepareNetwork(E = E.prep,
                                       topology = "metabolites",
                                       met.to.filter = met.to.filter,
                                       network = network.kegg,
                                       network.annotation = network.annotation.mm))
  
  expect_equal(colnames(r$result), expected = c("from", "to", "gene"))
  expect_false(any(rownames(r$result) %in% met.to.filter))
  expect_true(any(grepl("> Global metabolite network contains",
                        r$messages)))
})


test_that("prepareNetwork works for Mus musculus (Rhea), metabolites", {
  
  r <- evaluate_promise(prepareNetwork(E = E.prep,
                                       topology = "metabolites",
                                       met.to.filter = met.to.filter,
                                       network = network.rhea,
                                       network.annotation = network.annotation.mm))
  
  expect_equal(colnames(r$result), expected = c("from", "to", "gene"))
  expect_false(any(rownames(r$result) %in% met.to.filter))
  expect_true(any(grepl("> Largest connected component of this global network contains",
                    r$messages)))
})


test_that("prepareNetwork works for Mus musculus (KEGG), atoms", {
  
  r <- evaluate_promise(prepareNetwork(E = E.prep,
                                       topology = "atoms",
                                       met.to.filter = met.to.filter,
                                       network = network.kegg,
                                       network.annotation = network.annotation.mm))
  
  expect_equal(colnames(r$result), expected = c("from", "to", "gene"))
  expect_false(any(rownames(r$result) %in% met.to.filter))
  expect_true(any(grepl("> Global atom network contains",
                        r$messages)))
})


test_that("prepareNetwork works for Mus musculus (Rhea), atoms", {
  
  r <- evaluate_promise(prepareNetwork(E = E.prep,
                                       topology = "atoms",
                                       met.to.filter = met.to.filter,
                                       network = network.rhea,
                                       network.annotation = network.annotation.mm))
  
  expect_equal(colnames(r$result), expected = c("from", "to", "gene"))
  expect_false(any(rownames(r$result) %in% met.to.filter))
  expect_true(any(grepl("> Largest connected component of this global network contains",
                        r$messages)))
})


test_that("prepareNetwork works for Homo sapiens (KEGG), metabolites", {
  
  r <- evaluate_promise(prepareNetwork(E = E.prep.hs,
                                       topology = "metabolites",
                                       met.to.filter = met.to.filter,
                                       network = network.kegg,
                                       network.annotation = network.annotation.hs))
  
  expect_equal(colnames(r$result), expected = c("from", "to", "gene"))
  expect_true(nrow(network.prep) > 0)
  expect_false(any(rownames(r$result) %in% met.to.filter))
  expect_true(any(grepl("> Global metabolite network contains",
                        r$messages)))
})


test_that("prepareNetwork stops with no metabolic genes", {
  
  expect_error(prepareNetwork(E = E.prep.bad,
                                topology = "metabolites",
                                met.to.filter = met.to.filter,
                                network = network.kegg,
                                network.annotation = network.annotation.mm),
               "No metabolic genes from the analysed dataset mapped to the metabolic network.\n
      In this case GAM-clustering will not work. Please try another subset of genes if it is possible.")
})
