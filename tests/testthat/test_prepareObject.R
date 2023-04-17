data("bulkRNAdata_ex_good")
data("bulkRNAdata_ex_bad")
data("scRNAdata_small_ex")


met.to.filter <- data.table::fread(system.file("mets2mask_MashaE_now.lst", package="GAMclust"))$ID
# KEGG:
network.kegg <- readRDS(url("http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/network.kegg.rds"))
# Rhea
network.rhea <- readRDS(url("http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/network.rhea.rds"))

network.annotation.mm <- readRDS(url("http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/org.Mm.eg.gatom.anno.rds"))
network.annotation.hs <- readRDS(url("http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/org.Hs.eg.gatom.anno.rds"))


test_that("prepareObjects works for KEGG, metabolites (Mus musculus)", {

  metObjs <- prepareObjects(data = data.good,
                            network = network.kegg,
                            network.annotation = network.annotation.mm,
                            topology = "metabolites",
                            met.to.filter = met.to.filter)

  expect_false(any(rownames(metObjs$gene.exprs) %in% met.to.filter))
  expect_length(metObjs, n = 3)
  expect_true(igraph::is.igraph(metObjs$globalGraph))
  expect_message(prepareObjects(data = data.good,
                                network = network.kegg,
                                network.annotation = network.annotation.mm,
                                topology = "metabolites",
                                met.to.filter = met.to.filter),
                 "> Global metabolite network contains")
})



test_that("prepareObjects works for KEGG, atoms (Mus musculus)", {

  metObjs <- prepareObjects(data = data.good,
                            network = network.kegg,
                            network.annotation = network.annotation.mm,
                            topology = "atoms",
                            met.to.filter = met.to.filter)

  expect_false(any(rownames(metObjs$gene.exprs) %in% met.to.filter))
  expect_length(metObjs, n = 3)
  expect_true(igraph::is.igraph(metObjs$globalGraph))
})



test_that("prepareObjects works for Rhea, metabolites (Mus musculus)", {

  metObjs <- prepareObjects(data = data.good,
                            network = network.rhea,
                            network.annotation = network.annotation.mm,
                            topology = "metabolites",
                            met.to.filter = met.to.filter)

  expect_false(any(rownames(metObjs$gene.exprs) %in% met.to.filter))
  expect_length(metObjs, n = 3)
  expect_true(igraph::is.igraph(metObjs$globalGraph))
})



test_that("prepareObjects works for Rhea, atoms (Mus musculus)", {

  metObjs <- prepareObjects(data = data.good,
                            network = network.rhea,
                            network.annotation = network.annotation.mm,
                            topology = "atoms",
                            met.to.filter = met.to.filter)

  expect_false(any(rownames(metObjs$gene.exprs) %in% met.to.filter))
  expect_length(metObjs, n = 3)
  expect_true(igraph::is.igraph(metObjs$globalGraph))
})



test_that("prepareObjects stops with no metabolic genes", {

  expect_error(prepareObjects(data = data.bad,
                              network = network.kegg,
                              network.annotation = network.annotation.mm,
                              topology = "metabolites",
                              met.to.filter = met.to.filter),
               "No metabolic genes from the analysed dataset mapped to the metabolic network")

})



test_that("prepareObjects works for sc, KEGG, metabolites (Homo sapiens)", {

  metObjs <- prepareObjects(data = sc.data,
                            network = network.kegg,
                            network.annotation = network.annotation.hs,
                            topology = "metabolites",
                            met.to.filter = met.to.filter)

  expect_false(any(rownames(metObjs$gene.exprs) %in% met.to.filter))
  expect_length(metObjs, n = 3)
  expect_true(igraph::is.igraph(metObjs$globalGraph))
})
