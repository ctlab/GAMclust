data("scRNAseurat_ex")

met.to.filter <- data.table::fread(system.file("mets2mask_MashaE_now.lst", package="GAMclust"))$ID
# KEGG:
network <- readRDS(url("http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/network.kegg.rds"))
network.annotation <- readRDS(url("http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/org.Hs.eg.gatom.anno.rds"))
metabolites.annotation <- readRDS(url("http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/met.kegg.db.rds"))


test_that("Pipeline works for single cell data", {

  work.dir <- file.path(tempdir(), "gamTest_pipelineSc")
  dir.create(work.dir, showWarnings = F)

  cplex.dir <- "/opt/ibm/ILOG/CPLEX_Studio1271"
  solver <- mwcsr::virgo_solver(cplex_dir = cplex.dir)

  clust.dir <- file.path(tempdir(), "chR")
  sc.data <- prepareSingleCell(object = seurat_object,
                               cluster.data = TRUE,
                               return.seurat = TRUE,
                               reduction = "pca",
                               ndims = 10,
                               resolutions = c( seq(1, 2, 0.5), 2.5, 3, 5),
                               chooseR.cluster.n.times = 3,
                               clustering.folder = clust.dir,
                               assay = "RNA",
                               slot = "scale.data",
                               organism = "hsa")
  expect_length(sc.data, 2)
  expect_true(all(names(sc.data) == c("es", "seurat")))
  expect_true("ChooseR_kmeans" %in% colnames(sc.data$seurat@meta.data))
  expect_true(class(sc.data$es) == "ExpressionSet")

  data <- sc.data$es

  repeats <- rownames(Biobase::pData(data))

  metObjs <- prepareObjects(data = data,
                            network = network,
                            network.annotation = network.annotation,
                            topology = "metabolites",
                            met.to.filter = met.to.filter)
  expect_false(any(rownames(metObjs$gene.exprs) %in% met.to.filter))
  expect_length(metObjs, n = 3)

  exprsMetab.centers <- preClustering(gene.exprs = metObjs$exprsMetab,
                                      repeats = repeats,
                                      initial.number.of.clusters = 8,
                                      show.initial.clustering = TRUE)
  expect_gt(sum(dim(exprsMetab.centers)), expected = 2)

  results <- gamClustering(gene.exprs = metObjs$exprsMetab,
                           repeats = repeats,
                           cur.centers = exprsMetab.centers,
                           edge.table = metObjs$globalEdgeTable,
                           batch.solver = seq_batch_solver(solver),
                           show.intermediate.clustering = TRUE,
                           verbose = TRUE,
                           work.dir = work.dir,
                           collect.stats = TRUE)
  expect_true(all(sapply(results$modules, class) == "igraph"))

  r <- evaluate_promise(getGraphs(modules = results$modules,
                                  network.annotation = network.annotation,
                                  metabolites.annotation = metabolites.annotation,
                                  seed.for.layout = 22,
                                  work.dir = work.dir))
  expect_true(all(grepl("Graphs for module \\d are built", r$messages)))

  r <- evaluate_promise(getHeatmaps(modules = results$modules,
                                    patterns = results$patterns.pos,
                                    gene.exprs = metObjs$exprsMetab,
                                    order = 1:ncol(metObjs$exprsMetab),
                                    data.annotation = 1:ncol(results$patterns.pos),
                                    network.annotation = network.annotation,
                                    width = 25, height = 10,
                                    work.dir = work.dir))
  expect_true(is.null(r$result))

  r <- evaluate_promise(getGeneTables(modules = results$modules,
                                      nets = results$nets,
                                      patterns = results$patterns.pos,
                                      gene.exprs = metObjs$exprsMetab,
                                      network.annotation = network.annotation,
                                      work.dir = work.dir))
  expect_length(r$warnings, 0)

  r <- evaluate_promise(getModulesAnnotation(network.annotation = network.annotation,
                                             nets = results$nets,
                                             work.dir = work.dir))
  expect_true(all(grepl("Pathway annotation for module \\d is produced", r$messages)))

  r <- evaluate_promise(heatmapModulesAnnotation(work.dir = work.dir,
                                                 threshold = 0,
                                                 file_name = "Modules_heatmap.png"))
  expect_true(class(r$result) == "pheatmap")

  unlink(work.dir, recursive = T)
  unlink(clust.dir, recursive = T)

})
