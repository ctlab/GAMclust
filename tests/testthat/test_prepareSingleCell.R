data("scRNAseurat_ex")


test_that("prepareSingleCell works without clustering", {

  seurat_object <- SeuratObject::pbmc_small

  sc.data <- prepareSingleCell(seurat_object,
                               cluster.data = FALSE,
                               cluster.identity = "RNA_snn_res.1",
                               assay = "RNA",
                               slot = "scale.data",
                               organism = "hsa")
  expect_length(sc.data, 1)
  expect_true(names(sc.data) == "es")
  expect_true(all(
    dim(sc.data[["es"]]@assayData[["exprs"]]) == c(18, 3)))

})


test_that("prepareSingleCell works with clustering", {

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

  unlink(clust.dir, recursive = T)

})

