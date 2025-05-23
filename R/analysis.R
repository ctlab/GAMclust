#' Prepare gene expression matrix
#'
#' @param E Expression matrix with rownames as gene symbols.
#' @param gene.id.type Gene ID type.
#' @param keep.top.genes Which top of the most expressed genes to keep for the further analysis.
#' @param use.PCA Whether to reduce matrix dimensionality by PCA or not.
#' @param repeats Here you may collapse biological replicas by providing vector with repeated sample names
#' @param network.annotation Metabolic network annotation.
#' @return Expression matrix prepared for the analysis.
#' @import data.table
#' @export
prepareData <- function(
    E,
    gene.id.type = NULL,
    keep.top.genes = 12000,
    use.PCA = TRUE,
    use.PCA.n = 50,
    repeats = seq_len(ncol(E)),
    network.annotation){
  
  if(any(duplicated(repeats))){
    colnames(E) <- repeats
    E <- t(apply(E, 1, function(x) tapply(x, colnames(E), mean)))
  }
  
  new2old <- rownames(E)
  
  if(is.null(gene.id.type) || gene.id.type == network.annotation$baseId){
    flog.info("No gene annotation was performed.", name = "stats.logger")
  } else {
    
    if(gene.id.type %in% names(network.annotation$mapFrom)){
      
      rownames.dubl <- network.annotation$mapFrom[[gene.id.type]][rownames(E)]
      rownames(E) <- rownames.dubl[!duplicated(rownames.dubl[[gene.id.type]]), ]$gene
      
    } else {
      valid.types.string <- paste(c(names(network.annotation$mapFrom), "Entrez"), collapse = ", ")
      message.string <- sprintf("Invalid `gene.id.type`: %s. Must be one of: %s.", gene.id.type, valid.types.string)
      flog.error(message.string, name = "stats.logger")
      stop(message.string)
    }
  }
  
  names(new2old) <- rownames(E)

  E <- E[!is.na(rownames(E)), ]
  E <- E[order(rowMeans(E), decreasing = T), ]
  E <- E[!duplicated(rownames(E)), ]
  E <- E[head(order(rowMeans(E), decreasing = T), keep.top.genes), ]
  
  E <- t(base::scale(t(E), center = TRUE, scale = TRUE))
  E[is.na(E)] <- 0 # E <- E[rowSums(is.na(E)) == 0, ]
  
  if(use.PCA){
    if(use.PCA.n > ncol(E)){
      message.string <- sprintf("Please provide value of `use.PCA.n` smaller than `ncol(E)` = %s.", ncol(E))
      flog.error(message.string, name = "stats.logger")
      stop(message.string)
    }
    pcaRev <- irlba::prcomp_irlba(E, n = use.PCA.n, center = FALSE, scale. = FALSE, retx = TRUE)
    E.red <- pcaRev$x # E.red <- pcaRev$x %*% t(pcaRev$rotation)
    rownames(E.red) <- rownames(E)
    
    E <- E.red
    E <- t(base::scale(t(E), center = FALSE, scale = TRUE))
    E[is.na(E)] <- 0 # E <- E[rowSums(is.na(E)) == 0, ]
  } 
  
  attributes(E)$original.gene.names <- new2old
  
  E
}

#' Prepare network
#'
#' @param E Expression matrix after the `prepareData()` function.
#' @param network Metabolic network.
#' @param topology Vertices can be represented either as `metabolites`, either as `atoms`.
#' @param met.to.filter Metabolites that should not be used as connections in the module.
#' @param network.annotation Metabolic network annotation.
#' @param gene2reaction.extra For a combined network: supplementary file with genes that either come from proteome or are not linked to a specific enzyme.
#' @return Edges of the final network.
#' @import data.table
#' @export
prepareNetwork <- function(
    E,
    network,
    topology = c("metabolites", "atoms"),
    met.to.filter = data.table::fread(system.file("mets2mask.lst", package="GAMclust"))$ID,
    network.annotation,
    gene2reaction.extra = NULL){
  
  topology <- match.arg(topology)

  globalEdgeTable_pre <- as.data.frame(network$reaction2align)
  globalEdgeTable_pre <- merge(globalEdgeTable_pre, network$enzyme2reaction)
  globalEdgeTable_pre <- merge(globalEdgeTable_pre, network.annotation$gene2enzyme)
  
  if(!is.null(gene2reaction.extra)){
    globalEdgeTable_pre.extra <- merge(gene2reaction.extra, network$reaction2align)
    globalEdgeTable_pre.extra <- cbind(enzyme = "-.-.-.-", globalEdgeTable_pre.extra)
    globalEdgeTable_pre <- rbind(globalEdgeTable_pre, 
                                 globalEdgeTable_pre.extra) 
    }
  colnames(globalEdgeTable_pre)[which(colnames(globalEdgeTable_pre) == "atom.x")] <- "from"
  colnames(globalEdgeTable_pre)[which(colnames(globalEdgeTable_pre) == "atom.y")] <- "to"
  globalEdgeTable_pre$from.m <- network$atoms$metabolite[match(globalEdgeTable_pre$from, network$atoms$atom)]
  globalEdgeTable_pre$to.m <- network$atoms$metabolite[match(globalEdgeTable_pre$to, network$atoms$atom)]
  globalEdgeTable_pre <- globalEdgeTable_pre[which(!globalEdgeTable_pre$from.m %in% met.to.filter), ]
  globalEdgeTable_pre <- globalEdgeTable_pre[which(!globalEdgeTable_pre$to.m %in% met.to.filter), ]
  globalEdgeTable_pre <- globalEdgeTable_pre[which(globalEdgeTable_pre$gene %in% rownames(E)), ]
  
  if(topology == "atoms"){
    
    globalEdgeTable_pre <- globalEdgeTable_pre[, c("from", "to", "gene")]
    globalEdgeTable_pre <- globalEdgeTable_pre[!duplicated(globalEdgeTable_pre), ]
    
    flog.info("Global atom network contains %s edges.", dim(globalEdgeTable_pre)[1], 
              name = "stats.logger")
  }
  if(topology == "metabolites"){
    
    globalEdgeTable_pre <- globalEdgeTable_pre[, c("from.m", "to.m", "gene")]
    globalEdgeTable_pre <- globalEdgeTable_pre[!duplicated(globalEdgeTable_pre), ]
    colnames(globalEdgeTable_pre)[which(colnames(globalEdgeTable_pre) == "from.m")] <- "from"
    colnames(globalEdgeTable_pre)[which(colnames(globalEdgeTable_pre) == "to.m")] <- "to"
    
    flog.info("Global metabolite network contains %s edges.", dim(globalEdgeTable_pre)[1], 
              name = "stats.logger")
  }

  if (dim(globalEdgeTable_pre)[1] == 0) {
    message.string <- "No metabolic genes from the analysed dataset mapped to the metabolic network.\n
      In this case GAM-clustering will not work. Please try another subset of genes if it is possible."
    flog.error(message.string, name = "stats.logger")
    stop(message.string)
  }
  
  globalEdgeTable_pre_graph <- igraph::graph_from_data_frame(globalEdgeTable_pre, directed=FALSE)
  globalEdgeTable_pre_graph_cc <- igraph::decompose.graph(globalEdgeTable_pre_graph)
  globalGraph <- globalEdgeTable_pre_graph_cc[[which.max(sapply(globalEdgeTable_pre_graph_cc, igraph::vcount))]]
  # multi-edges, loops
  
  flog.info("Largest connected component of this global network contains %s nodes and %s edges.", 
            igraph::vcount(globalGraph), igraph::ecount(globalGraph), 
            name = "stats.logger")

  x.1p <- paste(globalEdgeTable_pre$from, globalEdgeTable_pre$to)
  x.2p <- with(igraph::as_data_frame(globalGraph), paste(c(from, to), c(to, from)))
  globalEdgeTable <- globalEdgeTable_pre[x.1p %in% x.2p, ]
  
  globalEdgeTable
}

#' Defining initial patterns
#'
#' @param E.prep Expression matrix after the `prepareData()` function.
#' @param network.prep Network edge table driven from `prepareNetwork()` function.
#' @param initial.number.of.clusters The number of clusters for the initial approximation of modules.
#' @param network.annotation Metabolic network annotation.
#' @return Initial patterns.
#' @export
preClustering <- function(E.prep,
                          network.prep,
                          initial.number.of.clusters = 32,
                          network.annotation,
                          use.ICA = FALSE
                          ){
  
  E.prep <- E.prep[rownames(E.prep) %in% network.prep$gene, , drop = F]
  flog.info("%d metabolic genes from the analysed dataset mapped to this component.",
            dim(E.prep)[1],
            name = "stats.logger")

  ### gene.cor <- cor(t(E.prep), use="pairwise.complete.obs")
  # gene.cor <- (E.prep %*% t(E.prep)) / max(rowSums(E.prep**2)) # max(rowSums(E.prep**2)) = x, while x+1 samples
  # gene.cor.dist <- as.dist(1 - gene.cor)
  # gene.pam <- cluster::pam(gene.cor.dist, k=initial.number.of.clusters)
  # cur.centers <- E.prep[gene.pam$medoids,]
  # OR
  gene.kmeans <- kmeans(E.prep, centers=initial.number.of.clusters)
  cur.centers <- gene.kmeans$centers
  
  if(use.ICA == TRUE){
    if(all(grepl("PC", colnames(E.prep)))) {
      ica_result <- fastICA::fastICA(t(E.prep), n.comp = initial.number.of.clusters) 
      cur.centers <- t(ica_result$S) } else {
        message.string <- "To perform ICA, set `use.PCA = TRUE` in `prepareData()` function."
        flog.error(message.string, name = "stats.logger")
        stop(message.string)
      }
  }
  
  cur.centers
}

#' GAM-clustering analysis
#'
#' @param E.prep Expression matrix after the `prepareData()` function.
#' @param network.prep Network edge table driven from `prepareNetwork()` function.
#' @param cur.centers Initial patterns produced by `preClustering()` function.
#' @param start.base The parameter which influences modules sizes.
#' @param base.dec The value by which `base` parameter should be reduced if some module's size is bigger that `max.module.size`.
#' @param max.module.size Maximal number of unique genes in the final module.
#' @param cor.threshold Threshold for correlation between module patterns.
#' @param p.adj.val.threshold Padj threshold of geseca score for final modules.
#' @param batch.solver Solver for SGMWCS problem.
#' @param work.dir Working directory where results should be saved.
#' @param show.intermediate.clustering Whether to show or not heatmap of intermediate clusters.
#' @param verbose Verbose running.
#' @param collect.stats Whether to save or not running statistics.
#' @return results$modules -- Metabolic modules.
#' @return results$nets -- Scored networks.
#' @return results$patterns.pos -- Modules' patterns (genes with positive score only considered).
#' @return results$patterns.all -- Modules' patterns (all genes considered).
#' @return results$iter.stats -- Statistics from iterations.
#' @export
gamClustering <- function(E.prep,
                          network.prep,
                          cur.centers,
                          
                          start.base = 0.5,
                          base.dec = 0.05,
                          max.module.size = 50,
                          
                          cor.threshold = 0.8,
                          p.adj.val.threshold = 0.001,
                          
                          batch.solver = seq_batch_solver(solver),
                          work.dir,
                          
                          show.intermediate.clustering = TRUE,
                          verbose = TRUE,
                          collect.stats = TRUE
                          ){
  
  flog.info("GAM-CLUSTERING starts here.", name = "stats.logger")
  
  iteration <- 1
  base <- start.base
  iter.stats <- list()
  
  while (T) {
    
    k <- 1
    revs <- list()
    
    while (T) {
      
      flog.info("[*] Iteration %s", iteration, name = "stats.logger")
      
      # 0. PREPARE ENVIRONMENT

      gK1 <- nrow(cur.centers)
      
      rev <- new.env()
      rev$modules <- list()
      rev$centers.pos <- matrix(nrow=gK1,
                                ncol=ncol(E.prep),
                                dimnames = list(
                                  paste0("c.pos", seq_len(gK1)),
                                  colnames(E.prep)))
      rev$centers.all <- matrix(nrow=gK1,
                                ncol=ncol(E.prep),
                                dimnames = list(
                                  paste0("c.all", seq_len(gK1)),
                                  colnames(E.prep)))
      
      # 1. CALCULATE CORRELATIONS -> DISTANCES -> SCORES

      # the projection of genes onto the centroids (measures of similarity between each gene and each centroid, cosine similarity between the two vectors):
      m <- cur.centers %*% t(E.prep) # 32 x samples * samples x genes = 32 x genes
      # ensure that the similarities between genes & centroids in m are not biased by differences in the magnitudes of ...
      # ...gene expression values (normalizes the columns of m):
      m <- m / sqrt(max(rowSums(E.prep**2))) # scales the rows of m by their Euclidean lengths -> m is [-1, 1]
      # ...centroid values (normalizes the rows of m):
      m <- sweep(x = m, MARGIN = 1, FUN = '/', STATS = sqrt(rowSums(cur.centers**2)))

      dist.to.centers <- 1-m
      dist.to.centers[dist.to.centers < 1e-10] <- 0

      idxs <- seq_len(gK1)
      
      posScores_keeping_var <- c()
      
      # TODO: replace base to `correlation.threshold`
      nets <- lapply(idxs, function(i) {

        if(nrow(dist.to.centers) > 1) { 
          minOther <- pmin(apply(dist.to.centers[-i, , drop=F], 2, min), base) } else { 
            minOther <- base }
        
        score <- log2(minOther) - log2(dist.to.centers[i, ])
        score[score == Inf] <- 0
        score <- pmax(score, -1000)
        posScores_keeping_var <<- c(posScores_keeping_var, length(which(score>0)))

        EdgeTable <- data.table::as.data.table(data.table::copy(network.prep))
        EdgeTable[, score := score[gene]]
        EdgeTable[from > to, c("from", "to") := list(to, from)]
        EdgeTable <- EdgeTable[order(score, decreasing = T)]
        EdgeTable <- unique(EdgeTable, by=c("from", "to"))
        # we still keep loops here

        scored_graph <- igraph::graph_from_data_frame(EdgeTable, directed = F)
        igraph::V(scored_graph)$score <- 0
        scored_graph
      })

      nets_attr <- lapply(nets, mwcsr::normalize_sgmwcs_instance,
                          edges.weight.column = "score",
                          nodes.weight.column = "score",
                          edges.group.by = "gene",
                          nodes.group.by = NULL,
                          group.only.positive = TRUE)
      
      # 2. SOLVE SGMWCS TO GET MODULES
      
      # cat("Calling: batch.solver(nets)\n")
      ms <- batch.solver(nets_attr)
      # cat("Done: batch.solver(nets)\n")
      ms_mods <- lapply(ms, `[[`, "graph")
      
      # 2.a. COLLECT CORRESPONDING LOGS
      
      m.size.unique <- unlist(lapply(ms_mods, function(x) ulength(igraph::edge_attr(x)$gene)))

      if (collect.stats) {
        iter.stats_add <- data.frame(
          genes.n = dim(E.prep)[1],
          genes.pos.scored = posScores_keeping_var,
          base = base,
          m.size = unlist(lapply(ms_mods, igraph::gsize)),
          m.size.unique = m.size.unique,
          m.pos = unlist(lapply(ms_mods, function(x) sum(igraph::edge_attr(x)$score > 0))),
          m.non.neg = unlist(lapply(ms_mods, function(x) sum(igraph::edge_attr(x)$score >= 0)))
        )
        iter.stats[[iteration]] <- iter.stats_add
      }
      if (verbose) {
        flog.info(">> base was equal to: %s", base, 
                  name = "stats.logger")
        flog.info(">> number of modules was equal to: %s", length(ms_mods), 
                  name = "stats.logger")
        flog.info(">> sizes of modules (unique genes) were in range: %s-%s", 
                  min(m.size.unique), max(m.size.unique), 
                  name = "stats.logger")
      }

      # 2.b. RECORD MODULES AND CENTERS
      
      rev$modules <- ms_mods
      
      for (i in idxs) {

        module <- ms_mods[[i]]

        center.pos <- if (ulength(igraph::E(module)[score > 0]$gene) >= 3) {
          getCenter(E.prep, unique(igraph::E(module)[score > 0]$gene))
        } else {
          cur.centers[i, ]
        }
        center.all <- if (ulength(igraph::E(module)$gene) >= 3) {
          getCenter(E.prep, unique(igraph::E(module)$gene))
        } else {
          cur.centers[i, ]
        }
        rev$centers.pos[i, ] <- center.pos
        rev$centers.all[i, ] <- center.all
      }

      if (show.intermediate.clustering) {
        heatmapTable <- rbind(cur.centers, rev$centers.pos)[rbind(
          seq_len(gK1),
          seq_len(gK1) + gK1), ]
        pheatmap::pheatmap(
          normalize.rows(heatmapTable),
          cluster_rows=F, cluster_cols=F,
          show_rownames=T, show_colnames=F)
      }
      
      revs[[k]] <- rev
      
      # 3. DID MODULES CONVERGE, i.e. CAN WE LEAVE THE SECOND LOOP
      
      # previous iterations, in which there was the same number of modules:
      revsToCheck <- revs[sapply(revs[seq_len(k-1)], function(rev) nrow(rev$centers.pos)) 
                          == nrow(rev$centers.pos)]
      
      diff <- max(abs(rev$centers.pos - cur.centers))
      
      if (length(revsToCheck) > 0) {
        diff <- min(sapply(revsToCheck,
                           function(prevRev) max(abs(rev$centers.pos - prevRev$centers.pos))))
      }
      
      # 4. UPDATE PARAMETERS

      cur.centers <- rev$centers.pos
      
      iteration <- iteration + 1
      
      if (verbose) {flog.info(">> max diff: %s", round(diff, 2), name = "stats.logger")}
      
      if (diff < 0.01) {break}

      k <- k + 1

    } # -------------------------------------------------------------------------------------- SECOND LOOP

    # 5. IF MODULES CONVERGED, WE CHECK THEM FOR PRESENCE OF 
    
    # (i) TOO BIG ONES:
    
    biggest.one <- max(sapply(ms_mods, function(m) ulength(igraph::E(m)$gene))) 
    
    if (biggest.one > max.module.size) {
      base <- base - base.dec
    }
    
    # (ii) CORRELATED ONES: 
    
    centers.cors <- cor(t(cur.centers))
    diag(centers.cors) <- 0
    correlation.max <- apply(centers.cors, 1, max, na.rm=T)
    
    if (any(correlation.max > cor.threshold)) {
      
      flog.info(">> max cor exceeded %s: %s", cor.threshold, round(max(correlation.max), 2),
                name = "stats.logger")
      max.cor.mod1 <- which.max(correlation.max) 
      max.cor.mod2 <- which.max(centers.cors[max.cor.mod1, ])
      cur.centers <- updCenters(cur.centers = cur.centers, 
                                m1 = max.cor.mod1, m2 = max.cor.mod2, 
                                E.prep = E.prep, ms_mods = ms_mods)
    } else {
      
      # (iii) or UNINFORMATIVE ONES:
      
      gesecaRes <- doGeseca(E.prep = E.prep,
                            network.prep = network.prep,
                            network.annotation = network.annotation,
                            modules = rev$modules,
                            scale = FALSE,
                            center = FALSE,
                            verbose = verbose)
      
      good <- gesecaRes$pathway[which(gesecaRes$padj < p.adj.val.threshold)]
      bad <- rownames(cur.centers)[!rownames(cur.centers) %in% good]
      
      if (length(bad) == 0 & biggest.one > max.module.size) {next} 
      if (length(bad) == 0 & biggest.one <= max.module.size) {break} 
      if (length(bad) != 0 & nrow(cur.centers) > 1) {
        
        max.cor.mod1 <- as.integer(gsub("c.pos", "", bad[which.max(apply(centers.cors, 1, max, na.rm=T)[bad])]))
        max.cor.mod2 <- which.max(centers.cors[max.cor.mod1, ])
        cur.centers <- updCenters(cur.centers = cur.centers, 
                                  m1 = max.cor.mod1, m2 = max.cor.mod2, 
                                  E.prep = E.prep, ms_mods = ms_mods)
      } else {
        saveStats(work.dir, rev, gesecaRes, iter.stats)
        message.string <- "[Attention!] No modules found.\n
          Try to tune method's parameters.\n
          Check '/stats' folder for the statistics of the run."
        flog.warn(message.string, name = "stats.logger")
        warning(message.string, call. = FALSE)
        return()
      }
    }

    # keep expressions devoted to sizes of modules:
    # m.sizes <- sapply(modules, function(m) ulength(igraph::E(m)$gene))
    # modules <- modules[m.sizes >= min.module.size] # add as param

  } # ---------------------------------------------------------------------------------------- FIRST LOOP
  
  # 9. FINAL ADJUSTMENTS OF MODULES 
  
  # (i) compactise
  modules_pre <- lapply(rev$modules, function(x) {
    igraph::graph.attributes(x)$signals[which(names(igraph::graph.attributes(x)$signals) %in%
                                                igraph::vertex_attr(x)$signal)] <- -0.001
    x
  })
  modules_set <- batch.solver(modules_pre)
  modules <- lapply(modules_set, `[[`, "graph")

  # (ii) recalc geseca & sort modules
  gesecaRes <- doGeseca(E.prep = E.prep,
                        network.prep = network.prep,
                        network.annotation = network.annotation,
                        modules = modules,
                        scale = FALSE,
                        center = FALSE,
                        verbose = verbose)
  modules <- modules[as.numeric(gsub("c.pos", "", gesecaRes$pathway))]
  gesecaRes$pathway <- paste0("m", 1:nrow(gesecaRes)) 

  saveStats(work.dir, rev, gesecaRes, iter.stats)
  
  flog.info("GAM-CLUSTERING ends here.", name = "stats.logger")

  return(list(
    modules = modules,
    nets = nets_attr,
    patterns.pos = rev$centers.pos,
    patterns.all = rev$centers.all,
    iter.stats = iter.stats
  ))
}
