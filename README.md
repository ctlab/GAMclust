# GAM-clustering

GAMclust is a tool for metabolic GAM-clustering analysis of

-   bulk RNA-seq data ([tutorial](https://rpubs.com/anastasiiaNG/GAMclust_BULK)),
-   single cell RNA-seq data ([tutorial](https://rpubs.com/anastasiiaNG/GAMclust_SC)),
-   spatial RNA-seq data.

The GAM-clustering method reveals metabolic variability within dataset using a novel network-based computational approach. This method utilizes cellular transcriptional profiles as proxies. The metabolic network of reactions is presented as a graph. This graph has vertices corresponding to metabolites and the edges corresponding to the reactions with the expressed genes. In the graph the method tries to find a set of connected subgraphs. Genes of each connected subgraph have highly correlated expression pattern.

Interactive results of GAM-clustering analysis of ImmGen Open Source and Tabula Muris data can be explored [here]((http://artyomovlab.wustl.edu/immgen-met/)), details are in [Gainullina et al, 2023](https://www.cell.com/cell-reports/fulltext/S2211-1247(23)00057-8).
