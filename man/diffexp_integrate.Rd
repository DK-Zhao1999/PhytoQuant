\name{diffexp_integrate}
\alias{diffexp_integrate}
\title{Integrated Differential Expression Analysis and Gene-Set Construction}
\usage{
diffexp_integrate(
  vector,
  geneSets_list,
  data_dir = getwd(),
  importance_threshold = 10,
  rra_p_threshold = 0.00001,
  sorenson_threshold = 0.6,
  min_geneset_size = 5,
  max_na_ratio = 0.3,
  save_intermediate = FALSE,
  methods = c("DESeq2", "edgeR", "limma"),
  linkage_method = "complete",
  verbose = TRUE
)
}
\arguments{
  \item{vector}{Character vector of dataset names to process.}

  \item{geneSets_list}{Named list of prior knowledge gene sets. Each element
  must be a character vector of gene identifiers.}

  \item{data_dir}{Directory containing the expression-count and group-annotation
  RData files. Defaults to the current working directory.}

  \item{importance_threshold}{Numeric threshold on the absolute importance score
  used to select de novo up- and down-regulated genes.}

  \item{rra_p_threshold}{Numeric threshold applied to the robust-rank-aggregation
  score when selecting genes across datasets.}

  \item{sorenson_threshold}{Numeric Sorensen-Dice similarity threshold used to
  merge redundant gene sets.}

  \item{min_geneset_size}{Minimum number of genes required for a final gene set.}

  \item{max_na_ratio}{Maximum allowed proportion of missing importance values
  for a gene across datasets.}

  \item{save_intermediate}{Logical. If \code{TRUE}, saves the per-dataset merged
  differential-expression tables.}

  \item{methods}{Character vector of differential-expression methods to run.
  May include \code{"DESeq2"}, \code{"edgeR"}, and \code{"limma"}.}

  \item{linkage_method}{Hierarchical clustering linkage method used by the
  Sorensen-Dice merging step.}

  \item{verbose}{Logical. If \code{TRUE}, prints progress messages.}
}
\value{
A list with components:
\item{activation_gene_sets}{Named list of activation gene sets.}
\item{inhibition_gene_sets}{Named list of inhibition gene sets.}
\item{importance_all}{Data frame of importance scores for all genes.}
\item{importance_up_gene}{Robust-rank-aggregation results for up-regulated genes.}
\item{importance_dn_gene}{Robust-rank-aggregation results for down-regulated genes.}
}
\description{
Combines DESeq2, edgeR, and limma with robust rank aggregation,
knowledge-guided refinement, and Sorensen-Dice-based redundancy removal to
discover direction-resolved activation and inhibition gene sets from plant
transcriptome cohorts.
}
\examples{
\dontrun{
# See the package README for the complete simulated-data workflow.
library(PhytoQuant)

# Simulate and save two count matrices and group tables.
set.seed(1)
genes <- paste0("PHYTO", sprintf("\%03d", 1:100))
make_counts <- function(seed, genes) {
  set.seed(seed)
  n <- 12
  mu <- matrix(100, nrow = length(genes), ncol = n)
  mu[1:10, 1:6] <- 1000
  mu[11:20, 7:12] <- 1000
  matrix(rnbinom(length(genes) * n, size = 20, mu = mu), nrow = length(genes),
         ncol = n, dimnames = list(genes, paste0("S", 1:n)))
}
dir.create(tmp <- tempfile())
for (ds in c("PhytoA", "PhytoB")) {
  expr <- make_counts(if (ds == "PhytoA") 1 else 2, genes)
  group <- data.frame(Tag = colnames(expr),
                      group = rep(c("H", "L"), each = 6),
                      stringsAsFactors = FALSE)
  assign(ds, expr)
  assign(paste0(ds, "_G"), group)
  save(list = ds, file = file.path(tmp, paste0(ds, ".Rdata")))
  save(list = paste0(ds, "_G"), file = file.path(tmp, paste0(ds, "_G.Rdata")))
}

knowledge <- list(auxin = genes[1:10], ethylene = genes[11:20])
result <- diffexp_integrate(
  vector = c("PhytoA", "PhytoB"),
  geneSets_list = knowledge,
  data_dir = tmp,
  importance_threshold = 1,
  rra_p_threshold = 0.01,
  sorenson_threshold = 0.7,
  min_geneset_size = 3
)
str(result)
}
}
