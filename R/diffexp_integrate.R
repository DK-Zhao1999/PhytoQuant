#' Integrated differential expression analysis and gene-set construction
#'
#' Runs DESeq2, edgeR, and limma on one or more plant transcriptome datasets,
#' synthesizes the per-gene evidence into a signed importance score, integrates
#' the resulting de novo gene sets with prior knowledge gene sets using robust
#' rank aggregation, and returns direction-resolved activation and inhibition
#' gene sets after Sorensen-Dice-based redundancy removal.
#'
#' @param vector Character vector of dataset names to process. Each dataset must
#'   have a corresponding \code{<name>.Rdata} expression-count file and a
#'   \code{<name>_G.Rdata} group-annotation file.
#' @param geneSets_list Named list of known gene sets used as prior knowledge.
#'   Each element must be a character vector of gene identifiers.
#' @param data_dir Directory containing the RData files. Defaults to the current
#'   working directory.
#' @param importance_threshold Numeric threshold on the absolute importance
#'   score used to select de novo up- and down-regulated genes. Default is 10.
#' @param rra_p_threshold Numeric threshold applied to the robust-rank-aggregation
#'   score when selecting genes across datasets. Default is 0.00001.
#' @param sorenson_threshold Numeric Sorensen-Dice similarity threshold used to
#'   merge redundant gene sets. Default is 0.6.
#' @param min_geneset_size Minimum number of genes required for a final gene set.
#'   Default is 5.
#' @param max_na_ratio Maximum allowed proportion of missing importance values
#'   for a gene across datasets. Default is 0.3.
#' @param save_intermediate Logical. If \code{TRUE}, saves the per-dataset merged
#'   differential-expression tables. Default is \code{FALSE}.
#' @param methods Character vector of differential-expression methods to run.
#'   May include \code{"DESeq2"}, \code{"edgeR"}, and \code{"limma"}. Defaults to
#'   all three.
#' @param linkage_method Hierarchical clustering linkage method used by the
#'   Sorensen-Dice merging step. Default is \code{"complete"}.
#' @param verbose Logical. If \code{TRUE}, prints progress messages. Default is
#'   \code{TRUE}.
#'
#' @return A list containing:
#'   \item{activation_gene_sets}{A named list of activation (up-regulated) gene
#'     sets.}
#'   \item{inhibition_gene_sets}{A named list of inhibition (down-regulated)
#'     gene sets.}
#'   \item{importance_all}{A data frame of importance scores for all genes across
#'     all input datasets.}
#'   \item{importance_up_gene}{Robust-rank-aggregation results for up-regulated
#'     genes.}
#'   \item{importance_dn_gene}{Robust-rank-aggregation results for down-regulated
#'     genes.}
#'
#' @export
diffexp_integrate <- function(
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
) {

  # Move to the requested data directory unless it is already active.
  if (!is.null(data_dir) && data_dir != getwd()) {
    setwd(data_dir)
    if (verbose) cat("Working directory set to:", data_dir, "\n")
  }

  # Load the packages used by the differential-expression backends.
  required_packages <- c("DESeq2", "edgeR", "limma", "RobustRankAggreg")
  for (pkg in required_packages) {
    if (!require(pkg, character.only = TRUE)) {
      stop(paste("Package", pkg, "is required but not installed."))
    }
  }

  # Initialize the container that will hold one importance table per dataset.
  if (verbose) cat("Processing", length(vector), "datasets...\n")

  # Store the per-dataset importance tables before the cross-study merge.
  all_importance <- list()

  # Analyze each supplied dataset independently.
  for (i in seq_along(vector)) {
    dataset_name <- vector[i]
    if (verbose) cat("\n=== Processing dataset:", dataset_name, "===\n")

    # Locate and load the expression-count and group-annotation objects.
    exp_file <- paste0(dataset_name, ".Rdata")
    g_file <- paste0(dataset_name, "_G.Rdata")

    if (!file.exists(exp_file) || !file.exists(g_file)) {
      warning(paste("Data files for", dataset_name, "not found. Skipping."))
      next
    }

    load(exp_file, envir = environment())
    load(g_file, envir = environment())

    exp <- get(dataset_name, envir = environment())
    exp_G <- get(paste0(dataset_name, "_G"), envir = environment())

    # Build a count matrix and a compatible sample-level design table.
    countData <- as.matrix(exp)
    colData <- data.frame(
      row.names = exp_G$Tag,
      condition = factor(exp_G$group, levels = c("L", "H"))
    )
    countData <- countData[, exp_G$Tag, drop = FALSE]

    # Initialize the list that collects successful method-specific tables.
    res_list <- list()

    # 1. DESeq2: model raw counts with a negative-binomial GLM.
    if ("DESeq2" %in% methods) {
      if (verbose) cat("Running DESeq2...\n")
      tryCatch({
        dds <- DESeq2::DESeqDataSetFromMatrix(
          countData = countData,
          colData = colData,
          design = ~ condition
        )
        dds$condition <- relevel(dds$condition, ref = "L")
        dds <- DESeq2::DESeq(dds)
        res_DESeq2 <- DESeq2::results(dds, contrast = c("condition", "H", "L"))
        res_DESeq2 <- data.frame(
          gene = rownames(res_DESeq2),
          logFC.DESeq2 = res_DESeq2$log2FoldChange,
          p.value.DESeq2 = res_DESeq2$pvalue
        )
        res_DESeq2$logFC.DESeq2[is.na(res_DESeq2$logFC.DESeq2)] <- 0
        res_DESeq2$p.value.DESeq2[is.na(res_DESeq2$p.value.DESeq2)] <- 1
        res_list$DESeq2 <- res_DESeq2
      }, error = function(e) {
        warning(paste("DESeq2 failed for", dataset_name, ":", e$message))
      })
    }

    # 2. edgeR: estimate dispersions and run a quasi-likelihood test.
    if ("edgeR" %in% methods) {
      if (verbose) cat("Running edgeR...\n")
      tryCatch({
        dge <- edgeR::DGEList(counts = countData, group = colData$condition)
        dge <- edgeR::calcNormFactors(dge)
        design <- stats::model.matrix(~ condition, data = colData)
        dge <- edgeR::estimateDisp(dge, design)
        fit <- edgeR::glmQLFit(dge, design)
        qlf <- edgeR::glmQLFTest(fit, coef = 2)
        res_edgeR <- edgeR::topTags(qlf, n = Inf)$table
        res_edgeR <- data.frame(
          gene = rownames(res_edgeR),
          logFC.edgeR = res_edgeR$logFC,
          p.value.edgeR = res_edgeR$PValue
        )
        res_edgeR$logFC.edgeR[is.na(res_edgeR$logFC.edgeR)] <- 0
        res_edgeR$p.value.edgeR[is.na(res_edgeR$p.value.edgeR)] <- 1
        res_list$edgeR <- res_edgeR
      }, error = function(e) {
        warning(paste("edgeR failed for", dataset_name, ":", e$message))
      })
    }

    # 3. limma-voom: transform counts and fit a weighted linear model.
    if ("limma" %in% methods) {
      if (verbose) cat("Running limma...\n")
      tryCatch({
        design <- stats::model.matrix(~ condition, data = colData)
        v <- limma::voomWithQualityWeights(countData, design = design, plot = FALSE)
        fit <- limma::lmFit(v, design = design)
        fit <- limma::eBayes(fit)
        res_limma <- limma::topTable(fit, coef = 2, number = Inf, sort.by = "none")
        res_limma <- data.frame(
          gene = rownames(res_limma),
          logFC.limma = res_limma$logFC,
          p.value.limma = res_limma$P.Value
        )
        res_limma$logFC.limma[is.na(res_limma$logFC.limma)] <- 0
        res_limma$p.value.limma[is.na(res_limma$p.value.limma)] <- 1
        res_list$limma <- res_limma
      }, error = function(e) {
        warning(paste("limma failed for", dataset_name, ":", e$message))
      })
    }

    # Skip this dataset if none of the requested methods succeeded.
    if (length(res_list) == 0) {
      warning(paste("No method succeeded for", dataset_name, ". Skipping."))
      next
    }

    # Merge all successful method-specific tables by gene identifier.
    if (verbose) cat("Calculating importance scores...\n")
    res_all <- Reduce(function(x, y) merge(x, y, by = "gene"), res_list)

    # Locate the log-fold-change and p-value columns produced by the methods.
    logFC_cols <- grep("logFC", colnames(res_all), value = TRUE)
    p_cols <- grep("p.value", colnames(res_all), value = TRUE)

    # Combine the evidence into a signed importance score.
    res_all$logFC <- rowMeans(res_all[, logFC_cols, drop = FALSE], na.rm = TRUE)
    res_all$p.value <- exp(rowMeans(log(res_all[, p_cols, drop = FALSE]), na.rm = TRUE))
    res_all$importance <- res_all$logFC * (-log10(res_all$p.value))

    # Construct direction-resolved de novo gene sets.
    res_all_imp <- res_all[, c("gene", "importance")]
    res_all_imp <- res_all_imp[order(res_all_imp$importance), ]

    res_all_up <- res_all_imp[res_all_imp$importance > importance_threshold, ]
    res_all_dn <- res_all_imp[res_all_imp$importance < -importance_threshold, ]

    assign(paste0(dataset_name, "_up_denovo"), res_all_up$gene, envir = environment())
    assign(paste0(dataset_name, "_dn_denovo"), res_all_dn$gene, envir = environment())

    # Store the per-dataset importance table for later cross-study integration.
    colnames(res_all_imp) <- c("gene", dataset_name)
    all_importance[[dataset_name]] <- res_all_imp

    # Optionally save the merged table for this dataset.
    if (save_intermediate) {
      save(res_all, file = paste0(dataset_name, "_results.Rdata"))
    }
  }

  # Merge the importance tables from all datasets into one gene-by-study matrix.
  if (verbose) cat("\n=== Merging importance scores across datasets ===\n")

  importance_all <- NULL
  for (i in seq_along(all_importance)) {
    imp <- all_importance[[i]]
    if (i == 1) {
      importance_all <- imp
    } else {
      importance_all <- merge(importance_all, imp, by = "gene", all = TRUE)
    }
  }

  if (is.null(importance_all) || nrow(importance_all) == 0) {
    stop("No data available. All datasets failed to process.")
  }

  rownames(importance_all) <- importance_all[, 1]
  importance_all <- importance_all[, -1, drop = FALSE]

  # Remove genes whose importance values are missing in too many datasets.
  na_ratio <- rowSums(is.na(importance_all)) / ncol(importance_all)
  importance_clean <- importance_all[na_ratio <= max_na_ratio, ]

  if (verbose) cat("Keeping", nrow(importance_clean), "genes after NA filtering\n")

  # RRA integration for down-regulated genes: low importance ranks first.
  if (verbose) cat("Running RRA for down-regulated genes...\n")
  importance_dn <- lapply(importance_clean, function(x) {
    if (is.numeric(x)) {
      ranks <- rank(x, na.last = "keep", ties.method = "min")
      rownames(importance_clean)[order(ranks, na.last = NA)]
    }
  })
  importance_dn <- RobustRankAggreg::aggregateRanks(importance_dn)
  importance_dn <- importance_dn[, "Score", drop = FALSE]
  importance_dn_gene <- importance_dn[importance_dn$Score < rra_p_threshold, , drop = FALSE]

  # RRA integration for up-regulated genes: high importance ranks first.
  if (verbose) cat("Running RRA for up-regulated genes...\n")
  importance_up <- lapply(importance_clean, function(x) {
    if (is.numeric(x)) {
      ranks <- rank(x, na.last = "keep", ties.method = "min")
      rownames(importance_clean)[order(ranks, na.last = NA, decreasing = TRUE)]
    }
  })
  importance_up <- RobustRankAggreg::aggregateRanks(importance_up)
  importance_up <- importance_up[, "Score", drop = FALSE]
  importance_up_gene <- importance_up[importance_up$Score < rra_p_threshold, , drop = FALSE]

  # Integrate prior knowledge gene sets and de novo gene sets.
  if (verbose) cat("Integrating knowledge gene sets and de novo gene sets...\n")

  # Add the per-dataset de novo gene sets to the knowledge list.
  for (dataset_name in names(all_importance)) {
    gs_up <- get(paste0(dataset_name, "_up_denovo"), envir = environment())
    gs_dn <- get(paste0(dataset_name, "_dn_denovo"), envir = environment())
    geneSets_list[[paste0(dataset_name, "_up")]] <- gs_up
    geneSets_list[[paste0(dataset_name, "_dn")]] <- gs_dn
  }

  # Retain only genes that also passed cross-dataset robust rank aggregation.
  geneSets_up_filtered <- lapply(geneSets_list, function(gs) {
    intersect(gs, rownames(importance_up_gene))
  })
  geneSets_up_filtered <- geneSets_up_filtered[sapply(geneSets_up_filtered, length) > 0]
  names(geneSets_up_filtered) <- paste0("ags", seq_along(geneSets_up_filtered))

  # Retain only genes that also passed cross-dataset robust rank aggregation.
  geneSets_dn_filtered <- lapply(geneSets_list, function(gs) {
    intersect(gs, rownames(importance_dn_gene))
  })
  geneSets_dn_filtered <- geneSets_dn_filtered[sapply(geneSets_dn_filtered, length) > 0]
  names(geneSets_dn_filtered) <- paste0("igs", seq_along(geneSets_dn_filtered))

  # Calculate the Sorensen-Dice coefficient between two gene sets.
  sorensen_dice <- function(set1, set2) {
    intersection <- length(intersect(set1, set2))
    total <- length(set1) + length(set2)
    if (total == 0) return(NA)
    (2 * intersection) / total
  }

  # Merge gene sets whose pairwise similarity is above the requested threshold.
  aux.merge <- function(gsets, minsds, linkage_method) {
    if (length(gsets) <= 1) {
      return(gsets)
    }

    integration_matrix <- outer(seq_along(gsets), seq_along(gsets),
                                Vectorize(function(x, y) sorensen_dice(gsets[[x]], gsets[[y]])))
    dimnames(integration_matrix) <- list(names(gsets), names(gsets))
    integration_matrix[is.na(integration_matrix)] <- 0
    integration_matrix[is.infinite(integration_matrix)] <- 0
    integration_matrix <- pmax(pmin(integration_matrix, 1), 0)

    dist_matrix <- as.dist(1 - integration_matrix)
    hc <- hclust(dist_matrix, method = linkage_method)

    # Ensure that cluster heights are monotonically increasing.
    if (is.unsorted(hc$height)) {
      ord <- order(hc$height)
      hc$height <- hc$height[ord]
      hc$merge <- hc$merge[ord, ]
    }

    clust <- cutree(hc, h = 1 - minsds)
    tapply(seq_along(gsets), clust, function(i) {
      sort(unique(unlist(gsets[i])))
    })
  }

  # Merge up-regulated gene sets.
  if (verbose) cat("Merging up-regulated gene sets...\n")
  update_geneSets_up_filtered <- aux.merge(
    gsets = geneSets_up_filtered,
    minsds = sorenson_threshold,
    linkage_method = linkage_method
  )
  update_geneSets_up_filtered <- update_geneSets_up_filtered[!duplicated(lapply(update_geneSets_up_filtered, sort))]
  update_geneSets_up_filtered <- update_geneSets_up_filtered[sapply(update_geneSets_up_filtered, length) >= min_geneset_size]
  names(update_geneSets_up_filtered) <- paste0("ags", seq_along(update_geneSets_up_filtered))

  # Merge down-regulated gene sets.
  if (verbose) cat("Merging down-regulated gene sets...\n")
  update_geneSets_dn_filtered <- aux.merge(
    gsets = geneSets_dn_filtered,
    minsds = sorenson_threshold,
    linkage_method = linkage_method
  )
  update_geneSets_dn_filtered <- update_geneSets_dn_filtered[!duplicated(lapply(update_geneSets_dn_filtered, sort))]
  update_geneSets_dn_filtered <- update_geneSets_dn_filtered[sapply(update_geneSets_dn_filtered, length) >= min_geneset_size]
  names(update_geneSets_dn_filtered) <- paste0("igs", seq_along(update_geneSets_dn_filtered))

  if (verbose) {
    cat("\n=== Summary ===\n")
    cat("Activation gene sets:", length(update_geneSets_up_filtered), "\n")
    cat("Inhibition gene sets:", length(update_geneSets_dn_filtered), "\n")
  }

  # Return the direction-resolved gene sets and the supporting evidence tables.
  return(list(
    activation_gene_sets = update_geneSets_up_filtered,
    inhibition_gene_sets = update_geneSets_dn_filtered,
    importance_all = importance_all,
    importance_up_gene = importance_up_gene,
    importance_dn_gene = importance_dn_gene
  ))
}
