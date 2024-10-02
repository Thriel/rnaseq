#' DESeq2 analysis
#'
#' @param txi The txi object returned by the import_kallisto function.
#' @param design The experimental design (see ?DESeqDataSetFromTximport).
#' @param formula The design formula in data.frame format (see
#'                ?DESeqDataSetFromTximport).
#' @param filter The minimum number of reads detected for a feature across all
#'               samples. Default: 2
#' @param count_matrix Use an alternative count matrix to use for the
#' differential analysis instead of \code{txi$counts}. Will use the
#' \code{DESeq2::DESeqDataSetFromMatrix} instead of the
#' \code{DESeq2::DESeqDataSetFromTximport} function, so will work even if txi
#' object is incomplete (i.e.: length matrix is missing).  Default: \code{NA}.
#' @param ... Extra param for the DESeq2::DESeq function
#'
#' @return A DESeqDataSet object.
#'
#' @examples
#' txi <- get_demo_txi()
#' design <- get_demo_design()
#' dds <- deseq2_analysis(txi, design, ~ group)
#' de <- DESeq2::results(dds, contrast = c("group", "A", "B"))
#'
#' @importFrom DESeq2 DESeqDataSetFromMatrix
#' @importFrom DESeq2 DESeqDataSetFromTximport
#' @importFrom DESeq2 counts
#' @importFrom DESeq2 DESeq
#'
#' @export
deseq2_analysis <- function(txi, design, formula, filter = 2,
                            count_matrix = NULL, ...) {
    validate_txi(txi)
    stopifnot(all(c("sample") %in% colnames(design)))
    stopifnot(identical(colnames(txi$counts), as.character(design$sample)))
    if (!is.null(count_matrix)) {
        stopifnot(is(count_matrix, "character"))
        stopifnot(count_matrix %in% names(txi))
        stopifnot(is(txi[[count_matrix]], "matrix"))
    }
    if (!is.null(txi$dummy)) {
        stopifnot(is(txi$dummy, "character"))
        if (is.null(count_matrix)) {
            stopifnot("count" %in% txi$dummy)
            stopifnot("length" %in% txi$dummy)
        } else {
            stopifnot(count_matrix %in% txi$dummy)
        }
    }

    if (!is.null(count_matrix)) {
        dds <- DESeq2::DESeqDataSetFromMatrix(txi[[count_matrix]], design,
                                              formula)
    } else {
        dds <- DESeq2::DESeqDataSetFromTximport(txi, design, formula)
    }
    dds <- dds[rowSums(DESeq2::counts(dds)) >= filter]
    dds <- DESeq2::DESeq(dds, ...)
    dds
}

#' Split differential expression results
#'
#' This function will split the DE results table into 4 elements:
#'     1) de: The original DE table
#'     2) signif: All the genes considered statistically differentially
#'                expressed
#'     3) up: The significant genes that are up-regulated
#'     4) down: The significant genes that are down-regulated
#'
#' To be considered significant, a gene must have a padj (qV) value lower or
#' equal to the \code{p_threshold} param, an absolute fold change greater or
#' equal to the \code{fc_threshold} param. Also, if \code{tpm_threshold} is not
#' \code{NULL}, the average TPM across all samples in the comparison must be
#' greater or equal to the \code{tpm_threshold} value.
#'
#' It is important to use the results produced with the
#' \code{format_de_results} function.
#'
#' @param de_res the \code{data.frame} object returned by the
#' \code{format_de_results} function.
#' @param fc_threshold The threshold of FC to be considered as significant.
#' Default: 0.05
#' @param p_threshold The threshold of p stat to be considered as significant.
#' Default: 1.5
#' @param tpm_threshold The threshold of the mean TPM of the current samples to
#' be considered as significant
#'
#' @return A list of DE tables.
#'
#' @importFrom dplyr filter
#'
#' @examples
#' \dontrun{
#' txi <- get_demo_txi()
#' design <- get_demo_design()
#' dds <- deseq2_analysis(txi, design, ~ group)
#' de <- format_de_results(dds, txi, c("group", "A", "B"))
#' split_de <- split_de_results(de)
#' }
#'
#' @export
split_de_results <- function(de_res, p_threshold = 0.05, fc_threshold = 1.5,
                             tpm_threshold = NULL) {
    print(colnames(de_res))
    stopifnot(is(de_res, "data.frame"))
    expected_cols <- c("padj", "log2FoldChange")
    print(expected_cols)
    stopifnot(all(expected_cols %in% colnames(de_res)))
    stopifnot(is(p_threshold, "numeric"))
    stopifnot(p_threshold >= 0 & p_threshold <= 1)
    stopifnot(is(fc_threshold, "numeric"))
    stopifnot(fc_threshold >= 0)
    if (!is.null(tpm_threshold)) {
      stopifnot(is(tpm_threshold, "numeric"))
      stopifnot(tpm_threshold >= 0)
      }
    
    idx_p <- de_res$padj <= p_threshold
    idx_fc <- abs(de_res$log2FoldChange) >= fc_threshold
    de_res$p_threshold <- FALSE
    de_res$p_threshold[idx_p] <- TRUE
    de_res$fc_threshold <- FALSE
    de_res$fc_threshold[idx_fc] <- TRUE
    
    if (!is.null(tpm_threshold)) {
      mean_tpm_1 <- de_res$mean_TPM_grp1
      mean_tpm_2 <- de_res$mean_TPM_grp2
      idx_tpm <- mean(mean_tpm_1, mean_tpm_2) >= tpm_threshold
      de_res$tpm_threshold <- FALSE
      de_res$tpm_threshold[idx_tpm] <- TRUE
      signif <- dplyr::filter(de_res, p_threshold, fc_threshold,
                              tpm_threshold)
    } else {
      signif <- dplyr::filter(de_res, p_threshold, fc_threshold)
    }
    
    up <- dplyr::filter(signif, log2FoldChange > 0)
    down <- dplyr::filter(signif, log2FoldChange < 0)
    
    list(de_res = de_res, signif = signif, up = up, down = down)
}
#' Minimal formatting of de results
#'
#' This function will call DESeq2::results on the dds object and add the
#' informations from txi$anno using a full join. This means that the id that
#' were removed in the dds production step by filtering the rows with a small
#' number of counts will be re-included in the results. In this case, the
#' log2FoldChange, the pvalue and the padj columns will all be NA.
#'
#' @param dds The DESeqDataSet object returned by deseq2_analysis.
#' @param txi The txi object returned by the import_kallisto function.
#' @param contrast A vector describing the contrasts in the c("<comp_group>",
#' "<comp1>", "<comp2>") format.
#' @param keep_stats Keep baseMean, lfcSE and stat values in the results?
#' Default: \code{TRUE}.
#' @param add_mean_dds Add the mean DESeq normalization value for each group of
#' the comparison. Default: \code{FALSE}
#'
#' @return A data.frame with the id, ensembl_gene, symbol, entrez_id,
#' transcript_type, log2FoldChange, pvalue, padj columns.
#'
#' @examples
#' txi <- get_demo_txi()
#' design <- get_demo_design()
#' dds <- deseq2_analysis(txi, design, ~ group)
#' de_res <- format_de_results(dds, txi,  c("group", "A", "B"))
#'
#' @importFrom DESeq2 results
#' @importFrom DESeq2 counts
#' @importFrom SummarizedExperiment colData
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr full_join
#' @importFrom dplyr relocate
#' @importFrom dplyr last_col
#' @importFrom dplyr select
#'
#' @export
format_de_results <- function(dds, txi, contrast, keep_stats = TRUE, add_mean_dds = FALSE) {
    stopifnot(is(dds, "DESeqDataSet"))
    validate_txi(txi)
    stopifnot(is(contrast, "character"))
    stopifnot(length(contrast) == 3)

    de_res <- DESeq2::results(dds, contrast = contrast) %>%
        as.data.frame %>%
        tibble::rownames_to_column("id")
    stopifnot(all(de_res$id %in% txi$anno$id))

    de_res <- de_res %>%
        dplyr::full_join(txi$anno, by = "id")
    if (keep_stats) {
        de_res <- de_res %>%
            dplyr::relocate(baseMean, lfcSE, stat, log2FoldChange, pvalue,
                            padj, .after = dplyr::last_col())
    } else {
        de_res <- de_res %>%
            dplyr::select(-baseMean, -lfcSE, -stat) %>%
            dplyr::relocate(log2FoldChange, pvalue, padj,
                            .after = dplyr::last_col())
    }
    if (add_mean_dds) {
        dds_counts <- DESeq2::counts(dds, normalized = TRUE)
        design <- SummarizedExperiment::colData(dds)
        i <- design[,contrast[1], drop = TRUE] == contrast[2]
        j <- design[,contrast[1], drop = TRUE] == contrast[3]
        mean_dds_grp1 <- rowMeans(dds_counts[,i])
        mean_dds_grp2 <- rowMeans(dds_counts[,j])
        stopifnot(all(names(mean_dds_grp1) %in% de_res$id))
        stopifnot(all(names(mean_dds_grp2) %in% de_res$id))
        mean_dds_grp1 <- mean_dds_grp1[de_res$id]
        mean_dds_grp2 <- mean_dds_grp2[de_res$id]
        de_res$mean_dds_grp1 <- mean_dds_grp1
        de_res$mean_dds_grp2 <- mean_dds_grp2
    }

    validate_de(de_res, txi, keep_stats, add_mean_dds)
    de_res
}

validate_de <- function(de, txi, keep_stats, add_mean_dds = FALSE) {
    stopifnot(is(de, "data.frame"))

    expected_cols <- c("id", "ensembl_gene", "symbol", "entrez_id",
                       "transcript_type", "log2FoldChange", "pvalue", "padj")
    if (keep_stats) {
        expected_cols <- c(expected_cols, "baseMean", "lfcSE", "stat")
    }
    if (add_mean_dds) {
        expected_cols <- c(expected_cols, "mean_dds_grp1", "mean_dds_grp2")
    }
    stopifnot(all(expected_cols %in% colnames(de)))

    stopifnot(is(de$id, "character"))
    stopifnot(is(de$ensembl_gene, "character"))
    stopifnot(is(de$symbol, "character"))
    stopifnot(is(de$entrez_id, "character") | is(de$entrez_id, "numeric"))
    stopifnot(is(de$log2FoldChange, "numeric"))
    stopifnot(is(de$pvalue, "numeric"))
    stopifnot(is(de$padj, "numeric"))
    if (keep_stats) {
        stopifnot(is(de$baseMean, "numeric"))
        stopifnot(is(de$lfcSE, "numeric"))
        stopifnot(is(de$stat, "numeric"))
    }
    if (add_mean_dds) {
        stopifnot(is(de$mean_dds_grp1, "numeric"))
        stopifnot(is(de$mean_dds_grp2, "numeric"))
    }

    stopifnot(all(txi$anno$id %in% de$id))
    stopifnot(all(de$id %in% txi$anno$id))

    invisible(TRUE)
}
