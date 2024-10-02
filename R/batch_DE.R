#' Produce differential expression analysis in batch
#'
#' The goal of this function is to take a csv file or a \code{data.frame}
#' describing all the differential expression (DE) analysis to produce and
#' launch them in batches.
#'
#' If the \code{r_objects} param is used, the dds object will be saved and
#' won't be recalculated if the function is called again, unless the
#' \code{force} param is set to \code{TRUE}.
#'
#' The table may contain the following columns, in no specific order:
#'   * id_de: unique identifier for this DE analysis. Mandatory.
#'   * group: column from the design file to define the group for the current
#'   DE analysis. Mandatory.
#'   * contrast_1: The first group for the DE analysis. Must be present in the
#'   group column. Mandatory.
#'   * contrast_2: The first group for the DE analysis. Must be present in the
#'   group column. Mandatory.
#'   * formula: The formula for the de_analysis. Default: ~ group
#'   * filter: Minimum mean number of reads across all samples to be part of
#'   the DE analysis
#'   * count_matrix: The count matrix to use for the differential analysis.
#'   Will use the \code{DESeq2::DESeqDataSetFromMatrix} instead of the
#'   \code{DESeq2::DESeqDataSetFromTximport} function, so will work even if txi
#'   object is incomplete (i.e.: length matrix is missing). Default: NA
#'
#' @param de_infos A csv file or a \code{data.frame} describing the DE analysis
#' to perform.
#' @param txi The txi object returned by the import_kallisto function.
#' @param design A csv file of a \code{data.frame} describing the groups for
#' the comparisons.
#' @param outdir The directory where to save DE analysis in csv format. If
#' \code{NULL}, the results won't be saved as csv The files will be saved as
#' <outdir>/<id_de>.csv. Default: \code{NULL}.
#' @param r_objects The directory where to save dds in rds format. If
#' \code{NULL}, the dds won't be saved as rds files. The files will be saved as
#' <r_objects>/<id_de>.rds. Default: \code{NULL}.
#' @param force Should the files be re-created if they already exists? Default:
#' \code{FALSE}.
#' @param cores Number of cores for the DE analysis. Default: 1
#'
#' @return Invisibly returns a \code{list} of all the DE results.
#'
#' @examples
#' \dontrun{
#' de_infos <- get_demo_de_infos_file()
#' txi <- get_demo_txi()
#' design <- get_demo_design()
#' de_list <- batch_de(de_infos, txi, design)
#' }
#'
#' @importFrom readr read_csv
#' @importFrom stringr str_detect
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr full_join
#' @importFrom dplyr select
#' @importFrom dplyr everything
#' @importFrom readr write_csv
#' @importFrom parallel mclapply
#' @importFrom dplyr filter
#' @importFrom tibble column_to_rownames
#' @importFrom DESeq2 results
#' @importFrom dplyr left_join
#'
#' @export
batch_de <- function(de_infos, txi, design, outdir = NULL, r_objects = NULL,
                     force = FALSE, cores = 1, return_normalized_counts = FALSE) {
  
  # 1. Data validation
  ## de_infos
  stopifnot(is(de_infos, "data.frame") | is(de_infos, "character"))
  if (is.character(de_infos)) {
    stopifnot(file.exists(de_infos))
    de_infos <- readr::read_csv(de_infos, show_col_types = FALSE)
  }
  stopifnot(nrow(de_infos) > 0)
  stopifnot("id_de" %in% colnames(de_infos))
  stopifnot(!any(is.na(de_infos$id_de)))
  stopifnot(!any(duplicated(de_infos$id_de)))
  expected_cols <- c("group", "contrast_1", "contrast_2")
  stopifnot(all(expected_cols %in% colnames(de_infos)))
  
  ## txi
  validate_txi(txi)
  
  ## design
  stopifnot(is(design, "data.frame") | is(design, "character"))
  if (is.character(design)) {
    stopifnot(file.exists(design))
    design <- readr::read_csv(design, show_col_types = FALSE)
  }
  stopifnot(nrow(design) > 0)
  
  ## Directories
  if (!is.null(outdir)) {
    stopifnot(is(outdir, "character"))
    if (!dir.exists(outdir)) {
      dir.create(outdir)
    }
  }
  if (!is.null(r_objects)) {
    stopifnot(is(r_objects, "character"))
    if (!dir.exists(r_objects)) {
      dir.create(r_objects)
    }
  }
  
  ## force
  stopifnot(is(force, "logical"))
  stopifnot(!is.na(force))
  
  ## cores
  stopifnot(is(cores, "numeric"))
  stopifnot(cores == round(cores))
  stopifnot(cores > 0)
  
  # Complete DE infos
  de_infos <- complete_de_infos(de_infos)
  stopifnot(length(validate_de_infos(de_infos, design, txi)) == 0)
  
  # Produce the DE
  res <- list()
  #Initier le dataframe histo_data
  
  shared_env <- new.env()
  shared_env$histo_data <- data.frame(sample = character(), sum_counts = numeric())
  
  
  de_analysis <- function(i) {
    current_id <- de_infos$id_de[i]
    output_rds <- paste0(r_objects, "/", current_id, ".rds")
    output_rds_de <- paste0(r_objects, "/", current_id, "_de.rds")
    if (!force & !is.null(r_objects) & file.exists(output_rds)) {
      dds <- readRDS(output_rds)
    } else {
      dds <- NULL
    }
    if (!force & !is.null(r_objects) & file.exists(output_rds_de)) {
      de <- readRDS(output_rds_de)
    } else {
      de <- NULL
    }
    res_de <- produce_single_de_batch(de_infos[i,,drop=FALSE],
                                      txi, design, dds, de)
    
    if (!is.null(outdir)) {
      output_csv <- paste0(outdir, "/", current_id, ".csv")
      if (!file.exists(output_csv) | force) {
        tmp <- res_de$de %>% as.data.frame()
        
        # replace qV and pV by padj and pvalue
        i <- stringr::str_detect(colnames(tmp), "qV")
        stopifnot(sum(i) %in% c(0,1))
        if (sum(i) == 1) {
          colnames(tmp)[i] <- "padj"
        }
        i <- stringr::str_detect(colnames(tmp), "pV")
        stopifnot(sum(i) %in% c(0,1))
        if (sum(i) == 1) {
          colnames(tmp)[i] <- "pvalue"
        }
        
        # Ajouter les comptages normalisés si demandé
        if (return_normalized_counts) {
          normalized_counts <- counts(res_de$dds, normalized = TRUE)
          print("TEST")
          histo_data <- shared_env$histo_data
          
          
          # Vérifier et ajouter chaque colonne de normalized_counts à histo_data
          for (col_name in colnames(normalized_counts)) {
            if (!any(histo_data$sample == col_name)) {
              col_sum <- sum(normalized_counts[, col_name], na.rm = TRUE)
              new_row <- data.frame(sample = col_name, sum_counts = col_sum)
              histo_data <- rbind(histo_data, new_row)
            }
          }
          print(histo_data)
          norm_counts_df <- as.data.frame(normalized_counts)
          
          # Assurez-vous que tmp a une colonne id
          if(!("id" %in% colnames(tmp))){
            tmp <- tmp %>% tibble::rownames_to_column("id")
          }
          norm_counts_df$id <- rownames(norm_counts_df)
          
          # Fusionner tmp et norm_counts_df
          combined_results <- dplyr::full_join(tmp, norm_counts_df, by = "id")
          
          # Réorganiser les colonnes par groupe
          sample_data <- colData(res_de$dds)
          group_order <- order(sample_data$group)
          sample_columns <- colnames(norm_counts_df)[colnames(norm_counts_df) != "id"]
          reordered_sample_columns <- sample_columns[group_order]
          
          # Créer le nouvel ordre des colonnes
          new_column_order <- c("id", colnames(tmp)[colnames(tmp) != "id"], reordered_sample_columns)
          
          # Réorganiser les colonnes
          combined_results <- combined_results[, new_column_order]
        } else {
          combined_results <- tmp
        }
        
        # Call split_de_results function
        split_results <- split_de_results(tmp, p_threshold = 0.05, fc_threshold = 1.5, tpm_threshold = NULL)
        
        
        # Écrire le fichier CSV combiné
        readr::write_csv(combined_results, output_csv)
        readr::write_csv(split_results$up, paste0(outdir, "/", current_id, "_up.csv"))
        readr::write_csv(split_results$down, paste0(outdir, "/", current_id, "_down.csv"))
      }
    }
    if (!is.null(r_objects)) {
      if (!file.exists(output_rds) | force) {
        saveRDS(res_de$dds, output_rds)
      }
    }
    if (!is.null(r_objects)) {
      if (!file.exists(output_rds_de) | force) {
        saveRDS(res_de$de, output_rds_de)
      }
    }
    # Update the shared environment
    shared_env$histo_data <- histo_data
    shared_env$split_results[[current_id]] <- split_results
    res[[current_id]] <- res_de$de
    
    list(de = res_de$de, histo_data = histo_data, split_results = split_results)
    
  }
  results <- parallel::mclapply(1:nrow(de_infos), de_analysis, mc.cores = cores)
  names(results) <- de_infos$id_de
  
  # Extraire les résultats DE
  de_list <- lapply(results, function(x) x$de)
  
  # Combiner les histo_data
  all_histo_data <- do.call(rbind, lapply(results, function(x) x$histo_data))
  histo_data <- aggregate(sum_counts ~ sample, data = all_histo_data, FUN = sum)
  
  # Extraire les split_results
  split_results <- lapply(results, function(x) x$split_results)
  
  return(list(de_list = de_list, histo_data = histo_data, split_results = split_results))
}

complete_de_infos <- function(de_infos) {
  # Fill missing
  if (!"formula" %in% colnames(de_infos))
    de_infos[["formula"]] <- "~ group"
  if (!"filter" %in% colnames(de_infos))
    de_infos[["filter"]] <- 2
  if (!"count_matrix" %in% colnames(de_infos))
    de_infos[["count_matrix"]] <- NA
  de_infos
}

validate_de_infos <- function(de_infos, design, txi) {
  errors <- list()
  for (i in seq_along(de_infos$id_de)) {
    current_id  <- de_infos$id_de[i]
    
    # Groups checks
    current_group <- de_infos$group[i]
    
    ## group
    if (!is(current_group, "character")) {
      msg <- "group column should be in character format"
      errors[[current_id]] <- c(errors[[current_id]], msg)
    } else {
      if (!current_group %in% colnames(design)) {
        msg <- "group column must be present in design"
        errors[[current_id]] <- c(errors[[current_id]], msg)
      }
    }
    
    ## contrasts
    current_contrast_1 <- de_infos$contrast_1[i]
    current_contrast_2 <- de_infos$contrast_2[i]
    
    if (!is(current_contrast_1, "character")) {
      msg <- "contrast_1 column should be in character format"
      errors[[current_id]] <- c(errors[[current_id]], msg)
    } else {
      if (!current_contrast_1 %in% design[[current_group]]) {
        msg <- "contrast_1 column not present in design group column"
        errors[[current_id]] <- c(errors[[current_id]], msg)
      }
    }
    if (!is(current_contrast_2, "character")) {
      msg <- "contrast_2 column should be in character format"
      errors[[current_id]] <- c(errors[[current_id]], msg)
    } else {
      if (!current_contrast_2 %in% design[[current_group]]) {
        msg <- "contrast_2 column not present in design group column"
        errors[[current_id]] <- c(errors[[current_id]], msg)
      }
    }
    
    ## formula
    current_formula <- de_infos$formula[1]
    
    if (!is(current_formula, "character") &
        !is(current_formula, "formula")) {
      msg <- "formula column should be in character or formula format"
      errors[[current_id]] <- c(errors[[current_id]], msg)
    }
    
    ## count_matrix
    current_count_matrix <- de_infos$count_matrix[1]
    
    if (!is.na(current_count_matrix)) {
      if (!is(current_count_matrix, "character")) {
        msg <- "count_matrix column should be in character format"
        errors[[current_id]] <- c(errors[[current_id]], msg)
      } else {
        if (!current_count_matrix %in% names(txi)) {
          # Validate counts, abundance, ruvg_counts, combat_counts or extra_matrix
          msg <- "count_matrix not found in txi"
          errors[[current_id]] <- c(errors[[current_id]], msg)
        }
      }
    }
  }
  
  if (length(errors) > 0) {
    for (i in seq_along(errors)) {
      message(paste0(names(errors)[i], ":"))
      for (j in seq_along(errors[[i]])) {
        message(paste0("    ", errors[[i]][j]))
      }
    }
  }
  errors
}

produce_single_de_batch <- function(current_de_info, txi, design, dds, de, return_normalized_counts = FALSE) {
  cdi <- current_de_info
  cd <- design[[cdi$group]]
  current_contrasts <- c(cdi$contrast_1, cdi$contrast_2)
  current_samples <- design$sample[cd %in% current_contrasts]
  txi <- filter_txi(txi, current_samples)
  design <- dplyr::filter(design, sample %in% current_samples) %>%
    tibble::column_to_rownames("sample") %>%
    .[colnames(txi$counts),,drop=FALSE] %>%
    tibble::rownames_to_column("sample")
  
  if (is.null(dds)) {
    de <- NULL
    
    if (is(cdi$formula, "character")) {
      formula <- eval(parse(text=cdi$formula))
    }
    
    design <- design[design[[cdi$group]] %in% c(current_contrasts),]
    
    if (is.na(cdi$count_matrix)) {
      count_matrix <- NULL
    } else {
      count_matrix <- cdi$count_matrix
    }
    
    dds <- deseq2_analysis(txi = txi,
                           design = design,
                           formula = formula,
                           filter = cdi$filter,
                           count_matrix = count_matrix)
  }
  
  if (is.null(de)) {
    contrast <- c(cdi$group, cdi$contrast_1, cdi$contrast_2)
    de <- format_de_results(dds, txi, contrast)
  }
  normalized_counts <- counts(dds, normalized = TRUE)
  return(list(dds = dds, de = de, normalized_counts = normalized_counts))  
}