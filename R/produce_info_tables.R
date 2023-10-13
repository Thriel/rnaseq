
#' Produce tables required for batch analysis from minimal information
#'
#' @param de_infos Path to a csv file describing the DE analyses to perform.
#' @param design_infos Path to a csv file describing the groups for the
#' comparisons.
#' @param metadata_infos Path to a csv file of sample metadata.
#' @param volcano_fc_threshold Numeric vector of fold-change thresholds to use
#' in volcano plots. Default: \code{c(1.5, 3)}.
#' @param outdir The directory where to save produced tables in csv format. The
#' files will be saved as <outdir>/<name>.csv. If \code{NULL}, the results
#' won't be saved as csv. Default: \code{NULL}.
#'
#' @return Invisibly returns a named \code{list} of all the tables produced.
#'
#' @examples
#' \dontrun{
#' }
#'
#' @importFrom readr read_csv
#' @importFrom readr write_csv
#' @importFrom dplyr left_join
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom dplyr mutate
#' @importFrom dplyr bind_rows
#' @importFrom purrr map
#' @importFrom purrr map_dfr
#' @importFrom purrr reduce
#' @importFrom tidyselect one_of
#' @importFrom stringr str_replace
#' @importFrom stringr str_replace_all
#' @importFrom stringr str_to_title
#' @importFrom stringr str_detect
#' @importFrom tibble tibble
#'
#'
#' @export
produce_info_tables <- function(de_infos, design_infos, metadata_infos,
                            volcano_fc_threshold = c(1.5, 3), outdir = NULL){

    # Import input files
    de <- read_csv(de_infos)
    meta <- read_csv(metadata_infos, col_types = "cc")
    design <- read_csv(design_infos, col_types = "cc")

    # 1 Prepare metadata file
    mcols_meta <- map(de$id_de, \(x) prepare_metadata_mcols(x, de, design)) %>%
    purrr::reduce(left_join, by = "sample")
    stopifnot(all(mcols_meta$sample %in% meta$sample))
    stopifnot(all(meta$sample %in% mcols_meta$sample))
    meta_complete <- left_join(meta, mcols_meta, by = "sample")

    # 2 Prepare volcano file
    volcano <- prepare_volcano_info(de, volcano_fc_threshold)

    # 3 Prepare PCA file
    pca_infos <- prepare_pca_info(meta, de)

    # 4 Report infos
    report_infos <- prepare_report_info(
        pca_infos, de, volcano, meta, header_file = "header_experiment_1.qmd",
        dir_r_obj = "r_objects_experiment_1"
    )

    if (!is.null(outdir)) {
        dir.create(outdir, recursive = TRUE)
        write_csv(meta_complete, file.path(outdir, "metadata_infos.csv"))
        write_csv(volcano, file.path(outdir, "volcano_infos.csv"))
        write_csv(pca_infos, file.path(outdir, "pca_infos.csv"))
        write_csv(report_infos, file.path(outdir, "report_infos.csv"))
    }

    invisible(
        list(
            metadata = meta_complete, volcano = volcano, pca = pca_infos,
            report = report_infos
        )
    )
}

prepare_metadata_mcols <- function(current_id_de, de, design) {
  print(current_id_de)
  current_de <- filter(de, id_de == current_id_de)
  stopifnot(nrow(current_de) == 1)

  current_grp_col <- colnames(design)[colnames(design) == current_de$group]
  stopifnot(length(current_grp_col) == 1)

  current_grp1 <- current_de$contrast_1
  current_grp2 <- current_de$contrast_2

  current_samples_grp1 <- design$sample[pull(design, one_of(current_grp_col)) %in% current_grp1]
  current_samples_grp2 <- design$sample[pull(design, one_of(current_grp_col)) %in% current_grp2]
  current_samples <- c(current_samples_grp1, current_samples_grp2)

  res <- dplyr::select(design, sample) %>%
    mutate(tmp = if_else(sample %in% current_samples, "yes", "no"))
  colnames(res)[2] <- current_id_de
  res
}

prepare_volcano_info <- function(de_infos, fc_threshold) {
  id_plot <- paste0(rep(de_infos$id_de, each = length(fc_threshold)),
                    "_FC_",
                    rep(fc_threshold, length(de_infos$id_de))) %>%
    str_replace_all("\\.", "_")
  tibble(id_plot = id_plot,
         id_de = rep(de_infos$id_de, each = length(fc_threshold)),
         y_axis = "padj",
         p_threshold = 0.05,
         fc_threshold = rep(fc_threshold, length(de_infos$id_de)),
         title = str_replace(id_de, "_vs_", " versus ") %>% str_replace_all("_", " "),
         show_signif_counts = TRUE,
         show_signif_lines = "vertical",
         show_signif_colors = TRUE,
         col_up = "#E73426",
         col_down = "#0020F5",
         size = 3)
}

prepare_pca_info <- function(metadata, de_infos) {
  col <- dplyr::select(metadata, -sample) %>% colnames

  ## Global
  titles <- paste0("Global PCA ", c(" ", rep("- ", length(col))), c("", col)) %>%
    str_replace_all("_", " ") %>%
    str_replace_all(" - $", "") %>%
    str_to_title %>%
    str_replace("Pca", "PCA")
  ids_plot <- paste0("global_pca", c("", paste0("_", col)))
  global_pcas <- tibble(id_plot = ids_plot,
                        group = NA,
                        group_val = NA,
                        use_normalization = "none",
                        min_counts = 5,
                        id_metadata = "sample",
                        color = c(NA, col),
                        title = titles)
  ## DE
  ids_de <- de_infos$id_de
  ids_pca <- paste0("pca_",
                    rep(ids_de, length(col)+1),
                    "_",
                    rep(c("", col), each = length(ids_de))) %>%
    str_replace_all("_$", "")
  titles <- paste0(rep(str_replace(ids_de, " vs ", " versus "), length(col)+1),
                   " - ",
                   rep(c("", col), each = length(ids_de))) %>%
    str_replace_all(" - $", "")
  de_pcas <- tibble(id_plot = ids_pca,
                    group = rep(ids_de, length(col)+1),
                    group_val = "yes",
                    use_normalization = "none",
                    min_counts = 5,
                    id_metadata = "sample",
                    color = rep(c(NA, col), each = length(ids_de)),
                    title = titles)

  bind_rows(global_pcas, de_pcas)
}

prepare_report_info <- function(pca_infos, de_infos, volcano_infos, metadata, header_file = "header.qmd", dir_r_obj = "r_objects") {
  ## Header
  header_report <- tibble(add = "file", value = header_file)

  ## Global PCA
  global_pca <- filter(pca_infos, str_detect(id_plot, "^global_pca"))

  global_pca_add <- c("text", rep("plot", nrow(global_pca)))
  global_pca_values <- c("## Global PCA",
                         paste0(dir_r_obj, "/", global_pca$id_plot, ".rds"))
  global_pca_report <- tibble(add = global_pca_add, value = global_pca_values)

  ## DE results
  col <- dplyr::select(metadata, -sample) %>% colnames

  de_results_header <- tibble(add = "text", value = "## Differential expression analysis results")

  produce_de_res_report <- function(current_id_de) {
    current_header_report <- tibble(add = "text", value = paste0("### ", current_id_de))

    ### PCA
    current_pca_values <- paste0(dir_r_obj, "/pca_", current_id_de, "_", c("", col), ".rds") %>%
      str_replace("_\\.rds$", "\\.rds")
    current_pca_report <- tibble(add = c("text", rep("plot", length(col) + 1)),
                                 value = c("#### PCA", current_pca_values))

    ### Volcano
    current_volcano <- filter(volcano_infos, id_de == current_id_de)
    current_volcano_add <- rep("plot", nrow(current_volcano))
    current_volcano_values <- paste0(dir_r_obj, "/", current_volcano$id_plot, ".rds")
    current_volcano_report <- tibble(add = c("text", current_volcano_add),
                                     value = c("#### Volcano plots", current_volcano_values))

    ### Merge everything
    bind_rows(current_header_report, current_pca_report, current_volcano_report)
  }
  de_res_report <- map_dfr(de_infos$id_de, produce_de_res_report)

  bind_rows(header_report, global_pca_report, de_res_report)
}

