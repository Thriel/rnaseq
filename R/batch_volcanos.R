#' Produce volcano plots in batch
#'
#' The goal of this function is to take a csv file or a \code{data.frame}
#' describing all the volcano plots to produce and launch them in batches.
#'
#' The table may contain the following columns, in no specific order:
#'   * id_plot: unique identifier for each volcano plot. Mandatory.
#'   * id_de: unique identifier for each DE analysis results. Mandatory.
#'   Default: \code{TRUE}.
#'   * y_axis: Column in the DE results to use for the y-axis. "padj" or
#'   "pvalue" Default: "padj"
#'   * p_threshold: Maximal pvalue or padj to be considered significative.
#'   Default: 0.05
#'   * fc_threshold: Minimal fold-change to be considered significative.
#'   Default: 1.5
#'   * title: The title of the volcano. Default: NA
#'   * show_signif_counts: show the number of up- and down-regulated genes?
#'   * show_signif_lines: show lines at the threshold for significance? "none",
#'   "vertical","horizontal" or "both". Default: "vertical"
#'   * show_signif_color: Show color for significant genes? Default: TRUE
#'   * col_up: Color of the up-regulated genes. Default: "#E73426"
#'   * col_down: Color of the down-regulated genes. Default: "#0020F5"
#'   * size: The point size. Default: 3
#'
#' @param volcano_infos A csv file or a \code{data.frame} describing the
#' volcanos to produce.
#' @param de_results A \code{list} of DE results where the name of each DE
#' corresponds to the id_de value. Or the \code{outdir} directory where
#' the csv files of the \code{batch_de} function where saved.
#' @param add_labels A vector of the symbols to show on the volcano plot. If
#' NULL, no symbols will be shown. Default: NULL
#' @param outdir The directory where to save the volcano plots in pdf format.
#' If \code{NULL}, the results won't be saved as pdf. The files will be saved
#' as <outdir>/<id_plot>.pdf. Default: \code{NULL}.
#' @param r_objects The directory where to save volcano plots in rds format. If
#' \code{NULL}, the volcano plots won't be saved as rds files. The files will
#' be saved as <r_objects>/<id_plot>.rds. Default: \code{NULL}.
#' @param force Should the files be re-created if they already exists? Default:
#' \code{FALSE}.
#' @param cores Number of cores for the volcano creation. Default: 1
#'
#' @return Invisibly returns a \code{list} of all the volcano plots.
#'
#' @examples
#' \dontrun{
#' # First launch batch_de
#' de_infos <- get_demo_de_infos_file()
#' txi <- get_demo_txi()
#' design <- get_demo_design()
#' de_list <- batch_de(de_infos, txi, design)
#' # Then do the volcano_batch
#' volcano_infos <- get_demo_volcano_infos_file()
#' volcano_list <- batch_volcano(volcano_infos, de_list)
#' }
#'
#' @importFrom readr read_csv
#' @importFrom purrr map
#'
#' @importFrom parallel mclapply
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr full_join
#' @importFrom dplyr select
#' @importFrom dplyr everything
#' @importFrom readr write_csv
#'
#' @export
batch_volcano <- function(volcano_infos, de_results, add_labels = NULL,
                          outdir = NULL, r_objects = NULL, force = FALSE,
                          cores = 1) {

    # 1. Data validation
    ## volcano_infos
    stopifnot(is(volcano_infos, "data.frame") | is(volcano_infos, "character"))
    if (is.character(volcano_infos)) {
        stopifnot(file.exists(volcano_infos))
        volcano_infos <- readr::read_csv(volcano_infos, show_col_types = FALSE)
    }
    stopifnot(nrow(volcano_infos) > 0)
    stopifnot("id_plot" %in% colnames(volcano_infos))
    stopifnot("id_de" %in% colnames(volcano_infos))
    stopifnot(!any(is.na(volcano_infos$id_plot)))
    stopifnot(!any(is.na(volcano_infos$id_de)))
    stopifnot(!any(duplicated(volcano_infos$id_plot)))

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

    ## add_labels
    if (!is.null(add_labels)) {
        stopifnot(is(add_labels, "character"))
        stopifnot(all(nchar(add_labels) > 0))
    }

    ## force
    stopifnot(is(force, "logical"))
    stopifnot(!is.na(force))

    ## cores
    stopifnot(is(cores, "numeric"))
    stopifnot(cores == round(cores))
    stopifnot(cores > 0)

    ## de_results
    stopifnot(is(de_results, "list") | is(de_results, "character"))
    if (is.character(de_results)) {
        unique_de_ids <- unique(volcano_infos$id_de)
        de_files <- paste0(de_results, "/", unique_de_ids, ".csv")
        names(de_files) <- unique_de_ids
        stopifnot(all(file.exists(de_files)))
        de_results <- purrr::map(de_files, readr::read_csv,
                                 show_col_types = FALSE)
    }

    # Complete volcano infos
    volcano_infos <- complete_volcano_infos(volcano_infos)
    stopifnot(length(validate_volcano_infos(volcano_infos, de_results)) == 0)

    # Produce the volcano plots
    res <- list()
    volcano_analysis <- function(i) {
        current_id <- volcano_infos$id_plot[i]
        output_pdf <- paste0(outdir, "/", current_id, ".pdf")
        output_rds <- paste0(r_objects, "/", current_id, ".rds")

        if (!force & !is.null(r_objects) & file.exists(output_rds)) {
            current_volcano <- readRDS(output_rds)
        } else {
            current_infos <- volcano_infos[i,,drop=FALSE]
            current_volcano <- produce_single_volcano_batch(current_infos,
                                                            de_results,
                                                            add_labels, outdir)
        }
        if (!is.null(outdir)) {
            if (!file.exists(output_pdf) | force) {
                pdf(output_pdf)
                print(current_volcano$p)
                dev.off()
            }
        }
        if (!is.null(r_objects)) {
            if (!file.exists(output_rds) | force) {
                saveRDS(current_volcano$p, output_rds)
            }
        }
        res[[current_id]] <- current_volcano$p
    }
    res_volcano <- parallel::mclapply(1:nrow(volcano_infos), volcano_analysis, mc.cores = cores)
    names(res_volcano) <- volcano_infos$id_plot
    invisible(res_volcano)
}

complete_volcano_infos <- function(volcano_infos) {
    # Fill missing
    if (!"y_axis" %in% colnames(volcano_infos))
        volcano_infos[["y_axis"]] <- "padj"
    if (!"p_threshold" %in% colnames(volcano_infos))
        volcano_infos[["p_threshold"]] <- 0.05
    if (!"fc_threshold" %in% colnames(volcano_infos))
        volcano_infos[["fc_threshold"]] <- 1.5
    if (!"show_signif_counts" %in% colnames(volcano_infos))
        volcano_infos[["show_signif_counts"]] <- TRUE
    if (!"show_signif_lines" %in% colnames(volcano_infos))
        volcano_infos[["show_signif_lines"]] <- "vertical"
    if (!"show_signif_color" %in% colnames(volcano_infos))
        volcano_infos[["show_signif_color"]] <- TRUE
    if (!"col_up" %in% colnames(volcano_infos))
        volcano_infos[["col_up"]] <- "#E73426"
    if (!"col_down" %in% colnames(volcano_infos))
        volcano_infos[["col_down"]] <- "#0020F5"
    if (!"size" %in% colnames(volcano_infos))
        volcano_infos$size <- 3
    if (!"title" %in% colnames(volcano_infos))
        volcano_infos$title <- NA
    volcano_infos
}

validate_volcano_infos <- function(volcano_infos, design, txi) {
    errors <- list()
    for (i in seq_along(volcano_infos$id_de)) {
        current_id  <- volcano_infos$id_de[i]

        # y_axis checks
        current_y_axis <- volcano_infos$y_axis[i]

        if (!is(current_y_axis, "character")) {
            msg <- "y_axis value should be in character format"
            errors[[current_id]] <- c(errors[[current_id]], msg)
        } else {
            if (!current_y_axis %in% c("pvalue", "padj")) {
                msg <- 'Invalid y_axis value. Should be "pvalue" or "padj".'
                errors[[current_id]] <- c(errors[[current_id]], msg)
            }
        }

        # p_threshold checks
        current_p_threshold <- volcano_infos$p_threshold[i]
        if (!is(current_p_threshold, "numeric")) {
            msg <- "p_threshold value should be in numeric format"
            errors[[current_id]] <- c(errors[[current_id]], msg)
        } else {
            if (current_p_threshold < 0 | current_p_threshold > 1) {
                msg <- "p_threshold value should be between 0 and 1"
                errors[[current_id]] <- c(errors[[current_id]], msg)
            }
        }

        # fc_threshold checks
        current_fc_threshold <- volcano_infos$fc_threshold[i]
        if (!is(current_fc_threshold, "numeric")) {
            msg <- "fc_threshold value should be in numeric format"
            errors[[current_id]] <- c(errors[[current_id]], msg)
        } else {
            if (current_fc_threshold <= 0) {
                msg <- "fc_threshold value should be greater or equal to 0"
                errors[[current_id]] <- c(errors[[current_id]], msg)
            }
        }

        # title checks
        current_title <- volcano_infos$title[i]
        if (!is.na(current_title)) {
            if (!is(current_title, "character")) {
                msg <- "title must be a character"
                errors[[current_id]] <- c(errors[[current_id]], msg)
            }
        }

        # show_signif_counts checks
        current_show_signif_counts <- volcano_infos$show_signif_counts[i]
        if (!is(current_show_signif_counts, "logical")) {
            msg <- "show_signif_counts value should be in logical format"
            errors[[current_id]] <- c(errors[[current_id]], msg)
        }

        # show_signif_lines checks
        current_show_signif_lines <- volcano_infos$show_signif_lines[i]
        if (!is(current_show_signif_lines, "character")) {
            msg <- "show_signif_lines value should be in character format"
            errors[[current_id]] <- c(errors[[current_id]], msg)
        } else {
            expected_values <- c("none", "both", "vertical", "horizontal")
            if (!current_show_signif_lines %in% expected_values) {
                msg <- "show_signif_lines value show be one of: "
                msg <- paste0(msg, '"none", "both", "vertical", "horizontal"')
                errors[[current_id]] <- c(errors[[current_id]], msg)
            }
        }

        # show_signif_color checks
        current_show_signif_color <- volcano_infos$show_signif_color[i]
        if (!is(current_show_signif_color, "logical")) {
            msg <- "show_signif_color value should be in logical format"
            errors[[current_id]] <- c(errors[[current_id]], msg)
        }

        # color checks
        ## col_up checks
        current_col_up <- volcano_infos$col_up[i]
        if (!is(current_col_up, "character")) {
            msg <- "col_up value should be in character format"
            errors[[current_id]] <- c(errors[[current_id]], msg)
        } else {
            if (!is_color(current_col_up)) {
                msg <- "col_up value is not a valid color"
                errors[[current_id]] <- c(errors[[current_id]], msg)
            }
        }

        ## col_down checks
        current_col_down <- volcano_infos$col_down[i]
        if (!is(current_col_down, "character")) {
            msg <- "col_down value should be in character format"
            errors[[current_id]] <- c(errors[[current_id]], msg)
        } else {
            if (!is_color(current_col_down)) {
                msg <- "col_down value is not a valid color"
                errors[[current_id]] <- c(errors[[current_id]], msg)
            }
        }

        # size checks
        current_size <- volcano_infos$size[i]

        if (!is.na(current_size)) {
            if (!is(current_size, "numeric")) {
                msg <- "size must be a numeric value"
                errors[[current_id]] <- c(errors[[current_id]], msg)
            } else {
                if (!identical(current_size, round(current_size))) {
                    msg <- "size must be an integer"
                    errors[[current_id]] <- c(errors[[current_id]], msg)
                }
                if (current_size <= 0) {
                    msg <- "size must be greater than 0"
                    errors[[current_id]] <- c(errors[[current_id]], msg)
                }
            }
        } else {
            msg <- "size must not be NA"
            errors[[current_id]] <- c(errors[[current_id]], msg)
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

produce_single_volcano_batch <- function(current_volcano_info, de_results, add_labels, outdir) {
    cvi <- current_volcano_info
    de_res <- de_results[[cvi$id_de]]

    # TODO: add_labels
    produce_volcano(de_res,
                    fc_threshold = cvi$fc_threshold,
                    p_threshold = cvi$p_threshold,
                    y_axis = cvi$y_axis,
                    show_signif_counts = cvi$show_signif_counts,
                    show_signif_lines = cvi$show_signif_lines,
                    show_signif_color = cvi$show_signif_color,
                    col_up = cvi$col_up,
                    col_down = cvi$col_down,
                    size = cvi$size,
                    title = cvi$title,
                    graph = FALSE,
                    outdir = outdir)

    #add faceting, modify

}
