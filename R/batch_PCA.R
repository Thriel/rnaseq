#' Produce PCA in batch
#'
#' The goal of this function is to take a csv file or a \code{data.frame}
#' describing all the PCA to produce and launch them in batches.
#'
#' The table may contain the following columns, in no specific order:
#'   * id_plot: unique identifier for this graph. Mandatory
#'   * group: column from the metadata file to subset the txi. If value is
#'            NA, all the samples will be used in the PCA. If value is not
#'            NA, group_val value must not be NA either. If NA, group_val
#'            must also be NA. Default: NA
#'   * group_val: Value in the group column to use in current PCA. If group
#'                is NA, group_val must also be NA. If group is not NA,
#'                group_val must be a valid column name from the metadata
#'                table. Default: NA
#'   * use_normalisation: "none", "ruvg" or "combat" (txi must be produced in
#'              consequence). Default: "none"
#'   * min_counts: Mean values of counts (TPM, ruvg or combat) to keep gene for
#'                 PCA. Default: 5
#'   * id_metadata: Column in metadata table that contains the sample names.
#'                  Must be present in the colnames of the metadata table.
#'                  Default: NA
#'   * size: The point size. Default: 3
#'   * shape: The column in the metadata file to use to define shapes.
#'            Default: NA
#'   * color: The column in the metadata file to use to define colors.
#'             Default: NA
#'   * title: The title of the PCA. Default: NA
#'   * legend.position: "left", "right", "top" or "bottom". Default: "right"
#'   * legend.box: "horizontal" or "vertical". Default: "vertical"
#'   * show_names: Show sample names on PCA? Default: TRUE
#'
#' Only the id_plot is mandatory. This value is used to name the PCA plots that
#' are created by the batch_pca function and the output and r_objects filename.
#'
#' If the other columns are absent, they will be automatically created and
#' filled with their default values.
#'
#' @param pca_infos A csv file or a \code{data.frame} describing the PCA to
#' produce.
#' @param txi The \code{txi} object returned by the \code{import_kallisto}
#' function.
#' @param metatada The metadata to add to the coordinate table to add color or
#' shape with the results and the \code{plot_pca} function. If there is no
#' metadata available, keep the default value of \code{NULL}. Value can either
#' be a csv file or a \code{data.frame}. Default: \code{NULL}.
#' @param outdir The directory where to save PCA in pdf format. If \code{NULL},
#' the plots won't be saved as pdf. The files will be saved as
#' <outdir>/<id_plot>.pdf. Default: \code{NULL}.
#' @param r_objects The directory where to save PCA in rds format. If
#' \code{NULL}, the plots won't be saved as rds files. The files will be saved
#' as <r_objects>/<id_plot>.rds. Default: \code{NULL}.
#' @param force Should the files be re-created if they already exists? Default:
#' \code{FALSE}.
#' @param cores Number of cores for the DE analysis. Default: 1
#'
#' @return Invisibly returns a \code{list} of all the ggplots.
#'
#' @examples
#' pca_infos <- get_demo_pca_infos_file()
#' txi <- get_demo_txi(large = TRUE)
#' metadata <- get_demo_metadata_file()
#' gg_list <- batch_pca(pca_infos, txi, metadata)
#'
#' @importFrom readr read_csv
#' @importFrom dplyr select
#' @importFrom dplyr left_join
#'
#' @export
batch_pca <- function(pca_infos, txi, metadata = NULL, outdir = NULL,
                      r_objects = NULL, force = FALSE, cores = 1) {
    # 1. Data validation
    ## pca_infos
    stopifnot(is(pca_infos, "data.frame") | is(pca_infos, "character"))
    if (is.character(pca_infos)) {
        stopifnot(file.exists(pca_infos))
        pca_infos <- readr::read_csv(pca_infos, show_col_types = FALSE)
    }
    stopifnot(nrow(pca_infos) > 0)
    stopifnot("id_plot" %in% colnames(pca_infos))
    stopifnot(!any(is.na(pca_infos$id_plot)))
    stopifnot(nrow(pca_infos) > 0)
    stopifnot(!any(duplicated(pca_infos$id_plot)))

    ## txi
    validate_txi(txi)

    ## metadata
    if (!is.null(metadata)) {
        stopifnot(is(metadata, "data.frame") | is(metadata, "character"))
        if (is.character(metadata)) {
            stopifnot(file.exists(metadata))
            metadata <- readr::read_csv(metadata, show_col_types = FALSE)
        }
        stopifnot(nrow(metadata) > 0)
    }

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

    # Complete PCA infos
    pca_infos <- complete_pca_infos(pca_infos)
    stopifnot(length(validate_pca_infos(pca_infos, metadata, txi)) == 0)

    # Produce the pca
    
    print("Starting PCA")
    min_pca_infos <- dplyr::select(pca_infos, group, group_val, id_metadata,
                                   use_normalisation, min_counts) %>%
                         unique
    create_pca_df <- function(i) {
        cg <- min_pca_infos$group[i]
        cgv <- min_pca_infos$group_val[i]
        cim <- min_pca_infos$id_metadata[i]
        cun <- min_pca_infos$use_normalisation[i]
        cmc <- min_pca_infos$min_counts[i] %>% as.character

        if (is.na(cg)) {
            cg <- "NA"
            cgv <- "NA"
        }

        pca_df <- list()
        pca_df[[cg]] <- list()
        pca_df[[cg]][[cgv]] <- list()
        current_dir <- paste0(cg, "/", cgv)
        current_prefix <- paste(cim, cun, cmc, sep = "_")
        pca_df[[cg]][[cgv]][[current_prefix]] <- list()
        current_pca_df <- NULL
        if (!is.null(r_objects)) {
            current_rds <- paste0(r_objects, "/min_pca_infos/", current_dir, "/",
                                  current_prefix, ".rds")
            if (!dir.exists(dirname(current_rds))) {
                dir.create(dirname(current_rds), recursive = TRUE)
            }
            if (file.exists(current_rds) & !force) {
                current_pca_df <- readRDS(current_rds)
            }
        }
        if (is.null(current_pca_df)) {
            cpi <- min_pca_infos[i,,drop=FALSE]
            current_pca_df <- produce_single_pca_df(cpi, txi, metadata)
            if (!is.null(r_objects)) {
                saveRDS(current_pca_df, current_rds)
            }
        }
        pca_df[[cg]][[cgv]][[current_prefix]] <- current_pca_df
        pca_df
    }

    pca_df <- parallel::mclapply(1:nrow(min_pca_infos), create_pca_df,
                                 mc.cores = cores) %>%
        merge_pca_lists
    pca_list <- list()
    for (i in 1:nrow(pca_infos)) {
        current_id <- pca_infos$id_plot[i]
        cg <- pca_infos$group[i]
        cgv <- pca_infos$group_val[i]
        cim <- pca_infos$id_metadata[i]
        cun <- pca_infos$use_normalisation[i]
        cmc <- pca_infos$min_counts[i] %>% as.character
        current_prefix <- paste(cim, cun, cmc, sep = "_")
        cmc <- as.numeric(cmc)

        if (is.na(cg)) {
            cg <- "NA"
            cgv <- "NA"
        }

        current_pca <- pca_df[[cg]][[cgv]][[current_prefix]]
        gg <- produce_single_pca_batch(pca_infos[i,,drop=FALSE], current_pca)

        if (!is.null(r_objects)) {
            output_rds <- paste0(r_objects, "/", current_id, ".rds")
            if (!file.exists(output_rds) | force) {
                saveRDS(gg, output_rds)
            }
        }
        if (!is.null(outdir)) {
            output_pdf <- paste0(outdir, "/", current_id, ".pdf")
            if (!file.exists(output_pdf) | force) {
                pdf(output_pdf)
                print(gg)
                dev.off()
            }
        }
        pca_list[[current_id]] <- gg
    }
    invisible(pca_list)
}

merge_pca_lists <- function(pl) {
    res <- list()
    for (cpd in pl) {
        for (cg in names(cpd)) {
            ccg <- cpd[[cg]]
            if (is.null(res[[cg]]))
                res[[cg]] <- list()
            for (cgv in names(ccg)) {
                ccgv <- ccg[[cgv]]
                if (is.null(res[[cg]][[cgv]]))
                    res[[cg]][[cgv]] <- list()
                for (cim in names(ccgv)) {
                    ccim <- ccgv[[cim]]
                    if (is.null(res[[cg]][[cgv]][[cim]]))
                        res[[cg]][[cgv]][[cim]] <- list()
                    res[[cg]][[cgv]][[cim]] <- ccim
                }
            }
        }
    }
    res
}

complete_pca_infos <- function(pca_infos) {
    # Fill missing
    if (!"group" %in% colnames(pca_infos))
        pca_infos$group <- NA
    if (!"group_val" %in% colnames(pca_infos))
        pca_infos$group_val <- NA
    if (!"use_normalisation" %in% colnames(pca_infos))
        pca_infos$use_normalisation <- "none"
    if (!"min_counts" %in% colnames(pca_infos))
        pca_infos$min_counts <- 5
    if (!"id_metadata" %in% colnames(pca_infos))
        pca_infos$id_metadata <- NA
    if (!"size" %in% colnames(pca_infos))
        pca_infos$size <- 3
    if (!"shape" %in% colnames(pca_infos))
        pca_infos$shape <- NA
    if (!"color" %in% colnames(pca_infos))
        pca_infos$color <- NA
    if (!"title" %in% colnames(pca_infos))
        pca_infos$title <- NA
    if (!"legend.position" %in% colnames(pca_infos))
        pca_infos[["legend.position"]] <- "right"
    if (!"legend.box" %in% colnames(pca_infos))
        pca_infos[["legend.box"]] <- "vertical"
    if (!"show_names" %in% colnames(pca_infos))
        pca_infos[["show_names"]] <- TRUE
    pca_infos
}

validate_pca_infos <- function(pca_infos, metadata, txi) {
    errors <- list()
    for (i in seq_along(pca_infos$id_plot)) {
        current_id  <- pca_infos$id_plot[i]

        # Groups checks
        current_group <- pca_infos$group[i]
        current_group_val <- pca_infos$group_val[i]
        if (!is.na(current_group) & is.null(metadata)) {
            msg <- "group is not NA but there is no metadata table"
            errors[[current_id]] <- c(errors[[current_id]], msg)
        } else {
            if (is.na(current_group) & !is.na(current_group_val)) {
                msg <- "group_val should be NA when group is NA"
                errors[[current_id]] <- c(errors[[current_id]], msg)
            }
            if (!is.na(current_group)) {
                if (!current_group %in% colnames(metadata)) {
                    msg <- "group not found in metadata"
                    errors[[current_id]] <- c(errors[[current_id]], msg)
                } else {
                    if (!current_group_val %in% metadata[[current_group]]) {
                        msg <- "group_val not in group column in metadata"
                        errors[[current_id]] <- c(errors[[current_id]], msg)
                    }
                }
            }
        }

        # Normalisations checks
        current_use_normalisation <- pca_infos$use_normalisation[i]
        if (!is(current_use_normalisation, "character")) {
            msg <- "use_normalisation must be \"none\", \"ruvg\", or \"combat\""
            errors[[current_id]] <- c(errors[[current_id]], msg)
        } else {
            if (!current_use_normalisation %in% c("none", "ruvg", "combat")) {
                msg <- "use_normalisation must be \"none\", \"ruvg\", or \"combat\""
                errors[[current_id]] <- c(errors[[current_id]], msg)
            } else {
                if (current_use_normalisation == "ruvg") {
                    if (!"ruvg_counts" %in% names(txi)) {
                        msg <- "use_normalisation is ruvg but ruvg_counts is missing from txi"
                        errors[[current_id]] <- c(errors[[current_id]], msg)
                    }
                }
                if (current_use_normalisation == "combat") {
                    if (!"combat_counts" %in% names(txi)) {
                        msg <- "use_normalisation is combat but combat_counts is missing from txi"
                        errors[[current_id]] <- c(errors[[current_id]], msg)
                    }
                }
            }
        }

        ## min_counts
        current_min_counts <- pca_infos$min_counts[i]
        if (!is.numeric(current_min_counts)) {
            msg <- "min_counts must be a numeric value"
            errors[[current_id]] <- c(errors[[current_id]], msg)
        } else {
            if (!identical(current_min_counts, round(current_min_counts))) {
                msg <- "min_counts must be an integer"
                errors[[current_id]] <- c(errors[[current_id]], msg)
            }
            if (current_min_counts < 0) {
                msg <- "min_counts must be greater than 0"
                errors[[current_id]] <- c(errors[[current_id]], msg)
            }
        }

        ## id_metadata
        current_id_metadata <- pca_infos$id_metadata[i]
        if (!is.null(metadata)) {
            if (!is(current_id_metadata, "character")) {
                msg <- "id_metadata must be a character value"
                errors[[current_id]] <- c(errors[[current_id]], msg)
            } else {
                if (!current_id_metadata %in% colnames(metadata)) {
                    msg <- "id_metadata not found in metadata column names"
                    errors[[current_id]] <- c(errors[[current_id]], msg)
                }
            }
        } else {
            if (!is.na(current_id_metadata)) {
                msg <- "id_metadata is not NA but metadata table is missing"
                errors[[current_id]] <- c(errors[[current_id]], msg)
            }
        }

        ## size
        current_size <- pca_infos$size[i]

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

        ## Shape
        current_shape <- pca_infos$shape[i]

        if (!is.na(current_shape) & is.null(metadata)) {
            msg <- "shape is not NA but there is no metadata table"
            errors[[current_id]] <- c(errors[[current_id]], msg)
        } else {
            if (!is.na(current_shape)) {
                if (!is(current_shape, "character")) {
                    msg <- "shape must be in character format"
                    errors[[current_id]] <- c(errors[[current_id]], msg)
                } else {
                    if (!current_shape %in% colnames(metadata)) {
                        msg <- "shape value is not in metadata table"
                        errors[[current_id]] <- c(errors[[current_id]], msg)
                    }
                }
            }
        }

        ## Color
        current_color <- pca_infos$color[i]

        if (!is.na(current_color) & is.null(metadata)) {
            msg <- "color is not NA but there is no metadata table"
            errors[[current_id]] <- c(errors[[current_id]], msg)
        } else {
            if (!is.na(current_color)) {
                if (!is(current_color, "character")) {
                    msg <- "color must be in character format"
                    errors[[current_id]] <- c(errors[[current_id]], msg)
                } else {
                    if (!current_color %in% colnames(metadata)) {
                        msg <- "color value is not in metadata table"
                        errors[[current_id]] <- c(errors[[current_id]], msg)
                    }
                }
            }
        }

        ## Title
        current_title <- pca_infos$title[i]
        if (!is.na(current_title)) {
            if (!is(current_title, "character")) {
                msg <- "title must be a character"
                errors[[current_id]] <- c(errors[[current_id]], msg)
            }
        }

        ## Legend
        current_legend_position <- pca_infos[["legend.position"]][i]
        current_legend_box <- pca_infos[["legend.box"]][i]
        if (!current_legend_position %in% c("top", "bottom", "right", "left")) {
            msg <- "Invalid legend.position value (top, bottom, right or left)"
            errors[[current_id]] <- c(errors[[current_id]], msg)
        }
        if (!current_legend_box %in% c("vertical", "horizontal")) {
            msg <- "Invalid legend.box value (vertical or horizontal)"
            errors[[current_id]] <- c(errors[[current_id]], msg)
        }

        ## show_names
        current_show_names <- pca_infos[["show_names"]][i]
        if (!is(current_show_names, "logical")) {
            msg <- "Invalid show_names value (logical)"
            errors[[current_id]] <- c(errors[[current_id]], msg)
        }
    }

    if (length(errors) > 0) {

        message("Errors in batch_PCA:\n")
        for (i in seq_along(errors)) {
            message(paste0(names(errors)[i], ":"))
            for (j in seq_along(errors[[i]])) {
                message(paste0("    ", errors[[i]][j]))
            }
        }
    }
    errors
}

produce_single_pca_df <- function(current_pca_info, txi, metadata) {
    cpi <- current_pca_info
    
    use_normalisation <- cpi$use_normalisation
      
    if (!is.null(metadata)) {

        id_metadata <- validate_metadata(metadata, cpi$id_metadata, txi)
        if ("group" %in% colnames(cpi) & !is.na(cpi$group)) {
            i <- metadata[[cpi$group]] == cpi$group_val
            i[is.na(i)] <- FALSE
            current_samples <- metadata[[cpi$id_metadata]][i]
            txi <- filter_txi(txi, current_samples)
        }
    }
    res_pca <- produce_pca_df(txi = txi,
                              use_normalisation = use_normalisation,
                              min_counts = cpi$min_counts,
                              ncp = 2)
    if (!is.null(metadata)){
      
      
        res_pca$coord <- dplyr::left_join(res_pca$coord, metadata, by = c("sample" = id_metadata))

    }
    res_pca
}

produce_single_pca_batch <- function(current_pca_info, res_pca) {
    cpi <- current_pca_info
    shape <- cpi$shape
    if (is.na(shape)) {
        shape <- NULL
    }
    color <- cpi$color
    if (is.na(color)) {
        color <- NULL
    }
    title <- cpi$title
    if (is.na(title)) {
        title <- NULL
    }
    legend.position <- cpi[["legend.position"]]
    if (is.na(legend.position)) {
        legend.position <- "right"
    }
    legend.box <- cpi[["legend.box"]]
    if (is.na(legend.box)) {
        legend.box <- "vertical"
    }

    plot_pca(res_pca = res_pca,
             size = cpi$size,
             color = color,
             shape = shape,
             title = title,
             graph = FALSE,
             legend.position = legend.position,
             legend.box = legend.box)
}
