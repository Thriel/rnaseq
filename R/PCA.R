#' Produce a PCA plot from txi
#'
#' The PCA is produced using the TPM by default (can also be the ruvg counts;
#' see the \code{ruvg_normalization} param)
#'
#' @param txi The \code{txi} object returned by the \code{import_kallisto}
#' function.
#' @param graph Produce the graph. \code{TRUE} or \code{FALSE}. Default:
#' \code{TRUE}.
#' @param use_normalisation What kind of normalisation should be used instead
#' of TPM? Can be either "none", "ruvg" or "combat". Default: "none"
#' @param min_counts The minimal mean number of count (TPM by default; see also
#' the \code{ruvg_normalization} param) required to keep the feature (gene or
#' transcript) for the PCA.
#' @param metatada The metadata to add to the coordinate table to add color or
#' shape with the results and the \code{plot_pca} function. If there is no
#' metadata available, keep the default value of \code{NULL}. Value can either
#' be a csv file or a \code{data.frame}. Default: \code{NULL}.
#' @param id_metadata The colname to use to join the metadata with the
#' coordinate table. Must contains all the sample names found in the txi. If
#' \code{NULL} and a metadata table is provided, the first column will be used
#' to join the tables. Default: \code{NULL}.
#' @param ncp Number of component to include in the graph. Default: 2.
#'
#' @return Produce the PCA and silently returns a \code{list} with the
#' coordinates and the x and y labels.
#'
#' @examples
#' txi <- get_demo_txi()
#' res_pca <- produce_pca_df(txi)
#'
#' @importFrom magrittr %>%
#' @importFrom FactoMineR PCA
#' @importFrom dplyr mutate
#' @importFrom dplyr group_by
#' @importFrom dplyr summarize
#' @importFrom dplyr filter
#' @importFrom dplyr pull
#' @importFrom tidyr gather
#' @importFrom tidyr spread
#' @importFrom stringr str_extract
#' @importFrom utils tail
#'
#' @export

produce_pca_df <- function(txi, use_normalisation = "none", min_counts = 5,
                           metadata = NULL, id_metadata = NULL, ncp = 2, first_analysis = FALSE) {

    validate_txi(txi)
    stopifnot(use_normalisation %in% c("none", "ruvg", "combat"))
    id_metadata <- validate_metadata(metadata, id_metadata, txi)
    stopifnot(is(ncp, "numeric"))
    stopifnot(identical(ncp, round(ncp)))
    stopifnot(ncp > 1)

    if (use_normalisation == "none") {
        tpm <- as.data.frame(txi$abundance)
    } else if (use_normalisation == "ruvg") {
        stopifnot("ruvg_counts" %in% names(txi))
        stopifnot(is(txi$ruvg_counts, "matrix"))
        tpm <- as.data.frame(txi$ruvg_counts)
    } else if (use_normalisation == "combat") {
        stopifnot("combat_counts" %in% names(txi))
        stopifnot(is(txi$combat_counts, "matrix"))
        tpm <- as.data.frame(txi$combat_counts)
    }

    if (!is.null(txi$dummy)) {
        stopifnot(is(txi$dummy, "character"))
        if (use_normalisation == "none") {
            stopifnot("abundance" %in% txi$dummy)
        } else if (use_normalisation == "ruvg") {
            stopifnot("ruvg_counts" %in% txi$dummy)
        } else if (use_normalisation == "combat") {
            stopifnot("combat_counts" %in% txi$dummy)
        }
    }

    
    if (first_analysis == FALSE) {
        tpm <- tpm %>%
            dplyr::mutate(ensembl_gene = rownames(tpm)) %>%
            tidyr::gather(sample, tpm, -ensembl_gene)
    } else if (first_analysis == TRUE) {
        tpm <- tpm %>%
            dplyr::mutate(ensembl_gene = rownames(tpm)) %>%
            tidyr::gather(sample, tpm, -ensembl_gene) %>%
            dplyr::left_join(metadata %>% select(id, group), by = c("sample" = "id")) %>%
            dplyr::filter(group != "unassigned group") %>%
            dplyr::select(-group)
    }
    min_tpm <- dplyr::group_by(tpm, ensembl_gene) %>%
        dplyr::summarize(tpm = sum(tpm)) %>%
        dplyr::filter(tpm >= min_counts) %>%
        dplyr::pull(ensembl_gene)

    tpm_filter <- dplyr::filter(tpm, ensembl_gene %in% min_tpm)

    df <- tidyr::spread(tpm_filter, sample, tpm) %>%
        as.data.frame
    rownames(df) <- df$ensembl_gene
    m <- dplyr::select(df, -ensembl_gene) %>%
        as.matrix %>%
        t

    pca <- FactoMineR::PCA(m, graph = FALSE)
    coord <- pca$ind$coord
    dims <- stringr::str_extract(colnames(coord), "[0-9]*$") %>%
        as.numeric
    stopifnot(ncp <= tail(dims, 1))
    df <- data.frame(Dim1 = coord[,1])
    for (i in seq(2, ncp)) {
        df[[paste0("Dim", i)]] <- coord[,i]
    }
    df <- df %>%
        dplyr::mutate(sample = rownames(df))

    if (!is.null(metadata)) {
        if (is.null(id_metadata)) {
            id_metadata <- colnames(metadata)[1]
        }
      if (first_analysis == FALSE) {
        df <- left_join(df, metadata, by = c("sample" = id_metadata))
      } else if (first_analysis == TRUE) {
          df <- left_join(df, metadata, by = c("sample" = "id"))
      }
    }

    xlab <- paste0("Dim1 (", pca$eig[1,2] %>% round(2), "%)")
    ylab <- paste0("Dim2 (", pca$eig[2,2] %>% round(2), "%)")

    list(coord = df, xlab = xlab, ylab = ylab)
}

validate_metadata <- function(metadata, id_metadata, txi) {
    if (!is.null(metadata)) {
        stopifnot(is(metadata, "data.frame") | is(metadata, "character"))
        if (is.character(metadata)) {
            stopifnot(file.exists(metadata))
            metadata <- readr::read_csv(metadata, show_col_types = FALSE)
        }
        if (is.null(id_metadata)) {
            id_metadata <- colnames(metadata)[1]
        }
        stopifnot(is(id_metadata, "character"))
        stopifnot(id_metadata %in% colnames(metadata))
        stopifnot(all(colnames(txi$abundance) %in% metadata[[id_metadata]]))
    }
    id_metadata
}

#' Produce a PCA plot from produce_pca_df results
#'
#' @param res_pca The \code{list} returned from the \code{produce_pca} function
#' @param size The size of the points. Default: 3.
#' @param color The name of the column in \code{res_pca$coord} to use to color the
#' points in the PCA. If \code{NULL}, won't add color. Default: \code{NULL}.
#' @param shape The name of the column in \code{res_pca$coord} to use to define
#' the shape of the points in the PCA. If \code{NULL}, won't add change the shape.
#' Default: \code{NULL}.
#' @param show_names Try to add the labels of each point? Will only add the
#' labels if the points are not too closely clustered together in the PCA.
#' Default: \code{TRUE}.
#' @param title Add a title to the graph? If \code{NULL}, no title is added to
#' the graph. Otherwise this param must be a \code{character}. Default:
#' \code{NULL}.
#' @param graph Produce the graph. \code{TRUE} or \code{FALSE}. Default:
#' \code{TRUE}.
#' @param legend.position Value for the \code{legend.position} param from
#' ggplot2. Default: "right".
#' @param legend.box Value for the \code{legend.box} param from ggplot2.
#' Default: "vertical".
#'
#' @return Returns the \code{ggplot} object
#'
#' @examples
#' txi <- get_demo_txi()
#' res_pca <- produce_pca_df(txi)
#' p <- plot_pca(res_pca, graph = FALSE)
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 scale_shape_manual
#' @importFrom ggplot2 ggtitle
#' @importFrom ggrepel geom_text_repel
#'
#' @export
plot_pca <- function(res_pca, size = 3, color = NULL, shape = NULL,
                     show_names = TRUE, title = NULL, graph = TRUE,
                     legend.position = "right",
                     legend.box = "vertical", show_ellipse = FALSE ) {
    
    # Validate the params
    stopifnot(is(res_pca, "list"))
    stopifnot(all(c("coord", "xlab", "ylab") %in% names(res_pca)))
    stopifnot(is(res_pca$coord, "data.frame"))
    stopifnot(nrow(res_pca$coord) > 0)
    stopifnot(all(c("Dim1", "Dim2", "sample") %in% colnames(res_pca$coord)))

    stopifnot(is(res_pca$xlab, "character"))
    stopifnot(is(res_pca$ylab, "character"))

    stopifnot(is(size, "numeric"))
    stopifnot(identical(size, round(size)))
    stopifnot(size > 0)
    if (!is.null(color)) {
        stopifnot(is(color, "character"))
        stopifnot(color %in% colnames(res_pca$coord))
    }
    if (!is.null(shape)) {
        stopifnot(is(shape, "character"))
        stopifnot(shape %in% colnames(res_pca$coord))
    }

    stopifnot(is(show_names, "logical"))
    if (!is.null(title)) {
        stopifnot(is(title, "character"))
    }
    stopifnot(is(graph, "logical"))

    # To avoid 6 shapes limitation
    if (!is.null(shape)) {
        shape_vec <- c(15,16,17,18,19,20,1,2,3,4,5,6,7,8,9,10,11,12,13,14)
        shape_vec <- shape_vec[1:length(unique(res_pca$coord[[shape]]))]
    }

    # Produce the graph
    gg <- ggplot2::ggplot(res_pca$coord, ggplot2::aes(x = Dim1, y = Dim2)) +
        ggplot2::theme_bw() +
        ggplot2::labs(x = res_pca$xlab, y = res_pca$ylab)

    ## Add color and shape
    if (!is.null(color) & !is.null(shape)) {
        gg <- gg +
            ggplot2::geom_point(ggplot2::aes_string(col = color, shape = shape),
                                size = size) +
            ggplot2::scale_shape_manual(values = shape_vec)
    } else if (!is.null(color) & is.null(shape)) {
        gg <- gg +
            ggplot2::geom_point(ggplot2::aes_string(col = color), size = size)
    } else if (is.null(color) & !is.null(shape)) {
        gg <- gg +
            ggplot2::geom_point(ggplot2::aes_string(shape = shape), size = size) +
            ggplot2::scale_shape_manual(values = shape_vec)
    } else {
        gg <- gg + ggplot2::geom_point(size = size)
    }
    
    # Ajouter les ellipses par groupe si demandé
    if (show_ellipse) {
      gg <- gg + ggplot2::stat_ellipse(ggplot2::aes_string(col = color, group = color), 
                                       level = 0.95, linetype = "solid")
    }

  

    # legend
    if (!is.null(color) | !is.null(shape)) {
        gg <- gg + ggplot2::theme(legend.position = legend.position,
                                  legend.box = legend.box)
    }

    # show_names
    if(show_names){
        gg <- gg + ggrepel::geom_text_repel(ggplot2::aes(label = sample),
                                            color = "black", force = 10)
    }

    # title
    if (!is.null(title)) {
        gg <- gg + ggplot2::ggtitle(title)
    }

    # graph
    if (graph) {
        print(gg)
    }

    invisible(gg)
}

#' Deprecated
#'
#' @param txi The \code{txi} object returned by the \code{import_kallisto}
#' function.
#' @param graph Produce the graph. \code{TRUE} or \code{FALSE}. Default:
#' \code{TRUE}.
#' @param use_ruvg Use RUVg normalization? Needs to be pre-computed using the
#' \code{ruvg_normalization} function. Default: \code{FALSE}.
#'
#' @export
produce_pca <- function(txi, graph = TRUE, use_ruvg = FALSE) {
    msg <- paste0("This function is now deprecated, please use the ",
                  "`produce_pca_df` and the `plot_pca` functions to produce ",
                  "your PCA.")
    message(msg)
}
