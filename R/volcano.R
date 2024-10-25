#' produce a volcano plot from DESeq2 results
#'
#' @param de_res the \code{data.frame} object returned by the
#' \code{DESeq2::results} function.
#' @param fc_threshold The threshold of FC to be considered as significant.
#' Default: 0.05
#' @param p_threshold The threshold of p stat to be considered as significant.
#' Default: 1.5
#' @param y_axis Statistical results to show in volcano plot. "pvalue" or "padj".
#' Defauld: "padj"
#' @param show_signif_counts show the number of up- and down-regulated genes?
#' @param show_signif_lines show lines at the threshold for significance? "none",
#' "vertical","horizontal" or "both". Default: "vertical"
#' @param show_signif_color Show color for significant genes? Default: TRUE
#' @param col_up Color of the up-regulated genes. Default: "#E73426"
#' @param col_down Color of the down-regulated genes. Default: "#0020F5"
#' @param size The size of the points. Default: 3.
#' @param graph produce the graph. \code{true} or \code{false}. default:
#' \code{true}.
#' @param title The title of the graph. If \code{NA}, no title will be
#' displayed. Default: \code{NA}.
#'
#' @return produce the volcano plot and silently returns the \code{ggplot}
#' object and the data.frame used. # TODO: data.frame w anno, df signif
#'
#' @examples
#' txi <- get_demo_txi()
#' design <- get_demo_design()
#' dds <- deseq2_analysis(txi, design, ~ group)
#' de <- DESeq2::results(dds, contrast = c("group", "A", "B"))
#' volcano <- produce_volcano(de, graph = FALSE)
#'
#' @importFrom stringr str_detect
#' @importFrom dplyr mutate
#' @importFrom dplyr if_else
#' @importFrom dplyr filter
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_vline
#' @importFrom ggplot2 geom_hline
#' @importFrom ggplot2 annotate
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 scale_colour_identity
#' @importFrom ggplot2 theme_minimal
#' @importFrom ggplot2 annotate
#' @importFrom ggplot2 ylab
#' @importFrom utils head
#' @importFrom dplyr filter
#' @importFrom dplyr mutate
#' @importFrom dplyr arrange
#' @importFrom dplyr pull
#' @importFrom dplyr if_else
#' @importFrom stringr str_detect
#'
#' @export
    produce_volcano <- function(de_res, fc_threshold = 3, p_threshold = 0.05,
                            y_axis = "padj", show_signif_counts = TRUE,
                            show_signif_lines = "vertical",
                            show_signif_color = TRUE, col_up = "#E73426",
                            col_down = "#0020F5", size = 3, graph = TRUE,
                            title = NA, outdir = NULL) {
    stopifnot(is.numeric(fc_threshold))
    stopifnot(fc_threshold > 0)
    stopifnot(is.numeric(p_threshold))
    stopifnot(p_threshold >= 0 & p_threshold <= 1)
    stopifnot(is(show_signif_counts, "logical"))
    stopifnot(is(show_signif_lines, "character"))
    stopifnot(is(title, "character") | is.na(title))
    expected_values <- c("none", "both", "vertical", "horizontal")
    stopifnot(show_signif_lines %in% expected_values)
    stopifnot(is(show_signif_color, "logical"))
    if (show_signif_color) {
        stopifnot(is(col_up, "character"))
        stopifnot(is_color(col_up))
        stopifnot(is(col_down, "character"))
        stopifnot(is_color(col_down))
    }
    stopifnot(is(graph, "logical"))
    stopifnot(is(y_axis, "character"))
    stopifnot(y_axis %in% c("pvalue", "padj"))

    if (is(de_res, "DESeqResults")) {
        de_res <- as.data.frame(de_res)
    }

    # Rename qV to padj
    i <- stringr::str_detect(colnames(de_res), "qV")
    stopifnot(sum(i) %in% c(0,1))
    if (sum(i) == 1) {
        colnames(de_res)[i] <- "padj"
    }
    # Rename pV to pvalue
    i <- stringr::str_detect(colnames(de_res), "pV")
    stopifnot(sum(i) %in% c(0,1))
    if (sum(i) == 1) {
        colnames(de_res)[i] <- "pvalue"
    }
    stopifnot(y_axis %in% colnames(de_res))

    if (show_signif_color) {
        red <- col_up
        blue <- col_down
        grey <- "#7C7C7C"
    } else {
        grey <- "#7C7C7C"
        red <- "#E73426"
        blue <- "#0020F5"
    }

    de_res <- dplyr::mutate(de_res,
                            padj = as.numeric(padj),
                            pvalue = as.numeric(pvalue),
                            log2FoldChange = as.numeric(log2FoldChange))

    # Remove NA
    de_res <- de_res[!is.na(de_res[[y_axis]]),]

    # Remove inf -log10(padj), i.e.: padj == 0
    min_y_axis <- de_res[de_res[[y_axis]] != 0, y_axis, drop = TRUE] %>% min
    i <- de_res[[y_axis]] == 0
    de_res[[y_axis]][i] <- min_y_axis

    # Add color
    de_res$y_axis <- de_res[[y_axis]]
    de_res <- dplyr::mutate(de_res, color = grey) %>%
        dplyr::mutate(color = dplyr::if_else(log2FoldChange <= -log2(fc_threshold) &
                                             y_axis <= p_threshold, blue, color)) %>%

        dplyr::mutate(color = dplyr::if_else(log2FoldChange >= log2(fc_threshold) &
                                             y_axis <= p_threshold, red, color))

    count_blue <- dplyr::filter(de_res, color == blue) %>% nrow
    count_red <- dplyr::filter(de_res, color == red) %>% nrow
    lbl <- c(count_blue, count_red) %>% as.character
    count_y <- round(max(-log10(de_res$padj), na.rm = TRUE))
    count_y <- count_y * 0.925
    min_x <- round(min(de_res$log2FoldChange, na.rm = TRUE))
    min_x <- min_x * 0.925
    max_x <- round(max(de_res$log2FoldChange, na.rm = TRUE))
    max_x <- max_x * 0.925

    if (!show_signif_color) {
        de_res <- mutate(de_res, color = grey)
        p <- ggplot2::ggplot(de_res, ggplot2::aes(x = log2FoldChange,
                                                  y = -log10(y_axis))) +
            ggplot2::geom_point(size = size, alpha = 0.8)
    } else {
        p <- ggplot2::ggplot(de_res, ggplot2::aes(x = log2FoldChange,
                                                  y = -log10(y_axis),
                                                  color = color)) +
            ggplot2::geom_point(size = size, alpha = 0.8) +
            ggplot2::scale_colour_identity()
    }
    if (show_signif_lines %in% c("both", "vertical")) {
        p <- p + ggplot2::geom_vline(xintercept = c(-log2(fc_threshold),
                                                    log2(fc_threshold)),
                                     linetype = "dashed")
    }
    if (show_signif_lines %in% c("both", "horizontal")) {
        p <- p + ggplot2::geom_hline(yintercept = -log10(p_threshold),
                                     linetype = "dashed")
    }
    if (show_signif_counts) {
        if (show_signif_color) {
            p <- p + ggplot2::annotate("text",
                                       x = c(min_x, max_x),
                                       y = count_y,
                                       label = lbl,
                                       size = 8,
                                       fontface = 2,
                                       color = c(blue, red))
        } else {
            p <- p + ggplot2::annotate("text",
                                       x = c(min_x, max_x),
                                       y = count_y,
                                       label = lbl,
                                       size = 8,
                                       fontface = 2,
                                       color = c(grey, grey))
        }
    }
    if(!is.na(title)){
        p <- p + ggplot2::ggtitle(title)
    }
    p <- p + ggplot2::theme_minimal() +
        ggplot2::ylab(paste0("-log10(", y_axis, ")"))

    if (isTRUE(graph)) {
        print(p)
    }
    invisible(list(p = p, df = de_res))
}

    
is_color <- function(x) {
    res <- try(col2rgb(x),silent=TRUE)
    return(!"try-error"%in%class(res))
}
    