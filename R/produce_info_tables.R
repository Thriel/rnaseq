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
#' @param header_file Path to the header file that should be added at the top
#' of the report.
#' @param r_objects A \code{character} string corresponding to the directory
#' containing the rds files that will be used to create the report.
#' Default: "r_objects".
#'
#' @return Invisibly returns a named \code{list} of all the tables produced.
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom readr read_csv
#' @importFrom readr write_csv
#' @importFrom readr cols
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
produce_info_tables <-
    function(de_infos,
             design_infos,
             metadata_infos,
             volcano_fc_threshold = c(1.5, 3),
             outdir = NULL,
             header_file = "header.qmd",
             r_objects = "r_objects") {
        # Import input files
        de <- read_csv(de_infos, col_types = cols(.default = "c"))
        stopifnot(c("id_de", "group", "contrast_1", "contrast_2") %in% colnames(de))

        meta <-
            read_csv(metadata_infos, col_types = cols(.default = "c"))
        stopifnot("sample" %in% colnames(meta))

        design <-
            read_csv(design_infos, col_types = cols(.default = "c"))
        stopifnot(c("sample", "group") %in% colnames(design))

        # 1 Prepare metadata file
        mcols_meta <-
            map(de$id_de, \(x) prepare_metadata_mcols(x, de, design)) %>%
            purrr::reduce(left_join, by = "sample")
        stopifnot(all(mcols_meta$sample %in% meta$sample))
        stopifnot(all(meta$sample %in% mcols_meta$sample))
        meta_complete <- left_join(meta, mcols_meta, by = "sample")

        # 2 Prepare volcano file
        volcano <- prepare_volcano_info(de, volcano_fc_threshold)

        # 3 Prepare PCA file
        pca_infos <- prepare_pca_info(meta, de)

        # 4 Report infos
        report_infos <-
            prepare_report_info(pca_infos, de, volcano, meta, header_file, r_objects)

        if (!is.null(outdir)) {
            dir.create(outdir, recursive = TRUE)
            write_csv(
                meta_complete,
                file.path(outdir, "metadata_infos.csv")
            )
            write_csv(volcano, file.path(outdir, "volcano_infos.csv"))
            write_csv(pca_infos, file.path(outdir, "pca_infos.csv"))
            write_csv(report_infos, file.path(outdir, "report_infos.csv"))
        }

        invisible(
            list(
                metadata = meta_complete,
                volcano_infos = volcano,
                pca_infos = pca_infos,
                report_infos = report_infos,
                de_infos = de,
                design_infos = design
            )
        )
    }

#' Parse and check a dataframe, a delimited text file or an Excel file
#'
#' This function is meant to check a table's format against custom
#' specifications and change the types of its columns as desired. A mismatch
#' between expectations and the input table will result in an error.
#'
#' @param path_or_table Either a \code{data.frame} or a single string giving
#' the path to an existing delimited-text or Excel file.
#' @param delim Either \code{NULL} or a single string to use as a text
#' delimiter. Ignored if \code{path_or_table} is not a path to a text file.
#' If \code{NULL}, the delimiter will be guessed from the contents of the file.
#' Default: \code{NULL}
#' @param sheet Either \code{NULL}, a single string or a single integer, giving
#' the name or number of the sheet to read from an Excel file. Ignored if
#' \code{path_or_table} is not a path to an Excel file, or if the sheet is
#' specified using the \code{range} argument. If \code{NULL}, the first sheet is
#' read. Default: \code{NULL}
#' @param range Either \code{NULL} or a single string specifying the rectangle
#' of cells to read in an Excel file. Ignored if \code{path_or_table} is not a
#' path to an Excel file. The rectangle can be specified in any of two formats,
#' "D3:G5" and "R3C4:R5C7" (both examples specify rows 3 to 5 and
#' columns 4 to 7). The name of the sheet to read, e.g. "mysheet", can be
#' specified with "mysheet!D3:G5" or "mysheet!R3C4:R5C7". A sheet name specified
#' in this way overrides the \code{sheet} argument. If \code{range} is
#' \code{NULL}, the rectangle is guessed from the contents of the file.
#' Default: \code{NULL}
#' @param required_col_names Either \code{NULL} or a character vector of column
#' names that should appear in the table. The function will produce an error if
#' any of these column names is missing from the table. Default: \code{NULL}
#' @param required_col_patterns Either \code{NULL} or a character vector. If it
#' is a named vector, all its names should appear as column names in the table.
#' If it has no names, it should be the same length as \code{required_col_names}
#' and its elements are then assumed to correspond to the elements of
#' \code{required_col_names} in the same order. The elements of
#' \code{required_col_patterns} should be valid (ICU) regular expressions. If an
#' input column has a corresponding pattern, every entry in this column must
#' contain a match of the pattern.
#' Default: \code{NULL}
#' @param required_col_types A \code{cols()} specification giving the types to
#' which input columns should be converted, e.g.:
#' \code{cols(age = "i", weight = "d", .default = "c")}. Every column mentioned
#' in the specification should appear in the input table, and its entries should
#' be convertible to the specified type. If the specification omits a column and
#' has no \code{.default} argument, the type of the column will be guessed from
#' its contents. Default: \code{cols()}
#'
#'
#' @return A tibble resulting from the conversion of the input table using
#' the \code{required_col_types} specification.
#'
#' @examples
#' library("tibble")
#' library("readr")
#' library("magrittr")
#' input <- tibble(x = c("A_B", "C_D"), y = c("1", "2"))
#' expected <- tibble(x = c("A_B", "C_D"), y = as.integer(c(1, 2)))
#' output <-
#'     input %>%
#'     parse_table(
#'         required_col_names = c("x", "y"),
#'         required_col_patterns = c(x = "[A-Z]_[A-Z]"),
#'         required_col_types = cols(y = "i", .default = "c")
#'     )
#' stopifnot(identical(output, expected))
#'
#' @importFrom readr cols
#' @importFrom readr read_delim
#' @importFrom readr type_convert
#' @importFrom dplyr pull
#' @importFrom stringr str_detect
#' @importFrom tibble as_tibble
#' @importFrom tibble tibble
#' @importFrom readxl excel_format
#' @importFrom readxl read_excel
#'
#' @export
parse_table <- function(path_or_table,
                        delim = NULL,
                        sheet = NULL,
                        range = NULL,
                        required_col_names = NULL,
                        required_col_patterns = NULL,
                        required_col_types = cols()) {
    # Check required_col_names
    if (!class(required_col_names) %in% c("NULL", "character")) {
        stop("required_col_names must be a character vector.")
    }

    # Check required_col_patterns
    if (is.null(required_col_patterns)) {
        patterns <- NULL
    } else {
        if (!is.character(required_col_patterns)) {
            stop("required_col_patterns must be a character vector.")
        }
        if (is.null(names(required_col_patterns))) {
            if (length(required_col_patterns) != length(required_col_names)) {
                stop(
                    paste0(
                        "If required_col_patterns is not a named vector, ",
                        "it must have the same length as required_col_names"
                    )
                )
            }
            patterns <-
                setNames(required_col_patterns, required_col_names)
        } else {
            patterns <- required_col_patterns
        }
    }

    # Check required_col_types
    if (!is.null(required_col_types)) {
        if (!is(required_col_types, "col_spec")) {
            stop("required_col_types must be a cols() specification.")
        }
    }

    # If needed, create tibble from delimited text file or Excel file
    if (is(path_or_table, "data.frame")) {
        table <-
            as_tibble(path_or_table) %>%
            mutate(across(everything(), as.character))
    } else if (is.character(path_or_table) &
        length(path_or_table) == 1) {
        path <- path_or_table

        if (!file.exists(path)) {
            stop(paste0(
                "No file at path ",
                path,
                " from current directory ",
                getwd()
            ))
        }

        if (is.na(excel_format(path))) {
            table <- read_delim(path, delim, col_types = cols(.default = "c"))
        } else {
            table <- read_excel(path, sheet, range, col_types = "text")
        }
    } else {
        stop("path_or_table must be a single string or a dataframe.")
    }

    all_col_names <-
        unique(c(
            required_col_names,
            names(patterns),
            names(required_col_types$cols)
        ))

    missing_col_names <- setdiff(all_col_names, colnames(table))

    if (length(missing_col_names) > 0) {
        stop(
            paste0(
                "The following required column names are missing from the table: ",
                paste(missing_col_names, collapse = ", ")
            )
        )
    }

    # Check that table entries match required patterns
    for (col_name in names(patterns)) {
        values <- table %>% pull(col_name)
        pattern <- patterns[col_name]
        matches <- values %>% str_detect(pattern)
        bad_values <- values[!matches]
        if (length(bad_values) > 0) {
            stop(
                paste0(
                    "Some values in column ",
                    col_name,
                    " (e.g. ",
                    paste(head(bad_values, 3), collapse = ", "),
                    ") do not match pattern \"",
                    pattern,
                    "\""
                )
            )
        }
    }

    # Convert columns to required types
    for (col_name in names(table)) {
        error_func <- function(e) {
            x <- tibble()
            x[[col_name]] <- character(0)
            expected_type <- class(
                type_convert(x, required_col_types)[[col_name]]
            )
            stop(paste0(
                "Could not convert column ",
                col_name,
                " to type ",
                expected_type,
                ":\n\t",
                e
            ))
        }
        tryCatch(
            table[, col_name] <-
                table[, col_name] %>%
                type_convert(required_col_types),
            error = error_func,
            warning = error_func
        )
    }

    return(table)
}

prepare_metadata_mcols <- function(current_id_de, de, design) {
    print(current_id_de)
    current_de <- filter(de, id_de == current_id_de)
    stopifnot(nrow(current_de) == 1)

    current_grp_col <-
        colnames(design)[colnames(design) == current_de$group]
    stopifnot(length(current_grp_col) == 1)

    current_grp1 <- current_de$contrast_1
    current_grp2 <- current_de$contrast_2

    current_samples_grp1 <-
        design$sample[pull(design, one_of(current_grp_col)) %in% current_grp1]
    current_samples_grp2 <-
        design$sample[pull(design, one_of(current_grp_col)) %in% current_grp2]
    current_samples <- c(current_samples_grp1, current_samples_grp2)

    res <- dplyr::select(design, sample) %>%
        mutate(tmp = if_else(sample %in% current_samples, "yes", "no"))
    colnames(res)[2] <- current_id_de
    res
}

prepare_volcano_info <- function(de_infos, fc_threshold) {
    id_plot <- paste0(
        rep(de_infos$id_de, each = length(fc_threshold)),
        "_FC_",
        rep(fc_threshold, length(de_infos$id_de))
    ) %>%
        str_replace_all("\\.", "_")
    tibble(
        id_plot = id_plot,
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
        size = 3
    )
}

prepare_pca_info <- function(metadata, de_infos) {
    col <- dplyr::select(metadata, -sample) %>% colnames()

    ## Global
    titles <-
        paste0("Global PCA ", c(" ", rep("- ", length(col))), c("", col)) %>%
        str_replace_all("_", " ") %>%
        str_replace_all(" - $", "") %>%
        str_to_title() %>%
        str_replace("Pca", "PCA")
    ids_plot <- paste0("global_pca", c("", paste0("_", col)))
    global_pcas <- tibble(
        id_plot = ids_plot,
        group = NA,
        group_val = NA,
        use_normalization = "none",
        min_counts = 5,
        id_metadata = "sample",
        color = c(NA, col),
        title = titles
    )
    ## DE
    ids_de <- de_infos$id_de
    ids_pca <- paste0(
        "pca_",
        rep(ids_de, length(col) + 1),
        "_",
        rep(c("", col), each = length(ids_de))
    ) %>%
        str_replace_all("_$", "")
    titles <-
        paste0(
            rep(str_replace(ids_de, " vs ", " versus "), length(col) + 1),
            " - ",
            rep(c("", col), each = length(ids_de))
        ) %>%
        str_replace_all(" - $", "")
    de_pcas <- tibble(
        id_plot = ids_pca,
        group = rep(ids_de, length(col) + 1),
        group_val = "yes",
        use_normalization = "none",
        min_counts = 5,
        id_metadata = "sample",
        color = rep(c(NA, col), each = length(ids_de)),
        title = titles
    )

    bind_rows(global_pcas, de_pcas)
}

prepare_report_info <-
    function(pca_infos,
             de_infos,
             volcano_infos,
             metadata,
             header_file,
             r_objects) {
        ## Header
        header_report <- tibble(add = "file", value = header_file)

        ## Global PCA
        global_pca <-
            filter(pca_infos, str_detect(id_plot, "^global_pca"))

        global_pca_add <- c("text", rep("plot", nrow(global_pca)))
        global_pca_values <- c(
            "## Global PCA",
            paste0(r_objects, "/", global_pca$id_plot, ".rds")
        )
        global_pca_report <-
            tibble(add = global_pca_add, value = global_pca_values)

        ## DE results
        col <- dplyr::select(metadata, -sample) %>% colnames()

        de_results_header <-
            tibble(add = "text", value = "## Differential expression analysis results")

        produce_de_res_report <- function(current_id_de) {
            current_header_report <-
                tibble(
                    add = "text",
                    value = paste0("### ", current_id_de)
                )

            ### PCA
            current_pca_values <-
                paste0(
                    r_objects,
                    "/pca_",
                    current_id_de,
                    "_",
                    c("", col),
                    ".rds"
                ) %>%
                str_replace("_\\.rds$", "\\.rds")
            current_pca_report <-
                tibble(
                    add = c("text", rep("plot", length(col) + 1)),
                    value = c("#### PCA", current_pca_values)
                )

            ### Volcano
            current_volcano <-
                filter(volcano_infos, id_de == current_id_de)
            current_volcano_add <-
                rep("plot", nrow(current_volcano))
            current_volcano_values <-
                paste0(r_objects, "/", current_volcano$id_plot, ".rds")
            current_volcano_report <-
                tibble(
                    add = c("text", current_volcano_add),
                    value = c("#### Volcano plots", current_volcano_values)
                )

            ### Merge everything
            bind_rows(
                current_header_report,
                current_pca_report,
                current_volcano_report
            )
        }
        de_res_report <-
            map_dfr(de_infos$id_de, produce_de_res_report)

        bind_rows(header_report, global_pca_report, de_res_report)
    }
