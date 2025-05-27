#' Get Results for Multiple Cell Types
#'
#' This function retrieves results for multiple cell types from an RCTD object.
#'
#' @param rctd An RCTD object.
#' @param confident_only Logical, whether to include only confident results. Default is TRUE.
#' @param order_by_abundance Logical, whether to order results by abundance. Default is FALSE.
#' @return A data frame with cell types and their corresponding weights.
#' @export
rctd_get_results_multi <- function(rctd, confident_only = TRUE, order_by_abundance = FALSE) {
    stopifnot("RCTD" %in% class(rctd))

    spots <- colnames(rctd@spatialRNA@counts)

    max_multi_types <- rctd@results |>
        lapply(\(r) length(r$cell_type_list)) |>
        unlist() |>
        max()

    lapply(rctd@results, \(r) {
        cell_types <- weights <- rep(NA, max_multi_types)

        if (confident_only) {
            ct_conf <- names(r$conf_list)[r$conf_list == TRUE]
            cell_types[seq_along(ct_conf)] <- ct_conf
            weights[seq_along(ct_conf)] <- r$sub_weights[ct_conf]
        } else {
            cell_types[seq_along(r$cell_type_list)] <- r$cell_type_list
            weights[seq_along(r$sub_weights)] <- r$sub_weights[r$cell_type_list]
        }

        if (order_by_abundance) {
            o <- order(weights, decreasing = TRUE)
            cell_types <- cell_types[o]
            weights <- weights[o]
        }

        list(cell_types, weights)
    }) |>
        (\(.) {
            cbind.data.frame(
                do.call(rbind, lapply(., "[[", 1)),
                do.call(rbind, lapply(., "[[", 2))
            )
        })() |>
        `rownames<-`(spots) |>
        `colnames<-`(
            c(
                paste0("cell_type_", seq_len(max_multi_types)),
                paste0("weight_", seq_len(max_multi_types))
            )
        )
}

#' @title Extract RCTD Results for doublet mode 'full'
#'
#' @description
#' Extracts and combines cell type predictions and weights from an RCTD object.
#'
#' @param rctd An RCTD object containing the results of the deconvolution.
#' @param n_discrete The number of top cell types to return for each spot.
#'   Defaults to 3.
#'
#' @return A data frame containing the top `n_discrete` predicted cell types
#'   and the corresponding weights of all cell types for each spot.
#'   The data frame has columns `cell_type_1`, `cell_type_2`, ..., `cell_type_n_discrete`
#'   containing the names of the top cell types, followed by columns containing the weights
#'   for each cell type.
#'
#' @examples
#' \dontrun{
#' # Assuming 'rctd_object' is an RCTD object
#' results <- rctd_get_results_full(rctd_object, n_discrete = 5)
#' head(results)
#' }
rctd_get_results_full <- function(rctd, n_discrete = 3) {
    stopifnot("RCTD" %in% class(rctd))

    # ensure weights sum to 1
    weights <- rctd@results$weights |>
        spacexr::normalize_weights() |>
        as.matrix()

    # colnames(weights)[max.col(weights)]
    cell_types <- apply(weights, 1, \(r) {
        o <- order(r, decreasing = TRUE)
        colnames(weights[, o])[seq_len(n_discrete)]
    }) |>
        t() |>
        as.data.frame() |>
        `colnames<-`(c(paste0("cell_type_", seq_len(n_discrete))))

    top_n_weights_names <- c(paste0("weight_", seq_len(n_discrete)))
    top_n_cell_types_weights <- map_df(rownames(cell_types), \(bc) {
        weights[bc, as.character(cell_types[bc, ])] |>
            `names<-`(top_n_weights_names)
    })

    cbind.data.frame(
        cell_types,
        top_n_cell_types_weights,
        weights
    )
}

#' Fix Names
#'
#' This function fixes names so that RCTD can handle them by replacing certain special
#' characters with underscores.
#'
#' @param names A character vector of names to be fixed.
#' @return A character vector with fixed names.
#' @example fix_names(c("T cells (T.CD8.96H)", "DC (DC.103-11B+F4-80LO.KD)"))
#' @export
fix_names <- function(names) {
    names |>
        # remove a subset of special characters
        (\(.) gsub("[\\/+.*()&%]", "_", .))() |>
        # remove excessive gaps
        (\(.) gsub("[_* *]+", "_", .))() |>
        # remove trailing underscore
        (\(.) gsub("_$", "", .))()
}

#' Get Weights for Multiple Cell Types
#'
#' This function retrieves weights for multiple cell types from a results data frame.
#'
#' @param results_df A data frame containing results.
#' @param cell_types A character vector of cell types.
#' @return A data frame with weights for each cell type.
#' @export
rctd_get_weights_multi3 <- function(results_df, cell_types) {
    rownames(results_df) |>
        sapply(\(spot) {
            row_vals <- results_df[spot, ]
            cell_types |>
                sapply(\(cell_type) {
                    ifelse(
                        cell_type %in% row_vals$cell_type_1,
                        row_vals$weight_1,
                        ifelse(
                            cell_type %in% row_vals$cell_type_2,
                            row_vals$weight_2,
                            ifelse(
                                cell_type %in% row_vals$cell_type_3,
                                row_vals$weight_3,
                                0.00
                            )
                        )
                    )
                })
        }) |>
        t() |>
        as.data.frame()
}

#' Create a Spatial Plot of RCTD Weights
#'
#' This function plots the corresponding weights for a given cell type variable in a Seurat object.
#'
#' @param so A Seurat object.
#' @param var A character string specifying the cell type to plot.
#' @param save Logical, whether to save the plot. Default is TRUE.
#' @param save_path A character string specifying the path to save the plot. Default is ".".
#' @param return Logical, whether to return the plot object. Default is FALSE.
#' @param ... Additional arguments passed to Seurat::SpatialFeaturePlot.
#' @return A ggplot object if return is TRUE.
#' @export
rctd_plot_weights <- function(
    so,
    var,
    name = NULL,
    save = TRUE,
    save_path = ".",
    return = FALSE,
    plot_aspect_ratio = 16 / 9,
    ...
) {
    stopifnot("Seurat" %in% class(so))
    stopifnot(var %in% colnames(so[[]]))

    sfp_defaults <- list(
        alpha = 1,
        ncol = 4,
        stroke = 0,
        image.alpha = 0.25
    )
    sfp_args <- dots_list(..., !!!sfp_defaults, .homonyms = "first")

    var_name_short <- paste(tail(strsplit(var, "_")[[1]], 3), collapse = "_")
    var_weight <- sub("cell_type", "weight", var)
    var_weight_short <- paste(tail(strsplit(var_weight, "_")[[1]], 2), collapse = "_")

    # so[[]][[var_weight]] <- so[[]][[var_weight]] |> (\(.) replace(., is.na(.), -1))()
    na_idx <- is.na(so[[]][[var_weight]])
    message(sprintf("Removing %s NAs", sum(na_idx)))

    p <- do.call(
        SpatialFeaturePlot,
        c(
            list(subset(so, cells = rownames(so[[]])[which(!na_idx)]), features = var_weight),
            sfp_args
        )
    ) +
        plot_layout(guides = "collect") &
        scale_fill_gradientn(
            (if (is.null(name)) var_weight_short else name),
            limits = c(0, 1),
            colours = SpatialColors(n = 100)
        ) &
        theme_legend_tight(legend_size_npc = 0.025, legend_width_npc = 0.0375) &
        theme(
            plot.margin = margin(0, 0, 0, 0),
            legend.position = "bottom",
            legend.title = element_text(margin = margin(15, 15, 15, 15, unit = "pt"))
        )

    if (save) {
        ggsave(
            sprintf("Spatial_%s.png", gsub("\\s", "_", var_name_short)),
            path = save_path,
            width = 12,
            height = 12 / plot_aspect_ratio,
            dpi = 300,
            device = ragg::agg_png
        )
    }

    if (return) {
        p
    }
}

#' Plot Cell Types One by One
#'
#' This function plots cell types one by one for a given variable in a Seurat object.
#'
#' @param so A Seurat object.
#' @param var A character string specifying the variable to plot.
#' @param cell_types A character vector of cell types to plot.
#' @param save Logical, whether to save the individual plots. Default is TRUE.
#' @param save_path A character string specifying the path to save the plot. Default is ".".
#' @param return Logical, whether to return the plot object. Default is FALSE.
#' @param ... Additional arguments passed to the plotting functions.
#' @return A list of ggplot objects if return is TRUE.
#' @export
rctd_plot_cell_types_1by1 <- function(
    so,
    var,
    cell_types,
    groups_from = "group",
    save = TRUE,
    save_path = ".",
    return = FALSE,
    plot_aspect_ratio = 16 / 9,
    ...
) {
    stopifnot("Seurat" %in% class(so))
    stopifnot(all(c(var, groups_from) %in% colnames(so[[]])))

    sdp_defaults <- list(
        label = FALSE,
        ncol = 4,
        stroke = 0,
        image.alpha = 0.25
    )
    sdp_args <- dots_list(..., !!!sdp_defaults, .homonyms = "first")

    so$cluster_group <- paste(so[[]][[var]], so[[]][[groups_from]], sep = "_")
    so$cluster_group <- factor(
        so[[]]$cluster_group,
        levels = unique(stringr::str_sort(so[[]]$cluster_group, numeric = TRUE))
    )
    Idents(so) <- "cluster_group"
    pl <- setNames(cell_types, cell_types) |>
        map(\(cell_type) {
            cluster_groups <- paste(cell_type, unique(so[[]][[groups_from]]), sep = "_")
            colors <- if (length(cluster_groups) > 2) {
                grDevices::hcl.colors(length(cluster_groups), "Red-Blue")
            } else {
                c("red", "blue")
            }
            color_mapping <- c("Unselected" = "grey70", setNames(colors, cluster_groups))

            do.call(
                possibly(SpatialDimPlot),
                c(
                    list(
                        so,
                        cells.highlight = possibly(CellsByIdentities)(so, idents = cluster_groups),
                        cols.highlight = color_mapping
                    ),
                    sdp_args
                )
            ) +
                plot_layout(guides = "collect") &
                theme(
                    plot.margin = margin(0, 0, 0, 0),
                    legend.position = "bottom",
                    legend.key = element_rect(fill = "transparent")
                ) &
                scale_fill_manual(
                    cell_type,
                    values = color_mapping,
                    limits = names(color_mapping),
                    drop = FALSE,
                    guide = guide_legend(override.aes = list(size = 3))
                )

            if (save) {
                ggsave(
                    sprintf("Spatial_%s.png", gsub("\\s", "_", cell_type)),
                    path = save_path,
                    width = 12,
                    height = 12 / plot_aspect_ratio,
                    dpi = 300,
                    device = ragg::agg_png
                )
            }
        })

    if (return) {
        pl
    }
}

#' Create a Spatial Plot of RCTD Cell Type Assignments
#'
#' This function plots cell types for a given variable in a Seurat object.
#'
#' @param so A Seurat object.
#' @param var A character string specifying the variable to plot.
#' @param colors A character vector of colors for the cell types.
#' @param name A character string specifying the name of the plot. Default is NULL.
#' @param save Logical, whether to save or return the plot. Default is TRUE.
#' @param save_path A character string specifying the path to save the plot. Default is ".".
#' @param plot_aspect_ratio A number specifying the aspect ratio of the saved image.
#' @param ... Additional arguments passed to the plotting functions.
#' @return A ggplot object if return is TRUE.
#' @export
rctd_plot_cell_types <- function(
    so,
    var,
    colors,
    name = NULL,
    save = TRUE,
    save_path = ".",
    plot_aspect_ratio = 16 / 9,
    ...
) {
    stopifnot("Seurat" %in% class(so))
    stopifnot(var %in% colnames(so[[]]))

    sdp_defaults <- list(
        label = FALSE,
        ncol = 4,
        stroke = 0,
        image.alpha = 0.25
    )
    sdp_args <- dots_list(..., !!!sdp_defaults, .homonyms = "first")

    var_name_short <- paste(tail(strsplit(var, "_")[[1]], 3), collapse = "_")

    # get unique cell types ordered by frequency
    cell_types <- names(sort(table(so[[]][[var]]), decreasing = TRUE))
    n_cell_types <- length(cell_types)

    stopifnot(length(colors) >= n_cell_types)

    color_mapping <- if (is.null(names(colors))) {
        cell_types |>
            purrr::modify_if(is.na, `<-`, "NA") |>
            setNames(object = colors[seq_len(n_cell_types)])
    } else {
        colors
    }

    # Spatial plot
    p <- do.call(SpatialDimPlot, c(list(so, group.by = var), sdp_args)) +
        plot_layout(guides = "collect") &
        theme(
            legend.position = "bottom",
            legend.key = element_rect(fill = "transparent", colour = "transparent"),
            legend.key.size = unit(0.60, "lines"),
            legend.margin = margin(1, 1, 1, 1),
            legend.box.margin = margin(1, 1, 1, 1),
            plot.margin = margin(0, 0, 0, 0)
        ) &
        scale_fill_manual(
            name = {
                if (!is.null(name)) name else var_name_short
            },
            values = color_mapping,
            breaks = cell_types,
            # never drop to enable collecting guides with patchwork
            drop = FALSE,
            guide = guide_legend(
                nrow = legend_rows(n_cell_types, 6),
                byrow = TRUE,
                override.aes = list(size = 3)
            )
        ) &
        guides(color = "none")

    if (save) {
        ggsave(
            sprintf("Spatial_%s.png", var_name_short),
            path = save_path,
            width = 12,
            height = 12 / plot_aspect_ratio,
            dpi = 300,
            device = ragg::agg_png
        )
    }

    p
}


#' Create a UMAP Plot from a RCTD Variable
#'
#' This function plots a UMAP for a given variable in a Seurat object.
#'
#' @param so A Seurat object.
#' @param var A character string specifying the variable to plot.
#' @param colors A character vector of colors for the cell types.
#' @param name A character string specifying the name of the plot. Default is NULL.
#' @param save Logical, whether to save the plot. Default is TRUE.
#' @param save_path A character string specifying the path to save the plot. Default is ".".
#' @param return Logical, whether to return the plot object. Default is FALSE.
#' @param ... Additional arguments passed to the plotting functions.
#' @return A ggplot object if return is TRUE.
#' @export
rctd_plot_umap <- function(
    so,
    var,
    colors,
    name = NULL,
    save = TRUE,
    save_path = ".",
    return = FALSE,
    plot_aspect_ratio = 5 / 4,
    ...
) {
    stopifnot("Seurat" %in% class(so))
    stopifnot(var %in% colnames(so[[]]))

    var_name_short <- paste(tail(strsplit(var, "_")[[1]], 3), collapse = "_")

    # get unique cell types ordered by frequency
    cell_types <- names(sort(table(so[[]][[var]]), decreasing = TRUE))
    n_cell_types <- length(cell_types)

    stopifnot(length(colors) >= n_cell_types)

    color_mapping <- if (is.null(names(colors))) {
        cell_types |>
            purrr::modify_if(is.na, `<-`, "NA") |>
            setNames(object = colors[seq_len(n_cell_types)])
    } else {
        colors
    }

    # UMAP
    p <- DimPlot(
        so,
        reduction = "umap",
        group.by = var,
        ...
    ) &
        theme(legend.position = "bottom") &
        theme_legend_tight() &
        scale_color_manual(
            name = {
                if (!is.null(name)) name else var_name_short
            },
            values = color_mapping,
            drop = FALSE,
            guide = guide_legend(
                nrow = legend_rows(n_cell_types, 7),
                byrow = TRUE,
                override.aes = list(size = 3)
            )
        )

    if (save) {
        ggsave(
            sprintf("UMAP_%s.png", var_name_short),
            path = save_path,
            width = 10,
            height = 10 / plot_aspect_ratio,
            dpi = 300,
            device = ragg::agg_png
        )
    }

    if (return) {
        p
    }
}

#' Get Manual Color Mapping for Cell Types
#'
#' This function generates a manual color mapping for cell types based on the provided reference and annotation.
#' It uses predefined color palettes to ensure consistent and distinguishable colors for different cell types.
#'
#' @param reference A string specifying the reference dataset name. Currently supports "NBAtlas".
#' @param annotation A string specifying the annotation type within the reference dataset.
#'  Supports "Cell_type" and "annot_immunezoom_v5".
#' @return A named vector where the names are cell types and the values are their corresponding colors.
#' @examples
#' color_mapping <- rctd_get_manual_color_mapping("NBAtlas", "Cell_type")
#' print(color_mapping)
#' @export
rctd_get_manual_color_mapping <- function(reference, annotation, fix_names = TRUE) {
    nba_base_mapping <- setNames(
        as.character(
            # skip the dark grey which could be mistaken for the default NA color
            paletteer::paletteer_d("ggsci::category20_d3", n = 12 + 1)[-8]
        ),
        c(
            "Neuroendocrine",
            "Myeloid",
            "Fibroblast",
            "Schwann",
            "T cell",
            "NK cell",
            "Endothelial",
            "RBCs",
            "B cell",
            "Plasma",
            "pDC",
            "Stromal other"
        )
    )
    nba_immune_zoom_mapping <- c(
        nba_base_mapping[c("Neuroendocrine", "B cell", "Plasma")],
        # rename Myeloid to Mono/macro
        setNames(nba_base_mapping[c("Myeloid")], "Mono/macro"),
        # greens for NKs and TCs
        setNames(
            as.character(paletteer::paletteer_d("RColorBrewer::Greens", direction = -1)[seq_len(7)]),
            c(
                "TOX2+/KIT+ NK cell",
                "Resident NK cell",
                "Circulating NK cell",
                "NKT cell",
                "Treg",
                "CD4+ T cell",
                "CD8+ T cell"
            )
        ),
        # pinks for DC types
        setNames(
            # c("#f9c2b2", "#FB6A4AFF", "#EF3B2CFF"),
            as.character(paletteer::paletteer_d("PrettyCols::Pinks")[seq_len(4)]),
            c("pDC", "Migratory cDC", "cDC1", "cDC2/DC3")
        ),
        # manual color selection for others
        setNames(
            c("#E377C2FF", "#9467BDFF"),
            c("Cycling", "Neutrophil")
        )
    )
    mapping <- switch(
        reference,
        "NBAtlas" = {
            switch(
                annotation,
                "Cell_type" = nba_base_mapping,
                "annot_immunezoom_v5" = nba_immune_zoom_mapping,
                "immune_zoom_old" = {
                    map <- c(nba_base_mapping, nba_immune_zoom_mapping)
                    map <- map[unique(names(map))]
                    # rename Neuroendocrine to match old name 'NE'
                    # rename Myeloid color to Macrophage
                    names(map)[match(c("Neuroendocrine", "Mono/macro"), names(map))] <- c("NE", "Macrophage")
                    # rename others old names
                    names(map)[match("Cycling", names(map))] <- "Immune cycling"
                    map <- c(
                        map,
                        c(
                            "Patrolling monocyte" = "#cc650b",
                            "Classical monocyte" = "#994b08",
                            "Immune doublets" = "grey80",
                            "Immune low quality" = "grey80"
                        )
                    )
                    map
                },
                stop(sprintf("No colormapping defined for NBAtlas annotation '%s'", annotation))
            )
        },
        "ImmGen" = {
            switch(
                annotation,
                "label.main" = {
                    setNames(
                        as.character(
                            # skip the dark grey which could be mistaken for the default NA color
                            paletteer::paletteer_d("ggsci::category20_d3", n = 15 + 1)[-8]
                        ),
                        c(
                            "Stem cells",
                            "Macrophages",
                            "Fibroblasts",
                            "ILC",
                            "T cells",
                            "NKT",
                            "Endothelial cells",
                            "Epithelial cells",
                            "B cells",
                            "Plasma",
                            "DC",
                            "Neutrophils",
                            "Mast cells",
                            "Monocytes",
                            "Tgd"
                        )
                    )
                },
                stop(sprintf("No colormapping defined for ImmGen annotation '%s'", annotation))
            )
        },
        stop(sprintf("No colormapping defined for reference '%s'", reference))
    )

    if (fix_names) {
        names(mapping) <- fix_names(names(mapping))
    }

    mapping
}
