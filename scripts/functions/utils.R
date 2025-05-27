get_space_ranger_dir <- function(sample, base_path = "temp", check_exists = TRUE) {
    dir <- paste(base_path, sample, "outs", sep = "/")

    if (check_exists && !dir.exists(dir)) {
        stop(sprintf("Dir '%s' does not exist.", dir))
    }

    dir
}

get_space_ranger_files <- function(sample, base_path = "temp", check_exists = TRUE) {
    dir <- get_space_ranger_dir(sample, base_path, check_exists)

    file_ls <- list(
        "raw_feature_matrix" = paste(dir, "raw_feature_bc_matrix.h5", sep = "/"),
        "filtered_feature_matrix" = paste(dir, "filtered_feature_bc_matrix.h5", sep = "/"),
        "tissue_positions" = paste(dir, "spatial", "tissue_positions.csv", sep = "/"),
        "image" = paste(dir, "spatial", "tissue_hires_image.png", sep = "/"),
        "scale_factors" = paste(dir, "spatial", "scalefactors_json.json", sep = "/")
    )

    if (check_exists) {
        check_res <- sapply(file_ls, file.exists)
        if (!any(unlist(check_res))) {
            stop(sprintf(
                "%s file(s) do not exist.",
                paste(names(check_res[which(check_res == FALSE)]), collapse = ", ")
            ))
        }
    }

    setNames(list(file_ls), sample)
}

# see https://stackoverflow.com/a/27626007/8631547
chunks_of <- function(x, n) {
    sequence <- seq.int(from = 1, to = length(x), by = n)
    mapply(
        \(a, b) (x[a:b]),
        sequence,
        pmin(sequence + (n - 1), length(x)),
        SIMPLIFY = FALSE
    )
}

get_giotto_brain_cell_sig_mat <- function() {
    sign_matrix_path <- system.file("extdata", "sig_matrix.txt", package = "Giotto")
    brain_sc_markers <- data.table::fread(sign_matrix_path)
    mat <- as.matrix(brain_sc_markers[, -1])
    rownames(mat) <- brain_sc_markers$Event

    mat
}

get_mt_rp_genes <- function(so, return.data = TRUE) {
    stopifnot("so must be a SeuratObject" = "Seurat" %in% class(so))

    n_genes <- length(rownames(so))
    mt_genes <- grep("^mt-", rownames(so), value = TRUE)
    rp_genes <- grep("^Rpl|^Rps", rownames(so), value = TRUE)
    message(sprintf(
        "mitochondrial genes: %s (%.1f%%)\nribosomal genes: %s (%.1f%%)",
        length(mt_genes),
        length(mt_genes) / n_genes * 100,
        length(rp_genes),
        length(rp_genes) / n_genes * 100
    ))

    if (return.data) {
        list(
            "mt" = mt_genes,
            "rp" = rp_genes
        )
    }
}

# from https://www.r-bloggers.com/2016/10/converting-mouse-to-human-gene-names-with-biomart-package/
# as usual, contacting biomaRt is super slow or times out entirely :/
convert_genes <- function(genes, from, to, identifier = "symbol") {
    local({
        source("/home/peter/projects/alk_thinning/scripts/functions/organism_conversion.R")
        res <<- convert_genes_org(genes, identifier = "symbol", from, to)
    })

    res
}

prettify_geneset_label <- function(i, pref = NULL, max_len = 30) {
    # remove prefix/suffix , e.g. "HALLMARK_" or "_TARGET_GENES"
    if (!is.null(pref)) {
        i <- gsub(paste0("_?", pref, "_?"), "", i)
    }
    # replace underscore by space
    i <- gsub("_", " ", i)
    # enforce max length of characters
    sapply(i, \(j) {
        if (stringr::str_length(j) > max_len) {
            paste0(substring(j, 1, max_len - 3), "...")
        } else {
            j
        }
    })
}

mkdir <- function(dir) {
    if (!dir.exists(dir)) {
        dir.create(dir, recursive = TRUE)
    }
    dir
}

load_if_not_exist <- function(object_name, object_path) {
    if (exists(object_name)) {
        message(sprintf(
            "Using %s present in environment. class: %s; size: %s",
            object_name,
            paste(class(get(object_name)), collapse = ", "),
            format(object.size(get(object_name)), unit = "MB")
        ))
    } else {
        message(sprintf("Loading %s as %s", object_path, object_name))
        assign(object_name, readRDS(object_path), envir = parent.frame())
    }
}

ls_search <- function(haystack, needle, all = TRUE) {
    if (!needle %in% unlist(haystack)) {
        return(FALSE)
    }
    # case named list: return name
    res <- if (!is.null(attributes(haystack))) {
        lapply(haystack, \(values) needle %in% values) |>
            (\(.) Filter(isTRUE, .))() |>
            names()
    } else {
        # case unnamed list: return index
        which(sapply(haystack, \(values) needle %in% values))
    }

    if (!all && length(res) > 1) {
        warning(sprintf(
            "Multiple occurences of %s found in %s. Returning first.",
            needle,
            deparse(substitute(haystack))
        ))
        return(res[1])
    }

    return(res)
}

# fix the consecutive naming scheme, i.e. seurat_hallmark1, seurat_hallmark2, ...
rename_module_names <- function(so, gs_db, name) {
    stopifnot("so must be a SeuratObject" = "Seurat" %in% class(so))

    name_rgx <- sprintf("%s\\d", name)
    tmp <- so[[]][grep(name_rgx, colnames(so[[]]))]
    # fwo: remove "Seurat_" prefix
    # names(tmp) <- paste0("Seurat_", names(gs_db))
    names(tmp) <- names(gs_db)
    so <- AddMetaData(so, tmp)

    so@meta.data <- so@meta.data[, -grep(name_rgx, colnames(so@meta.data))]

    so
}

proc_time_to_period <- function(proc_time, attr = "elapsed") {
    if (!"proc_time" %in% class(proc_time)) {
        warning("Given time is not of class 'proc_time'")
        return(proc_time)
    }

    proc_time[attr] |>
        lubridate::seconds_to_period() |>
        round()
}

get_protein_coding_feat <- function(db, type = "genes", id = "symbol") {
    stopifnot("Type must be one of genes, transcripts" = type %in% c("genes", "transcripts"))
    stopifnot("id must be one of gene_id, gene_name, symbol" = id %in% c("gene_id", "gene_name", "symbol"))

    pc_filter <- ensembldb::filter(db, filter = ~ tx_biotype == "protein_coding")

    if (type == "genes") {
        return(
            as.data.frame(ensembldb::genes(pc_filter))[, id]
        )
    } else {
        return(
            as.data.frame(ensembldb::transcripts(pc_filter))[, id]
        )
    }
}

adjust_scale_factors <- function(so) {
    so@images <- so@images |>
        imap(\(img, img_name) {
            img@scale.factors$spot <- switch(
                img_name,
                "NB1Pre1" = 225,
                "NB1Pre2" = 285,
                "NB1Post1" = 240,
                "NB1Post2" = 265,
                "NB2Post1" = 58,
                "NB2Post2" = 48
            )
            # img@spot.radius <- 0.012
            img
        })
    so
}

# hgTables to GRanges
ucsc_hg_cytoband_to_gr <- function(hg) {
    stopifnot(c("chr", "start", "end", "name", "gieStain") %in% colnames(hg))
    gr <- GRanges(
        seqnames = hg$chr,
        ranges = IRanges(start = hg$start, end = hg$end),
        name = hg$name,
        gieStain = hg$gieStain
    )
    seqlengths(gr) <- as.vector(tapply(end(gr), seqnames(gr), max))

    gr
}

# from https://jokergoo.github.io/ComplexHeatmap-reference/book/genome-level-heatmap.html
average_in_window <- function(window, gr, v, method = "weighted", empty_v = NA) {
    if (missing(v)) v <- rep(1, length(gr))
    if (is.null(v)) v <- rep(1, length(gr))
    if (is.atomic(v) && is.vector(v)) v <- cbind(v)

    v <- as.matrix(v)
    if (is.character(v) && ncol(v) > 1) {
        stop("`v` can only be a character vector.")
    }

    if (length(empty_v) == 1) {
        empty_v <- rep(empty_v, ncol(v))
    }

    u <- matrix(rep(empty_v, each = length(window)), nrow = length(window), ncol = ncol(v))

    mtch <- as.matrix(findOverlaps(window, gr))
    intersect <- pintersect(window[mtch[, 1]], gr[mtch[, 2]])
    w <- width(intersect)
    v <- v[mtch[, 2], , drop = FALSE]
    n <- nrow(v)

    ind_list <- split(seq_len(n), mtch[, 1])
    window_index <- as.numeric(names(ind_list))
    window_w <- width(window)

    if (is.character(v)) {
        for (i in seq_along(ind_list)) {
            ind <- ind_list[[i]]
            if (is.function(method)) {
                u[window_index[i], ] <- method(v[ind], w[ind], window_w[i])
            } else {
                tb <- tapply(w[ind], v[ind], sum)
                u[window_index[i], ] <- names(tb[which.max(tb)])
            }
        }
    } else {
        if (method == "w0") {
            gr2 <- reduce(gr, min.gapwidth = 0)
            mtch2 <- as.matrix(findOverlaps(window, gr2))
            intersect2 <- pintersect(window[mtch2[, 1]], gr2[mtch2[, 2]])

            width_intersect <- tapply(width(intersect2), mtch2[, 1], sum)
            ind <- unique(mtch2[, 1])
            width_setdiff <- width(window[ind]) - width_intersect

            w2 <- width(window[ind])

            for (i in seq_along(ind_list)) {
                ind <- ind_list[[i]]
                x <- colSums(v[ind, , drop = FALSE] * w[ind]) / sum(w[ind])
                u[window_index[i], ] <- (x * width_intersect[i] + empty_v * width_setdiff[i]) / w2[i]
            }
        } else if (method == "absolute") {
            for (i in seq_along(ind_list)) {
                u[window_index[i], ] <- colMeans(v[ind_list[[i]], , drop = FALSE])
            }
        } else if (method == "weighted") {
            for (i in seq_along(ind_list)) {
                ind <- ind_list[[i]]
                u[window_index[i], ] <- colSums(v[ind, , drop = FALSE] * w[ind]) / sum(w[ind])
            }
        } else {
            if (is.function(method)) {
                for (i in seq_along(ind_list)) {
                    ind <- ind_list[[i]]
                    u[window_index[i], ] <- method(v[ind], w[ind], window_w[i])
                }
            } else {
                stop("wrong method.")
            }
        }
    }

    return(u)
}

#' Add meta data to SeuratObject using a join operation on arbitrary columns.
#'
#' This function replaces the meta.data slot of a SeuratObject with new meta data obtained through
#' a left join of old and new meta data. Of note, the new meta data must have the appropriate columns
#' to join by.
#'
#' @param so A SeuratObject to whose meta.data slot will be the left side of the join.
#' @param meta A data.frame containing the new meta data to be joined as the right side.
#' @param by A character vector specifying the columns to join by.
#'
#' @return A SeuratObject with the updated metadata.
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr separate
#' @importFrom dplyr left_join
#' @export
add_meta_data_by <- function(so, meta, by) {
    old_meta <- if ("barcode" %in% by && !"barcode" %in% colnames(so[[]])) {
        warning("Column barcode not found in SeuratObject. Assuming barcode from rownames.")
        so[[]] |>
            tibble::rownames_to_column("tmp") |>
            tidyr::separate("tmp", c("barcode", "sample_idx"), sep = "_")
    } else {
        so[[]]
    }

    stopifnot(all(by %in% colnames(old_meta)))

    new_meta <- dplyr::left_join(old_meta, meta, by = by)
    # preserve original rownames for general compatibility
    rownames(new_meta) <- rownames(so[[]])
    so@meta.data <- new_meta

    so
}

do_goa_fisher <- function(
    genes_deg,
    genes_bg,
    gsea_db,
    min_genes
) {
    stopifnot("GSEA DB must be a list" = class(gsea_db) == "list")
    n_genesets <- length(names(gsea_db))

    cols <- c(
        "n_genes_pw",
        "n_genes_pw_neg",
        "n_genes_pw_pos",
        "prop_neg",
        "prop_pos",
        "OR",
        "p",
        "q",
        "genes"
    )
    goa_table <- matrix(
        data = NA,
        nrow = n_genesets,
        ncol = length(cols),
        dimnames = list(names(gsea_db), cols)
    )

    for (i in seq_len(nrow(goa_table))) {
        genes_overlap_bg <- intersect(unlist(gsea_db[[i]]), genes_bg)
        goa_table[i, "n_genes_pw"] <- length(genes_overlap_bg)

        if (length(genes_overlap_bg) < min_genes) {
            warning(sprintf(
                "Geneset %s too small: n = %s but min_genes = %s. Skipping.",
                names(gsea_db[[i]]),
                length(gsea_db[[i]]),
                min_genes
            ))
            next
        } else {
            isRetrieved <- genes_bg %in% genes_deg
            inDB <- factor(genes_bg %in% genes_overlap_bg, levels = c(FALSE, TRUE))
            compare_pw_temp_t <- table(isRetrieved, inDB)

            # fix: fisher.test: 'x' must have at least 2 rows and columns
            if (any(dim(compare_pw_temp_t) < c(2, 2))) {
                warning(sprintf(
                    "Not enough hits in geneset %s. Skipping.",
                    names(gsea_db[[i]])
                ))
                next
            }

            goa_table[i, c("n_genes_pw_neg", "n_genes_pw_pos")] <- compare_pw_temp_t[, "TRUE"]
            goa_table[i, c("prop_neg", "prop_pos")] <- round(100 * prop.table(compare_pw_temp_t, 1)[, "TRUE"], 1)
            goa_table[i, "genes"] <- paste(genes_bg[isRetrieved & inDB == TRUE], collapse = ",")

            f_test <- fisher.test(compare_pw_temp_t, alternative = "greater")

            goa_table[i, "OR"] <- f_test[["estimate"]]
            goa_table[i, "p"] <- f_test[["p.value"]]
        }
    }

    goa_table <- goa_table[!is.na(goa_table[, "p"]) & !duplicated(rownames(goa_table)), ]
    goa_table <- as.data.frame(goa_table)
    goa_table[, 1:8] <- apply(goa_table[, 1:8], 2, as.numeric)

    goa_table$q <- p.adjust(goa_table$p, "fdr")

    return(goa_table[order(goa_table$p), ])
}
