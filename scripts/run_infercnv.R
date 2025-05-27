suppressPackageStartupMessages({
    library(furrr)
    library(future)
    library(infercnv)
    library(purrr)
    library(Seurat)
})

source("scripts/functions/utils.R")

memory_profiling <- FALSE

cell_type_column <- "cell_type_hr"
patient <- "NBLU02"
patient_short <- sub("LU", "", patient)
exclude_samples <- c("NB2Pre1", "NB2Pre2")
version <- "v4"
use_ref <- TRUE

# v2:
# - split chromosome by arms
# - inferCNV version ‘1.19.1’
# v3:
# - use Space Ranger's raw data matrix (not filtered matrix)
# - use gene order file based on reference genome (not reference probe set)
# - do not split chromosomes before running the analysis; this does not make any biological sense; split chromosomes
#   later for visualisation purposes
# - inferCNV version ‘1.20.0’
# - change `HMM_report_by` from the default `subclusters` to `cell` (i.e. ST spot); if we ever want to use that
#   info, it makes more sense to have this per spot
# v4:
# - like v3, but on chromsome arm level (cf. v2)
# v5:
# - like v4, but using GENCODE v32 (the same version used for the Visium probe set 2020-A)

message <- function(...) {
    base::message(sprintf(
        "[%s][%s][%s] %s",
        patient_short,
        cell_type_column,
        version,
        paste(...)
    ))
}

# TODO: if old versions should be rerun, switch to hg38_gencode_v27.txt in general
# using the probe set is probably not the best way to look at it because
# probe lengths != gene lengths
if (version %in% c("", "v2")) {
    message("DEPRECATED: Version v1 and v2 of this script are hard-deprecated. Check gene order file definition.")
    stop()
}
params <- list(
    # the gene order file input for inferCNV (either created by this script or pre-downloaded)
    gene_order_file_path = switch(
        version,
        "v5" = "temp/hg38_gencode_v32_chr_arms.txt",
        "v4" = "temp/hg38_gencode_v27_chr_arms.txt",
        "v3" = "temp/hg38_gencode_v27.txt",
        "v2" = "temp/inferCNV_Visium_Hs_ProbeSet_v1_gene_order_file_chr_arms.tsv",
        stop("Invalid version")
    ),
    # the genome BED used to create the gene order file
    genome_bed_path = switch(
        version,
        "v5" = "temp/hg38_gencode_v32.txt",
        "v4" = "temp/hg38_gencode_v27.txt",
        "v3" = NULL,
        "v2" = "temp/Visium_Human_Transcriptome_Probe_Set_v1.0_GRCh38-2020-A.bed",
        stop("Invalid version")
    ),
    space_ranger_matrix = switch(
        version,
        "v5" = "raw_feature_bc_matrix.h5",
        "v4" = "raw_feature_bc_matrix.h5",
        "v3" = "raw_feature_bc_matrix.h5",
        "v2" = "filtered_feature_bc_matrix.h5",
        stop("Invalid version")
    ),
    resolution = switch(
        version,
        "v5" = "chr_arm",
        "v4" = "chr_arm",
        "v3" = "chr",
        "v2" = "chr_arm",
        stop("Invalid version")
    ),
    chr_exclude = switch(
        version,
        "v5" = c("Xp", "Xq", "Yp", "Yq", "MNA"),
        "v4" = c("Xp", "Xq", "Yp", "Yq", "MNA"),
        "v2" = c("Xp", "Xq", "Yp", "Yq"),
        c("chrX", "chrY", "chrM")
    )
)

out_dir <- mkdir(sprintf("temp/inferCNV/%s/%s/%s", (if (use_ref) "" else "noRef"), cell_type_column, version))
data_dir <- "/home/joachim/projects/NB_spatialT/temp/data/"

sr_dirs <- sprintf("spaceranger_out_%s", patient) |>
    lapply(\(i) list.dirs(paste0(data_dir, i), recursive = FALSE)) |>
    unlist()

samples <- sr_dirs |>
    strsplit("/") |>
    lapply(\(i) i[length(i)]) |>
    unlist() |>
    (\(.) gsub("^NB_?LU_?0([1-9]+)_?(Pre|Post)_?(\\d+)", "NB\\1\\2\\3", .))()

names(sr_dirs) <- samples
sr_dirs <- sr_dirs[!names(sr_dirs) %in% exclude_samples]
samples <- samples[!samples %in% exclude_samples]

message("Loading merged data...")
# this contains data transformed with SpotClean;
# only used to obtain cell type annotation;
# raw data as input for inferCNV obtained later from Space Ranger outs
# NB01: 12703 features, 8415 spots, 4 samples (NB1Pre1, NB1Pre2, NB1Post1, NB1Post2)
# NB02: 14140 features, 4582 spots, 4 samples (NB2Pre1, NB2Pre2, NB2Post1, NB2Post2)
so_merged <- readRDS(sprintf("%s/seurat_mergedAll_list.rds", data_dir)) |>
    pluck(patient_short) |>
    subset(orig.ident %in% samples)

#### prepare data ####
# sanity check and subset
if (!all(samples %in% so_merged[[]]$orig.ident)) {
    exclude_idx <- which(!samples %in% so_merged[[]]$orig.ident)
    exclude_samples <- samples[exclude_idx]
    message(sprintf(
        "Excluding %s samples because they weren't found in so_merged (%s)",
        length(exclude_samples),
        paste(exclude_samples, collapse = ", ")
    ))
    samples <- samples[!samples %in% exclude_samples]
    sr_dirs <- sr_dirs[!seq_along(sr_dirs) %in% exclude_idx]
}

# gene order file (provided by inferCNV/Trinity CTAT)
# from: https://data.broadinstitute.org/Trinity/CTAT/cnv/hg38_gencode_v27.txt
message("Obtaining gene order file...")
if (!file.exists(params$gene_order_file_path)) {
    # adapted for v3, v4, deprecated versions not supported anymore
    bed <- read.csv(
        params$genome_bed_path,
        sep = "\t",
        header = FALSE
    )
    # bed_sym <- bed$V4 |>
    #     strsplit("\\|") |>
    #     lapply("[[", 2) |>
    #     unlist()
    genome_anno <- as.data.frame(bed)

    if (any(duplicated(genome_anno[, 1]))) {
        genome_anno <- genome_anno[!duplicated(genome_anno[, 1]), ]
    }

    if (params$resolution == "chr_arm") {
        hg38_cyto <- data.table::fread(
            "https://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBand.txt.gz",
            col.names = c("chr", "start", "end", "name", "gieStain")
        )
        # chromosome arms
        chr_arms <- hg38_cyto |>
            within(arm <- substring(name, 1, 1)) |>
            subset(arm != "") |>
            (\(.) {
                by(., list(.$chr, .$arm), \(x) {
                    cbind(x$chr, x$arm, min(x$start), max(x$end), sum(x$end - x$start))
                })
            })() |>
            lapply(unique) |>
            do.call(what = rbind.data.frame) |>
            `colnames<-`(c("chr", "arm", "start", "end", "length")) |>
            dplyr::mutate(
                chr = factor(
                    chr,
                    levels = stringr::str_sort(unique(chr), numeric = TRUE)
                ),
                start = as.numeric(start),
                end = as.numeric(end),
                length = as.numeric(length)
            ) |>
            dplyr::arrange(chr)

        genome_anno <- genome_anno |>
            `colnames<-`(c("gene", "chr", "start", "end")) |>
            dplyr::group_by(chr) |>
            dplyr::left_join(
                chr_arms,
                by = dplyr::join_by(chr, within(x$start, x$end, y$start, y$end)),
                suffix = c("_input", "_hg38_ref")
            ) |>
            dplyr::mutate(chr = paste0(sub("chr", "", chr), arm)) |>
            dplyr::select(gene, chr, start_input, end_input) |>
            dplyr::filter()
    }

    write.table(
        genome_anno,
        file = params$gene_order_file_path,
        sep = "\t",
        quote = FALSE,
        col.names = FALSE,
        row.names = FALSE
    )
}

# not more than 8 threads to run inferCNV
n_cores_lvl1 <- length(samples) + 1
n_cores_lvl2 <- max(1, (availableCores() - n_cores_lvl1) %/% n_cores_lvl1)
plan(list(
    tweak(multisession, workers = I(n_cores_lvl1)),
    tweak(multisession, workers = I(min(8, n_cores_lvl2)))
))

message("Obtaining raw data from Space Ranger output...")
# get raw data from Space Ranger output via Seurat
walk2(
    sr_dirs,
    samples,
    \(path, sample) {
        message("Loading data for", sample)
        mkdir(sprintf("%s/%s", out_dir, sample))

        so <- Load10X_Spatial(
            sprintf("%s/outs", path),
            filename = params$space_ranger_matrix,
            # always remove background spots
            filter.matrix = TRUE,
            slice = sample
        )

        # raw count matrix
        LayerData(so, assay = "Spatial", layer = "counts") |>
            as.data.frame() |>
            write.table(
                sprintf("%s/%s/raw_counts.tsv", out_dir, sample),
                sep = "\t",
                quote = FALSE
            )

        # annotations file
        so_anno <- subset(so_merged, orig.ident %in% sample)
        anno <- cbind.data.frame(rownames(so_anno[[]]), as.character(so_anno[[]][[cell_type_column]]))
        colnames(anno) <- c("id", cell_type_column)
        anno$simple <- ifelse(grepl("^NE", anno[, cell_type_column]), "tumour", "normal")
        anno$id <- gsub("^.+?_", "", anno$id)

        write.table(
            anno[, c("id", cell_type_column)],
            sprintf("%s/%s/annotations.tsv", out_dir, sample),
            sep = "\t",
            quote = FALSE,
            col.names = FALSE,
            row.names = FALSE
        )
    }
    # .options = furrr_options(seed = NULL)
)

#### run inferCNV ####
if (use_ref) {
    normal_cells <- so_merged[[]][[cell_type_column]] |>
        (\(.) .[!grepl("^NE", .)])() |>
        unique() |>
        as.character()

    message(sprintf("Using %s as normal cells.", paste(normal_cells, collapse = ", ")))
} else {
    message("Running inferCNV reference-free.")
}


future_walk(
    samples,
    \(sample) {
        if (memory_profiling) {
            Rprof(
                mem_prof <- sprintf(
                    "%s/%s/_rprof_%s_%s.log",
                    out_dir,
                    sample,
                    format(Sys.time(), "%Y%m%d_%H%M%S"),
                    strsplit(uuid::UUIDgenerate(), "-")[[1]][1]
                ),
                memory.profiling = TRUE
            )
        }
        # required for analysis_mode = "subclusters"
        # see https://github.com/broadinstitute/infercnv/issues/396
        options(scipen = 100)

        # save inferCNV version
        writeLines(
            infercnv_version <- as.character(packageVersion("infercnv")),
            file(sprintf("%s/%s/_infercnv_version.txt", out_dir, sample))
        )

        n_workers <- future::nbrOfWorkers()
        message(sprintf(
            "[%s] Running inferCNV (v%s) with %s threads and memory profiling %s.",
            sample,
            infercnv_version,
            n_workers,
            (if (memory_profiling) "enabled" else "disabled")
        ))

        # only pass normal cells present in the current samples or
        # object creation will fail because of wonky if statement
        if (use_ref) {
            normals_present <- unique(read.table(sprintf("%s/%s/annotations.tsv", out_dir, sample))[, 2])
            normals_present <- normal_cells[normal_cells %in% normals_present]
        }

        infercnv::CreateInfercnvObject(
            raw_counts_matrix = sprintf("%s/%s/raw_counts.tsv", out_dir, sample),
            gene_order_file = params$gene_order_file_path,
            annotations_file = sprintf("%s/%s/annotations.tsv", out_dir, sample),
            ref_group_names = (if (use_ref) normals_present else NULL),
            chr_exclude = params$chr_exclude
        ) |>
            infercnv::run(
                # use 1 for smart-seq, 0.1 for 10x-genomics
                cutoff = 0.1,
                # dir is auto-created for storing outputs
                out_dir = sprintf("%s/%s/out", out_dir, sample),
                cluster_by_groups = TRUE,
                tumor_subcluster_pval = 0.05,
                denoise = TRUE,
                HMM = TRUE,
                HMM_type = "i3",
                HMM_report_by = "cell",
                # leave this at subcluster-level (default)
                analysis_mode = "subclusters",
                num_threads = n_workers,
                write_expr_matrix = TRUE
            )

        if (memory_profiling) Rprof(NULL)
    },
    .options = furrr_options(seed = TRUE)
)

message("Done.")
