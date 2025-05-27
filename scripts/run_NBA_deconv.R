library(ggplot2)
library(future)
library(patchwork)
library(purrr)
library(rlang)
library(Seurat)
library(spacexr)

source("scripts/functions/utils.R")
source("scripts/functions/RCTD.R")

so_ls <- readRDS("data/ST_NB_seurat.rds") |>
    map(\(so) SplitObject(so, split.by = "orig.ident")) |>
    reduce(c)

rctd_name <- "NBAtlas"
rctd_anno_cols <- c("Cell_type", "Cell_type_wImmuneZoomAnnot")
rctd_doublet_mode <- "full"
rctd_max_multi_cell_types <- 3

out_dir <- mkdir(sprintf("results/data/deconvolution/RCTD/%s", rctd_name))
fig_dir <- mkdir(sub("data", "figs", out_dir))

# create RCTD SpatialRNA objects
message("Creating RCTD SpatialRNA objects...")
puck_ls <- imap(so_ls, \(so, sample) {
    SpatialRNA(
        coords = {
            GetTissueCoordinates(so, image = sample) |>
                `colnames<-`(c("x", "y"))
        },
        # RCTD requirements: "Counts should be untransformed count-level data."
        # Spatial assay contains SpotClean-transformed counts, which are not integers
        # therefore, round them
        counts = round(LayerData(so, assay = "Spatial", layer = "counts"), 0)
    )
})

# diagnostic nUMI plot
message("Saving diagnostic nUMI plot...")
pw <- imap(puck_ls, \(puck, sample) {
    plot_puck_continuous(
        puck = puck,
        barcodes = colnames(puck@counts),
        plot_val = puck@nUMI,
        ylimit = c(0, round(quantile(puck@nUMI, 0.9))),
        size = 0.8
    ) +
        labs(title = sample, x = NULL, y = NULL) +
        theme_void() +
        theme(
            plot.title = element_text(size = rel(0.85), hjust = 0.5),
            legend.box.margin = margin(0, 0, 0, 0),
            legend.margin = margin(0, 0, 0, 0),
            legend.key.width = unit(0.15, "cm"),
            legend.text = element_text(size = rel(0.55)),
            plot.margin = margin(0, 0, 0, 0)
        )
}) |>
    (\(.) {
        Reduce("+", .) +
            plot_layout(ncol = 3) +
            plot_annotation(title = "nUMI", caption = "max. nUMI = Q90")
    })()
ggsave(
    "nUMI_Spatial.png",
    path = fig_dir,
    width = 10,
    height = 6,
    dpi = 300,
    device = ragg::agg_png
)

# save some memory
rm(so_ls)
gc()

# run RCTD
plan(multisession, workers = min(8, future::availableCores()))

walk(rctd_anno_cols, \(anno_col) {
    ref <- qs::qread(sprintf("temp/%s/RCTD_%s_%s_human_reference_v2.qs", rctd_name, rctd_name, anno_col))
    cell_min_instance <- min(25, min(table(ref@cell_types)))

    message(sprintf(
        "[%s] Preparing RCTD replicates using %s cores and CELL_MIN_INSTANCE of %s...",
        anno_col,
        future::nbrOfWorkers(),
        cell_min_instance
    ))

    rctd <- create.RCTD.replicates(
        spatialRNA.replicates = puck_ls,
        reference = ref,
        replicate_names = names(puck_ls),
        group_ids = setNames(
            dplyr::case_when(
                grepl("NB1Pre", names(puck_ls)) ~ 1,
                grepl("NB1Post", names(puck_ls)) ~ 2,
                grepl("NB2Post", names(puck_ls)) ~ 3
            ),
            names(puck_ls)
        ),
        max_cores = future::nbrOfWorkers(),
        # lower below default (25) if a cell type from the reference has less cells
        CELL_MIN_INSTANCE = cell_min_instance,
        MAX_MULTI_TYPES = rctd_max_multi_cell_types
    )

    message(sprintf(
        "[%s] Running RCTD in '%s' mode on %s replicates...",
        anno_col,
        rctd_doublet_mode,
        length(puck_ls)
    ))
    exec_time <- system.time({
        rctd <- run.RCTD.replicates(rctd, doublet_mode = rctd_doublet_mode)
    })
    message(sprintf(
        "[%s] ... Execution finished after %s.",
        anno_col,
        proc_time_to_period(exec_time)
    ))

    print(names(rctd@RCTD.reps))
    names(rctd@RCTD.reps) <- sub("\\.", "_", names(puck_ls))

    rctd_path <- sprintf("%s/RCTD_%s_%s_%s.qs", out_dir, rctd_name, rctd_doublet_mode, anno_col)
    message(sprintf("Saving %s to %s", format(object.size(rctd), unit = "GB"), out_dir))
    qs::qsave(rctd, rctd_path)
})
