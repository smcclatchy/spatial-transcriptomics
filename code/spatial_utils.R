# Custom util functions

suppressPackageStartupMessages(require(Matrix))
suppressPackageStartupMessages(require(Seurat))

#' Add metadata to Seurat object.
#'
#' This function is a simple wrapper around AddMetaData
#'
#' @param obj A Seurat object
#' @param metadata.df A data.frame expected to have all of the rownames (corresponding to cells/spots) in existing metadata for obj (obj[[]])
#' @return obj modified to hold the metadata in metadata.df
add.metadata.to.seurat.obj <- function(obj, metadata.df) {
  # rows <- rownames(obj[[]])
  # print(head(rows))
  # print(table(rows %in% rownames(metadata.df)))
  # stopifnot(rows %in% rownames(metadata.df))
  # metadata.df <- metadata.df[rows,]
  colnames(metadata.df) <- make.names(colnames(metadata.df))
  metadata.df <- metadata.df[, !(colnames(metadata.df) %in% colnames(obj[[]])), drop=FALSE]
  cols.to.add <- colnames(metadata.df)
  all.meta <- merge(obj[[]], metadata.df, by = "row.names", all.x = TRUE)
  rownames(all.meta) <- all.meta$Row.names
  ids <- Cells(obj)
  all.meta <- all.meta[ids,]
  obj <- AddMetaData(obj, all.meta[, cols.to.add, drop=FALSE])
  obj
}


#' Get the position information for each spot in a tissue
#'
#' A function that extracts spot position information from the spaceranger output corresponding to
#' a single tissue. See https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/output/images
#' for a full description.
#'
#' @param spaceranger_dir The path to the spaceranger output.
#' @return A data.frame with columns:
#'  barcode: The sequence of the barcode associated to the spot.
#'  in_tissue: Binary, indicating if the spot falls inside (1) or outside (0) of tissue.
#'  array_row: The row coordinate of the spot in the array from 0 to 77. The array has 78 rows.
#'  array_col: The column coordinate of the spot in the array. In order to express the orange crate arrangement of the spots, this column index uses even numbers from 0 to 126 for even rows, and odd numbers from 1 to 127 for odd rows. Notice then that each row (even or odd) has 64 spots.
#'  pxl_row_in_fullres: The row pixel coordinate of the center of the spot in the full resolution image.
#'  pxl_col_in_fullres: The column pixel coordinate of the center of the spot in the full resolution image.
get.tissue.position.metadata <- function(spaceranger_dir) {
  image.dir <- paste0(spaceranger_dir, "/", "spatial/")

  file <- file.path(image.dir, "tissue_positions_list.csv")
  col.names <- c("barcodes", "tissue", "row", "col", "imagerow", "imagecol")
  if(file.exists(file)) {
    tissue.positions <-
      read.csv(file = file, col.names = col.names, header = FALSE, as.is = TRUE, row.names = 1)
    return(tissue.positions)
  }
  file <- file.path(image.dir, "tissue_positions.csv")
  if(!file.exists(file)) {
    stop("Found neither tissue_positions_list.csv nor tissue_positions.csv")
  }
  tissue.positions <-
    read.csv(file = file, header = TRUE, as.is = TRUE)
  colnames(tissue.positions) <- col.names
  rownames(tissue.positions) <- tissue.positions$barcodes
  tissue.positions <- tissue.positions[, !(colnames(tissue.positions) %in% c("barcodes"))]
  return(tissue.positions)
}

#' Save a ggplot to file with specified dimensions and resolution (optional)
#'
#' This function saves a ggplot object to a file. The user specifies the filename, and can optionally
#' specify the resolution, width, and height for the saved plot.
#'
#' @param plot The plot object to be saved, typically a ggplot.
#' @param filename The name of the file where the plot should be saved. This can include a path.
#' @param dpi The resolution in dots per inch for the saved plot. Default is 300 dpi.
#' @param width The width of the saved plot in inches. Default is 10 inches.
#' @param height The height of the saved plot in inches. Default is 8 inches.
#' @return None
plot_and_save <- function(plot, filename, dpi = 300, width = 10, height = 8) {
  # Use ggsave to save the plot with specified or default parameters
  ggsave(filename = file.path(filename), plot = plot, dpi = dpi, width = width, height = height)
}

#' Run RCTD deconvolution by first ensuring that single-cell RNA-seq reference only includes genes in ST data.
#'
#' @param st.obj A Seurat object for the ST data
#' @param sc.counts The reference scRNA-seq count matrix (with genes as rows and cells as columns).
#'                  Genes should be in the same namespace (e.g., symbols or Ensembl ids) as the ST count data in st.obj.
#' @param sc.cell.types The cell type assignments of each cell in sc.counts stored in a vector, named according to the columns/cells of sc.counts
#' @param rds.output.file If non-null, an output file name to store the RCTD results. Default: NULL.
#'                        If the file exists, the output will be read and returned by this function.
#' @param run.rctd Boolean indicating whether to actually deconvolve data (run.RCTD) or simply to create.RCTD.
#' @param intersect.sc.and.st.genes Boolean indicating whether to subset scRNA-seq genes to those within
#'        ST data.
#' @return An RCTD object
#'
#' Assuming st.obj is the filtered feature obj, the associated count matrix should _already_ be filtered to include only
#' targeted genes, as described here:
#' https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/output/matrices
rctd.wrapper <- function(st.obj, sc.counts, sc.cell.types, rds.output.file = NULL, run.rctd = TRUE, intersect.sc.and.st.genes = FALSE) {
  if(!is.null(rds.output.file) && file.exists(rds.output.file)) {
    myRCTD <- readRDS(rds.output.file)
    gc()
    return(myRCTD)
  }

  # Get raw ST coiunts
  st.counts <- GetAssayData(st.obj, assay="Spatial", slot="counts")

  # Get the spot coordinates
  st.coords <- st.obj[[]][, c("col", "row")]
  colnames(st.coords) <- c("x","y")

  # Get the # of UMI per spot
  nUMI <- colSums(st.counts) # In this case, total counts per pixel is nUMI

  # Create the RCTD 'puck', representing the ST data
  puck <- SpatialRNA(st.coords, st.counts, nUMI)

  # Create an RCTD scRNA-seq reference.
  # NB: we do this for each of potentially multiple ST samples because some of those
  # samples may be from FFPE assays and others from fresh frozen. These will have
  # different gene sets that should be reflected (i.e., filtered) in the scRNA-seq reference (sc.counts).

  if(intersect.sc.and.st.genes) {
    stop("Would not recommend restricting scRNA-seq genes to those in ST!")
    # spacexr/RCTD:get_de_genes has an expression threshold (expr_threshold)
    # that any marker gene must exceed. Restricting the set of genes (to those in ST),
    # impacts their normalization (by total UMI) in RCTD. I have observed that this
    # can lead to significant differences in the number of genes passing the
    # expression threshold (e.g., twice as many following subsetting the genes).
    # Note that get_de_genes will only consider as markers those genes in common
    # anyways between scRNA-seq and ST. As such, it intuitively makes sense
    # to consistently define the normalized expression values of the scRNA-seq data
    # irrespective of the ST data -- i.e., without subsetting.
    sc.counts <- sc.counts[rownames(sc.counts) %in% rownames(st.counts),]
  }
  sc.cell.types <- sc.cell.types[colnames(sc.counts)]
  nUMI <- colSums(sc.counts)
  reference <- Reference(sc.counts, sc.cell.types, nUMI, n_max_cells = ncol(sc.counts) + 1)

  max_cores <- min(5, detectCores() - 1)

  # Note that we have keep_reference = TRUE here. Without it, the default is
  # keep_reference = FALSE, which would ignore the Reference we have
  # intentionally created above.
  myRCTD <- create.RCTD(puck, reference, max_cores = max_cores, keep_reference = TRUE)

  # NB: running in full mode!
  if(run.rctd) {
    myRCTD <- suppressPackageStartupMessages(run.RCTD(myRCTD, doublet_mode = 'full'))
  }

  if(!is.null(rds.output.file)) {
    saveRDS(myRCTD, rds.output.file)
  }

  gc()
  myRCTD
}

plot_tissue_prc_merge <-function(brain.merge){
  metadata <- brain.merge@meta.data
  # Create a data frame with 'sample', 'cluster', and 'count'
  data <- metadata %>%
    group_by(sample = orig.ident, cluster = Idents(brain.merge)) %>%
    summarise(count = n(), .groups = 'drop')

  total_counts_per_cluster <- data %>%
    group_by(cluster) %>%
    summarise(total_cluster_count = sum(count), .groups = 'drop')

  data <- data %>%
    left_join(total_counts_per_cluster, by = "cluster") %>%
    mutate(proportion = count / total_cluster_count) %>%
    select(-total_cluster_count)  # remove the total_cluster_count column as it's no longer needed

  total_counts <- data %>%
    group_by(sample) %>%
    summarise(count = sum(count), .groups = 'drop')

  # Calculate the grand total of all counts
  grand_total <- sum(total_counts$count)

  # Add a 'Total' proportion for each sample (divided by grand total for all clusters)
  total_counts <- total_counts %>%
    mutate(cluster = "Total", proportion = count / grand_total)

  data <- bind_rows(total_counts, data)

  data$cluster <- factor(data$cluster, levels = c("Total", unique(data$cluster)[-1]))

  ari <- adjustedRandIndex(brain.merge@meta.data[["layer_guess"]], Idents(brain.merge))

  plot <- ggplot(data, aes(x = cluster, y = proportion, fill = sample)) +
    geom_bar(stat = "identity") +
    scale_fill_brewer(palette = "Set1") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          legend.position = "bottom") +
    labs(x = "Cluster", y = "Proportion", fill = "Sample") +
    ggtitle(paste0("Proportion of Clusters by Sample with ARI:",ari))

  plot
}

#' Use a color safe palette in SpatialDimPlot
SpatialDimPlotColorSafe <- function(obj, group.by, ...) {
  vals <- unique(obj[[]][,group.by])
  cols <- setNames(carto_pal(length(vals), "Safe"), vals)

  g <- SpatialDimPlot(obj, group.by, cols = cols, label = TRUE, repel = TRUE, combine = FALSE, ...)[[1]]
  g <- g + guides(fill = guide_legend(override.aes = list(size=10)))
  g <- g + theme(text = element_text(size = 20))
  g
}

#' Extract the deconvolved fractions from RCTD results.
#'
#' @param rctd An RCTD object, from the spacexr package.
#' @param normalize Boolean indicating whether weights should be normalized to sum to 1.
#' @return A data frame whose rows are spots, whose columns are deconvolved populations, and whose
#'                    entries are the (predicted) fraction of a population in a given spot.
#'                    Also adds columns x and y, holding spatial coordinates of spots.
format.rctd.output_ <- function(rctd, normalize = FALSE) {
  barcodes <- colnames(rctd@spatialRNA@counts)
  weights <- rctd@results$weights
  if(normalize) {
    weights <- normalize_weights(weights)
  }
  df <- as.data.frame(weights)
  df$x <- rctd@spatialRNA@coords$x
  df$y <- rctd@spatialRNA@coords$y
  df
}

plot_RCTD_results <- function(obj, myRCTD.all) {
  df <- format.rctd.output_(myRCTD.all, normalize=FALSE)
  colnames(df) <- make.names(colnames(df))
  colnames(df)[colnames(df) == "Oligodendrocytes"] <- "Oligo"
  #regions <- c("AST.FB", "L2_3", "L4", "L5_6", "Oligo")
  regions <- c("AST.FB", "L2.3", "L4", "L5.6", "Oligo")
  #regions <- c("AST.FB", "L2_3", "L4", "L5_6", "Oligo")


  rctd.prediction <- apply(df[, regions], 1, function(row)
    names(row)[which.max(row)])
  df$Prediction <- rctd.prediction
  obj <- add.metadata.to.seurat.obj(obj, df) # 'objs' is replaced with 'seurat_object'

  anno.col <- "Prediction"
  flag <- is.na(obj[[]][,anno.col])
  obj <- subset(obj, cells = Cells(obj)[!flag])
  gt.col <- "layer_guess"
  flag <- is.na(obj[[]][,gt.col])
  obj <- subset(obj, cells = Cells(obj)[!flag])

  gt.vals <- unique(obj[[]][, gt.col])
  gt.cols <- setNames(carto_pal(length(gt.vals), "Safe"), gt.vals)

  vals <- unique(obj[[]][,anno.col])
  cols <- setNames(carto_pal(length(vals), "Safe"), vals)

  g1 <- SpatialDimPlot(obj, c(anno.col), cols = cols, label = TRUE, repel = TRUE, combine = FALSE)[[1]]
  g1 <- g1 + guides(fill = guide_legend(override.aes = list(size=10)))
  g1 <- g1 + theme(text = element_text(size = 20))

  g2 <- SpatialDimPlot(obj, c(gt.col), cols = gt.cols, label = TRUE, repel = TRUE, combine = FALSE)[[1]]
  g2 <- g2 + guides(fill = guide_legend(title = "Annotation", override.aes = list(size=10)))
  g2 <- g2 + theme(text = element_text(size = 20))

  df <- as.data.frame(table(obj[[]][,gt.col], obj[[]]$Prediction))
  colnames(df) <- c("Annotation", "Prediction", "Freq")
  df$Annotation <- factor(df$Annotation)
  df$Prediction <- factor(df$Prediction)

  g <- ggplot(data = df, aes(x = Annotation, y = Prediction, fill = Freq)) + geom_tile()
  g <- g + theme(text = element_text(size = 20))

  return(list("prediction_plot" = g1, "ground_truth_plot" = g2, "confusion_matrix" = g))
}


print_RCTD_results <- function(obj, myRCTD.all, prefix) {
  df <- format.rctd.output_(myRCTD.all, normalize=FALSE)
  colnames(df) <- make.names(colnames(df))
  colnames(df)[colnames(df) == "Oligodendrocytes"] <- "Oligo"
  #regions <- c("AST.FB", "L2_3", "L4", "L5_6", "Oligo")
  regions <- c("AST.FB", "L2.3", "L4", "L5.6", "Oligo")
  #regions <- c("AST.FB", "L2_3", "L4", "L5_6", "Oligo")


  rctd.prediction <- apply(df[, regions], 1, function(row)
    names(row)[which.max(row)])
  df$Prediction <- rctd.prediction
  obj <- add.metadata.to.seurat.obj(obj, df) # 'objs' is replaced with 'seurat_object'

  gt.col <- "Maynard"
  gt.vals <- unique(obj[[gt.col]][,gt.col])
  all.unique.predictions <-gt.vals
  all.gt.cols <- carto_pal(length(all.unique.predictions), "Safe")
  names(all.gt.cols) <- all.unique.predictions

  anno.col <- "Prediction"
  flag <- is.na(obj[[]][,anno.col])
  obj <- subset(obj, cells = Cells(obj)[!flag])
  vals <- unique(obj[[anno.col]][,anno.col])
  cols <- c()
  if ("AST.FB" %in% vals){
    f_val="AST.FB"
  } else if ("AST-FB" %in% vals){
    f_val="AST-FB"
  }
  if(f_val  %in% vals) {
    if("1" %in% gt.vals) {
      cols[f_val] <- all.gt.cols["1"]
    } else if("1_6" %in% gt.vals) {
      cols[f_val] <- all.gt.cols["1_6"]
    } else {
      cols[f_val] <- all.gt.cols["1"]
    }
  }

  if ("L2_3" %in% vals){
    sec_val="L2_3"
  } else if ("L2.3" %in% vals){
    sec_val="L2.3"
  } else if ("L2-3" %in% vals){
    sec_val="L2-3"
  }
  if( sec_val%in% vals) {
    cols[sec_val] <- all.gt.cols["2_3"]
  }
  if("L4" %in% vals) {
    cols["L4"] <- all.gt.cols["4"]
  }

  if ("L5_6" %in% vals){
    third_val="L5_6"
  } else if ("L5.6" %in% vals){
    third_val="L5.6"
  } else if ("L5-6" %in% vals){
    third_val="L5-6"
  }
  if(third_val %in% vals) {
    if("5_6" %in% gt.vals) {
      cols[third_val] <- all.gt.cols["5_6"]
    } else if("6" %in% gt.vals) {
      cols[third_val] <- all.gt.cols["6"]
    } else if("5" %in% gt.vals) {
      cols[third_val] <- all.gt.cols["5"]
    } else {
      cols[third_val] <- all.gt.cols["5_6"]
    }
  }
  if("Oligo" %in% vals) {
    cols["Oligo"] <- all.gt.cols["WM"]
  }

  g1 <- SpatialDimPlot(obj, c(anno.col), cols = cols, label = TRUE, repel = TRUE, combine = FALSE)[[1]]
  g1 <- g1 + guides(fill = guide_legend(override.aes = list(size=10)))
  g1 <- g1 + theme(text = element_text(size = 20))
  print(g1)

  g2 <- SpatialDimPlot(obj, c(gt.col), cols = all.gt.cols, label = TRUE, repel = TRUE, combine = FALSE)[[1]]
  g2 <- g2 + guides(fill = guide_legend(title = "Annotation", override.aes = list(size=10)))
  g2 <- g2 + theme(text = element_text(size = 20))

  png(paste0(img_dir, "/",prefix,"_prediction-vs-annotation-", levels(Idents(obj)[[1]]), ".png"), width = 2 * 480)
  print(g1 + g2)
  d <- dev.off()

  df <- as.data.frame(table(obj[[]]$Maynard, obj[[]]$Prediction))
  colnames(df) <- c("Annotation", "Prediction", "Freq")
  df$Annotation <- factor(df$Annotation)
  df$Prediction <- factor(df$Prediction)

  g <- ggplot(data = df, aes(x = Annotation, y = Prediction, fill = Freq)) + geom_tile()
  g <- g + theme(text = element_text(size = 20))
  png(paste0(img_dir, "/",prefix,"_prediction-vs-annotation-heatmap-", levels(Idents(obj)[[1]]), ".png"))
  print(g)
  d <- dev.off()
}


# Function to apply the threshold and plot quality control vars
apply_qc_threshold <- function(seurat_object, metric_name, threshold, is_greater=FALSE) {
  if (is_greater) {
    qc_flag <- seurat_object@meta.data[[metric_name]] > threshold
  } else {
    qc_flag <- seurat_object@meta.data[[metric_name]] < threshold
  }

  seurat_object@meta.data[paste0("qc_", metric_name)] <- qc_flag
  plot1 <- SpatialDimPlot(seurat_object, group.by = paste0("qc_", metric_name)) + theme(legend.position = "right")
  print(plot1)
  return(qc_flag)
}

################################################################################
# Given a Seurat object, write out each counts matrix to separate files,
# then null out the counts slots and save the Seurat object.
# We will write out the counts for each Assay, then set them to an empty
# matrix before saving the object.
#
# We write out the files as file_prefix, Assay, Layer, so that we can
# assemble the Seurat object later.
#
# Arguments:
# obj: Seurat object to save.
# file_prefix: character string to append to files that we have.
# Returns: Nothing.
save_seurat_object <- function(obj, file_prefix) {

  # Save the object metadata.
  saveRDS(obj[[]], file = paste0(file_prefix, "_meta_data.rds"))

  # Get the Assays in this object.
  obj_assays <- Assays(obj)

  # For each Assay, ...
  for(a in obj_assays) {

    # Set the default assay and get the Layers.
    DefaultAssay(obj) <- a
    assay_layers      <- Layers(obj)

    # For each Layer, write it out and null out the layer.
    for(c in assay_layers) {

      # Save the current counts matrix.
      saveRDS(object = LayerData(obj, c),
              file   = paste0(file_prefix, "_", a, "_", c, ".rds"))

      # Null out the current counts matrix.
      if(class(LayerData(obj, c))[1] == "dgCMatrix") {

        # Sparse matrix.
        LayerData(obj, c) <- SparseEmptyMatrix(nrow     = nrow(obj),
                                               ncol     = ncol(obj),
                                               rownames = rownames(obj),
                                               colnames = colnames(obj))

      } else {

        # Normal matrix.
        LayerData(obj, c) <- matrix()

      } # else

    } # for(c)

  } # for(a)

  # Save the stripped-down Seurat object.
  saveRDS(obj, file = paste0(file_prefix, "_seurat_obj.rds"))

} # save_seurat_object()


################################################################################
# Given a file prefix, look for the files that start with that prefix,
# infer the Assays and Layers and try to reassemble the Seurat object.
#
# Arguments:
# file_prefix: character string to append to files that we have.
# Returns:
# Seurat object with all of the Assay Layers populated.
load_seurat_object <- function(file_prefix) {

  # Get the files that start with the prefix.
  files <- dir(path = '.', pattern = file_prefix)

  # Find the file containing the Seurat object and read it in.
  wh  <- grep("_seurat_obj.rds", files)
  obj <- readRDS(files[wh])

  # Remove the Seurat object filename for the next part.
  files <- files[-wh]

  # Read in the cspot metadata and assign it to the Seurat object metadata.
  wh      <- grep("_meta_data.rds", files)
  obj[[]] <- readRDS(files[wh])

  # Remove the metadata filename for the next part.
  files <- files[-wh]

  # Create a data.frame with the filenames, assay, and layer.
  file_parts <- strsplit(sub("\\.rds$", "", files), split = '_')
  file_info  <- data.frame(files = files,
                           assay = sapply(file_parts, '[', 2),
                           layer = sapply(file_parts, '[', 3))

  # Go through each row in the file data.frame, read in the file, and place
  # it's contents in the correct Assay and Layer.
  for(i in 1:nrow(file_info)) {

    current_assay <- file_info$assay[i]
    current_layer <- file_info$layer[i]

    #DefaultAssay(obj) <- current_assay

    LayerData(obj, layer = current_layer, assay = current_assay) <- readRDS(file_info$files[i])

  } # for(i)

  return(obj)

} # load_seurat_object()


################################################################################
# Given an RCTD object, zero out the Reference counts (by replacing them with
# a small, dummy matrix.
# Arguments:
# rctd.results: RCTD results by run.RCTD
# Returns: rctd.results with Reference counts replaced by small dummy matrix
remove.RCTD.reference.counts <- function(rctd.results) {
  dummy_cnts = matrix(data=10000,nrow=2,ncol=2)
  rownames(dummy_cnts) <- rep(1:nrow(dummy_cnts))
  colnames(dummy_cnts) <- rep(1:ncol(dummy_cnts))
  cell_types <- as.factor(c("a","b"))
  names(cell_types) <- colnames(dummy_cnts)
  rctd.results@reference <- Reference(counts=dummy_cnts,cell_types=cell_types)
  rctd.results
}

