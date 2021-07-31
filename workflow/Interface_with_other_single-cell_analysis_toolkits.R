# Interface with other single-cell analysis toolkits ----

## Create a CellChat object from Seurat or SingleCellExperiment object ----

### Converting between single-cell objects (Seurat, SingleCellExperiment and anndata objects) ----

## Create a CellChat object from a data matrix extracting from Seurat and Scanpy tools ----

### Data input required in CellChat ----

#### Data format ----

### Extract the CellChat input files from a Seurat V3 object ----

data.input <- GetAssayData(seurat_object, assay = "RNA", slot = "data") # normalized data matrix
labels <- Idents(seurat_object)
meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels

### Extract the CellChat input files from a Seurat V2 object ----

data.input <- seurat_object@data # normalized data matrix
labels <- seurat_object@idents
meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels

### Extract the CellChat input files from a Scanpy object ----

library(reticulate)
ad <- import("anndata", convert = FALSE)
ad_object <- ad$read_h5ad("scanpy_object.h5ad")
# access normalized data matrix
data.input <- t(py_to_r(ad_object$X))
rownames(data.input) <- rownames(py_to_r(ad_object$var))
colnames(data.input) <- rownames(py_to_r(ad_object$obs))
# access meta data
meta.data <- py_to_r(ad_object$obs)
meta <- meta.data

### Create a CellChat object using data matrix as input ----

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")

### Add cell information into meta slot of the object ----

cellchat <- addMeta(cellchat, meta = meta, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
