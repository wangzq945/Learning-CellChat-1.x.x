# Comparison analysis of multiple datasets with different cell type compositions ----

## Load the required libraries ----

library(CellChat)
library(NMF)
library(ggplot2)
library(ggalluvial)
library(patchwork)
library(ComplexHeatmap)

# custom
sub1 <- "E13"
sub2 <- "E14"
sub <- paste0(sub1, "_", sub2)

## Create a directory to save figures ----

data.dir <- str_c("CellChat/figure", ".", sub)
dir.create(data.dir)
setwd(data.dir)

## Load data ----

cellchat.E14 <- updateCellChat(cellchat.E14)

cellchat1 <- readRDS(str_c("CellChat/cellchat_embryonic_", sub1, ".rds"))
cellchat1 <- updateCellChat(cellchat.1)
cellchat2 <- readRDS(str_c("CellChat/cellchat_embryonic_", sub2, ".rds"))
cellchat2 <- updateCellChat(cellchat.2)
object.list <- list(cellchat1, cellchat2)
names(object.list) <- c(sub1, sub2)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

# Part I: Comparison analysis of multiple datasets with slightly different cell type compositions ----

## Lift up CellChat object and merge together ----

# Define the cell labels to lift up
group.new = levels(cellchat.E14@idents)
cellchat.E13 <- liftCellChat(cellchat.E13, group.new)
object.list <- list(E13 = cellchat.E13, E14 = cellchat.E14)
cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)

## Visualize the inferred signaling network using the lifted object ----

# Hierarchy plot
pathways.show <- c("WNT") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
vertex.receiver = seq(1,10) # Left portion of hierarchy plot the shows signaling to dermal cells and right portion shows signaling to epidermal cells
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, vertex.receiver = vertex.receiver, edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

# Circle plot
pathways.show <- c("WNT") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

# Chord diagram
pathways.show <- c("WNT") 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
}

# Chord diagram
group.merged <- c(rep("Dermal", 10), rep("Epidermal", 3)) # grouping cell clusters into dermal and epidermal cells to study the cell-cell communication between dermal and epidermal
names(group.merged) <- levels(object.list[[1]]@idents)
pathways.show <- c("WNT") 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_cell(object.list[[i]], signaling = pathways.show, group = group.merged, title.name = paste0(pathways.show, " signaling network - ", names(object.list)[i]))
}

# Part II: Comparison analysis of multiple datasets with vastly different cell type compositions ----
