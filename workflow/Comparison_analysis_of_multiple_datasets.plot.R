# Comparison analysis of multiple datasets using CellChat ----

# package
library(CellChat)
library(NMF)
library(ggplot2)
library(ggalluvial)
library(patchwork)
library(ComplexHeatmap)

# custom
sub1 <- "NL"
sub2 <- "LS"
sub <- paste0(sub1, "_", sub2)

# path
paste0("CellChat/", c("network", "pathway", "LR", "gene", "role", "pattern", "classification")) %>% walk(function(x){
  if(!dir.exists(x)) {dir.create(x)}
})
paste0("CellChat/", sub) %>% walk(function(x){
  if(!dir.exists(x)) {dir.create(x)}
})

# continue
cellchat <- readRDS(str_c("CellChat/cellchat.", sub, ".rds"))

## Create a directory to save figures ----

## Load data ----

# Part I: Predict general principles of cell-cell communication ----

## Compare the total number of interactions and interaction strength ----

if(F){
  gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), color.use = myColor(2,1))
  gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), color.use = myColor(2,1), measure = "weight")
  gg <- gg1 + gg2
  w <- 8/2.54 # inch
  h <- w
  ggsave(
    filename = "compare.pdf", 
    plot = gg, 
    device = pdf, 
    path = str_c("CellChat/", sub), 
    scale = 1,
    width = w * 2, height = h, units = "in"
  )
}

gg <- compareInteractions(cellchat, show.legend = F, group = c(1,2), color.use = myColor(2,1), measure = "weight", size.text = 12)
w <- 8/2.54 # inch
h <- w
ggsave(
  filename = "compare.pdf", 
  plot = gg, 
  device = pdf, 
  path = str_c("CellChat/", sub), 
  scale = 1,
  width = w, height = h, units = "in"
)

## Compare the number of interactions and interaction strength among different cell populations ----

### Differential number of interactions or interaction strength among different cell populations ----

if(F){
  w <- 9/2.54 # inch
  h <- w * 0.618 * 2
  CairoPDF(
    file = str_c("CellChat/", sub, "/", "diff.pdf"), 
    width = w * 2, height = h
  )
  par(mfrow = c(1,2), xpd = TRUE)
  netVisual_diffInteraction(cellchat, weight.scale = T)
  netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
  # red (or blue) colored edges represent increased (or decreased) signaling in the second dataset compared to the first one.
  graphics.off()
}

w <- 9/2.54 # inch
h <- w * 0.618 * 2
CairoPDF(
  file = str_c("CellChat/", sub, "/", "diff.pdf"), 
  width = w, height = h
)
par(mfrow = c(1,1), xpd = TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
# red (or blue) colored edges represent increased (or decreased) signaling in the second dataset compared to the first one.
graphics.off()

groupSize.list <- object.list %>% map(~as.numeric(table(.@idents)))

w <- 9/2.54 # inch
h <- w * 0.618 * 2
CairoPDF(
  file = str_c("CellChat/network/circle.", sub, ".pdf"), 
  width = w * 2, height = h
)
weight.max <- getMaxWeight(object.list, attribute = c("idents", "weight"))
par(mfrow = c(1,2), xpd = TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$weight, vertex.weight = groupSize.list[[i]], weight.scale = T, label.edge = F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = names(object.list)[i])
}
graphics.off()

if(F){
  ht1 <- netVisual_heatmap(cellchat)
  ht2 <- netVisual_heatmap(cellchat, measure = "weight")
  # The top colored bar plot represents the sum of column of values displayed in the heatmap (incoming signaling). 
  # The right colored bar plot represents the sum of row of values (outgoing signaling). 
  # red  (or blue) represents increased (or decreased) signaling in the second dataset compared to the first one.
  w <- 13.2/2.54 # inch
  h <- 6.8/2.54 # inch
  CairoPDF(
    file = str_c("CellChat/", sub, "heatmap.pdf"), 
    width = w, height = h
  )
  par(mfrow = c(1,1), xpd = TRUE)
  draw(ht1 + ht2, ht_gap = unit(0.75, "cm"))
  graphics.off()
}

ht <- netVisual_heatmap(cellchat, measure = "weight")
# The top colored bar plot represents the sum of column of values displayed in the heatmap (incoming signaling). 
# The right colored bar plot represents the sum of row of values (outgoing signaling). 
# red  (or blue) represents increased (or decreased) signaling in the second dataset compared to the first one.
w <- 8.8/2.54 # inch
h <- 6.8/2.54 # inch
CairoPDF(
  file = str_c("CellChat/", sub, "/heatmap.pdf"), 
  width = w, height = h
)
par(mfrow = c(1,1), xpd = TRUE)
ht
graphics.off()

ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], measure = "weight", color.heatmap = "Reds", title.name = names(object.list)[i])
}
w <- 13.2/2.54 # inch
h <- 6.8/2.54 # inch
CairoPDF(
  file = str_c("CellChat/network/heatmap.", sub, ".pdf"), 
  width = w, height = h
)
draw(ht[[1]] + ht[[2]], ht_gap = unit(0.75, "cm"))
graphics.off()

### Differential number of interactions or interaction strength among different cell types ----

## Compare the major sources and targets in 2D space ----

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], dot.size = c(1, 3), title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
gg <- gg[[1]] + gg[[2]] +
  plot_layout(guides = 'collect')
w <- 8/2.54 # inch
h <- w * 0.75
ggsave(
  filename = str_c("scatter.", sub, ".pdf"), 
  plot = gg, 
  device = "pdf", 
  path = str_c("CellChat/role"), 
  width = w * 2, height = h, units = "in"
)

# Part II: Identify the conserved and context-specific signaling pathways ----

## Identify signaling networks with larger (or less) difference as well as signaling groups based on their functional/structure similarity ----

### Identify signaling groups based on their functional similarity ----

# Visualization in 2D-space
if(F){
  netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
  netVisual_embeddingPairwiseZoomIn(cellchat, type = "functional", nCol = 2)
}

gg <- netVisual_embeddingPairwise(cellchat, type = "functional", color.use = myColor(4, 1), dot.size = c(1, 3), dot.alpha = 0.5, title = "Functional signaling relationships", do.label = F) + 
  theme(
    legend.position = "left"
  )
w <- 10.3/2.54 # inch
h <- 7.2/2.54 # inch
ggsave(
  filename = str_c("functional.", sub, ".pdf"), 
  plot = gg, 
  device = "pdf", 
  path = "CellChat/classification", 
  width = w, height = h, units = "in"
)

gg <- netVisual_embeddingPairwiseZoomIn(cellchat, type = "functional", nCol = 2, color.use = myColor(4, 1), dot.size = c(2, 4), label.size = 3.5, dot.alpha = 0.5)
w <- 9.6/2.54 # inch
h <- 10.8/2.54 # inch
ggsave(
  filename = str_c("functionalZoomIn.", sub, ".pdf"), 
  plot = gg, 
  device = "pdf", 
  path = "CellChat/classification", 
  width = w, height = h, units = "in"
)

### Identify signaling groups based on structure similarity ----

### Compute and visualize the pathway distance in the learned joint manifold ----

gg <- rankSimilarity(cellchat, type = "functional")
w <- 5/2.54 # inch
h <- 7/2.54
ggsave(
  filename = "rankSimilarity.pdf", 
  plot = gg, 
  device = pdf, 
  path = str_c("CellChat/", sub), 
  width = w, height = h, units = "in"
)

## Identify and visualize the conserved and context-specific signaling pathways ----

### Compare the overall information flow of each signaling pathway ----

gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE, color.use = myColor(2,1)) + 
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5), 
    legend.position = "top"
  )
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE, color.use = myColor(2,1)) + 
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5), 
    legend.position = "top"
  )
gg <- gg1 + gg2
w <- 7/2.54 # inch
h <- 12/2.54
ggsave(
  filename = "rankNet.pdf", 
  plot = gg, 
  device = pdf, 
  path = str_c("CellChat/", sub), 
  width = w * 2, height = h, units = "in"
)

### Compare outgoing (or incoming) signaling associated with each cell population ----

pathway.union <- union(cellchat1@netP$pathways, cellchat2@netP$pathways)

w <- 2.5
h <- 9
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = paste0("\n", names(object.list)[i]), width = w, height = h, color.heatmap = "BuGn")
}
w <- 12/2.54 # inch
h <- 13.5/2.54 # inch
CairoPDF(
  file = str_c("CellChat/role/heatmap.outgoing.", sub, ".pdf"), 
  width = w, height = h
)
draw(ht[[1]] + ht[[2]], ht_gap = unit(1.5, "cm"))
graphics.off()

w <- 2.5
h <- 9
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = paste0("\n", names(object.list)[i]), width = w, height = h, color.heatmap = "GnBu")
}
w <- 12/2.54 # inch
h <- 13.5/2.54 # inch
CairoPDF(
  file = str_c("CellChat/role/heatmap.incoming.", sub, ".pdf"), 
  width = w, height = h
)
draw(ht[[1]] + ht[[2]], ht_gap = unit(1.5, "cm"))
graphics.off()

# Part III: Identify the upgulated and down-regulated signaling ligand-receptor pairs ----

## Identify dysfunctional signaling by comparing the communication probabities ----

# Bubble plot

## Identify dysfunctional signaling by using differential expression analysis ----

# Bubble plot

# Part IV: Visually compare cell-cell communication using Hierarchy plot, Circle plot or Chord diagram ----

pathways.show <- pathway.union

# Circle plot
w <- 9/2.54 # inch
h <- w * 0.618 * 2
i = 1
j = 1
for(i in seq_along(pathways.show)){
  x <- pathways.show[i]
  tryCatch({
    CairoPDF(
      file = str_c("CellChat/pathway/", x, "/circle.", sub, ".pdf"), 
      width = w * 2, height = h
    )
    weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = x) # control the edge weights across different datasets
    par(mfrow = c(1,2), xpd = TRUE)
    for (j in 1:length(object.list)) {
      netVisual_aggregate(object.list[[j]], signaling = x, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste0(x, " - ", names(object.list)[j], "\n"))
    }
    graphics.off()
  }, error = function(e){
    cat("ERROR :",conditionMessage(e), "\n")
    graphics.off()
    file.remove(str_c("CellChat/pathway/", x, "/circle.", sub, ".pdf"))
    file.create(str_c("CellChat/pathway/", x, "/circle.", sub, ".txt"))
  })
  rm(x)
}

# Heatmap
w <- 13.8/2.54 # inch
h <- 6.8/2.54 # inch
i = 1
j = 1
for(i in seq_along(pathways.show)){
  x <- pathways.show[i]
  ht.x <- list()
  tryCatch({
    for (j in 1:length(object.list)) {
      ht.x[[j]] <- netVisual_heatmap(object.list[[j]], signaling = x, color.heatmap = c("#fee5d9", "#fb6a4a", "#a50f15"),title.name = paste0(x, " - ",names(object.list)[j]))
    }
    CairoPDF(
      file = str_c("CellChat/pathway/", x, "/heatmap.", sub, ".pdf"), 
      width = w, height = h
    )
    draw(ht.x[[1]] + ht.x[[2]], ht_gap = unit(0.85, "cm"))
    graphics.off()
  }, error = function(e){
    cat("ERROR :",conditionMessage(e), "\n")
    graphics.off()
    file.remove(str_c("CellChat/pathway/", x, "/heatmap.", sub, ".pdf"))
    file.create(str_c("CellChat/pathway/", x, "/heatmap.", sub, ".txt"))
  })
  rm(ht.x)
  rm(x)
}

# Part V: Compare the signaling gene expression distribution between different datasets ----

cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = names(object.list)) # set factor level

gene.plot.list <- pathways.show %>% map(~plotGeneExpression(
  cellchat, 
  signaling = ., 
  split.by = "datasets", 
  colors.ggplot = T, 
  color.use = myColor(2,1)
  # enriched.only = FALSE
))
w <- 9/2.54 # inch
h <- w * 0.618
list(
  plot = gene.plot.list, 
  path = str_c("CellChat/pathway/", pathways.show)
) %>% 
  pwalk(
    ggsave, 
    filename = str_c("gene.", sub, ".pdf"), 
    device = "pdf", 
    width = w, height = h, units = "in"
  )

# Save the merged CellChat object ----
