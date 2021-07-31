# Inference and analysis of cell-cell communication using CellChat ----

# package
library(CellChat)
library(NMF)
library(ggplot2)
library(ggalluvial)
library(patchwork)
library(ComplexHeatmap)

# custom
sub <- "LS"

# path
paste0("CellChat/", c("network", "pathway", "LR", "gene", "role", "pattern", "classification")) %>% walk(function(x){
  if(!dir.exists(x)) {dir.create(x)}
})

# continue
cellchat <- readRDS(str_c("CellChat/cellchat.", sub, ".rds"))

## Create a directory to save figures ----

## Load data ----

# Part I: Data input & processing and initialization of CellChat object ----

## Create a CellChat object ----

## Set the ligand-receptor interaction database ----

## Preprocessing the expression data for cell-cell communication analysis ----

# Part II: Inference of cell-cell communication network ----

## Compute the communication probability and infer cellular communication network ----

## Extract the inferred cellular communication network as a data frame ----

## Infer the cell-cell communication at a signaling pathway level ----

## Calculate the aggregated cell-cell communication network ----

groupSize <- as.numeric(table(cellchat@idents))

if(F){
  w <- 9/2.54 # inch
  h <- w * 0.618 * 2
  CairoPDF(
    file = str_c("CellChat/network/circle.", sub, ".pdf"), 
    width = w * 2, height = h
  )
  par(mfrow = c(1,2), xpd = TRUE)
  netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge = F, title.name = "Number of interactions")
  netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge = F, title.name = "Interaction weights/strength")
  graphics.off()
}

w <- 9/2.54 # inch
h <- w * 0.618 * 2
CairoPDF(
  file = str_c("CellChat/network/circle.", sub, ".pdf"), 
  width = w, height = h
)
par(mfrow = c(1,1), xpd = TRUE)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge = F, title.name = "Interaction strength")
graphics.off()

mat <- cellchat@net$weight

if(T){
  w <- 9/2.54 # inch
  h <- w * 0.618 * 2
  CairoPDF(
    file = str_c("CellChat/network/circle.", sub, ".split.pdf"), 
    width = w * 4/1.5, height = h * 2/1.5
  )
  par(mfrow = c(2,4), xpd = TRUE)
  for (i in 1:nrow(mat)) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[, i] <- mat[, i]
    netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
  }
  graphics.off()
}

ht <- netVisual_heatmap(cellchat, measure = "weight", color.heatmap = "Reds")
w <- 8.4/2.54 # inch
h <- 6.8/2.54 # inch
CairoPDF(
  file = str_c("CellChat/network/heatmap.", sub, ".pdf"), 
  width = w, height = h
)
ht
graphics.off()

# Part III: Visualization of cell-cell communication network ----

pathways.show <- cellchat@netP$pathways

write_lines(pathways.show, file = str_c("CellChat/pathway.", sub, ".txt"))

pathways.show %>% str_c("CellChat/pathway/", .) %>% walk(function(x){
  if(!dir.exists(x)) {dir.create(x)}
})

## Visualize each signaling pathway using Hierarchy plot, Circle plot or Chord diagram ----

# Circle plot
w <- 9/2.54 # inch
h <- w * 0.618 * 2
i = 1
for(i in seq_along(pathways.show)){
  x <- pathways.show[i]
  CairoPDF(
    file = str_c("CellChat/pathway/", x, "/circle.", sub, ".pdf"), 
    width = w, height = h
  )
  par(mfrow = c(1,1), xpd = TRUE)
  netVisual_aggregate(cellchat, signaling = x, layout = "circle", signaling.name = paste0(x, "\n"))
  graphics.off()
  rm(x)
}

# Heatmap
w <- 8.4/2.54 # inch
h <- 6.8/2.54 # inch
i = 1
for(i in seq_along(pathways.show)){
  x <- pathways.show[i]
  ht.x <- netVisual_heatmap(cellchat, signaling = x, measure = "weight", color.heatmap = c("#fee5d9", "#fb6a4a", "#a50f15"))
  tryCatch({
    CairoPDF(
      file = str_c("CellChat/network/pathway/", x, "/heatmap.", sub, ".pdf"), 
      width = w, height = h
    )
    ht.x
    graphics.off()
  }, error = function(e){
    cat("ERROR :",conditionMessage(e), "\n")
    file.remove(str_c("CellChat/pathway/", x, "/heatmap.", sub, ".pdf"))
    file.create(str_c("CellChat/pathway/", x, "/heatmap.", sub, ".txt"))
  })
  rm(ht.x)
  rm(x)
}

## Compute the contribution of each ligand-receptor pair to the overall signaling pathway ----

contribution.plot.list <- pathways.show %>% map(~netAnalysis_contribution(cellchat, signaling = .))
w <- 9/2.54 # inch
h <- w * 0.618
list(
  plot = contribution.plot.list, 
  path = str_c("CellChat/pathway/", pathways.show)
) %>% 
  pwalk(
    ggsave, 
    filename = str_c("contribution.", sub, ".pdf"), 
    device = "pdf", 
    width = w, height = h, units = "in"
  )

## Visualize cell-cell communication mediated by a single ligand-receptor pair ----

if(F){
  # Circle plot
  w <- 9/2.54 # inch
  h <- w * 0.618 * 2
  i = 1
  j = 1
  for(i in seq_along(pathways.show)){
    x <- pathways.show[i]
    pairLR <- extractEnrichedLR(cellchat, signaling = x, geneLR.return = FALSE)
    if(F){LR.show <- pairLR[1, ]}
    LR.show <- pairLR
    for(j in seq_along(LR.show)){
      y <- LR.show[j, ]
      CairoPDF(
        file = str_c("CellChat/pathway/", x, "/circle.", y, ".", sub, ".pdf"), 
        width = w, height = h
      )
      par(mfrow = c(1,1), xpd = TRUE)
      netVisual_individual(cellchat, signaling = x, pairLR.use = y, layout = "circle")
      graphics.off()
    }
    rm(x)
  }
}

if(F){
  # Heatmap
  w <- 8.4/2.54 # inch
  h <- 6.8/2.54 # inch
  i = 1
  j = 1
  for(i in seq_along(pathways.show)){
    x <- pathways.show[i]
    pairLR <- extractEnrichedLR(cellchat, signaling = x, geneLR.return = FALSE)
    if(F){LR.show <- pairLR[1, ]}
    LR.show <- pairLR
    for(j in seq_along(LR.show)){
      y <- LR.show[j, ]
      ht.y <- netVisual_heatmap(cellchat, signaling = x, pairLR.use = y, measure = "weight", color.heatmap = c("#fee5d9", "#fb6a4a", "#a50f15"))
      tryCatch({
        CairoPDF(
          file = str_c("CellChat/network/pathway/", x, "/heatmap.", y, ".", sub, ".pdf"), 
          width = w, height = h
        )
        ht.y
        graphics.off()
      }, error = function(e){
        cat("ERROR :",conditionMessage(e), "\n")
        file.remove(str_c("CellChat/pathway/", x, "/heatmap.", y, ".", sub, ".pdf"))
        file.create(str_c("CellChat/pathway/", x, "/heatmap.", y, ".", sub, ".txt"))
      })
      rm(ht.y)
      rm(y)
    }
    rm(x)
  }
}

## Visualize cell-cell communication mediated by multiple ligand-receptors or signaling pathways ----

# Bubble plot 

## Plot the signaling gene expression distribution using violin/dot plot ----

gene.plot.list <- pathways.show %>% map(~plotGeneExpression(
  cellchat, 
  signaling = . # , 
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

# Part IV: Systems analysis of cell-cell communication network----

## Identify signaling roles (e.g., dominant senders, receivers) of cell groups as well as the major contributing signaling ----

### Compute and visualize the network centrality scores ----

# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
w <- 12/2.54 # inch
h <- 5/2.54 # inch
i = 1
for(i in seq_along(pathways.show)){
  x <- pathways.show[i]
  CairoPDF(
    file = str_c("CellChat/pathway/", x, "/role.", sub, ".pdf"), 
    width = w , height = h
  )
  netAnalysis_signalingRole_network(cellchat, signaling = x, width = w, height = h, font.size = 10)
  graphics.off()
  rm(x)
}

### Visualize the dominant senders (sources) and receivers (targets) in a 2D space ----

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg <- netAnalysis_signalingRole_scatter(cellchat, dot.size = c(1, 3))
w <- 8/2.54 # inch
h <- w * 0.75
ggsave(
  filename = str_c("scatter.", sub, ".pdf"), 
  plot = gg, 
  device = "pdf", 
  path = str_c("CellChat/role"), 
  width = w, height = h, units = "in"
)

### Identify signals contributing most to outgoing or incoming signaling of certain cell groups ----

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
w <- 6/2.54 # inch
h <- 18/2.54 # inch
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", width = w, height = h)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", width = w, height = h)
w <- 12/2.54 # inch
h <- 11/2.54 # inch
CairoPDF(
  file = str_c("CellChat/role/heatmap.", sub, ".pdf"), 
  width = w, height = h
)
par(mfrow = c(1,1), xpd = TRUE)
draw(ht1 + ht2, ht_gap = unit(1.5, "cm"))
graphics.off()

## Identify global communication patterns to explore how multiple cell types and signaling pathways coordinate together ----

### Identify and visualize outgoing communication pattern of secreting cells ----

# river plot
w <- 18/2.54 # inch
h <- 10/2.54
ggsave(
  filename = str_c("river.outgoing.", sub, ".pdf"), 
  plot = netAnalysis_river(cellchat, pattern = "outgoing"), 
  device = "pdf", 
  path = str_c("CellChat/pattern"), 
  width = w, height = h, units = "in"
)

### Identify and visualize incoming communication pattern of target cells ----

# river plot
w <- 18/2.54 # inch
h <- 10/2.54
ggsave(
  filename = str_c("river.incoming.", sub, ".pdf"), 
  plot = netAnalysis_river(cellchat, pattern = "incoming"), 
  device = "pdf", 
  path = str_c("CellChat/pattern"), 
  width = w, height = h, units = "in"
)

## Manifold and classification learning analysis of signaling networks ----

### Identify signaling groups based on their functional similarity ----

# Visualization in 2D-space
if(F){
  netVisual_embedding(cellchat, type = "functional", label.size = 3.5)
  netVisual_embeddingZoomIn(cellchat, type = "functional", nCol = 2)
}

gg <- netVisual_embedding(cellchat, type = "functional", color.use = myColor(4, 1), dot.size = c(1, 3), dot.alpha = 0.5, title = "Functional signaling relationships", do.label = F) + 
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

gg <- netVisual_embeddingZoomIn(cellchat, type = "functional", nCol = 2, color.use = myColor(4, 1), dot.size = c(2, 4), label.size = 3.5, dot.alpha = 0.5)
w <- 9.6/2.54 # inch
h <- 10.8/2.54 # inch
ggsave(
  filename = str_c("functionalZoomIn.", sub, ".pdf"), 
  plot = gg, 
  device = "pdf", 
  path = "CellChat/classification", 
  width = w, height = h, units = "in"
)

### Identify signaling groups based on their structure similarity ----

# Part V: Save the CellChat object ----

## Path ----
