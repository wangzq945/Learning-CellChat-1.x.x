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
c("CellChat") %>% walk(function(x){if(!dir.exists(x)) {dir.create(x)}})

## Create a directory to save figures ----

data.dir <- str_c("CellChat/figure", ".", sub)
dir.create(data.dir)
setwd(data.dir)

## Load data ----

load("../../data/data_humanSkin_CellChat.rda")

# Part I: Data input & processing and initialization of CellChat object ----

data.input = data_humanSkin$data # normalized data matrix
meta = data_humanSkin$meta # a dataframe with rownames containing cell mata data

if(T){
  # Prepare input data for CelChat analysis
  cell.use = rownames(meta)[meta$condition == sub] # extract the cell names
  data.input = data.input[, cell.use]
  meta = meta[cell.use, ]
}

unique(meta$labels) # check the cell labels

## Create a CellChat object ----

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")

if(F){
  # Add cell information into meta slot of the object
  cellchat <- addMeta(cellchat, meta = meta)
  # Set default cell identity
  cellchat <- setIdent(cellchat, ident.use = "labels")
}

levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

## Set the ligand-receptor interaction database ----

CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
if(F){
  # use all CellChatDB for cell-cell communication analysis
  CellChatDB.use <- CellChatDB # simply use the default CellChatDB
}

# set the used database in the object
cellchat@DB <- CellChatDB.use

## Preprocessing the expression data for cell-cell communication analysis ----

cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
future::plan("multiprocess", workers = 4) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)

# Part II: Inference of cell-cell communication network ----

## Compute the communication probability and infer cellular communication network ----

cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
# When analyzing unsorted single-cell transcriptomes, under the assumption that abundant cell populations tend to send collectively stronger signals than the rare cell populations, CellChat can also consider the effect of cell proportion in each cell group in the probability calculation. 
# USER can set population.size = TRUE.

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

## Extract the inferred cellular communication network as a data frame ----

if(F){
  df.net <- subsetCommunication(cellchat, slot.name = "net")
  df.net <- subsetCommunication(cellchat, slot.name = "netP")
  df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5))
  df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb"))
}

## Infer the cell-cell communication at a signaling pathway level ----

cellchat <- computeCommunProbPathway(cellchat)

## Calculate the aggregated cell-cell communication network ----

cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))

par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

mat <- cellchat@net$weight

par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

# Part III: Visualization of cell-cell communication network ----

pathways.show <- c("CXCL") 

## Visualize each signaling pathway using Hierarchy plot, Circle plot or Chord diagram ----

# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)

# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
#> Note: The first link end is drawn out of sector 'Inflam. FIB'.

# Chord diagram
group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4)) # grouping cell clusters into fibroblast, DC and TC cells
names(group.cellType) <- levels(cellchat@idents)
netVisual_chord_cell(cellchat, signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network"))
#> Plot the aggregated cell-cell communication network at the signaling pathway level
#> Note: The first link end is drawn out of sector 'Inflam. FIB'.

## Compute the contribution of each ligand-receptor pair to the overall signaling pathway ----

netAnalysis_contribution(cellchat, signaling = pathways.show)

## Visualize cell-cell communication mediated by a single ligand-receptor pair ----

pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)

LR.show <- pairLR.CXCL[1,] # show one ligand-receptor pair

# Hierarchy plot
vertex.receiver = seq(1,4) # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)

# Circle plot
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# Chord diagram
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")

## Visualize cell-cell communication mediated by multiple ligand-receptors or signaling pathways ----

### Bubble plot ----

# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), remove.isolate = FALSE)

# show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), signaling = c("CCL","CXCL"), remove.isolate = FALSE)

# show all the significant interactions (L-R pairs) based on user's input (defined by `pairLR.use`)
pairLR.use <- extractEnrichedLR(cellchat, signaling = c("CCL","CXCL","FGF"))
netVisual_bubble(cellchat, sources.use = c(3,4), targets.use = c(5:8), pairLR.use = pairLR.use, remove.isolate = TRUE)

### Chord diagram ----

# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')

# show all the interactions sending from Inflam.FIB
netVisual_chord_gene(cellchat, sources.use = 4, targets.use = c(5:11), lab.cex = 0.5,legend.pos.y = 30)

# show all the interactions received by Inflam.DC
netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4), targets.use = 8, legend.pos.x = 15)

# show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4), targets.use = c(5:11), signaling = c("CCL","CXCL"),legend.pos.x = 8)

# show all the significant signaling pathways from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4), targets.use = c(5:11), slot.name = "netP", legend.pos.x = 10)

## Plot the signaling gene expression distribution using violin/dot plot ----

plotGeneExpression(cellchat, signaling = "CXCL")

plotGeneExpression(cellchat, signaling = "CXCL", enriched.only = FALSE)

# Part IV: Systems analysis of cell-cell communication network----

## Identify signaling roles (e.g., dominant senders, receivers) of cell groups as well as the major contributing signaling ----

### Compute and visualize the network centrality scores ----

# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways

# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

### Visualize the dominant senders (sources) and receivers (targets) in a 2D space ----

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat)

# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("CXCL", "CCL"))
gg1 + gg2

### Identify signals contributing most to outgoing or incoming signaling of certain cell groups ----

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2

# Signaling role analysis on the cell-cell communication networks of interest
ht <- netAnalysis_signalingRole_heatmap(cellchat, signaling = c("CXCL", "CCL"))

## Identify global communication patterns to explore how multiple cell types and signaling pathways coordinate together ----

### Identify and visualize outgoing communication pattern of secreting cells ----

selectK(cellchat, pattern = "outgoing")

nPatterns = 3
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)

# river plot
netAnalysis_river(cellchat, pattern = "outgoing")

# dot plot
netAnalysis_dot(cellchat, pattern = "outgoing")

### Identify and visualize incoming communication pattern of target cells ----

selectK(cellchat, pattern = "incoming")

nPatterns = 4
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)

# river plot
netAnalysis_river(cellchat, pattern = "incoming")

# dot plot
netAnalysis_dot(cellchat, pattern = "incoming")

## Manifold and classification learning analysis of signaling networks ----

### Identify signaling groups based on their functional similarity ----

cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")

# Visualization in 2D-space
netVisual_embedding(cellchat, type = "functional", label.size = 3.5)
netVisual_embeddingZoomIn(cellchat, type = "functional", nCol = 2)

### Identify signaling groups based on their structure similarity ----

cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural")

# Visualization in 2D-space
netVisual_embedding(cellchat, type = "structural", label.size = 3.5)
netVisual_embeddingZoomIn(cellchat, type = "structural", nCol = 2)

# Part V: Save the CellChat object ----

saveRDS(cellchat, file = str_c("../cellchat.", sub, ".rds"))

## Path ----

setwd("../../")
