# CellCharDB ----

library(CellChat)

## Create a directory to save data ----

data.dir <- 'CellChatDB'
dir.create(data.dir)
setwd(data.dir)

## Extract ----

write_tsv(CellChatDB.mouse$interaction, file = "CellChatDB.mouse.interaction.tsv")
write_tsv(CellChatDB.human$interaction, file = "CellChatDB.human.interaction.tsv")
