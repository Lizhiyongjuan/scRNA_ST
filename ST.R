library(Seurat)
library(tidyverse)
library(Matrix)
library(spacexr)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(reshape2)
library(readr)
library(Seurat)
library(config)
library(ggpubr)
library(gridExtra)
library(reshape2)
library(png)
library(patchwork)
library(SingleR)
library(celldex)
Rcpp::sourceCpp(code='
                
                #include <Rcpp.h>
                
                using namespace Rcpp;
                
                // [[Rcpp::export]]
                
                IntegerMatrix asMatrix(NumericVector rp,
                                       
                                       NumericVector cp,
                                       
                                       NumericVector z,
                                       
                                       int nrows,
                                       
                                       int ncols){
                  
                  int k = z.size() ;
                  
                  IntegerMatrix  mat(nrows, ncols);
                  
                  for (int i = 0; i < k; i++){
                    
                    mat(rp[i],cp[i]) = z[i];
                    
                  }
                  
                  return mat;
                  
                }
                
                ' )
as_matrix <- function(mat){
  row_pos <- mat@i
  col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])
  tmp <- asMatrix(rp = row_pos, cp = col_pos, z = mat@x,
                  nrows =  mat@Dim[1], ncols = mat@Dim[2])
  row.names(tmp) <- mat@Dimnames[[1]]
  colnames(tmp) <- mat@Dimnames[[2]]
  return(tmp)
}

library(spacexr)
library(tidyverse)
library(Matrix)
setwd("XXX")
scRNA=readRDS('XXX.rds')
table(scRNA$celltype)

sc_counts <- as.matrix(scRNA[['RNA']]@counts)

sc_nUMI = colSums(sc_counts)

head(sc_nUMI)

cellType=data.frame(barcode=colnames(scRNA), celltype=scRNA$celltype)

names(cellType) = c('barcode', 'cell_type')

cell_types = cellType$cell_type; names(cell_types) <- cellType$barcode # create cell_types named list

cell_types =as.factor(cell_types) # convert to factor data type
cell_types

## construct referrence scRNA
reference = Reference(sc_counts, cell_types, sc_nUMI)

rm(sc_counts)
rm(scRNA)
gc()


# spatial --------------------------
stRNA=readRDS('XXX.RDS')
coords <- GetTissueCoordinates(stRNA,
                               cols = c("row", "col"),
                               scale = NULL) 
colnames(coords)=c('xcoord','ycoord')
sp_counts <- as.matrix(stRNA[['Spatial']]@counts)
#nUMI
sp_nUMI <- colSums(sp_counts)
# construct spatial 
puck <- SpatialRNA(coords, sp_counts, sp_nUMI)
rm(stRNA)
rm(sp_counts)
gc()
### attention: each celltype should contain more than 25 cells
## if error, go back to delete 
myRCTD <- create.RCTD(puck, reference)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet') #run.RCTD
save(myRCTD,file ='RCTD_p4.Rdata')
load('RCTD_p4.Rdata')
anno=myRCTD@results$results_df
anno=anno[,'first_type',drop=F]
anno2=as.matrix(myRCTD@results$weights)
anno2=as.data.frame(anno2)
anno2$celltype=anno$first_type

#stlearn
write.csv(anno2,file ='RCTD_decon_mat.csv')
table(myRCTD@internal_vars$gene_list_bulk)
results <- myRCTD@results
results_df <- results$results_df

weight=myRCTD@results$weights
weight=as.data.frame(weight)

#
stRNA=readRDS('xxx.RDS')
SpatialFeaturePlot(stRNA, features = c("Bmp10","S100a10","Abcc9","Lars2","Sfrp1","Ube2c"))
SpatialFeaturePlot(stRNA, features = c("Bmp10","Gja5","Myh6","Nppa","Tbx5"))
SpatialFeaturePlot(stRNA, features = c('Nr2f1','Irx4','Vsnl1'))
SpatialFeaturePlot(stRNA, features = c('Col1a1','Pdgfra','Dcn',
                                       'Serpinb2','Cd5l','Timd4',
                                       'Lrrn4','Bnc1','Msln',    
                                       'Kcna1','Plp1','Gfra3'   
))
SpatialFeaturePlot(stRNA, features = c('Uts2b'))

library(Seurat)
stRNA=stRNA[,rownames(results_df)]
stRNA$celltype1=results_df$first_type
stRNA$celltype2=results_df$second_type
Idents(stRNA)=stRNA$celltype1
coor <- read.csv('XXX.csv')
stRNA@meta.data$num <- colnames(stRNA)
stRNA@meta.data <- merge(stRNA@meta.data,coor,by='num')
rownames(stRNA@meta.data) <- stRNA@meta.data[,1]
saveRDS(stRNA,'XXX.rds')
table(stRNA$imagecol)
min(stRNA$imagecol)
max(stRNA$imagecol)
table(stRNA$imagerow)
min(stRNA$imagerow)
max(stRNA$imagerow)
stRNA1 <- stRNA
stRNA <- stRNA1
stRNA <- subset(stRNA, subset = imagecol > 0 &imagecol < 6500 & 
                  imagerow>0 & imagerow<7000)

plot <- SpatialDimPlot(stRNA,pt.size.factor = 4)
plot
stRNA <- stRNA1
min(stRNA$imagecol)
max(stRNA$imagecol)
table(stRNA$imagerow)
min(stRNA$imagerow)
max(stRNA$imagerow)
stRNA <- subset(stRNA, subset = imagecol > 6500 &imagecol < 12363 & 
                  imagerow>0 & imagerow<7000)

plot <- SpatialDimPlot(stRNA,pt.size.factor = 4)
plot
stRNA <- stRNA1
min(stRNA$imagecol)
max(stRNA$imagecol)
table(stRNA$imagerow)
min(stRNA$imagerow)
max(stRNA$imagerow)
stRNA <- subset(stRNA, subset = imagecol > 0 &imagecol < 6500 & 
                  imagerow>7000 & imagerow<13216)

plot <- SpatialDimPlot(stRNA,pt.size.factor = 4)
plot
table(stRNA$celltype1)
stRNA <- stRNA1
min(stRNA$imagecol)
max(stRNA$imagecol)
table(stRNA$imagerow)
min(stRNA$imagerow)
max(stRNA$imagerow)
stRNA <- subset(stRNA, subset = imagecol > 6500 &imagecol < 12363 & 
                  imagerow>7000 & imagerow<13216)

plot <- SpatialDimPlot(stRNA,pt.size.factor = 4)
plot
####
stRNA=readRDS('stRNA_x34.RDS')
library(Seurat)
stRNA=stRNA[,rownames(results_df)]
stRNA$celltype1=results_df$first_type
stRNA$celltype2=results_df$second_type
Idents(stRNA)=stRNA$celltype1
coor <- read.csv('XXX.csv')
stRNA@meta.data$num <- colnames(stRNA)
stRNA@meta.data <- merge(stRNA@meta.data,coor,by='num')
rownames(stRNA@meta.data) <- stRNA@meta.data[,1]
saveRDS(stRNA,'XXX.rds')
table(stRNA$imagecol)
min(stRNA$imagecol)
max(stRNA$imagecol)
table(stRNA$imagerow)
min(stRNA$imagerow)
max(stRNA$imagerow)
stRNA1 <- stRNA
stRNA <- stRNA1
stRNA <- subset(stRNA, subset = imagecol > 0 &imagecol < 6500 & 
                  imagerow>0 & imagerow<7000)

plot <- SpatialDimPlot(stRNA,pt.size.factor = 4)
plot
stRNA <- stRNA1
min(stRNA$imagecol)
max(stRNA$imagecol)
table(stRNA$imagerow)
min(stRNA$imagerow)
max(stRNA$imagerow)
stRNA <- subset(stRNA, subset = imagecol > 6500 &imagecol < 12363 & 
                  imagerow>0 & imagerow<7000)

plot <- SpatialDimPlot(stRNA,pt.size.factor = 4)
stRNA <- stRNA1
min(stRNA$imagecol)
max(stRNA$imagecol)
table(stRNA$imagerow)
min(stRNA$imagerow)
max(stRNA$imagerow)
stRNA <- subset(stRNA, subset = imagecol > 0 &imagecol < 6500 & 
                  imagerow>7000 & imagerow<13216)

plot <- SpatialDimPlot(stRNA,pt.size.factor = 4)
plot
table(stRNA$celltype1)
stRNA <- stRNA1
min(stRNA$imagecol)
max(stRNA$imagecol)
table(stRNA$imagerow)
min(stRNA$imagerow)
max(stRNA$imagerow)
stRNA <- subset(stRNA, subset = imagecol > 6500 &imagecol < 12363 & 
                  imagerow>7000 & imagerow<13216)

plot <- SpatialDimPlot(stRNA,pt.size.factor = 4)
plot
####
setwd("XXX")
scRNA <- readRDS('XXX.rds')
library(CellChat)
gc()
table(scRNA@meta.data$orig.ident)
scRNA@meta.data$tissue_type <- scRNA@meta.data$celltype1
table(scRNA@meta.data$celltype1)
scRNA_chat <- scRNA
scRNA_chat@meta.data$celltype1 <- paste(scRNA_chat@meta.data$celltype1,"1",sep="-")
meta =scRNA_chat@meta.data # a dataframe with rownames containing cell mata data
data_input <- as.matrix(scRNA_chat@assays$Spatial$data)
#data_input <- as.matrix(scRNA_chat@assays$RNA@data)
#data_input=data_input[,rownames(meta)]
identical(colnames(data_input),rownames(meta))

library(CellChat)
cellchat <- createCellChat(object = data_input, meta = meta, group.by = "celltype1")

CellChatDB <- CellChatDB.mouse  
groupSize <- as.numeric(table(cellchat@idents))
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
cellchat@DB <- CellChatDB.use 

dplyr::glimpse(CellChatDB$interaction)
cellchat <- subsetData(cellchat)

cellchat <- identifyOverExpressedGenes(cellchat)

cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- projectData(cellchat, PPI.human)
unique(cellchat@idents)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 1)
cellchat <- computeCommunProbPathway(cellchat)
df.net<- subsetCommunication(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
dev.off()
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= T, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()
table(scRNA$celltype1)          
p_bubble= netVisual_bubble(cellchat,
                            source.use = c('Lars2*-Cardiomyocyte-1'),
                           remove.isolate = FALSE)+coord_flip()
p_bubble
p_bubble= netVisual_bubble(cellchat,
                           #targets.use = c('IS_Mono','Other_Mono'),
                           #targets.use = c('Astrocytes','Epithelial cells', 'Microglia','Monocytes','Oligodendrocytes'),
                           targets.use = c('Abcc9*-Cardiomyocyte-1'),
                           remove.isolate = FALSE)+coord_flip()
p_bubble
p_bubble= netVisual_bubble(cellchat,
                           #targets.use = c('IS_Mono','Other_Mono'),
                           #targets.use = c('Astrocytes','Epithelial cells', 'Microglia','Monocytes','Oligodendrocytes'),
                           targets.use = c('Lars2*-Cardiomyocyte'),
                           remove.isolate = FALSE)+coord_flip()
p_bubble
ggsave("Lars2t.png", plot = p_bubble, width = 9, height = 7)
ggsave("Lars2t.pdf", plot = p_bubble, width = 9, height = 7)

p_bubble= netVisual_bubble(cellchat,
                           #targets.use = c('IS_Mono','Other_Mono'),
                           #targets.use = c('Astrocytes','Epithelial cells', 'Microglia','Monocytes','Oligodendrocytes'),
                           targets.use = c('S100a10*-Cardiomyocyte'),
                           remove.isolate = FALSE)+coord_flip()
p_bubble
ggsave("S100a10t.png", plot = p_bubble, width = 9, height = 7)
ggsave("S100a10t.pdf", plot = p_bubble, width = 9, height = 7)

p_bubble= netVisual_bubble(cellchat,
                           #targets.use = c('IS_Mono','Other_Mono'),
                           #targets.use = c('Astrocytes','Epithelial cells', 'Microglia','Monocytes','Oligodendrocytes'),
                           targets.use = c('Ube2c*-Cardiomyocyte'),
                           remove.isolate = FALSE)+coord_flip()
p_bubble
ggsave("Ube2ct.png", plot = p_bubble, width = 9, height = 7)
ggsave("Ube2ct.pdf", plot = p_bubble, width = 9, height = 7)


p_bubble= netVisual_bubble(cellchat,
                           #targets.use = c('IS_Mono','Other_Mono'),
                           #targets.use = c('Astrocytes','Epithelial cells', 'Microglia','Monocytes','Oligodendrocytes'),
                           targets.use = c('Sfrp1*-Cardiomyocyte'),
                           remove.isolate = FALSE)+coord_flip()
p_bubble
ggsave("Sfrp1t.png", plot = p_bubble, width = 9, height = 7)
ggsave("Sfrp1t.pdf", plot = p_bubble, width = 9, height = 7)

p_bubble= netVisual_bubble(cellchat,
                           #targets.use = c('IS_Mono','Other_Mono'),
                           #targets.use = c('Astrocytes','Epithelial cells', 'Microglia','Monocytes','Oligodendrocytes'),
                           targets.use = c('Abcc9*-Cardiomyocyte'),
                           remove.isolate = FALSE)+coord_flip()
p_bubble
ggsave("Abcc9t.png", plot = p_bubble, width = 9, height = 7)
ggsave("Abcc9t.pdf", plot = p_bubble, width = 9, height = 7)

p_bubble= netVisual_bubble(cellchat,
                           #targets.use = c('IS_Mono','Other_Mono'),
                           #targets.use = c('Astrocytes','Epithelial cells', 'Microglia','Monocytes','Oligodendrocytes'),
                           targets.use = c('Bmp10*-Cardiomyocyte'),
                           remove.isolate = FALSE)+coord_flip()
p_bubble
ggsave("Bmp10t.png", plot = p_bubble, width = 9, height = 7)
ggsave("Bmp10t.pdf", plot = p_bubble, width = 9, height = 7)



p_bubble= netVisual_bubble(cellchat,
                           #targets.use = c('IS_Mono','Other_Mono'),
                           #targets.use = c('Astrocytes','Epithelial cells', 'Microglia','Monocytes','Oligodendrocytes'),
                           sources.use = c('Lars2*-Cardiomyocyte'),
                           remove.isolate = FALSE)+coord_flip()
p_bubble
ggsave("Lars2.png", plot = p_bubble, width = 9, height = 7)
ggsave("Lars2.pdf", plot = p_bubble, width = 9, height = 7)

p_bubble= netVisual_bubble(cellchat,
                           #targets.use = c('IS_Mono','Other_Mono'),
                           #targets.use = c('Astrocytes','Epithelial cells', 'Microglia','Monocytes','Oligodendrocytes'),
                           sources.use = c('S100a10*-Cardiomyocyte'),
                           remove.isolate = FALSE)+coord_flip()
p_bubble
ggsave("S100a10.png", plot = p_bubble, width = 9, height = 7)
ggsave("S100a10.pdf", plot = p_bubble, width = 9, height = 7)

p_bubble= netVisual_bubble(cellchat,
                           #targets.use = c('IS_Mono','Other_Mono'),
                           #targets.use = c('Astrocytes','Epithelial cells', 'Microglia','Monocytes','Oligodendrocytes'),
                           sources.use = c('Ube2c*-Cardiomyocyte'),
                           remove.isolate = FALSE)+coord_flip()
p_bubble
ggsave("Ube2c.png", plot = p_bubble, width = 9, height = 7)
ggsave("Ube2c.pdf", plot = p_bubble, width = 9, height = 7)


p_bubble= netVisual_bubble(cellchat,
                           #targets.use = c('IS_Mono','Other_Mono'),
                           #targets.use = c('Astrocytes','Epithelial cells', 'Microglia','Monocytes','Oligodendrocytes'),
                           sources.use = c('Sfrp1*-Cardiomyocyte'),
                           remove.isolate = FALSE)+coord_flip()
p_bubble
ggsave("Sfrp1.png", plot = p_bubble, width = 9, height = 7)
ggsave("Sfrp1.pdf", plot = p_bubble, width = 9, height = 7)