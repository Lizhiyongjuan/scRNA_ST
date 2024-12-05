library(RColorBrewer) 
library(viridis)
library(wesanderson)
library(Seurat)
library(ggplot2)
setwd("XXXX")
scRNA <- readRDS('XXXX.RDS')
scRNA_CM1=subset(scRNA,celltype=='Ventricular cell')
scRNA_CM2=subset(scRNA,celltype=='Atrial cell')
scRNA_CM3=subset(scRNA,celltype=='Unknow')
scRNA_CM  <-  merge(scRNA_CM1, y = c(scRNA_CM2,scRNA_CM3),  project = "SC")
FeaturePlot(scRNA_CM,features = 'Abcc9',order = T,split.by = 'orig.ident')
gc()
setwd('XXXX')
samples=list.files()
samples=samples[1:1]
library(Seurat)
for (i in seq_along(samples)[1:1]){
  assign(paste0("scRNA", samples[i]), Read10X(data.dir = samples[i]))
}
for (i in seq_along(samples)[1:1]){
  assign(paste0("scRNA", samples[i]), CreateSeuratObject(counts = eval(parse(text = paste0("scRNA", samples[i]))), project = samples[i]))
}
scRNA <- scRNAD1
scRNA  <-  merge(scRNAD1, y = c(scRNAD2,scRNAD4,scRNAD6), add.cell.ids = samples[1:4], project = "SC")
dim(scRNA)
save(scRNA,file='XXXX.Rdata')
load('XXXX.Rdata')
library(stringr)
library(Seurat)
table(scRNA@meta.data$orig.ident)
head(scRNA@meta.data)
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^mt-")
HB.genes <- c("Hba1","Hba2","Hbb","Hbd","Hbe1","Hbg1","Hbg2","Hbm","Hbq1","Hbz")
HB_m <- match(HB.genes, rownames(scRNA@assays$RNA)) 
HB.genes <- rownames(scRNA@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
scRNA[["percent.HB"]]<-PercentageFeatureSet(scRNA, features=HB.genes) 
col.num <- length(levels(as.factor(scRNA@meta.data$orig.ident)))
dim(scRNA)
library(ggplot2)
violin <-VlnPlot(scRNA, group.by = "orig.ident",  
                 features = c("nFeature_RNA", "nCount_RNA","percent.mt","percent.HB"), 
                 cols =rainbow(col.num), 
                 pt.size = 0.01, 
                 ncol = 4) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
violin
ggsave("vlnplot_before_qc.pdf", plot = violin, width = 12, height = 6) 
ggsave("vlnplot_before_qc.png", plot = violin, width = 12, height = 6)  
plot1 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.HB")
pearplot <- CombinePlots(plots = list(plot1, plot2, plot3), nrow=1, legend="none") 
pearplot
ggsave("pearplot_before_qc.pdf", plot = pearplot, width = 12, height = 5) 
ggsave("pearplot_before_qc.png", plot = pearplot, width = 12, height = 5)
sc=FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
ggsave("XXXX.png", plot = sc, width = 12, height = 5)
minGene=200
maxGene=6000
pctMT=10
scRNA <- subset(scRNA, subset = nFeature_RNA > minGene & nFeature_RNA < maxGene & percent.mt < pctMT)
table(scRNA@meta.data$orig.ident)
col.num <- length(levels(as.factor(scRNA@meta.data$orig.ident)))
violin <-VlnPlot(scRNA, group.by = "orig.ident",
                 features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), 
                 cols =rainbow(col.num), 
                 pt.size = 0.1, 
                 ncol = 4) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
violin
ggsave("vlnplot_after_qc.pdf", plot = violin, width = 12, height = 6) 
ggsave("vlnplot_after_qc.png", plot = violin, width = 12, height = 6)
dim(scRNA)
save(scRNA,file='XXXX.Rdata')
load('XXXX.Rdata')
table(scRNA@meta.data$orig.ident)
scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)
scRNA<- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)
scale.genes <-  VariableFeatures(scRNA)
scRNA <- ScaleData(scRNA, features = scale.genes)
scRNA<- RunPCA(scRNA, features = VariableFeatures(scRNA))
plot <- DimPlot(scRNA, reduction = "pca", group.by = "orig.ident")
plot
ggsave("XXXX.png", plot = plot, width = 12, height = 6)
ggsave("XXXX.pdf", plot = plot, width = 12, height = 6)
library(harmony)
scRNA_harmony <- RunHarmony(scRNA, group.by.vars = "orig.ident")
plot <- DimPlot(scRNA_harmony, reduction = "harmony", group.by = "orig.ident")
plot
ggsave("XXXX.png", plot = plot, width = 12, height = 6)
ggsave("XXX.pdf", plot = plot, width = 12, height = 6)
saveRDS(scRNA_harmony,file='XXXX.rds')
scRNA_harmony=readRDS('XXXX.rds')
plot <- ElbowPlot(scRNA)
plot
ggsave("XXXX.png", plot = plot, width = 12, height = 6)
ggsave("XXXX.pdf", plot = plot, width = 12, height = 6)
scRNA <- FindNeighbors(scRNA, dims = 1:18)
scRNA <- FindClusters(scRNA)
scRNA <- RunUMAP(scRNA, dims = 1:18)
plot <- DimPlot(scRNA,split.by = 'orig.ident')
plot
ggsave("XXXX.png", plot = plot, width = 24, height = 6)
ggsave("XXXX.pdf", plot = plot, width = 24, height = 6)

scRNA <- RunPCA(scRNA, features = VariableFeatures(scRNA))
plot1 <- DimPlot(scRNA, group.by="orig.ident")
plot2 <- ElbowPlot(scRNA, ndims=30) 
plotc <- plot1+plot2
plotc
ggsave("XXXX.png", plot = plotc, width = 8, height = 4)
ggsave("XXXX.pdf", plot = plotc, width = 8, height = 4)
pc.num=1:24
scRNA <- FindNeighbors(scRNA, dims = 1:18)
scRNA <- FindClusters(scRNA,resolution = 0.8)
table(scRNA@meta.data$seurat_clusters)
metadata <- scRNA@meta.data
cell_cluster <- data.frame(cell_ID=rownames(metadata), cluster_ID=metadata$seurat_clusters)
write.csv(cell_cluster,'XXXX.csv',row.names = F)
scRNA = RunTSNE(scRNA, dims = pc.num)

embed_tsne <- Embeddings(scRNA, 'tsne')   
write.csv(embed_tsne,'XXXX.csv')
plot1 = DimPlot(scRNA, reduction = "tsne", label=T) 
plot1
library(ggplot2)
ggsave("XXXX.png", plot = plot1, width = 9, height = 7)
ggsave("XXXX.pdf", plot = plot1, width = 9, height = 7)
plot2 = DimPlot(scRNA, reduction = "tsne", group.by='orig.ident') 
ggsave("XXXX.png", plot = plot2, width = 9, height = 7)
ggsave("XXXX.pdf", plot = plot2, width = 9, height = 7)
plotc <- plot1+plot2
plotc
ggsave("XXXX.png", plot = plotc, width = 15, height = 6)
ggsave("XXXX.pdf", plot = plotc, width = 15, height = 6)
scRNA <- RunUMAP(scRNA, dims = pc.num)
scRNA <- RunUMAP(scRNA, dims = pc.num)
embed_umap <- Embeddings(scRNA, 'umap') 
write.csv(embed_umap,'XXXX.csv') 
plot3 = DimPlot(scRNA, reduction = "umap", label=T) 
ggsave("XXXX.png", plot = plot3, width = 9, height = 7)
ggsave("XXXX.pdf", plot = plot3, width = 9, height = 7)
plot4 = DimPlot(scRNA, reduction = "umap", group.by='orig.ident')
ggsave("XXXX.png", plot = plot4, width = 9, height = 7)
ggsave("XXXX.pdf", plot = plot4, width = 9, height = 7)
plotc <- plot3+plot4
plotc
ggsave("XXXX.png", plot = plotc, width = 12, height = 6)
ggsave("XXXX.pdf", plot = plotc, width = 12, height = 6)
library(patchwork)
plotc <- plot2+plot4+ plot_layout(guides = 'collect')
plotc
ggsave("XXXX.png", plot = plotc, width = 12, height = 6)
ggsave("XXXX.pdf", plot = plotc, width = 12, height = 6)
library(limma)
saveRDS(scRNA,file='XXXX.rds')
gc()
scRNA <- readRDS('XXXX.rds')
cluster5.markers <- FindMarkers(scRNA, ident.1 = c(11, 32), ident.2 = c(0,1,2,3,4,5,6,7,8,9,10,12,13,14,15,16,17,18,19,
                                                                        20,21,22,23,24,25,26,27,28,29,30,31), 
                                min.pct = 0.25)

write.csv(cluster5.markers,file ='XXXX.csv',quote = F)
Idents(scRNA) <- scRNA@meta.data$seurat_clusters
pbmc.markers <- FindAllMarkers(scRNA, only.pos = TRUE, min.pct = 0.25, 
                               logfc.threshold = 0.25)
write.csv(pbmc.markers,file ='XXXX.csv',quote = F)
setwd("XXXX")
scRNA <- readRDS("XXXX.RDS")
library(Seurat)
library(ggplot2)
plot3 = DimPlot(scRNA, reduction = "umap", label=T) 
plot3
ggsave("XXXX.pdf", plot = plot3, width = 9, height = 7)
plot4 = DimPlot(scRNA, reduction = "tsne", label=T)
plot4
ggsave("XXXX.pdf", plot = plot4, width = 9, height = 7)
plotc <- plot3+plot4
plotc
ggsave("XXXX.pdf", plot = plotc, width = 12, height = 7)


ZONG <- c('Slc4a1','Gypa','Bpgm',#红细胞
          'Rag1','Cd8a','Il7r',#淋巴
          'Aqp3','Sox2','S100a14',#增殖
          'Kcna1','Plp1','Gfra3',	#神经							
          'Cd79a','Ms4a1','Cd79b',#B
          'Ly6g','Cxcr2','Cd177',#Neutrophil
          'C1qa','Cd68','Cd14',#Myeloid
          'Pdgfrb','Kcnj8','Rgs5',	#Pericyte
          'Lrrn4','Bnc1','Msln',#MS	
          'Lmod1','Myh11','Acan',#SMC
          'Pecam1','Cdh5','Vwf',#FB
          'Col1a1','Pdgfra','Dcn',#Endo
          'Mybpc3','Ryr2','Ttn',#CM
          'Sln',#A心房
          'Myl2')#V心室


ZONG <- c('Aqp3','Sox2','S100a14',#增殖
          'Kcna1','Plp1','Gfra3',	#神经							
          'Cd79a','Ms4a1','Cd79b',#B
          'Ly6g','Cxcr2','Cd177',#Neutrophil
          'Lmod1','Myh11','Acan',#SMC
          'Pdgfrb','Kcnj8','Rgs5',#Pericyte
          'Rag1','Cd8a','Il7r',#淋巴
          'C1qa','Cd68','Cd14',#Myeloid
          'Lrrn4','Bnc1','Msln',#MS	
          'Mybpc3','Ryr2','Ttn',#CM
          'Col1a1','Pdgfra','Dcn',#Endo
          'Pecam1','Cdh5','Vwf')#FB

library(ggplot2)
library(viridis)
table(scRNA$orig.ident)
plot <- DotPlot(scRNA, features = ZONG,assay='RNA') + coord_flip()
plot <- DotPlot(scRNA, features = T_marker,assay='RNA') + coord_flip()
plot <- DotPlot(scRNA, features = ZONG,assay='RNA',group.by = 'orig.ident')+ coord_flip()
plot
ggsave("XXXX.png", plot = plot, width = 18, height = 12)
ggsave("XXXX.pdf", plot = plot, width = 18, height = 12)
VlnPlot(scRNA, features = c('Id1','Id3','S100a10','Anxa2','Ahnak'),
        group.by = 'orig.ident',slot = 'data', pt.size=0)
scRNA$celltype <- factor(x=scRNA$celltype,
                              levels = c('Proliferating cells','Neuronal','B cells','Neutrophil',
                                         'Smooth muscle cells','Pericyte','Lymphoid','Myeloid',
                                         'Mesothelial','Cardiomyocyte','Fibroblast','Endothelial cells'))
table(scRNA$celltype)
plot <- DoHeatmap(scRNA,
          features = ZONG,
          group.by = "orig.ident",
          slot = 'data')+
  scale_fill_gradientn(colors = c("#282A62","#692F7C","#B43970","#D96558","#EFA143","#F6C63C"))
plot
DotPlot(scRNA, features = ZONG,assay='RNA')
T_marker <- c('Irx4','Nr2f1')
T_marker <- c('Irx4','Nr2f1','Nppa','Myl7','Myl4','Myh7','Myl2','Fhl2')
FeaturePlot(scRNA,features = 'Myl2',order = T,reduction = 'umap')

library(readxl)

annotation_curated_imm <- read_excel("XXXX")

Idents(scRNA) <- scRNA@meta.data$seurat_clusters
new_ids_imm <- annotation_curated_imm$cell_type_imm
Idents(scRNA) <- scRNA@meta.data$seurat_clusters
names(new_ids_imm) <- levels(scRNA)
scRNA <- RenameIdents(scRNA, new_ids_imm)
scRNA@meta.data$celltype <- Idents(scRNA)

library(patchwork)
p1 = DimPlot(scRNA, group.by="celltype", label=T, label.size=4, reduction='tsne')
p2 = DimPlot(scRNA, group.by="celltype", label=T, label.size=4, reduction='umap')
p4 = DimPlot(scRNA, reduction = "umap", group.by='orig.ident')
p5 <-  p2+p4
p5
p2
p3 = plotc <- p1+p2+ plot_layout(guides = 'collect')  
p3  
scRNA@meta.data$patient=stringr::str_replace(scRNA@meta.data$orig.ident,"[1-9]",'')
scRNA_CM <- subset(scRNA,celltype=="Cardiomyocyte")
scRNA_FB <- subset(scRNA,celltype=="Fibroblast")
scRNA_Endo <- subset(scRNA,celltype=="Endothelial cells")

plot <- DotPlot(scRNA_FB, features = ZONG,assay='RNA',
                group.by = 'orig.ident')+ coord_flip()
plot

VlnPlot(scRNA, features = c('Id1','Id3','S100a10','Anxa2','Ahnak'),
        group.by = 'orig.ident',slot = 'data', pt.size=0)
scRNA$celltype <- factor(x=scRNA$celltype,
                         levels = c('Proliferating cells','Neuronal','B cells','Neutrophil',
                                    'Smooth muscle cells','Pericyte','Lymphoid','Myeloid',
                                    'Mesothelial','Cardiomyocyte','Fibroblast','Endothelial cells'))
table(scRNA$celltype)
plot <- DoHeatmap(scRNA_Endo,
                  features = ZONG,
                  group.by = "orig.ident",
                  slot = 'data')+
  scale_fill_gradientn(colors = c("pink","white","blue"))
plot
library(readxl)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(viridis)
library(reshape2)
library(ggplot2)
prop_df <- table(scRNA@meta.data$celltype,scRNA@meta.data$orig.ident) %>% melt()
colnames(prop_df) <- c("Cluster","Sample","Number")
prop_df$Cluster <- factor(prop_df$Cluster)
library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))) 
sample_color <- col_vector[1:20] 
prop <- ggplot(data = prop_df, aes(x =Number, y = Sample, fill =  Cluster)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values=col_vector[1:20]) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Ratio")+
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)
  ) 
prop
table(scRNA$celltype,scRNA$orig.ident)
FeaturePlot(scRNA,features = 'Myl2',order = T,reduction = 'umap',split.by = 'orig.ident',slot = 'data')
FeaturePlot(scRNA,features = 'Sln',order = T,reduction = 'umap',split.by = 'orig.ident',slot = 'data')


scRNA_CM=subset(scRNA,celltype=='Cardiomyocyte')



cluster5.markers <- FindMarkers(scRNA_CM, ident.1 = 'D1', ident.2 = 'D6', 
                                min.pct = 0.25)
Idents(scRNA_CM) <- scRNA_CM$orig.ident
pbmc.markers <- FindAllMarkers(scRNA_CM, only.pos = TRUE,min.pct = 0.25, 
                               logfc.threshold = 0.25)
FeaturePlot(scRNA_CM,features = 'Gpx3',order = T,reduction = 'tsne',split.by = 'orig.ident',slot = 'data')
VlnPlot(scRNA,features = "ADH1A",pt.size = 1,ncol = 2,
        group.by = 'orig.ident') 
scRNA1=subset(scRNA,celltype=='vCM')
scRNA2=subset(scRNA,celltype=='aCM')
saveRDS(scRNA1,file='vCM.rds')
saveRDS(scRNA2,file='aCM.rds')
setwd("XXXX")
scRNA <- readRDS("XXXX.rds")
table(scRNA@meta.data$seurat_clusters)
scRNA_CM=scRNA

table(scRNA_CM@meta.data$celltype)
library(Seurat)
library(readxl)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(viridis)
scRNAsub<- FindVariableFeatures(scRNA_CM, selection.method = "vst", nfeatures = 2000)
scale.genes <-  VariableFeatures(scRNAsub)
scRNAsub <- ScaleData(scRNAsub, features = scale.genes)
scRNAsub<- RunPCA(scRNAsub, features = VariableFeatures(scRNAsub))
plot <- DimPlot(scRNAsub, reduction = "pca", group.by = "orig.ident")
plot
plot <- ElbowPlot(scRNAsub,ndim=30)
plot
set.seed(1000)
scRNAsub <- RunHarmony(scRNAsub, group.by.vars = "orig.ident")
plot <- DimPlot(scRNAsub, reduction = "harmony", group.by = "orig.ident")
plot
plot <- ElbowPlot(scRNAsub,reduction = 'harmony',ndims = 30)
plot
scRNAsub <- FindNeighbors(scRNAsub,dims = 1:14)
scRNAsub <- FindClusters(scRNAsub,resolution = 0.05)
scRNAsub <- RunUMAP(scRNAsub,dims = 1:14)
table(scRNA@meta.data$orig.ident)
plot <- DimPlot(scRNAsub, label=TRUE,split.by = 'orig.ident') 
plot
plot <- DimPlot(scRNAsub, label=TRUE) 
plot
scRNAsub <- RunTSNE(scRNAsub,dims = 1:14)
table(scRNA@meta.data$orig.ident)
plot <- DimPlot(scRNAsub, label=TRUE,split.by = 'orig.ident',reduction = 'tsne') 
plot
plot <- DimPlot(scRNAsub, label=TRUE,reduction = 'tsne') 
plot
dev.off()
table(scRNAsub@meta.data$seurat_clusters)
library(reshape2)
library(ggplot2)
prop_df <- table(scRNAsub@meta.data$seurat_clusters,scRNAsub@meta.data$orig.ident) %>% melt()
colnames(prop_df) <- c("Cluster","Sample","Number")
prop_df$Cluster <- factor(prop_df$Cluster)
library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))) 
sample_color <- col_vector[1:32] 
prop <- ggplot(data = prop_df, aes(x =Number, y = Sample, fill =  Cluster)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values=col_vector[1:32]) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Ratio")+
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)
  ) 
prop
pbmc.markers <- FindAllMarkers(scRNAsub, only.pos = TRUE, min.pct = 0.25, 
                               logfc.threshold = 0.25)
FeaturePlot(scRNAsub,features = 'Bmp10',order = T,reduction = 'umap')
FeaturePlot(scRNAsub,features = 'Adm',order = T,reduction = 'umap')
FeaturePlot(scRNAsub,features = 'Id2',order = T,reduction = 'umap')
FeaturePlot(scRNAsub,features = 'Fabp5',order = T,reduction = 'umap')
FeaturePlot(scRNAsub,features = 'Col4a1',order = T,reduction = 'umap')
FeaturePlot(scRNAsub,features = 'Kcne1',order = T,reduction = 'tsne')
FeaturePlot(scRNAsub,features = 'Fam220a',order = T,reduction = 'tsne')
FeaturePlot(scRNAsub,features = 'Lrrfip1',order = T,reduction = 'tsne')
FeaturePlot(scRNAsub,features = 'Bmp2',order = T,reduction = 'tsne')
FeaturePlot(scRNAsub,features = 'Ttn',order = T,reduction = 'tsne')
FeaturePlot(scRNAsub,features = 'Pln',order = T,reduction = 'tsne')
FeaturePlot(scRNAsub,features = 'Myom2',order = T,reduction = 'tsne')
FeaturePlot(scRNAsub,features = 'Atp2a2',order = T,reduction = 'tsne')
FeaturePlot(scRNAsub,features = 'Cox7b',order = T,reduction = 'tsne')
FeaturePlot(scRNAsub,features = 'Sfrp1',order = T,reduction = 'tsne')
FeaturePlot(scRNAsub,features = 'Tcf21',order = T,reduction = 'tsne')
FeaturePlot(scRNAsub,features = 'Nrp1',order = T,reduction = 'tsne')
FeaturePlot(scRNAsub,features = 'Bmp4',order = T,reduction = 'tsne')
FeaturePlot(scRNAsub,features = 'Pdgfra',order = T,reduction = 'tsne')
T_marker <- c('Irx4','Nr2f1','Bmp10','Ube2c','S100a10','Lars2','Sfrp1','Abcc9')
T_marker <- c('Ttn','Irx4','Hey2','Slc1a3','Prrx1','Irx3','Mk167','Irx1',
              'Irx2','Cnn1','Crabp2','Bmp2','Tbx3','Rspo3','Nr2f1','Trpm3',
              'Angpt1','Shox2','Isl1','Hcn4')
T_marker <- c('Irx4','Nr2f1','Nppa','Myl7','Myl4','Myh7','Myl2','Fhl2')
FeaturePlot(scRNAsub,features = 'Nr2f1',order = T,reduction = 'tsne')
plot <- DotPlot(scRNAsub, 
                features = T_marker,
                group.by = "seurat_clusters") + coord_flip()
plot

annotation_curated_imm <- read_excel("XXXX.xlsx")
T_celltype <- annotation_curated_imm$cell_type_imm
T_celltype=c('vCM-1','vCM-2','vCM-3')

Idents(scRNAsub) <- scRNAsub@meta.data$seurat_clusters
names(T_celltype) <- levels(scRNAsub)
scRNAsub<- RenameIdents(scRNAsub, T_celltype)
table(scRNA@meta.data$celltype)
scRNAsub@meta.data$T_celltype <- Idents(scRNAsub)

FeaturePlot(scRNAsub,features = 'Abcc9',order = T,reduction = 'tsne',split.by = 'orig.ident')
VlnPlot(scRNA, 
        features = c("Abcc9"),
        pt.size = 1,
        ncol = 2,group.by = 'orig.ident') 
Idents(scRNAsub)=scRNAsub@meta.data$T_celltype
table(scRNA$orig.ident)
library(reshape2)
library(ggplot2)
prop_df <- table(scRNAsub@meta.data$T_celltype,scRNAsub@meta.data$tissue.type) %>% melt()
colnames(prop_df) <- c("Cluster","Sample","Number")
prop_df$Cluster <- factor(prop_df$Cluster)
library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))) 

sample_color <- col_vector[1:10] 

prop <- ggplot(data = prop_df, aes(x =Number, y = Sample, fill =  Cluster)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values=col_vector[1:20]) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Ratio")+
 
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)
  ) 
prop
Idents(scRNA)=scRNA$celltype
library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))) 
plot <- DimPlot(scRNA,split.by = 'orig.ident',cols = col_vector[19:33])
plot
library(Seurat)
library(reshape2)
library(ggplot2)
library(dplyr)
pB2_df <- table(scRNA@meta.data$celltype,scRNA@meta.data$orig.ident) %>% melt()
colnames(pB2_df) <- c("Cluster","Sample","Number")
pB2_df$Cluster <- factor(pB2_df$Cluster)

library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))) 
sample_color <- col_vector[1:15] 
pB4 <- ggplot(data = pB2_df, aes(x =Number, y = Sample, fill =  Cluster)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values=col_vector[1:20]) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Ratio")+
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)
  ) 
pB4

#########
library(Seurat)
scRNA <- readRDS('all.RDS')
library(CellChat)
gc()
table(scRNA@meta.data$orig.ident)
scRNA@meta.data$tissue_type <- scRNA@meta.data$orig.ident
table(scRNA@meta.data$seurat_clusters)
scRNA_chat <-subset(scRNA, tissue_type =='D1')
scRNA_chat <- scRNAsub
table(scRNA_chat@meta.data$seurat_clusters)
meta =scRNA_chat@meta.data # a dataframe with rownames containing cell mata data

data_input <- as.matrix(scRNA_chat@assays$RNA@data)
identical(colnames(data_input),rownames(meta))

library(CellChat)
cellchat <- createCellChat(object = data_input, meta = meta, group.by = "orig.ident")

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
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)

df.net<- subsetCommunication(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
dev.off()
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= T, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()

p_bubble= netVisual_bubble(cellchat,
                    
                           sources.use = c('Proliferating cell'),
                           remove.isolate = FALSE)+coord_flip()
p_bubble
ggsave("D2_Proliferating.png", plot = p_bubble, width = 9, height = 7)
p_bubble= netVisual_bubble(cellchat,
                           
                           sources.use = c('T'),
                           remove.isolate = FALSE)+coord_flip()
p_bubble
ggsave("D4_T.png", plot = p_bubble, width = 9, height = 7)









library(Seurat)
library(BiocManager)
library(stringr)
library(monocle)
library(tidyverse)
library(patchwork)

##monocle
library(monocle)
table(scRNA@meta.data$celltype)
Cells.sub <- subset(scRNA@meta.data, celltype=='Ventricular cell')
Cells.sub <- scRNAsub

data <- as(as.matrix(scRNAsub@assays$RNA@counts), 'sparseMatrix')

pd <- new('AnnotatedDataFrame', data = scRNAsub@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

mycds <- newCellDataSet(data,
                        phenoData = pd,
                        featureData = fd,
                        expressionFamily = negbinomial.size())


mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds, cores=4, relative_expr = TRUE)



disp_table <- dispersionTable(mycds)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
mycds <- setOrderingFilter(mycds, disp.genes)
p3 <- plot_ordering_genes(mycds)
p3
dev.off()
mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree')
mycds <- orderCells(mycds)
plot1 <- plot_cell_trajectory(mycds, color_by = "State")
plot1
table(scRNA@meta.data$seurat_clusters)
plot2 <- plot_cell_trajectory(mycds, color_by = "seurat_clusters")
plot2
plot3 <- plot_cell_trajectory(mycds, color_by = "Pseudotime")
plot3
plot4 <-  plot_cell_trajectory(mycds, color_by = "orig.ident")
plotc <- plot3|plot1|plot4|plot3
plotc
state=pData(mycds)
identical(rownames(scRNAsub@meta.data),rownames(state))
scRNAsub@meta.data$monocle_state=state$State
Idents(scRNAsub)=Idents(scRNAsub)=scRNAsub@meta.data$monocle_state