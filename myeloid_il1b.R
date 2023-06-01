library(Seurat)
library(dplyr)
library(ggplot2)
library(scCustomize)
library(ggrepel)
options(future.globals.maxSize = 40000 * 1024^2)
set.seed(1234)
gc()

count1.data <- Read10X_h5('~/Kropski_Lab Dropbox/Jon Kropski/Lab Data/Agrawal/8697_cb/8697_1_output_filtered.h5')
count1 <- CreateSeuratObject(count1.data$`Gene Expression`, project='Control')
count1_hto <- CreateAssayObject(counts = count1.data$`Antibody Capture`)
count1[['HTO']] <- count1_hto
count1 <- NormalizeData(count1, assay = "HTO", normalization.method = "CLR")
count1 <- HTODemux(count1, assay = "HTO", positive.quantile = 0.99)

#demultiplex and assign ID
Idents(count1) <- 'HTO_maxID'
Idents(count1, cells = WhichCells(count1, idents = c("TotalSeqAHashtag1"))) <- "Control_1"
Idents(count1, cells = WhichCells(count1, idents = c("TotalSeqAHashtag2"))) <- "Control_2"
Idents(count1, cells = WhichCells(count1, idents = c("TotalSeqAHashtag3"))) <- "Control_3"
Idents(count1, cells = WhichCells(count1, idents = c("TotalSeqAHashtag4"))) <- "Control_4"
count1$sample <- Idents(count1)


#read in L-Name/HFD
count2.data <- Read10X_h5('~/Kropski_Lab Dropbox/Jon Kropski/Lab Data/Agrawal/8697_cb/8697_2_output_filtered.h5')
count2 <- CreateSeuratObject(count2.data$`Gene Expression`, project='L-NAME/HFD')
count2_hto <- CreateAssayObject(counts = count2.data$`Antibody Capture`)
count2[['HTO']] <- count2_hto
count2 <- NormalizeData(count2, assay = "HTO", normalization.method = "CLR")
count2 <- HTODemux(count2, assay = "HTO", positive.quantile = 0.99)

#demultiplex and assign ID
Idents(count2) <- 'HTO_maxID'
Idents(count2, cells = WhichCells(count2, idents = c("TotalSeqAHashtag1"))) <- "L-NAME/HFD_1"
Idents(count2, cells = WhichCells(count2, idents = c("TotalSeqAHashtag2"))) <- "L-NAME/HFD_2"
Idents(count2, cells = WhichCells(count2, idents = c("TotalSeqAHashtag3"))) <- "L-NAME/HFD_3"
Idents(count2, cells = WhichCells(count2, idents = c("TotalSeqAHashtag4"))) <- "L-NAME/HFD_4"
count2$sample <- Idents(count2)

adata <- merge(count1, y = c(count2))
Idents(adata) <- 'HTO_classification.global'

#remove doublets and low quality
adata <- subset(adata, idents = c('Singlet', 'Negative'))
adata <- PercentageFeatureSet(adata, pattern = "^mt-", col.name = "percent.mt")
adata <- PercentageFeatureSet(adata, pattern = "^Rp-", col.name = "percent.rb")
adata <- subset(adata, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)

#normalization, scaling and clustering
adata <- SCTransform(adata, variable.features.n = 2000, vars.to.regress = c("percent.mt", 'percent.rb'),return.only.var.genes = F)
adata <- RunPCA(adata, verbose = F)
adata <- RunUMAP(adata, dims = 1:45,)
adata <- FindNeighbors(adata, dims = 1:45)
adata <- FindClusters(adata, resolution = 0.5)


Idents(adata) <- 'seurat_clusters'
DimPlot(adata, label = T, cols = 'polychrome')
DimPlot(adata, group.by = "orig.ident")
DotPlot(adata, features = c("Ptprc", "Epcam", "Pecam1", 'Vwf', 'F8', 'Car4','Ackr1', 'Hey1', 'Plvap', 'Apln', 'Aplnr', 'Ccl21a', "Lyz1", 'S100a8', 'S100a9', 'Itgax', 'Itgam', 'Siglecf', 'Siglec1', 'Ccr5', 'Sell', 'Cd14', 'Cd68', 'Adgre1', 'Mrc1', 'Pparg', 'Cd86', "Cd3e", 'Cd8a', 'Cd4', 'Foxp3', 'Nkg7',  'Ms4a1', 'Iglc1', 'Iglc2', 'Ighd', 'Jchain', 'Irf7', "Col1a1", "Pdgfra", 'Wnt2', 'Wnt5a', 'Pi16', 'Sfrp4', 'Cthrc1', 'Acta2', 'Cspg4', 'Sftpc', 'Foxj1', 'Scgb1a1', 'Hopx', 'Rtkn2', 'Trp63', 'Calca',  'Muc5b', 'Mki67', 'Cdkn1a', 'Cldn4', 'Krt8', 'Krt19', 'Wt1'), cluster.idents = TRUE, cols = c('grey', 'dark blue'))+ theme(axis.text.x = element_text(angle = 45, hjust=1))
DimPlot(adata, label = T, cols = 'polychrome', split.by = 'orig.ident')
VlnPlot(adata, features = 'nFeature_RNA')

#remove heterotypic doublets (cluster 24), low quality (6,16,22)
adata <- subset(adata, idents = c(0:5,7:15, 17:21,23,25:27))
adata <- RunPCA(adata, verbose = F)
adata <- RunUMAP(adata, dims = 1:45,)
adata <- FindNeighbors(adata, dims = 1:45)
adata <- FindClusters(adata, resolution = 0.5)

Idents(adata) <- 'seurat_clusters'
DimPlot(adata, label = T, cols = 'polychrome')
DimPlot(adata, group.by = "orig.ident")
DotPlot(adata, features = c("Ptprc", "Epcam", "Pecam1", 'Vwf', 'F8', 'Car4','Ackr1', 'Hey1', 'Plvap', 'Apln', 'Aplnr', 'Ccl21a', "Lyz1", 'S100a8', 'S100a9', 'Itgax', 'Itgam', 'Siglecf', 'Siglec1', 'Ccr5', 'Sell', 'Cd14', 'Cd68', 'Adgre1', 'Mrc1', 'Pparg', 'Cd86', "Cd3e", 'Cd8a', 'Cd4', 'Foxp3', 'Nkg7',  'Ms4a1', 'Iglc1', 'Iglc2', 'Ighd', 'Jchain', 'Irf7', "Col1a1", "Pdgfra", 'Wnt2', 'Wnt5a', 'Pi16', 'Sfrp4', 'Cthrc1', 'Acta2', 'Cspg4', 'Sftpc', 'Foxj1', 'Scgb1a1', 'Hopx', 'Rtkn2', 'Trp63', 'Calca',  'Muc5b', 'Mki67', 'Cdkn1a', 'Cldn4', 'Krt8', 'Krt19', 'Wt1'), cluster.idents = TRUE, cols = c('grey', 'dark blue'))+ theme(axis.text.x = element_text(angle = 45, hjust=1))

adata <- FindSubCluster(adata, 5, graph.name = 'SCT_snn',  subcluster.name = 'subcluster1', resolution = 0.2)
Idents(adata) <- 'subcluster1'
adata <- FindSubCluster(adata, 15, graph.name = 'SCT_snn',  subcluster.name = 'subcluster2', resolution = 0.2)
Idents(adata) <- 'subcluster2'
adata <- FindSubCluster(adata, 9, graph.name = 'SCT_snn',  subcluster.name = 'subcluster3', resolution = 0.2)
Idents(adata) <- 'subcluster3'
adata <- FindSubCluster(adata, 14, graph.name = 'SCT_snn',  subcluster.name = 'subcluster4', resolution = 0.3)
Idents(adata) <- 'subcluster4'
adata <- FindSubCluster(adata, 19, graph.name = 'SCT_snn',  subcluster.name = 'subcluster5', resolution = 0.4)
Idents(adata) <- 'subcluster5'
adata <- FindSubCluster(adata, 4, graph.name = 'SCT_snn',  subcluster.name = 'subcluster6', resolution = 0.3)
Idents(adata) <- 'subcluster6'
adata <- FindSubCluster(adata, 8, graph.name = 'SCT_snn',  subcluster.name = 'subcluster7', resolution = 0.3)
Idents(adata) <- 'subcluster7'
adata <- FindSubCluster(adata, 11, graph.name = 'SCT_snn',  subcluster.name = 'subcluster8', resolution = 0.2)
Idents(adata) <- 'subcluster8'
adata <- FindSubCluster(adata, 10, graph.name = 'SCT_snn',  subcluster.name = 'subcluster9', resolution = 0.2)
Idents(adata) <- 'subcluster9'
adata <- FindSubCluster(adata, 18, graph.name = 'SCT_snn',  subcluster.name = 'subcluster10', resolution = 0.2)
Idents(adata) <- 'subcluster10'

DimPlot(adata, cols = 'polychrome', label = TRUE, repel = TRUE)
DotPlot(adata, features = c("Ptprc", "Epcam", "Pecam1", 'Vwf', 'F8', 'Car4','Ackr1', 'Hey1', 'Plvap', 'Apln', 'Aplnr', 'Ccl21a', "Lyz1", 'S100a8', 'S100a9', 'Itgax', 'Itgam', 'Siglecf', 'Siglec1', 'Ccr5', 'Sell', 'Cd14', 'Cd68', 'Adgre1', 'Mrc1', 'Pparg', 'Cd86', 'Il1b', "Cd3e", 'Cd8a', 'Cd4', 'Foxp3', 'Nkg7',  'Ms4a1', 'Iglc1', 'Iglc2', 'Ighd', 'Jchain', 'Irf7', "Col1a1", "Pdgfra", 'Wnt2', 'Wnt5a', 'Pi16', 'Sfrp4', 'Cthrc1', 'Acta2', 'Cspg4', 'Sftpc', 'Foxj1', 'Scgb1a1', 'Hopx', 'Rtkn2', 'Trp63', 'Calca',  'Muc5b', 'Mki67', 'Cdkn1a', 'Cldn4', 'Krt8', 'Krt19', 'Wt1'), cluster.idents = TRUE, cols = c('grey', 'dark blue'))+ theme(axis.text.x = element_text(angle = 45, hjust=1))
VlnPlot(adata, features = 'nFeature_RNA')


#remove 14_0, 4_5, 9_3, 4_4, 5_2, 8_3, 10_2, 10_3, 18_0 doublets
adata <- subset(adata, idents = c(0,'11_1','11_0',2,22,'18_1', '9_1', '9_0', '9_2','14_1', '14_2', '15_0','5_0', 13,'5_1', 16,3,'8_1', '8_0', '8_2', '4_3', '4_2', '10_0', 6,'4_1', '4_0','10_1',7,1,'19_0', '19_1','10_4',12,21,17,20,23))
adata <- RunPCA(adata, verbose = F)
adata <- RunUMAP(adata, dims = 1:45)
adata <- FindNeighbors(adata, dims = 1:45)

Idents(adata) <- 'subcluster10'
DimPlot(adata, label = T, cols = 'polychrome')
DimPlot(adata, group.by = "orig.ident")
DotPlot(adata, features = c("Ptprc", "Epcam", "Pecam1", 'Vwf', 'F8', 'Car4','Ackr1', 'Hey1', 'Plvap', 'Apln', 'Aplnr', 'Ccl21a', "Lyz1", 'S100a8', 'S100a9', 'Itgax', 'Itgam', 'Siglecf', 'Siglec1', 'Ccr5', 'Sell', 'Cd14', 'Cd68', 'Adgre1', 'Mrc1', 'Pparg', 'Cd86', "Cd3e", 'Cd8a', 'Cd4', 'Foxp3', 'Nkg7',  'Ms4a1', 'Iglc1', 'Iglc2', 'Ighd', 'Jchain', 'Irf7', "Col1a1", "Pdgfra", 'Wnt2', 'Wnt5a', 'Pi16', 'Sfrp4', 'Cthrc1', 'Acta2', 'Cspg4', 'Sftpc', 'Foxj1', 'Scgb1a1', 'Hopx', 'Rtkn2', 'Trp63', 'Calca',  'Muc5b', 'Mki67', 'Cdkn1a', 'Cldn4', 'Krt8', 'Krt19', 'Wt1'), cluster.idents = TRUE, cols = c('grey', 'dark blue'))+ theme(axis.text.x = element_text(angle = 45, hjust=1))
DotPlot_scCustom(adata, features = c("Ptprc", "Epcam", "Pecam1", 'Vwf', 'F8', 'Car4','Ackr1', 'Hey1', 'Plvap', 'Apln', 'Aplnr', 'Ccl21a', "Lyz1", 'S100a8', 'S100a9', 'Ly6g', 'Elane', 'Ly6c', 'Il1b', 'Itgax', 'Itgam', 'Siglecf', 'Siglec1', 'Ccr5', 'Sell', 'Cd14', 'Cd68', 'Adgre1', 'Mrc1', 'Pparg', 'Cd86', "Cd3e", 'Cd8a', 'Cd4', 'Foxp3', 'Nkg7',  'Ms4a1', 'Iglc1', 'Iglc2', 'Ighd', 'Jchain', 'Irf7', "Col1a1", "Pdgfra", 'Wnt2', 'Wnt5a', 'Pi16', 'Sfrp4', 'Cthrc1', 'Acta2', 'Cspg4', 'Sftpc', 'Foxj1', 'Scgb1a1', 'Hopx', 'Rtkn2', 'Trp63', 'Calca',  'Muc5b', 'Mki67', 'Cdkn1a', 'Cldn4', 'Krt8', 'Krt19', 'Wt1'), cluster.idents = T)+ theme(axis.text.x = element_text(angle = 45, hjust=1))
marker_genes <- c('Rtkn2', 'Ager', 'Sftpd', 'Lamp3', 'Scgb3a2', 'Cyp2f2', 'Foxj1', 'Trp73', 'Col1a1', 'Pdgfra', 'Tcf21', 'Aspn', 'Pi16', 'Acta2', 'Pdgfrb', 'Wt1', 'Msln', 'Pecam1', 'Hey1', 'Aplnr', 'Car4', 'Vwf', 'Plvap', 'Ccl21a', 'Ptprc','Itgam', 'Itgax', 'Cd68', 'Pparg', 'Cd14', 'Cd86', 'S100a8', 'S100a9', 'Il1b', 'Cpa3', 'Cd19', 'Bank1', 'Cd3e', 'Cd4', 'Cd8a', 'Foxp3', 'Nkg7', 'Pf4', 'Mki67')
DotPlot_scCustom(adata, features = marker_genes, cluster.idents = TRUE)+ theme(axis.text.x = element_text(angle = 45, hjust=1))

#assign cell types
Idents(adata) <- 'subcluster10'
Idents(adata, cells = WhichCells(adata, idents = c('10_4'))) <- "Cycling"
Idents(adata, cells = WhichCells(adata, idents = c(22))) <- "MK/Platelet"

Idents(adata, cells = WhichCells(adata, idents = c(6))) <- "NK"
Idents(adata, cells = WhichCells(adata, idents = c('10_0'))) <- "NKT"
Idents(adata, cells = WhichCells(adata, idents = c('4_2', '4_3'))) <- "CD8"
Idents(adata, cells = WhichCells(adata, idents = c('4_0', '4_1', '10_1'))) <- "CD4"
Idents(adata, cells = WhichCells(adata, idents = c(1,7))) <- "B cell"

Idents(adata, cells = WhichCells(adata, idents = c(23))) <- "Mast"
Idents(adata, cells = WhichCells(adata, idents = c('8_0', '8_1'))) <- "PMN"
Idents(adata, cells = WhichCells(adata, idents = c('8_2'))) <- "Activated monocyte"
Idents(adata, cells = WhichCells(adata, idents = c('5_0'))) <- "Monocyte"
Idents(adata, cells = WhichCells(adata, idents = c(13))) <- "cDC"
Idents(adata, cells = WhichCells(adata, idents = c('5_1'))) <- "MDM"
Idents(adata, cells = WhichCells(adata, idents = c(16))) <- "moDC"
Idents(adata, cells = WhichCells(adata, idents = c(3))) <- "AMO"

Idents(adata, cells = WhichCells(adata, idents = c('19_0', '19_1'))) <- "Lymphatic"
Idents(adata, cells = WhichCells(adata, idents = c('11_0'))) <- "Venule"
Idents(adata, cells = WhichCells(adata, idents = c(2))) <- "aCap"
Idents(adata, cells = WhichCells(adata, idents = c(0))) <- "gCap"
Idents(adata, cells = WhichCells(adata, idents = c('11_1'))) <- "Arteriole"

Idents(adata, cells = WhichCells(adata, idents = c(17))) <- "Mesothelial"
Idents(adata, cells = WhichCells(adata, idents = c('14_1'))) <- "Pericyte"
Idents(adata, cells = WhichCells(adata, idents = c('14_2'))) <- "SMC"

Idents(adata, cells = WhichCells(adata, idents = c('15_0'))) <- "MyoFB"
Idents(adata, cells = WhichCells(adata, idents = c('9_0', '9_1'))) <- "Alveolar FB"
Idents(adata, cells = WhichCells(adata, idents = c('9_2'))) <- "Adventitial FB"

Idents(adata, cells = WhichCells(adata, idents = c(20))) <- "Multiciliated"
Idents(adata, cells = WhichCells(adata, idents = c(21))) <- "Secretory"
Idents(adata, cells = WhichCells(adata, idents = c(12))) <- "AT2"
Idents(adata, cells = WhichCells(adata, idents = c('18_1'))) <- "AT1"
adata$celltype <- Idents(adata)

DotPlot_scCustom(adata, features = marker_genes, cluster.idents = FALSE)+ theme(axis.text.x = element_text(angle = 45, hjust=1))

#reverse celltypes
Idents(adata) <- 'subcluster10'
Idents(adata, cells = WhichCells(adata, idents = c('18_1'))) <- "AT1"
Idents(adata, cells = WhichCells(adata, idents = c(12))) <- "AT2"
Idents(adata, cells = WhichCells(adata, idents = c(21))) <- "Secretory"
Idents(adata, cells = WhichCells(adata, idents = c(20))) <- "Multiciliated"
Idents(adata, cells = WhichCells(adata, idents = c('9_2'))) <- "Adventitial FB"
Idents(adata, cells = WhichCells(adata, idents = c('9_0', '9_1'))) <- "Alveolar FB"
Idents(adata, cells = WhichCells(adata, idents = c('15_0'))) <- "MyoFB"
Idents(adata, cells = WhichCells(adata, idents = c('14_2'))) <- "SMC"
Idents(adata, cells = WhichCells(adata, idents = c('14_1'))) <- "Pericyte"
Idents(adata, cells = WhichCells(adata, idents = c(17))) <- "Mesothelial"
Idents(adata, cells = WhichCells(adata, idents = c('11_1'))) <- "Arteriole"
Idents(adata, cells = WhichCells(adata, idents = c(0))) <- "gCap"
Idents(adata, cells = WhichCells(adata, idents = c(2))) <- "aCap"
Idents(adata, cells = WhichCells(adata, idents = c('11_0'))) <- "Venule"
Idents(adata, cells = WhichCells(adata, idents = c('19_0', '19_1'))) <- "Lymphatic"
Idents(adata, cells = WhichCells(adata, idents = c(3))) <- "AMO"
Idents(adata, cells = WhichCells(adata, idents = c(16))) <- "moDC"
Idents(adata, cells = WhichCells(adata, idents = c('5_1'))) <- "MDM"
Idents(adata, cells = WhichCells(adata, idents = c(13))) <- "cDC"
Idents(adata, cells = WhichCells(adata, idents = c('5_0'))) <- "Monocyte"
Idents(adata, cells = WhichCells(adata, idents = c('8_2'))) <- "Activated monocyte"
Idents(adata, cells = WhichCells(adata, idents = c('8_0', '8_1'))) <- "PMN"
Idents(adata, cells = WhichCells(adata, idents = c(23))) <- "Mast"
Idents(adata, cells = WhichCells(adata, idents = c(1,7))) <- "B cell"
Idents(adata, cells = WhichCells(adata, idents = c('4_0', '4_1', '10_1'))) <- "CD4"
Idents(adata, cells = WhichCells(adata, idents = c('4_2', '4_3'))) <- "CD8"
Idents(adata, cells = WhichCells(adata, idents = c('10_0'))) <- "NKT"
Idents(adata, cells = WhichCells(adata, idents = c(6))) <- "NK"
Idents(adata, cells = WhichCells(adata, idents = c(22))) <- "MK/Platelet"
Idents(adata, cells = WhichCells(adata, idents = c('10_4'))) <- "Cycling"
adata$celltype_rev <- Idents(adata)

Idents(adata) <- 'celltype'
glasbey_pal <- DiscretePalette_scCustomize(num_colors = 32, palette = "glasbey")
ditto_pal <- DiscretePalette_scCustomize(num_colors = 30, palette = "ditto_seq")
Idents(adata) <- 'celltype'
DimPlot_scCustom(adata, label = F, figure_plot = T, colors_use = DiscretePalette_scCustomize(num_colors = 30, palette = "ditto_seq"))
DimPlot_scCustom(adata, label = F, figure_plot = T, colors_use = DiscretePalette_scCustomize(num_colors = 30, palette = 'glasbey'))

Idents(adata) <- 'orig.ident'
DimPlot_scCustom(adata, label = F, figure_plot = T, colors_use = DiscretePalette_scCustomize(num_colors = 2, palette = "glasbey"), pt.size = 0.1) + theme(plot.title = element_blank())
Idents(adata) <- 'celltype'

Idents(adata) <- 'celltype_rev'
DotPlot_scCustom(adata, features = marker_genes, cluster.idents = FALSE)+ theme(axis.text.x = element_text(angle = 45, hjust=1))

#cell proportion plot
Idents(adata) <- 'celltype'
ct_prop_table <- as.data.frame(prop.table(table(Idents(adata), adata$sample), margin = 2))
ct_prop_table <- ct_prop_table %>% 
  rename(
    Cell_type = Var1,
    Condition = Var2
  )

#ct_prop_table$Condition <- factor(ct_prop_table$Condition, levels = c('Control', 'L-NAME/HFD'))

p1<- ggplot(ct_prop_table, aes(fill=Cell_type, y=Freq, x=Condition)) + 
  geom_bar(position="stack", stat="identity")+
  ggtitle("Cell Types by Condition") +
  xlab("")

glasbey_pal <- DiscretePalette_scCustomize(num_colors = 32, palette = "glasbey")

p1 + scale_fill_manual(values = glasbey_pal)


#save file
saveRDS(adata, file = '~/Kropski_Lab Dropbox/Jon Kropski/Lab Data/Agrawal/8697_cb/lname_control_20230426.rds')


#calculate all markers
all_markers <- FindAllMarkers(adata)
write.csv(all_markers, file ='~/Kropski_Lab Dropbox/Jon Kropski/Lab Data/Agrawal/8697_cb/all_markers_202304026.csv')

#calculate cell proportions
proportions <- prop.table(table(Idents(adata), adata$orig.ident), margin = 2)
proportions_sample <- prop.table(table(Idents(adata), adata$sample), margin = 2)

cell_counts <- table(Idents(adata), adata$orig.ident)
write.csv(cell_counts, file ='~/Kropski_Lab Dropbox/Jon Kropski/Lab Data/Agrawal/8697_cb/cell_count_20230426.csv')


write.csv(proportions, file ='~/Kropski_Lab Dropbox/Jon Kropski/Lab Data/Agrawal/8697_cb/proportions_20230426.csv')
write.csv(proportions_sample, file ='~/Kropski_Lab Dropbox/Jon Kropski/Lab Data/Agrawal/8697_cb/proportions_sample_20230426.csv')


#cell_type_deg
main_dir <- "~/Kropski_Lab Dropbox/Jon Kropski/Lab data/Agrawal/8697_cb/"
date <- gsub("-", "", Sys.Date())


dir.create(file.path(main_dir, date), showWarnings = FALSE)
setwd(file.path(main_dir, date))
getwd()

#DefaultAssay(adata) <- 'RNA'
#adata <- NormalizeData(adata)
#adata <- ScaleData(adata, vars.to.regress = c('percent.mt', 'percent.rb', 'nCount_RNA'))


DefaultAssay(adata) <- 'SCT'
cell_list <- SplitObject(adata, split.by = 'celltype')

lname_deg <- lapply(cell_list, function(xx){
  print(unique(xx@meta.data$celltype))
  if(length(unique(xx@meta.data$orig.ident)) > 1) {
    FindMarkers(xx, group.by = "orig.ident", ident.1 = "L-NAME/HFD", ident.2 = "Control", test.use = "wilcox")
  } 
  else{
    return(NULL)
  } 
})

for(i in 1:length(cell_list)){
  write.csv(lname_deg[[i]], paste(gsub("/", "", unique(cell_list[[i]]@meta.data$celltype)), "_deg", ".csv"), sep =",", quote = F)
}

gc()




Idents(adata) <- 'celltype_rev'
DotPlot_scCustom(adata, features = "Il1b")

Idents(adata) <- 'celltype_rev'
myeloid <- subset(adata, idents = c('Monocyte', 'AMO', 'moDC', 'MDM', 'cDC', 'Monocyte'))
Idents(myeloid) <- 'celltype_rev'
VlnPlot(myeloid, features = "Il1b", split.by = 'orig.ident')


#Volcano plot
amo_deg <- read.csv('~/Kropski_Lab Dropbox/Jon Kropski/Lab Data/Agrawal/8697_cb/20230329/AMO _deg2.csv')
amo_deg$diffexpressed <- "NO"
amo_deg$diffexpressed[amo_deg$avg_log2FC > 0.25 & amo_deg$p_val_adj < 0.05] <- "UP"
amo_deg$diffexpressed[amo_deg$avg_log2FC < -0.25 & amo_deg$p_val_adj < 0.05] <- "DOWN"

amo_deg$delabel <- NA
amo_deg$delabel[amo_deg$diffexpressed != "NO"] <- amo_deg$X[amo_deg$diffexpressed != "NO"]

p <- ggplot(data=amo_deg, aes(x=avg_log2FC, y=-log10(p_val_adj), col=diffexpressed, label=delabel)) + geom_point()
p2 <- p + geom_vline(xintercept=c(-0.25, 0.25), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")+ theme_minimal() +
  geom_text_repel()
p2

mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")
p3 <- p2 + scale_colour_manual(values = mycolors)
p3

gcap_deg <- read.csv('~/Kropski_Lab Dropbox/Jon Kropski/Lab Data/Agrawal/8697_cb/20230329/gcap _deg .csv')
gcap_deg$diffexpressed <- "NO"
gcap_deg$diffexpressed[gcap_deg$avg_log2FC > 0.25 & gcap_deg$p_val_adj < 0.05] <- "UP"
gcap_deg$diffexpressed[gcap_deg$avg_log2FC < -0.25 & gcap_deg$p_val_adj < 0.05] <- "DOWN"

gcap_deg$delabel <- NA
gcap_deg$delabel[gcap_deg$diffexpressed != "NO"] <- gcap_deg$X[gcap_deg$diffexpressed != "NO"]

p4 <- ggplot(data=gcap_deg, aes(x=avg_log2FC, y=-log10(p_val_adj), col=diffexpressed, label=delabel)) + geom_point()
p5 <- p4 + geom_vline(xintercept=c(-0.25, 0.25), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")+ theme_minimal() +
  geom_text_repel()
p5

mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")
p6 <- p5 + scale_colour_manual(values = mycolors)
p6

acap_deg <- read.csv('~/Kropski_Lab Dropbox/Jon Kropski/Lab Data/Agrawal/8697_cb/20230329/acap _deg .csv')
acap_deg$diffexpressed <- "NO"
acap_deg$diffexpressed[acap_deg$avg_log2FC > 0.25 & acap_deg$p_val_adj < 0.05] <- "UP"
acap_deg$diffexpressed[acap_deg$avg_log2FC < -0.25 & acap_deg$p_val_adj < 0.05] <- "DOWN"

acap_deg$delabel <- NA
acap_deg$delabel[acap_deg$diffexpressed != "NO"] <- acap_deg$X[acap_deg$diffexpressed != "NO"]

p7 <- ggplot(data=acap_deg, aes(x=avg_log2FC, y=-log10(p_val_adj), col=diffexpressed, label=delabel)) + geom_point()
p8 <- p7 + geom_vline(xintercept=c(-0.25, 0.25), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")+ theme_minimal() +
  geom_text_repel()
p8

mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")
p9 <- p8 + scale_colour_manual(values = mycolors)
p9

#cell_number_deg
cell_deg <- read.csv('~/Kropski_Lab Dropbox/Jon Kropski/Lab Data/Agrawal/8697_cb/deg_number.csv')

cell_deg$Cell_type <- factor(cell_deg$Cell_type, levels = c(
"AT1",
"AT2",
"Secretory",
"Multiciliated",
"Adventitial FB",
"Alveolar FB",
"MyoFB",
"SMC",
"Pericyte",
"Mesothelial",
"Arteriole",
"gCap",
"aCap",
"Venule",
"Lymphatic",
"AMO",
"moDC",
"MDM",
"cDC",
"PMN",
"B cell",
"CD4",
"CD8",
"NK",
"MK/Platelet",
"Cycling")
)

p11<- ggplot(cell_deg, aes(y=Total_deg, x=Cell_type, fill=Cell_type)) + 
  geom_bar(stat="identity")+
  ggtitle("Cell Types by Condition") +
  xlab("")

glasbey_pal <- DiscretePalette_scCustomize(num_colors = 32, palette = "ditto_seq")

p11 + scale_fill_manual(values = glasbey_pal)+ coord_flip()

p12<- ggplot(cell_deg, aes(y=Total_deg, x=factor(Cell_type,levels = rev(levels(factor(Cell_type)))), fill=Cell_type)) + 
  geom_bar(stat="identity")+
  ggtitle("Differentially Expressed Genes") +
  xlab("")
p12 + scale_fill_manual(values = glasbey_pal)+ coord_flip()

#misc
DotPlot(adata, features = 'Il1b', split.by = 'orig.ident')

Idents(adata) <- 'celltype_rev'
myeloid <- subset(adata, idents = c('Activated monocyte', 'Monocyte', 'AMO', 'moDC', 'MDM', 'cDC', 'PMN'))
Idents(myeloid) <- 'celltype_rev'
VlnPlot(myeloid, features = "Il1b", split.by = 'orig.ident')
DotPlot_scCustom(myeloid, features = 'Il1b', cluster.idents = FALSE)+ theme(axis.text.x = element_text(angle = 45, hjust=1))

Idents(adata) <- 'celltype_rev'
il1b <- subset(adata, idents = c('Activated monocyte','PMN'))
VlnPlot_scCustom(il1b, features = "Il1b", split.by = 'orig.ident')

Idents(il1b) <- 'orig.ident'
FindMarkers(il1b, features = 'Il1b', ident.1 = 'Control', logfc.threshold = 0.01 )

mono <- subset(adata, idents = c('Activated monocyte'))
Idents(mono) <- 'orig.ident'
VlnPlot_scCustom(mono, features = "Il1b", colors_use = DiscretePalette_scCustomize(num_colors = 2, palette = "glasbey" ))

il1b <- subset(adata, idents = c('Activated monocyte','PMN'))
Idents(il1b) <- 'orig.ident'
FindMarkers(il1b, features = 'Il1b', ident.1 = 'Control', logfc.threshold = 0.01 )


DotPlot_scCustom(adata, features = c('Ccl2', 'Cxcl1', 'Cxcl2'), cluster.idents = FALSE)+ theme(axis.text.x = element_text(angle = 45, hjust=1))
VlnPlot_scCustom(adata, features = 'Ccl2', split.by = 'orig.ident')

amo_deg = c("Rbfox1",
            "Fabp1",
            "Spag11b",
            "Cd36",
            "Acaa1b",
            "Cdc42ep3",
            "Acot1",
            "Syne2",
            "Rgl1",
            "Ccnd2",
            "Fcgr4",
            "Adgrl2",
            "Abcd2",
            "Card11",
            "Dipk1a",
            "Snx24",
            "Nabp1",
            "Pvt1",
            "Itpripl2",
            "Runx2",
            "Tnfrsf21",
            "mt-Nd3",
            "Plcb1",
            "Gns",
            "Hip1",
            "Peli1",
            "Bzw2",
            "Lima1",
            "Ubash3b",
            "mt-Nd1",
            "Cotl1",
            "Cd164",
            "Hsd17b4",
            "Gmds",
            "AI504432",
            "C5ar1",
            "Dhrs3",
            "BC005537",
            "Atp10a",
            "Nceh1",
            "Hnrnpf",
            "Frmd4a",
            "Fam129a",
            "Gnas",
            "Lpl",
            "Shtn1",
            "Atp6v0d2",
            "Fam49b",
            "mt-Nd2",
            "mt-Atp6",
            "S100a10",
            "Prdx5",
            "Cyba",
            "Cox5b",
            "Cd302",
            "Ubb",
            "Atpif1",
            "Mrps21",
            "Cd52",
            "Dennd4a",
            "Slfn2",
            "Syk",
            "Lyz2",
            "Mgst3",
            "Cd63",
            "Srgn",
            "Tspo",
            "Rsrp1",
            "St7",
            "Napsa",
            "Pde4b",
            "Clec4d",
            "Ccrl2",
            "Mrpl52",
            "Tnfaip3",
            "Osm",
            "Ccnd3",
            "Gm47283",
            "Wfdc21",
            "S100a6",
            "H2-Ab1",
            "Crip1",
            "H2-Aa",
            "Il1b",
            "Fkbp5",
            "Fn1",
            "Cd74",
            "Ier3",
            "Wfdc17",
            "Nfkbia",
            "Cxcl2",
            "Cfb")

amo <- subset(adata, idents = 'AMO')
Idents(amo) <- 'orig.ident'
amo <- ScaleData(amo, features = amo_deg)
DoHeatmap(amo, features = amo_deg, group.colors = glasbey_pal)

#Volcano plot
library(ggrepel)
amo_deg <- read.csv('~/Kropski_Lab Dropbox/Jon Kropski/Lab Data/Agrawal/20230426/amo_deg2.csv')
amo_deg$diffexpressed <- "NO"
amo_deg$diffexpressed[amo_deg$avg_log2FC > 0.25 & amo_deg$p_val_adj < 0.05] <- "UP"
amo_deg$diffexpressed[amo_deg$avg_log2FC < -0.25 & amo_deg$p_val_adj < 0.05] <- "DOWN"

amo_deg$delabel <- NA
amo_deg$delabel[amo_deg$diffexpressed != "NO"] <- amo_deg$X[amo_deg$diffexpressed != "NO"]

p <- ggplot(data=amo_deg, aes(x=avg_log2FC, y=-log10(p_val_adj), col=diffexpressed, label=delabel)) + geom_point()
p2 <- p + geom_vline(xintercept=c(-0.25, 0.25), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")+ theme_minimal() +geom_point(size = 0.5)+ scale_x_continuous("avg_log2FC", limits = c(-1.7, 1.7))
p2

mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")
p3 <- p2 + scale_colour_manual(values = mycolors) + geom_text_repel()
p3

