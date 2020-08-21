#extract t/nk cells
seu[,seu$celltype %in% c("T cell","NK cell")] ->seu_T_NK

#remove doubelt because driving clustering
seu_T_NK <- seu_T_NK[,seu_T_NK$doublet==FALSE]

seu_T_NK$previous_cluster <- seu$seurat_clusters[seu$celltype %in% c("T cell","NK cell")]
#remove unwanted genes
seu_T_NK = seu_T_NK[!(grepl("^MTRNR",rownames(seu_T_NK)) | grepl("MT-",rownames(seu_T_NK)) | grepl("ACTB",rownames(seu_T_NK))| grepl("RPS",rownames(seu_T_NK)) |grepl("RPL",rownames(seu_T_NK)) ),]

library(Seurat)
library(SeuratData)
library(SeuratWrappers)
seu_T_NK <- RunFastMNN(object.list = SplitObject(seu_T_NK, split.by = "Sample"),features = 2000)

seu_T_NK <- RunUMAP(seu_T_NK, reduction = "mnn", dims = 1:30,seed.use=22L)
seu_T_NK <- FindNeighbors(seu_T_NK, reduction = "mnn", dims = 1:30 ,k.param = 10)
seu_T_NK <- FindClusters(seu_T_NK)
#seu_T_NK <- FindClusters(seu_T_NK,resolution = 0.4)


library(SingleR)
pred_cell_fine_tnk <- pred_cell_fine[rownames(pred_cell_fine) %in% colnames(seu_T_NK),]
tab <- table(Assigned = pred_cell_fine_tnk$pruned.labels, Cluster = seu_T_NK$seurat_clusters)

library(pheatmap)
pdf("T_NK/seurat_pheatmap_cluster_labels_fine.pdf",height = 7,width=8)
p <-pheatmap(log2(tab + 10), color = colorRampPalette(c("white", "blue"))(101))
print(p)
dev.off()

#Idents(seu_T_NK) <- seu_T_NK$RNA_snn_res.0.4
Idents(seu_T_NK) <- seu_T_NK$seurat_clusters
FindAllMarkers(seu_T_NK, only.pos = TRUE) -> seu_nk_markers

saveRDS(seu_T_NK,"T_NK/seu_T_NK.rds")
saveRDS(seu_nk_markers,"T_NK/seu_nk_markers.rds")

top10 <- seu_nk_markers %>% group_by(cluster) %>% top_n(n = 30, wt = -p_val_adj) %>% top_n(n = 15, wt = avg_logFC)
top10[!duplicated(top10$gene),]->top10

for (i in c(3,8,16)) {
  seu_nk_markers$gene[seu_nk_markers$cluster==i] -> features
  if (length(features)>180){features<- features[1:180]}
  pdf(paste0("T_NK/individual/dotplot-",i,".pdf"),height =20,width = 7)
  p <- DotPlot(seu_T_NK,features = features,cols="Spectral",cluster.idents = TRUE) + coord_flip()
  print(p)
  dev.off()
}

for (i in unique(seu_T_NK$seurat_clusters)){
  seu_nk_markers$gene[seu_nk_markers$cluster==i][1:12] -> features
  png(paste0("T_NK/individual/umap-",i,".png"),width = 10, height = 8,units="in",res=300)
  p<- FeaturePlot(seu_T_NK, features = features)
  print(p)
  dev.off()
}



Idents(seu_T_NK) <- seu_T_NK$seurat_clusters


### assign cell types
seu_T_NK$celltype <- 1
seu_T_NK$celltype[seu_T_NK$seurat_clusters=="0"] <- "Naive CD4+ T"
seu_T_NK$celltype[seu_T_NK$seurat_clusters=="1"] <- "Effector memory CD4+ T"

seu_T_NK$celltype[seu_T_NK$seurat_clusters=="2"|seu_T_NK$seurat_clusters=="6"] <- "Effector CD8+ T"
seu_T_NK$celltype[seu_T_NK$seurat_clusters=="3"] <- "Naive CD8+ T"
seu_T_NK$celltype[seu_T_NK$seurat_clusters=="4"|seu_T_NK$seurat_clusters=="8"] <- "CD56dim CD16+ NK"
seu_T_NK$celltype[seu_T_NK$seurat_clusters=="5"] <- "Central memory CD8+ T"

seu_T_NK$celltype[seu_T_NK$seurat_clusters=="7"] <- "Effector memory CD4+ T"
seu_T_NK$celltype[seu_T_NK$seurat_clusters=="9"] <- "γδT"

seu_T_NK$celltype[seu_T_NK$seurat_clusters=="10"] <- "Effector CD8+ T"

seu_T_NK$celltype[seu_T_NK$seurat_clusters=="11"] <- "Effector CD4+ T"
seu_T_NK$celltype[seu_T_NK$seurat_clusters=="12"] <- "MAIT"
seu_T_NK$celltype[seu_T_NK$seurat_clusters=="13"] <- "Central memory CD4+ T"
seu_T_NK$celltype[seu_T_NK$seurat_clusters=="14"] <- "CD56bri CD16- NK"

seu_T_NK$celltype[seu_T_NK$seurat_clusters=="15"] <- "Naive CD4+ T"
seu_T_NK$celltype[seu_T_NK$seurat_clusters=="16"] <- "Proliferative T/NK"
seu_T_NK$celltype[seu_T_NK$seurat_clusters=="17"&seu_T_NK$previous_cluster=="22"] <- "iNK"
seu_T_NK$celltype[seu_T_NK$seurat_clusters=="17"&!(seu_T_NK$previous_cluster=="22")] <- "Th2-like"


###
seu_T_NK@meta.data ->metadata

tmp1 = metadata %>% group_by(patient,days) %>% summarise(cell_n=n())
tmp2 = metadata %>% group_by(patient,days,seurat_clusters) %>% summarise(cluster_n=n())
tmp2 = tmp2 %>% left_join(tmp1)
tmp2$cluster_n = tmp2$cluster_n/tmp2$cell_n
tmp2 = tmp2 %>% full_join(unique(metadata[,c("patient","condition")]),by=c("patient"="patient"))
tmp2 = tmp2 %>% full_join(unique(metadata[,c("patient","days","Sample")]))


Idents(seu_T_NK) <- seu_T_NK$celltype
UMAPPlot(seu_T_NK,cols=c(celltype_color,"orange","pink"))
ggsave("T_NK/celltype_umap.png")

seu_T_NK$celltype <- factor(seu_T_NK$celltype,levels=c("Naive CD4+ T","Central memory CD4+ T","Effector memory CD4+ T","Effector CD4+ T",
                                                       "Naive CD8+ T",
                                                       "Central memory CD8+ T","Effector CD8+ T",
                                                       "γδT","Proliferative T/NK","MAIT","CD56bri CD16- NK","CD56dim CD16+ NK",
                                                       "iNK","Th2-like"))
Idents(seu_T_NK) <- seu_T_NK$celltype
UMAPPlot(seu_T_NK,cols=c(celltype_color[-10],"orange","pink")) +theme(axis.ticks = element_blank(),axis.text = element_blank())
ggsave("fig3/new/celltype_umap3.png",width=10,height=6.5,units="in",dpi=300)

#assign stage
seu_T_NK$condition[seu_T_NK$condition=="Health control"]="HC"
seu_T_NK@meta.data ->metadata


metadata$stage = ">20"
metadata$stage[metadata$days<20] = "10-20"
metadata$stage[metadata$days<10] = "<10"
metadata$stage[metadata$condition=="Asymptomatic"]="Asymptomatic"
metadata$stage[metadata$condition=="HC"]="HC"
metadata$stage[metadata$condition=="Severe"]="Severe"
metadata$stage = factor(metadata$stage,levels = c("HC","Asymptomatic","Severe","<10","10-20",">20"))
#marker genes dotplots

#heatmap t_nk 
tmp1 = metadata %>% group_by(patient,days) %>% summarise(cell_n=n())
tmp2 = metadata %>% group_by(patient,days,celltype) %>% summarise(cluster_n=n())
tmp2 = tmp2 %>% left_join(tmp1)
tmp2$cluster_n = tmp2$cluster_n/tmp2$cell_n
tmp2 = tmp2 %>% full_join(unique(metadata[,c("patient","condition","stage","sex","days")]),by=c("patient"="patient","days"="days"))

tmp2 %>% dplyr::select(celltype,cluster_n) %>% distinct() %>% group_by(celltype) %>% summarise(n=sum(cluster_n)) -> celltype_per

tmp2 = tmp2 %>% full_join(celltype_per,by="celltype") %>% mutate(final_n=cluster_n/n)

tmp2$id <- paste0(tmp2$patient,"_",tmp2$days)


##use celltype_pct
library(reshape2)
dcast(tmp2 %>% dplyr::select(id,celltype, final_n),celltype ~ id,mean) -> tmp2_heatmap2
tmp2_heatmap2[tmp2_heatmap2=="NaN"] <-0 
rownames(tmp2_heatmap2) <- tmp2_heatmap2[,1]
tmp2_heatmap2[,-1] -> tmp2_heatmap2

tmp2 %>% dplyr::select(patient,days,sex,id,condition) %>% distinct() -> tmp_col

tmp2_heatmap2[,tmp_col$id] -> tmp2_heatmap2

tmp2_heatmap2[tmp2_heatmap2>0.13] <- 0.13
tmp2_heatmap2[tmp2_heatmap2<0.01] <- 0.01

patient_col <- celltype_color[1:13]

names(patient_col) <- unique(tmp_col$patient)
library(RColorBrewer)
sex_col <- c(brewer.pal(2,"Paired"))[1:2]
names(sex_col) <- unique(tmp_col$sex)
library(pheatmap)
pheatmap(tmp2_heatmap2,
         annotation_col = data.frame(
           Patient=tmp_col$patient,
           Sex= tmp_col$sex,
           Days= tmp_col$days,
           Condition=tmp_col$condition,
           row.names = colnames(tmp2_heatmap2)
         ),
         annotation_colors  = list(Patient=patient_col,Sex=sex_col,Condition=condition_color),
         max.labels = Inf,
         normalize = TRUE,
         show.labels = FALSE,
         fontsize = 12,
         show_colnames = FALSE,
         treeheight_row = 0,
         treeheight_col = 0,
         border_color=NA ,filename ="fig3/new/heatmap_celltypes_labels.pdf" ,height=5,width=13)
pdf("fig3/heatmap_celltypes_labels3.pdf",height=8,width=8)
print(p)
dev.off()


#stage percentage barplot
tmp2
library(ggpubr)
my_comparisons <- list( c("Asymptomatic", "HC"),
                        c("<10", "HC"),
                        c("Severe", "HC"),
                        c("<10", "Severe"),
                        c("Asymptomatic", "Severe"),
                        c("<10", "Asymptomatic"),c("<10",">20"),
                        c("<10","10-20"),
                        c("10-20",">20"),
                        c(">20","Asymptomatic"),
                        c(">20","HC"),
                        c(">20","Severe"),
                        c("10-20","Asymptomatic"),
                        c("10-20","HC"),
                        c("10-20","Severe")
)


ggplot(data=tmp2,aes(x=stage,cluster_n,fill=stage))+
  geom_boxplot(outlier.size = 0.1)+
  geom_jitter(size=0.1)+
  scale_fill_manual(values = cond_palette)+
  facet_wrap(~celltype,nrow=3,scales = "free_y")+
  stat_compare_means(label = "p.signif",method="wilcox.test",comparisons=my_comparisons ,hide.ns=TRUE)+
  theme_cowplot(font_size = 15)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +labs(y="Percentage")
ggsave("fig3/new/stage_allp_bar.pdf",width = 15,height = 12)


#reduce 
my_comparisons <- list( c("Asymptomatic", "HC"),
                        c("<10", "HC"),
                        c("Severe", "HC"),
                        c("<10", "Severe"),
                        c("Asymptomatic", "Severe"),
                        c("<10", "Asymptomatic"),
                        c(">20","Asymptomatic"),
                        c(">20","HC"),
                        c(">20","Severe"),
                        c("10-20","Asymptomatic"),
                        c("10-20","HC"),
                        c("10-20","Severe")
)

ggplot(data=tmp2 %>% filter(!celltype %in% c("Naive CD4+ T","MAIT","Naive CD8+ T","Proliferative T")),aes(x=stage,cluster_n,fill=stage))+
  geom_boxplot(outlier.size = 0.1)+
  geom_jitter(size=0.1)+
  scale_fill_manual(values = cond_palette)+
  facet_wrap(~celltype,nrow=2,scales = "free_y")+
  stat_compare_means(label = "p.signif",method="wilcox.test",comparisons=my_comparisons ,hide.ns=TRUE)+
  theme_cowplot(font_size = 15)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +labs(y="Percentage")
ggsave("fig3/stage_somep_bar.pdf",width = 15,height = 5)




##dotplot
Idents(seu_T_NK) <- seu_T_NK$celltype
featurstnk <- c("CD3D","CD4","GPR183","CCR7","SELL","S100A4","CD8A","GZMA","GZMB","NKG7","LTB","TRDV2","MKI67","TRAV1-2","NCAM1","FCGR3A","KIT","PTGDR2")
DotPlot(seu_T_NK,features=featurstnk,cols="Spectral") +theme(axis.text.x = element_text(hjust = 1,angle = 60))
ggsave("fig3/new/dotplot2.pdf",width = 8,height = 5)
##umap for those genes
FeaturePlot(seu_T_NK,featurstnk,order = TRUE)
ggsave("fig3/umaps_features.png",width=16,height=18)



#correspondence between previous clusters and new clusters
library(tidyverse)
tmp <- seu_T_NK@meta.data
tmp %>% dplyr::select(celltype) -> tmp
tmp$cluster_cell <- seu$cell_type[match(rownames(tmp),colnames(seu))]

tmp %>% group_by(celltype,cluster_cell) %>% summarise(n=n()) -> tmp2
tmp2 %>% full_join(tmp %>% group_by(celltype) %>% summarise(cn=n()),by=c("celltype"="celltype")) %>% mutate(per=n/cn) -> tmp3

library()
reshape2::dcast(tmp3,celltype~cluster_cell,mean)-> tmp4
tmp4[tmp4=="NaN"] <-NA 
rownames(tmp4) <- tmp4[,1]
tmp4[,-1] -> tmp4
#tmp4 <- t(scale(t(tmp4)))
library(pheatmap)
pheatmap(tmp4,
         #max.labels = Inf,
         #normalize = TRUE,
         #scale="column",
         color = colorRampPalette(colors = c("white","blue"))(100),
         fontsize = 10,
         na_col="white",
         cluster_rows = FALSE,
         cluster_cols = FALSE,filename ="fig3/new/celltypes_cor.pdf" ,width = 6,height =5)










