#seu[,seu$celltype %in% c("B cell")] ->seu_b
#remove doubelt because driving clustering
#seu_b <- seu_b[,seu_b$doublet==FALSE]

#seu_b$previous_cluster <- seu$seurat_clusters[seu$celltype %in% c("B cell")]
setwd("~/Dropbox/research/COVID19_MjM/analysis")
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(ggpubr)
library(RColorBrewer)
library(pheatmap)
##### color palette
brewer.pal(n = 4, name ="Dark2")
cond_palette = c("HC"="#1B9E77",
                 "Asymptomatic"="#D95F02",
                 "<10"="#6763A0",
                 "10-20"="#9F8BE0",
                 ">20"="#C5C9FF",
                 "Severe"="#E7298A")

########## 

seu_b = readRDS("seu_bcell.rds")
seu_b$condition[seu_b$condition=="Health control"]="HC"
seu_b <- seu_b[,seu_b$doublet==FALSE]
seu_b <- RunFastMNN(object.list = SplitObject(seu_b, split.by = "Sample"),features = 1000)
seu_b <- RunUMAP(seu_b, reduction = "mnn", dims = 1:10)
seu_b <- FindNeighbors(seu_b, reduction = "mnn", dims = 1:10)

seu_b <- FindClusters(seu_b,resolution = 0.2)
seu_b = ScaleData(seu_b)


DimPlot(seu_b,label = T)
ggsave("figs/UMAP_B_cluster.png",dpi=500,width = 6,height = 5)
seu_b$ct = "Naive B"
seu_b$ct[seu_b$seurat_clusters==1] = "Memory B"
seu_b$ct[seu_b$seurat_clusters==2] = "Plasma B"
DimPlot(seu_b,label = T,label.size = 4,group.by ="ct",cols="alphabet2")+theme(legend.position = "none",axis.ticks = element_blank(),axis.text = element_blank())
ggsave("figs/UMAP_B_celltype.pdf", width = 4.5,height = 4)

seu_b.marker = FindAllMarkers(seu_b,only.pos = T)
#
top10 <- seu_b.marker %>% group_by(cluster) %>% top_n(n = 17, wt = -p_val_adj) %>%  top_n(n = 5, wt = avg_logFC)
ident.1.marker01 = FindMarkers(seu_b,ident.1=2,ident.2=3)
DotPlot(seu_b, features = c(unique(top10$gene),rownames( head(ident.1.marker01[ident.1.marker01$avg_logFC<0,]))),cluster.idents=T)+coord_flip()
ggsave("figs/dotplot_B_marker.pdf",width = 4,height = 7)



saveRDS(seu_b,file = "seu_bcell_cluster.rds")
#seu_b_bc = seu_b
#seu_b = seu_b_bc
# seu_b = readRDS(file="seu_bcell_cluster.rds")

DotPlot(seu_b,features = c("CD79A","MS4A1","CD74","IGHD","CD27","CD38","MZB1","JCHAIN"),group.by = "ct",cols="Spectral",cluster.idents = TRUE) + 
  coord_flip()+theme(axis.text.x = element_text(hjust = 1,angle = 30))
ggsave("figs/dotplot_Bcell_marker.pdf",width = 4.3,height = 4)

get_bcr_clone_df = function(ix){
  filtered_contig_annotations <- read.csv(paste0("~/Dropbox/research/COVID19_MjM/analysis/bcr_filtered/BCR",ix,".filtered_contig_annotations.csv"))
  filtered_contig_annotations = filtered_contig_annotations[filtered_contig_annotations$productive=="True", ]
  filtered_contig_annotations = filtered_contig_annotations[!(filtered_contig_annotations$chain=="Multi"), ]
  filtered_contig_annotations$raw_clonotype_id = paste(ix,filtered_contig_annotations$raw_clonotype_id,sep="_")
  print(table(unique(filtered_contig_annotations$barcode) %in% seu_b$Barcode))
  tmp = filtered_contig_annotations[,c("raw_clonotype_id","chain","v_gene","j_gene")]
  tmp1 = tmp[tmp$chain=="IGK",]
  colnames(tmp1) = c("raw_clonotype_id", "chain", "v_IGK", "j_IGK")
  tmp1 = tmp1[!duplicated(tmp1$raw_clonotype_id),]
  tmp3 = tmp[tmp$chain=="IGL",]
  tmp3 = tmp3[!duplicated(tmp3$raw_clonotype_id),]
  colnames(tmp3) = c("raw_clonotype_id", "chain", "v_IGL", "j_IGL")
  tmp = filtered_contig_annotations[,c("raw_clonotype_id","chain","v_gene","j_gene","c_gene")]
  tmp2 = tmp[tmp$chain=="IGH",]
  tmp2 = tmp2[!duplicated(tmp2$raw_clonotype_id),]
  colnames(tmp2) = c("raw_clonotype_id", "chain", "v_IGH", "j_IGH","c_gene")
  tmp2 = tmp2 %>% full_join(tmp1,by=c("raw_clonotype_id"="raw_clonotype_id"))
  tmp2 = tmp2 %>% full_join(tmp3,by=c("raw_clonotype_id"="raw_clonotype_id"))
  tmp2 = tmp2[,c("raw_clonotype_id","v_IGH","j_IGH","v_IGK","j_IGK","v_IGL","j_IGL","c_gene")]
  
  filtered_contig_annotations = unique(filtered_contig_annotations[,c("barcode","raw_clonotype_id")])
  filtered_contig_annotations = filtered_contig_annotations[!duplicated(filtered_contig_annotations$barcode),]
  filtered_contig_annotations$Sample=ix
  filtered_contig_annotations = filtered_contig_annotations %>% left_join(tmp2,by=c("raw_clonotype_id"="raw_clonotype_id"))
  filtered_contig_annotations = filtered_contig_annotations[!(is.na(filtered_contig_annotations$c_gene)),]
  not_igk = !is.na(filtered_contig_annotations$v_IGK)
  not_igl = !is.na(filtered_contig_annotations$v_IGL)
  filtered_contig_annotations = filtered_contig_annotations[(not_igk | not_igl),]
  filtered_contig_annotations
}





bcr_clone_combined = Reduce(rbind,lapply(2:21,get_bcr_clone_df))
bcr_clone_combined = bcr_clone_combined[!(bcr_clone_combined$Sample==5),] # sample 5 is removed.
bcr_clone_combined = bcr_clone_combined %>% group_by(raw_clonotype_id) %>% mutate(clone_freq=n())
bcr_clone_combined = as.data.frame(bcr_clone_combined)
rownames(bcr_clone_combined) = paste(bcr_clone_combined$Sample,bcr_clone_combined$barcode,sep="-")
bcr_clone_combined = bcr_clone_combined[,!(colnames(bcr_clone_combined) %in% c("Sample","barcode"))]

seu_b = AddMetaData(seu_b, bcr_clone_combined)
seu_b$has_BCR = as.character(!is.na(seu_b$raw_clonotype_id))
seu_b$global_clone_id =  seu_b@meta.data %>% group_indices(v_IGH,v_IGK,v_IGL) 
seu_b$c_gene = substring(seu_b$c_gene,1,4)

seu_b$clonal_type = "Unique"
seu_b$clonal_type[seu_b$clone_freq>1] = "Clonal"
seu_b$clonal_type[is.na(seu_b$clone_freq)] = "Not detected"
seu_b$clonal_type = factor(seu_b$clonal_type,levels = c("Clonal","Unique","Not detected"))
  
seu_b$stage = ">20"
seu_b$stage[seu_b$days<20] = "10-20"
seu_b$stage[seu_b$days<10] = "<10"
seu_b$stage[seu_b$condition=="Asymptomatic"]="Asymptomatic"
seu_b$stage[seu_b$condition=="HC"]="HC"
seu_b$stage[seu_b$condition=="Severe"]="Severe"
seu_b$stage = factor(seu_b$stage,levels = c("HC","Asymptomatic","Severe","<10","10-20",">20"))


######
### B cell abundance
tmp_pct = seu_b@meta.data %>% group_by(patient,days,stage) %>% mutate(sum_c=n()) %>% group_by(patient,days,stage,ct) %>% summarise(ct_pct=n()/min(sum_c))
  
  
  
p1=ggplot(data=tmp_pct[tmp_pct$ct=="Naive B",],aes(x=stage,y=ct_pct,fill=stage))+
  geom_boxplot(outlier.colour = NA,show.legend = F)+
  labs(y="Percentage",title="Naive B")+
  stat_compare_means(label = "..p.signif..",method="wilcox.test",ref.group="HC" ,hide.ns=T,size=7)+
  geom_boxplot(outlier.colour = NA,show.legend = F)+
  scale_fill_manual(values = cond_palette)+
  geom_jitter(show.legend = F,width = 0.1,size=1,color="gray20")+
  theme_bw()+theme(axis.text.x = element_text(hjust = 1,angle = 30))+  theme(legend.position="top",axis.title.x=element_blank(),
                                                                             axis.text.x = element_text(angle = 30, hjust = 1),
                                                                             panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                             panel.background = element_blank(), axis.line = element_line(colour = "black"))
p2=ggplot(data=tmp_pct[tmp_pct$ct=="Memory B",],aes(x=stage,y=ct_pct,fill=stage))+
  geom_boxplot(outlier.colour = NA,show.legend = F)+
  labs(y="",title="Memory B")+
  stat_compare_means(label = "..p.signif..",method="wilcox.test",ref.group="HC" ,hide.ns=T,size=7)+
  geom_boxplot(outlier.colour = NA,show.legend = F)+
  scale_fill_manual(values = cond_palette)+
  geom_jitter(show.legend = F,width = 0.1,size=1,color="gray20")+
  theme_bw()+theme(axis.text.x = element_text(hjust = 1,angle = 30))+  theme(legend.position="top",axis.title.x=element_blank(),
                                                                             axis.text.x = element_text(angle = 30, hjust = 1),
                                                                             panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                             panel.background = element_blank(), axis.line = element_line(colour = "black"))
p3=ggplot(data=tmp_pct[tmp_pct$ct=="Plasma B",],aes(x=stage,y=ct_pct,fill=stage))+
  geom_boxplot(outlier.colour = NA,show.legend = F)+
  stat_compare_means(label = "..p.signif..",method="wilcox.test",ref.group="HC" ,hide.ns=T,size=7,label.y = 0.345)+
  geom_boxplot(outlier.colour = NA,show.legend = F)+
  labs(y="",title="Plasma B")+
  ylim(0,0.368)+
  scale_fill_manual(values = cond_palette)+
  geom_jitter(show.legend = F,width = 0.1,size=1,color="gray20")+
  theme_bw()+theme(axis.text.x = element_text(hjust = 1,angle = 30))+  theme(legend.position="top",axis.title.x=element_blank(),
                                                                             axis.text.x = element_text(angle = 30, hjust = 1),
                                                                             panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                             panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggarrange(p1,p2,p3,ncol=3)
ggsave("figs/B_cell_abundance.pdf",width = 7,height = 3)


######
tmp0 =  seu_b@meta.data %>% group_by(stage,clonal_type) %>% summarise(p_cnt=n())

ggplot(data=tmp0,aes(x=stage,y=p_cnt,fill=clonal_type))+
  geom_bar(stat="identity")

tmp0_m = tmp0[tmp0$stage %in% c("<10","10-20",">20"),]
tmp0_m %<>% group_by(stage) %>% mutate(fraction = p_cnt / sum(p_cnt),
                                   ymax = cumsum(fraction),
                                   ymin = c(0,ymax[1:length(ymax)-1]))
baseNum <- 2
#numCat <- length(unique(dat$ring))
tmp0_m$xmax <- as.numeric(tmp0_m$stage) + baseNum
tmp0_m$xmin = tmp0_m$xmax -0.9


ggplot(tmp0_m, aes(fill=clonal_type,
                     alpha = stage,
                     ymax=ymax, 
                     ymin=ymin, 
                     xmax=xmax, 
                     xmin=xmin)) +
  geom_rect() +
  coord_polar(theta="y")+
  labs(fill="",title="Moderate")+
  scale_alpha_discrete(range = c(0.4,1))+
  scale_fill_manual(values = c("#CF3B34","#379FB8","#BFBFBF"))+
  theme_void()+theme(plot.title = element_text(hjust = 0.5))
ggsave("figs/BCR_clonal_Moderate.pdf",width = 4,height = 3.2)
##

tmp0_m = tmp0[!(tmp0$stage %in% c("<10","10-20",">20")),]
tmp0_m %<>% group_by(stage) %>% mutate(fraction = p_cnt / sum(p_cnt),
                                       ymax = cumsum(fraction),
                                       ymin = c(0,ymax[1:length(ymax)-1]))
baseNum <- 2
#numCat <- length(unique(dat$ring))
tmp0_m$xmax <-  baseNum
tmp0_m$xmin = tmp0_m$xmax -0.9


ggplot(tmp0_m, aes(fill=clonal_type,
                   ymax=ymax, 
                   ymin=ymin, 
                   xmax=xmax, 
                   xmin=xmin)) +
  geom_rect() +
  coord_polar(theta="y")+
  labs(fill="")+
  facet_wrap(~stage)+
  xlim(c(0.7, 2)) +
  scale_fill_manual(values = c("#CF3B34","#379FB8","#BFBFBF"))+
  theme_void()+theme(plot.title = element_text(hjust = 0.5))
ggsave("figs/BCR_clonal_others.pdf",width = 6.5,height = 3.2)

##########
#ig type analysis

tmp = seu_b@meta.data[!is.na(seu_b$raw_clonotype_id),] %>% group_by(patient,days,c_gene) %>% mutate(freq=n())
tmp = tmp %>% group_by(patient,days) %>% mutate(ncnt=n())
tmp = tmp[,c("condition","stage","patient","days","c_gene","freq","ncnt")]
tmp = tmp %>% group_by(condition,stage,patient,days,c_gene) %>% summarise(igh_freq=unique(freq/ncnt))

pdf("figs/BCR_temporal_analysis.pdf")
ggplot(data=tmp,aes(x=condition,y=igh_freq,fill=condition))+
  geom_boxplot(outlier.colour = NA,show.legend = F)+
  geom_jitter(show.legend = F)+
  facet_wrap(~c_gene,ncol = 2)+theme_bw()+theme(axis.text.x = element_text(hjust = 1,angle = 30))

ggplot(data=tmp[tmp$c_gene=="IGHM",],aes(x=days,y=igh_freq,group=patient))+
  geom_point()+
  ggtitle("IGHM")+
  geom_line(linetype="dashed", color="blue")+theme_bw()+facet_wrap(~condition,ncol = 2)+theme_bw()

ggplot(data=tmp[tmp$c_gene=="IGHG",],aes(x=days,y=igh_freq,group=patient))+
  geom_point()+
  ggtitle("IGHG")+
  geom_line(linetype="dashed", color="blue")+theme_bw()+facet_wrap(~condition,ncol = 2)+theme_bw()

ggplot(data=tmp[tmp$c_gene=="IGHM",],aes(x=days,y=igh_freq))+
  geom_point()+
  ggtitle("IGHM")+
  facet_wrap(~patient,ncol = 2)+theme_bw()
ggplot(data=tmp[tmp$c_gene=="IGHG",],aes(x=days,y=igh_freq))+
  geom_point()+
  ggtitle("IGHG")+
  facet_wrap(~patient,ncol = 2)+theme_bw()
ggplot(data=tmp[tmp$c_gene=="IGHD",],aes(x=days,y=igh_freq))+
  geom_point()+
  ggtitle("IGHD")+
  facet_wrap(~patient,ncol = 2)+theme_bw()
ggplot(data=tmp[tmp$c_gene=="IGHA",],aes(x=days,y=igh_freq))+
  geom_point()+
  ggtitle("IGHA")+
  facet_wrap(~patient,ncol = 2)+theme_bw()

ggplot(data=tmp,aes(x=factor(days),y=igh_freq,fill=c_gene))+
  geom_bar(stat = "identity")+
  facet_wrap(~patient,ncol = 2,scales = "free_x")+theme_bw()
dev.off()




tmp_m = tmp


ggplot(data=tmp_m[tmp_m$c_gene=="IGHG",],aes(x=stage,y=igh_freq,fill=stage))+
  geom_boxplot(outlier.colour = NA,show.legend = F)+
  scale_fill_manual(values = cond_palette)+
  geom_jitter(show.legend = F,width = 0.1,size=1,color="gray20")+
  labs(title="IGHG",y="Percentage")+
  stat_compare_means(label = "..p.signif..",method="wilcox.test",ref.group="HC" ,hide.ns=T,size=7,label.y=0.26)+
  theme_bw()+theme(axis.text.x = element_text(hjust = 1,angle = 30))+  theme(legend.position="top",axis.title.x=element_blank(),
                                                                             axis.text.x = element_text(angle = 30, hjust = 1),
                                                                             panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                             panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave("figs/IGHG_summary_by_stage.pdf",width = 4,height = 3)

ggplot(data=tmp_m[tmp_m$c_gene=="IGHM",],aes(x=stage,y=igh_freq,fill=stage))+
  geom_boxplot(outlier.colour = NA,show.legend = F)+
  scale_fill_manual(values = cond_palette)+
  geom_jitter(show.legend = F,width = 0.1,size=1,color="gray20")+
  labs(title="IGHM",y="Percentage")+
  ylim(0.39,0.83)+
  stat_compare_means(label = "..p.signif..",method="wilcox.test",ref.group="HC" ,hide.ns=T,size=7,label.y=0.79)+
  theme_bw()+theme(axis.text.x = element_text(hjust = 1,angle = 30))+  theme(legend.position="top",axis.title.x=element_blank(),
                                                                             axis.text.x = element_text(angle = 30, hjust = 1),
                                                                             panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                             panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave("figs/IGHM_summary_by_stage.pdf",width = 4,height = 3)

ggplot(data=tmp_m[tmp_m$c_gene=="IGHD",],aes(x=stage,y=igh_freq,fill=stage))+
  geom_boxplot(outlier.colour = NA,show.legend = F)+
  scale_fill_manual(values = cond_palette)+
  geom_jitter(show.legend = F,width = 0.1,size=1,color="gray20")+
  ylim(0.,0.32)+
  labs(title="IGHD",y="Percentage")+
  stat_compare_means(label = "..p.signif..",method="wilcox.test",ref.group="HC" ,hide.ns=T,size=7,label.y=0.29)+
  theme_bw()+theme(axis.text.x = element_text(hjust = 1,angle = 30))+  theme(legend.position="top",axis.title.x=element_blank(),
                                                                             axis.text.x = element_text(angle = 30, hjust = 1),
                                                                             panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                             panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave("figs/IGHD_summary_by_stage.pdf",width = 4,height = 3)

ggplot(data=tmp_m[tmp_m$c_gene=="IGHA",],aes(x=stage,y=igh_freq,fill=stage))+
  geom_boxplot(outlier.colour = NA,show.legend = F)+
  scale_fill_manual(values = cond_palette)+
  geom_jitter(show.legend = F,width = 0.1,size=1,color="gray20")+
  labs(title="IGHA",y="Percentage")+
  stat_compare_means(label = "..p.signif..",method="wilcox.test",ref.group="HC" ,hide.ns=T,size=7)+
  theme_bw()+theme(axis.text.x = element_text(hjust = 1,angle = 30))+  theme(legend.position="top",axis.title.x=element_blank(),
                                                                             axis.text.x = element_text(angle = 30, hjust = 1),
                                                                             panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                             panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave("figs/IGHA_summary_by_stage.pdf",width = 4,height = 3)


##########

tmp = seu_b@meta.data[!is.na(seu_b$raw_clonotype_id),]


tmp_ighv = tmp %>% group_by(patient,days) %>% mutate(igh_tot=n()) %>% group_by(condition,patient,days,v_IGH) %>% summarise(igh_prop=n()/min(igh_tot) )

pp = ggplot(data=tmp_ighv,aes(x=condition,y=igh_prop))+
  geom_boxplot()+
  facet_wrap(~v_IGH,ncol=4,scales = "free_y")+
  stat_compare_means(label = "p.signif",method="wilcox.test",ref.group = "HC",hide.ns=T)+
  theme(axis.text.x = element_text(hjust = 1,angle = 30))
ggsave("figs/IGHV_dist.pdf",plot=pp,width = 10,height = 20)

tmp_ighv = tmp[!is.na(tmp$v_IGK),] %>% group_by(patient,days) %>% mutate(igh_tot=n()) %>% group_by(condition,patient,days,v_IGK) %>% summarise(igh_prop=n()/min(igh_tot) )

pp = ggplot(data=tmp_ighv,aes(x=condition,y=igh_prop))+
  geom_boxplot()+
  facet_wrap(~v_IGK,ncol=4,scales = "free_y")+
  stat_compare_means(label = "p.signif",method="wilcox.test",ref.group = "HC",hide.ns=T)+
  theme(axis.text.x = element_text(hjust = 1,angle = 30))
ggsave("figs/v_IGKV_dist.pdf",plot=pp,width = 10,height = 20)


tmp_ighv = tmp[!is.na(tmp$v_IGL),] %>% group_by(patient,days) %>% mutate(igh_tot=n()) %>% group_by(condition,patient,days,v_IGL) %>% summarise(igh_prop=n()/min(igh_tot) )

pp = ggplot(data=tmp_ighv,aes(x=condition,y=igh_prop))+
  geom_boxplot()+
  facet_wrap(~v_IGL,ncol=4,scales = "free_y")+
  stat_compare_means(label = "p.signif",method="wilcox.test",ref.group = "HC",hide.ns=T)+
  theme(axis.text.x = element_text(hjust = 1,angle = 30))
ggsave("figs/v_IGLV_dist.pdf",plot=pp,width = 10,height = 20)



tmp_ighv1 = tmp[!is.na(tmp$v_IGK),] %>% group_by(condition) %>% mutate(igh_tot=n()) %>% group_by(condition,v_IGK,j_IGK) %>% summarise(igh_prop=n()/min(igh_tot) )
colnames(tmp_ighv1) = c("condition","V","J","prop")
#tmp_ighv1$prop = scale(tmp_ighv1$prop)

tmp_ighv2 = tmp[!is.na(tmp$v_IGH),] %>% group_by(condition) %>% mutate(igh_tot=n()) %>% group_by(condition,v_IGH,j_IGH) %>% summarise(igh_prop=n()/min(igh_tot) )
colnames(tmp_ighv2) = c("condition","V","J","prop")
#tmp_ighv2$prop = scale(tmp_ighv2$prop)
#tmp_ighv2$prop[tmp_ighv2$prop>6] = 6

tmp_ighv3 = tmp[!is.na(tmp$v_IGL),] %>% group_by(condition) %>% mutate(igh_tot=n()) %>% group_by(condition,v_IGL,j_IGL) %>% summarise(igh_prop=n()/min(igh_tot) )
colnames(tmp_ighv3) = c("condition","V","J","prop")
#tmp_ighv3$prop = scale(tmp_ighv3$prop)
#tmp_ighv3$prop[tmp_ighv3$prop>4]=4
tmp_vj_comb = rbind(as.data.frame(tmp_ighv1),
                    as.data.frame(tmp_ighv2),
                    as.data.frame(tmp_ighv3))


tmp_a= tmp_vj_comb[tmp_vj_comb$condition=="Asymptomatic",] %>% pivot_wider(id_cols = 'V',names_from="J",values_from="prop",values_fill=0)
tmp_a = as.data.frame(tmp_a)
rownames(tmp_a) = tmp_a$V
tmp_a = tmp_a[,-1]
tmp_a = tmp_a[rowSums(tmp_a)>0.04,colSums(tmp_a)>0.03]
pheatmap(t(tmp_a),cluster_cols = F,cluster_rows = F,
         color = colorRampPalette((brewer.pal(n = 7, name ="Blues")))(100),
         filename = "figs/BCR_VJ_Asy.pdf",
         main="Asymptomatic",
         width=4.7,height = 3)


tmp_a= tmp_vj_comb[tmp_vj_comb$condition=="HC",] %>% pivot_wider(id_cols = 'V',names_from="J",values_from="prop",values_fill=0)
tmp_a = as.data.frame(tmp_a)
rownames(tmp_a) = tmp_a$V
tmp_a = tmp_a[,-1]
tmp_a[tmp_a<0]=0
tmp_a = tmp_a[rowSums(tmp_a)>0.03,colSums(tmp_a)>0.01]
pheatmap(t(tmp_a),cluster_cols = F,cluster_rows = F,
         color = colorRampPalette((brewer.pal(n = 7, name ="Blues")))(100),
         filename = "figs/BCR_VJ_HC.pdf",
         main="Health",
         width=5.1,height = 3)


tmp_a= tmp_vj_comb[tmp_vj_comb$condition=="Moderate",] %>% pivot_wider(id_cols = 'V',names_from="J",values_from="prop",values_fill=0)
tmp_a = as.data.frame(tmp_a)
rownames(tmp_a) = tmp_a$V
tmp_a = tmp_a[,-1]
tmp_a = tmp_a[rowSums(tmp_a)>0.03,colSums(tmp_a)>0.01]
pheatmap(t(tmp_a),cluster_cols = F,cluster_rows = F,
         color = colorRampPalette((brewer.pal(n = 7, name ="Blues")))(100),
         filename = "figs/BCR_VJ_Mod.pdf",
         main="Moderate",
         width=5.1,height = 3)


tmp_a= tmp_vj_comb[tmp_vj_comb$condition=="Severe",] %>% pivot_wider(id_cols = 'V',names_from="J",values_from="prop",values_fill=0)
tmp_a = as.data.frame(tmp_a)
rownames(tmp_a) = tmp_a$V
tmp_a = tmp_a[,-1]
tmp_a = tmp_a[rowSums(tmp_a)>0.04,colSums(tmp_a)>0.02]
pheatmap(t(tmp_a),cluster_cols = F,cluster_rows = F,
         color = colorRampPalette((brewer.pal(n = 7, name ="Blues")))(100),
         filename = "figs/BCR_VJ_Sev.pdf",
         main="Severe",
         width=4.7,height = 3)



tmp = seu_b@meta.data[!is.na(seu_b$raw_clonotype_id),]

tmp_ighv2 = tmp[!is.na(tmp$v_IGH),] %>% group_by(condition,patient,days) %>% mutate(igh_tot=n()) %>% group_by(condition,patient,days,v_IGH,j_IGH) %>% summarise(igh_prop=n()/min(igh_tot) )
#tmp_ighv2_j4 = tmp_ighv2[tmp_ighv2$j_IGH=="IGHJ4",]
tmp_sum = tmp_ighv2_j4 %>% group_by(v_IGH) %>% summarise(prop_sum=sum(igh_prop)) %>% top_n(10,wt=prop_sum)
tmp_ighv2_j4 = tmp_ighv2_j4[tmp_ighv2_j4$v_IGH %in% tmp_sum$v_IGH,]

my_comparisons <- list( c("Asymptomatic", "HC"),
                        c("Moderate", "HC"),
                        c("Severe", "HC"),
                        c("Moderate", "Severe"),
                        c("Asymptomatic", "Severe"),
                        c("Moderate", "Asymptomatic"))

ggplot(data=tmp_ighv2_j4,aes(x=condition,y=igh_prop,fill=condition))+
  geom_boxplot()+
  geom_jitter()+
  scale_fill_brewer(palette = "Dark2")+
  facet_wrap(~v_IGH,ncol=3,scales = "free_y")+
  stat_compare_means(label = "p.signif",method="t.test",comparisons=my_comparisons ,hide.ns=T)+
  theme(axis.text.x = element_text(hjust = 1,angle = 30))
ggsave("figs/BCR_VJ_IGHJ4_withtest.pdf")

ggplot(data=tmp_ighv2_j4,aes(x=condition,y=igh_prop,fill=condition))+
  geom_boxplot()+
  geom_jitter(width=0.1)+
  facet_wrap(~v_IGH,ncol=3)+
  theme_bw()+
  theme(axis.text.x = element_text(hjust = 1,angle = 30))
ggsave("figs/BCR_VJ_IGHJ4_notest.pdf",width = 6,height = 6)

ggplot(data=tmp_ighv2_j4,aes(x=days,y=igh_prop,col=patient,group=patient))+
  geom_point()+
  #ggtitle("IGHM")+
  geom_line(linetype="dashed", color="blue")+theme_bw()+facet_wrap(~v_IGH,ncol = 2)+theme_bw()
ggsave("figs/BCR_VJ_IGHJ4_byday.pdf")


####################   
### DE gene analysis.
library(topGO)
library(edgeR)
library(limma)
library(Glimma)
library(ggrepel)
library(org.Hs.eg.db)
library(AnnotationHub)
ah <- AnnotationHub()
EnsDb.Hsapiens.v98 <- query(ah, c("EnsDb", "Homo Sapiens", 98))[[1]]

# Add chromosome location so we can filter on mitochondrial genes.

rowData(summed)$gene_biotype =biotype
##
get_go_table = function(results,direction){
  gene_vec = as.vector(results@.Data)
  names(gene_vec) = rownames(results@.Data)

  sampleGOdata <- new("topGOdata",
                      description = "test", ontology = "BP",
                      allGenes = gene_vec, geneSelectionFun = function(a_g){a_g == direction},
                      nodeSize = 5,
                      annot=annFUN.org, mapping="org.Hs.eg.db", ID = "SYMBOL")

  res = runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")

  tabl = GenTable(sampleGOdata, classic = res,
                  orderBy = "weight", ranksOf = "classic", topNodes = 20)

  ann.genes = function(x){genesInTerm(sampleGOdata, whichGO=tabl$GO.ID[x])[[1]][scoresInTerm(sampleGOdata, whichGO=tabl$GO.ID[x])[[1]] != 0]}
  tabl$genes = unlist(lapply(lapply(1:20,ann.genes),function(x){paste(x[x %in% names(gene_vec)[gene_vec == direction]],collapse = ',')}))
  return(tabl)
}
##
#QC_genes = c("HBB","DDX3X","AL031123.2","AC103591.3","AC105285.1","EIF1AY","EEF1A1","AL391121.1","AC008735.2","AD000090.1","EEF1G","XIST","AC233755.1","AC136616.2","AC245060.6")
############
## single cell DE


seu_b_sel = seu_b[,seu_b$condition=="Severe"]
seu_b_sel = seu_b_sel[rowSums(as.matrix(seu_b_sel@assays$RNA@counts))>30,]
seu_b_sel = seu_b_sel[!(rownames(seu_b_sel) %in% QC_genes),]

d = seu_b_sel$days
ct = seu_b_sel$ct

design_mat = model.matrix(~d+ct)

allcounts = DGEList(counts=as.matrix(seu_b_sel@assays$RNA@counts)) # may need a lot of RAM
allcounts <- calcNormFactors(allcounts)
allcounts = estimateDisp(allcounts, design_mat,robust=TRUE)

fit = glmQLFit(allcounts, design_mat)

lrt = glmQLFTest(fit, coef=2)
results <- decideTests(lrt,p.value=0.05)

glMDPlot(lrt,counts=cpm(allcounts,log=T,prior.count = 1),groups=as.factor(d),status=results[,1],html="B_days_Severe",launch=FALSE)

tabl=get_go_table(results,1)
print(head(tabl))

tabl2=get_go_table(results,-1)
print(head(tabl2))


seu_b_sel = seu_b[,seu_b$condition=="Moderate"]
seu_b_sel = seu_b_sel[,seu_b_sel$patient != "G"]
seu_b_sel = seu_b_sel[rowSums(as.matrix(seu_b_sel@assays$RNA@counts))>30,]
seu_b_sel = seu_b_sel[!(rownames(seu_b_sel) %in% QC_genes),]

d = seu_b_sel$days
ct = seu_b_sel$ct
pt = seu_b_sel$patient
m = seu_b_sel@assays$RNA@scale.data["MALAT1",]
c = seu_b_sel@assays$RNA@scale.data["CD74",]
design_mat = model.matrix(~d+m+c+ct+pt)

allcounts = DGEList(counts=as.matrix(seu_b_sel@assays$RNA@counts)) # may need a lot of RAM
allcounts <- calcNormFactors(allcounts)
allcounts = estimateDisp(allcounts, design_mat,robust=TRUE)

fit = glmQLFit(allcounts, design_mat)

lrt = glmQLFTest(fit, coef=2)
results <- decideTests(lrt,p.value=0.01,lfc=0.015)
cpm_cnt = cpm(allcounts,log=T,prior.count = 1)

glMDPlot(lrt,counts=,groups=as.factor(d),status=results[,1],html="B_days_Moderate",launch=FALSE)

tp = topTags(lrt,n=Inf)@.Data[[1]]

ggplot(data=tp,aes(x=logFC,y=-log10(FDR)))+
  geom_point(alpha=0.7,size=1)+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))


tabl=get_go_table(results,1)
print(head(tabl))

tabl2=get_go_table(results,-1)
print(head(tabl2))

plot_df = data.frame(days=d,gene_expr=cpm_cnt["IL6",],patient=pt,celltype=ct)
ggplot(data=plot_df,aes(x=days,y=gene_expr))+
  geom_point()+
  geom_smooth()+
  facet_wrap(~celltype,ncol=1)


############
# pseudo-bulk DE

hist(log2(rowSums(seu_b@assays$RNA@counts)))

seu_b_filter = seu_b[rowSums(as.matrix(seu_b@assays$RNA@counts))>30,]
tmp_bigmat = as.matrix(seu_b_filter@assays$RNA@counts)
tmp_bigmat = as.data.frame(tmp_bigmat)
all_cell_id = colnames(tmp_bigmat)
tmp_bigmat$genename=rownames(tmp_bigmat)
all_meta = seu_b_filter@meta.data
all_meta$cell_id = rownames(all_meta)
all_meta$ct = make.names(all_meta$ct)
#rm(seu_b)
rm(seu_b_filter)  # dont have to run this. it is just to save RAM

pt_cond = unique(all_meta[,c("patient","condition")])

pseudobulk_mat = tmp_bigmat %>% pivot_longer(all_cell_id,names_to = "cell_id", values_to = "cnt")
pseudobulk_mat %<>% left_join(all_meta[,c("cell_id","patient","days","ct","condition")],by=c("cell_id"="cell_id")) %>% group_by(genename,patient,days,ct,condition) %>% summarise(gene_cnt=sum(cnt))
pseudobulk_mat = pseudobulk_mat %>% pivot_wider(id_cols = genename,names_from=c(patient,days,condition,ct),values_from=gene_cnt,values_fill=0)
pseudobulk_mat = as.data.frame(pseudobulk_mat)
rownames(pseudobulk_mat) = pseudobulk_mat$genename
pseudobulk_mat = pseudobulk_mat[,-1]

# remove genes related to sex or quality of cell
biotype<- mapIds(
  x = EnsDb.Hsapiens.v98, 
  # NOTE: Need to remove gene version number prior to lookup.
  keys = rownames(pseudobulk_mat),
  keytype = "SYMBOL",
  column = "GENEBIOTYPE")
pseudobulk_mat = pseudobulk_mat[biotype=="protein_coding" & !is.na(biotype),]

sel_ct="Memory.B"
#sel_ct="Naive.B"
#sel_ct="Plasma.B"
#QC_genes = c("MALAT1","FOSB","FOS","JUNB","UTY","HBB","DDX3X","DDX3Y","KDM5D","ZFY","AL031123.2","AC103591.3","AC105285.1","EIF1AY","EEF1A1","AL391121.1","AC008735.2","AD000090.1","EEF1G","XIST","AC233755.1","AC136616.2","AC245060.6")
#pseudobulk_sel = pseudobulk_mat[!(rownames(pseudobulk_mat) %in% QC_genes),grepl(sel_ct,colnames(pseudobulk_mat))]
pseudobulk_sel = pseudobulk_sel[rowSums(pseudobulk_sel)>15,]
pt = sapply(strsplit(colnames(pseudobulk_sel),"_"), function(x){x[1]})
cond = sapply(strsplit(colnames(pseudobulk_sel),"_"), function(x){x[3]})
dy = as.numeric(sapply(strsplit(colnames(pseudobulk_sel),"_"), function(x){x[2]}))


### by condition
design_mat = model.matrix(~0+cond)
contr.matrix <- makeContrasts(
  MvsA = condModerate-condAsymptomatic,
  MvsS = condModerate-condSevere,
  AvsS = condAsymptomatic-condSevere,
  MvsH = condModerate-condHC,
  AvsH = condAsymptomatic-condHC,
  SvsH = condSevere-condHC,
  levels = colnames(design_mat))


allcounts = DGEList(counts=pseudobulk_sel) # may need a lot of RAM
allcounts <- calcNormFactors(allcounts)
allcounts = estimateDisp(allcounts, design=design_mat, robust=TRUE)

fit = glmQLFit(allcounts, design_mat)

FDR_cutoff=0.05
print_top_genes = 10

for (i in 1:ncol(contr.matrix)){
  lrt = glmQLFTest(fit, contrast=contr.matrix[,i])
  results <- decideTests(lrt,p.value=FDR_cutoff,lfc = 1)
  comp_name=colnames(contr.matrix)[i]
  print(comp_name)
  print(table(results))
  # if(sum(results==1)>10){ # only do GO if there are more than 10 sig DE genes
  #   tabl=get_go_table(results,1)
  #   print(head(tabl))
  #   write.csv(tabl,file = paste0("processed_files/",sel_ct,"_",
  #                                paste(comp_name,"upin",substr(comp_name,1,1),sep = "_"),
  #                                "_GOtable.csv"))
  # }
  # if(sum(results==(-1))>10){ # only do GO if there are more than 10 sig DE genes
  #   tabl=get_go_table(results,-1)
  #   print(head(tabl))
  #   write.csv(tabl,file = paste0("processed_files/",sel_ct,"_",
  #                                paste(comp_name,"upin",substr(comp_name,4,4),sep = "_"),
  #                                "_GOtable.csv"))
  # }
  tp = topTags(lrt,n=Inf)@.Data[[1]]
  write.csv(tp,file=paste0("processed_files/",sel_ct,"_",comp_name,"_DGEtable.csv"))
  tp$genename=rownames(tp)
      if(sum(results!=0)>3){  ## make plot if more than 3 DE genes       
        tp$sig="no"
        tp$sig[tp$FDR<FDR_cutoff & tp$logFC>1] = "pos"
        tp$sig[tp$FDR<FDR_cutoff & tp$logFC<(-1)] = "neg"
        if(sum(results!=0)<print_top_genes){
          tp$genename[tp$FDR>FDR_cutoff] = NA
        }else{
          tp$genename[order(tp$FDR)[(print_top_genes+1):nrow(tp)]] = NA
        }
        p = ggplot(data=tp,aes(x=logFC,y=-log10(FDR),col=sig,label=genename))+
          geom_point(alpha=0.6,size=0.5,show.legend = F)+
          ggtitle(comp_name)+
          geom_text_repel(show.legend = F)+
          scale_color_manual(values=c("no"="grey30","pos"="#FF3E00","neg"="#6796FF"))+
          theme_bw()+
          theme(plot.title = element_text(hjust = 0.5),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black"))
        ggsave(paste0("figs/DEfig/",sel_ct,"_",comp_name,"_volcano.pdf"),plot=p,width = 3,height =4)
      }
}
groups = c("AvsH"="Asymptomatic-HC",
           "AvsS"="Asymptomatic-Severe",
           "MvsA"="Moderate-Asymptomatic",
           "MvsH"="Moderate-HC",
           "MvsS"="Moderate-Severe",
           "SvsH"="Severe-HC")

get_toptb = function(grp,sel_ct){
  B_MvsA_DGEtable <- read_csv(paste0("processed_files/",sel_ct,"_",grp,"_DGEtable.csv"))
  B_MvsA_DGEtable = B_MvsA_DGEtable[B_MvsA_DGEtable$logCPM>5 & B_MvsA_DGEtable$FDR<0.05,]
  B_MvsA_DGEtable_pos = B_MvsA_DGEtable[B_MvsA_DGEtable$logFC>0,]
  B_MvsA_DGEtable_pos=head(B_MvsA_DGEtable_pos[order(B_MvsA_DGEtable_pos$logFC,decreasing = T),],n=10)
  B_MvsA_DGEtable_neg = B_MvsA_DGEtable[B_MvsA_DGEtable$logFC<0,]
  B_MvsA_DGEtable_neg=head(B_MvsA_DGEtable_neg[order(B_MvsA_DGEtable_neg$logFC,decreasing = F),],n=10)
  rbind(B_MvsA_DGEtable_pos,B_MvsA_DGEtable_neg)
}

combined_df1 = Reduce(rbind,lapply(names(groups),get_toptb,sel_ct=sel_ct))

combined_gene = unique(combined_df1$X1)
cpm_cnt = cpm(allcounts,log=T,prior.count = 1)
ph_anno = data.frame(conditions=cond,days=dy)
rownames(ph_anno) = colnames(cpm_cnt)
tmp_cnt = t(scale(t(cpm_cnt[combined_gene,])))
tmp_cnt[tmp_cnt>2.5]=2.5
tmp_cnt[tmp_cnt<(-2.5)]=-2.5

anno_color=list(conditions=brewer.pal(n = 4, name ="Dark2"))
names(anno_color$conditions) = c( "HC","Asymptomatic","Moderate", "Severe")
pheatmap::pheatmap(tmp_cnt[,],cluster_cols = T,cluster_rows = T,
                   show_colnames = F,annotation_col = ph_anno,
                   annotation_colors = anno_color,treeheight_row=0,treeheight_col=0,
                   filename = paste0("figs/",sel_ct,"_DE_heatmap.pdf"),width = 3.5,height = 4.5,
                   fontsize = 5,border_color=NA)


#######
# time course analysis
# 
# sel_col = (pt %in% c("B","C","I","K","L","M"))
# 
# pseudobulk_mat_sel = pseudobulk_mat[,!sel_col]
# 
# timep = as.numeric( sapply(strsplit(colnames(pseudobulk_mat_sel),"_"), function(x){x[2]}))
# pt = sapply(strsplit(colnames(pseudobulk_mat_sel),"_"), function(x){x[1]})
# ct = sapply(strsplit(colnames(pseudobulk_mat_sel),"_"), function(x){x[3]})
# 
# 
# design_mat = model.matrix(~timep+pt+ct)
# 
# allcounts = DGEList(counts=pseudobulk_mat_sel) # may need a lot of RAM
# allcounts <- calcNormFactors(allcounts)
# allcounts = estimateDisp(allcounts, design=design_mat, robust=TRUE)
# 
# fit = glmQLFit(allcounts, design_mat)
# 
# lrt = glmQLFTest(fit, coef=2)
# results <- decideTests(lrt,p.value=FDR_cutoff)
# 
# glMDPlot(lrt,counts=cpm(allcounts,log=T,prior.count = 1),groups=timep,status=results[,1],html="B_timecouse",launch=FALSE)



#p1 = VlnPlot(seu_b,features=c("IFITM1","ISG15","IFI6","XAF1"),group.by = "stage",cols=cond_palette,pt.size=0.1,ncol=2)+theme(axis.title = element_blank())

p1 = VlnPlot(seu_b,features=c("CD69"),group.by = "stage",cols=cond_palette,pt.size=0.)+theme(axis.text.x = element_blank(),axis.title = element_blank())
p2 = VlnPlot(seu_b,features=c("ID3"),group.by = "stage",cols=cond_palette,pt.size=0)+theme(axis.text.x = element_blank(),axis.title = element_blank())
p2 = VlnPlot(seu_b,features=c("IFITM1"),group.by = "stage",cols=cond_palette,pt.size=0)+theme(axis.text.x = element_blank(),axis.title = element_blank())
p3 = VlnPlot(seu_b,features=c("ISG15"),group.by = "stage",cols=cond_palette,pt.size=0)+theme(axis.text.x = element_blank(),axis.title = element_blank())
p4 = VlnPlot(seu_b,features=c("XAF1"),group.by = "stage",cols=cond_palette,pt.size=0)+theme(axis.text.x = element_blank(),axis.title = element_blank())
p5 = VlnPlot(seu_b,features=c("STAT1"),group.by = "stage",cols=cond_palette,pt.size=0.0)+theme(axis.text.x = element_blank(),axis.title = element_blank())
p6 = VlnPlot(seu_b,features=c("IFI44L"),group.by = "stage",cols=cond_palette,pt.size=0.0)+theme(axis.text.x = element_blank(),axis.title = element_blank())
p7 = VlnPlot(seu_b,features=c("TSC22D3"),group.by = "stage",cols=cond_palette,pt.size=0.0)+theme(axis.text.x = element_blank(),axis.title = element_blank())
p8 = VlnPlot(seu_b,features=c("PMAIP1"),group.by = "stage",cols=cond_palette,pt.size=0.0)+theme(axis.text.x = element_blank(),axis.title = element_blank())
ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,ncol=4,nrow=2,common.legend = T,legend="right",align="h")
ggsave("figs/B_IFN1gene_vlnplot.pdf",width = 9,height = 4)

## make GO plot

GOtable <- read.csv("~/Dropbox/research/COVID19_MjM/analysis/processed_files/B_SvsH_upin_S_GOtable.csv", row.names=1)
GOtable=GOtable[order(GOtable$classic)[1:10],]
GOtable$Term = factor(GOtable$Term,levels = rev(GOtable$Term))

ggplot(data=GOtable,aes(x=Term,y=-log10(classic),fill=-log10(classic)))+
  geom_bar(stat="identity",show.legend=F)+
  scale_fill_continuous(type = "viridis")+
  labs(y="-log10(P value)",x="GO Biological Processes")+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                   panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  coord_flip()
ggsave("figs/SvsH_GO_barplot.pdf",width = 4.5,height =2.5)


groups = c("AvsH"="Asymptomatic-HC",
           "AvsS"="Asymptomatic-Severe",
           "MvsA"="Moderate-Asymptomatic",
           "MvsH"="Moderate-HC",
           "MvsS"="Moderate-Severe")


for (i in 1:length(groups)) {
  g =  names(groups[i])
  go_f = paste0("~/Dropbox/research/COVID19_MjM/analysis/processed_files/B_",g,"_upin_",substr(g,1,1),"_GOtable.csv")
  if(file.exists(go_f)){
    print(go_f)
    GOtable <- read.csv(go_f, row.names=1)
    GOtable = GOtable[!duplicated(GOtable$Term),]
    GOtable=GOtable[order(GOtable$classic)[1:6],]
    GOtable$Term = factor(GOtable$Term,levels = rev(GOtable$Term))
    ggplot(data=GOtable,aes(x=Term,y=-log10(classic),fill=-log10(classic)))+
      geom_bar(stat="identity",show.legend=F)+
      scale_fill_continuous(type = "viridis")+
      labs(y="-log10(P value)",x="GO Biological Processes",title=paste0(groups[i],"upin",substr(g,1,1) ))+
      theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                       panel.background = element_blank(), axis.line = element_line(colour = "black"))+
      coord_flip()
    ggsave(paste0("figs/DEfig/",g,"upin",substr(g,1,1),"_GO_barplot.pdf"),width = 4.5,height =2.5)
  }
  go_f = paste0("~/Dropbox/research/COVID19_MjM/analysis/processed_files/B_",g,"_upin_",substr(g,4,4),"_GOtable.csv")
  if(file.exists(go_f)){
    print(go_f)
    GOtable <- read.csv(go_f, row.names=1)
    GOtable = GOtable[!duplicated(GOtable$Term),]
    GOtable=GOtable[order(GOtable$classic)[1:6],]
    GOtable$Term = factor(GOtable$Term,levels = rev(GOtable$Term))
    ggplot(data=GOtable,aes(x=Term,y=-log10(classic),fill=-log10(classic)))+
      geom_bar(stat="identity",show.legend=F)+
      scale_fill_continuous(type = "viridis")+
      labs(y="-log10(P value)",x="GO Biological Processes",title=paste0(groups[i],"upin",substr(g,4,4) ))+
      theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                       panel.background = element_blank(), axis.line = element_line(colour = "black"))+
      coord_flip()
    ggsave(paste0("figs/DEfig/",g,"upin",substr(g,4,4),"_GO_barplot.pdf"),width = 4.5,height =2.8)
  }
}





