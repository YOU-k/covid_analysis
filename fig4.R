#vdj analysis
#t cells
setwd("/stornext/General/data/user_managed/grpu_mritchie_1/Yue/COVID19_MjM/report/vdj")

dir()
fileid <- c(2:4,6:21)
tra_tib <- tibble()
trb_tib <- tibble()
for (i in 1:length(fileid)) {
  read.csv(paste0("TCR",fileid[i],".filtered_contig_annotations.csv")) -> tcr
  tcr$barcode <- paste0(fileid[i],"-",tcr$barcode)
  
  tra <- tcr[tcr$chain == "TRA", ]
  tra_tib <- bind_rows(tra_tib,tra)
  trb <- tcr[tcr$chain == "TRB", ]
  trb_tib <- bind_rows(trb_tib,trb)
}

tra_tib <- tra_tib[!duplicated(tra_tib$barcode),]
rownames(tra_tib) <- tra_tib$barcode

trb_tib <- trb_tib[!duplicated(trb_tib$barcode),]
rownames(trb_tib) <- trb_tib$barcode

####
tra_tbl <- tibble(
  cluster = factor(levels(seu$seurat_clusters), levels(seu$seurat_clusters)),
  n_cells = as.vector(table(seu$seurat_clusters)),
  prop = as.vector(table(seu$seurat_clusters[lengths(TRA) > 0]) / n_cells),
  chain = "a")

trb_tbl <- tibble(
  cluster = factor(levels(seu$seurat_clusters), levels(seu$seurat_clusters)),
  n_cells = as.vector(table(seu$seurat_clusters)),
  prop = as.vector(table(seu$seurat_clusters[lengths(TRB) > 0]) / n_cells),
  chain = "b")

setwd("/stornext/General/data/user_managed/grpu_mritchie_1/Yue/COVID19_MjM/result/vdj_tcells")
pdf("percentage_ab_cluster.pdf")
ggplot(rbind(tra_tbl, trb_tbl)) + 
  geom_col(aes(x = cluster, y = prop, fill = cluster)) + 
  facet_grid(~ chain) +
  #scale_fill_manual(values = cluster_colours, name = "cluster") +
  theme_cowplot(font_size = 8) +
  ylim(0, 1) +
  ylab("Proportion") +
  guides(fill = FALSE) +
  ggtitle("TCR component")
dev.off()
#####
saveRDS(trb_tib,"trb_tib.rds")
saveRDS(tra_tib,"tra_tib.rds")


##only keep unique tra-trb pair

trb_tib[match(rownames(tra_tib),rownames(trb_tib)),] ->trb_tib # contain a but not b
trb_tib <- trb_tib[!is.na(trb_tib$barcode),]
rownames(trb_tib) <- trb_tib$barcode
#AddMetaData(seu,trb_tib) ->seu2

AddMetaData(seu_T_NK,trb_tib) ->seutnk

tra_tib[match(rownames(trb_tib),rownames(tra_tib)),] ->tra_tib 
rownames(tra_tib) <- tra_tib$barcode
library(tidyverse)
tra_tib %>% mutate(v_gene2=v_gene,d_gene2=d_gene,j_gene2=j_gene) %>% dplyr::select(v_gene2,d_gene2,j_gene2) ->tra_genes
rownames(tra_genes) <- rownames(tra_tib)

AddMetaData(seutnk,tra_genes) ->seutnk

png("trb_umap_tnk.png",width = 10, height = 8,units="in",res=200)
Idents(seutnk) <- seutnk$chain
UMAPPlot(seutnk,cols="pink")
dev.off()


#check clonetypes
seutnk$raw_clonotype_id <- paste0(seutnk$Sample,"-",seutnk$raw_clonotype_id)
seutnk@meta.data %>% group_by(v_gene,d_gene,j_gene,v_gene2,d_gene2,j_gene2) %>% summarise(n=n())-> clonoinfo



clonoinfo[clonoinfo$n>200,] 
#check where are they

clonoinfo %>% dplyr::filter(!is.na(v_gene)) -> clonoinfo

clonoinfo$assigned_clonotype <- paste0("clonotype-",c(1:nrow(clonoinfo)))
clonoinfo[clonoinfo$n>170,]-> highclono

left_join(seutnk@meta.data,clonoinfo,by = c("v_gene", "d_gene", "j_gene", "v_gene2", "d_gene2", "j_gene2")) ->assigned_clonotype

seutnk$all_clonotype <- assigned_clonotype$assigned_clonotype

seutnk$high_clonal <- seutnk$all_clonotype 
seutnk$high_clonal[!seutnk$high_clonal %in% highclono$assigned_clonotype] <- NA


library(wesanderson)
plot_grid(p1,p2)


#percentage plots 
seutnk@meta.data -> meta2

#noclonal
assigned_clonotype$assigned_clonotype[assigned_clonotype$n==1] -> noClonal
noClonal[!is.na(noClonal)] -> noClonal
meta2$clonal_status[meta2$all_clonotype %in% noClonal] <- "NoClonal"
meta2$clonal_status[is.na(meta2$all_clonotype)] <- "NoVDJ"
meta2$clonal_status[is.na(meta2$clonal_status)] <- "Clonal"

#distribution of clone status
library(RColorBrewer)
colors <- rev(brewer.pal(n=3,"Blues"))
names(colors) <- c("Clonal","NoClonal", "NoVDJ" )

setwd("/stornext/General/data/user_managed/grpu_mritchie_1/Yue/COVID19_MjM/result/vdj_tcells/")
pdf("percentage_clonal_status.pdf")
ggplot(meta2) +
  geom_bar(aes(seurat_clusters,fill=clonal_status),position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values=colors) +theme_bw()
dev.off()


n <- length(unique(metavdj$patient))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector[4] <- "lightblue"
patient_color <- col_vector 


#compare trbv genes

left_join(metavdj %>% group_by(condition,days,v_gene) %>% summarise(n=n()), metavdj %>% group_by(condition,days) %>% summarise(con_daysn=n())) %>% 
  mutate(percent=n/con_daysn) %>% 
  ggplot() +geom_boxplot(aes(x=condition,y=percent),col="grey") + 
  geom_point(aes(x=condition,y=percent,col=days))  +facet_wrap(.~v_gene) +theme_bw()

#incorporate stage
metavdj$stage[metavdj$days<=0] <- "stage1"
metavdj$stage[metavdj$days>0 & metavdj$days<=5] <- "stage2"
metavdj$stage[metavdj$days>5 & metavdj$days<=10] <- "stage3"
metavdj$stage[metavdj$days>10 & metavdj$days<=20] <- "stage4"
metavdj$stage[metavdj$days>20] <- "stage5"


colors <- brewer.pal(n=6,"Greens")[-1]

pdf("TRBV_stage.pdf",width=15,height = 15)
left_join(metavdj %>% group_by(condition,stage,v_gene) %>% summarise(n=n()), metavdj %>% group_by(condition,stage) %>% summarise(con_stagen=n())) %>% 
  mutate(percent=n/con_stagen) %>% filter(!v_gene=="None") %>% 
  ggplot() +geom_boxplot(aes(x=condition,y=percent)) + 
  geom_point(aes(x=condition,y=percent,col=stage))  +facet_wrap(.~v_gene) +theme_bw() + scale_color_manual(values = cond_palette) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        text = element_text(size=15))
dev.off()


#do not care time
left_join(metavdj %>% group_by(condition,v_gene) %>% summarise(n=n()), metavdj %>% group_by(condition) %>% summarise(con=n())) %>% 
  mutate(percent=n/con) %>% filter(!v_gene=="None") %>% 
  ggplot()  + 
  geom_col(aes(x=condition,y=percent,fill=condition))  +facet_wrap(.~v_gene) +theme_bw() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        text = element_text(size=15))


pdf("TRBJ_stage.pdf",width=8,height = 8)
left_join(metavdj %>% group_by(condition,stage,j_gene) %>% summarise(n=n()), metavdj %>% group_by(condition,stage) %>% summarise(con_stagen=n())) %>% 
  mutate(percent=n/con_stagen) %>% filter(!j_gene=="None") %>% 
  ggplot() +geom_boxplot(aes(x=condition,y=percent)) + 
  geom_point(aes(x=condition,y=percent,col=stage))  +facet_wrap(.~j_gene) +theme_bw() + scale_color_manual(values = cond_palette) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        text = element_text(size=15))
dev.off()

saveRDS(meta2,"vdj_tcells/meta2.rds")
saveRDS(clonoinfo,"vdj_tcells/clonoinfo.rds")


#barplots for each cell
assigned_clonotype$assigned_clonotype[assigned_clonotype$n==1] -> noClonal
noClonal[!is.na(noClonal)] -> noClonal
meta2$clonal_status[meta2$all_clonotype %in% noClonal] <- "Unique"
meta2$clonal_status[is.na(meta2$all_clonotype)] <- "Not detected"
meta2$clonal_status[is.na(meta2$clonal_status)] <- "Clonal"

colors <- rev(brewer.pal(n=3,"Greens"))
names(colors) <- c("Clonal","Unique","Not detected")
####
#meta2$clonal_status <- factor(meta2$clonal_status,levels = c("Clonal","Unique","Not detected"))

meta2$clonal_status <- factor(meta2$clonal_status,levels = c("Clonal","Unique","Not detected"))
meta2 %>% full_join(meta2 %>% group_by(celltype) %>% summarise(n=n()),by="celltype") %>%
  ggplot() +
  geom_bar(aes(reorder(celltype , n),fill=clonal_status),position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values=colors) +theme_bw() + coord_flip()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(hjust = 1,angle = 60),
        text = element_text(size=12))
ggsave("fig4/distribution_of_clonal_status_celltype.pdf",width = 6,height = 4)


meta2[!meta2$clonal_status=="Not detected",] ->metavdj


tmp1 = metavdj %>% group_by(condition) %>% mutate(igh_tot=n()) %>% group_by(condition,v_gene,j_gene) %>% summarise(igh_prop=n()/min(igh_tot) )
colnames(tmp1) = c("condition","V","J","prop")

tmp2 = metavdj %>% group_by(condition) %>% mutate(igh_tot=n()) %>% group_by(condition,v_gene2,j_gene2) %>% summarise(igh_prop=n()/min(igh_tot) )
colnames(tmp2) = c("condition","V","J","prop")


tmp_vj_comb = rbind(as.data.frame(tmp1),
                    as.data.frame(tmp2))

#heatmap_as
tmp_a= tmp_vj_comb[tmp_vj_comb$condition=="Asymptomatic",] %>% filter(!(V=="None"|J=="None")) %>% pivot_wider(id_cols = 'V',names_from="J",values_from="prop",values_fill = c(prop=0))
tmp_a = as.data.frame(tmp_a)
rownames(tmp_a) = tmp_a$V
tmp_a = tmp_a[,-1]
tmp_a = tmp_a[rowSums(tmp_a)>0.04,colSums(tmp_a)>0.03]
library(pheatmap)

pheatmap(tmp_a,cluster_cols = F,cluster_rows = F,
         color = colorRampPalette((brewer.pal(n = 7, name ="Blues")))(100),
         filename = "fig4/TCR_VJ_Asy.pdf",
         main="Asymptomatic",
         width=6,height =3)

#heatmap_hc

tmp_a= tmp_vj_comb[tmp_vj_comb$condition=="HC",] %>% filter(!(V=="None"|J=="None")) %>% pivot_wider(id_cols = 'V',names_from="J",values_from="prop",values_fill = c(prop=0))
tmp_a = as.data.frame(tmp_a)
rownames(tmp_a) = tmp_a$V
tmp_a = tmp_a[,-1]
tmp_a = tmp_a[rowSums(tmp_a)>0.04,colSums(tmp_a)>0.03]
#library(pheatmap)

pheatmap(tmp_a,cluster_cols = F,cluster_rows = F,
         color = colorRampPalette((brewer.pal(n = 7, name ="Blues")))(100),
         filename = "fig4/TCR_VJ_HC.pdf",
         main="Health",
         width=6,height =3)


#heatmap_moderate

tmp_a= tmp_vj_comb[tmp_vj_comb$condition=="Moderate",] %>% filter(!(V=="None"|J=="None")) %>% pivot_wider(id_cols = 'V',names_from="J",values_from="prop",values_fill = c(prop=0))
tmp_a = as.data.frame(tmp_a)
rownames(tmp_a) = tmp_a$V
tmp_a = tmp_a[,-1]
tmp_a = tmp_a[rowSums(tmp_a)>0.04,colSums(tmp_a)>0.03]
#library(pheatmap)

pheatmap(tmp_a,cluster_cols = F,cluster_rows = F,
         color = colorRampPalette((brewer.pal(n = 7, name ="Blues")))(100),
         filename = "fig4/TCR_VJ_Mod.pdf",
         main="Moderate",
         width=6,height =3)


#heatmap_severe

tmp_a= tmp_vj_comb[tmp_vj_comb$condition=="Severe",] %>% filter(!(V=="None"|J=="None")) %>% pivot_wider(id_cols = 'V',names_from="J",values_from="prop",values_fill = c(prop=0))
tmp_a = as.data.frame(tmp_a)
rownames(tmp_a) = tmp_a$V
tmp_a = tmp_a[,-1]
tmp_a = tmp_a[rowSums(tmp_a)>0.04,colSums(tmp_a)>0.03]
#library(pheatmap)

pheatmap(tmp_a,cluster_cols = F,cluster_rows = F,
         color = colorRampPalette((brewer.pal(n = 7, name ="Blues")))(100),
         filename = "fig4/TCR_VJ_Sev.pdf",
         main="Severe",
         width=6,height =3)



## abundance and diversity plots

library(alakazam)


productive_clonotype_tbl <-
  seutnk@meta.data %>%
  filter(productive == "True") %>%
  select(barcode, raw_clonotype_id) %>%
  distinct() %>%
  inner_join(
    seutnk@meta.data[, c("Sample", "patient", "barcode","ids","condition","stage")],
    by = c("barcode" = "barcode"))


condition_curve <- estimateAbundance(
  productive_clonotype_tbl,
  group = "condition",
  ci = 0.95,
  nboot = 200,
  clone = "raw_clonotype_id")
p1 <- plot(condition_curve, colors=condition_color, silent = TRUE)

p1

#ggsave("vdj_tcells/test.pdf",width = 5,height = 5)
ggsave("vdjplots/condition_abundance.pdf",width = 4,height = 3)


condition_alpha_curve <- alphaDiversity(
  condition_curve,
  min_q = 0,
  max_q = 4,
  step_q = 0.1,
  ci = 0.95,
  nboot = 200)
p2 <- plot(condition_alpha_curve, colors=condition_color, silent = TRUE)

p2

ggsave("vdjplots/condition_diversity.pdf",width = 4,height = 3)

stage_curve <- estimateAbundance(
  productive_clonotype_tbl,
  group = "stage",
  ci = 0.95,
  nboot = 200,
  clone = "raw_clonotype_id")
p1 <- plot(stage_curve, colors=cond_palette, silent = TRUE)

p1

#ggsave("vdj_tcells/test.pdf",width = 5,height = 5)
ggsave("vdjplots/stage_abundance.pdf",width = 4,height = 3)


stage_alpha_curve <- alphaDiversity(
  stage_curve,
  min_q = 0,
  max_q = 4,
  step_q = 0.1,
  ci = 0.95,
  nboot = 200)
p2 <- plot(stage_alpha_curve, colors=cond_palette, silent = TRUE)

p2

ggsave("vdjplots/stage_diversity.pdf",width = 4,height = 3)

