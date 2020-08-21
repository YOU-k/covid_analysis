# scripts for fig2
seu2$condition[seu2$condition=="Health control"]="HC"
seu2@meta.data ->metadata
metadata %>% dplyr::filter(!patient %in% c("0/1","1/0")) ->metadata

metadata$stage = ">20"
metadata$stage[metadata$days<20] = "10-20"
metadata$stage[metadata$days<10] = "<10"
metadata$stage[metadata$condition=="Asymptomatic"]="Asymptomatic"
metadata$stage[metadata$condition=="HC"]="HC"
metadata$stage[metadata$condition=="Severe"]="Severe"
metadata$stage = factor(metadata$stage,levels = c("HC","Asymptomatic","Severe","<10","10-20",">20"))



library(tidyverse)
tmp1 = metadata %>% group_by(patient,days) %>% summarise(cell_n=n())
tmp2 = metadata %>% group_by(patient,days,cell_type) %>% summarise(cluster_n=n())
tmp2 = tmp2 %>% left_join(tmp1)
tmp2$cluster_n = tmp2$cluster_n/tmp2$cell_n
tmp2 = tmp2 %>% full_join(unique(metadata[,c("patient","condition","stage","sex","days")]),by=c("patient"="patient","days"="days"))

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
  facet_wrap(~cell_type,nrow=3,scales = "free_y")+
  #stat_compare_means(label = "p.signif",method="wilcox.test",comparisons=my_comparisons ,hide.ns=TRUE)+
  theme_cowplot(font_size = 15)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +labs(y="Proportion",fill="Stage")
ggsave("fig2/stage_all_bar.pdf",width = 12,height = 6)



tmp_sig <- tmp2[tmp2$cell_type %in% c("Activated CD4+ T","Memory CD8+ T","Proliferative T/NK","Platelet",
                                      "Plasma B","NK"),]
ggplot(data=tmp_sig,aes(x=stage,cluster_n,fill=stage))+
  geom_boxplot(outlier.size = 0.1)+
  geom_jitter(size=1,width=0.1)+
  scale_fill_manual(values = cond_palette)+
  facet_wrap(~cell_type,nrow=1,scales = "free_y")+
  stat_compare_means(label = "p.signif",method="wilcox.test",comparisons=my_comparisons ,hide.ns=TRUE)+
  theme_cowplot(font_size = 15)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +labs(y="Proportion",fill="Stage")
ggsave("fig2/stage_nop_bar2.pdf",width = 14,height = 2)


#draw a heatmap
library(pheatmap)
library(reshape2)
tmp2$id <- paste0(tmp2$patient,"_",tmp2$days)

dcast(tmp2 %>% select(id,cell_type, cluster_n),cell_type ~ id,mean) -> tmp2_heatmap
tmp2_heatmap[tmp2_heatmap=="NaN"] <-0 
rownames(tmp2_heatmap) <- tmp2_heatmap[,1]
tmp2_heatmap[,-1] -> tmp2_heatmap

##give colors

patient_col <- celltype_color[1:13]

names(patient_col) <- unique(tmp_col$patient)

#patient_col["J"] <-"yellowgreen"
#patient_col["K"] <-"orange"
#patient_col["B"] <-"#33CC66"

sex_col <- c(brewer.pal(2,"Paired"))[1:2]
names(sex_col) <- unique(tmp_col$sex)
tmp2$condition[tmp2$condition=="Health control"] <- "HC"

tmp2 %>% select(patient,days,sex,id,condition) %>% distinct() -> tmp_col

tmp2_heatmap[,tmp_col$id] -> tmp2_heatmap

apply(tmp2_heatmap, 1, scale) -> tmp2_scale

t(tmp2_scale) -> tmp2_scale
colnames(tmp2_scale) <- colnames(tmp2_heatmap)

tmp2_scale[tmp2_scale>3]=3
tmp2_scale[tmp2_scale<(-1.5)]=-1.5

p <- pheatmap(tmp2_scale,
              annotation_col = data.frame(
                Patient=tmp_col$patient,
                Sex= tmp_col$sex,
                Days= tmp_col$days,
                Condition=tmp_col$condition,
                row.names = colnames(tmp2_heatmap)
              ),
              annotation_colors  = list(Patient=patient_col,Sex=sex_col,Condition=condition_color),
              max.labels = Inf,
              normalize = TRUE,
              show.labels = FALSE,
              fontsize = 15,
              show_colnames = FALSE,treeheight_row = 0, treeheight_col = 0 ,border_color=NA ,clustering_method = )
pdf("fig2/heatmap_celltypes3.pdf",height=5,width=10)
print(p)
dev.off()

