suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(psych)
  library(pheatmap)
  library(RColorBrewer)
  library(ggrepel)
  library(ggdendro)
  library(ggalluvial)
  library(gridExtra)
  library(ggpubr)
  library(cowplot)
  library(magrittr)
  library(circlize)
  library(xlsx)
  library(clusterProfiler)
  library(enrichplot)
  library(ggalt)
  library(tidygraph)
  library(ggraph)
  library(org.Hs.eg.db)
  library(ggsci)
  library(NOISeq)
  library(viridis)
  library(readxl)
  library(polycor)
  library(conflicted)
  library(tsne)
  library(aplot)
  library(dendextend)
})

source("functions.R")
source("theme.R")
# conflict_scout()
conflict_prefer("filter", "dplyr")
conflict_prefer("desc", "dplyr")
conflict_prefer("dotplot", "enrichplot")
conflict_prefer("select", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("spectrum", "igraph")
conflict_prefer("summarize", "dplyr")
conflict_prefer("pca","pcaMethods")
conflict_prefer("set_names","rlang")
setdiff = base::setdiff
melt = reshape2::melt
intersect = base::intersect
mutate = dplyr::mutate
arrange = dplyr::arrange

input.path="./data"
output.path="./results"
plt.path="./results/figures"


# figure 1a
cluster_des = readTXTfile(file.path(input.path,"rna_single_cell_cluster_description.tsv"),0)
cluster_des = cluster_des %>% 
  mutate(Cell.type=str_to_sentence(Cell.type)) %>% 
  mutate(Tissue=str_to_sentence(Tissue)) %>% 
  mutate(Cell.type.group=str_to_sentence(Cell.type.group)) %>% 
  filter(!(paste(Tissue, Cluster,sep="_") %in% exclude_cluster$label))


tissue_des.new = rbind( cbind(Tissue = c( "Colon", "Eye", "Heart muscle", "Kidney", "Pancreas", "Pbmc", "Placenta", "Rectum", "Skin", "Small intestine", "Testis"), version="v1"), 
                        cbind(Tissue = c("Adipose tissue" , "Bone marrow", "Brain", "Breast", "Bronchus", "Endometrium", "Esophagus", "Fallopian tube","Liver", "Lung", "Lymph node", "Ovary", "Prostate", "Salivary gland" , "Skeletal muscle", "Spleen", "Stomach", "Thymus", "Tongue", "Vascular"), version = "v2"))

cluster_num <- cluster_des %>% 
  filter(Cell.type.group!="Mixed cell types") %>% 
  group_by(Tissue, Cell.type.group) %>%
  summarise(Cell.count2 = sum(Cell.count)) %>% 
  group_by(Tissue) %>% 
  summarise(Cell.frac = 100*Cell.count2/sum(Cell.count2),
            Cell.type.group = Cell.type.group) %>% 
  dplyr::left_join(tissue_des.new %>% as.data.frame())


cluster_num %>% 
  ggplot(aes(y=Tissue, x=Cell.type.group))+
  geom_point(aes(size=as.numeric(Cell.frac)), color="firebrick4", alpha=0.7)+ 
  # scale_color_gradient(low="#FF9966",high="#CC3300")+
  theme_bw()+ 
  theme(plot.background = element_blank(),
        panel.border=element_rect(size=1, color="black"), 
        axis.text = element_text(size=13, face="bold"),
        axis.text.x = element_text(angle = 45, hjust =0),
        legend.position = "left")+ 
  scale_size(range = c(1,4))+ 
  xlab("")+ylab("")+ 
  guides(size=guide_legend(title="Percentage of cells (%)", color="none")) +
  ggplot2::facet_grid(rows = vars(version),space="free", scales = "free")+
  scale_x_discrete(expand = c(0.1, 0.1), position="top")
ggsave(file.path(plt.path,"tissue_cell_type_group_bubble_plot_pct.pdf"), height=8, width=10)


# figure 1b
sc.cluster2 = readTXTfile(file.path(input.path,"rna_single_cell_type_tissue.tsv"),0)
sc.cluster2 = sc.cluster2 %>% 
  mutate(Cell.type=str_to_sentence(Cell.type)) %>% 
  mutate(Tissue=str_to_sentence(Tissue))

exclude_cluster = readTXTfile(file.path(input.path,"excluded_clusters_v23.tsv"),0)
exclude_cluster = exclude_cluster %>% 
  mutate(tissue_name = str_to_sentence(tissue_name),
         cluster_id = paste0("c-",cluster_id),
         cell_type_name = str_to_sentence(cell_type_name),
         label = paste(tissue_name, cluster_id,sep="_"))

sc.cluster2 = sc.cluster2 %>% 
  filter(!(paste(Tissue, Cluster,sep="_") %in% exclude_cluster$label))

sc.cell_type.m = sc.cluster2 %>% 
  filter(Cell.type!="Mixed cell types") %>% 
  mutate(label = paste(Tissue, Cell.type,Cluster, sep="_")) %>% 
  select(-Tissue, -Cell.type, -Cluster, -Read.count, -Gene.name) %>% 
  pivot_wider(names_from = label, values_from = nTPM)

sc.cluster.des = cluster_des %>% 
  mutate(label = paste(Tissue, Cell.type,Cluster, sep="_")) %>% 
  distinct(label, .keep_all = T) %>% 
  mutate(cell.type.group.unique = replace(Cell.type.group, duplicated(Cell.type.group), NA))

sc.cell_type.m.pca_data = sc.cell_type.m %>% 
  # select(-Gene.name) %>% 
  column_to_rownames("Gene") %>% 
  {log10(.+1)} %>% 
  t() 

sc.cell_type.m.pca = pca_calc(sc.cell_type.m.pca_data, npcs=100)
print(sc.cell_type.m.pca$scores[,"PC1"] %>% {c(min(.), max(.))})
print(sc.cell_type.m.pca$scores[,"PC2"] %>% {c(min(.), max(.))})


tissue.pal = tissue.color
set.seed(123)
sc.cell_type.m.tsne <- 
  sc.cell_type.m.pca$scores%>% 
  tsne(initial_dims=35)
rownames(sc.cell_type.m.tsne) = rownames(sc.cell_type.m.pca$scores)

sc.cell_type.m.tsne %>% 
  as.tibble(rownames="label") %>% 
  left_join(sc.cluster.des) %>% 
  ggplot(aes(V1, V2, color=Cell.type.group)) +
  scale_color_manual(values=all_color)+
  geom_point(size=1.5, alpha=1)+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        text = element_text(size=25),
        axis.title = element_text(size=15),
        axis.text = element_text(size=13))+
  geom_text_repel(aes(label=cell.type.group.unique, color=Cell.type.group),  max.overlaps = 20)+
  xlab("tSNE 1")+ylab("tSNE 2")

ggsave(file.path(plt.path,"tSNE plot cell types.pdf"), width=5, height=5)


