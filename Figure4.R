
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
})

source("functions.R")
source("theme.R")
source("Deconvolution_functions.R")

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

# run Deconvolution_ImmuneCellType.R and Deconvolution_BulkTissue.R to get deconvolution results

pct_cutoff=5
method="DWLS"

dwls_res1 = readTXTfile(file.path(dev_output.path,"results_hpa",paste0("allCounts_",method,".txt")),1)
colnames(dwls_res1) = str_replace_all(colnames(dwls_res1),"_"," ")

dwls_res2 = readTXTfile(file.path(dev_output.path,"results_hpa",paste0("blood_allCounts_",method,".txt")),1)[,-19]
blood_celltype_temp = data.frame(Cell.type=names(blood.color)) %>% 
  mutate(cell.type.syn2 = str_replace_all(Cell.type, " ","_"),
         cell.type.syn2 = str_replace_all(cell.type.syn2, "-","_")) %>% 
  column_to_rownames("cell.type.syn2")
colnames(dwls_res2) = blood_celltype_temp[str_to_sentence(str_replace_all(colnames(dwls_res2),"\\.","_")),"Cell.type"]

dwls_res = merge(dwls_res1, dwls_res2, by="row.names")
dwls_res = dwls_res %>% column_to_rownames("Row.names")

cell_type_anno_temp = cell_type_anno %>% 
  mutate(cell.type.syn2 = str_replace_all(Cell.type, " ","\\."),
         cell.type.syn2 = str_replace_all(cell.type.syn2, "-","\\.")) %>% 
  column_to_rownames("cell.type.syn2")

rownames(dwls_res) = cell_type_anno_temp[rownames(dwls_res),"Cell.type"]


plot_order = dwls_res %>% 
  cluster_wide_data(distance_method = "euclidean",
                    clustering_method = "average", 
                    cluster_rows = T,
                    cluster_cols = T) %>% 
  {list(celltype = rownames(.),
        tissue = colnames(.))}


plot_data_ordered = dwls_res %>% 
  as.tibble(rownames="celltype") %>% 
  pivot_longer(!celltype, names_to = "tissue", values_to = "pct") %>% 
  mutate(pct = pct*100) %>% 
  mutate(pct = case_when(pct<10e-7 ~ 0,  T~ pct),
         dataset = case_when(tissue %in% colnames(dwls_res1) ~ "Tissue", 
                             tissue %in% colnames(dwls_res2) ~ "Immune cells",
                             T ~"")) %>% 
  mutate(dataset = factor(dataset, levels=c("Tissue","Immune cells"))) %>% 
  filter(pct>pct_cutoff) %>% 
  group_by(tissue) %>% 
  top_n(5, pct) %>% 
  arrange(tissue, desc(pct)) %>% 
  mutate(tissue = factor(tissue, levels=unique(names(sort(all_color2)))))

a  = plot_data_ordered %>% 
  left_join(data.frame(tissue =names(all_color2),color=all_color2), by="tissue") %>% 
  mutate(color=factor(color, levels=c("#1280C4",setdiff(unique(color),"#1280C4")))) %>% 
  arrange(dataset, color, tissue, desc(pct))

immune_order=c("Plasma cells","B-cells","Nk-cells","T-cells","Monocytes","Dendritic cells","Granulocytes","Erythroid cells","Hofbauer cells","Kupffer cells","Langerhans cells","Macrophages","Mixed immune cells", "Neutrophils")
immune_tissue_order  =c("Thymus", "Lymph node","Spleen", "Bone marrow","Memory b-cell","Naive b-cell","Gdt-cell","Nk-cell","Mait t-cell","Memory cd4 t-cell","Memory cd8 t-cell","Naive cd4 t-cell","Naive cd8 t-cell","T-reg")
a = a %>%
  mutate(dataset=factor(dataset, levels=unique(a$dataset)),
         tissue = factor(tissue, levels=c(immune_tissue_order,setdiff(unique(a$tissue),immune_tissue_order))),
         celltype = factor(celltype, levels=c(immune_order,setdiff(unique(a$celltype),immune_order))))



a %>% 
  ggplot(aes(y=celltype, x=tissue, size=pct, color=pct))+
  geom_point() +
  geom_tile(data = . %>% 
              select(celltype) %>% 
              distinct(),
            aes(0, celltype, fill = celltype),
            show.legend = F,
            inherit.aes = F,
            color = "black")  +
  geom_tile(data = . %>% 
              select(tissue, dataset) %>% 
              distinct(),
            aes(tissue, 0, fill = tissue),
            show.legend = F,
            inherit.aes = F,
            color = "black")+
  facet_grid(~dataset,
             space = "free",
             scales = "free_x") +
  scale_color_gradient(high = "orangered", 
                       low = "lightgray", guide="legend") +
  scale_fill_manual(values = all_color2) +
  scale_size_continuous(range = c(1,5)) +
  labs(y = "Single cell type") +
  # coord_fixed() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, 
                                   vjust = 0.5),
        text = element_text(size=20))
ggsave(file.path(dev_output.path,"results_hpa","aDWLS_res_bubble_plot_all_ordered.pdf"), width=15, height=15)

length(unique(a$tissue))


