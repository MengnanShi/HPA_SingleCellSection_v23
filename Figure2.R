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

# figure 2a
enriched_gene = readTXTfile(file.path(input.path,"proteinatlas.tsv"),0)
enriched_gene = enriched_gene %>% mutate(RNA.single.cell.type.specificity = case_when(RNA.single.cell.type.specificity == "" ~ "Not available", T~RNA.single.cell.type.specificity))

enriched_gene.v20 = readTXTfile(file.path(input.path,"proteinatlas_v20.tsv"),0)
enriched_gene.num = readTXTfile(file.path(input.path,"enriched_gene_num.txt"),0)
enriched_gene.num = na.omit(enriched_gene.num)
enriched_gene.num.v20 = readTXTfile(file.path(input.path,"enriched_gene_num.v20.txt"),0)

enriched_gene_uvial <- enriched_gene %>% 
  select(Ensembl, RNA.single.cell.type.specificity) %>% 
  left_join(enriched_gene.v20 %>% select(Ensembl, RNA.single.cell.type.specificity), by="Ensembl") %>% 
  dplyr::rename(v2_single_cell = RNA.single.cell.type.specificity.x) %>% 
  dplyr::rename(v1_single_cell = RNA.single.cell.type.specificity.y)

enriched_gene_uvial[is.na(enriched_gene_uvial)]="Not available"
enriched_gene_uvial[enriched_gene_uvial==""]="Not available"


enriched_gene_uvial_long <-  to_lodes_form(enriched_gene_uvial , key = "version", axes = c(2,3)) %>% 
  mutate(label = paste(stratum, version, sep="."))
unique(enriched_gene_uvial_long$stratum) 

enriched_gene_uvial_long$stratum <- 
  factor(enriched_gene_uvial_long$stratum, 
         levels = c("Cell type enriched", "Group enriched", "Cell type enhanced", "Low cell type specificity", "Not detected", "Not available"))

enriched_gene_uvial_long <- enriched_gene_uvial_long %>% 
  left_join(enriched_gene_uvial_long %>% group_by(label) %>% summarise(n=n()), by="label")


enriched_gene_uvial_long = enriched_gene_uvial_long %>% 
  mutate(version = case_when(version =="v2_single_cell" ~"SC_V2", T~"SC_V1")) %>% 
  mutate(version = factor(version, levels= c("SC_V1","SC_V2"))) %>% 
  rename("Stratum" = "stratum")

# dev.new()

ggplot(enriched_gene_uvial_long, aes(x = version, stratum = Stratum, alluvium = alluvium, y = alluvium, fill = Stratum, label=n)) +
  geom_lode() + 
  geom_flow(curve_type = "cubic", alpha = 2/3, width = 0.5) + 
  geom_stratum(alpha = 1, width = 0.5, color=NA) + 
  theme_minimal() + 
  scale_fill_manual(values = all_color) + 
  guides(fill=guide_legend(ncol=1)) +
  theme(legend.position = "bottom", 
        panel.grid = element_blank(), 
        axis.text.x = element_text(colour = "black", size = 15),
        axis.title.x = element_blank(), 
        axis.text.y = element_blank(), 
        axis.title.y = element_blank(),
        legend.text = element_text(color="black",size=15))+
  geom_text(stat = "Stratum", aes(label = n), size =7, color="white")

ggsave(file.path(plt.path,"uvial_plot_v20_v23.pdf"))

# figure 2b
celltype_order = enriched_gene.num %>%  
  arrange(desc(Total.elevated)) %>%
  pull(Cell.type) %>% str_to_sentence()

right_plot = enriched_gene.num %>% 
  select(-Cell.type.group, -Total.elevated) %>% 
  pivot_longer(!Cell.type, names_to = "Specificity", values_to = "value") %>% 
  mutate(Specificity = case_when(Specificity=="Enriched" ~ "Cell type enriched",
                                 Specificity=="Enhanced" ~ "Cell type enhanced",
                                 Specificity=="Group.enriched" ~ "Group enriched")) %>% 
  mutate(Cell.type = str_to_sentence(Cell.type),
         Specificity = factor(Specificity, levels=c("Cell type enriched","Group enriched", "Cell type enhanced")),
         Cell.type = factor(Cell.type, celltype_order))

p1 = right_plot %>%
  ggplot(aes(fill=Specificity, y=value, x=Cell.type)) + 
  geom_bar(stat="identity") + 
  theme_bw() + 
  scale_fill_manual(values = all_color) +
  guides(fill = guide_legend(override.aes = list(size = 5))) +
  scale_x_discrete(expand = c(0,0), position = "bottom")+
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 10, colour = "black"), 
        axis.title.x = element_blank(), 
        panel.background = element_blank(), 
        axis.text = element_text(size = 13, colour = "black"), 
        axis.text.x = element_blank(),
        # axis.text.x = element_text(angle = 90, hjust = 1, vjust=0),
        axis.line = element_blank(), 
        axis.ticks = element_blank(), 
        legend.position = "bottom", 
        legend.title = element_text(size = 11, colour = "black"),
        legend.text = element_text(size = 10, colour = "black")) + 
  ylab("Number of genes")

ggsave(file.path(plt.path,"elevated genes right plot bar plot.pdf"), width=15,height=3)

p2 = right_plot %>%
  ggplot(aes(fill=Cell.type, y=1, x=Cell.type)) + 
  geom_bar(stat="identity") + 
  theme_bw() + 
  scale_fill_manual(values = all_color) +
  guides(fill = guide_legend(override.aes = list(size = 5))) +
  scale_x_discrete(expand = c(0,0), position = "bottom")+
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(), 
        panel.background = element_blank(), 
        axis.text = element_text(size = 13, colour = "black"), 
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.5),
        axis.line = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.length = unit(.1, "cm"),
        legend.position = "none") 

ggsave(file.path(plt.path,"elevated genes right plot bar plot legend.pdf"), width=15, height=2.95)

p4 = right_plot %>%
  ggplot(aes(fill=Cell.type, y=1, x=Cell.type)) + 
  geom_bar(stat="identity") + 
  theme_bw() + 
  scale_fill_manual(values = all_color) +
  guides(fill = guide_legend(override.aes = list(size = 5))) +
  scale_x_discrete(expand = c(0,0), position = "bottom")+
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(), 
        panel.background = element_blank(), 
        axis.text = element_text(size = 13, colour = "black"), 
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.line = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.length = unit(.1, "cm"),
        legend.position = "none") 

celltype_order = enriched_gene.num %>%  
  arrange(desc(Total.elevated)) %>%
  pull(Cell.type) %>% str_to_sentence()

gene_anno_temp =   enriched_gene %>% 
  filter(RNA.single.cell.type.specificity != "Not detected" & RNA.single.cell.type.specificity != "Low cell type specificity" & RNA.single.cell.type.specificity != "Not available") %>% 
  mutate(Cell.type= str_split(RNA.single.cell.type.specific.nTPM, ";")) %>% 
  unnest() %>% 
  mutate( Cell.type = sapply(strsplit(Cell.type, ":"),"[[",1),
          Cell.type = str_to_sentence(Cell.type)) %>% 
  dplyr::select(gene=Ensembl,RNA.single.cell.type.specificity, Cell.type) 


enriched_gene_count = c()
for(celltype in unique(sc.cell_type$Cell.type)){
  temp = sc.cell_type %>% filter(Cell.type == celltype)
  total_count = sum(temp$nTPM)
  Enriched_count = temp %>% 
    filter(Gene %in% (gene_anno_temp %>% filter(Cell.type==celltype, RNA.single.cell.type.specificity=="Cell type enriched") %>% pull(gene))) %>% pull(nTPM) %>% sum()
  
  Group.enriched_count = temp %>% 
    filter(Gene %in% (gene_anno_temp %>% filter(Cell.type==celltype, RNA.single.cell.type.specificity=="Group enriched") %>% pull(gene)))%>% 
    pull(nTPM) %>% sum()
  Enhanced_count = temp %>% 
    filter(Gene %in% (gene_anno_temp %>% filter(Cell.type==celltype, RNA.single.cell.type.specificity=="Cell type enhanced") %>% pull(gene)))%>%    pull(nTPM) %>% sum()
  
  x=data.frame(Cell.type=celltype, 
               Enriched=as.numeric(Enriched_count), 
               Group.enriched=as.numeric(Group.enriched_count),
               Enhanced=as.numeric(Enhanced_count),
               total_count = as.numeric(total_count))
  enriched_gene_count = rbind(enriched_gene_count, x)
}

enriched_gene_count = enriched_gene_count %>% 
  mutate(Enriched.pct = Enriched/total_count*100,
         Group.enriched.pct = Group.enriched/total_count*100,
         Enhanced.pct = Enhanced/total_count*100)
writeTXTfile(enriched_gene_count, file.path(output.path,"gene_count_partial_celltype_elevated.txt"))


left_plot = enriched_gene_count %>% 
  select(Cell.type, Enriched.pct,Group.enriched.pct, Enhanced.pct) %>% 
  pivot_longer(!Cell.type, names_to = "Specificity", values_to = "value") %>% 
  mutate(Specificity = case_when(Specificity=="Enriched.pct" ~ "Cell type enriched",
                                 Specificity=="Enhanced.pct" ~ "Cell type enhanced",
                                 Specificity=="Group.enriched.pct" ~ "Group enriched")) %>% 
  mutate(Specificity = factor(Specificity, levels=c("Cell type enriched","Group enriched", "Cell type enhanced")),
         Cell.type = factor(Cell.type, celltype_order))

p3 = left_plot %>%
  ggplot(aes(fill=Specificity, y=value, x=Cell.type)) + 
  geom_bar(stat="identity") + 
  theme_bw() + 
  scale_fill_manual(values = all_color) +
  guides(fill = guide_legend(override.aes = list(size = 5))) +
  scale_x_discrete(expand = c(0,0), position = "bottom")+
  scale_y_reverse()+
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 10, colour = "black"), 
        axis.title.x = element_blank(), 
        panel.background = element_blank(), 
        axis.text = element_text(size = 13, colour = "black"), 
        axis.text.x = element_blank(),
        # axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),
        axis.line = element_blank(), 
        axis.ticks = element_blank(), 
        legend.position = "none", 
        legend.title = element_text(size = 11, colour = "black"),
        legend.text = element_text(size = 10, colour = "black")) + 
  ylab("Percentage of elevated gene counts")

p3%>% insert_top(p4, height=0.05)  %>% insert_top(p2, height=0.05) %>% insert_top(p1)

ggsave(file.path(plt.path,"elevated genes plot merge.pdf"), width=20, height=10)

# figure 2c
sc.cell_type = readTXTfile(file.path(input.path,"rna_single_cell_type.tsv"),0)
sc.cell_type = sc.cell_type %>% mutate(Cell.type=str_to_sentence(Cell.type))

i="CLCA1"
sc.cell_type %>% filter(Gene.name==i) %>% 
  select(Cell.type, nTPM) %>% 
  mutate(Cell.type = factor(Cell.type,levels=celltype_order_hpa)) %>% 
  ggplot() +
  geom_bar(aes(x=Cell.type, y=nTPM, fill=Cell.type),stat="identity")+
  scale_fill_manual(values=all_color)+
  theme(panel.border = element_rect(color="#78A493", fill=NA,linewidth = 2), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(), 
        panel.background = element_blank(), 
        axis.text = element_blank(),
        axis.text.x = element_blank(),
        # axis.text.x = element_text(angle = 90, hjust = 1, vjust=0),
        axis.line = element_blank(), 
        axis.ticks = element_blank(), 
        legend.position = "none")

ggsave(file.path(plt.path,paste0("sc_nTPM_bar_",i,".pdf")), width=1,height=1)

# figure 2f
enriched_gene_uvial <- enriched_gene %>% 
  dplyr::select(Ensembl, RNA.single.cell.type.specificity) %>% 
  left_join(enriched_gene.v20 %>% dplyr::select(Ensembl, RNA.single.cell.type.specificity), by="Ensembl") %>% 
  dplyr::rename(v2_single_cell = RNA.single.cell.type.specificity.x) %>% 
  dplyr::rename(v1_single_cell = RNA.single.cell.type.specificity.y)

enriched_gene_uvial[is.na(enriched_gene_uvial)]="Not detected"

enriched_gene_sub <- enriched_gene %>% dplyr::select(Gene, Ensembl, RNA.single.cell.type.specific.nTPM, RNA.single.cell.type.specificity.score) 

enriched_gene_sel = enriched_gene_uvial %>% 
  filter(v2_single_cell=="Cell type enriched" & v1_single_cell!="Cell type enriched") %>% 
  left_join(enriched_gene_sub, by="Ensembl") %>% 
  group_by(v1_single_cell) %>% 
  slice_max(order_by = RNA.single.cell.type.specificity.score, n = 5)


cell_group.pal = cell_type.color %>% 
  dplyr::select(-Cell.type) %>% 
  add_row(color="red",Cell.type.group="v2" ) %>% 
  add_row(color="green", Cell.type.group="v1") %$%
  set_names(color, Cell.type.group)

cell_type.v23 = unique(str_to_sentence(enriched_gene.num$Cell.type))
cell_type.v20 = unique(str_to_sentence(enriched_gene.num.v20$Cell.type))
cell_type.diff = setdiff(cell_type.v23, cell_type.v20)

cell_type_anno_temp = cell_type_anno %>% 
  mutate(version = case_when(Cell.type %in% cell_type.diff ~ "v2", T ~ "v1")) %>% 
  dplyr::select(Cell.type, Cell.type.group, version) 

cell_type_order = cell_type_anno_temp[order(cell_type_anno_temp[,3], cell_type_anno_temp[,2], decreasing=T),"Cell.type"]

cell_type_anno_temp_long = cell_type_anno_temp %>% 
  pivot_longer(!Cell.type,names_to = "anno", values_to = "label")
cell_type_anno_temp_long$anno = factor(cell_type_anno_temp_long$anno, 
                                       levels=c("Cell.type.group","version"))

cell_type_anno_temp_long$Cell.type = factor(cell_type_anno_temp_long$Cell.type, 
                                            levels=cell_type_order)

temp = sc.cell_type
temp$Cell.type = factor(temp$Cell.type, levels=cell_type_order)

for (i in c("IRF7", "NRAP","GPIHBP1")){
  temp %>% 
    filter(Gene.name==i) %>% 
    dplyr::select(Cell.type, nTPM) %>% arrange(desc(nTPM)) %>% 
    ggplot(aes(x=Cell.type, y=nTPM, fill=Cell.type))+
    geom_bar(stat="identity")+
    theme_bw()+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.position = "none",
          # axis.text.x = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1),
          axis.title.x = element_blank(),
          plot.title=element_text(hjust=0.5))+
    scale_fill_manual(values=all_color)+
    ggtitle(i)
  ggsave(file.path(plt.path,paste0("enriched_",i,".pdf")), height=3, width=14, limitsize = F)
  
}


legend_bar = cell_type_anno_temp_long %>% 
  mutate(anno1 = 1) %>% 
  # mutate(anno1 = ifelse(annot == "version", 1, ))
  ggplot(aes(x = Cell.type, y = anno1, fill=label))+ 
  geom_bar(stat='identity')+
  scale_fill_manual(values=cell_group.pal)+  
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        # axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        legend.position = "bottom")
ggsave(file.path(plt.path,"enriched_glegend.pdf"), height=3.6, width=13)































