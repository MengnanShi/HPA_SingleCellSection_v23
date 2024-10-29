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


# figure 3a
enriched_gene = readTXTfile(file.path(input.path,"proteinatlas.tsv"),0)
enriched_gene = enriched_gene %>% mutate(RNA.single.cell.type.specificity = case_when(RNA.single.cell.type.specificity == "" ~ "Not available", T~RNA.single.cell.type.specificity))


enriched_gene_uvial <- enriched_gene %>% 
  dplyr::select(Ensembl, RNA.single.cell.type.specificity, RNA.tissue.specificity) %>% 
  dplyr::rename(sc = RNA.single.cell.type.specificity) %>% 
  dplyr::rename(bulk = RNA.tissue.specificity) 

enriched_gene_uvial$sc = gsub("Cell type","Cell type/Tissue", enriched_gene_uvial$sc)
enriched_gene_uvial$sc = gsub("cell type","Cell type/Tissue", enriched_gene_uvial$sc)
enriched_gene_uvial$bulk = gsub("Tissue","Cell type/Tissue", enriched_gene_uvial$bulk)
enriched_gene_uvial$bulk = gsub("tissue","Cell type/Tissue", enriched_gene_uvial$bulk)


enriched_gene_uvial[is.na(enriched_gene_uvial)]="Not available"


enriched_gene_uvial_long <-  to_lodes_form(enriched_gene_uvial , key = "version", axes = c(2,3)) %>% 
  mutate(label = paste(stratum, version, sep="."))
unique(enriched_gene_uvial_long$stratum) 

enriched_gene_uvial_long$stratum <- 
  factor(enriched_gene_uvial_long$stratum, 
         levels = c("Cell type/Tissue enriched", "Group enriched", "Cell type/Tissue enhanced", "Low Cell type/Tissue specificity", "Not detected","Not available"))

enriched_gene_uvial_long <- enriched_gene_uvial_long %>% 
  left_join(enriched_gene_uvial_long %>% group_by(label) %>% summarise(n=n()), by="label")


dev.new()

ggplot(enriched_gene_uvial_long, aes(x = version, stratum = stratum, alluvium = alluvium, y = alluvium, fill = stratum, label=n)) +
  geom_lode() + 
  geom_flow(curve_type = "cubic", alpha = 2/3, width = 0.5) + 
  geom_stratum(alpha = 1, width = 0.5, color=NA) + 
  theme_minimal() + 
  scale_fill_manual(values = all_color) + 
  guides(fill=guide_legend(ncol=1)) +
  theme(legend.position = "right", panel.grid = element_blank(), axis.text.x = element_text(colour = "black", size = 10),
        axis.text.y = element_blank(), axis.title.y = element_blank())+
  geom_text(stat = "stratum", aes(label = n), size =5, color="white")

ggsave(file.path(plt.path,"uvial_plot_sc_bulk.pdf"))
dev.off()

# figure 3b
gene_class_sc =   
  enriched_gene %>% 
  mutate(Tissue= str_split(RNA.single.cell.type.specific.nTPM, ";")) %>% 
  unnest(cols = c(Tissue)) %>% 
  mutate(elevated_tissue = str_to_sentence(sapply(str_split(Tissue,":"),"[[",1))) %>% 
  dplyr::select(Ensembl,
                specificity_category = RNA.single.cell.type.specificity, 
                elevated_tissue) %>% 
  add_column(dataset = "singlecell") %>% 
  mutate(specificity_category = case_when(specificity_category== ""~"Not available", T~specificity_category))




gene_class_bulk =   
  enriched_gene %>% 
  mutate(Tissue= str_split(RNA.tissue.specific.nTPM, ";")) %>% 
  unnest(cols = c(Tissue)) %>% 
  mutate(elevated_tissue = str_to_sentence(sapply(str_split(Tissue,":"),"[[",1))) %>% 
  mutate(elevated_tissue = case_when(elevated_tissue == "Skin 1" ~ "Skin",
                                     elevated_tissue == "Stomach 1" ~ "Stomach",
                                     elevated_tissue == "Endometrium 1" ~ "Endometrium",
                                     TRUE ~ str_to_sentence(elevated_tissue))) %>% 
  dplyr::select(Ensembl,
                specificity_category = RNA.tissue.specificity, 
                elevated_tissue) %>% 
  add_column(dataset = "tissue")%>% 
  mutate(specificity_category = case_when(specificity_category=="" ~ "Not available", T~ specificity_category))


gene_class_long = as.data.frame(rbind(gene_class_sc, gene_class_bulk))

i="Liver"
temp_genes_bulk = gene_class_long %>% filter(dataset=="tissue" &
                                               specificity_category =="Tissue enriched"&
                                               elevated_tissue==i ) %>% pull(Ensembl)


temp_gene_sc_spec = gene_class_long %>% filter(dataset=="singlecell" & Ensembl %in% temp_genes_bulk)

# pie chart for the top cell type
top_celltype = temp_gene_sc_spec %>% group_by(elevated_tissue) %>% summarize(n=n()) %>% arrange(desc(n)) %>% pull(elevated_tissue)
top_celltype = top_celltype[1]
top_celltype_genes = temp_gene_sc_spec %>% filter(elevated_tissue==top_celltype) %>% distinct(Ensembl) %>% pull

top2_celltype = temp_gene_sc_spec %>% filter(!(Ensembl %in% top_celltype_genes) & !(specificity_category %in% c("Low cell type specificity","Not detected","Not available")) ) %>% 
  group_by(elevated_tissue) %>% summarize(n=n()) %>% arrange(desc(n)) %>% pull(elevated_tissue)
top2_celltype = top2_celltype[1]
top2_celltype_genes = temp_gene_sc_spec %>% 
  filter(!(Ensembl %in% top_celltype_genes) & elevated_tissue==top2_celltype) %>% 
  distinct(Ensembl) %>% pull


temp_color = c("#db7e8f","#00A087B2","#7999DC","grey40","grey57","grey")
names(temp_color) = c(top_celltype,top2_celltype,"Others","Low cell type specificity","Not detected","Not available")

# pie chart for category
p1 = temp_gene_sc_spec %>% distinct(Ensembl, .keep_all = T) %>% 
  group_by(specificity_category) %>% 
  summarize(n=n()) %>% 
  mutate(pct = n/sum(n)*100) %>% 
  arrange(desc(pct)) %>% 
  mutate(specificity_category = factor(specificity_category, levels = c("Cell type enriched","Group enriched", "Cell type enhanced","Low cell type specificity","Not detected","Not available"))) %>% 
  ggplot(aes(x="", y=pct, fill=specificity_category))+
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  theme_void() +
  theme(legend.position = "none",
        text = element_text(size=20))+
  scale_fill_manual(values=all_color)+
  geom_text_repel(aes(label=paste0(round(pct,2), "%")), 
                  position = position_stack(vjust=0.5),
                  size=7)+
  labs(title = paste0("Enriched genes of ",i,": ",length(unique(temp_gene_sc_spec$Ensembl))," genes"))


p2 = temp_gene_sc_spec %>% distinct(Ensembl, .keep_all = T) %>% 
  mutate(label = case_when(Ensembl %in% top_celltype_genes ~ top_celltype, 
                           Ensembl %in% top2_celltype_genes ~ top2_celltype, 
                           specificity_category %in% c("Low cell type specificity","Not detected","Not available") ~ specificity_category,
                           T~"Others")) %>% 
  group_by(label) %>% 
  summarize(n=n()) %>% 
  mutate(pct = n/sum(n)*100) %>% 
  select(label, pct) %>% 
  arrange(desc(pct)) %>% 
  mutate(label = factor(label, levels = c(top_celltype,top2_celltype,"Others","Low cell type specificity","Not detected","Not available"))) %>% 
  ggplot(aes(x="", y=pct, fill=label))+
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  theme_void() +
  theme(legend.position = "right",
        text = element_text(size=15),
        legend.text = element_text(size = 20))+
  # scale_fill_discrete(name="Genes elevated in")+
  scale_fill_manual(values=temp_color)+
  geom_text_repel(aes(label=paste0(round(pct,2), "%")), 
                  position = position_stack(vjust=0.5),
                  size=7)

p3 = ggarrange(p1, p2, ncol=1, nrow=2, legend="bottom")
ggsave(paste0(plt.path,"/Enriched genes of ",i,".pdf"), width=13, height=10)

# figure 3c
sc.cell_type = readTXTfile(file.path(input.path,"rna_single_cell_type.tsv"),0)
sc.cell_type = sc.cell_type %>% mutate(Cell.type=str_to_sentence(Cell.type))

i="FTCD"
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

# figure 3d
temp_gene = enriched_gene %>% 
  filter(RNA.tissue.specificity == "Low tissue specificity" & RNA.single.cell.type.specificity =="Low cell type specificity") 

low_specificity_gsoa_go = do_gsoa_go(temp_gene %>% pull(Ensembl), enriched_gene %>% pull(Ensembl) )

gotree_plot(low_specificity_gsoa_go, plot_title = "")
ggsave(file.path(plt.path,"tissue_cell_type_housekeeping_gotree.pdf"), width=14, height=8)

# figure 3e
bulk.ntpm = readTXTfile(file.path(input.path,"rna_tissue_consensus.tsv"),0)
bulk.ntpm = bulk.ntpm %>% mutate(Tissue = str_to_sentence(Tissue))

nodect_enriched_gene = enriched_gene %>% 
  filter(RNA.tissue.specificity != "Not detected" & RNA.single.cell.type.specificity %in% c("Not detected")) 

bulk.ntpm %>% filter(Gene %in% (nodect_enriched_gene %>% pull(Ensembl))) %>% 
  arrange(desc(nTPM)) %>% 
  distinct(Gene, .keep_all=T) %>% 
  mutate(label = case_when(Tissue %in% overlapped_tissues ~ "Overlapped tissues",
                           T ~ "Bulk specific tissues")) %>% 
  mutate(group="nTPM") %>% 
  ggplot(aes(x=group, y=log10(nTPM)))+
  geom_violin(fill="grey57", alpha=0.5)+
  geom_jitter(aes(color = label),width = 0.2, height = 0,  size = 1)+
  scale_color_manual(values=c("#E64B35B2","#4DBBD5B2" ))+
  xlab("")+
  ylab("Max Log10(nTPM)")+
  theme( panel.border = element_rect(fill=NA), 
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         axis.title.y = element_text(size = 10, colour = "black"), 
         panel.background = element_blank(), 
         # axis.text = element_text(size = 8, colour = "black"), 
         axis.text.x = element_blank(),
         axis.line = element_blank(), 
         axis.ticks.length.y = unit(.1, "cm"),
         axis.ticks.x = element_blank(),
         text = element_text(size=10),
         legend.text = element_text(size=10),
         legend.title = element_blank())
ggsave(file.path(plt.path,"tissue_cell_type_not_detected_sc_gene_violin.pdf"), height=3, width=5)

# figure 3f
sc.cluster = readTXTfile(file.path(input.path,"tissue_cell_type_exp_v23_exclude.tsv"),0)
sc.cluster = sc.cluster %>% 
  dplyr::rename(
    Gene = ensg_id,
    Tissue = tissue_name,
    Cell.type = cell_type_name,
    ave_nTPM = weighted_avg_exp) %>% 
  mutate(
    Cell.type = str_to_sentence(Cell.type),
    Tissue = case_when(Tissue=="pbmc" ~ "Pbmc",
                       TRUE ~ str_to_sentence(Tissue))
  )

nodect_enriched_gene = enriched_gene %>% 
  filter(RNA.tissue.specificity == "Not detected" & RNA.single.cell.type.specificity %in% c("Cell type enriched", "Cell type enhanced", "Low cell type specificity", "Group enriched")) 


temp = sc.cluster %>% 
  filter(Gene %in% (nodect_enriched_gene %>% pull(Ensembl))) %>% 
  group_by(Gene) %>% 
  summarize(max=max(ave_nTPM)) %>% 
  arrange(desc(max)) %>% 
  pull(Gene)

temp_color_pal = c("Detected in single"="#0A4DA1", "Detected in many"="#9FE3DB","Detected in some"="#6C8EEB","Detected in all"="#B0CCFF")

nodect_enriched_gene %>% 
  select(Ensembl,RNA.single.cell.type.distribution) %>% 
  group_by(RNA.single.cell.type.distribution) %>% summarize(n=n())%>% 
  mutate(RNA.single.cell.type.distribution = factor(RNA.single.cell.type.distribution, levels=c("Detected in single", "Detected in some","Detected in many","Detected in all"))) %>%
  ggplot(aes(x=RNA.single.cell.type.distribution, y=n, fill=RNA.single.cell.type.distribution)) +
  geom_bar(stat="identity", width=1)+
  xlab("Gene distribution in single cell")+ylab("Number of genes")+
  theme( panel.border = element_blank(), 
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         axis.title.y = element_text(size = 10, colour = "black"), 
         panel.background = element_blank(), 
         # axis.text = element_text(size = 8, colour = "black"), 
         axis.line = element_blank(), 
         axis.ticks.length = unit(.1, "cm"),
         axis.text.x = element_text(angle = 45, hjust = 1),
         # axis.title.x = element_blank(),
         legend.position = "right")+
  scale_fill_manual(values=temp_color_pal)
ggsave(file.path(plt.path,"tissue_cell_type_not_detected_bulk_gene_celltype_distribution.pdf"), height=5, width=5)




























