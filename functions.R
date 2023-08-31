
writeTXTfile <- function(variable, filename){
  write.table(cbind(rownames(variable),variable), file=paste0(filename), row.names=FALSE, quote=FALSE, sep="\t")
}

readTXTfile = function(filename,rownames){
  if(rownames==0){read.csv(file=filename, sep="\t", header=T, stringsAsFactors = F)}else{
    read.csv(file=filename, sep="\t", header=T, stringsAsFactors = F, row.names=rownames)
  }}


writeXLSXfile_loop <- function(variable, sheetname, filename){
  if(file.exists(filename)){
    exist=T
  }else{exist=F}
  write.xlsx(cbind(rownames(variable),variable), file= paste0(filename),
             sheetName = sheetname,
             append=exist,row.names=F,col.names=T)
}

makeMatrix <- function( rows,cols, value=0){
  data = matrix(value, ncol = length(cols), nrow=length(rows))
  colnames(data) = cols
  rownames(data) = rows
  return(data)
}

simple_theme <- 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

tmm_normalize <- function(x, ...) {
  median_column <-
    apply(x,
          MARGIN = 1,
          median)
  x[, dim(x)[2] + 1] <- median_column
  norm_data <- NOISeq::tmm(x, refColumn = dim(x)[2], ...)
  norm_data[, -(dim(x)[2])] 
}


pca_calc <- function(data, npcs) {
  require(pcaMethods)
  pca_res <-
    data %>%
    pca(nPcs = npcs)
  
  pca_stats <-
    tibble(PC = 1:npcs,
           R2cum = R2cum(pca_res))
  
  informative_pcs <- pca_stats$PC[which(pca_stats$R2cum > 0.95)[1]]
  
  pca_stats_plot <-
    pca_stats %>%
    dplyr::select(PC, R2cum) %>%
    ggplot(aes(PC, R2cum)) +
    geom_point() +
    geom_line() +
    simple_theme +
    geom_vline(xintercept = informative_pcs, linetype = "dashed") +
    annotate("text",
             x = informative_pcs,
             y = 0.55,
             label = paste0("PC ", informative_pcs,
                            "\nR2 = ", round(pca_stats[informative_pcs, ]$R2cum, 3)),
             hjust = 1,
             vjust = 0)
  
  list(pca = pca_res,
       scores = pcaMethods::scores(pca_res),
       loadings = pcaMethods::loadings(pca_res),
       stats = pca_stats,
       pc_95 = informative_pcs,
       plot = pca_stats_plot)
}

umap_calc <- function(data, npcs = 2, n_epochs = 400, n_neighbors = round(sqrt(dim(data)[1]))) {
  require(uwot)
  data %>%
    umap(n_components = npcs,
         n_epochs = n_epochs,
         n_neighbors = n_neighbors)
  
}


hpa_gene_classification <- 
  function(data, expression_col, tissue_col, gene_col, enr_fold, max_group_n, det_lim = 1) {
    data_ <- 
      data %>% 
      dplyr::select(gene = gene_col,
                    expression = expression_col,
                    tissue = tissue_col) %>% 
      mutate(expression = round(expression, 4)) 
    
    if(any(is.na(data_$expression))) stop("NAs in expression column")
    if(any(is.na(data_$gene))) stop("NAs in gene column")
    if(any(is.na(data_$tissue))) stop("NAs in tissue column")
    
    n_groups <- length(unique(data_$tissue))
    
    gene_class_info <- 
      data_ %>%
      group_by(gene) %>%
      summarise(
        
        # Gene expression distribution metrics
        mean_exp = mean(expression, na.rm = T),
        min_exp = min(expression, na.rm = T),
        max_exp = max(expression, na.rm = T), 
        max_2nd = sort(expression)[length(expression)-1],
        
        # Expression frequency metrics
        n_exp = length(which(expression >= det_lim)),
        frac_exp = n_exp/length(expression[!is.na(expression)])*100,
        
        # Limit of enhancement metrics
        lim = max_exp/enr_fold, 
        
        exps_over_lim = list(expression[which(expression >= lim & expression >= det_lim)]),
        n_over = length(exps_over_lim[[1]]), 
        mean_over = mean(exps_over_lim[[1]]),
        min_over = ifelse(n_over == 0, NA,
                          min(exps_over_lim[[1]])),
        
        max_under_lim = max(expression[which(expression < min_over)], det_lim*0.1),
        
        
        exps_enhanced = list(which(expression/mean_exp >= enr_fold & expression >= det_lim)),
        
        
        
        
        # Expression patterns
        enrichment_group = paste(sort(tissue[which(expression >= lim & expression >= det_lim)]), collapse=";"),
        
        n_enriched = length(tissue[which(expression >= lim & expression >= det_lim)]),
        n_enhanced = length(exps_enhanced[[1]]), 
        enhanced_in = paste(sort(tissue[exps_enhanced[[1]]]), collapse=";"),
        n_na = n_groups - length(expression),
        max_2nd_or_lim = max(max_2nd, det_lim*0.1),
        tissues_not_detected = paste(sort(tissue[which(expression < det_lim)]), collapse=";"),
        tissues_detected = paste(sort(tissue[which(expression >= det_lim)]), collapse=";")) 
    
    
    gene_categories <- 
      gene_class_info %>%
      
      mutate(
        spec_category = case_when(n_exp == 0 ~ "not detected", 
                                  
                                  # Genes with expression fold times more than anything else are tissue enriched
                                  max_exp/max_2nd_or_lim >= enr_fold ~ "tissue enriched", 
                                  
                                  # Genes with expression fold times more than other tissues in groups of max group_n - 1 are group enriched
                                  max_exp >= lim &
                                    n_over <= max_group_n & n_over > 1 &
                                    mean_over/max_under_lim >= enr_fold ~ "group enriched", 
                                  
                                  # Genes with expression in tissues fold times more than the mean are tissue enhance
                                  n_enhanced > 0 ~ "tissue enhanced", 
                                  
                                  # Genes expressed with low tissue specificity
                                  T ~ "low tissue specificity"), 
        
        
        dist_category = case_when(frac_exp == 100 ~ "detected in all",
                                  frac_exp >= 33 ~ "detected in many",
                                  n_exp > 1 ~ "detected in some",
                                  n_exp == 1 ~ "detected in single",
                                  n_exp == 0 ~ "not detected"),
        
        spec_score = case_when(spec_category == "tissue enriched" ~ max_exp/max_2nd_or_lim,
                               spec_category == "group enriched" ~ mean_over/max_under_lim, 
                               spec_category == "tissue enhanced" ~ max_exp/mean_exp)) 
    
    
    
    
    ##### Rename and format
    gene_categories %>%
      mutate(enriched_tissues = case_when(spec_category %in% c("tissue enriched", "group enriched") ~ enrichment_group,
                                          spec_category == "tissue enhanced" ~ enhanced_in),
             n_enriched = case_when(spec_category %in% c("tissue enriched", "group enriched") ~ n_enriched,
                                    spec_category == "tissue enhanced" ~ n_enhanced)) %>%
      dplyr::select(gene, 
                    spec_category, 
                    dist_category, 
                    spec_score,
                    n_expressed = n_exp, 
                    fraction_expressed = frac_exp,
                    max_exp = max_exp,
                    enriched_tissues,
                    n_enriched,
                    n_na = n_na,
                    tissues_not_detected,
                    tissues_detected) 
    
    
    
  }


do_phyper_test <- function(geneList1, geneList2, interGene=NA, backgroudGeneLength){
  if (is.na(interGene)){interGene= intersect(geneList1,geneList2)}
  if(class(interGene)=="matrix"){q = nrow(interGene)}else if(class(interGene)=="character"){q = length(interGene)}
  
  m = length(geneList1)
  n = length(geneList2)
  if(q==0){p=1}else{
    p = phyper(q-1, m, backgroudGeneLength-m, n , lower.tail = FALSE)}
  
  x = t(as.matrix(c(q,m,n,p,backgroudGeneLength)))
  colnames(x) = c("interGene", "list1", "list2","phyper","backgroundGene")
  return(x)
  
}

pheatmap.hygeo<- function(mat,path,file="pheatmap.pdf", anno=TRUE, save=TRUE){
  if(anno==TRUE){
    mat.anno = mat
    mat.anno[mat.anno<0.05]="*"
    mat.anno[mat.anno>0.05]=""
  }else{mat.anno=FALSE}
  
  my.break=c(seq(0,1,by=0.001))
  my.color = c(colorRampPalette(colors = c("tomato","white"))(length(my.break)))
  # pdf(paste0(path,"coexp_fuzz_phyper_heatmap.pdf"))
  p=pheatmap(mat,
             color=my.color,
             show_rownames = TRUE,
             cluster_cols = TRUE,
             cluster_rows = TRUE,
             display_numbers = mat.anno,
             number_color="gray0",
             breaks = my.break,
             annotation_names_row = FALSE, 
             annotation_names_col = FALSE,
             cellwidth = 15,
             cellheight = 15,
             border_color=NA,
             treeheight_col=7,
             treeheight_row = 0,
             fontsize = 8.5)
  if(save==TRUE){
    if (file.exists(paste0(path,file_name))){print("File exists. Be aware of overwriting.")}
    ggsave(paste0(path,file_name),plot=p,limitsize = FALSE)
    dev.off()
  }else{
    return(p)
  }
  # dev.off()
}


calculate_tau_score <- 
  function(wide_data) {
    max_exp <- 
      apply(wide_data,
            MARGIN = 1,
            function(x) max(x, na.rm = T))
    
    N <- 
      apply(wide_data,
            MARGIN = 1,
            function(x) length(which(!is.na(x))))
    
    expression_sum <- 
      wide_data %>% 
      sweep(MARGIN = 1, 
            STATS = max_exp, 
            FUN = `/`) %>% 
      {1 - .} %>% 
      apply(MARGIN = 1,
            function(x) sum(x, na.rm = T))
    
    
    tau_score <- 
      (expression_sum / (N - 1)) %>% 
      enframe("gene", "tau_score")
    
    tau_score
  }

plot_network <- function(all_color, node, corNet, fc, nodesize=1, labelsize=1){
  set.seed(123)
  node = node %>% left_join(all_color) %>% 
    mutate(color=case_when(nodename %in% (tissue_des %>% filter(Organ=="Brain") %>% select(Tissue) %>% pull) ~ "#FFDD00",T~color)) %>% 
    column_to_rownames("nodename")
  node = node[V(corNet),]
  
  V(corNet)$color= node %>%select(color) %>% pull
  V(corNet)$shape= node %>% mutate(shape=case_when(label=="Tissue"~"square", 
                                                   label=="Cell.type"~"circle")) %>% 
    select(shape) %>% pull
  
  # add_shape("nil")
  if(labelsize!=0){
    plot(fc, corNet, 
         vertex.size=nodesize, 
         layout=layout_with_graphopt, 
         vertex.label.dist = labelsize,
         cex=3)
  }else{
    plot(fc, corNet, 
         vertex.size=nodesize, 
         layout=layout_with_graphopt, 
         vertex.label=NA,
         cex=3)
  }
  
}



do_gsoa_go <- function(genelist, bggenes, keyType="ENSEMBL"){
  gsoa_results <- 
    clusterProfiler::enrichGO(gene = genelist, 
                              OrgDb = org.Hs.eg.db, 
                              keyType = keyType, 
                              universe = bggenes,
                              ont = "BP", 
                              pAdjustMethod = "fdr", 
                              readable = TRUE)
  if(is.null(gsoa_results)){return(NA)}else{
    return(gsoa_results %>% filter(p.adjust < 0.05))
  }
  
}

do_gsoa_kegg <- function(genelist, bggenes, keyType="ENSEMBL"){
  x=org.Hs.egENSEMBL2EG
  mapped_genes <- mappedkeys(x)
  xx <- as.list(x[mapped_genes])
  
  if(keyType=="ENSEMBL"){
    genelist=unlist(xx)[genelist]
    genelist = na.omit(genelist)
    bggenes=unlist(xx)[bggenes]
    bggenes = na.omit(bggenes)}
  gsoa_results <- 
    clusterProfiler::enrichKEGG(gene=genelist,
                                organism = "hsa",
                                keyType = "kegg",
                                pvalueCutoff = 0.05,
                                pAdjustMethod = "fdr",
                                universe=bggenes,
                                minGSSize = 5,
                                maxGSSize = 500,
                                qvalueCutoff = 0.2,
                                use_internal_data = FALSE
    )
  if (is.null(gsoa_results)){res=""}else{
    res=gsoa_results %>% filter(p.adjust < 0.05)
  }
  return(res)
}


perform_ORA <-
  function(gene_associations,
           database,
           universe,
           n_clusters = 5,
           minGSSize = 10,
           maxGSSize = Inf,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.2) {
    require(clusterProfiler)
    require(multidplyr)
    
    if(n_clusters != 1) {
      worker_cluster <- new_cluster(n = n_clusters)
      cluster_library(worker_cluster, c("dplyr",
                                        "tidyverse"))
      cluster_copy(worker_cluster, c("enricher",
                                     "universe",
                                     "database",
                                     "minGSSize",
                                     "maxGSSize",
                                     "pvalueCutoff",
                                     "qvalueCutoff" ))
      
      pre_out <- 
        gene_associations %>%
        group_by(partition) %>%
        partition(worker_cluster) 
    } else {
      pre_out <- 
        gene_associations %>%
        group_by(partition)
    }
    
    outdata <-
      pre_out %>% 
      do({
        g_data <- .
        pull(g_data, gene) %>%
          enricher(universe = universe,
                   TERM2GENE = database, 
                   minGSSize = minGSSize,
                   maxGSSize = maxGSSize,
                   pvalueCutoff = pvalueCutoff,
                   qvalueCutoff = qvalueCutoff) %>%
          as_tibble()
      }) %>%
      ungroup() %>%
      collect()
    
    if(n_clusters != 1) rm(worker_cluster)
    outdata
  }


cluster_wide_data <-  
  function(wide_data,
           distance_method = "euclidean",
           clustering_method = "ward.D2", 
           cluster_rows = T,
           cluster_cols = T, 
           fill = NA) {
    suppressMessages(require(tidyverse))
    
    order_row <- 
      rownames(wide_data)
    order_col <- 
      colnames(wide_data)
    
    
    if(cluster_rows) {
      order1 <- 
        wide_data %>% 
        dist(method = distance_method) %>% 
        hclust(method = clustering_method) %>% 
        with(labels[order])
    }
    
    if(cluster_cols) {
      order2 <- 
        wide_data %>% 
        t() %>% 
        dist(method = distance_method) %>% 
        hclust(method = clustering_method) %>% 
        with(labels[order])
    }
    wide_data[order1, order2]
    
  }

gotree_plot <- function(go_res, plot_title=NA){
  
  library(rrvgo)
  library(treemapify)
  library(pbapply)
  
  plot_data_all <-
    go_res@result %>%
    filter(p.adjust < 0.05)
  
  
  plot_mat <-
    calculateSimMatrix(plot_data_all$ID,
                       orgdb="org.Hs.eg.db",
                       ont="BP",
                       method="Rel")
  plot_mat2 <-
    reduceSimMatrix(plot_mat,
                    setNames(-log10(plot_data_all$p.adjust), plot_data_all$ID),
                    threshold=0.7,
                    orgdb="org.Hs.eg.db")
  
  plot_mat2 %>%
    as_tibble() %>%
    group_by(parentTerm) %>%
    mutate(n = row_number()) %>%
    ggplot(aes(area = size, subgroup = parentTerm)) +
    
    geom_treemap(aes(fill = parentTerm,
                     alpha = n),
                 color = "black",
                 show.legend = F) +
    geom_treemap_subgroup_border(color = "black") +
    
    geom_treemap_text(aes(label = term),
                      colour = "black",
                      place = "centre",
                      alpha = 0.4,
                      grow = TRUE) +
    
    geom_treemap_subgroup_text(place = "centre",
                               grow = T,
                               reflow = T,
                               alpha = 1,
                               colour = "white",
                               fontface = "bold",
                               min.size = 0) +
    scale_alpha_continuous(range = c(0.8, 1)) +
    scale_fill_manual(values = rep(ggthemes::tableau_color_pal(palette = "Color Blind",
                                                               type = c("regular"),
                                                               direction = 1)(10),
                                   10)) +
    theme_void()+ggtitle(plot_title)
  
} 

