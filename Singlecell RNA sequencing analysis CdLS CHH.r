library(Seurat)
library(ggplot2)

#IPS, NPCs and Neurons
seurat_obj <- readRDS("./seurat_obj.rds")
DimPlot(seurat_obj, reduction = "umap", group.by = "sample")
DimPlot(seurat_obj, reduction = "umap", group.by = "origin")
DimPlot(seurat_obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE, label.size = 5)


#########################Modules######################## 
path_to_genes <- "./all_genes.csv"
genes <- read.csv(path_to_genes, header = TRUE)
genes$gene_Name_ID <- paste(genes$gene_name, genes$gene_id, sep = ",")
genes <- genes %>%
  mutate(unique_gene_names = ifelse(duplicated(gene_name) | duplicated(gene_name, fromLast = TRUE), 
                                    paste(gene_id, gene_name, sep = "_"),
                                    gene_name))
all_markers <- read.csv("./module_markers.csv")
all_markers <- all_markers %>% left_join(genes, by = c("gene_name" = "gene_name"))

modules <- all_markers %>%
  group_by(cell_type) %>%
  summarise(gene_list = list(gene_name)) %>%
  deframe()

seurat_obj <- AddModuleScore(seurat_obj, features = modules, name = "Module")
module_score_columns <- colnames(seurat_obj@meta.data)[grepl("Module", colnames(seurat_obj@meta.data))]
original_module_names <- names(modules)  
generated_module_names <- module_score_columns  
module_name_mapping <- setNames(generated_module_names, original_module_names)

for (i in seq_along(module_score_columns)) {
  colnames(seurat_obj@meta.data)[colnames(seurat_obj@meta.data) == module_score_columns[i]] <- original_module_names[i]
}

module_order <- c("Stemness", "Proliferation", "NPC", "Neural crest", "Astrocytes","Mature neurons", "Glut", "GnRH-GABA","OSN") 
module_order <- rev(module_order)
valid_modules <- intersect(module_order, colnames(seurat_obj@meta.data))

data_to_plot <- seurat_obj@meta.data %>%
  dplyr::select(seurat_clusters, all_of(valid_modules)) %>%
  pivot_longer(
    cols = -seurat_clusters, 
    names_to = "module", 
    values_to = "expression"
  ) %>%
  mutate(
    module = factor(module, levels = module_order), 
    seurat_clusters = factor(seurat_clusters, levels = c("0", "1", "2", "3", "4", "5", "6", "7")) 
  )

dot_plot_data <- data_to_plot %>%
  group_by(seurat_clusters, module) %>%
  summarise(
    avg_expression = mean(expression, na.rm = TRUE),
    pct_expression = mean(expression > 0, na.rm = TRUE),
    .groups = 'drop'
  )

print(head(dot_plot_data))

p <- ggplot(dot_plot_data, aes(x = seurat_clusters, y = module, size = pct_expression, fill = avg_expression)) +
  geom_point(pch = 21, stroke = 0.5) +
  scale_size_area() +
  scale_fill_viridis_c(option = "C") +
  theme_minimal() +
  labs(x = "    ", y = "      ", size = "Percentage of Cells", fill = "Average Expression") +
  theme(
    axis.text.x = element_text(hjust = 1, size = 15, color = "black", face = "bold"), 
    axis.text.y = element_text(size = 15, color = "black", face = "bold"),
    axis.title.y = element_text(color = "black", face = "bold")
  )
print(p)


#####################fig:Expression of the genes from the modules
genes_df <- read.csv("./all_genes.csv", header = TRUE, stringsAsFactors = FALSE)
genes_df <- genes_df %>%
  mutate(unique_gene_names = ifelse(duplicated(gene_name) | duplicated(gene_name, fromLast = TRUE),
                                    paste(gene_id, gene_name, sep = "_"),
                                    gene_name))
all_markers2 <- read.csv("./module_markers.csv", header = TRUE, stringsAsFactors = FALSE) %>%
  merge(genes_df[, c("gene_name", "unique_gene_names")], by = "gene_name", all.x = TRUE)
modules_of_interest <- c(
  "Stemness", "Proliferation", "NPC", "Neural crest",
  "Astrocytes", "Mature neurons", "Glut", "GABAergic neurons", "GnRH-GABA", "OSN"
)
gene_module_map <- all_markers2 %>%
  filter(cell_type %in% modules_of_interest) %>%
  mutate(cell_type = factor(cell_type, levels = modules_of_interest)) %>%
  arrange(cell_type) %>%
  distinct(unique_gene_names, .keep_all = TRUE) %>%
  filter(unique_gene_names %in% rownames(seurat_obj))

genes_to_plot <- gene_module_map$unique_gene_names
if (length(genes_to_plot) == 0) {
  stop("check unique_gene_names and seurat_obj rownames.")
}
seurat_obj_subset <- subset(seurat_obj, subset =
                              grepl("H0", origin, ignore.case = TRUE) |
                              origin == "P02" |
                              origin == "P01")

expr_df <- FetchData(seurat_obj_subset, vars = c("seurat_clusters", "origin", genes_to_plot))
expr_long <- expr_df %>%
  rownames_to_column("cell") %>%
  pivot_longer(cols = -c(cell, seurat_clusters, origin), names_to = "gene", values_to = "expression") %>%
  mutate(
    gene = as.character(gene),
    seurat_clusters = factor(seurat_clusters, levels = sort(unique(seurat_clusters))),
    origin_group = case_when(
      grepl("H0", origin, ignore.case = TRUE) ~ "Control",
      origin == "P02" ~ "P02 (NIPBL)",
      origin == "P01" ~ "P01 (SMC3)"
    ),
    origin_group = factor(origin_group, levels = c("Control", "P01 (SMC3)", "P02 (NIPBL)"))
  ) %>%
  left_join(gene_module_map[, c("unique_gene_names", "cell_type")],
            by = c("gene" = "unique_gene_names")) %>%
  mutate(
    gene = factor(gene, levels = genes_to_plot),
    module = recode(cell_type,
                    "Stemness" = "Stemness",
                    "Proliferation" = "Prolif",
                    "NPC" = "NPC",
                    "Neural crest" = "Neural\ncrest",
                    "Astrocytes" = "Astro\ncytes",
                    "Mature neurons" = "Mature\nneurons",
                    "Glut" = "Glut",
                    "GnRH-GABA" = "GnRH-GABA",
                    "OSN" = "OSN"),
    gene_label = gene
  )

modules_of_interest_short <- c("Stemness", "Prolif","NPC", "Neural\ncrest", "Astro\ncytes", "Mature\nneurons",
                               "Glut","GnRH-GABA",  "OSN")
expr_long$module <- factor(expr_long$module, levels = modules_of_interest_short)
expr_long$gene_label <- factor(expr_long$gene_label,
                               levels = unique(expr_long$gene_label[order(expr_long$gene)]))
expr_long <- expr_long %>% mutate(origin_group = as.character(origin_group))
p <- ggplot(expr_long, aes(x = gene_label, y = expression))
if (nrow(filter(expr_long, origin_group == "Control")) > 0) {
  p <- p + gghalves::geom_half_violin(
    data = filter(expr_long, origin_group == "Control"),
    aes(fill = origin_group),
    side = "l", trim = TRUE, scale = "width", color = NA, width = 0.8
  )
}
if (nrow(filter(expr_long, origin_group == "P01 (SMC3)")) > 0) {
  p <- p + gghalves::geom_half_violin(
    data = filter(expr_long, origin_group == "P01 (SMC3)"),
    aes(fill = origin_group),
    side = "r", trim = TRUE, scale = "width", color = NA, width = 0.8
  )
}
if (nrow(filter(expr_long, origin_group == "P02 (NIPBL)")) > 0) {
  p <- p + gghalves::geom_half_violin(
    data = filter(expr_long, origin_group == "P02 (NIPBL)"),
    aes(fill = origin_group),
    side = "r", trim = TRUE, scale = "width", color = NA, width = 0.8, alpha = 0.5
  )
}

p <- p +
  facet_grid(seurat_clusters ~ module, scales = "free_x", space = "free_x", switch = "y") +
  labs(x = "Gene", y = "Expression Level", fill = NULL) +
  scale_fill_manual(
    values = c("Control" = "#66c2a5", "P01 (SMC3)" = "blue", "P02 (NIPBL)" = "red"),
    labels = c(
      "Control" = "H01 & H02",
      "P01 (SMC3)" = "P01",
      "P02 (NIPBL)" = "P02"
    )
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_markdown(angle = 60, hjust = 1, size = 15),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    strip.text.x = element_text(size = 15, face = "bold"),
    strip.background = element_blank(),
    strip.text.y = element_text(size = 15, face = "bold"),
    #strip.text.y.left = element_text(size = 15, face = "bold"),
    strip.placement = "outside",
    panel.spacing = unit(0.4, "lines"),
    plot.margin = margin(10, 10, 20, 10),
    legend.text = element_markdown(size = 15),
    legend.key.size = unit(0.8, "cm"),
    legend.position = "bottom",
    legend.justification = "center"
  )

print(p)


#IPS
seurat_ips  <- readRDS(file.path(indir, "seurat_ips.rds"))
DimPlot(seurat_ips, reduction = "umap", group.by = "seurat_clusters") + ggtitle("IPS_clusters")
DimPlot(seurat_ips, reduction = "umap", group.by = "origin") + ggtitle("IPS_origin")
genes_to_plot <- c( "TDGF1","POU5F1","NANOG","SIX3","OTX2","PAX6")
feature_plots <- lapply(genes_to_plot, function(gene) {
  FeaturePlot(seurat_ips, features = gene, reduction = "umap") + 
    ggtitle(paste(gene)) +
    theme(plot.title = element_text(hjust = 0.5, size = 14))
})
Fig_ips_markergenes <- wrap_plots(feature_plots, ncol = 3)  
Fig_ips_markergenes


#NPCs
seurat_npcs <- readRDS(file.path(indir, "seurat_npcs.rds"))
DimPlot(seurat_npcs, reduction = "umap", group.by = "origin") 
genes_to_plot <- c( "SOX2","PAX6","PAX3","MSX1")
feature_plots <- lapply(genes_to_plot, function(gene) {
  FeaturePlot(seurat_npcs, features = gene, reduction = "umap") + 
    ggtitle(paste(gene)) +
    theme(plot.title = element_text(hjust = 0.5, size = 14))
})
Fig_npcs_markergenes <- wrap_plots(feature_plots, ncol = 2)  
Fig_npcs_markergenes


#Neur
seurat_neur <- readRDS(file.path(indir, "seurat_neur.rds"))
DimPlot(seurat_neur, reduction = "umap", group.by = "origin") 
genes_to_plot <- c( "GAD1","GAD2","FOXG1","ISL1","DLX5")
feature_plots <- lapply(genes_to_plot, function(gene) {
  FeaturePlot(seurat_neur, features = gene, reduction = "umap") + 
    ggtitle(paste(gene)) +
    theme(plot.title = element_text(hjust = 0.5, size = 14))
})
Fig_neur_markergenes <- wrap_plots(feature_plots, ncol = 3)  
Fig_neur_markergenes







