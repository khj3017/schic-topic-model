suppressPackageStartupMessages(suppressWarnings({
    library(data.table)
    library(Matrix)
    library(ggplot2)
    library(cisTopic)
    library(dplyr)
    library(tibble)
    library(cluster)
}))

# load cisTopicObject
cisTopicObject = readRDS("cisTopicObject_500kb_10Mb.rds")
# select optimal topic number
cisTopicObject@selected.model = cisTopicObject@models[['45']]
# run UMAP
cisTopicObject = runUmap(cisTopicObject, target='cell', seed=999, method = 'Probability', n_components=2)

# visualize
umap_coord = as.data.frame(cisTopicObject@dr[['cell']]$Umap) %>% rownames_to_column()
cell_data = as.data.frame(cisTopicObject@cell.data) %>% rownames_to_column()
embedding_joined = left_join(umap_coord, cell_data, by="rowname")

ggplot(embedding_joined, aes(UMAP1, UMAP2, colour = as.factor(`Cell type`))) + 
    geom_point(alpha = 0.8, size = 0.1) + 
    labs(colour = "") +
    theme_minimal()

# silhouette calculation
silhouette_idx <- silhouette(as.double(cisTopicObject@cell.data$`Cell type`), 
                             dist(cisTopicObject@dr[['cell']]$Umap))
silhouette_indexes = data.frame(cell_type = cisTopicObject@cell.data$`Cell type`,
                                 sil.val = silhouette_idx[,3])

# plot silhouette violin plot
ggplot(silhouette_indexes, aes(x = cell_type, y = sil.val, fill = cell_type)) +
    geom_violin(scale = "width") +
    labs(x = "Cell type", y = "Silhouette Index") + 
    theme_minimal()

