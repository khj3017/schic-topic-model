suppressPackageStartupMessages(suppressWarnings({
    library(data.table)
    library(Matrix)
    library(ggplot2)
    library(cisTopic)
    library(dplyr)
    library(tibble)
    library(cluster)
    library(scales)
}))

cisTopicObject = readRDS("cisTopicObject_500kb_10Mb_final.rds")
umap_coord = as.data.frame(cisTopicObject@dr[['cell']]$Umap) %>% rownames_to_column()
cell_data = as.data.frame(cisTopicObject@cell.data) %>% rownames_to_column()
embedding_joined = left_join(umap_coord, cell_data, by="rowname")
embedding_joined$Color = ifelse(embedding_joined$Cluster == 5, "red", "grey80")

# plot mitotic cells in the embedding
ggplot(embedding_joined, aes(UMAP1, UMAP2)) + 
    geom_point(size = 0.1, alpha = 0.4, color = embedding_joined$Color) +  
    labs(colour = "") +
    theme_minimal() + theme(legend.position = "none")

# compare mitotic and short range contact counts between mitotic and other cells

gg = list()

gg[[1]] = ggplot(embedding_joined$Color, aes(x = Color, y = Mitotic)) +
    geom_violin(aes(fill = Color), adjust=1) +
    labs(x = "", y = "Proportion of \nmitotic contacts") + 
    theme_minimal() + theme(legend.position = "none")

gg[[2]] = ggplot(embedding_joined$Color, aes(x = Color, y = Short)) +
    geom_violin(aes(fill = Color), adjust=1) +
    labs(x = "", y = "Proportion of \nshort range contacts") + 
    theme_minimal() + theme(legend.position = "none")

do.call("grid.arrange", c(gg, ncol=2))


# plot short range contact vs mitotic contact plot
ggplot(embedding_joined, aes(Mitotic, Short)) +
    geom_point(alpha= 0.6, size = 1, color = df3$Color) +
    theme_minimal() + theme(legend.position = "none") 
