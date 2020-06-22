suppressPackageStartupMessages(suppressWarnings({
    library(ggplot2)
    library(cisTopic)
    library(dplyr)
    library(tibble)
    library(pheatmap)
}))

# load cisTopicObject
cisTopicObject = readRDS("cisTopicObject_500kb_10Mb_final.rds")


# get topic-cell matrix
alpha <- cisTopicObject@calc.params$runCGSModels$alpha / length(cisTopicObject@selected.model$topic_sums)
topic_cell_mat <- apply(cisTopicObject@selected.model$document_sums, 2, function(x) {(x + alpha)/sum(x + alpha)})
rownames(topic_cell_mat) <- paste('Topic', seq(1,nrow(topic_cell_mat)))
colnames(topic_cell_mat) <- cisTopicObject@cell.data$`Cell type`

cell_types = unique(cisTopicObject@cell.data$`Cell type`)
num_topics = nrow(topic_cell_mat)

# plot topic contributions
plot_distribution = function(topic_cell_mat, topic, xlims = NULL) {
    cell_types = colnames(topic_cell_mat)
    plot_df = data.frame(values = topic_cell_mat[topic,],
                         cell_type = cell_types)
    
    if (is.null(xlims)) {
        xlims = c(0, max(plot_df$values))
    } else {
        xlims = xlims
    }
    
    g = ggplot(plot_df, aes(x = values, colour = cell_type)) + 
            geom_density() +
            xlim(xlims) + 
            labs(x = "Normalized topic assignment", y = "Density", colour = "Cell type") +
            theme_minimal()
    return(g)
    
}

plot_distribution(topic_cell_mat, 10)
plot_distribution(topic_cell_mat, 5)


# Wilcoxson test
ks_list = list()
pval_list = list()

for (cell_type in cell_types) {
  idx = grepl(cell_type, colnames(topic_cell_mat))
  red_mat = topic_cell_mat[,idx]
  nx = length(idx)

  ks_stats_f = NULL
  pval_f = NULL  
    
  for (cell_type2 in cell_types) {
    idx2 = grepl(cell_type2, colnames(topic_cell_mat))
    if (cell_type == cell_type2) next
    rest_mat = topic_cell_mat[,idx2]
    num_topics = dim(topic_cell_mat)[1]
    ny = length(idx2)
    dim(red_mat)[2]
      
      
    ks_stats <- vector(, num_topics)
    p_values <- vector(, num_topics)

    for (topic in 1:(num_topics)) {
      topic_vec = red_mat[topic,]
      topic_vec2 = rest_mat[topic,]
      ks_stats[topic] = wilcox.test(topic_vec, topic_vec2, alternative="g")$statistic
      p_values[topic] = wilcox.test(topic_vec, topic_vec2, alternative="g")$p.value
    }
    
    p_values = p.adjust(p_values, method = "BH")
      
    if (is.null(ks_stats_f)) {
        ks_stats_f = ks_stats
        pval_f = p_values
    } else {
        ks_stats_f = rbind(ks_stats_f, ks_stats)
        pval_f = rbind(pval_f, p_values)
    } 
  }
  rownames(ks_stats_f) = cell_types[-which(cell_types == cell_type)]
  colnames(ks_stats_f) = paste0("Topic ", 1:ncol(ks_stats_f))
  
  rownames(pval_f) = rownames(ks_stats_f)
  colnames(pval_f) = colnames(ks_stats_f)
    
  ks_list[[cell_type]] = ks_stats_f
  pval_list[[cell_type]] = pval_f
}      

# assign topics to cell types
filtered_ks_list = list()
filtered_pval_list = list()

for (cell_type in names(pval_list)) {
    cell_types_filtered = cell_types[-which(cell_types == cell_type)]
    filtered_pval_list[[cell_type]] = as.data.frame(t(pval_list[[cell_type]])) %>%
                                rownames_to_column() %>%
                                rowwise() %>%
                                filter(max(get(cell_types_filtered[1]), 
                                           get(cell_types_filtered[2]), 
                                           get(cell_types_filtered[3]), 
                                           get(cell_types_filtered[4])) < 1e-2) %>%
                                column_to_rownames() %>%
                                as.data.frame()
    filtered_ks_list[[cell_type]] = as.data.frame(t(ks_list[[cell_type]])) %>%
                                rownames_to_column() %>%
                                rowwise() %>%
                                filter(max(get(cell_types_filtered[1]), 
                                           get(cell_types_filtered[2]), 
                                           get(cell_types_filtered[3]), 
                                           get(cell_types_filtered[4])) < 1e-2) %>%
                                column_to_rownames() %>%
                                as.data.frame()
}


# get aggregated cell type-topic matrix
cell_types = cisTopicObject@cell.data$`Cell type`
cell_type_unique = c("GM12878", "H1Esc", "HAP1", "HFF", "IMR90")

res = apply(topic_cell_mat, 1, function(col) {
    res = NULL
    for (cell_type in cell_type_unique) {
        idx = cell_types == cell_type
        res = c(res, mean(col[idx]))
    }
    return(res)
})

res = t(res)
colnames(res) = cell_type_unique
rownames(res) = 1:nrow(res)

cutoff = 1.5
topic_list = list()

# get cell_type specific topic list 
# and apply the second filter
for (cell_type in cell_type_unique) {
    topic_list[[cell_type]] = sapply(sort(rownames(filtered_pval_list[[cell_type]])), function(x) {
        as.double(strsplit(x, " ")[[1]][2])
    })
    filtered = topic_list[[cell_type]] %in% rownames(res[apply(res, 1, function(row) {
        max_val = max(row)
        second_max = sort(row, decreasing = T)[2]
        return(max_val / second_max)
    }) > cutoff, ])
    topic_list[[cell_type]] = topic_list[[cell_type]][filtered]
}

# reorder cell type topic mat
res = res[topic_anno$topics,]
rownames(topic_anno) = rownames(res)


## plot cell type-topic heatmap
topic_indices = topic_anno %>% dplyr::count(type) %>% pull(n)

for (idx in 2:length(topic_indices)) {
    topic_indices[idx] = topic_indices[idx] + topic_indices[idx-1]
}

pheatmap(res, 
         scale = "row",
         gaps_row = topic_indices[-length(topic_indices)],
         gaps_col = c(1,2,3,4),
         annotation_row = topic_anno %>% dplyr::select(type),
         annotation_names_row=F,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_colnames = FALSE,
         show_rownames = FALSE)


##### topic to lp annotation
cisTopicObject <- getRegionsScores(cisTopicObject, method='NormTop', scale=TRUE)
cisTopicObject <- binarizecisTopics(cisTopicObject, thrP=0.9975, plot=F)


out_dir = "Bed_500000_10Mb_NormTop_0.9975"
if (!file.exists(out_dir)) {
    dir.create(out_dir)
}
getBedFiles(cisTopicObject, path=out_dir)
