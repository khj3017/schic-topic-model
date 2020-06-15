suppressPackageStartupMessages(suppressWarnings({
    library(data.table)
    library(Matrix)
    library(ggplot2)
    library(cisTopic)
    library(dplyr)
    library(tibble)
    library(cluster)
    library(pheatmap)
    library(compare)
}))

cisTopicObject = readRDS("cisTopicObject_Nagano_final.rds")
cisTopicObject@selected.model = cisTopicObject@models[['25']]

# to remove coverage effect, run PCA first, remove first component
# and embed in UMAP
pca_coord = cisTopicObject@dr[['cell']]$PCA$ind.coord

custom.config = umap.defaults
custom.config$random_state = 444
custom.config$n_neighbors = 15
custom.config$min_dist = 0.3

Umap <- umap::umap(pca_coord[,-1], config = custom.config)
rownames(Umap$layout) <- rownames(pca_coord)
colnames(Umap$layout) <- paste0('UMAP', 1:ncol(Umap$layout))

umap_coord = as.data.frame(cisTopicObject@dr[['cell']]$Umap) %>% rownames_to_column()
umap_coord[,c(2,3)] = Umap$layout[,c(1,2)]
cell_data = as.data.frame(cisTopicObject@cell.data) %>% rownames_to_column()
embedding_joined = left_join(umap_coord, cell_data, by="rowname")

ggplot(embedding_joined, aes(UMAP1, UMAP2, colour = as.factor(`Cell type`))) +
    geom_point(alpha = 0.8, size = 0.1) +
    labs(colour = "") +
    theme_minimal()

cisTopicObject@dr[['cell']]$Umap = embedding_joined[,c(2,3)]
rownames(cisTopicObject@dr[['cell']]$Umap) = embedding_joined$rowname


# get topic-cell matrix
alpha <- cisTopicObject@calc.params$runCGSModels$alpha / length(cisTopicObject@selected.model$topic_sums)
topic_cell_mat <- apply(cisTopicObject@selected.model$document_sums, 2, function(x) {(x + alpha)/sum(x + alpha)})
rownames(topic_cell_mat) <- paste('Topic', seq(1,nrow(topic_cell_mat)))
colnames(topic_cell_mat) <- cisTopicObject@cell.data$Cycle

cell_types = unique(cisTopicObject@cell.data$Cycle)
num_topics = nrow(topic_cell_mat)

ks_list = list()
pval_list = list()

# Wilcoxson test
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
cell_types = cisTopicObject@cell.data$Cycle
cell_type_unique = c("G1", "Early-S", "Late-S/G2", "Pre-M", "Post-M")

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

topic_list = list()

# get cell_type specific topic list
for (cell_type in cell_type_unique) {
    topic_list[[cell_type]] = sapply(sort(rownames(filtered_pval_list[[cell_type]])), function(x) {
        as.double(strsplit(x, " ")[[1]][2])
    })
    topic_list[[cell_type]] = topic_list[[cell_type]] 
}

all_topics = 1:nrow(res)
general_topics = all_topics[-unlist(topic_list)]

topic_anno = NULL
for (cell_type in cell_type_unique) {
    if (length(topic_list[[cell_type]]) > 0) {
        topic_anno = rbind(topic_anno, 
                       data.frame(topics = topic_list[[cell_type]],
                                  type = cell_type))
    }
}

topic_anno = rbind(topic_anno, data.frame(topics = general_topics,
                        type = "General"))

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
         show_rownames = FALSE, 
         width = 4.5,
         height = 4)


##### topic to lp annotation
cisTopicObject <- getRegionsScores(cisTopicObject, method='NormTop', scale=TRUE)
cisTopicObject <- binarizecisTopics(cisTopicObject, thrP=0.99, plot=FALSE)


### plot LP distance by cell cycle stage
topic_lp_count_df = list()
plot_df = list()
num_topics = 25
resolution = 500000


for (topic in 1:num_topics) {
    binarized_data = cisTopicObject@binarized.cisTopics[[paste0("Topic", topic)]]
    bed_file = cisTopicObject@region.data[rownames(binarized_data), c('seqnames', 'start', 'end')] %>%
                dplyr::rename(chr = "seqnames")
    
    target_type = NULL
    for (cell_type in cell_types) {
        if (topic %in% topic_list[[cell_type]]) {
            target_type = cell_type
        }
    }
    
    if (is.null(target_type)) {
        target_type = "General"
    }
    
    null_df = data.frame(increment = seq(resolution, resolution * 50, by = resolution),
                         n = rep(0,50))
    
    proc_bed_file = bed_file %>% 
            mutate(increment = (end-start) + resolution) %>%
            arrange(chr, increment, start, end) %>%
            dplyr::count(increment)
    
    if (nrow(null_df) > nrow(proc_bed_file)) {
        temp = null_df
        temp[which(temp$increment %in% proc_bed_file$increment),2] = proc_bed_file$n
        proc_bed_file = temp
    }
    
    
    if (is.null(topic_lp_count_df[[target_type]])) {
        topic_lp_count_df[[target_type]] = proc_bed_file
        
    } else {
        target_df = topic_lp_count_df[[target_type]]
        
        if (nrow(null_df) > nrow(target_df)) {
            temp = null_df
            temp[which(temp$increment %in% target_df$increment),2] = target_df$n
            target_df = temp
        }
    
        target_df$n = target_df$n + proc_bed_file$n
        topic_lp_count_df[[target_type]] = target_df
    }
}



## plot cell cycle stage
plot_df = NULL
cell_types_plot = c("General", "Early-S", "Late-S/G2", "Pre-M", "Post-M", "G1")

for (cell_type in cell_types_plot) {
    plot_df = rbind(plot_df, topic_lp_count_df[[cell_type]] %>% mutate(cell_type = cell_type)) 
}

plot_df$cell_type = factor(plot_df$cell_type, levels = cell_types_plot)
ggplot(plot_df, aes(x = increment, y = n)) + 
    geom_bar(stat = "identity", colour = "black", fill = "white") + 
    scale_x_continuous(limits = c(0, 26000000),
                       breaks = seq(0, 26000000, by = 2000000),
                       labels = as.character(seq(0, 26, by = 2))) + 
    labs(x = "Distance between contact (Mb)", y = "Proportion of contacts") + 
    facet_wrap(~cell_type, ncol=2, scales = "free") +
    theme_minimal()
 
