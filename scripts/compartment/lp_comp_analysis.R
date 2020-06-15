suppressPackageStartupMessages(suppressWarnings({
    library(data.table)
    library(Matrix)
    library(ggplot2)
    library(irlba)
    library(dplyr)
    library(tibble)
    library(grid)
    library(gridExtra)
    library(reshape2)
    library(gplots)
}))

# define parameters
topicDir = '/net/trapnell/vol1/home/khj3017/proj/scripts/cisTopic_analyses/cisTopic_beds/500000_10Mb_NormTop_0.999'
outDir = '/net/trapnell/vol1/home/khj3017/proj/scripts/cisTopic_analyses/cisTopic_beds/500000_10Mb_NormTop_0.999'
LPDir = '/net/noble/vol5/user/khj3017/sparse_matrices_LP_names/ESR7_20_LPnames.txt'

comp_dir = "/net/noble/vol5/user/khj3017/compartment/aggregated/bedgraphs_new"
resolution = 500000
cell_types = c("H1Esc", "GM12878", "HFF", "IMR90", "HAP1")

num_topics = 45
topics = 1:num_topics

topic_H1Esc = c(13,15,16,17,35,40,5)
topic_GM = c(11,12,14,18,3,39,42)
topic_HAP1 = c(10,21,24,37)
topic_HFF = c(20,28,45,9)
topic_IMR = c(31,32)

topic_list = list()
topic_list[['HAP1']] = topic_HAP1
topic_list[['H1Esc']] = topic_H1Esc
topic_list[['IMR90']] = topic_IMR
topic_list[['HFF']] = topic_HFF
topic_list[['GM12878']] = topic_GM



# Parse compartment bedgraph
read_comp_bedgraph = function(file_name) {
    comp = read.table(file_name, skip = 1, 
                      col.names = c("chr", "start", "end", "comp_value"), 
                      colClasses = c("character", "integer", "integer", "double")) %>% 
        mutate(midpoint = (end + start - 1) / 2) %>% 
        select(-c(start, end)) %>% 
        select(chr, midpoint, comp_value)
    return(comp)
}

# read compartments.bedgraph files for each cell type
comp_dfs = list()

for (cell_type in cell_types) {
    filename = file.path(comp_dir, paste0(cell_type, ".500kb.compartments.bedgraph"))
    comp_dfs[[cell_type]] = read_comp_bedgraph(filename)
}


# read LP file that is used for cisTopic
LP_file = read.table(LPDir, sep = ":", 
           col.names = c("chr", "lp"), 
           colClasses = c("character", "character")) %>%
    rowwise() %>%
    mutate(start = as.double(strsplit(lp, "-")[[1]][1]),
           end = as.double(strsplit(lp, "-")[[1]][2])) %>%
    as.data.frame() %>%
    select(-lp)

## Assign compartment to LPs
join_LP_comp_value = function(LP_file, comp) {
    comp_start = comp %>% rename(start = midpoint)
    comp_end = comp %>% rename(end = midpoint)
    LP_file = LP_file %>%
        left_join(., comp_start, by = c("chr", "start")) %>%
        rename(start_comp_value = comp_value) %>%
        left_join(., comp_end, by = c("chr", "end")) %>%
        rename(end_comp_value = comp_value) %>%
        rowwise() %>%
        filter(start_comp_value != 0 & end_comp_value != 0) %>%
        mutate(compartment = 
              if (sign(start_comp_value) == sign(end_comp_value)) {
                  if (start_comp_value > 0) {
                      "AA"
                  } else {
                      "BB"
                  }
              } else {
                  "AB"
              }
        ) %>%
        as.data.frame()
    return(LP_file)
}

## Assign compartment to LPs
comp_assigned_dfs = list()

for (cell_type in cell_types) {
    comp_assigned_dfs[[cell_type]] = join_LP_comp_value(LP_file, comp_dfs[[cell_type]])
}


## Annotate topics to cell types

out_df = NULL

for (topic in topics) {
    is_cell_type_specific = F
    cell_list = NULL
    anno = "general"
    
    for (cell_type in names(topic_list)) {
        if (topic %in% topic_list[[cell_type]]) {
            cell_list = c(cell_list, cell_type)
        }
    }
    
    if (length(cell_list) == 1) {
        is_cell_type_specific = T
    }
    
    if (is_cell_type_specific) {
        anno = cell_list
    } 
    
    filename = file.path(topicDir, paste0("Topic_", topic, ".bed"))
    bed_file = read.table(filename, 
                    col.names = c("chr", "start", "end"), 
                    colClasses = c("character", "integer", "integer")) %>%
                    mutate(topic = topic, 
                           cell_type = anno)
    
    
    if (is.null(out_df)) {
        out_df = bed_file
    } else {
        out_df = rbind(out_df, bed_file)
    }
        
}

# Identify LPs that are shared among cell types
duplicate_df = out_df %>%
    group_by(chr, start, end, cell_type) %>%
    add_tally(name = "lp_count") %>% filter(lp_count > 1) %>%
    as.data.frame() %>%
    distinct(chr, start, end, .keep_all = T) %>%
    mutate(LP = paste0(chr, ":", start, "-", end))



## Get LP comp assignment for each topic for each cell type
topic_res = list()
plot_df = list()
background_df = list()

for (cell_type in cell_types) {
    topic_res[[cell_type]] = list()
    comp_file = comp_dfs[[cell_type]]
    out_df = NULL
    
    for (topic in topics) {
        filename = file.path(topicDir, paste0("Topic_", topic, ".bed"))
        bed_file = read.table(filename, 
                        col.names = c("chr", "start", "end"), 
                        colClasses = c("character", "integer", "integer")) 
        
        if (nrow(bed_file) > 0) {
            joined_res = join_LP_comp_value(bed_file, comp_file) %>%
                                mutate(LP = paste0(chr, ":", start, "-", end)) %>%
                                filter(!LP %in% duplicate_df$LP)
            
            topic_res[[cell_type]][[as.character(topic)]] = joined_res
            
            count_joined_res = joined_res %>% 
                                count(compartment)    
            
            if (is.null(out_df)) {
                if (nrow(count_joined_res) < 3) {
                    temp = data.frame(compartment = c("AA", "AB", "BB"),
                                      n = rep(0,3))
                    temp[temp$compartment %in% count_joined_res$compartment,]$n = count_joined_res$n
                    count_joined_res = temp
                }
                out_df = count_joined_res
            } else {
                out_df = left_join(out_df, count_joined_res, by = "compartment") %>% 
                                mutate_each(funs(replace(., which(is.na(.)), 0)))       
            }
            
        }
    }
    
    background_df[[cell_type]] = data.frame(compartment = out_df$compartment,
           n = rowSums(out_df[,-1]))
}




# converts topic data to count data
normalize_comp_df = function(comp_assigned_dfs, isTopic = F, normalize = T) {
    if (isTopic) {
        for (cell_type in names(comp_assigned_dfs)) {
               for (topic in names(comp_assigned_dfs[[cell_type]])) {
                   comp_assigned_dfs[[cell_type]][[topic]] = comp_assigned_dfs[[cell_type]][[topic]] %>% 
                                count(compartment)
                   
                   if (normalize) {
                        comp_assigned_dfs[[cell_type]][[topic]] = comp_assigned_dfs[[cell_type]][[topic]] %>% 
                                mutate(n = n / sum(n))
                   }
               }
        }
    } else {
        for (cell_type in names(topic_res)) {
            comp_assigned_dfs[[cell_type]] = comp_assigned_dfs[[cell_type]] %>% 
                count(compartment) 
            
            if (normalize) {
                comp_assigned_dfs[[cell_type]] = comp_assigned_dfs[[cell_type]] %>% 
                    mutate(n = n / sum(n))
            }
        }
    }
    
    return(comp_assigned_dfs)
}


perform_chi_square_test = function(topic_res, background_df, signif = 0.01) {
    normalized_topic_res = normalize_comp_df(topic_res, isTopic = T, normalize = F)
    normalized_background_df = normalize_comp_df(background_df, isTopic = F, normalize = T)
    
    res_df = list()
    
    for (cell_type in names(topic_res)) {
        display(cell_type)
        res = sapply(names(topic_res[[cell_type]]), function(topic) {
            
            target_df = normalized_topic_res[[cell_type]][[topic]]
            background_df = normalized_background_df[[cell_type]]
            
            if (nrow(target_df) < 3) {
                comp_to_fill = background_df[!( background_df$compartment %in% target_df$compartment),]$compartment
                zero_df = data.frame(compartment = comp_to_fill,
                                     n = rep(0, length(comp_to_fill)))
                target_df = rbind(target_df, zero_df)
                target_df = target_df %>% arrange(as.character(compartment))
            }
            
            if (nrow(target_df) == 0) target_df = data.frame(AA = 0, AB = 0, BB = 0)
            
            chisq_res = chisq.test(target_df$n, 
                                   p = background_df$n)
            res = c(chisq_res$statistic, chisq_res$p.value)
            return(res)
        })
        
        res = as.data.frame(t(res))
        colnames(res) = c("statistic", "pvalue")
        rownames(res) = paste0("Topic ", 1:length(topic_res[[cell_type]]))
        res$qvalue = p.adjust(res$pvalue, method = "BH")
        
        res_df[[cell_type]] = res
    }
    
    return(res_df)
}

get_log2_obs_exp = function(topic_res, background_df, diff = F) {
    normalized_topic_res = normalize_comp_df(topic_res, isTopic = T, normalize = T)
    
    if (diff) {
        normalized_background_df = background_df
    } else {
        normalized_background_df = normalize_comp_df(background_df, isTopic = F, normalize = T)
    }
    
    
    out_df = list()
    for (cell_type in names(topic_res)) {
        out_df[[cell_type]] = list()
        background_df_cell = normalized_background_df[[cell_type]]
        
        if (nrow(background_df_cell) < 3) {
            comp_to_fill = c("AA", "AB", "BB")
            zero_df = data.frame(compartment = comp_to_fill,
                                 n = rep(0, length(comp_to_fill)))
            zero_df[zero_df$compartment %in% background_df_cell$compartment,]$n = background_df_cell$n 
            background_df_cell = zero_df %>% arrange(as.character(compartment))
        }
        
        
        if (diff) {
            background_df_cell$n = background_df_cell$n / sum(background_df_cell$n)
        }
        
        for (topic in names(topic_res[[cell_type]])) {
            target_df = normalized_topic_res[[cell_type]][[topic]]
            
            if (nrow(target_df) < 3) {
                comp_to_fill = background_df_cell[!(background_df_cell$compartment %in% target_df$compartment),]$compartment
                zero_df = data.frame(compartment = comp_to_fill,
                                     n = rep(0, length(comp_to_fill)))
                target_df = rbind(target_df, zero_df)
                target_df = target_df %>% arrange(as.character(compartment))
            }
            
            
            out_df[[cell_type]][[topic]] = target_df
            target_df$n = 
                log2((target_df$n + 1e-5) / (background_df_cell$n+ 1e-5))
            target_df$n[target_df$n < -3] = -3
            target_df$n[target_df$n > 3] = 3
            
            out_df[[cell_type]][[topic]]$n = target_df$n
        }
    }
    
    return(out_df)
}


chisq_df = perform_chi_square_test(topic_res, background_df)
log2_obs_exp = get_log2_obs_exp(topic_res, background_df)


# plot lp topic compartment heatmaps for each cell type
plot_heatmap_topics = function(log2_obs_exp, chisq_df) {
    for (cell_type in cell_types) {
        mat = matrix(, nrow = num_topics, ncol = 3)
        for (topic in 1:length(names(log2_obs_exp[[cell_type]]))) {
            mat[topic,] = log2_obs_exp[[cell_type]][[topic]]$n
        }

        rownames(mat) = paste0("Topic ", topics)
        colnames(mat) = c("AA", "AB", "BB")

        signif_topics = chisq_df[[cell_type]] %>%
                            rownames_to_column("topic") %>%
                            arrange(topic) %>%
                            filter(qvalue < 0.01) %>%
                            pull(topic)

        cell_type_topics = intersect(signif_topics, paste0("Topic ", topic_list[[cell_type]]))

        mat = as.matrix(mat[cell_type_topics, , drop=F])
        if (nrow(mat) == 1) mat = rbind(mat, mat)
	n = 500
        color_list = rev(colorRampPalette(RColorBrewer::brewer.pal(7, "YlGnBu"))(n))
        
        
        pdf_name = paste0("comp_heatmap_", cell_type, ".pdf")
        pdf(file = pdf_name, width = 3, height = 3)
        heatmap.2(mat, Colv = F, key.par = list(mar=c(1,1,1,1)),
                  col=color_list, #bluered(100),
                  dendrogram = "row", 
                  density.info = "none",
                  trace = "none",
                  srtCol = -45,
                  adjCol = c(.1, .5),
                  scale = "none",
                  xlab = "",
                  ylab = "",
                  cexRow = 0.8,
                  cexCol = 1,
                  key=F,
                  breaks = seq(-1.5,1.5, length.out = n + 1))
        
        dev.off()
    }
}

plot_heatmap_topics(log2_obs_exp, chisq_df)

## see what lp compartments are enriched
get_num_topics_enriched = function(log2_obs_exp) {
    num_topics_enriched = list()
    
    for (cell_type in cell_types) {
        mat = matrix(, nrow = num_topics, ncol = 3)
        for (topic in 1:length(names(log2_obs_exp[[cell_type]]))) {
            mat[topic,] = log2_obs_exp[[cell_type]][[topic]]$n
        }

        rownames(mat) = paste0("Topic ", 1:num_topics)
        colnames(mat) = c("AA", "AB", "BB")

        mat = mat[topic_list[[cell_type]], ]
        
        mat_filtered = apply(mat, 1, function(row) {
            which(row > log2(1.25))
        })
        
        if (is.list(mat_filtered)) {
            temp = NULL
            for (i in mat_filtered) {
                if (length(i) == 2) {
                    temp = c(temp, 2)
                } else {
                    temp = c(temp, i)
                }
            }
        } else {
            temp = mat_filtered
        }
        
        num_topics_enriched[[cell_type]] = temp
        
    }
    return(num_topics_enriched)
}

LP_comp_count = get_num_topics_enriched(log2_obs_exp)
comp_count_concat = NULL

# concatenate lp comp counts 
for (type in names(LP_comp_count)) {
    lp_count = LP_comp_count[[type]]
    out_temp = data.frame(temp = lp_count) %>%
                    dplyr::count(temp, name = "comp_count") %>%
                    rowwise() %>%
                    mutate(comp = if (temp == 1) {
                        "AA only"
                    } else if (temp == 2) {
                        "Both AA and BB"
                    } else {
                        "BB only"
                    }) %>%
                    select(-temp) %>%
                    select(2,1) %>%
                    as.data.frame() %>%
                    mutate(cell_type = type,
                           n = sum(comp_count)) %>%
                    mutate(comp_count = comp_count / n)
    
    comp_count_concat = rbind(comp_count_concat, out_temp)
}

## plot comp_count_concat
ggplot(comp_count_concat, aes(x = cell_type, y = comp_count, fill = comp)) +
    geom_bar(stat = "identity", color = "black") +
    labs(x = "", y = "") + 
    theme_minimal()

