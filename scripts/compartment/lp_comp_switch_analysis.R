suppressPackageStartupMessages(suppressWarnings({
    library(data.table)
    library(Matrix)
    library(ggplot2)
    library(dplyr)
    library(tibble)
    library(grid)
    library(gridExtra)
    library(reshape2)
    library(gplots)
    library(pheatmap)
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

# filter locus pairs
filtered_out_df = out_df %>% 
    arrange(chr, start, end) %>%
    group_by(chr, start, end, cell_type) %>% add_tally() %>% 
    filter(n == 1) %>% select(-n) %>%
    distinct(chr, start, end, cell_type, .keep_all=T) %>%
    as.data.frame()

# locus pairs in all topics
all_lps = out_df %>%
    arrange(chr, start, end) %>%
    group_by(chr, start, end) %>% add_tally() %>% 
    filter(n == 1) %>% select(-n) %>%
    distinct(chr, start, end, .keep_all=T) %>%
    as.data.frame()

# concatenate compartment calls
concat_comp_dfs = function(comp_dfs) {
    out = NULL
    for (cell_type in names(comp_dfs)) {
        temp = comp_dfs[[cell_type]] %>%
                dplyr::rename(!!cell_type := comp_value)
        if (is.null(out)) {
            out = temp
        } else {
            out = left_join(out, temp, by = c("chr", "midpoint"))
        }
    }
    return(out)
}

concatenated_comp_df = concat_comp_dfs(comp_dfs)


# get compartment switching regions
get_comp_switches = function(concatenated_comp_df) {
    cell_types = c("GM12878", "H1Esc", "HFF", "IMR90", "HAP1")
    out = list()
    for (type in cell_types) {
        cell_types_filtered = cell_types[-which(cell_types == type)]
        out[[type]] = concatenated_comp_df %>%
            filter(get(type) * get(cell_types_filtered[1]) < 0 & 
                   get(type) * get(cell_types_filtered[2]) < 0 & 
                   get(type) * get(cell_types_filtered[3]) < 0 & 
                   get(type) * get(cell_types_filtered[4]) < 0)
    }
    return(out)
}

comp_switch_df = get_comp_switches(concatenated_comp_df)


## annotate lps whether they contain compartment switching regions
# CS = compartment switching
annotated_df = NULL
annotated_general_df = NULL
annotated_df_all = NULL

annotate_lps = function(lp_df, type) {
    comp_temp = comp_switch_df[[type]] %>%
                    rowwise() %>%
                    mutate(comp = paste0(chr, "-", midpoint)) %>% as.data.frame()
    out = lp_df %>%
        rowwise() %>%
        mutate(Start_CS = ifelse(paste0(chr, "-", start) %in% comp_temp$comp, "Y", "N"),
               End_CS = ifelse(paste0(chr, "-", end) %in% comp_temp$comp, "Y", "N")) %>%
        mutate(CS = ifelse(Start_CS == "Y" | End_CS == "Y", "Y", "N")) %>%
        as.data.frame()
    return(out)
}

for (type in names(comp_switch_df)) { 
    annotated_df[[type]] = annotate_lps(filtered_out_df %>% filter(cell_type == type), type)
    annotated_general_df[[type]] = annotate_lps(filtered_out_df %>% filter(cell_type == "general"), type)
    annotated_df_all[[type]] = annotate_lps(all_lps, type)
}



# compute enrichment of compartment switching regions using lps from all topics as baseline
comp_switch_enrichment = list()
temp = all_lps

for (type in names(comp_switch_df)) { 
    
    backgd = annotated_df_all[[type]] %>%
        count(CS) %>%
        mutate(n = n / sum(n)) %>%
        as.data.frame()
    
    
    obs = annotated_df[[type]] %>% 
                count(CS) %>%
                mutate(n = n / sum(n)) %>%
                as.data.frame()
     
    ddf = data.frame(CS = backgd$CS,
                      n = rep(0, length(backgd$CS)))
    
    ddf[which(ddf$CS %in% obs$CS),]$n = obs$n
    
    print(chisq.test(ddf$n, p = backgd$n))
    
    ddf$n = log2(ddf$n / backgd$n)
    ddf = ddf %>% 
            mutate_if(is.numeric, function(x) ifelse(is.infinite(x), 0, x)) %>%
            mutate(sample_size = annotated_df[[type]] %>% count(CS) %>% pull(n))

    comp_switch_enrichment[[type]] = ddf
}

# concatenate enrichment results
comp_switch_enrichment_concatenated = NULL
for (type in names(comp_switch_enrichment)) {
    comp_switch_enrichment_concatenated = rbind(comp_switch_enrichment_concatenated,
                                                comp_switch_enrichment[[type]] %>%
                                                    mutate(cell_type = type))
}


# plot concatenated enrichment results
offset = 0.15

ggplot(comp_switch_enrichment_concatenated, 
       aes(x = cell_type, y = n, fill = CS)) + 
        geom_bar(stat = "identity", position = "dodge", colour = "black") + 
        geom_text(aes(y = ifelse( sign(n) > 0, n + offset, n - offset ),
                      label = sample_size), 
             position = position_dodge(width=0.9)) + 
        labs(x = "", y = "log2 ratio") + 
        scale_fill_manual(values = c("mediumaquamarine", "tan1")) + 
	theme_minimal()



## compute enrichment of A/B compartments in identified comp. switching regions
comp_switch_enrichment = list()

for (type in names(comp_switch_df)) { 
    
    comp_temp = comp_switch_df[[type]] %>% 
                        rowwise() %>%
                        mutate(comp = paste0(chr, "-", midpoint)) %>% as.data.frame() %>%
                        mutate(assign = ifelse(.[[type]] > 0, "A", "B"))

    

    backgd = annotated_df_all[[type]] %>%
        filter(CS == "Y") %>%
        rowwise() %>%
        mutate(Start_Assign = ifelse(Start_CS == "Y", 
                               comp_temp[comp_temp$midpoint == start,]$assign, "")) %>%
        mutate(End_Assign = ifelse(End_CS == "Y", 
                               comp_temp[comp_temp$midpoint == end,]$assign, "")) %>%
        mutate(assign = paste0(Start_Assign, End_Assign)) %>%
        mutate(assign = ifelse(assign == "BA", "AB", assign)) %>%
        as.data.frame() %>%
        count(assign) %>%
        tidyr::spread(assign, n)

    obs = annotated_df[[type]] %>% 
        filter(CS == "Y") %>%
        rowwise() %>%
        mutate(Start_Assign = ifelse(Start_CS == "Y", 
                               comp_temp[comp_temp$midpoint == start,]$assign, "")) %>%
        mutate(End_Assign = ifelse(End_CS == "Y", 
                               comp_temp[comp_temp$midpoint == end,]$assign, "")) %>%
        mutate(assign = paste0(Start_Assign, End_Assign)) %>%
        mutate(assign = ifelse(assign == "BA", "AB", assign)) %>%
        as.data.frame() %>%
        count(assign) %>%
        tidyr::spread(assign, n)
        mutate(n = n / sum(n))) 
    
    
    ddf = data.frame(A = 0, AA = 0, AB = 0,
                      B = 0, BB = 0)
    
    ddf[1,which(colnames(ddf) %in% colnames(obs))] = obs
    
    df = data.frame(A = ddf$A + 2*ddf$AA + ddf$AB,
                    B = ddf$B + 2*ddf$BB + ddf$AB)
    
    
    
    ddf = data.frame(A = 0, AA = 0, AB = 0,
                      B = 0, BB = 0)
    
    ddf[1,which(colnames(ddf) %in% colnames(backgd))] = backgd

    dff = data.frame(A = ddf$A + 2*ddf$AA + ddf$AB,
                    B = ddf$B + 2*ddf$BB + ddf$AB)
    
    
    
    df = df %>%
        tidyr::gather() %>%
        mutate(value = value / sum(value))
    
    dff = dff %>%
        tidyr::gather() %>%
        mutate(value = value / sum(value)) %>%
        dplyr::rename(back = value) %>%
        mutate(obs = df$value) %>%
        mutate(log2_enrich = log2(obs / back))
    
    display(dff)
    display(chisq.test(dff$obs, p = dff$back))
    
    comp_switch_enrichment[[type]] = dff
}

## plot
p = matrix(0L, nrow = 2, ncol = length(cell_types))

cell_types_plot = c("GM12878", "H1Esc", "HAP1", "HFF", "IMR90")

for (i in 1:length(cell_types_plot)) {
    p[,i] = comp_switch_enrichment[[cell_types_plot[i]]]$log2_enrich
}

rownames(p) = c("A", "B")
colnames(p) = cell_types_plot

paletteLength = 500
myColor <- rev(colorRampPalette(RColorBrewer::brewer.pal(7, "YlGnBu"))(paletteLength))


options(repr.plot.width=5, repr.plot.height=3)
pheatmap(p, breaks = c(seq(-0.6, 0, length.out=ceiling(paletteLength/2) + 1),
            seq(0.3/paletteLength, 0.3, length.out=floor(paletteLength/2))),
         color = myColor,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         width = 5,
         height = 3,
         filename = "comp_switch_heatmap.pdf")



##### Enrichment of compartment switching regions in cell-type specific genes
options(scipen=10)

genes_df = readRDS("de_res_combined_GM_H1_IMR_HFFc6_HAP1.rds")

ensemble_gene_table_file = '/net/trapnell/vol1/home/khj3017/proj/scripts/figures/EnsemblGeneTable.hg19.withStrand.txt'
ensemble_gene_table = read.table(ensemble_gene_table_file, skip = 1, 
                                 col.names = c("id", "chr", "start", "end", "strand", "gene_short_name", "type"), 
                                 colClasses = c("character", "character", "integer", "integer", "integer", "character", "character")) %>%
                        filter(type == "protein_coding")


# map genes to appropriate 500kb bin
get_gene_df = function(ensemble_genes_de, ensemble_gene_table) {
    df = ensemble_gene_table %>%
        filter(gene_short_name %in% ensemble_genes_de$gene_short_name) %>%
        select(chr, start, end, gene_short_name)
    
    df <- left_join(df, ensemble_genes_de, by = "gene_short_name")
    
    return(df)
}

mapped_gene_df = lapply(genes_df, get_gene_df, ensemble_gene_table)


filtered_mapped_gene_df = list()

# filter mapped gene df
for (type in names(mapped_gene_df)) {
    filtered_mapped_gene_df[[type]] = mapped_gene_df[[type]] %>%
                    mutate(chr = paste0("chr", chr),
                           midpoint = floor((start + end) / 2 / 500000) * 500000 + 250000) %>%
                    mutate(LP = paste0(chr, ":", midpoint)) %>%
                    filter(chr %in% paste0("chr", 1:22)) %>%
                    select(chr:LP, gene_short_name,starts_with('log2'), starts_with('padj')) %>%
                    arrange(-(.[[7]] + .[[8]] + .[[9]] + .[[10]]))
}

## Compute frequency of comp switching regions around upregulated
## cell-type specific DE genes

cell_types = names(filtered_mapped_gene_df)

count_anno_de = list()
genome_anno = list()

for (type in cell_types) {
    count_df = comp_switch_df[[type]] %>%
            mutate(LP = paste0(chr, ":", midpoint),
                   label = ifelse(!!sym(type) > 0, "A", "B")) %>%
            left_join(., filtered_mapped_gene_df[[type]], by = "LP") %>% na.omit()
    
    # observed
    count_anno_de[[type]] = 
            count_df %>%
            count(label) %>% 
            mutate(n = n / sum(n))
    
    # expected (genome-wide)
    genome_anno[[type]] = comp_switch_df[[type]] %>%
            mutate(LP = paste0(chr, ":", midpoint),
                   label = ifelse(!!sym(type) > 0, "A", "B")) %>%
            count(label) %>% 
            mutate(n = n / sum(n))
    
}

## compute enrichment
out = lapply(names(count_anno_de), function(type) {
    a = count_anno_de[[type]] %>% pull(n)
    b = genome_anno[[type]] %>% pull(n)
    out = genome_anno[[type]]
    
    if (length(a) == 1) {
        if (count_anno_de[[type]]$label == "A") {
            a = c(a, 0)
        } else {
            a = c(0, a)
        }
    }
    
    out$n = log2(a / b)
    out
})

names(out) = names(count_anno_de)

gene_comp_switch_enrich = NULL
for (type in names(out)) {
    gene_comp_switch_enrich = rbind(gene_comp_switch_enrich,
                   out[[type]] %>%
                      mutate(type = type))
}

infinite_filter = is.infinite(gene_comp_switch_enrich$n)
gene_comp_switch_enrich[infinite_filter,]$n = 
		ifelse(gene_comp_switch_enrich[infinite_filter,]$n > 0, 3, -4)

## plot results
options(repr.plot.width=4, repr.plot.height=4)
ggplot(gene_comp_switch_enrich, 
       aes(x = label, y = n, fill = type)) + 
    geom_bar(stat = "identity", position = "dodge", color = "black") + 
    geom_hline(yintercept = 0, color = "black") + 
    theme_minimal()
