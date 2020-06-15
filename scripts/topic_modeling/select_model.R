suppressPackageStartupMessages(suppressWarnings({
    library(data.table)
    library(Matrix)
    library(cisTopic)
}))

args <- commandArgs(trailingOnly = TRUE)
cisTopicFile = as.character(args[1])
mat_file = as.character(args[2])
lp_label_file = as.character(args[3])
cell_file = as.character(args[4])
out_dir = as.character(args[5])

## load data


df = fread(
    mat_file,
    col.names = c("cell.idx", "lp.idx", "count"),
    colClasses = c("integer", "integer", "integer"))

lp.annotations = read.table(
    lp_label_file,
    col.names = c("lp"),
    colClasses = c("character"))

cell.annotations = read.table(
    cell_file,
    col.names = c("cell", "cell.type"),
    colClasses = c("character", "factor"))

rownames(lp.annotations) = lp.annotations$lp
rownames(cell.annotations) = cell.annotations$cell

df$lp.idx = df$lp.idx + 1

# add a dummy cell to ensure that all genes are included in the matrix
# even if a gene isn't expressed in any cell
df = rbind(df, data.frame(
    cell.idx = rep(nrow(cell.annotations)+1, 2),
    lp.idx =  c(1, nrow(lp.annotations)),
    count = c(1, 1)))

## make sparse matrix
mat = sparseMatrix(i = df$lp.idx, j = df$cell.idx, x = df$count)
mat = mat[, 1:(ncol(mat)-1)]

rownames(mat) = lp.annotations$lp
colnames(mat) = cell.annotations$cell

cell_list = cell.annotations$cell.type
idx = order(cell_list)
mat <- mat[,idx]
cell_list <- as.matrix(cell_list[idx])
rownames(cell_list) <- colnames(mat)
colnames(cell_list) = "Cell type"

cell_list_ordered <- cell_list
LIBIDs <- c("H1Esc", "HFF", "GM12878", "IMR90", "HAP1")
for (LIBID in LIBIDs) {
    cell_list_ordered[grepl(LIBID, cell_list)] <- LIBID
}

cell_list_ordered <- as.matrix(cell_list_ordered)
rownames(cell_list_ordered) <- colnames(mat)

dim(mat)

# Adapted from ldatuning
JS_divergence <- function(models, beta) {
  metrics <- sapply(models, function(model) {
    # topic-word matrix
    m1 = t(apply(model$topics, 1, function(x) {(x + beta)/sum(x + beta)}))
    # prevent NaN
    if (any(m1 == 0)) { m1 <- m1 + .Machine$double.xmin }
    # pair-wise Jensen-Shannon divergence
    pairs  <- utils::combn(nrow(m1), 2)
    jsd <- apply(pairs, 2, function(pair) {
      x <- m1[pair[1], ]
      y <- m1[pair[2], ]
      jsd <- 0.5 * sum(x*log(x/y)) + 0.5 * sum(y*log(y/x))
      return(jsd)
    })

    topic = nrow(model$topics)
    metric <- sum(jsd) / (topic*(topic-1)/2)
    return(metric)
  })
  return(metrics)
}

# from cisTopic
get_log_like = function(models) {
    loglikelihood <- sapply(seq_along(models), FUN=function(i) models[[i]]$log.likelihood[2,ncol(models[[i]]$log.likelihood)])
    # Topics
    topics <-  sapply(seq_along(models), FUN=function(i) nrow(models[[i]]$topics))
    object.log.lik <- data.frame(topics=topics, LL=loglikelihood)
    object.log.lik <- object.log.lik[order(object.log.lik$topics),]
    ll <- object.log.lik$LL
    return(ll)
}


cisTopicObject = readRDS(cisTopicFile)
models = cisTopicObject@models
beta = cisTopicObject@calc.params[['runCGSModels']]$beta
alpha = cisTopicObject@calc.params[['runCGSModels']]$alpha
dtm = t(mat[rownames(cisTopicObject@region.data),])

JS = JS_divergence(models, beta = beta)
LL = get_log_like(models)


out_df = data.frame(JS = JS, LL = LL)

outfilename = "JS_LL.rds"
saveRDS(out_df, file.path(out_dir, outfilename))
