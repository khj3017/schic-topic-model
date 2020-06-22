#############################################
# Run cisTopic
#############################################

library(cisTopic)
library(data.table)
library(Matrix)

args <- commandArgs(trailingOnly = TRUE)
mat_file <- as.character(args[1]) # cell-LP matrix in sparse matrix format
lp_label_file <- as.character(args[2]) # LP annotations (i.e. LP names)
cell_file <- as.character(args[3]) # cell annotations (i.e. cell type and replicate info)
out_dir <- as.character(args[4]) # out directory
low_b <- as.double(args[5]) # lowest topic number to try
up_b <- as.double(args[6]) # highest topic number to try
increment <- as.double(args[7]) # increment for topics
resolution <- as.double(args[8]) # resolution of contact matrix for name purposes

options(scipen = 10) 
topic_list = seq(low_b,up_b, by = increment)
outname = paste0("cisTopic_", resolution)

outDir <- file.path(out_dir, outname, paste(low_b, "-", up_b, sep=""))

if (!file.exists(file.path(out_dir, outname))) {
  dir.create(file.path(out_dir, outname))
}
if (!file.exists(outDir)) {
  dir.create(outDir)
}

## load data
ptm = proc.time()

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
    col.names = c("cell", "cell.type", "replicate"),
    colClasses = c("character", "factor", "replicate"))

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



###
cisTopicObject <- createcisTopicObject(mat, project.name="full_cisTopic", keepCountsMatrix=FALSE)
cisTopicObject <- addCellMetadata(cisTopicObject, cell.data = cell_list_ordered)
rm(mat)


PDFfilename = file.path(outDir, paste0("cisTopic_Reports_", low_b, "-", up_b, ".pdf"))
pdf(PDFfilename)

# runModel
cisTopicObject <- runCGSModels(cisTopicObject, topic=topic_list, seed=999, nCores=length(topic_list))
print(proc.time() - ptm)

saveRDS(cisTopicObject, file=file.path(outDir, paste0("cisTopicObject_CGS_", resolution, "_", low_b, "-", up_b, ".rds")))
dev.off()

