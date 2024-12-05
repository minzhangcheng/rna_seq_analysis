library(dplyr)
library(readr)
library(limma)
library(edgeR)


dispersion=0.1
work_dir <- "/home/minzhang/workspace/RNA-Seq"
expr_dir <- file.path(work_dir, "expr")
deg_dir <- file.path(work_dir, "deg")

if(!dir.exists(deg_dir)){
  dir.create(deg_dir)
}

expr_files <- list.files(expr_dir, pattern = "*.genes.results$", full.names = TRUE)
sample_names <- sapply(basename(expr_files), function(x) sub(".genes.results$", "", x))

expr_list <- lapply(expr_files, function(file){
  df <- read_table(file, col_types = cols())
  df %>% select(gene_id, expected_count)
})
count_matrix <- Reduce(function(x, y) full_join(x, y, by = "gene_id"), expr_list)
count_matrix[is.na(count_matrix)] <- 0
count_data <- as.data.frame(count_matrix)
rownames(count_data) <- count_data$gene_id
count_data <- count_data[ , -which(names(count_data) == "gene_id")]
colnames(count_data) <- sample_names

samples <- colnames(count_data)
pairwise_combinations <- combn(samples, 2, simplify = FALSE)

# 定义差异表达分析函数
perform_de_analysis <- function(pair){
  group <- factor(c(pair[1], pair[2]))
  counts = count_data[, pair]
  dge <- DGEList(counts = counts, group = group)
  dge$common.dispersion <- dispersion
  et <- exactTest(dge)
  result <- topTags(et, n = nrow(dge))
  result_table <- result$table
  result_table <- cbind(gene_id = rownames(result_table), result_table)
  out_file <- file.path(deg_dir, paste0(pair[1], "_vs_", pair[2], "_DEG.csv"))
  write_csv(result_table, out_file)
}

# 执行所有两两比较
for(pair in pairwise_combinations){
  perform_de_analysis(pair)
}

