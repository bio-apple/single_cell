import os
import re
import subprocess
import tempfile
import argparse

# R脚本内容，直接复用之前的测序饱和度分析代码
r_script_content = '''
library(Matrix)
library(ggplot2)
library(gridExtra)
library(Rsamtools)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
raw_prefix <- args[1]
filtered_prefix <- args[2]
bam_path <- args[3]
output <- ifelse(length(args) >= 4, args[4], "summary_with_cells.png")

load_matrix <- function(prefix) {
  mtx <- readMM(paste0(prefix, ".matrix.mtx.gz"))
  features <- read.delim(paste0(prefix, ".features.tsv.gz"), header = FALSE)
  barcodes <- read.delim(paste0(prefix, ".barcodes.tsv.gz"), header = FALSE)
  rownames(mtx) <- features$V1
  colnames(mtx) <- barcodes$V1
  return(mtx)
}

downsample_matrix <- function(mat, depth_ratio) {
  mat@x <- rbinom(length(mat@x), mat@x, prob = depth_ratio)
  return(mat)
}

get_bam_coverage <- function(bam_path, barcodes) {
  cat("Scanning BAM for CB tag coverage...\\n")
  param <- ScanBamParam(tag="CB", what="qname")
  bam <- scanBam(bam_path, param=param)[[1]]
  cb_tags <- bam$tag$CB
  cb_tags <- cb_tags[!is.na(cb_tags)]
  tab <- data.table(table(cb_tags))
  setnames(tab, c("barcode", "count"))
  tab <- tab[barcode %in% barcodes]
  return(tab)
}

plot_saturation_summary <- function(raw_prefix, filtered_prefix, bam_path = NULL, output = "summary_with_cells.png") {
  raw_mtx <- load_matrix(raw_prefix)
  filtered_mtx <- load_matrix(filtered_prefix)

  barcodes <- intersect(colnames(raw_mtx), colnames(filtered_mtx))
  raw_mtx <- raw_mtx[, barcodes]

  depths <- seq(0.05, 1, 0.1)
  results <- data.frame()

  for (d in depths) {
    mtx_d <- downsample_matrix(raw_mtx, d)

    umi_per_cell <- Matrix::colSums(mtx_d)
    mean_reads <- mean(umi_per_cell)
    median_umi <- median(umi_per_cell)

    gene_counts <- Matrix::colSums(mtx_d > 0)
    median_genes <- median(gene_counts)
    cells_detected <- sum(gene_counts >= 200)

    total_genes <- sum(Matrix::rowSums(mtx_d > 0) > 0)

    umi_triplets <- sum(umi_per_cell)
    unique_triplets <- Matrix::nnzero(mtx_d)
    saturation <- ifelse(umi_triplets == 0, 0, 1 - (unique_triplets / umi_triplets))

    results <- rbind(results, data.frame(
      depth_ratio = d,
      mean_reads_per_cell = mean_reads,
      total_genes_detected = total_genes,
      median_genes = median_genes,
      median_umis = median_umi,
      saturation = saturation * 100,
      cells_detected = cells_detected
    ))
  }

  p5 <- NULL
  if (!is.null(bam_path) && file.exists(bam_path)) {
    bam_cov <- get_bam_coverage(bam_path, colnames(raw_mtx))
    p5 <- ggplot(bam_cov, aes(x = count)) +
      geom_histogram(bins = 50, fill = "#1f78b4") +
      labs(title = "Reads per Cell from BAM", x = "Reads per Cell", y = "Cell Count") +
      theme_bw()
  }

  p1 <- ggplot(results, aes(mean_reads_per_cell, saturation)) +
    geom_point() + geom_line() +
    labs(title = "Sequencing Saturation", x = "Mean Reads per Cell", y = "Percent") +
    theme_bw()

  p2 <- ggplot(results, aes(mean_reads_per_cell, total_genes_detected)) +
    geom_point() + geom_line() +
    labs(title = "Total Genes Detected", x = "Mean Reads per Cell", y = "Genes") +
    theme_bw()

  p3 <- ggplot(results, aes(mean_reads_per_cell, median_genes)) +
    geom_point() + geom_line() +
    labs(title = "Median Genes per Cell", x = "Mean Reads per Cell", y = "Genes") +
    theme_bw()

  p4 <- ggplot(results, aes(mean_reads_per_cell, median_umis)) +
    geom_point() + geom_line() +
    labs(title = "Median UMI Counts per Cell", x = "Mean Reads per Cell", y = "UMI Counts") +
    theme_bw()

  p6 <- ggplot(results, aes(mean_reads_per_cell, cells_detected)) +
    geom_point(color = "#d95f02") + geom_line(color = "#d95f02") +
    labs(title = "Cells Detected (≥200 genes)", x = "Mean Reads per Cell", y = "Number of Cells") +
    theme_bw()

  if (!is.null(p5)) {
    ggsave(output, arrangeGrob(p1, p2, p3, p4, p5, p6, ncol = 2), width = 14, height = 12)
  } else {
    ggsave(output, arrangeGrob(p1, p2, p3, p4, p6, ncol = 2), width = 12, height = 10)
  }
}

plot_saturation_summary(raw_prefix, filtered_prefix, bam_path, output)
'''

def find_prefixes(input_dir):
    # 查找所有matrix.mtx.gz文件，根据文件名前缀分组
    files = os.listdir(input_dir)
    matrix_files = [f for f in files if f.endswith(".matrix.mtx.gz")]
    prefixes = set()

    for f in matrix_files:
        # 移除后缀得到前缀
        prefix = f[:-len(".matrix.mtx.gz")]
        prefixes.add(prefix)
    return list(prefixes)

def select_raw_filtered_prefixes(prefixes):
    # 按照含filtered关键字区分原始和过滤
    raw_prefix = None
    filtered_prefix = None
    for p in prefixes:
        if "filtered" in p:
            filtered_prefix = p
        else:
            raw_prefix = p
    if raw_prefix is None or filtered_prefix is None:
        raise ValueError("找不到原始和过滤矩阵前缀，请确认文件命名包含filtered区分")
    return raw_prefix, filtered_prefix

def find_bam_file(input_dir, prefix):
    # 试着找跟prefix匹配的bam文件
    # 例如 prefix=PIP-seq-lib-1-Ctrl.scRNA.filtered
    # 先去掉.filtered再找bam
    base = prefix.replace(".filtered", "")
    files = os.listdir(input_dir)
    for f in files:
        if f.endswith(".bam") and base in f:
            return os.path.join(input_dir, f)
    # 找不到就报错
    raise FileNotFoundError(f"找不到匹配 BAM 文件，基于前缀 {base}")

def run_saturation(input_dir, output):
    prefixes = find_prefixes(input_dir)
    raw_prefix, filtered_prefix = select_raw_filtered_prefixes(prefixes)

    bam_path = find_bam_file(input_dir, filtered_prefix)

    raw_prefix_path = os.path.join(input_dir, raw_prefix)
    filtered_prefix_path = os.path.join(input_dir, filtered_prefix)

    with tempfile.NamedTemporaryFile("w", suffix=".R", delete=False) as tmp_rfile:
        tmp_rfile.write(r_script_content)
        r_script_path = tmp_rfile.name

    try:
        cmd = [
            "Rscript",
            r_script_path,
            raw_prefix_path,
            filtered_prefix_path,
            bam_path,
            output
        ]
        print("运行R脚本做测序饱和度和细胞数稀疏曲线分析...")
        subprocess.run(cmd, check=True)
        print(f"分析完成，结果保存在 {output}")
    finally:
        os.remove(r_script_path)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="根据文件名前缀调用R做测序饱和度分析")
    parser.add_argument("-i", "--input_dir", required=True, help="包含matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz和bam的文件夹")
    parser.add_argument("-o", "--output", default="summary_with_cells.png", help="输出PNG文件路径")
    args = parser.parse_args()

    if not os.path.isdir(args.input_dir):
        raise FileNotFoundError(f"输入目录不存在: {args.input_dir}")

    run_saturation(args.input_dir, args.output)
