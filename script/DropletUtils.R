#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(DropletUtils)
  library(SingleCellExperiment)
})

# ================================
# Step 1: 定义命令行参数
# ================================
option_list <- list(
  make_option(c("-i", "--input"), type="character", help="Input 10X folder"),
  make_option(c("-f", "--fdr"), type="double", default=0.001,
              help="FDR threshold [default: %default]"),
  make_option(c("-l", "--lower"), type="integer", default=100,
              help="A numeric scalar specifying the lower bound on the total UMI count, at or below which all barcodes are assumed to correspond to empty droplets. [default: %default]")
)

opt_parser <- OptionParser(option_list=option_list)

# 如果没有参数，显示帮助并退出
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  print_help(opt_parser)
  quit(status = 0)
}

opt <- parse_args(opt_parser)

# 如果输入为空，显示帮助
if (is.null(opt$input)) {
  print_help(opt_parser)
  quit(status = 1)
}

# ================================
# 脚本逻辑
# ================================
input_dir <- opt$input
fdr_cutoff <- opt$fdr
lower_bound <- opt$lower

plot_file <- file.path(input_dir, "knee_plot.png")
output_dir <- file.path(input_dir, "filtered_matrix")

if (!dir.exists(input_dir)) {
  stop("输入目录不存在: ", input_dir)
}

if (dir.exists(output_dir)) {
  stop("输出文件夹已经存在，请重新定义输入目录或删除: ", output_dir)
}

message("Reading 10X data from: ", input_dir)
sce <- read10xCounts(input_dir)

br.out <- barcodeRanks(counts(sce))
knee_val <- metadata(br.out)$knee
inflection_val <- metadata(br.out)$inflection

png(plot_file, width=1600, height=1600, res=200)
plot(br.out$rank, br.out$total, log="xy",
     xlab="Rank", ylab="Total UMI Count", main="Barcode Rank Plot")
o <- order(br.out$rank)
lines(br.out$rank[o], br.out$fitted[o], col="red")
abline(h=knee_val, col="dodgerblue", lty=2)
abline(h=inflection_val, col="forestgreen", lty=2)
text(x=max(br.out$rank)*0.5, y=knee_val, labels=paste0("knee=", round(knee_val, 0)),
     col="dodgerblue", pos=3, cex=0.9)
text(x=max(br.out$rank)*0.5, y=inflection_val, labels=paste0("inflection=", round(inflection_val, 0)),
     col="forestgreen", pos=3, cex=0.9)
legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"),
       legend=c("knee", "inflection"))
dev.off()

set.seed(100)
message("Running emptyDrops with lower=", lower_bound)
e.out <- emptyDrops(counts(sce), lower=lower_bound, niters=10000)

retained <- sce[, which(e.out$FDR <= fdr_cutoff)]
message("Retained ", ncol(retained), " barcodes after filtering (FDR <= ", fdr_cutoff, ")")

write10xCounts(path=output_dir,
               x=counts(retained),
               gene.id=rowData(retained)$ID,
               gene.symbol=rowData(retained)$Symbol,
               barcode=colData(retained)$Barcode)

message("Filtered matrix saved to: ", output_dir)
