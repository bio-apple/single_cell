#!/usr/bin/env Rscript

# ----------------------------
# Check required packages
# ----------------------------
required_pkgs <- c("optparse", "DropletUtils", "SingleCellExperiment", "scDblFinder", "ggplot2")
for (p in required_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    stop("Package ", p, " is required but not installed.")
  }
}

suppressPackageStartupMessages({
  library(optparse)
  library(DropletUtils)
  library(SingleCellExperiment)
  library(scDblFinder)
  library(ggplot2)
})

# ----------------------------
# Define command-line options
# ----------------------------
option_list <- list(
  make_option(c("-i", "--input"), type="character", help="Input 10X folder (required)"),
  make_option(c("-s", "--sample"), type="character", help="Sample name (required)"),
  make_option(c("-o", "--output"), type="character", help="Output directory (required)"),
  make_option(c("-f", "--fdr"), type="double", default=0.001,
              help="FDR threshold for emptyDrops [default: %default]"),
  make_option(c("-l", "--lower"), type="integer", default=100,
              help="Lower bound on UMI count for emptyDrops [default: %default]"),
  make_option(c("-r", "--rate"), type="double", default=0.075,
              help="Expected doublet rate for scDblFinder [default: %default]")
)

opt_parser <- OptionParser(option_list=option_list)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  print_help(opt_parser)
  quit(status = 0)
}

opt <- parse_args(opt_parser)

# ----------------------------
# Check required arguments
# ----------------------------
missing_params <- c()
if (is.null(opt$input)) missing_params <- c(missing_params, "-i/--input")
if (is.null(opt$sample)) missing_params <- c(missing_params, "-s/--sample")
if (is.null(opt$output)) missing_params <- c(missing_params, "-o/--output")

if (length(missing_params) > 0) {
  message("Error: The following arguments are required: ", paste(missing_params, collapse = ", "), "\n")
  print_help(opt_parser)
  quit(status = 1)
}

# ----------------------------
# Assign variables
# ----------------------------
input_dir <- opt$input
sample_name <- opt$sample
output_dir <- opt$output
fdr_cutoff <- opt$fdr
lower_bound <- opt$lower
dbrate <- opt$rate

if (!dir.exists(input_dir)) stop("Input directory does not exist: ", input_dir)

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  message("Output directory does not exist, created: ", output_dir)
} else {
  message("Output directory exists, continuing: ", output_dir)
}

# Define output files
plot_file_knee <- file.path(output_dir, paste0(sample_name, "_knee_plot.png"))
plot_file_doublet <- file.path(output_dir, paste0(sample_name, "_doublet_score.png"))
final_matrix_dir <- file.path(output_dir, paste0(sample_name, "_final_matrix"))
qc_file <- file.path(output_dir, paste0(sample_name, "_QC_summary.csv"))

# ----------------------------
# Step 1: Read input data
# ----------------------------
message("Reading 10X data from: ", input_dir)
sce <- read10xCounts(input_dir)

# ----------------------------
# Step 2: Barcode rank plot
# ----------------------------
br.out <- barcodeRanks(counts(sce))
knee_val <- metadata(br.out)$knee
inflection_val <- metadata(br.out)$inflection

if (file.exists(plot_file_knee)) message("Overwriting existing file: ", plot_file_knee)
png(plot_file_knee, width=1600, height=1600, res=200)
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

# ----------------------------
# Step 3: emptyDrops filtering
# ----------------------------
set.seed(100)
message("Running emptyDrops with lower=", lower_bound)
e.out <- emptyDrops(counts(sce), lower=lower_bound, niters=10000)

sce_filtered <- sce[, which(e.out$FDR <= fdr_cutoff)]
if (ncol(sce_filtered) == 0) {
  stop("No cells passed emptyDrops filtering. Check your input or adjust parameters.")
}
message("After emptyDrops: ", ncol(sce_filtered), " barcodes retained")

# ----------------------------
# Step 4: scDblFinder for doublets
# ----------------------------
sce_filtered <- scDblFinder(sce_filtered, dbr=dbrate)
saveRDS(sce_filtered, file=file.path(output_dir, paste0(sample_name, "_sce_filtered.rds")))

dbl_scores <- colData(sce_filtered)$scDblFinder.score
dbl_class <- colData(sce_filtered)$scDblFinder.class
df_plot <- data.frame(Score=dbl_scores, Class=dbl_class)

if (file.exists(plot_file_doublet)) message("Overwriting existing file: ", plot_file_doublet)
p <- ggplot(df_plot, aes(x=Score, fill=Class)) +
  geom_histogram(alpha=0.6, bins=50, position="identity") +
  scale_fill_manual(values=c("singlet"="steelblue", "doublet"="firebrick")) +
  labs(title=paste("Doublet Score Distribution -", sample_name),
       x="Doublet Score", y="Cell Count", fill="Class") +
  theme_minimal(base_size = 14)
ggsave(filename=plot_file_doublet, plot=p, width=8, height=6, dpi=300)

# ----------------------------
# Step 5: Filter singlets & export
# ----------------------------
sce_final <- sce_filtered[, dbl_class == "singlet"]
message("Final singlets: ", ncol(sce_final), " / ", ncol(sce_filtered))

# If final matrix directory exists, delete it
if (dir.exists(final_matrix_dir)) {
  message("Output final matrix directory already exists, removing: ", final_matrix_dir)
  unlink(final_matrix_dir, recursive = TRUE)
}

write10xCounts(path=final_matrix_dir,
               x=counts(sce_final),
               gene.id=rowData(sce_final)$ID,
               gene.symbol=rowData(sce_final)$Symbol,
               barcode=colData(sce_final)$Barcode)

# ----------------------------
# Step 6: Write summary
# ----------------------------
total_droplets <- ncol(sce)
after_emptyDrops <- ncol(sce_filtered)
dbl_count <- sum(dbl_class == "doublet")
final_cells <- ncol(sce_final)

summary_df <- data.frame(
  Sample=sample_name,
  Droplets=total_droplets,
  After_emptyDrops=after_emptyDrops,
  knee=round(knee_val, 0),
  inflection=round(inflection_val, 0),
  lower=lower_bound,
  FDR=fdr_cutoff,
  Doublets=dbl_count,
  Final_Singlets=final_cells
)

if (file.exists(qc_file)) message("Overwriting existing QC summary file: ", qc_file)
write.csv(summary_df, file=qc_file, row.names=FALSE)

message("QC completed. Results saved to: ", output_dir)
