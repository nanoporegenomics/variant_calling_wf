suppressMessages(library(dplyr))
suppressMessages(library(data.table))
library("optparse")

options(warn = -1)  # Suppress all warnings

# Define the option parser
parser <- OptionParser()
parser <- add_option(parser, c("-i", "--input"), type="character", help="Input methyl TSV from bedtools map")
parser <- add_option(parser, c("-m", "--mincpgs"), type="integer", default=10, help="Minimum number of CpGs required for regional methylation calculation [default %default]")
parser <- add_option(parser, c("-o", "--output_prefix"), type="character", default="sample", help="Output prefix for all output files (methyl_stats.bed, qc_stats.tsv and methyl_plots.pdf) [default %default]")

# Parse command-line arguments
args <- parse_args(parser)

# Print parsed arguments (for debugging)
print(args)

# Check if input file is provided
if (is.null(args$input)) {
  stop("Error: Input file is required. Use -i or --input to specify a file.", call. = FALSE)
}

# Read the TSV input file
bedmap.df <- fread(args$input, sep = "\t")
ncols=ncol(bedmap.df)
if (ncols == 11) {
  colnames(bedmap.df) <- c("chrom", "start", "end", "region_name", "n_total", "n_pass", "mean_methyl", "mean_cov", "pass_meth_values", "pass_cov_values", "pass_pos_values")
} else {
  colnames(bedmap.df) <- c("chrom", "start", "end", "region_name", paste0("col_", 5:(ncols-7)), "n_total", "n_pass", "mean_methyl", "mean_cov", "pass_meth_values", "pass_cov_values", "pass_pos_values")
}

# Compute standard deviations for methylation and coverage values
bedmap.df <- bedmap.df %>%
  mutate(
    mean_methyl = as.numeric(mean_methyl),
    mean_cov = as.numeric(mean_cov),

    # Convert pass_meth_values from comma-separated string to numeric list
    pass_meth_values = strsplit(pass_meth_values, ","),
    pass_meth_values = lapply(pass_meth_values, as.numeric),
    stdev_meth = sapply(pass_meth_values, function(x) sd(x, na.rm = TRUE)),

    # Convert pass_cov_values from comma-separated string to integer list
    pass_cov_values = strsplit(pass_cov_values, ","),
    pass_cov_values = lapply(pass_cov_values, as.integer),
    stdev_cov = sapply(pass_cov_values, function(x) sd(x, na.rm = TRUE)),

    # Convert pass_pos_values from comma-separated string to integer list
    pass_pos_values = strsplit(pass_pos_values, ","),
    pass_pos_values = lapply(pass_pos_values, as.integer),

    # Convert list columns back to comma-separated strings
    pass_meth_values = sapply(pass_meth_values, function(x) paste(x, collapse = ",")),
    pass_cov_values = sapply(pass_cov_values, function(x) paste(x, collapse = ",")),
    pass_pos_values = sapply(pass_pos_values, function(x) paste(x, collapse = ","))
  )
##### Save as output TSV
filename = paste0(args$output_prefix,".methyl_stats.bed")
write.table(bedmap.df, file = filename, sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)

##### Compute qc stats
qc.df <- bedmap.df %>%
  summarize(
    regions.noCpgs = sum(is.na(mean_methyl)),  # Count regions where mean_methyl is NA
    regions.lessCpgs = sum(n_pass < args$mincpgs, na.rm = TRUE),  # Count regions where n_pass < mincpgs
    regions.zero.stdevMeth = sum(stdev_meth == 0, na.rm = TRUE), # Count regions where methylation stdev is 0 (We can trust the mean methylation there)
    regions.zero.stdevCov = sum(stdev_cov == 0, na.rm = TRUE), # Count regions where coverage stdev is 0 (We can trust the mean coverage there)
  )
filename = paste0(args$output_prefix,".qc_stats.tsv")
write.table(qc.df, file = filename, sep = "\t", row.names = FALSE, quote = FALSE)

###### Make plots using Base R

# Open a PDF file to save multiple plots
filename = paste0(args$output_prefix, ".methyl_plots.pdf")
pdf(filename, width = 8, height = 6)

# Plot 1: Distribution of regional methylation
hist(bedmap.df$mean_methyl, 
     breaks = 20, 
     col = "#a7bec7", 
     border = "black", 
     main = "Distribution of regional methylation", 
     xlab = "Mean methylation", 
     ylab = "Number of regions")

# Plot 2: Distribution of mean coverage
hist(bedmap.df$mean_cov, 
     breaks = 20, 
     col = "#a7bec7", 
     border = "black", 
     main = "Distribution of mean coverage", 
     xlab = "Mean coverage", 
     ylab = "Number of regions")

# Plot 3: Distribution of methylation standard deviation
hist(bedmap.df$stdev_meth, 
     breaks = 20, 
     col = "#a7bec7", 
     border = "black", 
     main = "Distribution of methylation standard deviation", 
     xlab = "Methylation standard deviation", 
     ylab = "Number of regions")

# Plot 4: Distribution of coverage standard deviation
hist(bedmap.df$stdev_cov, 
     breaks = 20, 
     col = "#a7bec7", 
     border = "black", 
     main = "Distribution of coverage standard deviation", 
     xlab = "Coverage standard deviation", 
     ylab = "Number of regions")

# Plot 5: Mean methylation vs methylation standard deviation (Scatter plot)
plot(bedmap.df$mean_methyl, bedmap.df$stdev_meth, 
     col = "black", 
     pch = 16, 
     main = "Mean methylation vs Methylation SD", 
     xlab = "Mean regional methylation", 
     ylab = "Methylation standard deviation", 
     cex = 0.7)

# Plot 6: Mean coverage vs coverage standard deviation (Scatter plot)
plot(bedmap.df$mean_cov, bedmap.df$stdev_cov, 
     col = "black", 
     pch = 16, 
     main = "Mean coverage vs Coverage SD", 
     xlab = "Mean coverage", 
     ylab = "Coverage standard deviation", 
     cex = 0.7)

# Close the PDF file
dev.off()

cat("All plots saved in '", filename, "'\n")


# ##### Make plots
# # Open a PDF file to save multiple plots
# filename = paste0(args$output_prefix,".methyl_plots.pdf")
# pdf(filename, width = 8, height = 6)

# # Plot 1: Distribution of regional methylation
# ggplot(bedmap.df, aes(x = mean_methyl)) +
#   geom_histogram(binwidth = 5, fill = "#a7bec7", color = "black") +
#   theme_bw() +
#   labs(title = "Distribution of regional methylation", x = "mean methylation", y = "number of regions")

# # Plot 2: Distribution of mean coverage
# ggplot(bedmap.df, aes(x = mean_cov)) +
#   geom_histogram(binwidth = 2, fill = "#a7bec7", color = "black") +
#   theme_bw() +
#   labs(title = "Distribution of mean coverage", x = "mean coverage", y = "number of regions")

# # Plot 3: Distribution of methylation standard deviation
# ggplot(bedmap.df, aes(x = stdev_meth)) +
#   geom_histogram(binwidth = 2, fill = "#a7bec7", color = "black") +
#   theme_bw() +
#   labs(title = "Distribution of methylation standard deviation", x = "methylation standard deviation", y = "number of regions")

# # Plot 4: Distribution of coverage standard deviation
# ggplot(bedmap.df, aes(x = stdev_cov)) +
#   geom_histogram(binwidth = 1, fill = "#a7bec7", color = "black") +
#   theme_bw() +
#   labs(title = "Distribution of coverage standard deviation", x = "coverage standard deviation", y = "number of regions")

# # Plot 5: Mean methylation vs methylation standard deviation
# ggplot(bedmap.df, aes(x = mean_methyl, y = stdev_meth)) +
#   geom_point(color = "black", alpha=0.5, size=1.5) +
#   theme_bw() +
#   labs(title = "Mean methylation vs methylation standard deviation", x = "mean regional methylation", y = "methylation standard deviation")

# # Plot 6: Mean coverage vs coverage standard deviation
# ggplot(bedmap.df, aes(x = mean_cov, y = stdev_cov)) +
#   geom_point(color = "black", alpha=0.5, size=1.5) +
#   theme_bw() +
#   labs(title = "Mean coverage vs coverage standard deviation", x = "mean coverage", y = "coverage standard deviation")

# # Close the PDF file
# dev.off()

# cat("All plots saved in 'methyl_plots.pdf'\n")