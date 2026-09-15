library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1 || length(args) > 3) {
    stop("Usage: Rscript plot_pca_angsd.R <pcangsd_prefix> [population_file] [bamlist]")
}

prefix <- args[1]
popfile <- if (length(args) >= 2) args[2] else NULL
bamlist <- if (length(args) == 3) args[3] else NULL

covariance_file <- paste0(prefix, ".cov")
if (!file.exists(covariance_file)) {
    stop("Could not find PCAngsd covariance matrix: ", covariance_file)
}

covariance <- as.matrix(read.table(covariance_file))
pca <- eigen(covariance, symmetric = TRUE)
eigenval <- pca$values
variance_percent <- 100 * eigenval / sum(eigenval)
pc1_label <- paste0("PC1 (", round(variance_percent[1], 1), "%)")
pc2_label <- paste0("PC2 (", round(variance_percent[2], 1), "%)")
pc_coordinates <- sweep(pca$vectors, 2, sqrt(pmax(eigenval, 0)), "*")

if (is.null(bamlist)) {
    stop("A BAM/CRAM list is required to associate PCAngsd rows with sample IDs.")
}

sample_files <- readLines(bamlist)
sample_ids <- basename(sample_files)
sample_ids <- sub("\\.dedup\\.(cram|bam)$", "", sample_ids)

if (length(sample_ids) != nrow(pc_coordinates)) {
    stop("The number of samples in the BAM/CRAM list does not match the number of PCAngsd rows.")
}

pca_data <- as.data.frame(pc_coordinates)
colnames(pca_data) <- paste0("PC", seq_len(ncol(pca_data)))
pca_data$IID <- sample_ids

if (!is.null(popfile)) {
    pop <- read.table(popfile, header = FALSE, stringsAsFactors = FALSE)
    if (ncol(pop) < 2) {
        stop("The population file must contain sample ID and population columns.")
    }
    pop <- pop[, 1:2]
    colnames(pop) <- c("IID", "POP")
    if (!setequal(pca_data$IID, pop$IID)) {
        stop("Sample IDs in the ANGSD input list and population file do not match.")
    }
    pca_data <- merge(pca_data, pop, by = "IID", sort = FALSE)
}

if (!is.null(popfile)) {
    p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = as.factor(POP))) +
        geom_point(size = 3) +
        labs(x = pc1_label, y = pc2_label, color = "Population") +
        theme_minimal()
} else {
    p <- ggplot(pca_data, aes(x = PC1, y = PC2)) +
        geom_point(size = 3) +
        labs(x = pc1_label, y = pc2_label) +
        theme_minimal()
}

ggsave(filename = paste0(prefix, "_pca_plot.png"), plot = p, width = 8, height = 8)

eigenval_data <- data.frame(
    PC = factor(paste0("PC", seq_along(eigenval)), levels = paste0("PC", seq_along(eigenval))),
    Variance = eigenval
)

bar_plot <- ggplot(eigenval_data, aes(x = PC, y = Variance)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Principal Component", y = "Eigenvalue")

ggsave(filename = paste0(prefix, "_eigenval_bar_plot.png"), plot = bar_plot)
