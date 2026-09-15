library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

prefix <- args[1]

if (!length(args) > 2) {
    if (length(args) == 2){
        popfile <- args[2]
    } 
} else {
    stop("Too many arguments provided")
}

# find eigenval and eigenvecs
eigenval_file <- paste0(prefix, ".eigenval")
eigenvec_file <- paste0(prefix, ".eigenvec")

# parse input files
eigenval <- scan(eigenval_file)
eigenvec <- read.table(eigenvec_file, header = TRUE, comment.char = "", check.names = FALSE)
colnames(eigenvec)[1] <- "IID"
variance_percent <- 100 * eigenval / sum(eigenval)
pc1_label <- paste0("PC1 (", round(variance_percent[1], 1), "%)")
pc2_label <- paste0("PC2 (", round(variance_percent[2], 1), "%)")

if (exists("popfile")) {
    pop <- read.table(popfile, header = FALSE)
}

# prepare data for plotting
pca_data <- eigenvec
if (exists("pop")) {
    colnames(pop) <- c("IID", "POP")
    if (!setequal(pca_data$IID, pop$IID)) {
        stop("Sample IDs in the PLINK eigenvector file and population file do not match.")
    }
    pca_data <- merge(pca_data, pop, by = c("IID")) %>%
        # arrange by population and make factor
        arrange(POP) %>%
        mutate(POP = factor(POP, levels = unique(POP)))
}



# plot PCA

# if popfile, include pops with color:
if (exists("popfile")) {
    p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = POP)) +
        geom_point(size = 3) +
        labs(x = pc1_label, y = pc2_label, color = "Population") +
        theme_minimal()
} else {
    p <- ggplot(pca_data, aes(x = PC1, y = PC2)) +
        geom_point(size = 3) +
        labs(x = pc1_label, y = pc2_label) +
        theme_minimal()
}

# include eigenvalues as an inset bar plot 
eigenval_data <- data.frame(PC = paste0("PC", 1:length(eigenval)), Variance = eigenval)

# ensure that the PCs are ordered correctly
eigenval_data$PC <- factor(eigenval_data$PC, levels = paste0("PC", 1:length(eigenval)))

bar_plot <- ggplot(eigenval_data, aes(x = PC, y = Variance)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Principal Component", y = "Eigenvalue")

# save the bar plot separately
ggsave(filename = paste0(prefix, "_eigenval_bar_plot.png"), plot = bar_plot)


ggsave(filename = paste0(prefix, "_pca_plot.png"), plot = p, width = 8, height = 8)

