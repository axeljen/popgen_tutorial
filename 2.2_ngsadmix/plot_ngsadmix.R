args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
    stop("Usage: Rscript plot_ngsadmix.R <ngsadmix_prefix> <K> <bamlist> <population_file>")
}

prefix <- args[1]
K <- args[2]
bamlist <- args[3]
popfile <- args[4]

suppressPackageStartupMessages(library(tidyverse))

q_file <- paste0(prefix, ".qopt")
sample_files <- readLines(bamlist)
sample_ids <- basename(sample_files)
sample_ids <- sub("\\.dedup\\.(cram|bam)$", "", sample_ids)
q <- read.table(q_file, header = FALSE, stringsAsFactors = FALSE)
pop <- read.table(
    popfile,
    header = FALSE,
    col.names = c("IID", "POP"),
    stringsAsFactors = FALSE
)

if (nrow(q) != length(sample_ids)) {
    stop("The .qopt file and BAM list contain different numbers of individuals.")
}

if (!setequal(sample_ids, pop$IID)) {
    stop("Sample IDs in the BAM list and population file do not match.")
}

q_data <- as_tibble(q)
colnames(q_data) <- paste0("K", seq_len(ncol(q_data)))
q_data$IID <- sample_ids

plot_data <- q_data %>%
    left_join(pop, by = "IID") %>%
    arrange(POP, IID) %>%
    mutate(IID = factor(IID, levels = unique(IID))) %>%
    pivot_longer(
        starts_with("K"),
        names_to = "Component",
        values_to = "Ancestry"
    )

ngsadmix_plot <- ggplot(
    plot_data,
    aes(x = IID, y = Ancestry, fill = Component)
) +
    geom_col(width = 1, color = "black", linewidth = 0.2) +
    facet_grid(~ POP, scales = "free_x", space = "free_x") +
    labs(x = "Individuals", y = "Ancestry proportion") +
    theme_minimal() +
    theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.spacing.x = unit(0.1, "lines")
    )

ggsave(
    paste0(prefix, ".qopt.png"),
    ngsadmix_plot,
    width = 12,
    height = 5,
    dpi = 300
)
