args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
    stop("Usage: Rscript plot_admixture.R <admixture_prefix> <K> <population_file>")
}

prefix <- args[1]
K <- args[2]
popfile <- args[3]

suppressPackageStartupMessages(library(tidyverse))

fam_file <- paste0(prefix, ".fam")
q_file <- paste0(prefix, ".", K, ".Q")

fam <- read.table(fam_file, header = FALSE, stringsAsFactors = FALSE)
q <- read.table(q_file, header = FALSE, stringsAsFactors = FALSE)
pop <- read.table(
    popfile,
    header = FALSE,
    col.names = c("IID", "POP"),
    stringsAsFactors = FALSE
)

if (nrow(fam) != nrow(q)) {
    stop("The .fam and .Q files contain different numbers of individuals.")
}

if (!setequal(fam$V2, pop$IID)) {
    stop("Sample IDs in the .fam and population files do not match.")
}

q_data <- as_tibble(q)
colnames(q_data) <- paste0("K", seq_len(ncol(q_data)))
q_data$IID <- fam$V2

plot_data <- q_data %>%
    left_join(pop, by = "IID") %>%
    arrange(POP, IID) %>%
    mutate(IID = factor(IID, levels = unique(IID))) %>%
    pivot_longer(
        starts_with("K"),
        names_to = "Component",
        values_to = "Ancestry"
    )

admixture_plot <- ggplot(
    plot_data,
    aes(x = IID, y = Ancestry, fill = Component)
) +
    geom_col(width = 1, color = "black", linewidth = 0.2) +
    facet_grid(~ POP, scales = "free_x", space = "free_x") +
    labs(x = "Individuals", y = "Ancestry proportion", fill = "Component") +
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
    paste0(prefix, ".", K, ".Q.png"),
    admixture_plot,
    width = 12,
    height = 5,
    dpi = 300
)
