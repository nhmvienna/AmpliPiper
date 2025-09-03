# Load required libraries
library(tidyverse)

# Read from command line
args <- commandArgs(trailingOnly = TRUE)
input <- args[1]

# Read input file
DATA <- read.csv(input, header = TRUE, na.string = "NA")

# Function to generate an extended color palette
extend_palette <- function(n) {
    base_colors <- c("#7ac1dc", "#ff8e3e", "#21614a", "#C3E6DA") # Main 4 colors
    if (n <= length(base_colors)) {
        return(base_colors[1:n])
    }
    # Dynamically generate additional colors if needed
    extra_colors <- scales::hue_pal()(n - length(base_colors))
    return(c(base_colors, extra_colors))
}

# Calculate dimensions with a minimum value
min_width <- 5 # Minimum width for single-locus datasets
min_height <- 5 # Minimum height for single-sample datasets
width <- max(length(levels(as.factor(DATA$Locus))), min_width)
height <- max(length(levels(as.factor(DATA$ID))) / 4, min_height)

# Create a palette based on the number of unique NoAlleles
no_alleles_palette <- extend_palette(length(unique(DATA$NoAlleles)))

# Make heatmap for number of haplotypes
PLOT <- ggplot(DATA, aes(y = ID, x = Locus, fill = as.factor(NoAlleles))) +
    geom_tile(color = "black") +
    theme_bw() +
    ylab("Samples") +
    xlab("Loci") +
    scale_fill_manual(
        name = "Number of \nhaplotypes",
        values = no_alleles_palette
    ) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# Save file as PNG
OUTPUT <- paste0(input, ".png")
ggsave(PLOT, file = OUTPUT, width = width, height = height)

OUTPUT <- paste0(input, ".pdf")
ggsave(PLOT, file = OUTPUT, width = width, height = height)

# Create a palette based on the number of unique ExpectedPloidy
expected_ploidy_palette <- extend_palette(length(unique(DATA$ExpectedPloidy)))

# Make heatmap for expected ploidy
PLOT <- ggplot(DATA, aes(y = ID, x = Locus, fill = as.factor(ExpectedPloidy))) +
    geom_tile(color = "black") +
    theme_bw() +
    ylab("Samples") +
    xlab("Loci") +
    scale_fill_manual(
        name = "Modelled\nPloidy",
        values = expected_ploidy_palette
    ) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# Save file as PNG
OUTPUT <- paste0(input, "_ExpectedPloidy.png")
ggsave(PLOT, file = OUTPUT, width = width, height = height)

OUTPUT <- paste0(input, "_ExpectedPloidy.pdf")
ggsave(PLOT, file = OUTPUT, width = width, height = height)
