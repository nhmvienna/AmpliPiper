# Load necessary R libraries
library("ggtree")
library("tidyverse")
library(phytools) # For midpoint rooting and tree height calculation
library(treeio)

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

input <- args[1]
output <- args[2]
title <- args[3]
offset <- as.numeric(args[4])
width <- as.numeric(args[5])
height <- as.numeric(args[6])
outgroup <- args[7]
watermark <- args[8]

# Read tree
tree <- read.tree(input)

# Root the tree if an outgroup is specified
if (outgroup != "no") {
  outgroup_samples <- unlist(str_split(outgroup, ","))

  # Check if all outgroup samples are in the tree's tip labels
  if (all(outgroup_samples %in% tree$tip.label)) {
    tree <- root(tree, outgroup = outgroup_samples)
  } else {
    cat("Outgroup samples not found in tree's tip labels. No rooting performed.\n")
  }
}

# Handle missing or invalid edge lengths
tree$edge.length[is.nan(tree$edge.length)] <- 0

# Calculate tree height (x-axis)
Xmax <- max(nodeHeights(tree))

# Filter bootstrapping support: only retain > 75%
tree$node.label <- ifelse(as.numeric(tree$node.label) >= 0.75, tree$node.label, NA)

# Adjust tip label parameters
tip_label_size <- 3 # Size of tip labels
tip_label_align <- TRUE # Align tip labels horizontally
tip_label_linesize <- 0.5 # Connecting line size for tip labels

# Plot tree with or without watermark
PLOT.tree <- ggtree(tree, layout = "roundrect") +
  ggtitle("ASTRAL") +
  theme_tree2() +
  theme_bw() +
  ggplot2::xlim(0, Xmax + Xmax * offset) +
  xlab("coalescent units") +
  geom_nodelab(
    hjust = 1.25,
    vjust = -0.75,
    size = 3,
    color = "#21614a"
  ) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  coord_cartesian(clip = "off") +
  geom_tiplab(
    align = tip_label_align,
    size = tip_label_size,
    linesize = tip_label_linesize
  )

if (watermark == "YES") {
  PLOT.tree <- PLOT.tree +
    annotate(
      geom = "text",
      x = (Xmax + Xmax * offset) / 2,
      y = 1 + length(tree$node.label) / 2,
      label = "PREVIEW",
      color = "#ff8e3e",
      angle = 45,
      fontface = "bold",
      size = 25,
      alpha = 0.5
    )
}

# Export tree to PNG and PDF
PNG <- paste0(output, ".png")
PDF <- paste0(output, ".pdf")

ggsave(
  filename = PDF,
  plot = PLOT.tree,
  width = width,
  height = height,
  limitsize = FALSE
)
ggsave(
  filename = PNG,
  plot = PLOT.tree,
  width = width,
  height = height,
  limitsize = FALSE
)
