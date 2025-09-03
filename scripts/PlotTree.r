# Load necessary R libraries
library(ggtree)
library(tidyverse)
library(phytools) # For midpoint rooting and tree height calculation

args <- commandArgs(trailingOnly = TRUE)

input <- args[1]
output <- args[2]
title <- args[3]
offset <- as.numeric(args[4])
width <- as.numeric(args[5])
height <- as.numeric(args[6])
outgroup <- args[7]
watermark <- args[8]

# Load tree file and root midpoint
tree <- read.tree(input)

if (outgroup == "no") {
  tree <- midpoint.root(tree)
} else {
  outgroup_samples <- unlist(str_split(outgroup, ","))

  # Check if all outgroup samples are in the tree's leaves (tip labels)
  if (all(outgroup_samples %in% tree$tip.label)) {
    tree <- root(tree, outgroup = outgroup_samples)
  } else {
    cat("Outgroup samples not found in tree's tip labels. Midpoint rooting will be applied.\n")
    tree <- midpoint.root(tree)
  }
}

# Calculate tree height (x-axis range)
Xmax <- max(nodeHeights(tree))

# Only retain bootstrapping support > 75%
tree$node.label[as.numeric(tree$node.label) < 75] <- NA

# Set parameters for tip labels
tip_label_size <- 3 # Adjust size as needed
tip_label_align <- TRUE # Align tip labels horizontally
tip_label_linesize <- 0.5 # Line size between tips and labels

if (watermark == "YES") {
  # Plot tree with watermark
  PLOT.tree <- ggtree(tree, layout = "roundrect") +
    ggtitle(title) +
    theme_tree2() +
    theme_bw() +
    geom_rootpoint(size = 2, color = "#21614a") +
    geom_rootedge(size = 0.5, color = "#21614a") +
    ggplot2::xlim(0, Xmax + Xmax * offset) +
    xlab("av. subst./site") +
    geom_nodelab(
      hjust = 1.25,
      vjust = -0.75,
      size = 3,
      color = "#21614a"
    ) +
    annotate(
      geom = "text",
      x = (Xmax + Xmax * offset) / 2,
      y = 1 + length(tree$node.label) / 2,
      label = "PREVIEW",
      color = "#ff8e3e",
      angle = 45,
      fontface = "bold",
      size = 25,
      alpha = 0.5,
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
} else {
  # Plot tree without watermark
  PLOT.tree <- ggtree(tree, layout = "roundrect") +
    ggtitle(title) +
    theme_tree2() +
    geom_rootpoint(size = 2, color = "#21614a") +
    geom_rootedge(linewidth = 0.5, color = "#21614a") +
    theme_bw() +
    ggplot2::xlim(0, Xmax + Xmax * offset) +
    xlab("av. subst./site") +
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
