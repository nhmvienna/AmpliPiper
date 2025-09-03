library("TreeDist")
library(tidyverse)
library(reshape2)

args <- commandArgs(trailingOnly = TRUE)

input <- args[1]
output <- args[2]

# List the file names
file_names <- list.files(
    path = input, pattern = "\\.treefile$",
    full.names = TRUE,
    recursive = TRUE
)

# Create an empty list to store the trees
tree_list <- list()

# Read each tree file and store it in the list
for (i in seq_along(file_names)) {
    tree <- ape::read.tree(file_names[i])
    tree_list[[gsub(
        ".treefile",
        "",
        basename(file_names[i])
    )]] <- tree
}

# View the list of trees
# print(tree_list)

# Function to count common leaves between trees in the list
count_common_leaves <- function(tree_list) {
    n <- length(tree_list)
    common_leaves <- matrix(0,
        nrow = n,
        ncol = n,
        dimnames = list(
            names(tree_list),
            names(tree_list)
        )
    )

    no_common_leaves_count <- 0

    for (i in 1:(n - 1)) {
        for (j in (i + 1):n) {
            common <- intersect(
                tree_list[[i]]$tip.label,
                tree_list[[j]]$tip.label
            )
            common_leaves[
                names(tree_list)[i],
                names(tree_list)[j]
            ] <- length(common)
            common_leaves[
                names(tree_list)[j],
                names(tree_list)[i]
            ] <- length(common)
            if (length(common) == 0) {
                no_common_leaves_count <- no_common_leaves_count + 1
            }
        }
    }

    return(list(
        common_leaves = common_leaves,
        no_common_leaves_count = no_common_leaves_count
    ))
}

# Call the function to count common leaves
result <- count_common_leaves(tree_list)

common_leaves_matrix <- result$common_leaves
no_common_leaves_count <- result$no_common_leaves_count

# print(common_leaves_matrix)
print(paste("Number of cases where there are no common leaves among trees:", no_common_leaves_count))

## only proceed if there is overlap among the samples in the tres
if (no_common_leaves_count == 0) {
    raw <- MutualClusteringInfo(tree_list,
        normalize = T
    )

    raw_melt <- melt(raw)
    raw_melt$value <- ifelse(raw_melt$Var1 == raw_melt$Var2, NA, raw_melt$value)
    PLOT <- ggplot(raw_melt, aes(Var1, Var2, fill = value)) +
        geom_tile(color = "white") +
        scale_fill_gradientn(colours = c("blue", "red"), values = c(0, 0.5, 1), na.value = "grey", limits = c(0, 1)) +
        labs(x = "Trees", y = "Trees", title = "Normalized Robinson Foulds Distance") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        coord_fixed() +
        theme(legend.position = "right")

    ggsave(
        file = output,
        PLOT,
        width = 10,
        height = 8
    )
}
