# Load required libraries
library(seqinr)
library(ape)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)

input <- args[1]
outgroups_arg <- if (length(args) >= 2) args[2] else ""
TH <- as.numeric(ifelse(length(args) >= 3, args[3], 2))

# Parse outgroups from comma-separated list
outgroups <- if (outgroups_arg != "") {
    trimws(unlist(strsplit(outgroups_arg, ",")))
} else {
    character(0)
}

setwd(input)

# Function to identify outlier taxa based on average distance
identify_outliers <- function(file_path, threshold_sd = 0.5, exclude_taxa = character(0)) {
    tryCatch(
        {
            # Read the alignment file
            alignment <- read.alignment(file_path, format = "fasta")

            # Calculate genetic distances
            distances <- dist.alignment(alignment)

            # Convert distances to a matrix
            distance_matrix <- as.matrix(distances)

            # Calculate average distance for each taxon (row means, excluding diagonal)
            n_taxa <- nrow(distance_matrix)
            avg_distances <- numeric(n_taxa)
            taxa_names <- rownames(distance_matrix)

            for (i in 1:n_taxa) {
                # Calculate average distance excluding self (diagonal = 0)
                avg_distances[i] <- mean(distance_matrix[i, -i])
            }

            # Filter out excluded taxa (outgroups) from outlier analysis
            if (length(exclude_taxa) > 0) {
                cat(sprintf("    Outgroups to exclude: %s\n", paste(exclude_taxa, collapse = ", ")))
                cat(sprintf("    Available taxa: %s\n", paste(taxa_names, collapse = ", ")))

                # Remove suffixes from taxa names for matching (everything after the last "-")
                taxa_base_names <- sub("-[^-]*$", "", taxa_names)
                cat(sprintf("    Taxa base names: %s\n", paste(taxa_base_names, collapse = ", ")))

                # Find indices of taxa that are NOT in the exclude list (based on base names)
                included_indices <- which(!taxa_base_names %in% exclude_taxa)
                analysis_avg_distances <- avg_distances[included_indices]

                cat(sprintf("    Included taxa count: %d (out of %d total)\n", length(included_indices), n_taxa))
                if (length(included_indices) > 0) {
                    cat(sprintf("    Included taxa: %s\n", paste(taxa_names[included_indices], collapse = ", ")))
                }
            } else {
                included_indices <- 1:n_taxa
                analysis_avg_distances <- avg_distances
                cat(sprintf("    No outgroups to exclude, analyzing all %d taxa\n", n_taxa))
            }

            # Only proceed if we have taxa to analyze
            if (length(analysis_avg_distances) == 0) {
                return(NULL)
            }

            # Calculate mean and standard deviation of average distances (excluding outgroups)
            mean_avg_dist <- mean(analysis_avg_distances)
            sd_avg_dist <- sd(analysis_avg_distances)

            # Identify outliers among the included taxa only
            outlier_threshold <- mean_avg_dist + threshold_sd * sd_avg_dist
            outlier_indices_in_analysis <- which(analysis_avg_distances > outlier_threshold)

            # Map back to original indices
            outlier_indices <- included_indices[outlier_indices_in_analysis]

            # Extract gene/locus name from file path
            file_name <- basename(file_path)
            locus <- str_extract(file_name, ".*(?=_aln.fasta)")

            # Create results data frame
            if (length(outlier_indices) > 0) {
                outlier_data <- data.frame(
                    Locus = locus,
                    Taxon = taxa_names[outlier_indices],
                    Average_Distance = avg_distances[outlier_indices],
                    Mean_Distance = mean_avg_dist,
                    SD_Distance = sd_avg_dist,
                    Threshold = outlier_threshold,
                    Deviation_from_Mean = avg_distances[outlier_indices] - mean_avg_dist,
                    SD_Units = (avg_distances[outlier_indices] - mean_avg_dist) / sd_avg_dist,
                    stringsAsFactors = FALSE
                )
                return(outlier_data)
            } else {
                return(NULL)
            }
        },
        error = function(e) {
            message(sprintf("Error processing file %s: %s", file_path, e$message))
            return(NULL)
        }
    )
}

# Function to read alignment file and plot distance matrix
read_and_plot <- function(file_path, max_dist) {
    tryCatch(
        {
            # Read the alignment file
            alignment <- read.alignment(file_path, format = "fasta")

            # Calculate genetic distances
            distances <- dist.alignment(alignment)**2

            # Convert distances to a distance matrix
            distance_matrix <- as.matrix(distances)

            # Convert to upper triangle matrix
            distance_matrix[lower.tri(distance_matrix)] <- NA

            # Convert to data frame
            distance_df <- melt(distance_matrix, na.rm = TRUE)

            # Extract ID from file name
            file_name <- basename(file_path)
            id <- str_extract(file_name, ".*(?=_aln.fasta)")

            # Plot the distance matrix as a heatmap using ggplot2
            p <- ggplot(data = distance_df, aes(Var2, Var1, fill = value)) +
                geom_tile() +
                scale_fill_gradient(low = "#7ac1dc", high = "#ff8e3e", limits = c(0, max_dist)) +
                theme_minimal() +
                theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
                labs(title = paste("Pairwise Genetic Distances -", id), x = "Sequence", y = "Sequence")

            return(p)
        },
        error = function(e) {
            message(sprintf("Error processing file %s: %s", file_path, e$message))
            return(NULL)
        }
    )
}

# Get a list of all files with the extension "_aln.fasta" in the current directory and its subdirectories
file_paths <- list.files(path = ".", pattern = "_aln.fasta$", recursive = TRUE, full.names = TRUE)

if (length(file_paths) == 0) {
    stop("No '_aln.fasta' files found in the directory.")
}

# OUTLIER ANALYSIS: Identify taxa with average distances > 0.5 SD from mean
cat("\n=== OUTLIER ANALYSIS ===\n")
cat("Identifying taxa whose average distance is larger than 0.5 SD from the mean\n")

if (length(outgroups) > 0) {
    cat(sprintf("Excluding outgroups from analysis: %s\n", paste(outgroups, collapse = ", ")))
} else {
    cat("No outgroups specified for exclusion\n")
}
cat("\n")

all_outliers <- data.frame()

for (file_path in file_paths) {
    cat(sprintf("Processing file: %s\n", basename(file_path)))
    outliers <- identify_outliers(file_path, threshold_sd = TH, exclude_taxa = outgroups)
    if (!is.null(outliers)) {
        all_outliers <- rbind(all_outliers, outliers)
    }
}

if (nrow(all_outliers) > 0) {
    cat("OUTLIER TAXA DETECTED:\n")
    cat("======================\n")

    # Print detailed results
    if (nrow(all_outliers) > 0) {
        for (i in seq_len(nrow(all_outliers))) {
            cat(sprintf("Locus: %s\n", all_outliers$Locus[i]))
            cat(sprintf("  Taxon: %s\n", all_outliers$Taxon[i]))
            cat(sprintf("  Average Distance: %.4f\n", all_outliers$Average_Distance[i]))
            cat(sprintf("  Population Mean: %.4f\n", all_outliers$Mean_Distance[i]))
            cat(sprintf("  Population SD: %.4f\n", all_outliers$SD_Distance[i]))
            cat(sprintf("  Threshold (Mean + 0.5*SD): %.4f\n", all_outliers$Threshold[i]))
            cat(sprintf("  Deviation from Mean: %.4f\n", all_outliers$Deviation_from_Mean[i]))
            cat(sprintf("  Standard Deviations above Mean: %.2f\n", all_outliers$SD_Units[i]))
            cat("\n")
        }
    }

    # Summary by locus
    cat("SUMMARY BY LOCUS:\n")
    cat("=================\n")
    locus_summary <- aggregate(Taxon ~ Locus, data = all_outliers, FUN = function(x) paste(x, collapse = ", "))
    if (nrow(locus_summary) > 0) {
        for (i in seq_len(nrow(locus_summary))) {
            cat(sprintf("%s: %s\n", locus_summary$Locus[i], locus_summary$Taxon[i]))
        }
    }
    cat("\n")

    # Save outlier results to CSV
    write.csv(all_outliers, "outlier_taxa.csv", row.names = FALSE, quote = FALSE)
    cat("Detailed outlier analysis saved to 'outlier_taxa.csv'\n\n")
} else {
    cat("No outlier taxa detected across all loci.\n\n")

    # Create empty CSV with headers
    empty_outliers <- data.frame(
        Locus = character(0),
        Taxon = character(0),
        Average_Distance = numeric(0),
        Mean_Distance = numeric(0),
        SD_Distance = numeric(0),
        Threshold = numeric(0),
        Deviation_from_Mean = numeric(0),
        SD_Units = numeric(0),
        stringsAsFactors = FALSE
    )
    write.csv(empty_outliers, "outlier_taxa.csv", row.names = FALSE, quote = FALSE)
    cat("Empty outlier analysis file with headers saved to 'outlier_taxa.csv'\n\n")
}

# Calculate the maximum distance across all files for plotting
max_dist <- 0
for (file_path in file_paths) {
    tryCatch(
        {
            # Read the alignment file
            alignment <- read.alignment(file_path, format = "fasta")

            # Calculate genetic distances
            distances <- dist.alignment(alignment)

            # Find the maximum distance
            max_dist <- max(max_dist, max(distances))
        },
        error = function(e) {
            message(sprintf("Error reading file %s: %s", file_path, e$message))
        }
    )
}

# Calculate the optimal dimensions for the grid
num_files <- length(file_paths)
num_cols <- ceiling(sqrt(num_files))
num_rows <- ceiling(num_files / num_cols)

# Create a list to store plots
plots <- list()

# Read each alignment file and plot the distance matrix
for (file_path in file_paths) {
    p <- read_and_plot(file_path, max_dist)
    if (!is.null(p)) {
        plots[[file_path]] <- p
    }
}

if (length(plots) == 0) {
    stop("No valid plots generated from the files.")
}

# Save the plot as a PNG file
png_output <- "distance_matrices.png"
ggsave(
    filename = png_output,
    plot = grid.arrange(grobs = plots, ncol = num_cols),
    width = num_cols * 12,
    height = num_rows * 12
)

pdf_output <- "distance_matrices.pdf"
ggsave(
    filename = pdf_output,
    plot = grid.arrange(grobs = plots, ncol = num_cols),
    width = num_cols * 12,
    height = num_rows * 12
)

cat("=== ANALYSIS COMPLETE ===\n")
cat(sprintf("Distance matrix plots saved to: %s and %s\n", png_output, pdf_output))
if (nrow(all_outliers) > 0) {
    cat(sprintf("Outlier analysis results saved to: outlier_taxa.csv\n"))
    cat(sprintf(
        "Total outlier taxa identified: %d across %d loci\n",
        nrow(all_outliers), length(unique(all_outliers$Locus))
    ))
} else {
    cat("No outlier taxa detected.\n")
}
cat("Analysis complete!\n")
