library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)

samples <- c('HRR1795890', 'HRR1795892')
outdir = "/Users/shinyiouyang/Downloads/scRNAseq_analysis/output/"

for (samplename in samples) {
    curr_outdir = paste0(outdir, samplename, "/")
    xchrom_annotation = read.csv(paste0(curr_outdir, samplename, ".xchrom_annotation.csv"), header = TRUE)
    cell_types = read.csv(paste0(curr_outdir, samplename, "_CellTypes_Azimuth_simple.csv"), header = TRUE)

    xchrom_annotation_bycelltype = merge(xchrom_annotation, cell_types, by = "cell_barcode")

    agg_data <- xchrom_annotation_bycelltype %>%
    group_by(Annotation, Xa) %>%
    summarise(count = n()) %>%
    ungroup()

    write.csv(agg_data,paste0(curr_outdir, samplename,".xchrom_annotation_bycelltype.csv"))

    #Create Graph
    pbmc1_starting_cells_grouped = cell_types %>%
    group_by(Annotation) %>%
    summarise(count = n_distinct(cell_barcode))

    # Calculate total Xa counts per Annotation
    total_xa_counts <- agg_data %>%
    group_by(Annotation) %>%
    summarise(total_xa_count = sum(count), .groups = 'drop')

    # Calculate remaining counts
    remaining_counts <- pbmc1_starting_cells_grouped %>%
    left_join(total_xa_counts, by = "Annotation") %>%
    mutate(total_xa_count = ifelse(is.na(total_xa_count), 0, total_xa_count),
            remaining_count = count - total_xa_count) %>%
    select(Annotation, remaining_count)

    # Combine all counts into one data frame
    combined_data <- remaining_counts %>%
    mutate(Xa = "Uncategorized") %>%
    rename(xa_count = remaining_count) %>%
    select(Annotation, Xa, xa_count) %>%
    bind_rows(agg_data %>%
                rename(xa_count = count))

    # Ensure all Annotation categories from pbmc1_starting_cells_grouped are included
    combined_data <- pbmc1_starting_cells_grouped %>%
    select(Annotation) %>%
    left_join(combined_data, by = "Annotation") %>%
    mutate(xa_count = ifelse(is.na(xa_count), 0, xa_count),
            Xa = ifelse(is.na(Xa), "Remaining", Xa))

    # Remove rows with zero counts (optional)
    combined_data <- combined_data %>%
    filter(xa_count > 0)

    # Set factor levels for proper stacking order
    combined_data$Xa <- factor(combined_data$Xa, levels = c("Allele_A", "Allele_B", "Uncategorized"))

    # Print combined_data for verification
    print(combined_data)

    # Create the bar graph
    plot1 <- ggplot(combined_data, aes(x = Annotation, y = xa_count, fill = Xa)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = c("Allele_A" = "maroon", "Allele_B" = "navy", "Uncategorized" = "grey")) +
    labs(y = "Count", x = "Annotation", fill = "Xa") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

    ggsave(paste0(curr_outdir, samplename,".subsets.pdf"), plot = plot1, width = 8, height = 6)
    write.csv(combined_data[,c(1:3)],paste0(curr_outdir, samplename,".subset_stats.csv"))

    counts_df <- combined_data %>%
        filter(Xa %in% c("Allele_A", "Allele_B")) %>%
        group_by(Annotation, Xa) %>%
        summarise(xa_count = sum(xa_count)) %>%
        spread(Xa, xa_count, fill = 0) %>%
        mutate(total = Allele_A + Allele_B) %>%
        select(Annotation, Allele_A, Allele_B, total)

    plot_with_labels <- ggplot(counts_df, aes(x = Allele_A, y = Allele_B, color = Annotation, label = Annotation)) +
        geom_point() +
        geom_text(aes(label = Annotation), size = 3, vjust = -0.5, hjust = 0.5) +
        theme_bw() +
        labs(x = "Cell Count for Allele_A", y = "Cell Count for Allele_B", title = paste0(samplename, " Scatter Plot with All Labels"))

    plot(plot_with_labels)
    ggsave(paste0(curr_outdir, samplename,".scatter_with_labels.pdf"), plot = plot_with_labels, width = 8, height = 6)

    plot_with_labels_log <- ggplot(counts_df, aes(x = Allele_A, y = Allele_B, color = Annotation, label = Annotation)) +
        geom_point() +
        geom_text_repel(aes(label = Annotation), size = 3) +  # Use geom_text_repel to avoid overlapping
        theme_bw() +
        labs(x = "Cell Count for Allele_A (log scale)", y = "Cell Count for Allele_B (log scale)", title = paste0(samplename, " Scatter Plot with All Labels Log Scale")) +
        scale_x_log10() +
        scale_y_log10()

    print(plot_with_labels_log)

    ggsave(paste0(curr_outdir, samplename,".scatter_logscale.pdf"), plot = plot_with_labels_log, width = 12, height = 9)

}