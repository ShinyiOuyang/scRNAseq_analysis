library(data.table)
library(Rsamtools)
library(GenomicRanges)
library(ggplot2)
library(data.table)
library(dplyr)

# Load the SAM file
bam_file = "/Users/shinyiouyang/Downloads/to_scp/HRR1795856/HRR1795856_s162p_subset.bam"
annotation_file = "/Users/shinyiouyang/Downloads/scRNAseq_analysis/HRR1795856.xchrom_annotation.csv"

params = ScanBamParam(tag = c("CB"), what=c("cigar", "pos", "qual", "seq"))
bam_data = scanBam(bam_file, param = params)


#We get back a list of length 1; this is because scanBam() can return output from multiple genomic regions, 
# and here we have only one (everything). 
#We therefore subset the output; this again gives us a list and we show the information from the first alignment
bam_data = bam_data[[1]]

# Define the target position
target_pos <- 79171491
SNP = target_pos

# Initialize counts for each cell
cell_counts <- list()

# Function to parse the CIGAR string and map the target position to the read sequence
parse_cigar <- function(cigar, read_start, target_pos) {
    # Split CIGAR string into operations (e.g., 58M252485N32M)
    cigar_ops <- unlist(strsplit(cigar, "(?<=\\D)(?=\\d)|(?<=\\d)(?=\\D)", perl = TRUE))
    cigar_lengths <- as.numeric(cigar_ops[seq(1, length(cigar_ops), by = 2)])
    cigar_types <- cigar_ops[seq(2, length(cigar_ops), by = 2)]

    # Initialize position tracking
    ref_pos <- read_start
    read_pos <- 1

    # Iterate through CIGAR operations
    for (i in seq_along(cigar_types)) {
        op_len <- cigar_lengths[i]
        op_type <- cigar_types[i]

        if (op_type == "M") {
            # Match: target position corresponds to read sequence
            if (target_pos >= ref_pos && target_pos < ref_pos + op_len) {
                read_offset <- target_pos - ref_pos + read_pos
                return(read_offset)
            }
            ref_pos <- ref_pos + op_len
            read_pos <- read_pos + op_len
        } else if (op_type == "N") {
            # Skipped region: target position is in the reference but not in the read
            if (target_pos >= ref_pos && target_pos < ref_pos + op_len) {
                return(NA)  # Target position is in the skipped region
            }
            ref_pos <- ref_pos + op_len
        } else if (op_type == "I") {
            # Insertion: position exists in the read but not the reference
            read_pos <- read_pos + op_len
        } else if (op_type == "D") {
            # Deletion: position exists in the reference but not the read
            ref_pos <- ref_pos + op_len
        } else if (op_type %in% c("S", "H")) {
            # Clipping: soft/hard clipping, does not consume reference positions
            if (op_type == "S") {
                read_pos <- read_pos + op_len
            }
        }
    }

    # If target position is not found in the alignment
    return(NA)
}

# Iterate through rows of the SAM file
#print(names(sam_data))
for (i in seq_len(length(bam_data$cigar))) {
    # Extract necessary fields
    read_start <- bam_data$pos[i]  # Start position of the read
    cigar <- bam_data$cigar[i]       # CIGAR string
    seq <- bam_data$seq[i]        # SEQ field
    qual <- bam_data$qual[i]      # QUAL field
    #tags <- bam_data$tag  # Optional tags
    cell_barcode = bam_data$tag$CB[i]

    # Skip rows with missing or invalid fields
    if (is.na(read_start) || is.na(cigar) || is.na(seq)) next

    # Map the target position to the read sequence using the CIGAR string
    read_offset <- parse_cigar(cigar, read_start, target_pos)

    # Skip reads where the target position is not mapped (deletion or skipped region)
    if (is.na(read_offset) || read_offset > nchar(seq)) next

    # # Extract cell barcode (CB tag)
    # cell_barcode <- NA
    # for (tag in tags) {
    #     if (grepl("^CB:Z:", tag)) {
    #         cell_barcode <- sub("^CB:Z:", "", tag)
    #         break
    #     }
    # }

    # Skip rows without cell barcodes
    if (is.na(cell_barcode)) next

    # Extract base and its quality score at the mapped position
    base <- substr(seq, read_offset, read_offset)  # Base at the mapped position
    base_quality <- as.integer(charToRaw(substr(qual, read_offset, read_offset))) - 33  # Phred score

    # Apply base quality filter (e.g., minimum Phred score of 30)
    if (base_quality < 30) next

    # Count C and T bases for the cell
    if (base %in% c("C", "T")) {
        if (!cell_barcode %in% names(cell_counts)) {
            cell_counts[[cell_barcode]] <- c(C = 0, T = 0)
        }
        cell_counts[[cell_barcode]][base] <- cell_counts[[cell_barcode]][base] + 1
    }
    cell_counts[[cell_barcode]]["read_pos"] = read_start

}
print(length(bam_data$cigar))


# Convert list to a data frame
cell_counts_df <- do.call(rbind, lapply(names(cell_counts), function(cell) {
    data.frame(Cell = cell, C = cell_counts[[cell]]["C"], T = cell_counts[[cell]]["T"], Read_Start = cell_counts[[cell]]["read_pos"])
}))

# Print results
print(cell_counts_df)

# Save results to a CSV file
output_file <- "C2_R18_base_counts.csv"
write.csv(cell_counts_df, output_file, row.names = FALSE)

#for (i in c("C2","C3"))

#Read Bam
#region = GRanges("chrX", IRanges(start = SNP, end = SNP))

param <- ScanBamParam(
    what = c("qname", "seq", "pos","cigar"),  # Extract sequence and position
    tag = c("CB")  # Extract cell barcode tag
)
reads = scanBam(bam_file, param = param)

read_start = reads[[1]]$pos
read_end = read_start + nchar(reads[[1]]$seq) - 1
overlapping_reads = (SNP >= read_start & SNP <= read_end & reads[[1]]$cigar=="98M")


# Extract relevant fields
reads_df = cbind(cell_barcode = reads[[1]]$tag$CB, sequence = as.character(reads[[1]]$seq), position = reads[[1]]$pos)

write.csv(reads_df, file = "GPR_174_reads_df.csv", row.names = FALSE)
overlapping_reads_df = reads_df[overlapping_reads,]
overlapping_reads_df = cbind(overlapping_reads_df, R18 = substr(overlapping_reads_df[,2],SNP - as.numeric(overlapping_reads_df[,3])+1,SNP - as.numeric(overlapping_reads_df[,3])+1))
overlapping_reads_df = na.omit(overlapping_reads_df, na.action = "omit")
print(overlapping_reads_df)

write.csv(overlapping_reads_df, file = "GPR_174_overlapping_df_qced.csv", row.names = FALSE)

GPR174cells = overlapping_reads_df[,c(1,4)]

annotation = read.csv(annotation_file, header = TRUE)

GPR174cell_annotation = merge(GPR174cells, annotation, by = "cell_barcode")

counts <- GPR174cell_annotation %>%
  group_by(R18, Xa) %>%
  summarize(count = n(), .groups = "drop")

ggplot(data = counts, aes(x = R18, y = count, fill = Xa)) +
  geom_bar(stat = "identity") +  # Creates a bar graph with actual counts
  labs(
    title = paste0("Assignment of chrX:",SNP," to Allele A or B"),
    x = "SNP",
    y = "Count"
  ) +
  scale_fill_manual(
    values = c("Allele_A" = "navy", "Allele_B" = "darkred"),
    name = "Allele"
  ) +
  theme_minimal()

ggsave(paste0(i,"_",SNP,"_annotation.pdf"), width = 6, height = 4)

write.csv(GPR174cell_annotation, file = "HRR1795890_GPR_174_cell_annotation.csv", row.names = FALSE)

# cell_annotation = read.csv(paste0("/Users/yun-huang/Documents/GPR174/scRNAseq/20250712_resequence_analysis/",i,"_resequence.CellTypes_Azimuth_simple.csv"))

# GPR174celltype_annotation = merge(cell_annotation, GPR174cell_annotation,by = "cell_barcode")

# write.csv(GPR174celltype_annotation, file = paste0(i,"_",SNP,"_categorization.csv"))

