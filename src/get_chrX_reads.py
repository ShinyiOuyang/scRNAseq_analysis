# Python script
# Extract the X chromosome
# Replace the suffix of the barcode (-1) with the data modality (RNA or ATAC)
# Renaming the barcode is actually unnecessary in this case, but it will be necessary if you want to merge multiple runs of the same sample.
# Therefore, we are changing the barcode name here as an example.

import pysam
import pandas as pd
import subprocess
import time
import sys

bam_files = []
DIR = sys.argv[1]  # Get directory from command line argument
curr_file = "possorted_genome_bam.bam"
BAM_LIST = [
    "/possorted_genome_bam"
]
# BAM_LIST = [
#     "pbmc_granulocyte_sorted_10k_gex_possorted_bam",
#     "pbmc_granulocyte_sorted_10k_atac_possorted_bam",
# ]
start_time = time.time()
print(start_time)

for BAM in BAM_LIST:
    input_bam = DIR + BAM + ".bam"
    output_bam = DIR + BAM + ".barcode.renamed.X.bam"
    bam_files.append(output_bam)

    in_bam = pysam.AlignmentFile(input_bam, "rb")
    out_bam = pysam.AlignmentFile(output_bam, "wb", template=in_bam)

    for aln in in_bam.fetch("chrX"):
        out_bam.write(aln)
    #         if aln.has_tag("CB"):
    #             old_bc = aln.get_tag("CB")
    #             new_bc = old_bc.replace("-1", "-RNA")
    #             aln.set_tag("CB", new_bc)
    #

    # elif BAM == "pbmc_granulocyte_sorted_10k_atac_possorted_bam":
    #     for aln in in_bam.fetch("chrX"):
    #         if aln.has_tag("CB"):
    #             old_bc = aln.get_tag("CB")
    #             new_bc = old_bc.replace("-1", "-ATAC")
    #             aln.set_tag("CB", new_bc)
    #             out_bam.write(aln)

    in_bam.close()
    out_bam.close()

location = DIR + "/possorted_genome_X.sorted.bam"
subprocess.run(
    [
        "samtools",
        "sort",
        "-o",
        location,
        DIR + "/possorted_genome_bam.barcode.renamed.X.bam",
    ],
    check=True,
)
subprocess.run(["samtools", "index", location], check=True)

print(f"Total time taken: {time.time() - start_time} seconds")
