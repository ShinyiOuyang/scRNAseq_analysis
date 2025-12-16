library(scLinaX)
library(ggplot2)
library(readr)
library(dplyr)


samplename="sample_HRR1795888"
data("XCI_ref")
data("AIDA_QCREF")
#data("multiome_ASE_df")
multiome_ASE_df = read_tsv("HRR1795888_RNA_QC_passed_SNP_df.tsv")
multiome_ASE_df$Sample_ID = "sample_HRR1795896"


# Inactive_Gene_ratio_THR=0.00,SNP_DETECTION_DP=0,SNP_DETECTION_MAF=0,QC_total_allele_THR=10,
#                              HE_allele_cell_number_THR=0,REMOVE_ESCAPE=TRUE,PVAL_THR=0.01,RHO_THR=0.5)
scLinaX_res<-run_scLinaX(ASE_df=multiome_ASE_df,XCI_ref=XCI_ref,QCREF=AIDA_QCREF,
                             Inactive_Gene_ratio_THR=0.05,SNP_DETECTION_DP=30,SNP_DETECTION_MAF=0.1,QC_total_allele_THR=10,
                             HE_allele_cell_number_THR=50,REMOVE_ESCAPE=TRUE,PVAL_THR=0.01,RHO_THR=0.5)
                            

PBMC_summary<-summarize_scLinaX(scLinaX_res,QC_total_allele_THR=0,Annotation=NULL)


if (require("ggplot2")) {
    PBMC_summary$Gene_class<-factor(PBMC_summary$Gene_class,levels=c("PAR1","nonPAR_escape","nonPAR_variable","nonPAR_inactive","nonPAR_unknown","PAR2"))

    p<-ggplot(PBMC_summary,aes(x=Gene_class,y=minor_allele_ratio,fill=Gene_class))+
    geom_boxplot()+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),aspect.ratio=2/3,legend.position="none")+
    scale_fill_manual(
        values=c("PAR1"="#ae2d68","nonPAR_escape"="#ff1700","nonPAR_variable"="#ffa600",
        "nonPAR_inactive"="#9db300","nonPAR_unknown"="#666666","PAR2"="#729efd"))+
    xlab("")+ylab("ratio of the expression from Xi")
    p
}
ggsave("./sclinaX_test_boxplot.png",device="png", width=6,height=4)

write.csv(scLinaX_res$Max_Num_Table_result, paste0(samplename,".sclinax_cell_counts.csv"), row.names = FALSE)
write.csv(scLinaX_res$df_snp_summary, paste0(samplename,".sclinax_snp_summary.csv"), row.names = FALSE)
write.csv(scLinaX_res$Fail_list, paste0(samplename,".fail_list.csv"), row.names = FALSE)
write.csv(PBMC_summary, paste0(samplename,".summarize.csv"), row.names = FALSE)
write.csv(scLinaX_res$raw_exp_result, paste0(samplename,".sclinax_raw_exp.csv"), row.names = FALSE)

#Pick out uniquely labeled cells
xchrom_annotation = scLinaX_res$raw_exp_result[,c(7,18)]
#xchrom_annotation = scLinaX_res$raw_exp_result[,c(6,16)]
head(xchrom_annotation)


unique_assignments <- xchrom_annotation %>%
  group_by(cell_barcode) %>%
  summarise(unique_count = n_distinct(Xa))

unique_barcodes <- unique_assignments %>%
  filter(unique_count == 1) %>%
  select(cell_barcode)

unique_xchrom_annotation <- xchrom_annotation %>%
  inner_join(unique_barcodes, by = "cell_barcode")

unique_xchrom_annotation_deduplicated = unique_xchrom_annotation %>% distinct(cell_barcode, .keep_all = TRUE)

write.csv(xchrom_annotation, paste0(samplename,".xchrom_annotation.csv"), row.names = FALSE)

barplot = ggplot(unique_xchrom_annotation_deduplicated) + 
  geom_bar(aes(x=Xa, fill=Xa))

ggsave("./allele_boxplot.png",device="png", width=6,height=4)

