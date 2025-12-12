library(scLinaX)
library(ggplot2)
library(readr)
library(dplyr)


data("XCI_ref")
data("AIDA_QCREF")
multiome_ASE_df = read_tsv("Multiome_RNA_QC_passed_SNP_df")

scLinaX_res<-run_scLinaX(ASE_df=multiome_ASE_df,XCI_ref=XCI_ref,QCREF=AIDA_QCREF,
                             Inactive_Gene_ratio_THR=0.05,SNP_DETECTION_DP=30,SNP_DETECTION_MAF=0.1,QC_total_allele_THR=10,
                             HE_allele_cell_number_THR=50,REMOVE_ESCAPE=TRUE,PVAL_THR=0.01,RHO_THR=0.5)
                            
print(head(scLinaX_res$results))

print(head(scLinaX_res$raw_exp_result))

PBMC_summary<-summarize_scLinaX(scLinaX_res,QC_total_allele_THR=10,Annotation=NULL)

print(head(PBMC_summary))

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
ggsave("./sclinX_test_boxplot.pdf",device="png", width=6,height=4)