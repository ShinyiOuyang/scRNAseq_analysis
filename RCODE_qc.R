#RCODE_QC.r
library(readr)
library(dplyr)
library(stringr)
library(viridis)
library(ggpointdensity)
library(patchwork)

# argv<-commandArgs(T)

# RNA_SNPCELL=argv[1]
# RNA_ANNOVAR=argv[2]
# ATAC_SNPCELL=argv[3]
# ATAC_ANNOVAR=argv[4]

RNA_SNPCELL="RNA.chrX.snp.call.summary.tsv.gz"
RNA_ANNOVAR="RNA.chrX.snp.call.hg38_multianno.txt" 
# ATAC_SNPCELL="ATAC.chrX.snp.call.summary.tsv.gz"
# ATAC_ANNOVAR="ATAC.chrX.snp.call.hg38_multianno.txt"

#load raw data
rna<-read_tsv(RNA_SNPCELL)

#load annotation data
load_anno<-function(ANNOVAR){
anno<-read_tsv(ANNOVAR)%>%
mutate(SNP_ID=str_c(Chr,Start,Ref,Alt,sep=":"))%>%
dplyr::select(SNP_ID,Func.refGene,Gene.refGene,ALL.sites.2015_08)
colnames(anno)<-c("SNP_ID","Region","Gene","ALL_Freq")
anno$ALL_Freq[is.na(anno$ALL_Freq)]<-0
return(anno)
}

rna_anno<-load_anno(RNA_ANNOVAR)

 
#merge raw data + annotation data
merge_anno<-function(data,anno){
df<-left_join(data,anno,by="SNP_ID")%>%
filter(!str_detect(Gene,pattern=";"))%>%
filter(!is.na(Gene))%>%filter(Region %in% c("intronic","UTR5","UTR3","exonic","ncRNA_exonic","ncRNA_intronic","splicing"))
return(df)
}
rna_df<-merge_anno(rna,rna_anno)

rna_df_QC_pass<-filter(rna_df,ALL_Freq>0.01)
write_tsv(rna_df_QC_pass,"HRR1795889_RNA_QC_passed_SNP_df.tsv.gz") #This QC-passed data should be used for the subsequent analyses