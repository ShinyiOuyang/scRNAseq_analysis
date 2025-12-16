#RCODE_process_cellsnp.r
library(dplyr)
library(readr)
library(stringr)
library(tidyr)

argv<-commandArgs(T)

SAMPLE_NAME=argv[1]
ADmtx=argv[2]
DPmtx=argv[3]
OTHmtx=argv[4]
VCF=argv[5]
OUT=argv[6]

sample_name<-read_tsv(SAMPLE_NAME,col_names=FALSE)
AD<-Matrix::readMM(ADmtx)
DP<-Matrix::readMM(DPmtx)
OTH<-Matrix::readMM(OTHmtx)
RD<-DP-AD
vcf<-read_tsv(VCF,comment="#",col_names=FALSE)

SNP_INFO<-mutate(vcf,SNP_ID=str_c(X1,":",X2,":",X4,":",X5))%>%
dplyr::select(SNP_ID,X1,X2,X4,X5)
colnames(SNP_INFO)<-c("SNP_ID","CHR","POS","REF","ALT")

alt_res<-AD%>%as.matrix()%>%as_tibble()
colnames(alt_res)<-sample_name$X1
alt_res$SNP<-SNP_INFO$SNP_ID

ref_res<-RD%>%as.matrix()%>%as_tibble()
colnames(ref_res)<-sample_name$X1
ref_res$SNP<-SNP_INFO$SNP_ID

oth_res<-OTH%>%as.matrix()%>%as_tibble()
colnames(oth_res)<-sample_name$X1
oth_res$SNP<-SNP_INFO$SNP_ID

alt_res<-alt_res%>%
pivot_longer(names_to="cell_barcode",values_to="ALTcount",-SNP)%>%
filter(ALTcount>0)

ref_res<-ref_res%>%
pivot_longer(names_to="cell_barcode",values_to="REFcount",-SNP)%>%
filter(REFcount>0)

oth_res<-oth_res%>%
pivot_longer(names_to="cell_barcode",values_to="OTHcount",-SNP)%>%
filter(OTHcount>0)

result<-full_join(ref_res,alt_res,by=c("SNP","cell_barcode"))%>%
full_join(oth_res,by=c("SNP","cell_barcode"))
result<-left_join(SNP_INFO,result,alt_res,by=c("SNP_ID"="SNP"))
result[is.na(result)]<-0
result<-dplyr::distinct(result)

write_tsv(result,OUT)