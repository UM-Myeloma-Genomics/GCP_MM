##################################################################################
##################################################################################
##############################################################################1###
#                                                                                #
#        .     .     .                                               .     .     .
####     .     .     .   PredictionModel - Genomic Classification   .     .     .     ####
#        .     .     .                                               .     .     .
#                                                                                #
##################################################################################
##################################################################################
##################################################################################

library(hdp)
library(RColorBrewer)
library(pheatmap)
library(survival)
library(ggplot2)
library(survminer)

WF="~/Google Drive/My Drive/Bachisio/Documents/Project/MM/Prediction_Model/GCP_MM/genomic/" # /where/is/the/folder/genomic

INPUT=paste(WF,"/matrixDEF_test.txt",sep="")
#######
#### . . .    - Preparation Data -     . . . ####
#######

all<- read.delim(INPUT,
                 stringsAsFactors= F, header=T)

# . . . Change translocation values 
table(all$t_CCND1)
table(all$t_MAF)
table(all$t_MMSET)
for(i in 156:158){
  all[,i]=ifelse(all[,i]==1,2,all[,i])
}

#### . . . - Only Genomic Features 
all=all[,c(1,27:158)] # Subset with only genomic features
rownames(all)=all$sample
all=all[,-1]

all=as.data.frame(all)
dim(all) # 30 132

#### . . . - Remove features
# info=read.delim("~/Google Drive/My Drive/Bachisio/Documents/Project/MM/Prediction_Model/analysis/ANALYSIS_2023/features_NAs_info.txt",
#                 stringsAsFactors = F)
# 
# info$NA.percent=info$value.NA/length(rownames(all))
# info$value.2=ifelse(is.na(info$value.2),0,info$value.2)
# info$Frequence.alt=(info$value.1+info$value.2)/length(rownames(all))
# 
# median(info$NA.percent)
# info[info$NA.percent>median(info$NA.percent),]
# median(info$Frequence.alt)
# summary=summary(info$Frequence.alt)
# 
# remove=info[info$Frequence.alt<summary[2],] # Using the 1st quartile 0.012287
# dim(remove) # 35 8

remove$feature
# [1] "SNV_IDH1"      "SNV_IDH2"      "SNV_PIK3CA"    "SNV_ABCF1"     "SNV_BHLHE41"   "SNV_DDX3X"     "SNV_DTX1"      "SNV_HIST1H1B"  "SNV_HIST1H1D"  "SNV_HIST1H2BK" "SNV_INO80E"   
# [12] "SNV_IRF1"      "SNV_LCE1D"     "SNV_MAN2C1"    "SNV_MRGPRX4"   "SNV_NDUFC2"    "SNV_PANK3"     "SNV_RFTN1"     "SNV_RPRD1B"    "SNV_RPS3A"     "SNV_SLC35G5"   "SNV_SYT17"    
# [23] "SNV_TBC1D29"   "SNV_TNFRSF11B" "SNV_USP7"      "SNV_ZFP36L1"   "SNV_ATRX"      "SNV_KDM5C"     "SNV_KDM6A"     "SNV_RPL10"     "SNV_MAF"       "SNV_MAFB"      "SNV_MAML2"    
# [34] "SNV_PIM1"      "SNV_TCL1A"


#### . . . - Remove features and samples with too much NAs ####
all=all[,!colnames(all) %in% c("SNV_IDH1","SNV_IDH2","SNV_PIK3CA","SNV_ABCF1","SNV_BHLHE41","SNV_DDX3X","SNV_DTX1","SNV_HIST1H1B","SNV_HIST1H1D","SNV_HIST1H2BK","SNV_INO80E",   
                               "SNV_IRF1","SNV_LCE1D","SNV_MAN2C1","SNV_MRGPRX4","SNV_NDUFC2","SNV_PANK3","SNV_RFTN1","SNV_RPRD1B","SNV_RPS3A","SNV_SLC35G5","SNV_SYT17",    
                               "SNV_TBC1D29","SNV_TNFRSF11B","SNV_USP7","SNV_ZFP36L1","SNV_ATRX","SNV_KDM5C","SNV_KDM6A","SNV_RPL10","SNV_MAF","SNV_MAFB","SNV_MAML2",  
                               "SNV_PIM1","SNV_TCL1A")]
dim(all) # 30 97

all=all[complete.cases(all),]
dim(all) # 30


all_gen=all
all_gen$sample=rownames(all_gen)
all_gen=all_gen[,c(ncol(all_gen),1:(ncol(all_gen)-1))]
all_gen$HRD<-0
all_gen$HRD[all_gen$CNV_chr3.gain+ all_gen$CNV_chr5.gain+all_gen$CNV_chr7.gain+
              all_gen$CNV_chr9.gain + all_gen$CNV_chr11.gain + all_gen$CNV_chr15.gain +
              all_gen$CNV_chr19.gain + all_gen$CNV_chr21.gain>=2]<- 1

all_gen$add_gains<-0
all_gen$add_gains[all_gen$CNV_Amp_14q11.2+ all_gen$CNV_Amp_17q22+all_gen$CNV_Amp_2q24.3+
                    all_gen$CNV_Amp_2q32.1 + all_gen$CNV_Amp_4q13.2 + all_gen$CNV_Amp_6p24.3 +
                    all_gen$CNV_Amp_8q24.21 + all_gen$CNV_chr18.gain>=3]<- 1


#all_gen$deletions<-0
# all_gen$deletions[all_gen$CNV_Del_10p15.3+ all_gen$CNV_Del_10q24.32+all_gen$CNV_Del_10q26.3+
#                     all_gen$CNV_Del_11q22.1 + all_gen$CNV_Del_12p13.2 + all_gen$CNV_Del_12q24.31 +
#                     all_gen$CNV_Del_14q24.3 + all_gen$CNV_Del_16p11.2 + all_gen$CNV_Del_16q22.2 + 
#                     all_gen$CNV_Del_16q24.3 + all_gen$CNV_Del_20p13 + all_gen$CNV_Del_20q13.12 + 
#                     all_gen$CNV_Del_21p11.2 + all_gen$CNV_Del_22q13.32 + all_gen$CNV_Del_2q31.1 + 
#                     all_gen$CNV_Del_2q37.3 + all_gen$CNV_Del_6q26 + all_gen$CNV_Del_7p11.2 + 
#                     all_gen$CNV_Del_7p22.2 + all_gen$CNV_Del_8p21.3 + all_gen$CNV_Del_8p23.1 + 
#                     all_gen$CNV_Del_8q24.21 + all_gen$CNV_Del_9p21.3  + all_gen$CNV.SNV_ARID1A + 
#                     all_gen$CNV.SNV_ARID2 + all_gen$CNV.SNV_ATM + all_gen$CNV.SNV_CDKN2C + all_gen$CNV.SNV_CYLD + 
#                     all_gen$CNV.SNV_FAM46C + all_gen$CNV.SNV_MAX + all_gen$CNV.SNV_NCOR1 + all_gen$CNV.SNV_NF1 +
#                     all_gen$CNV.SNV_NFKBIA + all_gen$CNV.SNV_PRDM1 + all_gen$CNV.SNV_RASA2 + all_gen$CNV.SNV_TRAF2 + 
#                     all_gen$CNV.SNV_TRAF3 + all_gen$CNV.SNV_RB1 + all_gen$CNV.SNV_RPL5 + all_gen$CNV.SNV_TP53 + all_gen$CNV.SNV_SETD2 +
#                     all_gen$CNV.SNV_SP140
#                   >=6]<- 1
all_gen$deletions1=0
for(i in 1:length(all_gen$sample)){
  all_gen[i,101]=sum(all_gen[i,c(19:69)],na.rm = T)
}
# all_gen$gains1=0
# for(i in 1:length(all_gen$sample)){
#   all_gen[i,120]=sum(all_gen[i,c(13:20)])
# }
# summary(all_gen$gains1)
# summary(all_gen$deletions1)

all_gen$deletions=ifelse(all_gen$deletions1>=5,1,0) # The 10% of the number of the features analyzed (Losses and TSGs).
# They are 50, so the threshold is 5

all_gen=all_gen[,-101]

all_gen$x1q=all_gen$CNV_Gain_Amp1q        
#all_gen$x1q[all_gen$CNV_Amp_1q21.1+all_gen$CNV_Amp_1q21.3+all_gen$CNV_Amp_1q44>=1]=1  


######
# . . From dis3_domain.R
annot=read.delim(paste(WF,"/input/snv/CoMMpass_all_coding_mutations_annot_test.txt",sep=""),
                 stringsAsFactors = F)
annot$sampleID=gsub("_1_BM","",annot$sampleID)
info_samples=read.delim(paste(WF,"/input/commpass_info_test.txt",sep=""),
                        stringsAsFactors = F)
info_samples=info_samples[info_samples$Genomic=="Yes",]
info_samples.1=info_samples[,c("sample","Genomic")]

annot.1=merge(annot,info_samples.1,by.x="sampleID",by.y="sample",all.x=T)
annot.1=annot.1[complete.cases(annot.1$sample),]
dim(annot.1) # 1844 17
annot.1$sampleID=annot.1$sample
unique(annot.1$sampleID) # 30
annot=annot.1

annot=annot[,-ncol(annot)]

annot_dis<- annot[annot$gene =="DIS3",]
dim(annot_dis) # 4 16

annot_dis=annot_dis[!(is.na(annot_dis$sample)),]
# D479, D488, and R780

annot_dis$cristal<- 0  
annot_dis$cristal[grep("D479", annot_dis$aachange)]<- "D479"
annot_dis$cristal[grep("D488", annot_dis$aachange)]<- "D488"
annot_dis$cristal[grep("R780", annot_dis$aachange)]<- "R780"
annot_dis$cristal[annot_dis$impact=="Synonymous"]<-"Synonymous"
colnames(annot_dis)[1]<-"sample"
annot_dis2<- annot_dis[,c("sample","cristal", "impact","gene")]
annot_dis2$sample<- gsub("_1_BM","", annot_dis2$sample)
length(unique(annot_dis2$sample)) # 4
clinc=read.delim(INPUT,
                    stringsAsFactors= F, header=T)
length(unique(clinc$sample[!is.na(clinc$SNV_DIS3)])) # 30
length(unique(clinc$sample[clinc$SNV_DIS3==1])) # 4
kk<-clinc$sample[clinc$SNV_DIS3==1][! clinc$sample[clinc$SNV_DIS3==1] %in% annot_dis2$sample]
kk[!is.na(kk)]
annot_dis2$sample[!annot_dis2$sample %in% clinc$sample[clinc$SNV_DIS3==1]]
int_clin<- merge(annot_dis2, clinc,by="sample")
int_clin$cristal_code<-1
int_clin$cristal_code[int_clin$cristal==0]<- 0
hs.dis3.samples=unique(int_clin[int_clin$cristal_code==1,"sample"])
################################################################################
# . . . Add the DIS3 information in the matrix
all_gen$SNV.HS_DIS3=ifelse(all_gen$sample %in% hs.dis3.samples,1,
                           ifelse(is.na(all_gen$SNV_DIS3),NA,0))
table(all_gen$SNV_DIS3,all_gen$SNV.HS_DIS3)

all_gen$SNV.noHS_DIS3=ifelse(all_gen$sample %in% hs.dis3.samples,0,all_gen$SNV_DIS3)
table(all_gen$SNV_DIS3,all_gen$SNV.noHS_DIS3)

# 
# 
all_gen$clusters_fra<- 0
all_gen$clusters_fra[all_gen$HRD==1 & all_gen$SNV_KRAS>0 & all_gen$add_gains==0 & all_gen$deletions==0 & all_gen$t_MMSET==0]<- "A_HRD_RAS"
all_gen$clusters_fra[all_gen$HRD==1 & all_gen$SNV_NRAS>0& all_gen$add_gains==0 & all_gen$deletions==0 & all_gen$t_MMSET==0]<- "A_HRD_RAS"
all_gen$clusters_fra[all_gen$HRD==1 & all_gen$SNV_BRAF>0& all_gen$add_gains==0 & all_gen$deletions==0 & all_gen$t_MMSET==0]<- "A_HRD_RAS"
all_gen$clusters_fra[all_gen$HRD==1 & all_gen$SNV_FGFR3>0& all_gen$add_gains==0 & all_gen$deletions==0 & all_gen$t_MMSET==0]<- "A_HRD_RAS"

all_gen$clusters_fra[all_gen$HRD==1 &
                       all_gen$deletions==0 & all_gen$add_gains==1 & all_gen$clusters_fra==0 & all_gen$t_MMSET==0]<- "B_HRD_gains"
all_gen$clusters_fra[all_gen$HRD==1 &
                       all_gen$deletions==1 & all_gen$add_gains==0 & all_gen$clusters_fra==0 & all_gen$t_MMSET==0]<- "C_HRD_deletions"

all_gen$clusters_fra[all_gen$HRD==1 & all_gen$clusters_fra==0 & all_gen$t_MMSET==0]<- "C_HRD_deletions"
# 

table(all_gen$clusters_fra, all_gen$HRD)

all_gen$clusters_fra[all_gen$t_MAF>0 & all_gen$SNV.HS_DIS3==0 & all_gen$APOBEC>1]<-"I_MAF_APOBEC"
all_gen$clusters_fra[all_gen$t_MAF>0 & all_gen$SNV.HS_DIS3==0 & all_gen$APOBEC<2]<-"I_MAF_APOBEC"

all_gen$clusters_fra[all_gen$t_CCND1>0 & all_gen$HRD==1 & all_gen$deletions==0]<-"E_CCND1_HRD"
all_gen$clusters_fra[all_gen$t_CCND1>0 & all_gen$SNV_KRAS>0 & all_gen$deletions==0]<-"F_CCND1_RAS"
all_gen$clusters_fra[all_gen$t_CCND1>0 & all_gen$SNV_BRAF>0 & all_gen$deletions==0]<-"F_CCND1_RAS"
all_gen$clusters_fra[all_gen$t_CCND1>0 & all_gen$SNV_NRAS>0 & all_gen$deletions==0]<-"F_CCND1_RAS"
all_gen$clusters_fra[all_gen$t_CCND1>0 & all_gen$SNV_FGFR3>0 & all_gen$deletions==0]<-"F_CCND1_RAS"
all_gen$clusters_fra[all_gen$t_CCND1>0 & all_gen$deletions>0 ]<-"G_CCND1_complex"
all_gen$clusters_fra[all_gen$t_CCND1>0 & all_gen$clusters_fra==0 & all_gen$deletions==0]<-"FF_CCND1_RAS"
all_gen$clusters_fra[all_gen$t_CCND1>0 & all_gen$clusters_fra %in% c("E_CCND1_HRD","F_CCND1_RAS","FF_CCND1_RAS") & all_gen$x1q>=1 ]<-"G_CCND1_complex"
all_gen$clusters_fra[all_gen$t_CCND1>0 & all_gen$clusters_fra %in% c("E_CCND1_HRD","F_CCND1_RAS","FF_CCND1_RAS") & all_gen$CNV.SNV_TP53>1 ]<-"G_CCND1_complex"

table(all_gen$clusters_fra)

all_gen$clusters_fra[all_gen$t_MMSET>0 & all_gen$x1q>0 & (all_gen$CNV.SNV_RB1>0 | all_gen$CNV.SNV_TGDS>0)  & all_gen$SNV.HS_DIS3==0 & all_gen$clusters_fra==0]<-"N_NSD2_1q_del13q"
all_gen$clusters_fra[all_gen$t_MMSET>0 & all_gen$x1q==0 & (all_gen$CNV.SNV_RB1>0 | all_gen$CNV.SNV_TGDS>0)  & all_gen$SNV.HS_DIS3==0 & all_gen$clusters_fra==0]<-"P_NSD2_del13q"
all_gen$clusters_fra[all_gen$t_MMSET==0 & all_gen$x1q>0 & (all_gen$CNV.SNV_RB1>0 | all_gen$CNV.SNV_TGDS>0) & all_gen$SNV.HS_DIS3==0 & all_gen$clusters_fra==0]<-"Q_x1q_del13q"

#
all_gen$clusters_fra[all_gen$APOBEC>1 & all_gen$t_MAF>0 & all_gen$clusters_fra == 0]<-"L_MAF_APOBEC"
all_gen$clusters_fra[all_gen$APOBEC>1 & all_gen$t_MMSET==0 & all_gen$clusters_fra == 0]<-"L_MAF_APOBEC"
all_gen$clusters_fra[all_gen$APOBEC>1 & all_gen$t_CCND1==0 & all_gen$clusters_fra == 0]<-"L_MAF_APOBEC"

all_gen$clusters_fra[all_gen$APOBEC>1 & all_gen$clusters_fra == "Q_x1q_del13q"]<-"L_MAF_APOBEC"
all_gen$clusters_fra[all_gen$APOBEC>1 & 
                       (all_gen$clusters_fra == "A_HRD_RAS" | all_gen$clusters_fra %in% c("B_HRD_gains","BB_HRD_gains") | all_gen$clusters_fra %in% c("C_HRD_deletions","CC_HRD_deletions") | all_gen$clusters_fra %in% c("D_HRD_other","DD_HRD_other"))]<-"M_MAF_APOBEC"

all_gen$clusters_fra[all_gen$t_MMSET>0 & all_gen$clusters_fra == 0]<-"O_NSD2"

table(all_gen[all_gen$clusters_fra==0,"deletions"])
all_gen$clusters_fra[all_gen$clusters_fra==0 & all_gen$deletions==1]<-"R_Multiple_losses"
all_gen[all_gen$clusters_fra==0,]
all_gen$clusters_fra[all_gen$clusters_fra==0]<-"RR_Simple"
table(all_gen$clusters_fra)

all_gen$clusters_fra[all_gen$clusters_fra==0 & all_gen$t_CCND1>0]<-"HH_CCND1"

all_gen$clusters_fra[all_gen$clusters_fra %in% c("E_CCND1_HRD","F_CCND1_RAS","FF_CCND1_RAS") &
                       (all_gen$APOBEC>1 | all_gen$CNV.Sig==1)]<-"G_CCND1_complex"


all_gen$clusters_fra[all_gen$clusters_fra %in% c("M_MAF_APOBEC","I_MAF_APOBEC","L_MAF_APOBEC") &
                       (all_gen$t_MMSET>0 & all_gen$t_MAF>0)]<-"N_NSD2_1q_del13q"

all_gen$clusters_fra[all_gen$clusters_fra %in% c("RR_Simple") &
                       (all_gen$CNV.Sig==1 | all_gen$x1q==2)]<-"R_Multiple_losses"


table(all_gen$clusters_fra)
all_gen<- all_gen[order(all_gen$clusters_fra),]
table(all_gen$clusters_fra)

colnames(all_gen)=gsub("CNV_","",colnames(all_gen))
colnames(all_gen)=gsub("CNV.SNV_","",colnames(all_gen))
colnames(all_gen)=gsub("SNV_","",colnames(all_gen))
colnames(all_gen)=gsub("SNV.","",colnames(all_gen))


colnames(all_gen)[96]="CCND1-IGH"
colnames(all_gen)[97]="NSD2-IGH"
colnames(all_gen)[98]="MAF-IGH"


# tc6=as.data.frame(readxl::read_excel("~/Google Drive/My Drive/Bachisio/Share/MM_PredictionModel/Heatmap/files/MGP_TC6.xlsx"))
# tc6=tc6[tc6$Study=="MMRF",]
# tc6=tc6[,c("Sample_Name","TC6_calls","RNA_high_risk")]
# tc6$Sample_Name=gsub("_1_BM","",tc6$Sample_Name)
# colnames(tc6)[1]="sample"
# tc6$TC6_calls=ifelse(is.na(tc6$TC6_calls),"NA",tc6$TC6_calls)
# tc6$RNA_high_risk=ifelse(is.na(tc6$RNA_high_risk),"NA",tc6$RNA_high_risk)
# 
# uams.c=read.csv("~/Google Drive/My Drive/Bachisio/Documents/Project/MM/Prediction_Model/analysis/UAMS/Call Table.csv",
#                 stringsAsFactors = F)
# uams.c.1=uams.c[,c("SAMPLE","UAMS.CALL")]
# colnames(uams.c.1)=c("sample","UAMS.call")
# uams.c.1=uams.c.1[grep("_1_BM",uams.c.1$sample),]
# uams.c.1$sample=gsub("_1_BM","",uams.c.1$sample)
# 
# 
# all_gen1=merge(all_gen,tc6,by="sample",all.x=T)
# all_gen1=merge(all_gen1,uams.c.1,by="sample",all.x=T)
dim(all_gen) # 30 105
all_gen1=all_gen[,c("sample", #samples
                     "chr3.gain","chr5.gain","chr7.gain","chr9.gain","chr11.gain","chr15.gain","chr19.gain","chr21.gain", #HRD
                     "Gain_Amp1q", #1q
                     "chr18.gain","Amp_14q11.2","Amp_17q22","Amp_2q24.3","Amp_2q32.1","Amp_4q13.2","Amp_6p24.3","Amp_8q24.21", #chr18.gain + other gains
                     "ARID1A","CDKN2C","FAM46C","FUBP1","RPL5", #chr1
                     "Del_2q31.1","Del_2q37.3","DNMT3A","SP140", #chr2
                     "RASA2","SETD2", #chr3
                     "TET2", #chr4
                     "Del_6q26","PRDM1","ZNF292", #chr6
                     "Del_7p11.2","Del_7p22.2","KMT2C","POT1", #chr7
                     "chr8p.loss","Del_8q24.21","UBR5", # chr8
                     "Del_9p21.3","TRAF2", #chr9
                     "Del_10p15.3","Del_10q24.32","Del_10q26.3", #chr10
                     "Del_11q22.1","ATM", #chr11
                     "Del_12p13.2","Del_12q24.31","ARID2","BTG1","CDKN1B", #chr12
                     "RB1","TGDS", #chr13
                     "Del_14q24.3","MAX","NFKBIA","TRAF3", #chr14
                     "CREBBP","CYLD", #chr16
                     "NCOR1","NF1","TP53", #chr17
                     "KMT2B", #chr19
                     "Del_20p13","Del_20q13.12", #chr20
                     "Del_21p11.2", #chr21
                     "Del_22q13.32","EP300", #chr22
                     "BRAF","KRAS","NRAS", #BRAF,KRAS,NRAS # (MAPK pathway)
                     "CCND1","FGFR3","HIST1H1E","IRF4","NFKB2","PTPN11","SF3B1","ACTG1","DUSP2","HIST1H1C","HLA.C",
                     "HUWE1","KLHL6","LTB","PABPC1","PRKD2","SAMHD1","BCL7A","EGR1","XBP1",#other mutations  # N.B.: Remove DIS3
                     "HS_DIS3","noHS_DIS3", # DIS3 HS,DIS3 noHS
                     "CCND1-IGH","NSD2-IGH","MAF-IGH", #translocations
                     "APOBEC", #APOBEC
                     "CNV.Sig", #CNV.Sig
                     "HRD","clusters_fra"#"TC6_calls","RNA_high_risk","UAMS.call" [These classification are not present in this code]
)]
dim(all_gen1) # 30 101


annotation_col_genomic<- as.data.frame((all_gen1[,c("clusters_fra")]))
rownames(annotation_col_genomic)<-all_gen1$sample
colnames(annotation_col_genomic)<-c("clusters")
# annotation_col_genomic$TC6=ifelse(is.na(annotation_col_genomic$TC6),"NA",annotation_col_genomic$TC6)
# annotation_col_genomic$GEP70=ifelse(is.na(annotation_col_genomic$GEP70),"NA",annotation_col_genomic$GEP70)
# annotation_col_genomic$UAMS7.ms=ifelse(is.na(annotation_col_genomic$UAMS7.ms),"NA",annotation_col_genomic$UAMS7.ms)
# 
# annotation_col_genomic$UAMS7.ms=ifelse(annotation_col_genomic$TC6=="NA","NA",annotation_col_genomic$UAMS7.ms)

annotation_col_genomic$clusters.new=ifelse(annotation_col_genomic$clusters %in% c("A_HRD_RAS"),"HRD_RAS",
                                           ifelse(annotation_col_genomic$clusters %in% c("B_HRD_gains"),"HRD_Gains",
                                                  ifelse(annotation_col_genomic$clusters %in% c("C_HRD_deletions"),"HRD_Complex_Cytogenetic",
                                                         ifelse(annotation_col_genomic$clusters %in% c("E_CCND1_HRD","F_CCND1_RAS","FF_CCND1_RAS"),"CCND1_Simple",
                                                                ifelse(annotation_col_genomic$clusters %in% c("G_CCND1_complex"),"CCND1_Complex_Cytogenetic",
                                                                       ifelse(annotation_col_genomic$clusters %in% c("I_MAF_APOBEC","L_MAF_APOBEC","M_MAF_APOBEC"),"MAF_and_or_HyperAPOBEC",
                                                                              ifelse(annotation_col_genomic$clusters %in% c("N_NSD2_1q_del13q"),"NSD2_GainAmp1q_Del13q",
                                                                                     ifelse(annotation_col_genomic$clusters %in% c("O_NSD2"),"NSD2_HRD",
                                                                                            ifelse(annotation_col_genomic$clusters %in% c("P_NSD2_del13q"),"NSD2_Del13q",
                                                                                                   ifelse(annotation_col_genomic$clusters %in% c("Q_x1q_del13q"),"GainAmp1q_Del13q",
                                                                                                          ifelse(annotation_col_genomic$clusters %in% c("R_Multiple_losses"),"Multiple_Losses",
                                                                                                                 ifelse(annotation_col_genomic$clusters %in% c("RR_Simple"),"Simple",NA))))))))))))


mycol_plus<- c(brewer.pal(8,"Dark2"),brewer.pal(8,"Paired"), brewer.pal(8,"Set2"))
ann_colors = list(
              #clusters =c(
              # "A_HRD_RAS"="lightpink3",#mycol_plus[1],
              # "B_HRD_gains"="mediumpurple2",#mycol_plus[3],
              # "C_HRD_deletions"="darkorchid4",#mycol_plus[2],
              # "E_CCND1_HRD"="skyblue3",#mycol_plus[5],
              # "F_CCND1_RAS"="skyblue3",#mycol_plus[8],
              # "FF_CCND1_RAS"="skyblue3",#mycol_plus[6],
              # "G_CCND1_complex"="navy",#mycol_plus[9],
              # "I_MAF_APOBEC"="darkgreen",#mycol_plus[13],
              # "L_MAF_APOBEC"="darkgreen",#mycol_plus[14],
              # "M_MAF_APOBEC"="darkgreen",#mycol_plus[7],
              # "N_NSD2_1q_del13q"="lightgoldenrod",#mycol_plus[10],
              # "O_NSD2"="gold",#mycol_plus[16],
              # "P_NSD2_del13q"="goldenrod4",#mycol_plus[11],
              # "Q_x1q_del13q"="sienna3",#mycol_plus[12],
              # "R_Multiple_losses"="firebrick",#mycol_plus[15],
              # "RR_Simple"="slategray"#mycol_plus[22]
              clusters.new =c("HRD_RAS"="lightpink3",#mycol_plus[1],
                          "HRD_Gains"="mediumpurple2",#mycol_plus[3],
                          "HRD_Complex_Cytogenetic"="darkorchid4",#mycol_plus[2],
                          "CCND1_Simple"="skyblue3",#mycol_plus[5],
                          "CCND1_Complex_Cytogenetic"="navy",#mycol_plus[9],
                          "MAF_and_or_HyperAPOBEC"="darkgreen",#mycol_plus[13],
                          "NSD2_GainAmp1q_Del13q"="lightgoldenrod",#mycol_plus[10],
                          "NSD2_HRD"="gold",#mycol_plus[16],
                          "NSD2_Del13q"="goldenrod4",#mycol_plus[11],
                          "GainAmp1q_Del13q"="sienna3",#mycol_plus[12],
                          "Multiple_Losses"="firebrick",#mycol_plus[15],
                          "Simple"="slategray"#mycol_plus[22]
  )#,
  # TC6=c(   "MMSET"="tan1",
  #          "MAF"="olivedrab3",
  #          "MAFA"="limegreen",
  #          "MAFB"="springgreen4",
  #          "CCND1"="royalblue3",
  #          "CCND3"="turquoise2",
  #          "D1"="maroon4",
  #          "D2"="lightpink",
  #          "NA"="floralwhite"),
  # GEP70=c("TRUE"="lightcoral",
  #         "FALSE"="lightblue1",
  #         "NA"="floralwhite"),
  # UAMS7.ms=c("MS"="mediumseagreen",
  #            "MF"="mediumturquoise",
  #            "LB"="lightsalmon",
  #            "HY"="magenta2",
  #            "PR"="black",
  #            "CD1"="lightskyblue",
  #            "CD2"="lightslateblue",
  #            "NA"="floralwhite")
  )


space_heat<- as.numeric(table(all_gen1$clusters_fra))
length(space_heat)

# space_heat2<- c(space_heat[1],
#                 space_heat[1] + space_heat[2],
#                 space_heat[1] + space_heat[2] + space_heat[3],
#                 space_heat[1] + space_heat[2] + space_heat[3] + (space_heat[4]+space_heat[5]+space_heat[6]),
#                 space_heat[1] + space_heat[2] + space_heat[3] + (space_heat[4]+space_heat[5]+space_heat[6]) + space_heat[7],
#                 space_heat[1] + space_heat[2] + space_heat[3] + (space_heat[4]+space_heat[5]+space_heat[6]) + space_heat[7] + (space_heat[8]+space_heat[9]+space_heat[10]),
#                 space_heat[1] + space_heat[2] + space_heat[3] + (space_heat[4]+space_heat[5]+space_heat[6]) + space_heat[7] + (space_heat[8]+space_heat[9]+space_heat[10]) + space_heat[11],
#                 space_heat[1] + space_heat[2] + space_heat[3] + (space_heat[4]+space_heat[5]+space_heat[6]) + space_heat[7] + (space_heat[8]+space_heat[9]+space_heat[10]) + space_heat[11] + space_heat[12],
#                 space_heat[1] + space_heat[2] + space_heat[3] + (space_heat[4]+space_heat[5]+space_heat[6]) + space_heat[7] + (space_heat[8]+space_heat[9]+space_heat[10]) + space_heat[11] + space_heat[12] + space_heat[13],
#                 space_heat[1] + space_heat[2] + space_heat[3] + (space_heat[4]+space_heat[5]+space_heat[6]) + space_heat[7] + (space_heat[8]+space_heat[9]+space_heat[10]) + space_heat[11] + space_heat[12] + space_heat[13] + space_heat[14],
#                 space_heat[1] + space_heat[2] + space_heat[3] + (space_heat[4]+space_heat[5]+space_heat[6]) + space_heat[7] + (space_heat[8]+space_heat[9]+space_heat[10]) + space_heat[11] + space_heat[12] + space_heat[13] + space_heat[14] + space_heat[15],
#                 space_heat[1] + space_heat[2] + space_heat[3] + (space_heat[4]+space_heat[5]+space_heat[6]) + space_heat[7] + (space_heat[8]+space_heat[9]+space_heat[10]) + space_heat[11] + space_heat[12] + space_heat[13] + space_heat[14] + space_heat[15] + space_heat[16]
# )

space_heat2<- c(space_heat[1],
                space_heat[1] + space_heat[2],
                space_heat[1] + space_heat[2] + space_heat[3],
                space_heat[1] + space_heat[2] + space_heat[3] + space_heat[4],
                space_heat[1] + space_heat[2] + space_heat[3] + space_heat[4] + space_heat[5],
                space_heat[1] + space_heat[2] + space_heat[3] + space_heat[4] + space_heat[5] + c(space_heat[6] + space_heat[7]),
                space_heat[1] + space_heat[2] + space_heat[3] + space_heat[4] + space_heat[5] + c(space_heat[6] + space_heat[7]) + space_heat[8],
                space_heat[1] + space_heat[2] + space_heat[3] + space_heat[4] + space_heat[5] + c(space_heat[6] + space_heat[7]) + space_heat[8] + space_heat[9]
)

space.row=c(8,
            8 + 1,
            8 + 1 + 8,
            8 + 1 + 8 + 51,
            8 + 1 + 8 + 51 + 3,
            8 + 1 + 8 + 51 + 3 + 22,
            8 + 1 + 8 + 51 + 3 + 22 + 3,
            8 + 1 + 8 + 51 + 3 + 22 + 3 + 1,
            8 + 1 + 8 + 51 + 3 + 22 + 3 + 1 + 1)

annotation_col_genomic1=as.data.frame(annotation_col_genomic[,c("clusters.new")])
colnames(annotation_col_genomic1)="clusters.new"
rownames(annotation_col_genomic1)=rownames(annotation_col_genomic)
rownames(all_gen1)=all_gen1$sample
all_gen1=all_gen1[,-1]
mtx=all_gen1[,c(1:98)]
identical(rownames(mtx),rownames(annotation_col_genomic))
annotation_col_genomic=annotation_col_genomic[order(annotation_col_genomic$clusters),]
mtx=mtx[match(rownames(annotation_col_genomic),rownames(mtx)),]
identical(rownames(mtx),rownames(annotation_col_genomic))

# space_heat2<- space_heat2[!is.na(space_heat2)]
pheatmap(as.matrix(t(mtx)), annotation_col= annotation_col_genomic1, 
         annotation_colors = ann_colors,
         cluster_cols = FALSE, show_colnames = F, fontsize_row = 5,
         cluster_rows =  FALSE, border_color = F, legend = F, 
         # col=c("grey80","white","gold3","forestgreen","dodgerblue","darkorchid1","red"),
         col=c("grey85","salmon1","darkred"),
         gaps_col = space_heat2,
         gaps_row = space.row
)
