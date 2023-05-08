##################################################################################
##################################################################################
##############################################################################1###
#                                                                                #
#        .     .     .                                               .     .     .
####     .     .     .   PredictionModel - Matrix Generation   .     .     .     ####
#        .     .     .                                               .     .     .
#                                                                                #
##################################################################################
##################################################################################
##################################################################################



# * * * * * * * * * *         IMPORTANT         * * * * * * * * * * 
#
#                          Change the folder location (WD: working directory)
#
# * * * * * * * * * *         IMPORTANT         * * * * * * * * * * 



WF="~/Google Drive/My Drive/Bachisio/Documents/Project/MM/Prediction_Model/GCP_MM/genomic/" # /where/is/the/folder/genomic 

##################################################################################
##################################################################################
##################################################################################
##################################################################################

INPUT.FILE.FOLDER=paste(WF,"/output/",sep="")
CLINIC.FILE=paste(WF,"/input/clinic_info_test.txt",sep="")
CNV.SIG=paste(WF,"/input/commpass_CNV.Sig.txt",sep="")

OUTPUT=paste(WF,"/matrixDEF_test.txt",sep="")

##############################################################################################################################################################
######
# . . . . .                                    PredictionModel - Matrix Creation 
######
################################################################################################################################################################

library(dplyr)
library(stringr)


clinic=read.delim(CLINIC.FILE,stringsAsFactors = F)

clinic=clinic[,c(1:26)]

#######
#### . . .    - HRD Chromosomes -     . . . ####
#######
hrd=read.delim(paste(INPUT.FILE.FOLDER,"HRD_chromosomes.txt",sep=""),stringsAsFactors = F)
colnames(hrd)=paste(colnames(hrd),"gain",sep="")
unique(rownames(hrd))

rownames(hrd)=gsub("_1_BM","",rownames(hrd))

# . - Check samples
intersect(clinic$sample,rownames(hrd))
setdiff(clinic$sample,rownames(hrd))
setdiff(rownames(hrd),clinic$sample)
colnames(hrd)=paste("CNV_",colnames(hrd),sep="")
setdiff(clinic$sample,rownames(hrd))
# . - Merge CLINIC info with HRD info 
clinic1=merge(clinic,hrd,by.x="sample","by.y"="row.names",all.x=T)



#######
#### . . .    - 1q -     . . . ####
#######
x1q=read.delim(paste(INPUT.FILE.FOLDER,"x1q_amp_gain.txt",sep=""),stringsAsFactors = F)
unique(rownames(x1q))

rownames(x1q)=gsub("_1_BM","",rownames(x1q))

# . - Check samples
intersect(clinic$sample,rownames(x1q))
setdiff(clinic$sample,rownames(x1q))
setdiff(rownames(x1q),clinic$sample)
setdiff(clinic$sample,rownames(x1q))
head(x1q)
colnames(x1q)=paste("CNV_",colnames(x1q),sep="")
# . - Merge CLINIC info with 1q info 
clinic1=merge(clinic1,x1q,by.x="sample","by.y"="row.names",all.x=T)
clinic1=clinic1[,-ncol(clinic1)]


#######
#### . . .    - GISTIC peaks -     . . . ####
#######
gistic=read.delim(paste(INPUT.FILE.FOLDER,"gistic.other_All_peaks.txt",sep=""),stringsAsFactors = F)
unique(rownames(gistic))
# Remove the GISTIC peaks with TSGs
# Del_1p32.3,Del_1p12,Del_1p22.1,Del_13q14.2,Del_16p13.3,Del_16q12.1,Del_17p13.1
colnames(gistic)
gistic=gistic[,!colnames(gistic) %in% c("Del_1p32.3","Del_1p12","Del_1p22.1","Del_13q14.2","Del_16p13.3","Del_16q12.1","Del_17p13.1",
                                        "Amp_14q32.33","Del_14q32.33","Del_22q11.22","Del_2p11.2","Amp_2p11.2")]
colnames(gistic)
rownames(gistic)=gsub("_1_BM","",rownames(gistic))
gistic=gistic[,-ncol(gistic)]

for(i in 1:ncol(gistic)){
  gistic[,i]=ifelse(gistic[,i]>1,1,gistic[,i])
}

# . - Check samples
intersect(clinic$sample,rownames(gistic))
setdiff(clinic$sample,rownames(gistic))
setdiff(rownames(gistic),clinic$sample)
setdiff(clinic$sample,rownames(gistic))
colnames(gistic)=paste("CNV_",colnames(gistic),sep="")

# . - Merge CLINIC info with HRD info 
clinic2=merge(clinic1,gistic,by.x="sample","by.y"="row.names",all.x=T)
head(clinic2)



#######
#### . . .    - ONCOGENEs -     . . . ####
#######
oncogenes=read.delim(paste(INPUT.FILE.FOLDER,"oncogene_commpass.txt",sep=""),stringsAsFactors = F)


# . - Check samples
intersect(clinic2$sample,rownames(oncogenes))
setdiff(clinic2$sample,rownames(oncogenes))
setdiff(rownames(oncogenes),clinic2$sample)
setdiff(clinic2$sample,rownames(oncogenes))
for(i in 1:ncol(oncogenes)){
  oncogenes[,i]=ifelse(oncogenes[,i]>1,1,oncogenes[,i])
}
colnames(oncogenes)=paste("SNV_",colnames(oncogenes),sep="")

# . - Merge CLINIC info with HRD info 
clinic3=merge(clinic2,oncogenes,by.x="sample","by.y"="row.names",all.x=T)



#######
#### . . .    - TSGs -     . . . ####
#######
tsg=read.delim(paste(INPUT.FILE.FOLDER,"tsg_commpass_test.txt",sep=""),stringsAsFactors = F)


# . - Check samples
intersect(clinic3$sample,rownames(tsg))
setdiff(clinic3$sample,rownames(tsg))
setdiff(rownames(tsg),clinic3$sample)
setdiff(clinic3$sample,rownames(tsg))
colnames(tsg)=paste("CNV.SNV_",colnames(tsg),sep="")

# . - Merge CLINIC info with HRD info 
clinic4=merge(clinic3,tsg,by.x="sample","by.y"="row.names",all.x=T)


#######
#### . . .    - UNKNOWNs -     . . . ####
#######
unk=read.delim(paste(INPUT.FILE.FOLDER,"unk_mutated_genes_commpass_test.txt",sep=""),stringsAsFactors = F)

# . - Check samples
intersect(clinic3$sample,rownames(unk))
setdiff(clinic3$sample,rownames(unk))
setdiff(rownames(unk),clinic3$sample)
setdiff(clinic4$sample,rownames(unk))
colnames(unk)=paste("SNV_",colnames(unk),sep="")

# . - Merge CLINIC info with HRD info 
clinic5=merge(clinic4,unk,by.x="sample","by.y"="row.names",all.x=T)
#clinic3=clinic3[,-c(ncol(clinic3))]


#######
#### . . .    - APOBEC -     . . . ####
#######
apobec=read.delim(paste(INPUT.FILE.FOLDER,"APOBEC_relative_contribution_matrix.txt",sep=""),stringsAsFactors = F)
rownames(apobec)=gsub("_1_BM","",rownames(apobec))

setdiff(clinic5$sample,rownames(apobec))

# . - Merge CLINIC info with APOBEC info 
clinic6=merge(clinic5,apobec,by.x="sample","by.y"="row.names",all.x=T)
dim(clinic6)
# remove offset column
clinic6=clinic6[,-ncol(clinic6)]


#######
#### . . .    - APOBEC presence or absence -     . . . ####
#######
deciles=quantile(na.omit(clinic6$apobec), probs = seq(0, 1, 1/10),)

# . . . - These values are the deciles obtained using the whole cohort.
#         0%        10%        20%        30%        40%        50%        60%        70%        80%        90%       100% 
# 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.02891755 0.04139446 0.06430805 0.11203844 0.91384142 
table(is.na(clinic6$apobec))
clinic6$APOBEC=ifelse(clinic6$apobec<deciles[7],0,
                      ifelse(clinic6$apobec>=deciles[7] & clinic6$apobec<deciles[10],1,
                             ifelse(clinic6$apobec>=deciles[10],2,NA)))

table(clinic6$APOBEC,useNA="ifany")

head(clinic6)



#######
#### . . .    - TRANSLOCATION -     . . . ####
#######
tra=read.delim(paste(INPUT.FILE.FOLDER,"translocations.txt",sep=""),stringsAsFactors = F)
rownames(tra)=gsub("_1_BM","",rownames(tra))

setdiff(clinic6$sample,rownames(tra))

# . - Merge CLINIC info with translocation info 
clinic7=merge(clinic6,tra,by.x="sample","by.y"="row.names",all.x=T)
dim(clinic7)
head(clinic7)



#######
#### . . .    - ChrX genes -     . . . ####
#######
head(clinic7)
clinic7=clinic7[,!colnames(clinic7) %in% c("CNV.SNV_ATRX","CNV.SNV_KDM5C","CNV.SNV_KDM6A","CNV.SNV_RPL10")]

chrX=read.delim(paste(INPUT.FILE.FOLDER,"chrx_tsg_commpass_test.txt",sep=""),
                stringsAsFactors = F)
colnames(chrX)=paste("SNV_",colnames(chrX),sep="")

clinic8=merge(clinic7,chrX,by.x="sample",by.y="row.names",all.x=T)


#######
#### . . .    - Matrix -     . . . ####
#######

clinic8=clinic8[,-137] # IGLL5

clinic9=clinic8[,c(1:26, # Clinic
                   27:70, # copy number
                   93:124, # TSGs
                   71:92,125:157,163:166, # Oncogenes and Unknown
                   159, # APOBEC
                   160:162 # translocation
                   )]
#clinic9$t_MAF=ifelse(clinic9$t_MAF>0,2,clinic9$t_MAF)

#### Collapse 1q peaks - . . . ####
clinic9$Gain_Amp1q=paste(clinic9$CNV_Amp_1q21.1,clinic9$CNV_Amp_1q21.3,clinic9$CNV_Amp_1q44,sep="_")
table(clinic9$Gain_Amp1q)
clinic9$Gain_Amp1q=ifelse(clinic9$Gain_Amp1q=="0_0_0",0,
                 ifelse(clinic9$Gain_Amp1q=="NA_NA_NA",NA,
                        ifelse(clinic9$Gain_Amp1q %in% c("1_1_1","1_1_0","0_1_0","0_1_1","0_0_1","1_0_0","1_0_1"),1,
                               ifelse(clinic9$Gain_Amp1q %in% c("1_2_1","0_2_0","2_1_1","0_2_1","1_2_0","2_2_0","2_2_2","2_2_1","0_1_2",
                                                        "2_1_0","1_1_2","1_2_2","0_0_2","0_2_2","2_0_0"),2,NA))))
colnames(clinic9)[166]="CNV_Gain_Amp1q"
#### . . . - Collapse 8p - . . . ####
clinic9$CNV_chr8p.loss=paste(clinic9$CNV_Del_8p21.3,clinic9$CNV_Del_8p23.1,sep="_")
table(clinic9$CNV_chr8p.loss,useNA="ifany")

clinic9$CNV_chr8p.loss=ifelse(clinic9$CNV_chr8p.loss=="0_0",0,
                               ifelse(clinic9$CNV_chr8p.loss=="NA_NA",NA,
                                      ifelse(clinic9$CNV_chr8p.loss %in% c("1_1","1_0","0_1"),1,NA)))

#### . . . - CNVSig - . . . ####
cnv.sig=read.delim(CNV.SIG,stringsAsFactors = F)
clinic9=merge(clinic9,cnv.sig,by="sample")


#### . . . - Column selection - . . . ####
clinic10=clinic9[,c("sample","age","gender","ecog","ISS","SCT_first_line","time_SCT",      
                    "SCT_line","pfs_time","pfs_code","os_time","os_code","KAR","chemo",           
                    "BORT","LEN","THAL","DARA","ELO","duration","continuos_treat",
                    "combo","study","LDH_level","Race","phase","CNV_chr3.gain","CNV_chr5.gain",   
                    "CNV_chr7.gain","CNV_chr9.gain","CNV_chr11.gain","CNV_chr15.gain","CNV_chr19.gain","CNV_chr21.gain","CNV_chr18.gain",  
                    "CNV_Gain_Amp1q","CNV_Amp_14q11.2","CNV_Amp_17q22","CNV_Amp_2q24.3","CNV_Amp_2q32.1","CNV_Amp_4q13.2","CNV_Amp_6p24.3",  
                    "CNV_Amp_8q24.21","CNV_chr8p.loss","CNV_Del_10p15.3","CNV_Del_10q24.32","CNV_Del_10q26.3","CNV_Del_11q22.1","CNV_Del_12p13.2", 
                    "CNV_Del_12q24.31","CNV_Del_14q24.3","CNV_Del_20p13","CNV_Del_20q13.12","CNV_Del_21p11.2","CNV_Del_22q13.32","CNV_Del_2q31.1",  
                    "CNV_Del_2q37.3","CNV_Del_6q26","CNV_Del_7p11.2","CNV_Del_7p22.2","CNV_Del_8q24.21","CNV_Del_9p21.3","CNV.SNV_ARID1A",  
                    "CNV.SNV_ARID2","CNV.SNV_ATM","CNV.SNV_BTG1","CNV.SNV_CDKN1B","CNV.SNV_CDKN2C","CNV.SNV_CREBBP","CNV.SNV_CYLD",    
                    "CNV.SNV_DNMT3A","CNV.SNV_EP300","CNV.SNV_FAM46C","CNV.SNV_FUBP1","CNV.SNV_KMT2B","CNV.SNV_KMT2C","CNV.SNV_MAX",     
                    "CNV.SNV_NCOR1","CNV.SNV_NF1","CNV.SNV_NFKBIA","CNV.SNV_POT1","CNV.SNV_PRDM1","CNV.SNV_RASA2","CNV.SNV_RB1",     
                    "CNV.SNV_RPL5","CNV.SNV_SETD2","CNV.SNV_SP140","CNV.SNV_TET2","CNV.SNV_TGDS","CNV.SNV_TP53","CNV.SNV_TRAF2",   
                    "CNV.SNV_TRAF3","CNV.SNV_UBR5","CNV.SNV_ZNF292","SNV_BRAF","SNV_CCND1","SNV_DIS3","SNV_FGFR3",       
                    "SNV_HIST1H1E","SNV_IDH1","SNV_IDH2","SNV_IRF4","SNV_KRAS","SNV_NFKB2","SNV_NRAS",        
                    "SNV_PIK3CA","SNV_PTPN11","SNV_SF3B1","SNV_ABCF1","SNV_ACTG1","SNV_BHLHE41","SNV_DDX3X",       
                    "SNV_DTX1","SNV_DUSP2","SNV_HIST1H1B","SNV_HIST1H1C","SNV_HIST1H1D","SNV_HIST1H2BK","SNV_HLA.C",       
                    "SNV_HUWE1","SNV_INO80E","SNV_IRF1","SNV_KLHL6","SNV_LCE1D","SNV_LTB","SNV_MAN2C1",      
                    "SNV_MRGPRX4","SNV_NDUFC2","SNV_PABPC1","SNV_PANK3","SNV_PRKD2","SNV_RFTN1","SNV_RPRD1B",      
                    "SNV_RPS3A","SNV_SAMHD1","SNV_SLC35G5","SNV_SYT17","SNV_TBC1D29","SNV_TNFRSF11B","SNV_USP7",        
                    "SNV_ZFP36L1","SNV_ATRX","SNV_KDM5C","SNV_KDM6A","SNV_RPL10","SNV_BCL7A","SNV_EGR1",        
                    "SNV_MAF","SNV_MAFB","SNV_MAML2","SNV_PIM1","SNV_TCL1A","SNV_XBP1","CNV.Sig",         
                    "APOBEC","t_CCND1","t_MMSET","t_MAF")]



write.table(clinic10,
            OUTPUT,
            col.names = T,
            row.names = F,
            quote = F,
            sep = "\t")
