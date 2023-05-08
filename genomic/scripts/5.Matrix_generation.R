##################################################################################
##################################################################################
##############################################################################1###
#                                                                                #
#        .     .     .                                               .     .     .
####     .     .     .   PredictionModel - Matrix Generation   .     .     .     ####
#        .     .     .                                               .     .     .
#                                                                                #
#B###########################################################################Z####
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
CLINIC.FILE=paste(WF,"/input/commpass_info_test.txt",sep="")

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
#          0%        10%        20%        30%        40%        50%        60%        70%        80%        90%       100% 
#  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.02280722 0.04085643 0.06469131 0.11747495 0.98904383
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

clinic8=clinic8[,-36] # IGLL5

clinic9=clinic8[,c(1:26, # Clinic
                   27:68, # copy number
                   91:122, # TSGs
                   69:90,123:156,162:165, # Oncogenes and Unknown
                   158, # APOBEC
                   159:161 # translocation
                   )]
clinic9$t_MAF=ifelse(clinic9$t_MAF>0,2,clinic9$t_MAF)

write.table(clinic9,
            OUTPUT,
            col.names = T,
            row.names = F,
            quote = F,
            sep = "\t")


