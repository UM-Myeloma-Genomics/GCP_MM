##################################################################################
##################################################################################
##############################################################################1###
#                                                                                #
#        .     .     .                                               .     .     .
####     .     .     .   PredictionModel - Cytogenetic abnormalies   .     .     .     ####
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



###############################################################################################################################################################
######
# . . . . .                                    PredictionModel - Matrix HRD chromosome
######
################################################################################################################################################################
#####################
######## --- .                 
#####################

matrix=read.delim(paste(WF,"/output/CNV_matrix_focal_large_5Mb.txt",sep=""),
                  stringsAsFactors = F)
pres.abs=matrix
for(j in c(1:(ncol(matrix)-1))){
  pres.abs[,j]=ifelse(pres.abs[,j]>0,1,pres.abs[,j])
}
hyperdiploid=read.delim(paste(WF,"/output/60%_ARMS_hyperdiploid.txt",sep=""),
                        stringsAsFactors = F)
matrix1=merge(matrix,hyperdiploid,by.x="row.names",by.y="sample",all.x=T)
pres.abs1=merge(pres.abs,hyperdiploid,by.x="row.names",by.y="sample",all.x=T)
for(www in c(69:77)){
  pres.abs1[,www]=ifelse(is.na(pres.abs1[,www]),0,pres.abs1[,www])
}
pres.abs1$hyp.sample=ifelse(is.na(pres.abs1$hyp.sample),"No.HRD",pres.abs1$hyp.sample)
table(pres.abs1$hyp.sample)
matrix.hyper=pres.abs1[,c(1,69:76)]
rownames(matrix.hyper)=matrix.hyper$Row.names
matrix.hyper=matrix.hyper[,-1]
colnames(matrix.hyper)=c("chr3+","chr5+","chr7+","chr9+","chr11+","chr15+","chr19+","chr21+","chr18+)
write.table(matrix.hyper,
            paste(WF,"/output/HRD_chromosomes.txt",sep=""),
            col.names = T,
            row.names = T,
            quote = F,
            sep = "\t")
#

##############################################################################################################################################################
######
# . . . . .                                    PredictionModel - Matrix Other GISTIC peaks
######
################################################################################################################################################################
#####################
######## --- .                 
#####################

matrix=read.delim(paste(WF,"/output/CNV_matrix_focal_large_5Mb.txt",sep=""),
                  stringsAsFactors = F)

# gistic.amp=matrix[,c(4,5,8,16,17,18,20,24,25,27)]
# for(j in c(1:(ncol(gistic.amp)))){
#   gistic.amp[,j]=ifelse(gistic.amp[,j] %in% c(0.5,1),1,
#                         ifelse(gistic.amp[,j] %in% c(1.5,2),2,
#                                ifelse(gistic.amp[,j]==0,0,gistic.amp[,j])))
# }

gistic.amp=matrix[,c("Amp_14q11.2","Amp_14q32.33","Amp_17q22","Amp_2p11.2","Amp_2q24.3","Amp_2q32.1","Amp_4q13.2","Amp_6p24.3","Amp_6q11.1","Amp_8q24.21")]
for(j in c(1:(ncol(gistic.amp)))){
  gistic.amp[,j]=ifelse(gistic.amp[,j] %in% c(0.5,1),1,
                        ifelse(gistic.amp[,j] %in% c(1.5,2),2,
                               ifelse(gistic.amp[,j]==0,0,gistic.amp[,j])))
}



# gistic.del=matrix[,c(31:65)]
# for(j in c(1:(ncol(gistic.del)-1))){
#   gistic.del[,j]=ifelse(gistic.del[,j] %in% c(0.5,1),1,
#                         ifelse(gistic.del[,j] %in% c(1.5,2),2,
#                                ifelse(gistic.del[,j]==0,0,gistic.del[,j])))
# }
gistic.del=matrix[,c("Del_10p15.3","Del_10q24.32","Del_10q26.3","Del_11q22.1","Del_12p13.2","Del_12q24.31","Del_13q12.11","Del_13q14.2","Del_14q24.3","Del_14q32.33",
                     "Del_16p11.2","Del_16p13.3","Del_16q12.1","Del_16q22.2","Del_16q24.3","Del_17p13.1","Del_1p12","Del_1p22.1","Del_1p32.3","Del_20p13",   
                     "Del_20q13.12","Del_21p11.2","Del_22q11.22","Del_22q13.32","Del_2p11.2","Del_2q31.1","Del_2q37.3","Del_6q26","Del_7p11.2","Del_7p22.2",  
                     "Del_8p21.3","Del_8p23.1","Del_8q24.21","Del_9p21.3","database" )]
for(j in c(1:(ncol(gistic.del)-1))){
  gistic.del[,j]=ifelse(gistic.del[,j] %in% c(0.5,1),1,
                        ifelse(gistic.del[,j] %in% c(1.5,2),2,
                               ifelse(gistic.del[,j]==0,0,gistic.del[,j])))
}



amp=matrix[grep("Amp",colnames(matrix))]
del=matrix[grep("Del",colnames(matrix))]
same=intersect(gsub("Amp_","",colnames(amp)),gsub("Del_","",colnames(del)))
matrix$sample=rownames(matrix)
matrix[,c("sample",paste("Amp_",same,sep=""),paste("Del_",same,sep=""))]

gistic.other=merge(gistic.amp,gistic.del,by="row.names")
rownames(gistic.other)=gistic.other[,1]
gistic.other=gistic.other[,-1]

write.table(gistic.other,
            paste(WF,"/output/gistic.other_All_peaks.txt",sep=""),
            col.names = T,
            row.names = T,
            quote = F,
            sep = "\t")
#


###############################################################################################################################################################
######
# . . . . .                                    PredictionModel - Matrix 1q
######
################################################################################################################################################################
#####################
######## --- .                 
#####################

matrix=read.delim(paste(WF,"/output/CNV_matrix_focal_large_5Mb.txt",sep=""),
                  stringsAsFactors = F)

x1q=matrix[,c("Amp_1q21.1","Amp_1q21.3","Amp_1q44","database")]
for(j in c(1:(ncol(x1q)-1))){
  x1q[,j]=ifelse(x1q[,j] %in% c(0.5,1),1,
                 ifelse(x1q[,j] %in% c(1.5,2),2,
                        ifelse(x1q[,j]==0,0,x1q[,j])))
}

write.table(x1q,
            paste(WF,"/output/x1q_amp_gain.txt",sep=""),
            col.names = T,
            row.names = T,
            quote = F,
            sep = "\t")
#

###############################################################################################################################################################
######
# . . . . .                                    PredictionModel - Matrix Translocations
######
################################################################################################################################################################
#####################
######## --- .                 
#####################

matrix=read.delim(paste(WF,"/input/commpass_info_test.txt",sep=""),
                  stringsAsFactors = F)
translocation=matrix[,c(1,163:165)]

for(j in c(1:(ncol(translocation)-1))){
  translocation[,j]=ifelse(translocation[,j] %in% c(1),2,translocation[,j])
}

rownames(translocation)=translocation[,1]
translocation=translocation[,-1]

write.table(translocation,
            paste(WF,"/output/translocations.txt",sep=""),
            col.names = T,
            row.names = T,
            quote = F,
            sep = "\t")
#




