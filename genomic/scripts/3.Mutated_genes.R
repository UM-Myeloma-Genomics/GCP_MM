##################################################################################
##################################################################################
##############################################################################1###
#                                                                                #
#        .     .     .                                               .     .     .
####     .     .     .   PredictionModel - Mutated Genes (Oncogenes, TSGs and Unknown)   .     .     .     ####
#        .     .     .                                               .     .     .
#                                                                                #
#B###########################################################################Z####
##################################################################################
##################################################################################



# * * * * * * * * * *         IMPORTANT         * * * * * * * * * * 
#
#                          Change the folder location (WF: working folder)
#
# * * * * * * * * * *         IMPORTANT         * * * * * * * * * * 


WF="~/Google Drive/My Drive/Bachisio/Documents/Project/MM/Prediction_Model/GCP_MM/genomic/" # /where/is/the/folder/genomic
UCSC.HG19=paste(WF,"/Reference/UCSC_hg19_ref_genes.dms",sep="")



##############################################################################################################################################################
######
# . . . . .                                    PredictionModel - ONCOGENE
######
################################################################################################################################################################


library(VariantAnnotation)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Hsapiens.UCSC.hg19)
library(Biostrings)
library(annotate)
library(org.Hs.eg.db)
library(dplyr)
library(stringr)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

#####################
######## --- .       Check old matrix          
#####################

matrix=read.delim(paste(WF,"/input/commpass_info_test.txt",sep=""),
                        stringsAsFactors = F)

#####################
########
######## --- .       Selection drivers            ######
########
#####################
drivers=read.delim(paste(WF,"/Reference/Driver_summary.txt",sep=""),
                   stringsAsFactors = F)
drivers$n.CallBy.new=drivers$DNDS+drivers$MUTSIG+drivers$ONCODRIVER+drivers$FISHHOOK

drivers1=drivers[drivers$n.CallBy.new>1 | drivers$MM.driver=="MM",]

colnames(drivers1)[1]="Gene.Symbol"
unique(drivers1[order(drivers1$Gene.Symbol),"Gene.Symbol"])

drivers2=drivers1[,c("Gene.Symbol","Role","Role_in_Cancer","MM.driver")]

#####################
########
######## --- .       Oncogenes             ######        
########
#####################
oncogene=drivers2[drivers2$Role=="Oncogene",]
oncogene=oncogene[complete.cases(oncogene$Gene.Symbol),]

#########################################################################################################
########
######## --- .                                   CoMMpass Cohort                      ######
########
#########################################################################################################
# . . . - All mutation file by DNDS analysis
mut.dnds=read.delim(paste(WF,"/input/snv/CoMMpass_all_coding_mutations_annot_test.txt",sep=""),
                    stringsAsFactors = F)
mut.dnds=mut.dnds[mut.dnds$impact != "Synonymous",]
mut.dnds.gene=mut.dnds[mut.dnds$gene %in% oncogene$Gene.Symbol,]

mat=as.data.frame.matrix(table(mut.dnds.gene$sampleID,mut.dnds.gene$gene))
mat=as.matrix(mat)
mat=ifelse(mat>1,1,mat)
mat=as.data.frame(mat)
samples.NA=matrix[matrix$Genomic=="No","sample"]
# - Using the old matrix - NA matrix
matrix.NA=matrix(NA,ncol=ncol(mat),nrow=length(samples.NA))
colnames(matrix.NA)=colnames(mat)
rownames(matrix.NA)=samples.NA
# Matrix with the samples without mutations
samples.0=setdiff(matrix[matrix$Genomic=="Yes","sample"],c(rownames(mat),rownames(matrix.NA)))
matrix.0=matrix(0,ncol=ncol(mat),nrow=length(samples.0))
colnames(matrix.0)=colnames(mat)
rownames(matrix.0)=samples.0

matrix.oncogene=rbind(mat,matrix.0,matrix.NA)
dim(matrix.oncogene)


missing=setdiff(unique(oncogene$Gene.Symbol),colnames(matrix.oncogene))
for(i in 1:length(missing)){
  matrix.oncogene[,ncol(matrix.oncogene)+1]=0
  matrix.oncogene[,ncol(matrix.oncogene)]=ifelse(rownames(matrix.oncogene) %in% samples.NA,NA,matrix.oncogene[,ncol(matrix.oncogene)])
  colnames(matrix.oncogene)[ncol(matrix.oncogene)]=missing[i]
}
matrix.oncogene=matrix.oncogene[,order(colnames(matrix.oncogene))]

write.table(matrix.oncogene,
            paste(WF,"/output/oncogene_commpass.txt",sep=""),
            col.names = T,
            row.names = T,
            quote = F,
            sep="\t")





# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# . . . -> REMEMBER: - - - For ATRX,KDM5C,KDM6A,RPL10 we consider only the mutations because they are in the X chromosome.
#
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#####################
########
######## --- .       Tumor Suppressor Genes - TSG         . ---########        
########
#####################
tsg=drivers2[drivers2$Role=="TSG",]
tsg=tsg[complete.cases(tsg$Gene.Symbol),]

#########################################################################################################
########
######## --- .                                   CoMMpass Cohort                      ######
########
#########################################################################################################
# . . . - All mutation file by DNDS analysis
mut.dnds=read.delim(paste(WF,"/input/snv/CoMMpass_all_coding_mutations_annot_test.txt",sep=""),
                    stringsAsFactors = F)
mut.dnds=mut.dnds[mut.dnds$impact != "Synonymous",]
mut.dnds.gene=mut.dnds[mut.dnds$gene %in% tsg$Gene.Symbol,]

mat=as.data.frame.matrix(table(mut.dnds.gene$sampleID,mut.dnds.gene$gene))
mat=as.matrix(mat)
mat=ifelse(mat>1,1,mat)
mat=as.data.frame(mat)
samples.NA=matrix[matrix$Genomic=="No","sample"]
# -  NA matrix
matrix.NA=matrix(NA,ncol=ncol(mat),nrow=length(samples.NA))
colnames(matrix.NA)=colnames(mat)
rownames(matrix.NA)=samples.NA
# - Matrix with the samples without mutations
samples.0=setdiff(matrix[matrix$Genomic=="Yes","sample"],c(rownames(mat),rownames(matrix.NA)))
matrix.0=matrix(0,ncol=ncol(mat),nrow=length(samples.0))
colnames(matrix.0)=colnames(mat)
rownames(matrix.0)=samples.0

matrix.tsg=rbind(mat,matrix.0,matrix.NA)
dim(matrix.tsg)

missing=setdiff(unique(tsg$Gene.Symbol),colnames(matrix.tsg))
for(i in 1:length(missing)){
  matrix.tsg[,ncol(matrix.tsg)+1]=0
  matrix.tsg[,ncol(matrix.tsg)]=ifelse(rownames(matrix.tsg) %in% samples.NA,NA,matrix.tsg[,ncol(matrix.tsg)])
  colnames(matrix.tsg)[ncol(matrix.tsg)]=missing[i]
}
matrix.tsg=matrix.tsg[,order(colnames(matrix.tsg))]

write.table(matrix.tsg,
            paste(WF,"/output/tmp/tsg_mutated_genes_commpass_test.txt",sep=""),
            col.names = T,
            row.names = T,
            quote = F,
            sep="\t")



#####
# . . . - Upload GISTIC file to performe a comparison at the end of the analysis
gistic=read.delim(paste(WF,"/output/gistic.other_All_peaks.txt",sep=""),
                  stringsAsFactors = F)
gistic=gistic[,grep("Del_",colnames(gistic))]
rownames(gistic)=gsub("_1_BM","",rownames(gistic))
gistic.study=gistic[rownames(gistic) %in% matrix[matrix$study=="MMRF","sample"],]

gistic.study.1=gistic.study[,c("Del_1p12","Del_17p13.1","Del_1p32.3",
                               "Del_1p22.1"),]
colnames(gistic.study.1)=c("Del_FAM46C","Del_TP53","Del_CDKN2C",
                           "Del_RPL5")
#####

#######
# . . . - Upload study CNV file to check the gene region CNV values 
#######
cnv=read.delim(paste(WF,"/input/cnv/commpass_cnv_test.txt",sep=""),   
               stringsAsFactors = F )
cnv$sample=gsub("_1_BM","",cnv$sample)

tsg$Gene.Symbol

ref.genes=read.delim(UCSC.HG19,
                     stringsAsFactors = F)
ts.genes=ref.genes[ref.genes$UCSC.hg19.kgXref.geneSymbol %in% tsg$Gene.Symbol,]
ts.genes=ts.genes[,c("UCSC.hg19.knownGene.chrom","UCSC.hg19.knownGene.txStart","UCSC.hg19.knownGene.txEnd",
                     "UCSC.hg19.kgXref.geneSymbol")]
ts.genes=ts.genes %>%
  distinct(UCSC.hg19.kgXref.geneSymbol,.keep_all = T)
cnv.1=cnv[cnv$sample %in% unique(matrix$sample),]
length(unique(cnv.1$sample))
head(cnv.1)
head(ts.genes)
colnames(ts.genes)=c("chrom","start","end","Symbol")
ts.genes$chrom=gsub("chr","",ts.genes$chrom)


gr.cnv=with(cnv.1, GRanges(chr, IRanges(start=start-50000, end=end+50000)))
values(gr.cnv)=cnv.1[,c(1,5:ncol(cnv.1))]

gr.tsg=with(ts.genes,GRanges(chrom,IRanges(start=start,end=end)))
values(gr.tsg)=ts.genes[,c(4:ncol(ts.genes))]

ranges <- merge(as.data.frame(gr.cnv),
                as.data.frame(gr.tsg),
                by="seqnames",
                suffixes=c("A","B"))

ranges_1=ranges[with(ranges, startB <= startA & endB <= endA & endB > startA),]
ranges_2=ranges[with(ranges, startB >= startA & endB >= endA & startB < endA),]
ranges_3=ranges[with(ranges, startB >= startA & endB <= endA),]
ranges_4=ranges[with(ranges, startB < startA & endB > endA),]
ranges.all<- rbind.data.frame(ranges_1, ranges_2, ranges_3, ranges_4)

ranges.all=unique(ranges.all)
ranges.all.1=ranges.all %>%
  arrange(sample,seqnames,startA)
head(ranges.all.1)
ranges.all.2=ranges.all.1[,c("sample","seqnames","startA","endA","tot","min","X")]
unique(ranges.all.2$tot)
table(ranges.all.2$tot)
table(round(ranges.all.2$tot,0))

################################################################################
# TOTAL CNV -> 1.5 will be loss, 2.5 will be consider diploid, 3.5 a triploid and 4.5 tetraploid
# MINOR CNV -> 0.5 will be homozygous deletion (0)

ranges.all.2$tot=ifelse(ranges.all.2$tot %in% c(2.5),2,
                          ifelse(ranges.all.2$tot %in% c(1.5),1,
                                 ifelse(ranges.all.2$tot==3.5,3,
                                        ifelse(ranges.all.2$tot==4.5,4,ranges.all.2$tot))))


table(ranges.all.2$min)
ranges.all.2$min=ifelse(ranges.all.2$min==0.5,0,ranges.all.2$min)
head(ranges.all.2)
################################################################################


ranges.all.2$major=ranges.all.2$tot-ranges.all.2$min
ranges.all.2$code=ifelse(ranges.all.2$major+ranges.all.2$min>=2 & ranges.all.2$min!=0,0,
                         ifelse(ranges.all.2$major+ranges.all.2$min==1 | (ranges.all.2$min==0 & ranges.all.2$major>=2),1,
                                ifelse(ranges.all.2$tot==0,2,NA)))

ranges.all.2[is.na(ranges.all.2$code),]

ranges.all.2$code=ifelse(is.na(ranges.all.2$code) & ranges.all.2$tot==1,1,
                         ifelse(is.na(ranges.all.2$code) & ranges.all.2$tot==0,2,
                                ifelse(is.na(ranges.all.2$code) & ranges.all.2$tot==2,0,
                                       ifelse(is.na(ranges.all.2$code) & ranges.all.2$tot>=3,0,ranges.all.2$code)
                                )
                         )
)
ranges.all.2[is.na(ranges.all.2$min),]

matrix.cnv=data.table::dcast(unique(ranges.all.2[,c("X","sample","code")]),
                             formula=sample~X,
                             toString,
                             value.var = "code")
data.table::setDF(matrix.cnv)
matrix.cnv=as.data.frame(matrix.cnv)
rownames(matrix.cnv)=matrix.cnv$sample
matrix.cnv=matrix.cnv[,-1]
for(j in 1:ncol(matrix.cnv)){
  matrix.cnv[,j]=ifelse(matrix.cnv[j]=="",0,matrix.cnv[,j])
}
for(i in 1:length(colnames(matrix.cnv))){
  for(z in 1:length(rownames(matrix.cnv))){
    cel=matrix.cnv[z,i]
    c1=str_split_fixed(cel,",",Inf)
    c1=as.matrix(c1)
    c1=gsub(" ","",c1)
    c1=as.data.frame(c1)
    for(zz in 1:length(colnames(c1))){
      c1[,zz]=as.numeric(c1[,zz])
    }
    c1$max=max(c1)
    matrix.cnv[z,i]=c1$max
  }
}


samples.NA=matrix[matrix$Genomic=="No","sample"]
# -  NA matrix
matrix.NA=matrix(NA,ncol=ncol(mat),nrow=length(samples.NA))
colnames(matrix.NA)=colnames(mat)
rownames(matrix.NA)=samples.NA
# - Matrix with the samples without mutations
samples.0=setdiff(cnv$sample,c(rownames(matrix.cnv),rownames(matrix.NA)))
matrix.0=matrix(0,ncol=ncol(matrix.cnv),nrow=length(samples.0))
colnames(matrix.0)=colnames(matrix.cnv)
rownames(matrix.0)=samples.0

matrix.tsg=rbind(matrix.cnv,matrix.0,matrix.NA)
dim(matrix.tsg)

colnames(matrix.tsg)=paste("Del_",colnames(matrix.tsg),sep="")


write.table(matrix.tsg,
            paste(WF,"/output/tmp/tsg_deleted_genes_commpass_test.txt",sep=""),
            col.names = T,
            row.names = T,
            quote = F,
            sep="\t")

##################################################################################################################################
##################################################################################################################################

# . . . Merge two matrices
# Mutated gene matrix
mut=read.delim(paste(WF,"/output/tmp/tsg_mutated_genes_commpass_test.txt",sep=""),
               stringsAsFactors = F)
# missing1=setdiff(tsg$Gene.Symbol,colnames(mut))
# for(i in 1:length(missing1)){
#   mut[,ncol(mut)+1]=0
#   colnames(mut)[ncol(mut)]=missing1[i]
# }
mut=mut[,order(colnames(mut))]
 
colnames(mut)=paste("Mut_",colnames(mut),sep="")


# Deleted gene matrix
cnv=read.delim(paste(WF,"/output/tmp/tsg_deleted_genes_commpass_test.txt",sep=""),
               stringsAsFactors = F)

missing2=setdiff(tsg$Gene.Symbol,gsub("Del_","",colnames(cnv)))
# 0
# for(i in 1:length(missing2)){
#   cnv[,ncol(cnv)+1]=0
#   colnames(cnv)[ncol(cnv)]=paste("Del_",missing2[i],sep="")
# }
cnv=cnv[,order(colnames(cnv))]


study.matrix=merge(cnv,mut,by="row.names")

study.matrix.1=study.matrix
rownames(study.matrix.1)=study.matrix.1$Row.names
study.matrix.1=study.matrix.1[,-1]

head(study.matrix.1)

for(i in 1:ncol(study.matrix.1)){
  study.matrix.1[,i]=ifelse(is.na(study.matrix.1[,i]),999,study.matrix.1[,i])
}

for(j in 1:36){
  study.matrix.1[,72+j]=rowSums(study.matrix.1[,c(j,36+j)],na.rm = T)
  colnames(study.matrix.1)[ncol(study.matrix.1)]=gsub("Del_","",colnames(study.matrix.1)[j])
}


for(z in 1:72){
  study.matrix.1[,z]=ifelse(study.matrix.1[,z]==999,NA,study.matrix.1[,z])
}
for(z in 73:ncol(study.matrix.1)){
  study.matrix.1[,z]=ifelse(study.matrix.1[,z] %in% c(3,1001),2,
                            ifelse(study.matrix.1[,z]==1000,1,
                                   ifelse(study.matrix.1[,z]==999,0,
                                          ifelse(study.matrix.1[,z]==1998,NA,study.matrix.1[,z]))))
}


summary(study.matrix.1[,73:108]) # only new columns
head(study.matrix.1)
# ATRX,KDM5C,KDM6A,RPL10
colnames(study.matrix.1)
study.matrix.2=study.matrix.1[,c(73:108)]
colnames(study.matrix.2)=gsub("Mut_","",colnames(study.matrix.2))
study.matrix.2=study.matrix.2[,order(colnames(study.matrix.2))]

head(study.matrix.2)
dim(study.matrix.2)
write.table(study.matrix.2,
            paste(WF,"/output/tsg_commpass_test.txt",sep=""),
            col.names = T,
            row.names = T,
            quote = F,
            sep="\t")



#####################
########
######## --- .       Unknown         . ---########        
########
#####################
ukn=drivers2[drivers2$Role=="Unknown",]
ukn1=drivers2[is.na(drivers2$Role),]

ukn=rbind(ukn,ukn1)
ukn=ukn[complete.cases(ukn$Gene.Symbol),]


#########################################################################################################
########
######## --- .                                   CoMMpass Cohort                                   ######
########
#########################################################################################################
# . . . - All mutation file by DNDS analysis
mut.dnds=read.delim(paste(WF,"/input/snv/CoMMpass_all_coding_mutations_annot_test.txt",sep=""),
                    stringsAsFactors = F)
mut.dnds=mut.dnds[mut.dnds$impact != "Synonymous",]
mut.dnds.gene=mut.dnds[mut.dnds$gene %in% ukn$Gene.Symbol,]

mat=as.data.frame.matrix(table(mut.dnds.gene$sampleID,mut.dnds.gene$gene))
mat=as.matrix(mat)
mat=ifelse(mat>1,1,mat)
mat=as.data.frame(mat)
samples.NA=matrix[matrix$Genomic=="No","sample"]
# -  NA matrix
matrix.NA=matrix(NA,ncol=ncol(mat),nrow=length(samples.NA))
colnames(matrix.NA)=colnames(mat)
rownames(matrix.NA)=samples.NA
# - Matrix with the samples without mutations
samples.0=setdiff(matrix[matrix$Genomic=="Yes","sample"],c(rownames(mat),rownames(matrix.NA)))
matrix.0=matrix(0,ncol=ncol(mat),nrow=length(samples.0))
colnames(matrix.0)=colnames(mat)
rownames(matrix.0)=samples.0

matrix.ukn=rbind(mat,matrix.0,matrix.NA)
dim(matrix.ukn)

missing=setdiff(unique(ukn$Gene.Symbol),colnames(matrix.ukn))
for(i in 1:length(missing)){
  matrix.ukn[,ncol(matrix.ukn)+1]=0
  matrix.ukn[,ncol(matrix.ukn)]=ifelse(rownames(matrix.ukn) %in% samples.NA,NA,matrix.ukn[,ncol(matrix.ukn)])
  colnames(matrix.ukn)[ncol(matrix.ukn)]=missing[i]
}
matrix.ukn=matrix.ukn[,order(colnames(matrix.ukn))]


write.table(matrix.ukn,
            paste(WF,"/output/unk_mutated_genes_commpass_test.txt",sep=""),
            col.names = T,
            row.names = T,
            quote = F,
            sep="\t")



#########################################################################################################
########
######## --- .                            Chromosome X Genes                               . --- ########   
########
#########################################################################################################

# "ATRX","KDM5C","KDM6A","RPL10"

commpass=read.delim(paste(WF,"/output/tmp/tsg_mutated_genes_commpass_test.txt",sep=""),
                    stringsAsFactors = F)
commpass=commpass[,colnames(commpass) %in% c("ATRX","KDM5C","KDM6A","RPL10")]

write.table(commpass,
            paste(WF,"/output/chrx_tsg_commpass_test.txt",sep=""),
            col.names = T,
            row.names = T,
            quote = F,
            sep = "\t")
