##################################################################################
##################################################################################
##############################################################################1###
#                                                                                #
#        .     .     .                                               .     .     .
####     .     .     .   GISTIC - Identification of abnormal peaks   .     .     .     ####
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
WF="~/Google Drive/My Drive/Bachisio/Documents/Project/MM/Prediction_Model/genomic/" # /where/is/the/folder/genomic

##################################################################################
##################################################################################
##################################################################################
##################################################################################


# Variables
CYTOBAND.FILE=paste(WF,"/Reference/cytoBand.txt",sep="")
ARMS=paste(WF,"/Reference/arm_hg19.txt",sep="")
CHROM.LENGHT=paste(WF,"/Reference/hg19.chrom_sizes.txt",sep="")

GISTIC.FILE=paste(WF,"/Reference/Genomic_Coordinates_GISTICpeaks.txt",sep="")
GISTIC.ARMS=paste(WF,"/Reference/arm_level.txt",sep="")

INPUT.PATH=paste(WF,"/input/cnv/All_coMMpass/",sep="")
OUT.DIR=paste(WF,"/output/",sep="")
PERCENTAGE=0.60 #(percentage of chromosome with gain/amp [for hyperdiploid chromosomes])

##################################################################################
##################################################################################
##################################################################################
##################################################################################
library(ggplot2)
library(gtools)
library(pheatmap)
library(RColorBrewer)
library("GenomicRanges")
library(dplyr)
library(plyr)
library(ComplexHeatmap)
library(stringr)

################################################################################################################################################################
######
# . . . . .                                              Matrix with distinction focal and large
######
################################################################################################################################################################
#####################
######## --- .                 Creation matrix with abnormal GISTIC peaks & CNV information
#####################
cytoband=read.delim(CYTOBAND.FILE,stringsAsFactors = F,header=F)
colnames(cytoband)=c("X.chrom","chromStart","chromEnd","name","type")

cytoband$cytoband=paste(gsub("chr","",cytoband$X.chrom),cytoband$name,sep="")

gistic.peaks=read.delim(GISTIC.FILE,
                        stringsAsFactors = F)
gistic.peaks.1=merge(gistic.peaks,cytoband,by="cytoband",all.x=T)

table(gistic.peaks.1$type)
gistic.peaks.1[gistic.peaks.1$type=="acen",]


new.gistic=gistic.peaks
new.gistic$gistic.seg.length=new.gistic$end-new.gistic$start

setwd(INPUT.PATH)
list.txt=list.files(patter="txt",full.names=T)
list.txt<-str_split_fixed(list.txt,"/",2)[,2]
name.dataset=str_split_fixed(list.txt,"_",3)[,1]

all.cohort=NULL
for(j in 1:length(list.txt)){
  cnv=read.delim(list.txt[j],stringsAsFactors = F)
  cnv$type=rep(name.dataset[j],length(cnv$sample))
  cnv$seg.length=cnv$end-cnv$start
  cnv$seg.type=ifelse(cnv$seg.length<=5000000,"focal","large")
  new.gistic$chr=gsub("chr","",new.gistic$chr)
  
  gr.cnv<- with(cnv, GRanges(chr, IRanges(start=start, end=end)))
  values(gr.cnv) <- cnv[,c(1,5:ncol(cnv))]
  
  gr.gistic = with(new.gistic, GRanges(chr, IRanges(start=start, end=end)))
  values(gr.gistic) <- new.gistic[,c(4:ncol(new.gistic))]
  
  ranges <- merge(as.data.frame(gr.cnv),
                  as.data.frame(gr.gistic),
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
  
  amp=ranges.all.1[ranges.all.1$lesion=="Amp",]
  # Gain Large 0.5
  # Gain Focal 1
  # Amp Large 1.5
  # Amp Focal 2
  # Normal => 0
  amp$anomaly=ifelse(amp$tot>=4 & amp$seg.type=="focal",2,
                     ifelse(amp$tot>=4 & amp$seg.type=="large",1.5,
                            ifelse(amp$tot>=3 & amp$tot<=3.5 & amp$seg.type=="focal",1,
                                   ifelse(amp$tot>=3 & amp$tot<=3.5 & amp$seg.type=="large",0.5,0))))
  table(amp$anomaly)
  table(amp$anomaly,amp$seg.type)
  amp.mtx=data.table::dcast(unique(amp[,c("sample","cytoband","anomaly")]),
                            formula=cytoband~sample,
                            toString,
                            value.var = "anomaly")
  data.table::setDF(amp.mtx)
  rownames(amp.mtx)=amp.mtx$cytoband
  amp.mtx=amp.mtx[,-1]
  amp.mtx=as.matrix(amp.mtx)
  amp.mtx=ifelse(amp.mtx=="",NA,amp.mtx)
  amp.mtx=as.data.frame(amp.mtx)
  for(i in 1:length(colnames(amp.mtx))){
    for(z in 1:length(rownames(amp.mtx))){
      cel=amp.mtx[z,i]
      c1=str_split_fixed(cel,",",Inf)
      c1=as.matrix(c1)
      c1=gsub(" ","",c1)
      c1=as.data.frame(c1)
      for(zz in 1:length(colnames(c1))){
        c1[,zz]=as.numeric(c1[,zz])
      }
      c1$max=ifelse(is.na(c1),
                    NA,
                    max(c1,na.rm = T))
      amp.mtx[z,i]=c1$max
    }
  }
  amp.mtx[1:5,1:5]
  missing.samples.amp=setdiff(unique(cnv$sample),colnames(amp.mtx))
  if(length(missing.samples.amp)>0){
    mis.amp.mtx=matrix(NA,nrow=length(rownames(amp.mtx)),ncol=length(missing.samples.amp))
    mis.amp.mtx=ifelse(is.na(mis.amp.mtx),0,mis.amp.mtx)
    amp.mtx1=cbind(amp.mtx,mis.amp.mtx)
  }else{
    amp.mtx1=amp.mtx
  }
  amp.mtx1[1:5,1:5]
  rownames(amp.mtx1)=paste("Amp_",rownames(amp.mtx1),sep="")
  
  
  del=ranges.all.1[ranges.all.1$lesion=="Del",]
  # Biallelic Focal => 2
  # Biallelic Large => 1.5
  # Monoallelic Focal => 1
  # Monoallelic Large => 0.5
  # Normal => 0
  del$pr=ifelse(del$tot>=2 & is.na(del$min),0,"qqq")
  del$anomaly=ifelse(del$tot==0 & del$seg.type=="focal",2,
                     ifelse(del$tot==0 & del$seg.type=="large",1.5,
                            ifelse(del$tot>=2 & del$min==0 & del$seg.type=="focal",1,
                                   ifelse(del$tot<=1 & del$tot>=0.5 & del$seg.type=="focal",1,
                                          ifelse(del$tot>=2 & del$min==0 & del$seg.type=="large",0.5,
                                                 ifelse(del$tot<=1 & del$tot>=0.5 & del$seg.type=="large",0.5,0))))))
  
  table(del$anomaly,useNA="ifany")
  table(del$pr,useNA="ifany")
  del$anomaly=ifelse(is.na(del$anomaly),del$pr,del$anomaly)
  table(del$anomaly,useNA="ifany")
  
  table(del$anomaly,del$seg.type)
  del.mtx=data.table::dcast(unique(del[,c("sample","cytoband","anomaly")]),
                            formula=cytoband~sample,
                            toString,
                            value.var = "anomaly")
  data.table::setDF(del.mtx)
  rownames(del.mtx)=del.mtx$cytoband
  del.mtx=del.mtx[,-1]
  del.mtx=as.matrix(del.mtx)
  del.mtx=ifelse(del.mtx=="",NA,del.mtx)
  del.mtx=as.data.frame(del.mtx)
  for(i in 1:length(colnames(del.mtx))){
    for(z in 1:length(rownames(del.mtx))){
      cel=del.mtx[z,i]
      c1=str_split_fixed(cel,",",Inf)
      c1=as.matrix(c1)
      c1=gsub(" ","",c1)
      c1=as.data.frame(c1)
      for(zz in 1:length(colnames(c1))){
        c1[,zz]=as.numeric(c1[,zz])
      }
      c1$max=ifelse(is.na(c1),
                    NA,
                    max(c1,na.rm = T))
      del.mtx[z,i]=c1$max
    }
  }
  del.mtx[1:5,1:5]
  missing.samples.del=setdiff(unique(cnv$sample),colnames(del.mtx))
  if(length(missing.samples.del)>0){
    mis.del.mtx=matrix(NA,nrow=length(rownames(del.mtx)),ncol=length(missing.samples.del))
    mis.del.mtx=ifelse(is.na(mis.del.mtx),0,mis.del.mtx)
    del.mtx1=cbind(del.mtx,mis.del.mtx)
  }else{
    del.mtx1=del.mtx
  }
  del.mtx1[1:5,1:5]
  rownames(del.mtx1)=paste("Del_",rownames(del.mtx1),sep="")
  
  identical(colnames(amp.mtx1),colnames(del.mtx1))
  all=rbind(amp.mtx1,del.mtx1)
  t.all=t(all)
  t.all=as.data.frame(t.all)
  for(jj in 1:length(colnames(t.all))){
    t.all[,jj]=as.numeric(as.character(t.all[,jj]))
  }
  
  mis.peak=setdiff(new.gistic$code,colnames(t.all))
  if(length(mis.peak)!=0){
    for(jjj in length(mis.peak)){
      t.all[,paste(mis.peak[jjj],"",sep="")]=0
    }
  }else{
    t.all=t.all
  }
  t.all$database=rep(name.dataset[j],length(rownames(t.all)))
  all.cohort=rbind(all.cohort,t.all)
}

write.table(all.cohort,
            paste(OUT.DIR,"/CNV_matrix_focal_large_5Mb.txt",sep=""),
            col.names = T,
            row.names = T,
            quote = F,
            sep = "\t")

################################################################################################################################################################
##########################################################################################################################################################
##################################################################################
##################################################################################
#
##################################################################################1
#     .     .     .                                                   .     .     .
#     .     .     .   PredictionModel - Identification hyperdiploid   .     .     .
#     .     .     .                                                   .     .     .
#B###############################################################################Z#

library(ggplot2)
library(gtools)
library(pheatmap)
library(RColorBrewer)
library("GenomicRanges")
library(dplyr)
library(plyr)
library(ComplexHeatmap)
library(stringr)
library(stringi)
################################################################################################################################################################
######
# . . . . .                                              ARMS
######
################################################################################################################################################################



#####################
######## --- .                 Creation matrix with abnormal GISTIC arms & CNV information    
#####################
#####################
######## --- .                 CoMMpass, Bolli, MRCxi, MSK and UAMS   
#####################
cytoband=read.delim(CYTOBAND.FILE,stringsAsFactors = F,header=F)
colnames(cytoband)=c("X.chrom","chromStart","chromEnd","name","type")
cytoband$cytoband=paste(gsub("chr","",cytoband$X.chrom),cytoband$name,sep="")
gistic.arms=read.delim(GISTIC.ARMS,
                       stringsAsFactors = F)

new.gistic=gistic.arms

setwd(INPUT.PATH)
list.txt=list.files(patter="txt",full.names=T)
list.txt<-str_split_fixed(list.txt,"/",2)[,2]
name.dataset=str_split_fixed(list.txt,"_",3)[,1]

#####################
######## --- .                 Segment collapsing with the same tot and min (based on code variable)    
new=NULL
for(j in 1:length(list.txt)){
  cnv=read.delim(list.txt[j],stringsAsFactors = F)
  samples=unique(cnv$sample)
  for (i in 1:length(samples))
  {
    cnv.s=cnv[cnv$sample==samples[i],]
    chroms=unique(cnv.s$chr)
    cnv.s$tot=round(cnv.s$tot,0)
    cnv.s$min=round(cnv.s$min,0)
    for(z in 1:length(chroms)){
      chr<-cnv.s[cnv.s$chr==chroms[z],]
      chr$code=paste(chr$tot,chr$min,sep="_")
      
      a<-rle(chr$code)$lengths
      chr$first<-unlist(lapply(a, function(x) seq(1,x)))
      chr$last<-unlist(lapply(a, function (x) seq(x,1,-1)))
      head(chr)
      
      chr1<-chr[chr$first==1 | chr$last==1,]
      head(chr1)
      
      a1<-rle(chr1$code)$lengths
      chr1$first2<-unlist(lapply(a1, function(x) seq(1,x)))
      chr1$last2<-unlist(lapply(a1, function (x) seq(x,1,-1)))
      head(chr1)
      dim(chr1)
      
      ##_____ 
      chr1_start<-subset(chr1,chr1$first2==1)
      head(chr1_start)
      chr1_end<-subset(chr1,chr1$last2==1)
      head(chr1_end)
      colnames(chr1_end)<-c("sample","chr1","start1","end1","tot1","min1","code","first1","last1","first21","last21")
      head(chr1_end)
      dim(chr1_start)
      dim(chr1_end)
      ##_____ 
      df1<-data.frame(chr1_start$sample,chr1_start$chr,chr1_start$start,chr1_end$end1,chr1_end$tot1,chr1_end$min1,
                      chr1_end$code,chr1_end$first1)
      df1
      colnames(df1)<-c("sample","chr","start","end","tot","min","code","Num_Segments")
      df1
      head(df1)
      dim(df1)
      df1$type=rep(name.dataset[j],length(df1$sample))
      new=rbind(new,df1)
    }
  }
}
cnv.n=new
cnv.gain=cnv.n[cnv.n$tot>2,]
cnv.gain$seg.length=cnv.gain$end-cnv.gain$start

cnv.gain1=cnv.gain

#####################
######## --- .                 Add arms information
arms=read.delim(ARMS,stringsAsFactors = F)
arms$arm.length=arms$End_Position-arms$Start_Position

gr.cnv<- with(cnv.gain, GRanges(chr, IRanges(start=start, end=end)))
values(gr.cnv) <- cnv.gain[,c(1,5:ncol(cnv.gain))]

gr.arms = with(arms, GRanges(Chromosome, IRanges(start=Start_Position, end=End_Position)))
values(gr.arms) <- arms[,c(1,5)]

ranges <- merge(as.data.frame(gr.cnv),
                as.data.frame(gr.arms),
                by="seqnames",
                suffixes=c("A","B"))


ranges_1=ranges[with(ranges, startB <= startA & endB <= endA & endB > startA),]
ranges_2=ranges[with(ranges, startB >= startA & endB >= endA & startB < endA),]
ranges_3=ranges[with(ranges, startB >= startA & endB <= endA),]
ranges_4=ranges[with(ranges, startB < startA & endB > endA),]
ranges.all<- rbind.data.frame(ranges_1, ranges_2, ranges_3, ranges_4)

ranges.all=unique(ranges.all)
ranges.all.1=ranges.all %>%
  arrange(type,sample,seqnames,startA)
head(ranges.all.1)

cnv.gain=ranges.all.1[,c("sample","seqnames","startA","endA","tot","min","code","Num_Segments","type","startB","endB","Arm","arm.length")]
colnames(cnv.gain)=c(c("sample","chr","start","end","tot","min","code","Num_Segments","type","start.arm","end.arm","arm","arm.length"))

#cnv.gain$p.q=substring(cnv.gain$arm,2) # keep the second character
cnv.gain$p.q=stri_sub(cnv.gain$arm,-1)  # keep the last character

cnv.gain$seg.length=ifelse(cnv.gain$end <= cnv.gain$end.arm & cnv.gain$start > cnv.gain$start.arm,
                           cnv.gain$end-cnv.gain$start,
                           ifelse(cnv.gain$end > cnv.gain$end.arm & cnv.gain$start > cnv.gain$start.arm,
                                  cnv.gain$end.arm-cnv.gain$start,
                                  ifelse(cnv.gain$start < cnv.gain$start.arm,
                                         cnv.gain$end-cnv.gain$start.arm,NA)))
#save(cnv.gain,file = "~/Box/Bachisio/Documents/Project/MM/Prediction_Model/analysis/output/Rdata/cnv.gain.RData")                               
#load("~/Box/Bachisio/Documents/Project/MM/Prediction_Model/analysis/output/Rdata/cnv.gain.RData")
#####################
######## --- .                 Chromosome identification with more than CN 2 in >60% of the chromosome (and no.hyperdiploid)
chrom.length=read.delim(CHROM.LENGHT,
                        stringsAsFactors = F)
f=c(23,155270560)
m=c(24,59373566)
chrom.length=rbind(chrom.length,f,m)

cnv.gain$chr=as.character(cnv.gain$chr)

hyperdiploid=NULL
no.hyperdiploid=NULL
samples=unique(cnv.n$sample)
for(ii in 1:length(samples)){
  new.s=cnv.gain[cnv.gain$sample==samples[ii],]
  chroms=unique(new.s$chr)
  ARMS=unique(new.s$arm)
  if(length(new.s$sample)>0){
    for(zz in 1:length(chroms)){
      chr<-new.s[new.s$chr==chroms[zz],]
      for(zzz in 1:length(ARMS)){
        arm<-new.s[new.s$arm==ARMS[zzz],]
        sam=samples[ii]
        cro=unique(arm$chr)
        CN=median(as.numeric(as.character(arm$tot)),na.rm = T)
        max.cn=max(as.numeric(as.character(arm$tot)),na.rm = T)
        min.cn=min(as.numeric(as.character(arm$tot)),na.rm = T)  
        type=unique(arm$type)
        length.seg=sum(as.numeric(as.character(arm$seg.length)))
        p.q=unique(arm$p.q)
        seg.gain.prop=as.numeric(as.character(length.seg))/as.numeric(as.character((arms[arms$Arm==ARMS[zzz],"arm.length"])))
        line=c(sam,cro,CN,max.cn,min.cn,type,length.seg,p.q,seg.gain.prop)
        hyperdiploid=rbind(hyperdiploid,line)
      }
    }
  }else{
    no.hyperdiploid=rbind(no.hyperdiploid,samples[ii])
  }
}

hyper=as.data.frame(hyperdiploid)

#save(hyper,file="~/Box/Bachisio/Documents/Project/MM/Prediction_Model/analysis/output/hyper.arms.RData")
#save(no.hyperdiploid,file="~/Box/Bachisio/Documents/Project/MM/Prediction_Model/analysis/output/no.hyper1.arms.RData")
#load("~/Box/Bachisio/Documents/Project/MM/Prediction_Model/analysis/output/hyper.arms.RData")
#load("~/Box/Bachisio/Documents/Project/MM/Prediction_Model/analysis/output/no.hyper1.arms.RData")
rownames(hyper)=c(1:length(hyper$V1))
colnames(hyper)=c("sample","chromosome","CN.median","CN.max","CN.min","type","length.seg","arm","seg.gain.prop")

hyper$seg.gain.prop=as.numeric(hyper$seg.gain.prop)
hyper=unique(hyper)

samples=unique(hyper$sample)
hyper.new=NULL
for(ii in 1:length(samples)){
  new.s=hyper[hyper$sample==samples[ii],]
  chroms=unique(new.s$chr)
  for(zz in 1:length(chroms)){
    chr<-new.s[new.s$chr==chroms[zz],]
    sam=samples[ii]
    cro=unique(chr$chr)
    CN.p=chr[chr$arm=="p","CN.median"]
    CN.p[length(CN.p)==0]=NA
    max.cn.p=chr[chr$arm=="p","CN.max"]
    max.cn.p[length(max.cn.p)==0]=NA
    min.cn.p=chr[chr$arm=="p","CN.min"]
    min.cn.p[length(min.cn.p)==0]=NA
    length.seg.p=chr[chr$arm=="p","length.seg"]
    length.seg.p[length(length.seg.p)==0]=NA
    arm.p=chr[chr$arm=="p","arm"]
    arm.p[length(arm.p)==0]=NA
    seg.gain.prop.p=chr[chr$arm=="p","seg.gain.prop"]
    seg.gain.prop.p[length(seg.gain.prop.p)==0]=NA
    
    CN.q=chr[chr$arm=="q","CN.median"]
    CN.q[length(CN.q)==0]=NA
    max.cn.q=chr[chr$arm=="q","CN.max"]
    max.cn.q[length(max.cn.q)==0]=NA
    min.cn.q=chr[chr$arm=="q","CN.min"]
    min.cn.q[length(min.cn.q)==0]=NA
    length.seg.q=chr[chr$arm=="q","length.seg"]
    length.seg.q[length(length.seg.q)==0]=NA
    arm.q=chr[chr$arm=="q","arm"]
    arm.q[length(arm.q)==0]=NA
    seg.gain.prop.q=chr[chr$arm=="q","seg.gain.prop"]
    seg.gain.prop.q[length(seg.gain.prop.q)==0]=NA
    
    type=unique(chr$type)
    gain.length=sum(as.numeric(as.character(chr$length.seg)))
    perc.gain.chr=as.numeric(as.character(gain.length))/(as.numeric(as.character(chrom.length[chrom.length$Offset==chroms[zz],"Size"])))
    line.p=c(sam,cro,CN.p,max.cn.p,min.cn.p,type,length.seg.p,arm.p,seg.gain.prop.p,gain.length,perc.gain.chr)
    line.q=c(sam,cro,CN.q,max.cn.q,min.cn.q,type,length.seg.q,arm.q,seg.gain.prop.q,gain.length,perc.gain.chr)
    hyper.new=rbind(hyper.new,line.p,line.q)
  }
}


hyper.new.1=as.data.frame(hyper.new)
rownames(hyper.new.1)=c(1:length(hyper.new.1$V1))
colnames(hyper.new.1)=c("sample","chromosome","CN.median","CN.max","CN.min","type","length.seg","arm","seg.gain.prop","gain.length","perc.gain.chr")

hyper.new.2=hyper.new.1[complete.cases(hyper.new.1$CN.median),]
hyper.new.2$arm.gain=ifelse(hyper.new.2$seg.gain.prop>PERCENTAGE,1,0)
hyper.new.2$chr.gain=ifelse(hyper.new.2$perc.gain.chr>PERCENTAGE,1,0)
head(hyper.new.2)
# #########
# # . . # # Understand
# #     # #
# ####################
# head(hyper.new.2[hyper.new.2$arm.gain==1 & hyper.new.2$chr.gain==0,])
# hyper.new.2[hyper.new.2$sample=="PD5858a" & hyper.new.2$chromosome==1,]
# # non hanno nessun gain nell'altro arm
# head(hyper.new.2[hyper.new.2$arm.gain==0 & hyper.new.2$chr.gain==1,])
# table(hyper.new.2$arm.gain)
# #0     1 
# #5317 12606 
# table(hyper.new.2$chr.gain)
# #0     1 
# #6027 11896 
# table(hyper.new.2$chr.gain,hyper.new.2$arm.gain)
# #              ARM
# #            0     1
# #  C   0  4813  1214
# #  H   
# #  R   1   504 11392

hyper=hyper.new.2[hyper.new.2$chr.gain==1,]
#####################
######## --- .                 Hyperdiploid and DoubleGenome identification  (and no.hyperdiploid)
chroms.hyper=c(3,5,7,9,11,15,19,21)
samples=unique(hyper$sample)
hyper.1=NULL
no.hyperdiploid.1=NULL
for(iii in 1:length(samples)){
  hyper.s=hyper[hyper$sample==samples[iii],]
  if(sum(hyper.s[,"arm.gain"])>30){
    sam=unique(hyper.s$sample)
    type=unique(hyper.s$type)
    whole.median.CN=median(as.numeric(as.character(hyper.s$CN.median)),na.rm = T)
    
    # . p arm
    chr3.p=hyper.s[hyper.s$chromosome==3 & hyper.s$arm=="p",]
    seg.gain.pro.chr3.p=ifelse(chr3.p$chromosome==3,chr3.p$seg.gain.prop,NA)
    seg.gain.pro.chr3.p[length(seg.gain.pro.chr3.p)==0]=NA
    chr5.p=hyper.s[hyper.s$chromosome==5 & hyper.s$arm=="p",]
    seg.gain.pro.chr5.p=ifelse(chr5.p$chromosome==5,chr5.p$seg.gain.prop,NA)
    seg.gain.pro.chr5.p[length(seg.gain.pro.chr5.p)==0]=NA
    chr7.p=hyper.s[hyper.s$chromosome==7 & hyper.s$arm=="p",]
    seg.gain.pro.chr7.p=ifelse(chr7.p$chromosome==7,chr7.p$seg.gain.prop,NA)
    seg.gain.pro.chr7.p[length(seg.gain.pro.chr7.p)==0]=NA
    chr9.p=hyper.s[hyper.s$chromosome==9 & hyper.s$arm=="p",]
    seg.gain.pro.chr9.p=ifelse(chr9.p$chromosome==9,chr9.p$seg.gain.prop,NA)
    seg.gain.pro.chr9.p[length(seg.gain.pro.chr9.p)==0]=NA
    chr11.p=hyper.s[hyper.s$chromosome==11 & hyper.s$arm=="p",]
    seg.gain.pro.chr11.p=ifelse(chr11.p$chromosome==11,chr11.p$seg.gain.prop,NA)
    seg.gain.pro.chr11.p[length(seg.gain.pro.chr11.p)==0]=NA
    chr15.p=hyper.s[hyper.s$chromosome==15 & hyper.s$arm=="p",]
    seg.gain.pro.chr15.p=ifelse(chr15.p$chromosome==15,chr15.p$seg.gain.prop,NA)
    seg.gain.pro.chr15.p[length(seg.gain.pro.chr15.p)==0]=NA
    chr19.p=hyper.s[hyper.s$chromosome==19 & hyper.s$arm=="p",]
    seg.gain.pro.chr19.p=ifelse(chr19.p$chromosome==19,chr19.p$seg.gain.prop,NA)
    seg.gain.pro.chr19.p[length(seg.gain.pro.chr19.p)==0]=NA
    chr21.p=hyper.s[hyper.s$chromosome==21 & hyper.s$arm=="p",]
    seg.gain.pro.chr21.p=ifelse(chr21.p$chromosome==21,chr21.p$seg.gain.prop,NA)
    seg.gain.pro.chr21.p[length(seg.gain.pro.chr21.p)==0]=NA
    # . q arm  
    chr3.q=hyper.s[hyper.s$chromosome==3 & hyper.s$arm=="q",]
    seg.gain.pro.chr3.q=ifelse(chr3.q$chromosome==3,chr3.q$seg.gain.prop,NA)
    seg.gain.pro.chr3.q[length(seg.gain.pro.chr3.q)==0]=NA
    chr5.q=hyper.s[hyper.s$chromosome==5 & hyper.s$arm=="q",]
    seg.gain.pro.chr5.q=ifelse(chr5.q$chromosome==5,chr5.q$seg.gain.prop,NA)
    seg.gain.pro.chr5.q[length(seg.gain.pro.chr5.q)==0]=NA
    chr7.q=hyper.s[hyper.s$chromosome==7 & hyper.s$arm=="q",]
    seg.gain.pro.chr7.q=ifelse(chr7.q$chromosome==7,chr7.q$seg.gain.prop,NA)
    seg.gain.pro.chr7.q[length(seg.gain.pro.chr7.q)==0]=NA
    chr9.q=hyper.s[hyper.s$chromosome==9 & hyper.s$arm=="q",]
    seg.gain.pro.chr9.q=ifelse(chr9.q$chromosome==9,chr9.q$seg.gain.prop,NA)
    seg.gain.pro.chr9.q[length(seg.gain.pro.chr9.q)==0]=NA
    chr11.q=hyper.s[hyper.s$chromosome==11 & hyper.s$arm=="q",]
    seg.gain.pro.chr11.q=ifelse(chr11.q$chromosome==11,chr11.q$seg.gain.prop,NA)
    seg.gain.pro.chr11.q[length(seg.gain.pro.chr11.q)==0]=NA
    chr15.q=hyper.s[hyper.s$chromosome==15 & hyper.s$arm=="q",]
    seg.gain.pro.chr15.q=ifelse(chr15.q$chromosome==15,chr15.q$seg.gain.prop,NA)
    seg.gain.pro.chr15.q[length(seg.gain.pro.chr15.q)==0]=NA
    chr19.q=hyper.s[hyper.s$chromosome==19 & hyper.s$arm=="q",]
    seg.gain.pro.chr19.q=ifelse(chr19.q$chromosome==19,chr19.q$seg.gain.prop,NA)
    seg.gain.pro.chr19.q[length(seg.gain.pro.chr19.q)==0]=NA
    chr21.q=hyper.s[hyper.s$chromosome==21 & hyper.s$arm=="q",]
    seg.gain.pro.chr21.q=ifelse(chr21.q$chromosome==21,chr21.q$seg.gain.prop,NA)
    seg.gain.pro.chr21.q[length(seg.gain.pro.chr21.q)==0]=NA
    
    perc.gain.chr3=unique(hyper.s[hyper.s$chromosome==3,"perc.gain.chr"])
    perc.gain.chr3[length(perc.gain.chr3)==0]=NA
    perc.gain.chr5=unique(hyper.s[hyper.s$chromosome==5,"perc.gain.chr"])
    perc.gain.chr5[length(perc.gain.chr5)==0]=NA
    perc.gain.chr7=unique(hyper.s[hyper.s$chromosome==7,"perc.gain.chr"])
    perc.gain.chr7[length(perc.gain.chr7)==0]=NA
    perc.gain.chr9=unique(hyper.s[hyper.s$chromosome==9,"perc.gain.chr"])
    perc.gain.chr9[length(perc.gain.chr9)==0]=NA
    perc.gain.chr11=unique(hyper.s[hyper.s$chromosome==11,"perc.gain.chr"])
    perc.gain.chr11[length(perc.gain.chr11)==0]=NA
    perc.gain.chr15=unique(hyper.s[hyper.s$chromosome==15,"perc.gain.chr"])
    perc.gain.chr15[length(perc.gain.chr15)==0]=NA
    perc.gain.chr19=unique(hyper.s[hyper.s$chromosome==19,"perc.gain.chr"])
    perc.gain.chr19[length(perc.gain.chr19)==0]=NA
    perc.gain.chr21=unique(hyper.s[hyper.s$chromosome==21,"perc.gain.chr"])
    perc.gain.chr21[length(perc.gain.chr21)==0]=NA
    
    
    double.chrs=sum(hyper.s[,"arm.gain"]) # Number of gained/amplified arms in the sample
    Cytogenetic="Double_Genomes"
    line.double.p=c(sam,type,whole.median.CN,seg.gain.pro.chr3.p,seg.gain.pro.chr5.p,seg.gain.pro.chr7.p,seg.gain.pro.chr9.p,
                    seg.gain.pro.chr11.p,seg.gain.pro.chr15.p,seg.gain.pro.chr19.p,seg.gain.pro.chr21.p,double.chrs,Cytogenetic,"p",
                    perc.gain.chr3,perc.gain.chr5,perc.gain.chr7,perc.gain.chr9,perc.gain.chr11,perc.gain.chr15,perc.gain.chr19,perc.gain.chr21)
    line.double.q=c(sam,type,whole.median.CN,seg.gain.pro.chr3.q,seg.gain.pro.chr5.q,seg.gain.pro.chr7.q,seg.gain.pro.chr9.q,
                    seg.gain.pro.chr11.q,seg.gain.pro.chr15.q,seg.gain.pro.chr19.q,seg.gain.pro.chr21.q,double.chrs,Cytogenetic,"q",
                    perc.gain.chr3,perc.gain.chr5,perc.gain.chr7,perc.gain.chr9,perc.gain.chr11,perc.gain.chr15,perc.gain.chr19,perc.gain.chr21)
    hyper.1=rbind(hyper.1,line.double.p,line.double.q)
  }else if(sum(hyper.s[hyper.s$chromosome %in% chroms.hyper,"arm.gain"])>=2){
    sam=unique(hyper.s$sample)
    type=unique(hyper.s$type)
    whole.median.CN=median(as.numeric(as.character(hyper.s$CN.median)),na.rm = T)
    
    # . p arm
    chr3.p=hyper.s[hyper.s$chromosome==3 & hyper.s$arm=="p",]
    seg.gain.pro.chr3.p=ifelse(chr3.p$chromosome==3,chr3.p$seg.gain.prop,NA)
    seg.gain.pro.chr3.p[length(seg.gain.pro.chr3.p)==0]=NA
    chr5.p=hyper.s[hyper.s$chromosome==5 & hyper.s$arm=="p",]
    seg.gain.pro.chr5.p=ifelse(chr5.p$chromosome==5,chr5.p$seg.gain.prop,NA)
    seg.gain.pro.chr5.p[length(seg.gain.pro.chr5.p)==0]=NA
    chr7.p=hyper.s[hyper.s$chromosome==7 & hyper.s$arm=="p",]
    seg.gain.pro.chr7.p=ifelse(chr7.p$chromosome==7,chr7.p$seg.gain.prop,NA)
    seg.gain.pro.chr7.p[length(seg.gain.pro.chr7.p)==0]=NA
    chr9.p=hyper.s[hyper.s$chromosome==9 & hyper.s$arm=="p",]
    seg.gain.pro.chr9.p=ifelse(chr9.p$chromosome==9,chr9.p$seg.gain.prop,NA)
    seg.gain.pro.chr9.p[length(seg.gain.pro.chr9.p)==0]=NA
    chr11.p=hyper.s[hyper.s$chromosome==11 & hyper.s$arm=="p",]
    seg.gain.pro.chr11.p=ifelse(chr11.p$chromosome==11,chr11.p$seg.gain.prop,NA)
    seg.gain.pro.chr11.p[length(seg.gain.pro.chr11.p)==0]=NA
    chr15.p=hyper.s[hyper.s$chromosome==15 & hyper.s$arm=="p",]
    seg.gain.pro.chr15.p=ifelse(chr15.p$chromosome==15,chr15.p$seg.gain.prop,NA)
    seg.gain.pro.chr15.p[length(seg.gain.pro.chr15.p)==0]=NA
    chr19.p=hyper.s[hyper.s$chromosome==19 & hyper.s$arm=="p",]
    seg.gain.pro.chr19.p=ifelse(chr19.p$chromosome==19,chr19.p$seg.gain.prop,NA)
    seg.gain.pro.chr19.p[length(seg.gain.pro.chr19.p)==0]=NA
    chr21.p=hyper.s[hyper.s$chromosome==21 & hyper.s$arm=="p",]
    seg.gain.pro.chr21.p=ifelse(chr21.p$chromosome==21,chr21.p$seg.gain.prop,NA)
    seg.gain.pro.chr21.p[length(seg.gain.pro.chr21.p)==0]=NA
    # . q arm  
    chr3.q=hyper.s[hyper.s$chromosome==3 & hyper.s$arm=="q",]
    seg.gain.pro.chr3.q=ifelse(chr3.q$chromosome==3,chr3.q$seg.gain.prop,NA)
    seg.gain.pro.chr3.q[length(seg.gain.pro.chr3.q)==0]=NA
    chr5.q=hyper.s[hyper.s$chromosome==5 & hyper.s$arm=="q",]
    seg.gain.pro.chr5.q=ifelse(chr5.q$chromosome==5,chr5.q$seg.gain.prop,NA)
    seg.gain.pro.chr5.q[length(seg.gain.pro.chr5.q)==0]=NA
    chr7.q=hyper.s[hyper.s$chromosome==7 & hyper.s$arm=="q",]
    seg.gain.pro.chr7.q=ifelse(chr7.q$chromosome==7,chr7.q$seg.gain.prop,NA)
    seg.gain.pro.chr7.q[length(seg.gain.pro.chr7.q)==0]=NA
    chr9.q=hyper.s[hyper.s$chromosome==9 & hyper.s$arm=="q",]
    seg.gain.pro.chr9.q=ifelse(chr9.q$chromosome==9,chr9.q$seg.gain.prop,NA)
    seg.gain.pro.chr9.q[length(seg.gain.pro.chr9.q)==0]=NA
    chr11.q=hyper.s[hyper.s$chromosome==11 & hyper.s$arm=="q",]
    seg.gain.pro.chr11.q=ifelse(chr11.q$chromosome==11,chr11.q$seg.gain.prop,NA)
    seg.gain.pro.chr11.q[length(seg.gain.pro.chr11.q)==0]=NA
    chr15.q=hyper.s[hyper.s$chromosome==15 & hyper.s$arm=="q",]
    seg.gain.pro.chr15.q=ifelse(chr15.q$chromosome==15,chr15.q$seg.gain.prop,NA)
    seg.gain.pro.chr15.q[length(seg.gain.pro.chr15.q)==0]=NA
    chr19.q=hyper.s[hyper.s$chromosome==19 & hyper.s$arm=="q",]
    seg.gain.pro.chr19.q=ifelse(chr19.q$chromosome==19,chr19.q$seg.gain.prop,NA)
    seg.gain.pro.chr19.q[length(seg.gain.pro.chr19.q)==0]=NA
    chr21.q=hyper.s[hyper.s$chromosome==21 & hyper.s$arm=="q",]
    seg.gain.pro.chr21.q=ifelse(chr21.q$chromosome==21,chr21.q$seg.gain.prop,NA)
    seg.gain.pro.chr21.q[length(seg.gain.pro.chr21.q)==0]=NA
    
    perc.gain.chr3=unique(hyper.s[hyper.s$chromosome==3,"perc.gain.chr"])
    perc.gain.chr3[length(perc.gain.chr3)==0]=NA
    perc.gain.chr5=unique(hyper.s[hyper.s$chromosome==5,"perc.gain.chr"])
    perc.gain.chr5[length(perc.gain.chr5)==0]=NA
    perc.gain.chr7=unique(hyper.s[hyper.s$chromosome==7,"perc.gain.chr"])
    perc.gain.chr7[length(perc.gain.chr7)==0]=NA
    perc.gain.chr9=unique(hyper.s[hyper.s$chromosome==9,"perc.gain.chr"])
    perc.gain.chr9[length(perc.gain.chr9)==0]=NA
    perc.gain.chr11=unique(hyper.s[hyper.s$chromosome==11,"perc.gain.chr"])
    perc.gain.chr11[length(perc.gain.chr11)==0]=NA
    perc.gain.chr15=unique(hyper.s[hyper.s$chromosome==15,"perc.gain.chr"])
    perc.gain.chr15[length(perc.gain.chr15)==0]=NA
    perc.gain.chr19=unique(hyper.s[hyper.s$chromosome==19,"perc.gain.chr"])
    perc.gain.chr19[length(perc.gain.chr19)==0]=NA
    perc.gain.chr21=unique(hyper.s[hyper.s$chromosome==21,"perc.gain.chr"])
    perc.gain.chr21[length(perc.gain.chr21)==0]=NA
    
    double.chrs=sum(hyper.s[hyper.s$chromosome %in% chroms.hyper,"arm.gain"])
    Cytogenetic="Trisomie"
    line.double.p=c(sam,type,whole.median.CN,seg.gain.pro.chr3.p,seg.gain.pro.chr5.p,seg.gain.pro.chr7.p,seg.gain.pro.chr9.p,
                    seg.gain.pro.chr11.p,seg.gain.pro.chr15.p,seg.gain.pro.chr19.p,seg.gain.pro.chr21.p,double.chrs,Cytogenetic,"p",
                    perc.gain.chr3,perc.gain.chr5,perc.gain.chr7,perc.gain.chr9,perc.gain.chr11,perc.gain.chr15,perc.gain.chr19,perc.gain.chr21)
    line.double.q=c(sam,type,whole.median.CN,seg.gain.pro.chr3.q,seg.gain.pro.chr5.q,seg.gain.pro.chr7.q,seg.gain.pro.chr9.q,
                    seg.gain.pro.chr11.q,seg.gain.pro.chr15.q,seg.gain.pro.chr19.q,seg.gain.pro.chr21.q,double.chrs,Cytogenetic,"q",
                    perc.gain.chr3,perc.gain.chr5,perc.gain.chr7,perc.gain.chr9,perc.gain.chr11,perc.gain.chr15,perc.gain.chr19,perc.gain.chr21)
    hyper.1=rbind(hyper.1,line.double.p,line.double.q)
  }else{
    no.hyperdiploid.1=rbind(no.hyperdiploid.1,samples[iii])
  }
}

hyper1=as.data.frame(hyper.1)
colnames(hyper1)=c("sample","type","whole.median.CN","seg.gain.pro.chr3","seg.gain.pro.chr5","seg.gain.pro.chr7","seg.gain.pro.chr9",
                   "seg.gain.pro.chr11","seg.gain.pro.chr15","seg.gain.pro.chr19","seg.gain.pro.chr21","double.chrs","Cytogenetic","arm",
                   "perc.gain.chr3","perc.gain.chr5","perc.gain.chr7","perc.gain.chr9","perc.gain.chr11","perc.gain.chr15","perc.gain.chr19","perc.gain.chr21"
)
rownames(hyper1)=c(1:length(hyper1$sample))
table(hyper1$type,hyper1$Cytogenetic)
hyper1
for(q in c(4:11,15:22)){
  hyper1[,q]=as.numeric(as.character(hyper1[,q]))
}

SAMPLE=unique(hyper1$sample)
hyper2=NULL
for(qq in 1:length(SAMPLE)){
  hyper1.sample=hyper1[hyper1$sample==SAMPLE[qq],]
  sam=SAMPLE[qq]
  type=unique(hyper1.sample$type)
  whole.median.CN=unique(hyper1.sample$whole.median.CN)
  
  hyp.chr3=unique(ifelse(hyper1.sample$perc.gain.chr3>=0.6,1,0))
  hyp.chr5=unique(ifelse(hyper1.sample$perc.gain.chr5>=0.6,1,0))
  hyp.chr7=unique(ifelse(hyper1.sample$perc.gain.chr7>=0.6,1,0))
  hyp.chr9=unique(ifelse(hyper1.sample$perc.gain.chr9>=0.6,1,0))
  hyp.chr11=unique(ifelse(hyper1.sample$perc.gain.chr11>=0.6,1,0))
  hyp.chr15=unique(ifelse(hyper1.sample$perc.gain.chr15>=0.6,1,0))
  hyp.chr19=unique(ifelse(hyper1.sample$perc.gain.chr19>=0.6,1,0))
  hyp.chr21=unique(ifelse(hyper1.sample$perc.gain.chr21>=0.6,1,0))
  
  hyp.sample=ifelse(sum(hyp.chr3,hyp.chr5,hyp.chr7,hyp.chr9,hyp.chr11,hyp.chr15,hyp.chr19,hyp.chr21,na.rm = T)>=2,"Hyperdiploidy","No.HRD")
  line=c(sam,type,whole.median.CN,
         hyp.chr3,hyp.chr5,hyp.chr7,hyp.chr9,hyp.chr11,hyp.chr15,hyp.chr19,hyp.chr21,
         hyp.sample)
  hyper2=rbind(line,hyper2)
}
hyper.2=as.data.frame(hyper2)
colnames(hyper.2)=c("sample","type","whole.median.CN",
                    "hyp.chr3","hyp.chr5","hyp.chr7","hyp.chr9","hyp.chr11","hyp.chr15","hyp.chr19","hyp.chr21",
                    "hyp.sample")
rownames(hyper.2)=c(1:length(hyper.2$sample))

#plot(density(hyper1$seg.gain.median))
#plot(density(hyper1$seg.gain.sum))
#hyper1[hyper1$seg.gain.median<0.80,]
#hyper1[hyper1$seg.gain.sum<1.60,] # 0.80 per two chromosomes (hyperdiploid is when at least two chromosome are gained)
#hyper1=hyper1[hyper1$seg.gain.sum>=SUM,]

no.hyper=as.data.frame(no.hyperdiploid)
if(length(no.hyper)>0){
  colnames(no.hyper)="sample"  
}
no.hyper1=as.data.frame(no.hyperdiploid.1)
if(length(no.hyper1)>0){
  colnames(no.hyper1)="sample"  
}


no.hyper.1=rbind(no.hyper,no.hyper1)
no.hyper.1$Cytogenetic="No_HRD"

hyper.5cohorts=hyper.2
summary(hyper.5cohorts)
no.hyper.5cohorts=no.hyper.1


hyper.def=hyper.5cohorts
no.hyper.def=no.hyper.5cohorts


write.table(
  hyper.def,
  paste(
    OUT.DIR,
    PERCENTAGE*100,
    "%_ARMS_",
    "hyperdiploid.txt",
    sep=""
  ),
  col.names = T,
  row.names = F,
  quote = F,
  sep="\t"
)

write.table(
  no.hyper.def,
  paste(
    OUT.DIR,
    PERCENTAGE*100,
    "%_ARMS_",
    "no.hyperdiploid.txt",
    sep=""
  ),
  col.names = T,
  row.names = F,
  quote = F,
  sep="\t"
)













