##################################################################################
##################################################################################
##############################################################################1###
#                                                                                #
#        .     .     .                                               .     .     .
####     .     .     .   PredictionModel - Mutational Signatures   .     .     .     ####
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
ITER=10 # mmsig parameter: 1000 iterations recommended for stable results


##################################################################################
##################################################################################
library(GenomicRanges)
library(IRanges)
library(data.table)
library(dplyr)
library(tidyr)
library(plyr)
library(stringr)
library(tibble)
library(ggplot2)
library(MutationalPatterns)
library(deconstructSigs)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(mmsig)
library(BSgenome.Hsapiens.UCSC.hg19)
library(clonevol)
set.seed(1)

rotatedAxisElementText = function(angle,position='x'){
  angle     = angle[1];
  position  = position[1]
  positions = list(x=0,y=90,top=180,right=270)
  if(!position %in% names(positions))
    stop(sprintf("'position' must be one of [%s]",paste(names(positions),collapse=", ")),call.=FALSE)
  if(!is.numeric(angle))
    stop("'angle' must be numeric",call.=FALSE)
  rads  = (angle - positions[[ position ]])*pi/180
  hjust = 0.5*(1 - sin(rads))
  vjust = 0.5*(1 + cos(rads))
  element_text(angle=angle,vjust=vjust,hjust=hjust)
}

scale_fill_sigs <- function(...){
  ggplot2:::manual_scale('fill',
                         values = setNames(c(RColorBrewer::brewer.pal(8, "Dark2"),"lightpink",
                                             "firebrick","lightblue1","darkslategray4"),
                                           c("SBS1", "SBS2", "SBS5","SBS8", "SBS9","SBS13","SBS18","SBS-MM1",
                                             "SBS17b","SBS31","SBS35","E_SBS37")),
                         ...)
}

# . . . - Signatures definition
sig_ref <- signature_ref[c("class","sub", "tri", "SBS1", "SBS2", "SBS5", "SBS8", 
                           "SBS9", "SBS13", "SBS18")]
sig<-sig_ref[,-1]
sig<-sig[,c("sub", "tri","SBS1","SBS2","SBS13","SBS5", "SBS8", 
            "SBS9", "SBS18")]


########################################################################################################################

dnds.input=read.delim(paste(WF,"/input/snv/snv_for_mutationalsignatures_test.txt",sep=""),
                      stringsAsFactors = F)
sig_out100 <- mm_fit_signatures(muts.input=dnds.input,
                              sig.input=sig,
                              input.format = "vcf",
                              strandbias = T,
                              bootstrap = TRUE,
                              iterations = ITER,
                              refcheck=FALSE,
                              cos_sim_threshold = 0.01,
                              force_include = c("SBS1","SBS5"),
                              dbg=FALSE)


sig_order = c("SBS1", "SBS2", "SBS13", "SBS5", "SBS8", "SBS9","SBS17b", 
              "SBS18","SBS-MM1","SBS31","E_SBS25","SBS35","E_SBS37")
est<- sig_out100$estimate

pr=est[est$SBS2!=0 | est$SBS13!=0,]
bootSigsPlot(sig_out100$bootstrap)

apobec=sig_out100$bootstrap[sig_out100$bootstrap$signature %in% c("SBS2","SBS13"),]

samples=unique(apobec$sample)
APOBEC=NULL
for(i in 1:length(samples)){
  sub=apobec[apobec$sample==samples[i],]
  estimate=sub[1,6]+sub[2,6]
  row=c(samples[i],estimate)
  APOBEC=rbind(APOBEC,row)
}
APOBEC=as.data.frame(APOBEC)
rownames(APOBEC)=1:length(APOBEC$V1)
colnames(APOBEC)=c("sample","apobec")
APOBEC$apobec=as.numeric(as.character(APOBEC$apobec))

summary(APOBEC$apobec)


APOBEC$offset=0
rownames(APOBEC)=APOBEC$sample
APOBEC=APOBEC[,-1]

write.table(APOBEC,
            paste(WF,"/output/APOBEC_relative_contribution_matrix.txt",sep=""),
            col.names = T,
            row.names = T,
            quote = F,
            sep = "\t")
