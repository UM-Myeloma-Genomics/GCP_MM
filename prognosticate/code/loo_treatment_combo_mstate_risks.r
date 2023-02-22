rm(list=ls())

library(plyr)
library(readr)
suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()

parser$add_argument("--model", default='neural_cox_non_prop', type="character",
                    help="model name")

parser$add_argument("--sample_id", default='PD5886a', type="character",
                    help="model name")

parser$add_argument("--code_path", default='/Users/axr2376/Desktop/pred_1_0_paper/code/', type="character",
                    help="for code.r")

parser$add_argument("--base_path", default='/Users/axr2376/Desktop/pred_1_0_paper/data_out/expt_1/loo/model_probs/', type="character",
                    help="folds count")

parser$add_argument("--output_path", default='/Users/axr2376/Desktop/pred_1_0_paper/data_out/expt_1/loo/state_probs/', type="character",
                    help="folds count")

args <- parser$parse_args()

base_path = args$base_path

output_path = args$output_path

model = args$model

sample_id = args$sample_id

pattern_files= paste(sample_id, '.*~', model, '.csv', sep='')
model_indiv_pred_files = list.files(path=base_path, pattern=pattern_files)

for (x in 1:length(model_indiv_pred_files)) {
  Data = read.csv(  paste(base_path, model_indiv_pred_files[x], sep = ''))[-1,]

  Sample = Data[c("sample")]
  days = Data$Time..years.*365
  SP2 = Data$Move.to.phase.2
  SPODinIND = Data$Progress..P1. # Survival for move to POD in induction (2nd column of matrix)
  SDeathafterPODinIND = Data$Progress...deceased..P1. # Survival for death after POD in induction (4th column of matrix)
  
  source(paste(args$code_path, 'code.r', sep=''))
  
  PrPODinInduction =  calculatePODInduction(days,SPODinIND,SP2)
  PrmoveP2 =  calculateMToP2(days,SP2,SPODinIND)
  PrAliveInduction = claculateAliveInInduction(PrPODinInduction,PrmoveP2)
  PrDeathinPODinInduction =  calculateDeathInPODInduction(days,SDeathafterPODinIND,days,PrPODinInduction)
  L = length(PrPODinInduction)
  PrProgressioninP2X = rep(PrPODinInduction[L],length(PrDeathinPODinInduction))
  PrProgressioninP2X[1:L] = PrPODinInduction
  
  PrAliveinPODinInduction = PrProgressioninP2X-PrDeathinPODinInduction
  
  SDP2 = Data$Non.progress...deceased..P2.# Survival for death in P2
  
  SP2P2 = Data$Progress..P2.# Survival for progression from P2
  
  PrDeathinP2 =  calculateDeathInP2(days,SDP2,SP2P2,days,PrmoveP2)
  
  PrProgressioninP2 = calculateToProgression(days,SDP2,SP2P2,days,PrmoveP2)
  
  PAliveinP2 =  claculateAliveInP2(PrmoveP2,PrDeathinP2,PrProgressioninP2)
  
  SDpostP = Data$Progress...deceased..P2.# Survival for death in progression
  PDeathpostPr = calculateDeathfromProgression(days,SDpostP,days,SDP2,SP2P2,days,PrmoveP2)
  L = length(PrProgressioninP2)
  PrProgressioninP2X = rep(PrProgressioninP2[L],length(PDeathpostPr))
  PrProgressioninP2X[1:L] = PrProgressioninP2
  
  XX = PrProgressioninP2X-PDeathpostPr
  
  Days = 1:(length(days)*3)
  PrDeathinP2E = rep(PrDeathinP2[length(PrDeathinP2)],length(Days))
  PrDeathinP2E[1:length(PrDeathinP2)] = PrDeathinP2
  PAliveinP2E = rep(PAliveinP2[length(PAliveinP2)],length(Days))
  PAliveinP2E[1:length(PAliveinP2)] = PAliveinP2
  PrDeathinPODinInductionE = rep(PrDeathinPODinInduction[length(PrDeathinPODinInduction)],length(Days))
  PrDeathinPODinInductionE[1:length(PrDeathinPODinInduction)] =PrDeathinPODinInduction
  PrAliveinPODinInductionE = rep(PrAliveinPODinInduction[length(PrAliveinPODinInduction)],length(Days))
  PrAliveinPODinInductionE[1:length(PrAliveinPODinInduction)] =PrAliveinPODinInduction
  PrAliveInductionE = rep(PrAliveInduction[length(PrAliveInduction)],length(Days))
  PrAliveInductionE[1:length(PrAliveInduction)] = PrAliveInduction
  
  SaveData = cbind(Sample, Days,PrAliveInductionE,PrAliveinPODinInductionE,PAliveinP2E,XX,PrDeathinPODinInductionE,PrDeathinP2E,PDeathpostPr)
  SaveData2 = cbind(Days,rowSums(SaveData[,6:8]),1-rowSums(SaveData[,6:8]))
  colnames(SaveData2) =  c('Time (Days)','Death','Alive')
  colnames(SaveData) = c('sample', 'Time (Days)','Alive in Induction', 'Alive in POD (induction)', 'Alive in Phase 2','Alive in POD','Death in POD (induction)','Death in Phase 2','Death in POD (Phase 2)')
  write.csv(SaveData,row.names = FALSE,file=paste(output_path, model_indiv_pred_files[x], sep = ''))

}
  

