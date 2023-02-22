rm(list=ls())

library(plyr)
library(readr)
suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()

parser$add_argument("--model", default='neural_cox_non_prop', type="character",
                    help="model name")
parser$add_argument("--num_folds", default=5, type="integer",
                    help="folds count")
parser$add_argument("--group", default='7', type="character",
                    help="folds count")

parser$add_argument("--code_path", default='/Users/axr2376/Desktop/pred_1_0_paper/code/', type="character",
                    help="for code.r")

parser$add_argument("--base_path", default='/Users/axr2376/Desktop/pred_1_0_paper/data_out/expt_1/kfold/model_probs/', type="character",
                    help="folds count")

parser$add_argument("--output_path_os", default='/Users/axr2376/Desktop/pred_1_0_paper/data_out/expt_1/kfold/state_probs/os/', type="character",
                    help="output path for OS")

parser$add_argument("--output_path_efs", default='/Users/axr2376/Desktop/pred_1_0_paper/data_out/expt_1/kfold/state_probs/efs/', type="character",
                    help="output path for EFS")

parser$add_argument("--output_path_all", default='/Users/axr2376/Desktop/pred_1_0_paper/data_out/expt_1/kfold/predictions_test_set/', type="character",
                    help="all states")

args <- parser$parse_args()

base_path = args$base_path

output_path_os = args$output_path_os

output_path_efs = args$output_path_efs

output_path_all = args$output_path_all

num_folds = args$num_folds #indexed from 0

model = args$model

group = args$group

for (fold in 1: num_folds-1) {
    preds_groups_per_fold = paste('group_', group, '_fold_', fold, '_', model, sep='')
    model_indiv_pred_files = list.files(path=base_path, pattern=preds_groups_per_fold)

    prog_p1 = read.csv(paste(base_path, 'group_', group, '_fold_', fold, '_', model, '_Progress (P1).csv', sep = ''))[-1,]
    prog_dec_p1 = read.csv(paste(base_path, 'group_', group, '_fold_', fold, '_', model, '_Progress & deceased (P1).csv', sep = ''))[-1,]
    move_to_phase_2 = read.csv(paste(base_path, 'group_', group, '_fold_', fold, '_', model, '_Move to phase 2.csv', sep = ''))[-1,]
    prog_p2 = read.csv(paste(base_path, 'group_', group, '_fold_', fold, '_', model, '_Progress (P2).csv', sep = ''))[-1,]
    prog_dec_p2 = read.csv(paste(base_path, 'group_', group, '_fold_', fold, '_', model, '_Progress & deceased (P2).csv', sep = ''))[-1,]
    non_prog_p2 = read.csv(paste(base_path, 'group_', group, '_fold_', fold, '_', model, '_Non-progress & deceased (P2).csv', sep = ''))[-1,]

    Data = data.frame()
    adj_Data2OS = NULL
    adj_Data2EFS = NULL
    for (ind in 1: dim(prog_p1)[2]) {
      Data = prog_p1[,ind, drop=FALSE]
      Sample = colnames(Data)[1]
      colnames(Data)[1] = "Progress..P1."
      Data['Progress...deceased..P1.'] = prog_dec_p1[,ind, drop=FALSE]
      Data['Move.to.phase.2'] = move_to_phase_2[,ind, drop=FALSE]
      Data['Progress..P2.'] = prog_p2[,ind, drop=FALSE]
      Data['Progress...deceased..P2.'] = prog_dec_p2[,ind, drop=FALSE]
      Data['Non.progress...deceased..P2.'] = non_prog_p2[,ind, drop=FALSE]
      Data['Time..years.'] = seq(1, length.out=nrow(Data), by=1) / 365

      # start
      days = Data$Time..years. *365
      SP2 = Data$Move.to.phase.2
      SPODinIND = Data$Progress..P1. # Survival for move to POD in induction (2nd column of matrix)
      SDeathafterPODinIND = Data$Progress...deceased..P1. # Survival for death after POD in induction (4th column of matrix)
      SDP2 = Data$Non.progress...deceased..P2.# Survival for death in P2
      SP2P2 = Data$Progress..P2.# Survival for progression from P2
      SDpostP = Data$Progress...deceased..P2.# Survival for death in progression
      
      source(paste(args$code_path, 'code.r', sep=''))
      #P_POD(t): calculate probability of POD in induction by time t while adjusting for competing event of transitioning to phase 2 (i.e. responding to a treatment)
      PrPODinInduction =  calculatePODInduction(days,SPODinIND,SP2)
      #P_P2(t): calculate probability of transitioning to phase 2 (i.e. responding to a treatment)  by time t while adjusting for competing event of POD in induction
      PrmoveP2 =  calculateMToP2(days,SP2,SPODinIND)
      # calculate probability of being alive in induction (i.e. receiving induction treatment): P_I(t) = 1-P_POD(t) - P_P2(r)
      PrAliveInduction = claculateAliveInInduction(PrPODinInduction,PrmoveP2)
      # P_Death_in_POD(t):calculate probability of death in POD that happened in induction: function of P_POD(t) and dying afterwards (SDeathafterPODinIND)
      PrDeathinPODinInduction =  calculateDeathInPODInduction(days,SDeathafterPODinIND,days,PrPODinInduction)
      L = length(PrPODinInduction)
      PrProgressioninP2X = rep(PrPODinInduction[L],length(PrDeathinPODinInduction))
      PrProgressioninP2X[1:L] = PrPODinInduction
      # P_Alive_in_POD(t): probability of being alive in POD, equivalent to being in POD but not dead in POD
      PrAliveinPODinInduction = PrProgressioninP2X-PrDeathinPODinInduction
      
      # P_Death_in_P2(t): probability of dying without POD by time t in phase 2 part (P2). Competing event progression in P2, prerequisite is to transition to P2
      PrDeathinP2 =  calculateDeathInP2(days,SDP2,SP2P2,days,PrmoveP2)
      # P_POD_in_P2(t): probability of POD by time t in phase 2 part (P2). Competing event death in P2, prerequisite is to transition to P2
      PrProgressioninP2 = calculateToProgression(days,SDP2,SP2P2,days,PrmoveP2)
      #  P_Alive_in_P2(t): probability of alive by time t in phase 2 part (P2): equivalent not to POD or die without POD, prerequisite is to transition to P2
      PAliveinP2 =  claculateAliveInP2(PrmoveP2,PrDeathinP2,PrProgressioninP2)
      
      #  P_Death_in_POD_P2(t): probability of dying by time t in POD that happened after phase 2 part (P2), prerequisite is to transition to P2
      PDeathpostPr = calculateDeathfromProgression(days,SDpostP,days,SDP2,SP2P2,days,PrmoveP2)
      L = length(PrProgressioninP2)
      PrProgressioninP2X = rep(PrProgressioninP2[L],length(PDeathpostPr))
      PrProgressioninP2X[1:L] = PrProgressioninP2
      
      #  P_Death_in_POD_P2(t): probability of being alive by time t in POD that happened after phase 2 part (P2),  prerequisite is to transition to POD in P2
      PrAlive_POST_PROG_P2 = PrProgressioninP2X-PDeathpostPr
      
      #
      # This part is to impute states outside of the observation window
      #
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
      
      # All states
      SaveData = cbind(Days,PrAliveInductionE,PrAliveinPODinInductionE,PAliveinP2E,PrAlive_POST_PROG_P2,PrDeathinPODinInductionE,PrDeathinP2E,PDeathpostPr)
      colnames(SaveData) = c('Time (Days)','Alive in Induction', 'Alive in POD (induction)', 'Alive in Phase 2','Alive in POD','Death in POD (induction)','Death in Phase 2','Death in POD (Phase 2)')

      # OS
      SaveData2OS = cbind(Days,rowSums(SaveData[,6:8]),1-rowSums(SaveData[,6:8]))
      colnames(SaveData2OS) =  c('Time (Days)', 'Death', 'Alive')

      # EFS
      prr = PrPODinInduction[length(PrPODinInduction)]
      PrPODinInductionX = rep(prr,length(PrProgressioninP2))
      PrPODinInductionX[1:length(PrPODinInduction)] = PrPODinInduction
      DaysX = Days[1:length(PrProgressioninP2)]
    
      SaveDataEFS = cbind(DaysX,PrPODinInductionX,PrProgressioninP2)
      SaveData2EFS = cbind(DaysX,rowSums(SaveDataEFS),1-rowSums(SaveDataEFS[,2:3]))
      colnames(SaveData2EFS) =  c('Time (Days)', 'POD', 'NON_POD')
      
      # save
      adj_Data2EFS = cbind(adj_Data2EFS, c(Sample,SaveData2EFS[,'NON_POD']))
      adj_Data2OS = cbind(adj_Data2OS, c(Sample,SaveData2OS[,'Alive']))

      # write all states only for model with +Top Genomics (IRMM)
      # for calculation of mean/median risks across folds
      if (group == '9' && model == 'neural_cox_non_prop') {
        write.csv(SaveData[720:1830,],row.names = FALSE,file=paste(output_path_all,'group~', group, '~sample~', Sample, '~fold~', fold, '~', model,'.csv',sep=''))
      }

    }

    # OS
    colnames(adj_Data2OS) = adj_Data2OS[1,]
    adj_Data2OS = adj_Data2OS[-1, ]
    write.csv(adj_Data2OS,row.names = FALSE,file=paste(output_path_os,'group-', group, '-fold-', fold, '-', model,'.csv',sep=''))
    
    # EFS
    colnames(adj_Data2EFS) = adj_Data2EFS[1,]
    adj_Data2EFS = adj_Data2EFS[-1, ]
    write.csv(adj_Data2EFS,row.names = FALSE,file=paste(output_path_efs,'group-', group, '-fold-', fold, '-', model,'.csv',sep=''))
    
  }


