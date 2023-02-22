#!/bin/bash

# Description: incrementally add feature groups
# Out: 
#    model/multistate predictions, 
#    c-index of model predictions

# expt. settings
folds=5
repeats=10
total_splits=$(($folds * $repeats))

#### set paths/dirs
current_dir=$(pwd)
path_data_in="$current_dir/data_in"
path_data_legend="$path_data_in/legend_PMMM_2022.xlsx"

out_dir="/home/axr2376/pred_1_0_paper"
path_data_out="$out_dir/data_out/expt_1"
path_permut="$path_data_out/permut_feat_imps"
path_permut_feat_compose="$path_permut/compose/"
path_permut_feat_agg="$path_permut/aggregate/"
path_permut_feat_agg_reranked="$path_permut/aggregate_reranked/"
path_model_preds="$path_data_out/kfold/model_probs/"

path_os_multistate_preds="$path_data_out/kfold/state_probs/os/"
path_efs_multistate_preds="$path_data_out/kfold/state_probs/efs/"
path_all_states_multistate_preds="$path_data_out/kfold/state_probs/all_states/"

path_plots="$path_data_out/plots"

path_loo_model_preds="$path_data_out/loo/model_probs/"
path_loo_state_preds="$path_data_out/loo/state_probs/"
path_treat_preds="$path_data_out/preds_treatment"


mkdir -p $out_dir $path_data_out $path_plots $path_permut $path_permut_feat_agg $path_permut_feat_compose "$path_data_out/feat_matrix" $path_model_preds $path_os_multistate_preds $path_efs_multistate_preds "$path_data_out/kfold/metrics/model_feat_combo" "$path_data_out/kfold/metrics/model_genomics" $path_permut_feat_agg $path_permut_feat_compose "$path_data_out/shap_feat_imps" $path_loo_model_preds $path_loo_state_preds $path_treat_preds $path_permut_feat_agg_reranked $path_all_states_multistate_preds

# func
task(){
  current_dir=$(pwd)
  group_id=${1}
  folds=${2}
  repeats=${3}
  model=${4}
  total_splits=${5}
  ext_val=${6}

  if [ "$group_id" = "100" ]; then
    py_group_id="None"
  else
    py_group_id=$group_id
  fi
 
  # if predictions on external test set
  if [ "$ext_val" = "Yes" ]; then
    CUDA_VISIBLE_DEVICES="" python3 "$current_dir/code/kfold_model_preds_ext_study.py" --get_feats_group $py_group_id --surv_model $model --kfolds $folds --n_repeat_kfold $repeats --in_data_file $path_data_in --path $path_data_out --legend $path_data_legend
    Rscript "$current_dir/code/inference_mstate_risks_ext_study.r" --model $model --num_folds $total_splits --group $group_id --base_path $path_model_preds --code_path "$current_dir/code/" --output_path_os $path_os_multistate_preds --output_path_efs $path_efs_multistate_preds --output_path_all $path_all_states_multistate_preds 
  else
    CUDA_VISIBLE_DEVICES="" python3 "$current_dir/code/kfold_model_train_preds.py" --get_feats_group $py_group_id --surv_model $model --kfolds $folds --n_repeat_kfold $repeats --in_data_file $path_data_in --path $path_data_out --legend $path_data_legend 
    Rscript "$current_dir/code/kfold_mstate_risks.r" --model $model --num_folds $total_splits --group $group_id --base_path $path_model_preds --code_path "$current_dir/code/" --output_path_os $path_os_multistate_preds --output_path_efs $path_efs_multistate_preds --output_path_all $path_all_states_multistate_preds &

  fi

}

# func
get_best_gen_feats(){
  folds=$1
  repeats=$2
  model=$3
  top_n_gen_feats=$4
  
  # univariate + permutation feat imp ranking
  CUDA_VISIBLE_DEVICES="" python3 "$current_dir/code/kfold_model_train_preds.py" --surv_model $model --top_gen_feats --top_n_gen_feats $top_n_gen_feats --kfolds $folds --n_repeat_kfold $repeats --in_data_file $path_data_in --path $path_data_out --legend $path_data_legend --agg_reranked
}

###### All feats
# model for feat sel; gen add-on
model="neural_cox_non_prop"
ext_val="No"

task "100" "$folds" "$repeats" "$model" "$total_splits" "$ext_val" &
task "100" "$folds" "$repeats" "rsf" "$total_splits" "$ext_val" &
task "100" "$folds" "$repeats" "cph" "$total_splits" "$ext_val" &

###### Feat rankers: (1) Permutation Importance
CUDA_VISIBLE_DEVICES=0 python3 "$current_dir/code/kfold_model_train_preds.py" --surv_model $model --kfolds $folds --n_repeat_kfold $repeats --in_data_file $path_data_in --path $path_data_out --legend $path_data_legend --perm_feat_imp &
#perm_imp_pid=$! 

### Feature rankers : (2) SHAP
#CUDA_VISIBLE_DEVICES=1 python3 "$current_dir/code/main_kfold_train_preds.py" --surv_model $model --kfolds $folds --n_repeat_kfold $repeats --in_data_file $path_data_in --path $path_data_out --legend $path_data_legend --shap_feat_imp &

###### ISS/R-ISS/R2-ISS

# python3 "$current_dir/code/create_feature_combo_matrix.py" --path $path_data_out --legend $path_data_legend --permut_feat_path $path_permut_feat_agg --iss_feat_groups

task "0" "$folds" "$repeats" "cph" "$total_splits" "$ext_val" &
task "1" "$folds" "$repeats" "cph" "$total_splits" "$ext_val" &
task "2" "$folds" "$repeats" "cph" "$total_splits" "$ext_val" &

wait $perm_imp_pid

### Aggregate permutation importance rank scores

python3 "$current_dir/code/metrics.py" --surv_model $model --in_data_file $path_data_in --path $path_data_out --agg_perm_imp --re_rank_feats --legend $path_data_legend

### Top GEN feats
N=11 # total 136 gen features (without translocs)
(
for num_gen_feats in $(seq 0 10); do 
   ((i=i%N)); ((i++==0)) && wait
   get_best_gen_feats "$folds" "$repeats" "$model" "$num_gen_feats" &
done
)

python3 "$current_dir/code/metrics.py" --surv_model $model --in_data_file $path_data_in --path $path_data_out --get_top_gen_feats --heatmap_bygroup --re_rank_feats --legend $path_data_legend

###### Aggregate --> score --> rank
# create feature combinations to feed the right slices of features to train models

: '
python3 "$current_dir/code/create_feature_combo_matrix.py" --path $path_data_out --legend $path_data_legend --permut_feat_path_reranked $path_permut_feat_agg_reranked --iss_feat_groups --all_feat_groups --permut_feat_path $path_permut_feat_agg --reverse_order

task "10" "$folds" "$repeats" "cph" "$total_splits" "$ext_val" & 

master_file="$path_data_out/feat_matrix/main_keeper.csv"
total_jobs_count=$(wc -l $master_file | awk '{print $1}') # with header

# without header and python is zero indexed
total_jobs_count="$((total_jobs_count-2))"

N=11 # except group 10 (GEP)
(
for group_id in $(seq 3 $total_jobs_count); do 
   ((i=i%N)); ((i++==0)) && wait
   if [ "$group_id" != "10" ]; then
     task "$group_id" "$folds" "$repeats" "$model" "$total_splits" "$ext_val" &
   fi

done 
)
'

: '
# external test set
task "9" "$folds" "$repeats" "$model" "$total_splits" "Yes" &
task "2" "$folds" "$repeats" "cph" "$total_splits" "Yes" &
task "100" "$folds" "$repeats" "$model" "$total_splits" "Yes" &
task "100" "$folds" "$repeats" "rsf" "$total_splits" "Yes" &
task "100" "$folds" "$repeats" "cph" "$total_splits" "Yes" &
task "0" "$folds" "$repeats" "cph" "$total_splits" "Yes" &
task "1" "$folds" "$repeats" "cph" "$total_splits" "Yes" &
'

: '
###### Treatment predictions
# loo predictions
indiv_treat_preds(){
  CUDA_VISIBLE_DEVICES="" python3 "$current_dir/code/loo_treatment_combo_model_preds.py" --sample $1 --in_data_file $path_data_in --path $path_data_out --legend $path_data_legend --get_feats_group $2
  Rscript "$current_dir/code/loo_treatment_combo_mstate_risks.r" --sample_id $1 --base_path $path_loo_model_preds --code_path "$current_dir/code/" --output_path $path_loo_state_preds --model $model
}

# python3 "$current_dir/code/viz_diff_combo.py" --in_data_file $path_data_in --path $path_data_out  --get_treat_stats

# predictions on specified patient IDs
filename="$path_treat_preds/patients.txt"
IFS=$'\n' read -d '' -r -a pat_ids < $filename
printf "%s\n" "${pat_ids[@]}"

N=16
total_pats=${#pat_ids[@]}
echo $total_pats

(
for num_pat in $(seq 0 $total_pats); do 
   ((i=i%N)); ((i++==0)) && wait
   echo "${pat_ids[num_pat]}"
   indiv_treat_preds "${pat_ids[num_pat]}" "9" &
done
)
'
