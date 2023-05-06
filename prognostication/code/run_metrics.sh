#!/bin/bash

# Description: Compute  metrics 

# expt. settings
folds=5
repeats=10
total_splits=$(($folds * $repeats))

#### set paths/dirs
current_dir=$(pwd)
path_data_in="$current_dir/data_in"
path_data_legend="$path_data_in/legend_PMMM_2023.xlsx"

out_dir="/home/axr2376/pred_1_0_paper"
path_data_out="$out_dir/data_out/expt_1_matrix_04052023"

path_permut="$path_data_out/permut_feat_imps"
path_permut_feat_compose="$path_permut/compose/"
path_permut_feat_agg="$path_permut/aggregate/"
path_model_preds="$path_data_out/kfold/model_probs/"
path_plots="$path_data_out/plots"

time=1825

c_index_type="C_index_harrell"

: '
## section 1
##### c-index (cross-validated)
python3 $(pwd)/code/metrics.py --preds_time_in_days $time --in_data_file $path_data_in --path $path_data_out --outcome "os" --group_feat_add "0~1~2~3~4~5~6~7~8~9~10~11~12~13~14~15~16~17~18~19~20~21~22~23~24~100" --get_metric $c_index_type --c_index_group_addition & # --plot_sub_groups  

python3 $(pwd)/code/metrics.py --preds_time_in_days $time --in_data_file $path_data_in --path $path_data_out --outcome "efs" --group_feat_add "0~1~2~3~4~5~6~7~8~9~10~11~12~13~14~15~16~17~18~19~20~21~22~23~24~100" --get_metric $c_index_type --c_index_group_addition & # --plot_sub_groups  
'

: '
## section 2
##### average risk & mean frac of patients plot 
python3 $(pwd)/code/metrics.py --preds_time_in_days $time --in_data_file $path_data_in --path $path_data_out --legend $path_data_legend --get_frac_plot &
python3 $(pwd)/code/metrics.py --preds_time_in_days $time --in_data_file $path_data_in --path $path_data_out --legend $path_data_legend --get_risk_scores &

## comment lines as mentioned in plot_treat_cohort() when running for subsequent times due to computational costs
python3 $(pwd)/code/viz_diff_combo.py --pred_risk_at $time --outcome "risk_pfs" --in_data_file $path_data_in --path $path_data_out --legend $path_data_legend --plot_treat_combo_all_pats &
'

: '
## section 3
##### external test set
python3 $(pwd)/code/metrics.py --preds_time_in_days $time --in_data_file $path_data_in --path $path_data_out --outcome "os" --group_feat_add "0~1~2~9~100" --get_metric $c_index_type --ext_val --c_index_group_addition &

python3 $(pwd)/code/metrics.py --preds_time_in_days $time --in_data_file $path_data_in --path $path_data_out --outcome "efs" --group_feat_add "0~1~2~9~100" --get_metric $c_index_type --ext_val --c_index_group_addition &
'

: '
##### predictions on external test set
python3 $(pwd)/code/viz_diff_combo.py --pred_risk_at $time --outcome "risk_pfs" --in_data_file $path_data_in --path $path_data_out --legend $path_data_legend --plot_treat_combo_indiv_pat --ext_val &
'


