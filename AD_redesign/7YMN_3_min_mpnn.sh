#!/bin/bash
#SBATCH --job-name=MPNN_7YMN_3_min
#SBATCH --partition=256GBv1
#SBATCH --nodes=1
#SBATCH --time=12:20:00  # Correct time format
#SBATCH --output=/work/CAND/shared/Chens/AD_fibrils/AD_redesign/7YMN_3_min/log/submit_jobs_%j.out
#SBATCH --error=/work/CAND/shared/Chens/AD_fibrils/AD_redesign/7YMN_3_min/log/submit_jobs_%j.err

source activate /endosome/work/CAND/shared/Chens/miniforge3/envs/upside2
output_dir=/work/CAND/shared/Chens/AD_fibrils/AD_redesign/7YMN_3_min

if [ ! -d $output_dir ]
then
    mkdir -p $output_dir
fi

path_for_parsed_chains=$output_dir"/parsed_pdbs.jsonl"
path_for_assigned_chains=$output_dir"/assigned_pdbs.jsonl"
path_for_fixed_positions=$output_dir"/fixed_pdbs.jsonl"
path_for_tied_positions=$output_dir"/tied_pdbs.jsonl"
pdb_path=/work/CAND/shared/Chens/AD_fibrils/AD_redesign/7YMN_3_min

chains_to_design="A B"
fixed_positions="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76"
tied_positions="10 15 24 10,10 15 24 10"

python /work/CAND/shared/Chens/AD_fibrils/ProteinMPNN/helper_scripts/parse_multiple_chains.py --input_path=$pdb_path --output_path=$path_for_parsed_chains

python /work/CAND/shared/Chens/AD_fibrils/ProteinMPNN/helper_scripts/assign_fixed_chains.py --input_path=$path_for_parsed_chains --output_path=$path_for_assigned_chains --chain_list "$chains_to_design"

python /work/CAND/shared/Chens/AD_fibrils/ProteinMPNN/helper_scripts/make_fixed_positions_dict.py --input_path=$path_for_parsed_chains --output_path=$path_for_fixed_positions --chain_list "$chains_to_design" --position_list "$fixed_positions"

python /work/CAND/shared/Chens/AD_fibrils/ProteinMPNN/helper_scripts/make_tied_positions_dict.py --input_path=$path_for_parsed_chains --output_path=$path_for_tied_positions --chain_list "$chains_to_design" --position_list "$tied_positions"

python /work/CAND/shared/Chens/AD_fibrils/ProteinMPNN/protein_mpnn_run.py \
        --jsonl_path $path_for_parsed_chains \
        --chain_id_jsonl $path_for_assigned_chains \
        --fixed_positions_jsonl $path_for_fixed_positions \
        --tied_positions_jsonl $path_for_tied_positions \
        --out_folder $output_dir \
        --num_seq_per_target 50 \
        --sampling_temp 0.4 \
        --batch_size 1 \
        --suppress_print 0 \
        --seed 7 \
        --save_score 0 \
        --save_probs 0 \
        --score_only 0 \
        --conditional_probs_only 0 \
        --conditional_probs_only_backbone 0 \
        --unconditional_probs_only 0 \
        --backbone_noise 0.0 \
        --pssm_multi 0.0 \
        --pssm_threshold 0.0 \
        --pssm_log_odds_flag 0 \
        --pssm_bias_flag 0 \
                                                                                

    