import yaml
import os
import json
import numpy as np
import urllib.request
import shutil
import Bio
from helper_scripts.fibrilNormalization import extend_fibril, fast_relax

#remove any terminal warning
import warnings
warnings.filterwarnings("ignore")

def _retrieve_pdb_data_bank(pdb_names, output_path):
    print(pdb_names)
    for pdb_name in pdb_names:
        #fetch from pdb data bank
        try:
            urllib.request.urlretrieve(f'http://files.rcsb.org/download/{pdb_name}.pdb', f'{output_path}/{pdb_name}.pdb')
        except Exception as e:
            print(f"Error fetching {pdb_name}.pdb: {e}")

def _parse_pdb_file(pdb_file_path, output_path):
    for file in os.listdir(pdb_file_path):
        if file.endswith(".pdb"):
            try:    
                #copy it to the output path
                shutil.copy(os.path.join(pdb_file_path, file), os.path.join(output_path, file))
            except Exception as e:
                print(f"Error copying {file}: {e}")

def _parse_pdb_file_single(pdb_file_path, output_path):
    try:
        shutil.copy(pdb_file_path, output_path)
    except Exception as e:
        print(f"Error copying {pdb_file_path}: {e}")


def _clean_up_pdbs(input_folder, output_folder, num_total_layers):
    for file in os.listdir(input_folder):
        if file.endswith(".pdb"):
            pdb_path = os.path.join(input_folder, file)
            output_path = os.path.join(output_folder, file)
            extend_fibril(pdb_path, output_path, num_total_layers=num_total_layers)

def _fast_relax(input_folder, output_folder):
    for file in os.listdir(input_folder):
        #ignore all minimized pdbs
        if file.endswith("min.pdb"):
            continue
        else:
            pdb_path = os.path.join(input_folder, file)
            fast_relax(pdb_path)

def run_mpnn(pdb_name, config):

    #generate chains to design based on user input
    if config["mpnn_model_config"]["positions_to_design"] != "":
        positions_to_design = list(map(str, config["mpnn_model_config"]["positions_to_design"].split()))
        #open pdb file and get how many residues in each chain
        pdb_path = os.path.join(config["filepath"]["out_folder"],pdb_name)
        if not os.path.exists(pdb_path):
            os.makedirs(pdb_path)
            #create a symolink link from the cleaned pdb to the original pdb
            os.symlink(os.path.join(CLEANED_PDB_FOLDER, pdb_name + '.pdb'), pdb_path + f'/{pdb_name}.pdb')
        structure = Bio.PDB.PDBParser().get_structure(pdb_name, pdb_path + f'/{pdb_name}.pdb')
        chain = structure[0]['A']
        chain_length = len(chain)

        fixed_chain = [str(i) for i in range(1, chain_length + 1) if i not in positions_to_design]
        fixed_chain_position = ' '.join(list(fixed_chain))
        tied_position = ' '.join(list(positions_to_design))
        tied_chain_position = ""
        #tied all the chains you want to redesign together since its fibrils
        for i in range(len(config["mpnn_model_config"]["chains_to_design"].split())):
            tied_chain_position += (tied_position + ",")
        tied_chain_position = tied_chain_position[:-1]
    else:
        fixed_chain_position = config["mpnn_model_config"]["fixed_positions"]
        tied_chain_position = config["mpnn_model_config"]["tied_positions"]
    proteinMPNN_path = config["filepath"]["proteinMPNN_path"]
    output_dir = os.path.join(config["filepath"]["out_folder"], pdb_name)
    print(output_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    log_dir = os.path.join(output_dir, "log")
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)
    run_bash_file = f"""#!/bin/bash
#SBATCH --job-name=MPNN_{pdb_name}
#SBATCH --partition=256GBv2
#SBATCH --nodes=1
#SBATCH --time=5:00:00  # Correct time format
#SBATCH --output={output_dir}/log/submit_jobs_%j.out
#SBATCH --error={output_dir}/log/submit_jobs_%j.err

source activate {config["filepath"]["conda_env_path"]}
output_dir={output_dir}

if [ ! -d $output_dir ]
then
    mkdir -p $output_dir
fi

path_for_parsed_chains=$output_dir"/parsed_pdbs.jsonl"
path_for_assigned_chains=$output_dir"/assigned_pdbs.jsonl"
path_for_fixed_positions=$output_dir"/fixed_pdbs.jsonl"
path_for_tied_positions=$output_dir"/tied_pdbs.jsonl"
pdb_path={pdb_path}

chains_to_design="{config["mpnn_model_config"]["chains_to_design"]}"
fixed_positions="{fixed_chain_position}"
tied_positions="{tied_chain_position}"

python {proteinMPNN_path}/helper_scripts/parse_multiple_chains.py --input_path=$pdb_path --output_path=$path_for_parsed_chains

python {proteinMPNN_path}/helper_scripts/assign_fixed_chains.py --input_path=$path_for_parsed_chains --output_path=$path_for_assigned_chains --chain_list "$chains_to_design"

python {proteinMPNN_path}/helper_scripts/make_fixed_positions_dict.py --input_path=$path_for_parsed_chains --output_path=$path_for_fixed_positions --chain_list "$chains_to_design" --position_list "$fixed_positions"

python {proteinMPNN_path}/helper_scripts/make_tied_positions_dict.py --input_path=$path_for_parsed_chains --output_path=$path_for_tied_positions --chain_list "$chains_to_design" --position_list "$tied_positions"

python {proteinMPNN_path}/protein_mpnn_run.py \\
        --jsonl_path $path_for_parsed_chains \\
        --chain_id_jsonl $path_for_assigned_chains \\
        --fixed_positions_jsonl $path_for_fixed_positions \\
        --tied_positions_jsonl $path_for_tied_positions \\
        --out_folder $output_dir \\
        --num_seq_per_target {config["mpnn_model_config"]["num_seq_per_target"]} \\
        --sampling_temp {config["mpnn_model_config"]["sampling_temp"]} \\
        --batch_size {config["mpnn_model_config"]["batch_size"]} \\
        --suppress_print {config["mpnn_model_config"]["suppress_print"]} \\
        --seed {config["mpnn_model_config"]["seed"]} \\
        --save_score {config["mpnn_model_config"]["save_score"]} \\
        --save_probs {config["mpnn_model_config"]["save_probs"]} \\
        --score_only {config["mpnn_model_config"]["score_only"]} \\
        --conditional_probs_only {config["mpnn_model_config"]["conditional_probs_only"]} \\
        --conditional_probs_only_backbone {config["mpnn_model_config"]["conditional_probs_only_backbone"]} \\
        --unconditional_probs_only {config["mpnn_model_config"]["unconditional_probs_only"]} \\
        --backbone_noise {config["mpnn_model_config"]["backbone_noise"]} \\
        --pssm_multi {config["mpnn_model_config"]["pssm_multi"]} \\
        --pssm_threshold {config["mpnn_model_config"]["pssm_threshold"]} \\
        --pssm_log_odds_flag {config["mpnn_model_config"]["pssm_log_odds_flag"]} \\
        --pssm_bias_flag {config["mpnn_model_config"]["pssm_bias_flag"]} \\
        {f'--ca_only {str(config["mpnn_model_config"]["ca_only"])}' if config["mpnn_model_config"]["ca_only"] != False else ""} \
        {f'--use_soluble_model {str(config["mpnn_model_config"]["use_soluble_model"])}' if config["mpnn_model_config"]["use_soluble_model"] != False else ""} \
        {f'--omit_AAs {config["mpnn_model_config"]["omit_AAs"]}' if config["mpnn_model_config"]["omit_AAs"] != "" else ""} \
        {f'--path_to_fasta {config["mpnn_model_config"]["path_to_fasta"]}' if config["mpnn_model_config"]["path_to_fasta"] != "" else ""} \
        {f'--bias_AA_jsonl {config["mpnn_model_config"]["bias_AA_jsonl"]}' if config["mpnn_model_config"]["bias_AA_jsonl"] != "" else ""} \
        {f'--bias_by_res_jsonl {config["mpnn_model_config"]["bias_by_res_jsonl"]}' if config["mpnn_model_config"]["bias_by_res_jsonl"] != "" else ""} \
        {f'--omit_AA_jsonl {config["mpnn_model_config"]["omit_AA_jsonl"]}' if config["mpnn_model_config"]["omit_AA_jsonl"] != "" else ""} \
        {f'--pssm_jsonl {config["mpnn_model_config"]["pssm_jsonl"]}' if config["mpnn_model_config"]["pssm_jsonl"] != "" else ""} \
        {f'--tied_positions_jsonl {config["mpnn_model_config"]["tied_positions_jsonl"]}' if config["mpnn_model_config"]["tied_positions_jsonl"] != "" else ""}

    """
    print(f"Bash file path: {config['filepath']['out_folder']}/{pdb_name}_mpnn.sh")
    with open(f"{config['filepath']['out_folder']}/{pdb_name}_mpnn.sh", "w") as f:
        f.write(run_bash_file)
    os.system(f"sbatch {config['filepath']['out_folder']}/{pdb_name}_mpnn.sh")

if __name__ == "__main__":
    config_path = "mpnn_config.yaml"
    with open(config_path, "r") as f:
        config = yaml.safe_load(f)

    global CLEANED_PDB_FOLDER
    PDB_FOLDER = os.path.join(config["filepath"]["out_folder"], "original_pdbs")
    CLEANED_PDB_FOLDER = os.path.join(config["filepath"]["out_folder"], "cleaned_pdbs")
    
    if config["Phases"]["Phase1"]:
        #make output folders and subfolders
        os.makedirs(config["filepath"]["out_folder"], exist_ok=True)
        os.makedirs(config["filepath"]["out_folder"] + "/original_pdbs", exist_ok=True)
        os.makedirs(config["filepath"]["out_folder"] + "/cleaned_pdbs", exist_ok=True)

        #Get all the pdb fiels onto the sampe place
        if config["filepath"]["data_bank_pdb_name"] != "":
            _retrieve_pdb_data_bank(list(config["filepath"]["data_bank_pdb_name"].split(',')), PDB_FOLDER)
        if config["filepath"]["pdb_folder_path"] != "":
            _parse_pdb_file(config["filepath"]["pdb_folder_path"], PDB_FOLDER)
        if config["filepath"]["pdb_path"] != "":
            _parse_pdb_file_single(config["filepath"]["pdb_path"], PDB_FOLDER)

        print("PDB obtained")
        #Clean up the pdbs
        _clean_up_pdbs(PDB_FOLDER, CLEANED_PDB_FOLDER, config["pdb_clean_up"]["layers"])
        print("PDB cleaned")
        #Fast relax the cleaned pdbs
        if config["pdb_clean_up"]["fast_relax"]:
            _fast_relax(CLEANED_PDB_FOLDER, config["pdb_clean_up"]["layers"])

        print("Clean Up Complete")

    if config["Phases"]["Phase2"]:
        for pdb_name in os.listdir(CLEANED_PDB_FOLDER):
            if pdb_name.endswith("min.pdb"):
                run_mpnn(pdb_name.split('.')[0], config)
        print("MPNN Submitted")


