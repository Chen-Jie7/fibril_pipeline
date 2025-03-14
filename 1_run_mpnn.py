import yaml
import os
import json
import numpy as np
import urllib.request
import shutil
import Bio
from helper_scripts.fibrilNormalization import extend_fibril, fast_relax
from helper_scripts.mpnn_cleaning import combine_fasta_files
#remove any terminal warning
import warnings
warnings.filterwarnings("ignore")
import pandas as pd

def _retrieve_pdb_data_bank(pdb_names, output_path):
    print(pdb_names)
    successful_pdbs = []
    for pdb_name in pdb_names:
        #fetch from pdb data bank
        try:
            urllib.request.urlretrieve(f'http://files.rcsb.org/download/{pdb_name}.pdb', f'{output_path}/{pdb_name}.pdb')
            successful_pdbs.append(pdb_name)
        except Exception as e:
            print(f"Error fetching {pdb_name}.pdb: {e}")
    return successful_pdbs

def _parse_pdb_file(pdb_file_path, output_path):
    successful_pdbs = []
    for file in os.listdir(pdb_file_path):
        if file.endswith(".pdb"):
            try:    
                #copy it to the output path
                shutil.copy(os.path.join(pdb_file_path, file), os.path.join(output_path, file))
                # Extract just the base filename without extension
                pdb_name = os.path.splitext(file)[0]
                successful_pdbs.append(pdb_name)
            except Exception as e:
                print(f"Error copying {file}: {e}")
    return successful_pdbs

def _parse_pdb_file_single(pdb_file_path, output_path):
    try:
        # Extract just the base filename without path and extension
        filename = os.path.basename(pdb_file_path)
        pdb_name = os.path.splitext(filename)[0]
        
        shutil.copy(pdb_file_path, output_path)
        return [pdb_name]
    except Exception as e:
        print(f"Error copying {pdb_file_path}: {e}")
        return []

def _clean_up_pdbs(input_folder, output_folder, num_total_layers):
    for file in os.listdir(input_folder):
        if file.endswith(".pdb"):
            pdb_path = os.path.join(input_folder, file)
            output_path = os.path.join(output_folder, file)
            extend_fibril(pdb_path, output_path, num_total_layers=num_total_layers)

def _fast_relax(input_folder, output_folder):
    for file in os.listdir(input_folder):
        #ignore all minimized pdbs
        if file.split('.')[0] not in RUNNING_PDBS or file.endswith("min.pdb"):
            continue
        else:
            pdb_path = os.path.join(input_folder, file)
            fast_relax(pdb_path)

def run_mpnn(pdb_name, config, temp ):
    #generate chains to design based on user input
    if config["mpnn_model_config"]["positions_to_design"] != "":
        positions_to_design = list(map(int, config["mpnn_model_config"]["positions_to_design"].split()))
        #open pdb file and get how many residues in each chain
        pdb_path = os.path.join(config["filepath"]["out_folder"],pdb_name)
        if not os.path.exists(pdb_path):
            os.makedirs(pdb_path)
            #create a symolink link from the cleaned pdb to the original pdb
            os.symlink(os.path.join(CLEANED_PDB_FOLDER, pdb_name + '.pdb'), pdb_path + f'/{pdb_name}.pdb')
        structure = Bio.PDB.PDBParser().get_structure(pdb_name, pdb_path + f'/{pdb_name}.pdb')
        chain = structure[0]['A']
        chain_length = len(chain)

        # Generate fixed positions by excluding positions_to_design
        fixed_chain = [str(i) for i in range(1, chain_length + 1) if i not in positions_to_design]
        fixed_chain = ' '.join(fixed_chain)
        
        #repeat fix for number of chains to design
        fixed_chain_position = ""
        for i in range(len(config["mpnn_model_config"]["chains_to_design"].split())):
            fixed_chain_position += (fixed_chain + ", ")
        fixed_chain_position = fixed_chain_position[:-2]
        
        tied_position = ' '.join(map(str, positions_to_design))
        tied_chain_position = ""
        #tied all the chains you want to redesign together since its fibrils
        for i in range(len(config["mpnn_model_config"]["chains_to_design"].split())):
            tied_chain_position += (tied_position + ", ")
        tied_chain_position = tied_chain_position[:-2]
    else:
        fixed_chain_position = config["mpnn_model_config"]["fixed_positions"]
        tied_chain_position = config["mpnn_model_config"]["tied_positions"]
    
    proteinMPNN_path = config["filepath"]["proteinMPNN_path"]
    output_dir = os.path.join(config["filepath"]["out_folder"], pdb_name)
    if not os.path.exists(os.path.join(output_dir, temp)):
        os.makedirs(os.path.join(output_dir, temp))
    log_dir = os.path.join(output_dir, "logs")
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)
    run_bash_file = f"""#!/bin/bash
#SBATCH --job-name={temp}_MPNN_{pdb_name}
#SBATCH --partition=256GBv2
#SBATCH --nodes=1
#SBATCH --time=5:00:00  # Correct time format
#SBATCH --output={output_dir}/logs/submit_jobs_%j.out
#SBATCH --error={output_dir}/logs/submit_jobs_%j.err

source activate {config["filepath"]["conda_env_path"]}
output_dir={output_dir}/{temp}

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
        --sampling_temp {temp} \\
        --batch_size {config["mpnn_model_config"]["batch_size"]} \\
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
    print(f"Bash file path: {config['filepath']['out_folder']}/{pdb_name}/{pdb_name}_mpnn_{temp}.sh")
    with open(f"{config['filepath']['out_folder']}/{pdb_name}/{pdb_name}_mpnn_{temp}.sh", "w") as f:
        f.write(run_bash_file)
    os.system(f"sbatch {config['filepath']['out_folder']}/{pdb_name}/{pdb_name}_mpnn_{temp}.sh")

#create a bash file to script the pyrosetta energetic calculation

def create_energetics_bash_file(config):
    output_folder = config['filepath']['out_folder']
    log_folder = os.path.join(output_folder, "log")
    os.makedirs(log_folder, exist_ok=True)
    with open(f"{output_folder}/energetics_launch.sh", "w") as f:
        f.write(f"""#!/bin/bash
#SBATCH --partition=super
#SBATCH -c 4
#SBATCH --time=12:00:00
#SBATCH --output={log_folder}/energetics_launch_%j.out
#SBATCH --error={log_folder}/energetics_launch_%j.err
                
source {config['filepath']['conda_env_path']}
conda activate {config['filepath']['conda_env_path']}

python {os.path.join(CURRENT_FOLDER, "helper_scripts/energetics_calculation.py")} --config mpnn_config.yaml
""")
        
#Two files, one are the command line arguments, one is the bash file to launch the array job for running energetics
def optimize_energetics(config):
    output_folder = config['filepath']['out_folder']
    proteins = config['mpnn_cleaning']['proteins']
    for protein in proteins:
        log_folder = os.path.join(f"{output_folder}/", f"{protein}/logs")
        os.makedirs(log_folder, exist_ok=True)
        prev_line = None
        # Read all sequences and headers first
        sequences = []
        headers = []
        with open(f"{output_folder}/{protein}/combined_fasta.fa", "r") as fasta_file:
            for line in fasta_file:
                if line.startswith(">"):
                    parts = line[1:].strip().split(", ")
                    temp = parts[0].split("=")[1]
                    sample = parts[1].split("=")[1]
                    headers.append(f"T_{temp}_S_{sample}")
                else:
                    sequences.append(line.strip())

        # # Write commands in batches of 200
        batch_size = 20
        num_batches = (len(sequences)// batch_size) + 1
        os.makedirs(f"{output_folder}/{protein}/energetics_seqs", exist_ok=True)
        if os.path.exists(f"{output_folder}/{protein}/energetics_command_line.sh"):
            os.remove(f"{output_folder}/{protein}/energetics_command_line.sh")


        with open(f"{output_folder}/{protein}/energetics_command_line.sh", "a") as write_command_file:
            for batch in range(num_batches):
                with open(f"{output_folder}/{protein}/energetics_seqs/{protein}_seq_bash_{batch+1}.fa", "w") as write_file:
                    batch_start = batch * batch_size
                    batch_end = min((batch + 1) * batch_size, len(sequences))
                    #write the batch amount of sequences to a new file
                    for i in range(batch_start, batch_end):
                        write_file.write(f"{headers[i]}-{sequences[i]}\n")
                    #check if energetics_command_line.sh exists if it does delete it and create a new one
                    write_command_file.write(f"python {os.path.join(CURRENT_FOLDER, 'helper_scripts/energetics_calculation.py')} --sequence_file {output_folder}/{protein}/energetics_seqs/{protein}_seq_bash_{batch+1}.fa --pdb_path {output_folder}/{protein}/{protein}.pdb\n")
        #create a bash file to launch the array job for running energetics
        with open(f"{output_folder}/{protein}/energetics_launch.sh", "w") as f:
            f.write(f"""#!/bin/bash
#SBATCH --partition=super
#SBATCH --job-name=array_job
#SBATCH --output={log_folder}/energetics_%A_%a.out
#SBATCH --error={log_folder}/energetics_%A_%a.err
#SBATCH --array=0-{num_batches-1}   # This will create jobs numbered 0 to num_batches-1
#SBATCH --time=10:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=2

# Get the command from the command file based on array task ID
CMD=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" {output_folder}/{protein}/energetics_command_line.sh)
# Execute the command
$CMD
""")

#Two files, one are the command line arguments, one is the bash file to launch the array job for running energetics
def design_top_fibrils(protein, config):
    output_folder = os.path.join(config['fibril_design']['output_path'], protein)
    log_folder = os.path.join(f"{output_folder}/", f"logs")
    os.makedirs(log_folder, exist_ok=True)
    #find energetics csv file
    total_energetics = pd.DataFrame()
    energetics_csv = os.path.join(output_folder, "energetics")
    for file in os.listdir(energetics_csv):
        if file.endswith(".csv"):
            df = pd.read_csv(os.path.join(energetics_csv, file))
            total_energetics = pd.concat([total_energetics, df], ignore_index=True)
    total_energetics.sort_values(by="score", inplace=True, ignore_index=True)
    top_N_fibrils = total_energetics.head(config['fibril_design']['N_fibrils'])
    num_batches = len(top_N_fibrils)
    if os.path.exists(f"{output_folder}/top_N_command_line.sh"):
        os.remove(f"{output_folder}/top_N_command_line.sh")
    os.makedirs(os.path.join(output_folder, 'top_N_fibrils'), exist_ok=True)
    with open(f"{output_folder}/top_N_command_line.sh", "a") as write_command_file:
        print(top_N_fibrils)
        for index, row in top_N_fibrils.iterrows():
            write_command_file.write(f"python {CURRENT_FOLDER}/helper_scripts/create_fibril.py --pdb_path {os.path.join(output_folder, protein + '.pdb')} --sequence {row['sequence']} --name {row['protein']} --num {index+1} --save True --path {os.path.join(output_folder, 'top_N_fibrils')}\n")
    with open(f"{output_folder}/top_N_launch.sh", "w") as f:
        f.write(f"""#!/bin/bash
#SBATCH --partition=super
#SBATCH --job-name=array_job
#SBATCH --output={log_folder}/energetics_%A_%a.out
#SBATCH --error={log_folder}/energetics_%A_%a.err
#SBATCH --array=0-{num_batches-1}
#SBATCH --time=2:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=2

# Load Python environment if needed
source {config['filepath']['conda_env_path']}

# Get the command from the command file based on array task ID
CMD=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" {output_folder}/top_N_command_line.sh)

# Execute the command with python explicitly
python $CMD
""")
    print(f"Top N Launch File: {output_folder}/top_N_launch.sh")

if __name__ == "__main__":
    config_path = "mpnn_config.yaml"
    with open(config_path, "r") as f:
        config = yaml.safe_load(f)

    global CLEANED_PDB_FOLDER, PDB_FOLDER, CURRENT_FOLDER, RUNNING_PDBS
    PDB_FOLDER = os.path.join(config["filepath"]["out_folder"], "original_pdbs")
    CLEANED_PDB_FOLDER = os.path.join(config["filepath"]["out_folder"], "cleaned_pdbs")
    CURRENT_FOLDER = os.path.dirname(os.path.abspath(__file__))
    RUNNING_PDBS = []

        #Get all the pdb fiels onto the sampe place
    if config["filepath"]["data_bank_pdb_name"] != "":
        RUNNING_PDBS.extend(_retrieve_pdb_data_bank(list(config["filepath"]["data_bank_pdb_name"].split(',')), PDB_FOLDER))
    if config["filepath"]["pdb_folder_path"] != "":
        RUNNING_PDBS.extend(_parse_pdb_file(config["filepath"]["pdb_folder_path"], PDB_FOLDER))
    if config["filepath"]["pdb_path"] != "":
        RUNNING_PDBS.extend(_parse_pdb_file_single(config["filepath"]["pdb_path"], PDB_FOLDER))

    if config["Phases"]["Phase1"]:
        #make output folders and subfolders
        os.makedirs(config["filepath"]["out_folder"], exist_ok=True)
        os.makedirs(config["filepath"]["out_folder"] + "/original_pdbs", exist_ok=True)
        os.makedirs(config["filepath"]["out_folder"] + "/cleaned_pdbs", exist_ok=True)

        print("PDB obtained")
        #Clean up the pdbs
        _clean_up_pdbs(PDB_FOLDER, CLEANED_PDB_FOLDER, config["pdb_clean_up"]["layers"])
        print("PDB cleaned")
        #Fast relax the cleaned pdbs
        if config["pdb_clean_up"]["fast_relax"]:
            _fast_relax(CLEANED_PDB_FOLDER, config["pdb_clean_up"]["layers"])

        print("Clean Up Complete")

    if config["Phases"]["Phase2"]:
        #collect pdbs running rn

        for pdb_name in os.listdir(CLEANED_PDB_FOLDER):
            if pdb_name.split('.')[0].split('_')[0] in RUNNING_PDBS and "min" in pdb_name:
                for temp in config["mpnn_model_config"]["sampling_temp"]:
                    run_mpnn(pdb_name.split('.')[0], config, temp)
                print("MPNN Submitted")

    if config["Phases"]["Phase3"]:
        for protein in config['mpnn_cleaning']['proteins']:
            #cleans up the fasta files by merging them into one file and removing duplicates
            #outputs a final fasta file combining all the temperatures
            combine_fasta_files(f"{config['filepath']['out_folder']}/{protein}", config['mpnn_cleaning']['temperatures'])
            optimize_energetics(config)
            #calcualte the energetic core of each sequence by threading the sequence to the pdb
            #and then calculating the energetics of a fibril through pyrosetta
    if config["Phases"]["Phase4"]:
        for protein in config['fibril_design']['proteins']:
            design_top_fibrils(protein, config)