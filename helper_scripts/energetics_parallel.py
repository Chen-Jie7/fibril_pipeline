import pyrosetta.toolbox
import pyrosetta.toolbox.mutants
import pandas as pd
import pyrosetta; pyrosetta.init()
from pyrosetta import *
import os
import yaml
import argparse

init("-ignore_unrecognized_res 1 -ex1 -ex2aro -detect_disulf 0")
scorefxn = create_score_function("ref2015_cart.wts")

# Open the intended fasta file and read the sequence
def design_sequence(pose, sequence,pack_radius=4, name = "", save=False, path=None):
    # Replace the sequence in the pose with the given sequence
    for i, new_residue in enumerate(sequence, start=1):
        if pose.residue(i).name1() != new_residue:
            pyrosetta.toolbox.mutants.mutate_residue(pose, i, new_residue, pack_radius)
    #minimize through fastrelax
    # The task factory accepts all the task operations
    tf = pyrosetta.rosetta.core.pack.task.TaskFactory()

    # These are pretty standard
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.InitializeFromCommandline())
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.IncludeCurrent())
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.NoRepackDisulfides())       
    # Set up a MoveMapFactory
    mmf = pyrosetta.rosetta.core.select.movemap.MoveMapFactory()
    mmf.all_bb(setting=True)
    mmf.all_bondangles(setting=True)
    mmf.all_bondlengths(setting=True)
    mmf.all_chi(setting=True)
    mmf.all_jumps(setting=True)
    mmf.set_cartesian(setting=True)
    display_pose = pyrosetta.rosetta.protocols.fold_from_loops.movers.DisplayPoseLabelsMover()
    display_pose.tasks(tf)
    display_pose.movemap_factory(mmf)
    display_pose.apply(pose)

    fr = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn_in=scorefxn, standard_repeats=1)
    fr.cartesian(True)
    fr.set_task_factory(tf)
    fr.set_movemap_factory(mmf)
    fr.min_type("lbfgs_armijo_nonmonotone") # For non-Cartesian scorefunctions, use "dfpmin_armijo_nonmonotone"
    fr.apply(pose)
    
    if save:
        if path is None:
            #get the current path
            current_path = os.path.dirname(os.path.abspath(__file__))
            pose.dump_pdb(f"{current_path}/{name}.pdb")
        else:
            pose.dump_pdb(f"{path}/{name}.pdb")
    else:   
        return scorefxn(pose)

def calculate_energetics_from_fasta(sequences, temperature,config):
    #create a dataframe to store the energetic core of each sequence
    fasta_sequences = []
    
    # Get configuration values
    proteins = config['mpnn_cleaning']['proteins']
    output_path = config['filepath']['energetic_calculation_path']
    
    for protein in proteins:
        df = pd.DataFrame(columns=["protein", "sequence", "score"])
        protein_name = protein.split("/")[-1]
        fasta_file = f"{config['filepath']['out_folder']}/{protein}/seqs/{protein_name}.fa"
        #open the fasta file and read the sequence
        with open(fasta_file, "r") as f:
            for line in f:
                if not line.startswith(">"):
                    fasta_sequences.append(line.strip().split("/")[0].upper())
    
        #remove duplicates
        fasta_sequences = list(set(fasta_sequences))
        
        for sequence in fasta_sequences:
            #create a pose, get one chain, design the sequence, calculate the energetic core
            pose = Pose()
            pose = pose_from_pdb(f"{config['filepath']['out_folder']}/{protein}/{protein_name}.pdb")
            pose_1 = pose.split_by_chain(1)
            score = design_sequence(pose_1, sequence)
            df = pd.concat([df, pd.DataFrame({"protein": [protein_name], "sequence": [sequence], "score": [score]})], ignore_index=True)
            print(fasta_file.split("/")[-1])
        os.makedirs(output_path, exist_ok=True)
        df.to_csv(f"{output_path}/{protein_name}_{temperature}energetics.csv", index=False)

if __name__ == "__main__":
    args = argparse.ArgumentParser()
    args.add_argument("--config", type=str, required=True)
    args = args.parse_args()
    # Load config file
    with open(args.config, "r") as f:
        config = yaml.safe_load(f)
    temperatures = config['mpnn_cleaning']['energy_calculation']['temperatures']
    for temperature in temperatures:
        calculate_energetics_from_fasta(config, temperature)