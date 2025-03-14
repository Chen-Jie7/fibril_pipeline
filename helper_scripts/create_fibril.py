import pyrosetta.toolbox
import pyrosetta.toolbox.mutants
import pandas as pd
import pyrosetta; pyrosetta.init()
import argparse
from pyrosetta import *
import os
init("-ignore_unrecognized_res 1 -ex1 -ex2aro -detect_disulf 0")
scorefxn = create_score_function("ref2015_cart.wts")
# Open the intended fasta file and read the sequence

def design_sequence(pdb_path, name, sequence, num, save=True, path=None):
    pack_radius=4
    pose = Pose()
    pose = pose_from_pdb(pdb_path)
    # Get the number of chains in the pose
    num_chains = pose.num_chains()
    seq_length = len(sequence)
    
    for chain_num in range(1, num_chains + 1):
        # Get the start and end residue numbers for the current chain
        chain_begin = pose.chain_begin(chain_num)
        chain_end = pose.chain_end(chain_num)
        chain_length = chain_end - chain_begin + 1
        
        # Ensure the chain length matches the sequence length
        if chain_length != seq_length:
            print(f"Warning: Chain {chain_num} length ({chain_length}) does not match the sequence length ({seq_length}).")
            # You can choose to skip this chain, adjust the sequence, or handle as needed
            continue
        
        # Mutate the residues in the current chain
        for i, new_residue in zip(range(chain_begin, chain_end + 1), sequence):
            if pose.residue(i).name1() != new_residue:
                pyrosetta.toolbox.mutants.mutate_residue(pose, i, new_residue, pack_radius)
    
        #minimize through fastrelax
        # The task factory accepts all the task operations
        tf = pyrosetta.rosetta.core.pack.task.TaskFactory()

        # These are pretty standard
        tf = pyrosetta.rosetta.core.pack.task.TaskFactory()
        tf.push_back(pyrosetta.rosetta.core.pack.task.operation.InitializeFromCommandline())
        tf.push_back(pyrosetta.rosetta.core.pack.task.operation.IncludeCurrent())
        tf.push_back(pyrosetta.rosetta.core.pack.task.operation.NoRepackDisulfides())
        tf.push_back(pyrosetta.rosetta.core.pack.task.operation.RestrictToRepacking())  # Add this line

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
        pdb_naming = f"{name}_{num}" if num is not None else name

        if path is None:
            #get the current path
            current_path = os.path.dirname(os.path.abspath(__file__))
            pose.dump_pdb(f"{current_path}/{pdb_naming}.pdb")
        else:
            pose.dump_pdb(f"{path}/{pdb_naming}.pdb")
    else:   
        return scorefxn(pose)

def design_multiple_sequence(pose, sequence, starting_chain, ending_chain, pack_radius=2, name="", save=False, path=None):
    # Get the number of chains in the pose
    num_chains = pose.num_chains()
    seq_length = len(sequence)
    
    for chain_num in range(starting_chain, ending_chain + 1):
        # Get the start and end residue numbers for the current chain
        residue_begin = (chain_num-1)*seq_length + 1
        residue_end = residue_begin + seq_length*(ending_chain-starting_chain)
        # Mutate the residues in the current chain
        for i, new_residue in zip(range(residue_begin, residue_end + 1), sequence):
            if pose.residue(i).name1() != new_residue:
                pyrosetta.toolbox.mutants.mutate_residue(pose, i, new_residue, pack_radius)
    
        #minimize through fastrelax
        # The task factory accepts all the task operations
        tf = pyrosetta.rosetta.core.pack.task.TaskFactory()

        # These are pretty standard
        tf = pyrosetta.rosetta.core.pack.task.TaskFactory()
        tf.push_back(pyrosetta.rosetta.core.pack.task.operation.InitializeFromCommandline())
        tf.push_back(pyrosetta.rosetta.core.pack.task.operation.IncludeCurrent())
        tf.push_back(pyrosetta.rosetta.core.pack.task.operation.NoRepackDisulfides())
        tf.push_back(pyrosetta.rosetta.core.pack.task.operation.RestrictToRepacking())  # Add this line

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
        return pose

#by default, the function will design the top sequences for a single chain pdb file
def design_top_sequence(pdb_name, sequence, num, save=True, path=None):
    #make sure math exist
    os.makedirs(path, exist_ok=True)
    pose = Pose()
    pose = pose_from_pdb(f"/work/CAND/shared/Chens/deNovoFibrils/PDBNormalization_PDB/output/{pdb_name}/{pdb_name}.pdb")
    if single:
        pose_1 = pose.split_by_chain(1)
    else:
        pose_1 = pose.clone()
    
    design_single_sequence(pose_1, sequence, name = f"{pdb_name}_top_{num}", save=save, path=path)
def design_bot_sequence(pdb_name, sequence, num, single = True, save=True, path=None):
    os.makedirs(path, exist_ok=True)
    
    pose = Pose()
    pose = pose_from_pdb(f"/work/CAND/shared/Chens/deNovoFibrils/PDBNormalization_PDB/output/{pdb_name}/{pdb_name}.pdb")
    if single:
        pose_1 = pose.split_by_chain(1)
    else:
        pose_1 = pose.clone()
    
    design_single_sequence(pose_1, sequence, name = f"{pdb_name}_bot_{num}", save=save, path=path)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb_path", type=str)
    parser.add_argument("--sequence", type=str)
    parser.add_argument("--name", type=str)
    parser.add_argument("--num", type=int)
    parser.add_argument("--save", type=bool)
    parser.add_argument("--path", type=str)
    args = parser.parse_args()  
    
    design_sequence(args.pdb_path, args.name, args.sequence, args.num, save=args.save, path=args.path)
    