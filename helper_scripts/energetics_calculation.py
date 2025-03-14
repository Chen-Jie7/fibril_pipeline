import pyrosetta
import pandas as pd
import os
import yaml
import argparse
import pyrosetta
from pyrosetta import rosetta

pyrosetta.init("-ignore_unrecognized_res 1 -ex1 -ex2aro -detect_disulf 0 -load_PDB_components false")
scorefxn = pyrosetta.create_score_function("ref2015_cart.wts")

def design_sequence(pose, sequence, pack_radius=4, name="", save=False, path=None):
    # Replace the sequence in the pose with the given sequence
    # Loop over the target sequence (using 1-indexing for the pose)
    for i, new_residue in enumerate(sequence, start=1):
        # Only mutate if the current residue doesn't match the target.
        if pose.residue(i).name1() != new_residue:
            # Mutate the residue using the MutateResidue mover.
            mutator = rosetta.protocols.simple_moves.MutateResidue(i, new_residue)
            mutator.apply(pose)
            
            # Set up a task factory for repacking the neighborhood.
            packer = rosetta.protocols.minimization_packing.PackRotamersMover()
            tf = rosetta.core.pack.task.TaskFactory()
            # First, restrict the entire pose to repacking (i.e. no design).
            tf.push_back(rosetta.core.pack.task.operation.RestrictToRepacking())
            # Now, for the mutated residue’s neighborhood,
            # explicitly allow repacking using a RestrictToRepackingRLT operation.
            nbr_selector = rosetta.core.select.residue_selector.NeighborhoodResidueSelector(
                rosetta.core.select.residue_selector.ResidueIndexSelector(str(i)),
                pack_radius
            )
            op = rosetta.core.pack.task.operation.RestrictToRepackingRLT()
            tf.push_back(rosetta.core.pack.task.operation.OperateOnResidueSubset(op, nbr_selector))
            
            packer.task_factory(tf)
            packer.apply(pose)

    # Set up a task factory for the FastRelax protocol.
    tf_relax = rosetta.core.pack.task.TaskFactory()
    tf_relax.push_back(rosetta.core.pack.task.operation.InitializeFromCommandline())
    tf_relax.push_back(rosetta.core.pack.task.operation.IncludeCurrent())
    tf_relax.push_back(rosetta.core.pack.task.operation.NoRepackDisulfides())
    tf_relax.push_back(rosetta.core.pack.task.operation.RestrictToRepacking())

    # Create a MoveMapFactory to allow movements.
    mmf = rosetta.core.select.movemap.MoveMapFactory()
    mmf.all_bb(setting=True)
    mmf.all_bondangles(setting=True)
    mmf.all_bondlengths(setting=True)
    mmf.all_chi(setting=True)
    mmf.all_jumps(setting=True)
    mmf.set_cartesian(setting=True)

    # (Optional) Display pose labels for debugging.
    display_pose = rosetta.protocols.fold_from_loops.movers.DisplayPoseLabelsMover()
    display_pose.tasks(tf_relax)
    display_pose.movemap_factory(mmf)
    display_pose.apply(pose)

    # Set up and run the FastRelax protocol.
    fr = rosetta.protocols.relax.FastRelax(scorefxn_in=scorefxn, standard_repeats=1)
    fr.cartesian(True)
    fr.set_task_factory(tf_relax)
    fr.set_movemap_factory(mmf)
    fr.min_type("lbfgs_armijo_nonmonotone")
    fr.apply(pose)
    
    if save:
        if path is None:
            current_path = os.path.dirname(os.path.abspath(__file__))
            pose.dump_pdb(f"{current_path}/{name}.pdb")
        else:
            pose.dump_pdb(f"{path}/{name}.pdb")
    score = scorefxn(pose)
    print(score)
    return score

def design_top_sequences(pdb_path, sequences, single=True, save=True, path=None):
    pose = pyrosetta.Pose()
    pdb_name = pdb_path.split("/")[-1]
    pose = pyrosetta.pose_from_pdb(f"{pdb_path}/{pdb_name}.pdb")
    if single:
        pose_1 = pose.split_by_chain(1)
    else:
        pose_1 = pose.clone()
    
    for num, sequence in enumerate(sequences):
        pose_2 = pose_1.clone()
        design_sequence(pose_2, sequence, name=f"{pdb_name}_top_{num+1}", save=save, path=path)
    return pose_1

#functin to map single letter amino acid to three letter amino acid
def map_aa(aa):
    dict = {
        "A": "ALA",
        "C": "CYS",
        "D": "ASP",
        "E": "GLU",
        "F": "PHE",
        "G": "GLY",
        "H": "HIS",
        "I": "ILE",
        "K": "LYS",
        "L": "LEU",
        "M": "MET",
        "N": "ASN",
        "P": "PRO",
        "Q": "GLN",
        "R": "ARG",
        "S": "SER",
        "T": "THR",
        "V": "VAL",
        "W": "TRP",
        "Y": "TYR"
    }
    return dict[aa]
        
def calculate_energetics_from_fasta(sequence_file, pdb_path):

    out = pdb_path.split("/")
    out = "/".join(out[:-1])
    output_path =os.path.join(out, "energetics")
  

    os.makedirs(output_path, exist_ok=True)
    df = pd.DataFrame(columns=["protein", "sequence", "score"])
    #read the fasta file where each line is a sequence and separated by dashes is the protein name
    with open(sequence_file, "r") as f:
        for line in f:
            protein_name, sequence = line.strip().split("-")            
            pose = pyrosetta.Pose()
            pose = pyrosetta.pose_from_pdb(pdb_path)
            pose_1 = pose.split_by_chain(1)
            three_letter_sequence = [map_aa(aa) for aa in sequence]
            score = design_sequence(pose_1, three_letter_sequence, save=True, name=protein_name, path=output_path)
            df = pd.concat([df, pd.DataFrame({
                "protein": [protein_name], 
                "sequence": [sequence], 
                "score": [score]
            })], ignore_index=True)
    
    print(output_path)
    print("OUTPUT PATHHHHHHH")
    print(f"RESULTS SAVED TO {output_path}/{sequence_file.split('/')[-1].split('.')[0]}_energetics.csv")
    
    df.to_csv(f"{output_path}/{sequence_file.split('/')[-1].split('.')[0]}_energetics.csv", index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--sequence_file", type=str, required=True)
    parser.add_argument("--pdb_path", type=str, required=True)
    args = parser.parse_args()

    calculate_energetics_from_fasta(args.sequence_file, args.pdb_path)