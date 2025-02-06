from Bio.PDB import PDBParser, PDBIO, Structure, Model, Chain, MMCIFIO
import numpy as np
from copy import deepcopy
import string
from scipy.spatial.distance import pdist, squareform
from icecream import ic
from itertools import product

import pyrosetta
from pyrosetta import *
from pyrosetta.rosetta.protocols.minimization_packing import MinMover
from pyrosetta.rosetta.core.kinematics import MoveMap
import argparse
import sys




def _initialize_pyrosetta(memory_limit=None):
    """Initialize PyRosetta with optional memory limit"""
    init_options = []
    if memory_limit:
        init_options.append(f"-mem_gb {memory_limit}")
    pyrosetta.init(' '.join(init_options))

def _minimize_backbone(pose, scorefxn=None, max_iter=1000, minimize_sidechains=True):
    """
    Perform backbone minimization using PyRosetta
    
    Args:
        pose: PyRosetta Pose object
        scorefxn: Score function to use (default: ref2015)
        max_iter: Maximum number of iterations for minimization
    
    Returns:
        minimized_pose: Minimized pose object
        final_score: Final energy score
    """
    # Create a copy of the input pose
    minimized_pose = Pose()
    minimized_pose.assign(pose)
    
    # Set up score function if not provided
    if scorefxn is None:
        scorefxn = pyrosetta.get_fa_scorefxn()
    
    # Create MoveMap for backbone minimization
    movemap = MoveMap()
    movemap.set_bb(True)  # Enable backbone movement
    movemap.set_chi(minimize_sidechains)  # Control side-chain movement
    
    
    # Create MinMover
    min_mover = MinMover(movemap, scorefxn, 'lbfgs_armijo_nonmonotone', 0.001, True)
    min_mover.max_iter(max_iter)
    
    # Perform minimization
    min_mover.apply(minimized_pose)
    
    # Calculate final score
    final_score = scorefxn(minimized_pose)
    
    return minimized_pose, final_score

def fast_relax(input_file, max_iter=1000, sidechains=True):

    _initialize_pyrosetta()
    pose = pyrosetta.pose_from_pdb(input_file)
    # Perform minimization
    print(f"Starting minimization (max {max_iter} iterations)...")
    initial_score = pyrosetta.get_fa_scorefxn()(pose)
    print(f"Initial score: {initial_score:.2f}")
    
    minimized_pose, final_score = _minimize_backbone(pose, max_iter=max_iter, minimize_sidechains=sidechains)
    
    print(f"Final score: {final_score:.2f}")
    print(f"Score improvement: {initial_score - final_score:.2f}")
    
    # Save the minimized structure
    output_file = f"{input_file.split('.')[0]}_min.pdb"
    minimized_pose.dump_pdb(output_file)
    print(f"Minimized structure saved to {output_file}")
        

# Function to extract atom coordinates indexed by residue ID and atom name
def get_atom_coords(chain):
    atom_coords = {}
    for residue in chain.get_residues():
        res_id = residue.get_id()
        for atom in residue.get_atoms():
            atom_name = atom.get_name()
            key = (res_id, atom_name)
            atom_coords[key] = atom.get_coord()
    return atom_coords
    


def remove_non_residues(model):
    standard_residues = [
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS",
        "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
        "TYR", "VAL"
    ]
    chains = list(model.get_chains())
    for chain in chains:
        residues = list(chain.get_residues())
        for residue in residues:
            if residue.get_resname() not in standard_residues:
                chain.detach_child(residue.id)
    return model

def normalize_fibrils(structure, chain_id_list):    
    model = structure[0]
    fibril_chains = list(model.get_chains())
    
    chain_coords = []
    for chain in fibril_chains:
        avg_coords = np.zeros(3)
        for residue in chain:
            atom = residue['CA'] if 'CA' in residue else next(residue.get_atoms())
            avg_coords += atom.get_coord()
        avg_coords /= len(chain)
        chain_coords.append(avg_coords)
    
    distances = pdist(np.array(chain_coords))
    distances = squareform(distances)
    #cluster the distances into two groups based on the median
    median_distance = np.median(distances)
    ic(median_distance)
    # # Select fibril chains based on distance from the first chain
    fibrils = [chain.id for num, chain in enumerate(fibril_chains) if distances[0][num] < median_distance]
    ic(fibrils)

    # Remove non-fibril chains from the model
    for chain in fibril_chains:
        if chain.id not in fibrils:
            model.detach_child(chain.id)
        else:
            chain.id = f"Z{chain.id}"
 
    # Recalculate chain coordinates for remaining chains
    chain_coords = {chain.id: np.mean([atom.get_coord() for atom in chain.get_atoms()], axis=0) for chain in model.get_chains()}

    # Sort chains based on their Z-coordinate
    chain_coords = dict(sorted(chain_coords.items(), key=lambda item: item[1][2]))
    chain_ids_sorted = list(chain_coords.keys())

   
    chain_id_mapping = {}
    for i, old_chain_id in enumerate(chain_ids_sorted):
        chain_id_mapping[old_chain_id] = chain_id_list[i]
    
    # Update chain IDs and residue numbers
    for chain in model.get_chains():
        if chain.id in chain_id_mapping:
            chain.id = chain_id_mapping[chain.id]
        for i, residue in enumerate(chain.get_residues(), start=1):
            residue.id = (' ', i, ' ')  
       
    return structure
# Function to modify the label_asym_id instead of auth_asym_id
def rename_label_asym_id(structure, chain_id_list):
    model = structure[0]
    for chain, new_id in zip(model.get_chains(), chain_id_list):
        chain.id = new_id  # This will modify label_asym_id when saved to CIF

def extend_fibril(input_file, output_file, num_total_layers=6, specific_range=None):

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('fibril', input_file)    
    
    model = structure[0]
    if specific_range is not None:
        start_residue, end_residue = specific_range
        #remove from the structure
        for chain in model.get_chains():
            for residue in chain.get_residues():
                if residue.get_id()[1] < start_residue or residue.get_id()[1] > end_residue:
                    chain.detach_child(residue.id)
    

    chain_id_list = list(string.ascii_uppercase)
    
    model = remove_non_residues(model)
    normalize_fibrils(structure, chain_id_list)
    model = structure[0]  # Reassign model after modifications

    protein_chains = {}
    for chain in model.get_chains():
        protein_chains[chain.id] = chain
    # Calculate the translation vector using centroids
    if len(protein_chains) < 2:
        print("Error: Not enough chains to calculate translation vector.")
        return
    ic(protein_chains)
    last_chain = protein_chains[chain_id_list[len(protein_chains)-1]]
    penultimate_chain = protein_chains[chain_id_list[len(protein_chains)-2]]
    ic(protein_chains)
    ic(last_chain)
    ic(penultimate_chain)
    # Get centroids using all relevant atoms
    centroid_last = np.mean([atom.get_coord() for atom in last_chain.get_atoms()], axis=0)
    centroid_penultimate = np.mean([atom.get_coord() for atom in penultimate_chain.get_atoms()], axis=0)
    
    # Calculate the translation vector
    translation_vector = centroid_last - centroid_penultimate
    ic(translation_vector)
    # Obtain the rotational matrix between the last chain and the penultimate chain

    def get_rotation_matrix(chain1, chain2):
        coords1 = np.array([atom.get_coord() for atom in chain1.get_atoms()])
        coords2 = np.array([atom.get_coord() for atom in chain2.get_atoms()])
        
        # Center the coordinates
        centroid1 = np.mean(coords1, axis=0)
        centroid2 = np.mean(coords2, axis=0)
        coords1 -= centroid1
        coords2 -= centroid2
        
        # Calculate the covariance matrix
        H = np.dot(coords1.T, coords2)
        
        # Singular Value Decomposition
        U, S, Vt = np.linalg.svd(H)
        
        # Calculate the rotation matrix
        R = np.dot(Vt.T, U.T)
        
        # Ensure proper rotation
        if np.linalg.det(R) < 0:
            Vt[-1,:] *= -1
            R = np.dot(Vt.T, U.T)
        
        return R

    # Find the rotation matrix between the last chain and the penultimate chain
    rotation_matrix = get_rotation_matrix(last_chain, penultimate_chain)
    # Ensure all chains have the correct IDs (remove any slashes)

    
    # Extend the fibril
    for layer in range(len(protein_chains), num_total_layers):
        # Duplicate the chain using deepcopy to avoid reference issues
        new_chain = deepcopy(protein_chains[chain_id_list[layer-1]])
        new_chain.id = chain_id_list[layer]
        print(f"Adding chain {new_chain.id}")
        
        # Reset residue numbers
        # for i, residue in enumerate(new_chain.get_residues(), start=1):
        #     residue.id = (' ', i, ' ')

        # Apply the rotation and translation to the new chain
        centroid_new_chain = np.mean([atom.get_coord() for atom in new_chain.get_atoms()], axis=0)
        for atom in new_chain.get_atoms():
            coord = atom.get_coord()
            centered_coord = coord - centroid_new_chain
            rotated_coord = np.dot(rotation_matrix, centered_coord)
            new_coord = rotated_coord + centroid_new_chain + translation_vector
            atom.set_coord(new_coord)

        model.add(new_chain)
        protein_chains[new_chain.id] = new_chain
        #print all the chain
        for chain in model.get_chains():
            print(chain.id)
        # Update the translation vector for the next chain
        translation_vector = np.dot(rotation_matrix, translation_vector)

    #set model number to 1
    model.id = "1"

    # Save the structure with updated chain IDs into PDB format
    io = PDBIO()
    io.set_structure(structure)
    output_file = f"{output_file.split('.')[0]}_{num_total_layers}.pdb"
    print(output_file)
    io.save(output_file) 



