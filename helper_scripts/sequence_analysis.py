#import the config file
import yaml
import sys
import os
from mpnn_cleaning import combine_fasta_files
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import pairwise2
from Bio.Seq import Seq
import pandas as pd
import multiprocessing as mp
from functools import partial
import argparse

#calculate the percent identity between two sequences
def calculate_percent_identity(seq1, seq2):
    """Calculate percent identity between two sequences."""
    # Perform global alignment with simple scoring (match=1, mismatch=0, gap=-1)
    alignments = pairwise2.align.globalms(seq1, seq2, 1, 0, -1, -1, one_alignment_only=True)
    if not alignments:
        return 0
    
    alignment = alignments[0]
    seq1_aligned, seq2_aligned = alignment[0], alignment[1]
    
    # Count matches
    matches = sum(c1 == c2 for c1, c2 in zip(seq1_aligned, seq2_aligned) 
                 if c1 != '-' and c2 != '-')
    
    # Calculate percent identity based on alignment length excluding gaps
    alignment_length = sum(1 for c1, c2 in zip(seq1_aligned, seq2_aligned) 
                         if c1 != '-' or c2 != '-')
    
    if alignment_length == 0:
        return 0
    
    return (matches / alignment_length) * 100

# Helper function for parallel processing
def calculate_row(i, sequences, n):
    """Calculate percent identity for row i of the matrix."""
    row = np.zeros(n)
    for j in range(n):
        if i == j:
            row[j] = 100  # Self-comparison
        else:
            row[j] = calculate_percent_identity(sequences[i], sequences[j])
    return i, row

def create_percent_identity_matrix(sequences, names=None, path=None, num_processes=None):
    """Generate a percent identity matrix for a list of sequences using parallel processing."""
    n = len(sequences)
    identity_matrix = np.zeros((n, n))
    
    # Determine number of processes to use
    if num_processes is None:
        num_processes = mp.cpu_count()
    num_processes = min(num_processes, n)  # Don't use more processes than sequences
    
    print(f"Calculating percent identity matrix using {num_processes} processes...")
    
    # Create a pool of workers
    with mp.Pool(processes=num_processes) as pool:
        # Create a partial function with fixed parameters
        func = partial(calculate_row, sequences=sequences, n=n)
        # Map the function to all row indices and collect results
        results = pool.map(func, range(n))
    
    # Assemble the matrix from results
    for i, row in results:
        identity_matrix[i] = row
    
    # Create a DataFrame if names are provided
    if names:
        identity_matrix = pd.DataFrame(identity_matrix, index=names, columns=names)
    
    # Save the identity matrix to a csv file
    if path:
        identity_matrix.to_csv(f'{path}/identity_matrix.csv')
        print(f'Identity matrix saved to {path}/identity_matrix.csv')
    
    return identity_matrix

def visualize_percent_identity_matrix(identity_matrix, title="Percent Identity Matrix", path=None):
    """Visualize the percent identity matrix as a heatmap."""
    plt.figure(figsize=(10, 8))
    
    if isinstance(identity_matrix, pd.DataFrame):
        sns.heatmap(identity_matrix, annot=True, cmap='viridis', vmin=0, vmax=100)
    else:
        sns.heatmap(identity_matrix, annot=True, cmap='viridis', vmin=0, vmax=100)
    
    plt.title(title)
    plt.tight_layout()
    
    if path:
        plt.savefig(f'{path}/identity_matrix_heatmap.png')
        print(f'Heatmap saved to {path}/identity_matrix_heatmap.png')
    else:
        plt.show()

def create_slurm_script(config_file, num_processes=None):
    """Create a SLURM submission script for running the analysis."""
    script_content = f"""#!/bin/bash
#SBATCH --job-name=seq_analysis
#SBATCH --output=seq_analysis_%j.out
#SBATCH --error=seq_analysis_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task={num_processes if num_processes else 4}
#SBATCH --mem=16G

python {os.path.abspath(__file__)} --config {config_file} --processes {num_processes if num_processes else 4}
"""
    
    script_path = f"run_seq_analysis_{os.path.basename(config_file).split('.')[0]}.sh"
    with open(script_path, 'w') as f:
        f.write(script_content)
    
    print(f"SLURM script created: {script_path}")
    print(f"Submit with: sbatch {script_path}")
    
    return script_path

#run the script
if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Calculate sequence identity matrix in parallel')
    parser.add_argument('--config', help='Config file name (without .yaml extension)')
    parser.add_argument('--processes', type=int, default=None, 
                        help='Number of processes to use (default: number of CPU cores)')
    parser.add_argument('--submit', action='store_true', 
                        help='Create and submit SLURM job instead of running directly')
    
    # Parse command line arguments or use sys.argv
    if len(sys.argv) > 1 and sys.argv[1].startswith('--'):
        args = parser.parse_args()
        config_name = args.config
        num_processes = args.processes
        submit_job = args.submit
    else:
        # Legacy support for old command line format
        input_args = sys.argv[1:]
        config_name = input_args[0]
        num_processes = None
        submit_job = False
    
    # Create SLURM script if requested
    if submit_job:
        script_path = create_slurm_script(config_name, num_processes)
        # Submit the job
        subprocess.run(['sbatch', script_path])
        sys.exit(0)
    
    # Load the config file
    config_path = f'/work/CAND/shared/Chens/AD_fibrils/{config_name}.yaml'
    with open(config_path, 'r') as file:
        config = yaml.safe_load(file)
    
    # Get sequences
    path = os.path.join(config['filepath']['out_folder'], config['mpnn_cleaning']['proteins'][0])
    fasta_files = list(combine_fasta_files(path, config['mpnn_cleaning']['temperatures']))
    names = [f"Seq{i}" for i in range(len(fasta_files))]
    
    # Calculate identity matrix
    identity_matrix = create_percent_identity_matrix(fasta_files, names, path, num_processes)
