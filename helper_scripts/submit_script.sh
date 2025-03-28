#!/bin/bash
#SBATCH --job-name=sequence_analysis
#SBATCH --partition=super
#SBATCH --nodes=1
#SBATCH --time=5:00:00
#SBATCH --output=/work/CAND/shared/Chens/AD_fibrils/6HRE_5_min/logs/submit_jobs_%j.out
#SBATCH --error=/work/CAND/shared/Chens/AD_fibrils/6HRE_5_min/logs/submit_jobs_%j.err

source activate /endosome/work/CAND/shared/Chens/miniforge3/envs/upside2
# Get input arguments as a string
input_args="$1"

# Split the input string into an array
IFS=' ' read -ra args_array <<< "$input_args"

# Print array elements for verification
echo "Input arguments split into array:"
for i in "${!args_array[@]}"; do
    echo "args_array[$i] = ${args_array[$i]}"
done

python sequence_analysis.py "${args_array[@]}"