#Script to convert the current fasta file format into a new fasta file format where the titles is removed
import os
import yaml
existing_sequences = []

def convert_fasta(protein):
    # Check if seqs directory exists
    if not os.path.exists(f"{protein}/seqs"):
        print(f"Error: Directory {protein}/seqs does not exist")
        return

    # Check if input fasta file exists
    input_fasta = f"{protein}/seqs/{protein.split('/')[-1]}.fa"
    if not os.path.exists(input_fasta):
        print(f"Error: Input file {input_fasta} does not exist")
        return

    # Create fasta_files directory if it doesn't exist
    if not os.path.exists(f"{protein}/fasta_files"):
        os.makedirs(f"{protein}/fasta_files")

    # Create a new file in proteinsName
    with open(input_fasta, "r") as file:
        lines = file.readlines()
        starting = 2

        for i in range(starting, len(lines), 2):
            name = f"{protein.split('/')[-1]}_{i//2}"
            if lines[i+1] in existing_sequences:
                continue
            existing_sequences.append(lines[i+1])
            output_file = f"{protein}/fasta_files/{name}.fa"
            try:
                with open(output_file, "w") as file:
                    file.write(f">{name}\n")
                    file.write(lines[i+1])
            except IOError as e:
                print(f"Error writing to {output_file}: {e}")

#Find all fasta files in the output folder from various temps and combine them into one fasta file
def combine_fasta_files(path, temps):
    existing_sequences = set()  # Changed to set for faster lookups
    # Get all fasta files in the output folder from various temps
    fasta_files = []
    for temp in temps:
        for file in os.listdir(f"{path}/{temp}/seqs"):
            if file.endswith(".fa"):
                with open(f"{path}/{temp}/seqs/{file}", "r") as file:
                    prev_line = None
                    #skip the two lines of the fasta file
                    file = file.readlines()[2:]
                    for line in file:
                        if not line.startswith(">"):
                            sequence = line.strip().split("/")[0].upper()
                            if sequence not in existing_sequences:
                                existing_sequences.add(sequence)
                                fasta_files.append(prev_line.strip() + '\n' + sequence + '\n')
                        else:
                            prev_line = line
    #write to a new file
    print(path)
    with open(f"{path}/combined_fasta.fa", "w") as file:
        for fasta in fasta_files:
            file.write(fasta)
    print(f"Combined fasta files saved to {path}/combined_fasta.fa")
    print(f"Total unique sequences: {len(fasta_files)}")
    return fasta_files  

# if __name__ == "__main__":
#     path = "/work/CAND/shared/Chens/AD_fibrils/Testing"
#     temps = ["0.3", "0.5"]
#     fasta_files = combine_fasta_files(path, temps)

# if __name__ == "__main__":
#     with open("/work/CAND/shared/Chens/AD_fibrils/mpnn_config.yaml", "r") as file:
#         config = yaml.safe_load(file)
#         output_folder = config["filepath"]["out_folder"]
#     for protein in config["mpnn_cleaning"]["proteins"]:
#         convert_fasta(f"{output_folder}/{protein}")