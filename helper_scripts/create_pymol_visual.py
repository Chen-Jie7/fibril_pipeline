
#based on the distribution, generate a pymol command to color the residues based on selected feature
def generate_pymol_bfactor_commands(temps, distribution_list, protein_name="6HRE_5_min", feature="mutation_rate"):
    """
    Generates PyMOL commands to color residues based on probability values using B-factors
    for multiple temperature copies.
    
    Parameters:
    -----------
    temps : list
        List of temperature values
    distribution_list : list
        List of dictionaries containing residue probabilities for each temperature
        
    Returns:
    --------
    str
        Complete PyMOL command string
    """
    # Initialize command string
    commands = "set grid_mode,1\n"
    
    maximum = 0
    # Generate commands for each temperature
    for temp, distributions in zip(temps, distribution_list):
        # Create copy of protein for this temperature
        commands += f"copy {protein_name}_{temp}, {protein_name}\n"
        
        # Add B-factor commands for each residue
        commands += f"# PyMOL commands for coloring residues by probability using B-factors\n"
        for position, body_dict in distributions.items():
            bfactor = body_dict[feature] * 100
            if bfactor > maximum:
                maximum = bfactor
            commands += f"alter {protein_name}_{temp} and resi {position}, b={bfactor};\n"
            
    # Add coloring command
    commands += "\n# Color by B-factor\n"
    commands += f"spectrum b, rainbow, minimum=0, maximum={maximum+0.5};\n"
        
    # Add visualization settings
    commands += "\n# Basic visualization settings\n"
    commands += "as cartoon;\n"
    commands += "show sticks;\n"
    commands += "zoom;\n"

    with open(f"{protein_name}_bfactor_command.txt", "w") as f:
        f.write(commands)
    print(f"PyMOL command saved to {protein_name}_bfactor_command.txt")




