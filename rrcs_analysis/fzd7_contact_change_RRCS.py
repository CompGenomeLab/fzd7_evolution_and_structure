# -*- coding: utf-8 -*-
"""
Created on Sun Dec 31 15:08:04 2023

@author: selcuk.1
"""
def remove_hydrogens_from_pdb(file_path):
    """
    Removes hydrogen atoms from a PDB file.

    Parameters:
    file_path (str): The directory path to the PDB file.
    """
    try:
        # Open the original file in read mode and a temporary file in write mode
        with open(file_path, 'r') as file:
            lines = file.readlines()

        # Filter out lines that represent hydrogen atoms
        new_lines = [line for line in lines if not (line.startswith('ATOM') and line.strip().endswith('H'))]

        # Write the filtered lines back to the original file
        file_path=file_path[:-4]+"H_removed.pdb"
        with open(file_path, 'w') as file:
            file.writelines(new_lines)
        
        print("Hydrogen atoms successfully removed.")
    except Exception as e:
        print(f"An error occurred: {e}")


remove_hydrogens_from_pdb(r"C:\Users\selcuk.1\OneDrive - The Ohio State University\Desktop\ClassF_files\colab_structural\Mono432-coot-0.pdb")


#%%
def parse_line(line):
    """Parse a line and return residue information and contact score."""
    parts = line.split()
    chain_id1, residue_info1 = parts[0].split(':')
    chain_id2, residue_info2 = parts[1].split(':')
    residue_num1 = residue_info1.split('_')[0]
    residue_num2 = residue_info2.split('_')[0]
    contact_score = float(parts[2])
    return (chain_id1, residue_num1, chain_id2, residue_num2, contact_score)

def read_file(file_path):
    """Read a file and return a dictionary of contact information."""
    contacts = {}
    if "Mono" in file_path:
        chain="A"
    else:
        chain="R"
    with open(file_path, 'r') as file:
        for line in file:
            chain_id1, residue_num1, chain_id2, residue_num2, contact_score = parse_line(line)
            if chain_id1 == chain and chain_id2 == chain:
                contacts[(residue_num1, residue_num2)] = contact_score
    return contacts

def calculate_differences(contacts_active, contacts_inactive):
    """Calculate differences in contact scores between active and inactive states, considering all pairs."""
    differences = {}
    # Create a set of all unique pairs from both active and inactive states
    all_pairs = set(contacts_active.keys()) | set(contacts_inactive.keys())
    
    for pair in all_pairs:
        if pair in contacts_active and pair in contacts_inactive:
            # Calculate difference if pair is in both active and inactive
            difference = contacts_active[pair] - contacts_inactive[pair]
        elif pair in contacts_active:
            # If pair is only in active, use its score directly
            difference = contacts_active[pair]
        else:
            # If pair is only in inactive, multiply its score by -1
            difference = -contacts_inactive[pair]
        
        differences[pair] = (difference, abs(difference))
    
    return differences


def write_output(differences, output_file):
    """Write the calculated differences and their absolute values to an output file."""
    with open(output_file, 'w') as file:
        for pair, (diff, abs_diff) in differences.items():
            file.write(f"{pair[0]}\t{pair[1]}\t{diff}\t{abs_diff}\n")

def process_contacts(file_active, file_inactive, output_file):
    contacts_active = read_file(file_active)
    contacts_inactive = read_file(file_inactive)
    differences = calculate_differences(contacts_active, contacts_inactive)
    write_output(differences, output_file)
    
def calculate_total_change_per_residue(input_file, output_file):
    """
    Calculate the total absolute contact change score for each residue.
    The input file is expected to have columns: residue_number_1, residue_number_2, difference, absolute_difference.
    """
    residue_changes = {}

    # Read the input file and accumulate changes for each residue
    with open(input_file, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            residue1, residue2, _, abs_change = parts
            abs_change = float(abs_change)

            # Update the total change for each residue
            residue_changes[residue1] = residue_changes.get(residue1, 0) + abs_change
            residue_changes[residue2] = residue_changes.get(residue2, 0) + abs_change

    # Write the total changes to the output file
    with open(output_file, 'w') as file:
        for residue, total_change in residue_changes.items():
            file.write(f"{residue}\t{total_change}\n")


inactive=r"C:\Users\selcuk.1\OneDrive - The Ohio State University\Desktop\ClassF_files\colab_structural\Mono432-coot-0H_removed.pdb.cscore"
active=r"C:\Users\selcuk.1\OneDrive - The Ohio State University\Desktop\ClassF_files\colab_structural\7evwb-coot-5_real_space_refined_004.pdb.cscore"
output=r"C:\Users\selcuk.1\OneDrive - The Ohio State University\Desktop\ClassF_files\colab_structural\fzd7_activation_changes.cscore"
process_contacts(active,inactive,output)
change_per_residue=r"C:\Users\selcuk.1\OneDrive - The Ohio State University\Desktop\ClassF_files\colab_structural\fzd7_activation_change_per_residue.tsv"
calculate_total_change_per_residue(output, change_per_residue)
