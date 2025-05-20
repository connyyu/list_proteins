import csv
import re
import requests
import sys
from collections import defaultdict

def download_cif_file(pdb_code):
    pdb_code_dl = pdb_code.lower()
    url = f"https://www.ebi.ac.uk/pdbe/entry-files/download/{pdb_code_dl}_updated.cif"
    response = requests.get(url)
    if response.status_code == 200:
        with open(f"{pdb_code_dl}.cif", 'wb') as f:
            f.write(response.content)        
        return f"{pdb_code_dl}.cif"
    else:
        return None

def parse_resolution(cif_file):
    resolution = None
    method = ""
    symmetry_type_found = False

    with open(cif_file, 'r') as f:
        lines = f.readlines()

    for i, line in enumerate(lines):
        if line.startswith('_reflns.d_resolution_high'):
            parts = line.split()
            if len(parts) > 1:
                resolution = parts[1]
                method = "Xtal"
                break
        elif line.startswith('_em_3d_reconstruction.resolution'):
            parts = line.split()
            if len(parts) > 1:
                resolution = parts[1]
                method = "EM"
                break
        elif line.startswith("_exptl.method 'Solution NMR'"):
            method = "NMR"
            break
        elif line.startswith('_em_3d_reconstruction.symmetry_type'):
            symmetry_type_found = True
        elif symmetry_type_found:  # Look at the next line after symmetry_type
            parts = line.split()
            if len(parts) >= 7:
                resolution = parts[6]  # Take the 7th column as resolution
                method = "EM"
            break
    # 2 decimal places for resolution
    if resolution is not None:
        try:
            resolution = f"{float(resolution):.2f}"
        except ValueError:
            resolution = None

    return resolution, method

def parse_hetatms(cif_file):
    hetatms = set()
    with open(cif_file, 'r') as f:
        lines = f.readlines()
    for line in lines:
        if line.startswith('HETATM'):
            parts = line.split()
            hetatm_name = parts[5]  # HETATM name (e.g., GLC, ZN)
            hetatms.add(hetatm_name)  # Use a set to avoid duplicates
    return sorted(hetatms)  # Return a sorted list

def fetch_uniprot_id(uniprot_ac):
    url = f'https://rest.uniprot.org/uniprotkb/{uniprot_ac}.json'
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        return data.get('uniProtkbId', None)
    return None

def parse_uniprot_mapping(cif_file):
    uniprot_mapping = {}
    uniprot_ids = []
    chain_mapping = {}
    seen_chain_ids = {}  # To track which chains are associated with each UniProt accession
    
    with open(cif_file, 'r') as f:
        lines = f.readlines()

    i = 0
    process_sifts = False
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith('_struct_ref.pdbx_db_accession'):
            found_valid_line = False
            j = i + 1
            while j < len(lines):
                if 'UNP' in lines[j]:  # Check if 'UNP' is in the line
                    parts = lines[j].strip().split()
                    if 'UNP' in parts and parts[1] == 'UNP':  # Ensure 'UNP' is in the expected position
                        for part in parts:
                            if len(part) == 6:
                                uniprot_accession = part  # Get UniProt accession
                                uniprot_id = parts[2]  # Get UniProt ID
                                entity_id = parts[0]  # Get entity_id
                                if entity_id not in uniprot_mapping:
                                    uniprot_mapping[entity_id] = uniprot_accession
                                    uniprot_ids.append((entity_id, uniprot_id))
                                    found_valid_line = True
                j += 1
            if not found_valid_line:
                for k in range(len(lines)):
                    if lines[k].startswith('_struct_ref.entity_id'):
                        entity_id = lines[k].strip().split()[-1]
                    if lines[k].startswith('_struct_ref.pdbx_db_accession'):
                        uniprot_accession = lines[k].strip().split()[-1]
                    if lines[k].startswith('_struct_ref.db_code'):
                        uniprot_id = lines[k].strip().split()[-1]
                if entity_id and uniprot_accession:
                    if entity_id not in uniprot_mapping:
                        uniprot_mapping[entity_id] = uniprot_accession
                        uniprot_ids.append((entity_id, uniprot_id))
        elif line.startswith('_struct_ref_seq.align_id'):
            i += 1  # Skip the header line
            while i < len(lines) and not lines[i].strip().startswith('loop_'):
                parts = lines[i].strip().split()
                if len(parts) > 3:
                    entity_id = parts[1]
                    chain_id = parts[3]
                    if entity_id in uniprot_mapping:
                        if entity_id not in chain_mapping:
                            chain_mapping[entity_id] = set()  # Initialize as a set to avoid duplicates
                        chain_mapping[entity_id].add(chain_id)  # Add the chain ID to the set of chains for this entity ID
                i += 1
            continue
        # Update the UniProt ID and AC with the new mappings in SIFTS
        elif line.startswith('_pdbx_sifts_unp_segments.identity'):
            process_sifts = True
        elif process_sifts and line.startswith('loop'):
            process_sifts = False
        elif process_sifts:
            parts = line.strip().split()
            if len(parts) > 3:
                entity_id = parts[0]
                uniprot_accession = parts[2]
                if len(uniprot_accession) == 6 and entity_id in uniprot_mapping:
                    uniprot_mapping[entity_id] = uniprot_accession
                    new_uniprot_id = fetch_uniprot_id(uniprot_accession)
                    if new_uniprot_id:
                        # Update existing uniprot_id for the entity_id
                        for idx, (eid, uid) in enumerate(uniprot_ids):
                            if eid == entity_id:
                                uniprot_ids[idx] = (entity_id, new_uniprot_id)
                                break  # Stop after finding the match

        i += 1
    
    chain_mapping = {entity_id: sorted(list(chains)) for entity_id, chains in chain_mapping.items()}
    
    return uniprot_mapping, uniprot_ids, chain_mapping

def check_chimera(table_data):
    chain_uniprot_map = defaultdict(set)
    for row in table_data:
        pdb_code = row['PDB'].upper()
        chain_id = row['Chain ID']
        uniprot_accession = row['UniProt AC']
        chain_uniprot_map[chain_id].add(uniprot_accession)
    warnings = []
    for chain_id, uniprot_accessions in chain_uniprot_map.items():
        if len(uniprot_accessions) > 1:
            warnings.append(f"{pdb_code} chain {chain_id} is associated with mutiple UniProt ACs: {', '.join(uniprot_accessions)}.")
    return warnings
    
def process_cif_file(pdb_code):
    cif_file = download_cif_file(pdb_code)
    if not cif_file:
        print(f"\033[1mFailed to download CIF file for {pdb_code}. Skipping...\033[0m")
        return

    resolution, method = parse_resolution(cif_file)    
    hetatms = parse_hetatms(cif_file)
    uniprot_mapping, uniprot_ids, chain_mapping = parse_uniprot_mapping(cif_file)

    hetatm_names = ', '.join(hetatms)
    table_data = []
    
    printed_chains = set()
    for entity_id, uniprot_accession in uniprot_mapping.items():
        uniprot_id = next((id for ent_id, id in uniprot_ids if ent_id == entity_id), None)
        chain_ids = chain_mapping.get(entity_id, "A")
        for chain_id in chain_ids:
            if (uniprot_accession, chain_id) in printed_chains:
                continue
            printed_chains.add((uniprot_accession, chain_id))
            table_data.append({
                'PDB': pdb_code.lower(),
                'Entity': entity_id,
                'Chain ID': chain_id,
                'UniProt AC': uniprot_accession,
                'UniProt ID': uniprot_id,
                'HETATM': hetatm_names,
                'Resolution': resolution,
                'Method': method
            })

    # Check for warnings and print them
    warnings = check_chimera(table_data)
    if warnings:
        print("\n\033[1;31mChimera warning!!\033[0m")
        for warning in warnings:
            print(f"\033[1;31m{warning}\033[0m")
    
    return table_data

def pymol_command(csv_file):
    # List of predefined colors
    colors = [
        "carbon", "cyan", "lightmagenta", "yellow", "salmon", "slate", "orange",
        "deepteal", "violetpurple", "hydrogen", "marine", "olive", "smudge",
        "teal", "wheat", "lightpink", "skyblue"
    ]
    
    uniprot_chains = defaultdict(lambda: defaultdict(list))
    with open(csv_file, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            pdb_id = row['PDB']
            uniprot_id = row['UniProt ID']
            chain_id = row['Chain ID']
            # Add the chain to the respective UniProt ID and PDB ID
            uniprot_chains[uniprot_id][pdb_id].append(chain_id)
    
    # Open the output file for writing
    with open("pdb_pymol.txt", "w") as outfile:
        outfile.write("util.cbaw;\n")
        color_index = 0  # Start with the first color in the list
        
        reformatted_selections = {}
        for uniprot_id, pdb_dict in uniprot_chains.items():
            selections = []
            for pdb_id, chains in pdb_dict.items():
                chains_str = "+".join(chains)
                selections.append(f"(chain {chains_str} and {pdb_id})")
            
            reformatted_selections[uniprot_id] = " or ".join(selections)
        
        for uniprot_id, selection in reformatted_selections.items():
            color = colors[color_index % len(colors)]
            outfile.write(f"sele {uniprot_id}, {selection}; color {color}, {uniprot_id};\n")
            color_index += 1  # Move to the next color
        
        outfile.write("deselect; util.cnc")
    
    # Print out the reformatted content
    with open("pdb_pymol.txt", "r") as outfile:
        print(outfile.read())
    
def print_simplified_output(csv_file):
    # defaultdict with list to collect chain IDs
    data = defaultdict(list)
    
    # Read the parsed CSV file and group the chain IDs by UniProt accession and PDB ID
    with open(csv_file, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            pdb_id = row['PDB']
            entity_id = row['Entity']
            chain_id = row['Chain ID']
            uniprot_ac = row['UniProt AC']
            uniprot_id = row['UniProt ID']
            hetatm = row['HETATM']
            resolution = row['Resolution'] 
            method = row['Method']

            # Group chain IDs by PDB ID and UniProt accession
            data[(pdb_id, entity_id, uniprot_ac, uniprot_id, hetatm, resolution, method)].append(chain_id)
    
    # Print the simplified output
    print("")
    print(f"{'PDB':<5} {'Entity':<6} {'Chain ID':<20} {'UniProt AC':<11} {'UniProt ID':<15} {'HETATM':<30} {'Resolution':<12} {'Method':<8}")
    print("-" * 120)
    for (pdb_id, entity_id, uniprot_ac, uniprot_id, hetatm, resolution, method), chain_ids in data.items():
        # Combine the chain IDs into a single string
        chain_ids_str = '/'.join(sorted(chain_ids))
        print(f"{pdb_id:<5} {entity_id:<6} {chain_ids_str:<20} {uniprot_ac:<11} {uniprot_id:<15} {hetatm:<30} {resolution:<12} {method:<8}")

def main():
    # Clear pdb_codes from last run
    open('pdb_codes.txt', 'w').close()
    
    if len(sys.argv) > 1:
        pdb_codes = [code.strip().upper() for code in sys.argv[1:]]
    else:
        pdb_codes = re.split(r'[,\s]+', input("Enter PDB codes (comma or space-separated): ").strip().upper())
    
    pdb_codes = [code.strip() for code in pdb_codes]

    with open('pdb_codes.txt', 'w') as f:
        f.write('\n'.join(pdb_codes))

    # Stop execution if pdb_codes.txt is empty
    if not pdb_codes or all(code == "" for code in pdb_codes):
        print("\033[1mNo PDB codes provided.\n<End>\033[0m")
        return
    
    all_data = []
    successful_downloads = 0
    
    for pdb_code in pdb_codes:
        print(f"Processing {pdb_code}...")
        table_data = process_cif_file(pdb_code)
        if table_data:
            all_data.extend(table_data)
            successful_downloads += 1
    
    # Stop execution if no valid PDBs were processed
    if successful_downloads == 0:
        print("\033[1mNo valid PDB files were found.\n<End>\033[0m")
        return

    with open('parsed_cif.csv', 'w', newline='') as csvfile:
        fieldnames = ['PDB', 'Entity', 'Chain ID', 'UniProt AC', 'UniProt ID', 'HETATM', 'Resolution', 'Method']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_data)
    
    print_simplified_output('parsed_cif.csv')
    print()
    pymol_command('parsed_cif.csv')

if __name__ == "__main__":
    main()
