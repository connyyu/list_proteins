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

def parse_uniprot_mapping(cif_file, chains_to_check=None):
    uniprot_mapping = {}
    uniprot_ids = []
    chain_mapping = {}

    fallback_accession = None
    fallback_ID = None
    fallback_chain = None

    with open(cif_file, 'r') as f:
        lines = f.readlines()

    # Scan for fallbacks
    for i, line in enumerate(lines):
        line = line.strip()
        if line.startswith('_struct_ref.pdbx_db_accession'):
            fallback_accession = (line.split()[1] if len(line.split()) > 1
                                  else lines[i + 1].strip().split()[0])
        elif line.startswith('_struct_ref.db_code'):
            fallback_ID = (line.split()[1] if len(line.split()) > 1
                           else lines[i + 1].strip().split()[0])
        elif line.startswith('_struct_ref_seq.pdbx_strand_id'):
            # Try to get chain from current line first
            parts = line.split()
            if len(parts) > 1:
                fallback_chain_candidate = parts[1]
                if len(fallback_chain_candidate) == 1 and fallback_chain_candidate.isalnum():
                    fallback_chain = fallback_chain_candidate
            else:
                # If not found on current line, look at next lines
                k = i + 1
                while k < len(lines):
                    next_line = lines[k].strip()
                    if not next_line or next_line.startswith('_'):
                        break  # End of this data block or new header
                    fallback_chain_candidate = next_line.split()[0]
                    if len(fallback_chain_candidate) == 1 and fallback_chain_candidate.isalnum():
                        fallback_chain = fallback_chain_candidate
                        break
                    k += 1

    i = 0
    process_sifts = False
    sifts_mapping = {}
    sifts_chain_mapping = {}

    while i < len(lines):
        line = lines[i].strip()

        if line.startswith('_struct_ref.pdbx_db_accession'):
            j = i + 1
            while j < len(lines) and not lines[j].startswith('loop_') and lines[j].strip():
                parts = lines[j].strip().split()
                if 'UNP' in parts and len(parts) > 2:
                    entity_id = parts[0]
                    uniprot_id = parts[2]
                    uniprot_accession = None
                    for part in parts:
                        if len(part) >= 6 and '-' not in part:
                            uniprot_accession = part
                            break
                    if uniprot_accession:
                        uniprot_mapping[entity_id] = uniprot_accession
                        uniprot_ids.append((entity_id, uniprot_id))
                j += 1
            i = j
            continue

        elif line.startswith('_struct_ref_seq.align_id'):
            i += 1
            while i < len(lines):
                next_line = lines[i].strip()
                if not next_line or next_line.startswith(('loop_', '#')):
                    break
                parts = next_line.split()
                if len(parts) > 3:
                    entity_id, chain_id = parts[1], parts[3]
                    chain_mapping.setdefault(entity_id, set()).add(chain_id)
                i += 1
            continue

        elif line.startswith('_pdbx_sifts_unp_segments.identity'):
            process_sifts = True

        elif process_sifts:
            if line.startswith(('loop_', '#')):
                process_sifts = False
            elif line:
                parts = line.split()
                if len(parts) > 3:
                    entity_id, chain_id, sifts_accession = parts[0], parts[1], parts[2]
                    if len(sifts_accession) >= 6 and '-' not in sifts_accession:
                        sifts_mapping[entity_id] = sifts_accession
                        sifts_chain_mapping.setdefault(entity_id, set()).add(chain_id)
        i += 1

    # Use fallback if no mapping found
    if not uniprot_mapping:
        # Choose fallback_entity (prefer numeric from existing mappings or fallback to '1')
        if sifts_mapping:
            fallback_entity = sorted(sifts_mapping.keys())[0]
        elif chain_mapping:
            fallback_entity = sorted(chain_mapping.keys())[0]
        else:
            fallback_entity = "1"
    
        # Assign fallback_accession, fallback_ID as before
        uniprot_mapping[fallback_entity] = fallback_accession
        uniprot_ids.append((fallback_entity, fallback_ID))
    
        # Assign all chains for that entity, if any, else fallback_chain if available
        chains_for_entity = chain_mapping.get(fallback_entity)
        if chains_for_entity:
            chain_mapping[fallback_entity] = chains_for_entity
        else:
            chain_mapping[fallback_entity] = set([fallback_chain] if fallback_chain else [])

    # Final reconciliation:
    final_uniprot_mapping = {}
    final_uniprot_ids = []
    final_chain_mapping = {}

    for entity in uniprot_mapping:
        original_ac = uniprot_mapping[entity]
        uniprot_id = next((uid for eid, uid in uniprot_ids if eid == entity), None)
        sifts_ac = sifts_mapping.get(entity, None)

        # Use SIFTS AC if available; else original AC
        final_ac = sifts_ac if sifts_ac else original_ac

        # Use only original chains — NO merging with SIFTS chains
        chains_orig = chain_mapping.get(entity, set())
        final_chains = sorted(list(chains_orig))

        final_uniprot_mapping[entity] = final_ac
        final_uniprot_ids.append((entity, uniprot_id))
        final_chain_mapping[entity] = final_chains

    return final_uniprot_mapping, final_uniprot_ids, final_chain_mapping

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
