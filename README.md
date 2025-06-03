# List the protein chains and ligands in structural file(s).

This repository contains a Python script that retrieves and lists all the protein chains and bound ligands in the structural files.

The tool downloads updated mmCIF files from PDBe, parse the files and lists the UniProt accessions and names (UniProt ID) of each protein chain. Using this information, the tool generates a PyMOL script to label each protein chain by creating a selection with its UniProt ID and colouring each chain for visualisation.

## Features

- Provides a quick summary of multiple structural files.
- Generates a PyMOL script to allow easy navigation and visualisation of protein complexes.
- Limitation: not compatible with chimeric and non-protein constructs.
  
### Input
- **PDB code(s)**: separated by space or comma, in upper or lower cases e.g. `python haku_list_proteins.py 9eih 9eii 9eij`.

### Output
- Download the updated mmCIF files (pdb_list.txt)
- Parse the structural files into a CSV file (parsed_cif.csv)
- List the UniProt AC and ID for each protein chain, ligands, resolution and experimental method
- Generate a PyMOL script for labelling and visualisation (pdb_pymol.txt)

## Example Usage

<img src="https://github.com/user-attachments/assets/68eb1d63-4635-45cc-8a09-40ce0af048b9" alt="Example Screenshot" width="600"/>

## Prerequisites

- **Python 3.x**
- Python libraries: `requests`, `csv`, `re`, `sys`, `defaultdict` from `collections`.
    
## Author

- **Conny Yu** – [GitHub Profile](https://github.com/connyyu)  
  Haku series 1.1 _May 2025_
