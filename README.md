# List the protein chains and ligands in structural file(s).

This repository contains a Python script that retrieves all downloads and lists all the chains and bound ligands (HETATM) in the structural files (mmCIF format).

The tool downloads updated mmCIF files from PDBe, parse the files and list the UniProt accessions and names (UniProt ID) of each protein chain.

Using this information, the tool generates a PyMOL script to label each protein chain by creating a selection with its UniProt ID and colouring each protein for visualisation.

## Features

- Provides a quick summary of multiple structural files.
- Generates a PyMOL script to allow easy navigation and visualisation of protein complexes.
  
### Input
- **PDB code(s)**: separated by space or comma, in upper or lower cases e.g. `9eih 9eii 9eij`.

### Output
- List all the protein chains in the structural files (protein chains only!).
- List the ligands, resolution and experimental method.
- Generate a PyMOL script for labelling and visualisation.

## Example Usage

<img src="https://github.com/user-attachments/assets/68eb1d63-4635-45cc-8a09-40ce0af048b9" alt="Example Screenshot" width="600"/>

## Prerequisites

- **Python 3.x**
- Python libraries: `requests`, `csv`, `re`, `sys`, `defaultdict` from `collections`.
    
## Author

- **Conny Yu** – [GitHub Profile](https://github.com/connyyu)  
  Haku series 1.0 _April 2025_
