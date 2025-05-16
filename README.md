# Docking Pipeline

This repository contains a Python-based docking pipeline for molecular docking simulations using AutoDock Vina. The pipeline automates the process of docking ligands to protein receptors, applying Lipinski's Rule of Five to filter drug-like compounds, and generating complex structures. It outputs the best poses and corresponding ligand-receptor complexes for further analysis.

## Features

- **Ligand Preparation**: Filters ligands based on Lipinski's Rule of Five (Molecular Weight, LogP, Hydrogen Bond Donors/Acceptors).
- **Docking**: Uses AutoDock Vina to dock ligands to receptors.
- **Best Pose Extraction**: Extracts and saves the best docking pose for each ligand-receptor pair.
- **Complex Creation**: Combines the best ligand pose with the receptor structure to create a docking complex.
- **Logging**: Detailed logs for every docking process and successful steps.

## Requirements

- Python 3.x
- AutoDock Vina
- Open Babel (for file conversion)
- RDKit (for molecular descriptors)
