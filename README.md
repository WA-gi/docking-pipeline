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

The docking pipeline follows a systematic sequence of steps to dock ligands to protein targets and generate the corresponding ligand-receptor complexes. Below is a detailed breakdown of the steps involved:
1. Ligand Collection

    The first step involves gathering all ligand files in .pdbqt format from the ligands/ directory. These files represent the 3D structures of the small molecules (ligands) to be docked against protein receptors.

    If any ligands have corresponding .sdf files in the ligands/ directory, they are further processed to filter out any compounds that don't meet Lipinski's Rule of Five (drug-likeness criteria). This includes filtering based on molecular weight (MW), octanol-water partition coefficient (LogP), and hydrogen bond donors/acceptors.

2. Filtering Ligands

    A drug-likeness filter is applied to ensure the ligands meet the following criteria:

        Molecular Weight (MW): Less than or equal to 500 Da.

        LogP: Less than or equal to 5 (this represents the hydrophobicity of the molecule).

        Hydrogen Bond Donors (HDonors): Less than or equal to 5.

        Hydrogen Bond Acceptors (HAcceptors): Less than or equal to 10.

    Only ligands that pass these criteria are retained for docking unless the user chooses to skip the filter and dock all available ligands.

3. Grid Box Setup

    Before docking can begin, the docking region on the receptor (target) must be defined. This is achieved through a grid box configuration (gdf.txt), which specifies the size and center of the docking grid in three-dimensional space.

    If the gdf.txt file is missing or invalid, the user is prompted to manually input the grid box parameters:

        Center Coordinates (X, Y, Z): The coordinates that define the center of the grid box.

        Grid Box Dimensions (X, Y, Z): The size of the docking region in each dimension.

4. File Conversion (PDBQT to PDB)

    If needed, ligand and receptor files are converted from .pdbqt to .pdb format using Open Babel. This conversion allows for easier visualization and analysis of the docked poses and complexes.

5. Docking Process with AutoDock Vina

    The docking itself is performed using AutoDock Vina, a popular molecular docking software. For each ligand-target pair, the following process occurs:

        AutoDock Vina docks the ligand to the receptor (target) using the previously defined grid box and parameters.

        The command executed includes the ligand and receptor file paths, along with the grid box dimensions and center.

        The docking is performed, and the results are saved in .pdbqt format, including a log file (.txt) containing the docking information.

6. Pose Extraction and Best Pose Selection

    After docking, the results contain multiple possible poses (conformations of the ligand within the receptor’s binding site). The best pose is selected based on the docking score, which represents the binding affinity of the ligand-receptor complex.

    The pipeline extracts the best pose (lowest energy score) from the .pdbqt file and saves it for further use.

    The best pose is then converted to the .pdb format for easier analysis.

7. Building Ligand-Receptor Complexes

    Once the best pose is identified, the ligand is combined with the receptor to form a ligand-receptor complex. This complex is saved in .pdb format and can be used for further analysis, such as visualizing the interaction between the ligand and receptor or running additional simulations.

8. Logging and Results

    Throughout the docking process, detailed logs are generated to track the progress and identify any errors or issues encountered. These logs provide information such as:

        The number of ligands and receptors processed.

        Errors or issues with ligand or receptor files.

        Success or failure of each docking attempt.

    The final output includes:

        The docked ligand poses in the out/ directory.

        The best poses in the best_poses/ directory.

        The complex structures in the complexes/ directory.

        The logs in the Log/ directory for further review.

9. Final Output

    Ligand-Target Complexes: For each ligand-receptor pair, a complex file is saved in the complexes/ folder. This file represents the docked structure of the ligand bound to the receptor.

    Best Pose Files: The best docking pose for each ligand-receptor pair is saved in the best_poses/ directory in both .pdbqt and .pdb formats.

    Log Files: Each docking attempt generates a log file that contains the docking score and any relevant information or errors. These are saved in the Log/ directory.
