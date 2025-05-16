import os
import subprocess
from rdkit import Chem
from rdkit.Chem import Descriptors
from datetime import datetime

# === CONFIGURATION ===
ligand_dir = "ligands"
target_dir = "target"
out_dir = "out"
log_dir = "Log"
pdb_dir = "converted"
complex_dir = "complexes"
best_pose_dir = "best_poses"

# === MAKE DIRECTORIES ===
for d in [out_dir, log_dir, pdb_dir, complex_dir, best_pose_dir]:
    os.makedirs(d, exist_ok=True)

# === LOGGING FUNCTION ===
def log(message):
    time = datetime.now().strftime("[%Y-%m-%d %H:%M:%S]")
    print(time, message)
    with open("run_log.txt", "a") as f:
        f.write(time + " " + message + "\n")

# === LIGAND COLLECTION ===
ligands_all = [f for f in os.listdir(ligand_dir) if f.endswith(".pdbqt")]
targets = [f for f in os.listdir(target_dir) if f.endswith(".pdbqt")]

log(f"Found {len(ligands_all)} ligands and {len(targets)} targets.")

# === DRUG-LIKENESS FILTER ===
def is_druglike(mol):
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    h_donors = Descriptors.NumHDonors(mol)
    h_acceptors = Descriptors.NumHAcceptors(mol)
    return mw <= 500 and logp <= 5 and h_donors <= 5 and h_acceptors <= 10

filtered_ligands = []
for lig in ligands_all:
    sdf_path = os.path.join(ligand_dir, lig.replace(".pdbqt", ".sdf"))
    if not os.path.exists(sdf_path):
        continue
    mol = Chem.SDMolSupplier(sdf_path)[0]
    if mol and is_druglike(mol):
        filtered_ligands.append(lig)
    else:
        log(f"❌ Skipped {lig} (fails Lipinski filter)")

log(f"🧪 {len(filtered_ligands)} ligands passed Lipinski filter out of {len(ligands_all)} total.")

# === USER FILTER CHOICE ===
choice = input("👉 Dock only filtered ligands (Y) or all ligands (N)? [Y/N]: ").strip().upper()
if choice == "Y":
    ligands = filtered_ligands
    log("📦 Proceeding with filtered ligands only.")
else:
    ligands = ligands_all
    log("📦 Proceeding with all ligands (Lipinski filter ignored).")

# === GRID BOX PARAMETERS ===
try:
    with open("gdf.txt", "r") as sdf:
        lines = sdf.readlines()
        x_size, y_size, z_size = lines[2].split()[1:4]
        x_center, y_center, z_center = lines[3].split()[1:4]
        log("✅ Grid parameters loaded.")
except Exception as e:
    log(f"⚠️ Grid box config missing or invalid. Error: {e}")
    x_center = input("Center X: ")
    y_center = input("Center Y: ")
    z_center = input("Center Z: ")
    x_size = input("Size X: ")
    y_size = input("Size Y: ")
    z_size = input("Size Z: ")

# === FILE CONVERSION ===
def convert_to_pdb(pdbqt_path, pdb_path):
    result = subprocess.run(["obabel", "-ipdbqt", pdbqt_path, "-O", pdb_path],
                            capture_output=True, text=True)
    if result.returncode != 0:
        log(f"❌ Open Babel failed converting {pdbqt_path}:\n{result.stderr.strip()}")
    elif not os.path.exists(pdb_path):
        log(f"❌ PDB file not created for {pdbqt_path}")

# === BEST POSE EXTRACTION FROM PDBQT ===
def extract_best_pose_from_pdbqt(pdbqt_path):
    with open(pdbqt_path, "r") as file:
        lines = file.readlines()

    poses = []
    current_pose = []
    score = None

    for line in lines:
        if line.startswith("MODEL"):
            current_pose = [line]
            score = None
        elif "REMARK VINA RESULT:" in line:
            try:
                score = float(line.strip().split()[3])
            except Exception:
                score = None
            current_pose.append(line)
        elif line.startswith("ENDMDL"):
            current_pose.append(line)
            if score is not None:
                poses.append((score, list(current_pose)))
        else:
            current_pose.append(line)

    if not poses:
        return None, None

    best_score, best_pose_lines = sorted(poses, key=lambda x: x[0])[0]
    return best_score, best_pose_lines

# === COMPLEX BUILDER ===
def build_complex(ligand_pdb, receptor_pdb, output_path):
    with open(output_path, "w") as out_f:
        with open(receptor_pdb) as r:
            out_f.write(r.read())
        with open(ligand_pdb) as l:
            out_f.write("\n")
            out_f.write(l.read())

# === DOCKING AND COMPLEX CREATION ===
log("🚀 Starting sequential docking and complex creation...")

for ligand in ligands:
    for receptor in targets:
        ligand_base = ligand.replace(".pdbqt", "")
        receptor_base = receptor.replace(".pdbqt", "")

        out_pdbqt = os.path.join(out_dir, f"{ligand_base}_{receptor_base}_docked.pdbqt")
        log_file = os.path.join(log_dir, f"{ligand_base}_{receptor_base}_log.txt")

        log(f"🔄 Docking {ligand_base} with {receptor_base}...")

        cmd = [
            "vina", "--receptor", os.path.join(target_dir, receptor),
            "--ligand", os.path.join(ligand_dir, ligand),
            "--center_x", x_center, "--center_y", y_center, "--center_z", z_center,
            "--size_x", x_size, "--size_y", y_size, "--size_z", z_size,
            "--num_modes", "9", "--out", out_pdbqt, "--log", log_file
        ]
        subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        best_score, best_pose_lines = extract_best_pose_from_pdbqt(out_pdbqt)
        if best_score is None:
            log(f"⚠️ No valid poses found in {out_pdbqt}")
            continue

        best_pose_pdbqt = os.path.join(best_pose_dir, f"{ligand_base}_{receptor_base}_best.pdbqt")
        with open(best_pose_pdbqt, "w") as f:
            f.writelines(best_pose_lines)

        best_pose_pdb = os.path.join(best_pose_dir, f"{ligand_base}_{receptor_base}_best.pdb")
        convert_to_pdb(best_pose_pdbqt, best_pose_pdb)

        # === Ensure receptor PDB file exists ===
        receptor_pdb = os.path.join(target_dir, f"{receptor_base}.pdb")
        if not os.path.exists(receptor_pdb):
            log(f"❌ Receptor PDB file not found: {receptor_pdb}")
            continue

        # === Ensure ligand PDB file exists ===
        if not os.path.exists(best_pose_pdb):
            log(f"❌ Ligand PDB file not found: {best_pose_pdb}")
            continue

        complex_path = os.path.join(complex_dir, f"{ligand_base}_{receptor_base}_complex.pdb")
        build_complex(best_pose_pdb, receptor_pdb, complex_path)

        log(f"✅ Best pose and complex for {ligand_base} with {receptor_base} saved (score: {best_score:.2f})")

log("🎯 Docking and complex generation completed for all ligand-target pairs.")
