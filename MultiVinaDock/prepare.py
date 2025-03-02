import os
import subprocess
import pandas as pd

mgltools_install_dir = '/home/hy/Softwares/MGLTools/'
# Prepare ligand and receptor scripts paths
prepare_ligand4_path = os.path.join(mgltools_install_dir, 'MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py')
prepare_receptor4_path = os.path.join(mgltools_install_dir, 'MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py')
pythonsh_path = os.path.join(mgltools_install_dir, 'bin/pythonsh')
parent_dir = os.path.abspath(os.path.dirname(__name__))

# Directory to save the PDBQT files
def generate_pdbqt(df, pdbqt_save_dir, is_ligand=True):
    ls = []
    for index, row in df.iterrows():
        pdb_path = row['Ligand PDB Path'] if is_ligand else row['Protein PDB Path']
        pdb_file = os.path.basename(pdb_path)
        pdb_dir = os.path.dirname(pdb_path)
        pdbqt_name = pdb_file.replace('.pdb', '.pdbqt')
        pdbqt_path = os.path.join(pdbqt_save_dir, pdbqt_name)
        ls.append(pdbqt_path)
        if is_ligand:
            command = f'cd {pdb_dir}; {pythonsh_path} {prepare_ligand4_path} -l {pdb_path} -o {pdbqt_path}'
        else:
            command = f'cd {pdb_dir}; {pythonsh_path} {prepare_receptor4_path} -r {pdb_path} -o {pdbqt_path}'
        try:
            subprocess.run(command, shell=True, check=True)
            print(f"Successfully converted {pdb_path} to {pdbqt_path}")
        except subprocess.CalledProcessError as e:
            print(f"Error processing {pdb_path}: {e}")
    if is_ligand:
        df['Ligand PDBQT Path'] = ls
    else:
        df['Protein PDBQT Path'] = ls



def batch_pdb2pdbqt(task_name):
    protein_pdbqt_save_dir = parent_dir+f'/task/{task_name}/protein_pdbqt'
    ligand_pdbqt_save_dir = parent_dir+f'/task/{task_name}/ligand_pdbqt'
    if not os.path.exists(protein_pdbqt_save_dir): os.mkdir(protein_pdbqt_save_dir)
    if not os.path.exists(ligand_pdbqt_save_dir): os.mkdir(ligand_pdbqt_save_dir)
    # Read the CSV file
    df = pd.read_csv(f'task/{task_name}/{task_name}.csv', sep='\t')
    # Generate PDBQT for ligands
    generate_pdbqt(df, ligand_pdbqt_save_dir, is_ligand=True)
    # Generate PDBQT for proteins
    generate_pdbqt(df, protein_pdbqt_save_dir, is_ligand=False)
    df.to_csv(f'task/{task_name}/{task_name}.csv', sep='\t', index=False)


if __name__ == "__main__":
    batch_pdb2pdbqt('rec-lig')


