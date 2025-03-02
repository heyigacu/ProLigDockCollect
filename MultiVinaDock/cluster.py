
import os
import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from scipy.spatial.distance import squareform
import matplotlib.pyplot as plt
import rdkit.Chem as Chem
import shutil
import pandas as pd

def calculate_rmsd_simple(mol1, mol2):
    conf1 = mol1.GetConformer()
    conf2 = mol2.GetConformer()
    coords1 = np.array([conf1.GetAtomPosition(i) for i in range(mol1.GetNumAtoms())])
    coords2 = np.array([conf2.GetAtomPosition(i) for i in range(mol2.GetNumAtoms())])
    diff = coords1 - coords2
    rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
    return rmsd

def calculate_rmsd_matrix(conformers):
    num_conformers = len(conformers)
    rmsd_matrix = np.zeros((num_conformers, num_conformers))
    for i in range(num_conformers):
        for j in range(i + 1, num_conformers):
            rmsd = calculate_rmsd_simple(conformers[i],conformers[j])
            rmsd_matrix[i, j] = rmsd_matrix[j, i] = rmsd
    return rmsd_matrix

def plot(rmsd_matrix_path, out_png_path):
    rmsd_mat = np.loadtxt(rmsd_matrix_path)
    Z = linkage(squareform(rmsd_mat), method='single')
    new_labels = [str(int(label) + 1) for label in range(rmsd_mat.shape[0])]
    plt.figure()
    dendrogram(Z, labels=new_labels)
    plt.title('Hierarchical Clustering of Conformers')
    plt.xlabel('Conformer index')
    plt.ylabel('RMSD (Å)')
    plt.tight_layout()
    plt.savefig(f'{os.path.dirname(rmsd_matrix_path)}/cluster.png', dpi=600)
    plt.close()


def calculate_rmsd_matrix_mol2(mol2_dir):
    conformers = [Chem.MolFromMol2File(f'{mol2_dir}/run_{i}_best.mol2', removeHs=False, sanitize=False) for i in range(1,101)]
    rmsd_matrix = calculate_rmsd_matrix(conformers)
    rmsd_matrix_path = f'{os.path.dirname(mol2_dir)}/rmsd_matrix.txt'
    rmsd_cluster_png_path = f'{os.path.dirname(mol2_dir)}/rmsd_cluster.png'
    np.savetxt(rmsd_matrix_path, rmsd_matrix)
    plot(rmsd_matrix_path, rmsd_cluster_png_path)
    return rmsd_matrix_path

def find_representative_by_mean_rmsd(rmsd_matrix_path):
    rmsd_mat = np.loadtxt(rmsd_matrix_path)
    mean_rmsd = np.mean(rmsd_mat, axis=1)
    representative_idx = np.argmin(mean_rmsd)
    shutil.copy(os.path.dirname(rmsd_matrix_path)+f'/best_mol2/run_{representative_idx+1}_best.mol2', 
                os.path.dirname(rmsd_matrix_path)+f'/center.mol2')
    return representative_idx



def cluster_and_find_representative(rmsd_matrix_path, cutoff=30):
    rmsd_mat = np.loadtxt(rmsd_matrix_path)
    Z = linkage(squareform(rmsd_mat), method='single')
    clusters = fcluster(Z, cutoff, criterion='distance') 
    unique_clusters = np.unique(clusters)
    max_cluster = max(unique_clusters, key=lambda x: np.sum(clusters == x))
    indices = np.where(clusters == max_cluster)[0]
    intra_cluster_rmsd = np.sum(rmsd_mat[indices][:, indices], axis=0)
    representative_idx = indices[np.argmin(intra_cluster_rmsd)]
    print(f"Cluster {max_cluster} is the largest with representative molecule {representative_idx+1}.")
    shutil.copy(os.path.dirname(rmsd_matrix_path)+f'/best_mol2/run_{representative_idx+1}_best.mol2', os.path.dirname(rmsd_matrix_path)+f'/run_{representative_idx+1}_best.mol2')
    return representative_idx


def cluster_and_find_representative_auto(rmsd_matrix_path):
    rmsd_mat = np.loadtxt(rmsd_matrix_path)
    Z = linkage(squareform(rmsd_mat), method='single')
    new_labels = [str(int(label) + 1) for label in range(rmsd_mat.shape[0])]
    dendro_info = dendrogram(Z, labels=new_labels)
    color_dict = {}
    for i, label in enumerate(dendro_info['ivl']):
        cluster_color = dendro_info['leaves_color_list'][i]
        color_dict[label] = cluster_color
    color_groups = {}
    for label, color in color_dict.items():
        if color not in color_groups:
            color_groups[color] = []
        color_groups[color].append(label)
    max_cluster_color = max(color_groups, key=lambda k: len(color_groups[k]))
    max_cluster = color_groups[max_cluster_color]
    indices = [int(i)-1 for i in max_cluster]
    intra_cluster_rmsd = np.sum(rmsd_mat[indices][:, indices], axis=0)
    representative_idx = indices[np.argmin(intra_cluster_rmsd)]
    print(f"Cluster {max_cluster} is the largest with representative molecule {representative_idx+1}.")
    shutil.copy(os.path.dirname(rmsd_matrix_path)+f'/best_mol2/run_{representative_idx+1}_best.mol2', os.path.dirname(rmsd_matrix_path)+f'/run_{representative_idx+1}_best.mol2')
    # os.system(f'rm '+ os.path.dirname(rmsd_matrix_path)+f'/center.mol2')
    return representative_idx, os.path.dirname(rmsd_matrix_path)+f'/run_{representative_idx+1}_best.mol2'


if __name__ == '__main__':
    # mol2_dir = 'CXCR4-triptolide/best_mol2'
    # rsmd_matrix_path = calculate_rmsd_matrix_mol2(mol2_dir)
    # # rsmd_matrix_path = 'CXCR4-triptolide/rmsd_matrix.txt'
    # cutoff = float(input(f"please see cluster.png in {os.path.dirname(mol2_dir)}, and input the cutoff for cluster: ")) # 0.5
    # cluster_and_find_representative(rsmd_matrix_path, cutoff=cutoff) 

    def task_batch_center(task_name):
        df = pd.read_csv(f'task/{task_name}/{task_name}.csv', sep='\t')
        score_ls = []
        path_ls = []
        rmsd_ls = []
        for index,row in df.iterrows():
            code = row['Protein Name']
            lig = row['Ligand Name']
            mol2_dir = f'task/{task_name}/dock/{code}/best_mol2'
            scores_file = f'task/{task_name}/dock/{code}/best_scores.txt'
            scores = []
            with open(scores_file, 'r') as f:
                lines = f.readlines()
                for line in lines:
                    scores.append(float(line.split(':')[1].strip())) 
            rsmd_matrix_path = calculate_rmsd_matrix_mol2(mol2_dir)
            # rsmd_matrix_path = f'task/{task_name}/dock/{code}/rmsd_matrix.txt'
            representative_idx, representative_path = cluster_and_find_representative_auto(rsmd_matrix_path) 
            score_ls.append(scores[representative_idx])  
            path_ls.append(representative_path)  
            mol1 = Chem.MolFromMol2File(representative_path, removeHs=False, sanitize=False) 
            OBABEL = 'obabel'
            cry_pdbqt = f"task/{task_name}/ligand_pdbqt/{lig}.pdbqt"
            cry_path = f'task/{task_name}/ligand_mol2/{lig}.mol2'
            os.system(f'{OBABEL} -ipdbqt {cry_pdbqt} -omol2 -O {cry_path}')
            mol2 = Chem.MolFromMol2File(cry_path, removeHs=False, sanitize=False)
            rmsd = calculate_rmsd_simple(mol1, mol2)
            rmsd_ls.append(rmsd)
        df['representative affinity'] = score_ls
        df['representative mol2'] = path_ls
        df['dock RMSD'] = rmsd_ls
        df.to_csv(f'task/{task_name}/{task_name}.csv', sep='\t', index=False)
        
    task_batch_center('rec-lig')

