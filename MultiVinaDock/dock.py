import subprocess
import re
import os
import pandas as pd

VINA = 'vina'
OBABEL = 'obabel'

def extract_best_conformation(output_file):
    with open(output_file, 'r') as f:
        lines = f.readlines()
    best_score = float('inf')
    best_model = []
    current_model = []
    recording = False
    for line in lines:
        if line.startswith("MODEL"):
            recording = True
            current_model = []
        elif "ENDMDL" in line:
            recording = False
            score_line = [l for l in current_model if "REMARK VINA RESULT" in l]
            if score_line:
                score = float(score_line[0].split()[3])
                if score < best_score:
                    best_score = score
                    best_model = current_model
        if recording:
            current_model.append(line)
    return best_model, best_score

def run_vina_and_extract_best(task_dir, center_x, center_y, center_z, size_x, size_y, size_z, ligand_pdbqt_path, receptor_pdbqt_path, repeat_times):

    log_dir = task_dir+'/log'
    output_dir = task_dir+'/output'
    best_pdbqt_dir = task_dir+'/best_pdbqt'
    best_mol2_dir = task_dir+'/best_mol2'
    for dir in [task_dir, log_dir, output_dir, best_pdbqt_dir, best_mol2_dir]:
        if not os.path.exists(dir):
            os.mkdir(dir)

    best_scores = []
    for i in range(repeat_times):
        log_file = f"{log_dir}/run_{i + 1}.log"
        output_file = f"{output_dir}/run_{i + 1}.pdbqt"
        subprocess.run(f"{VINA} --center_x {center_x} --center_y {center_y} --center_z {center_z} --size_x {size_x} --size_y {size_y} --size_z {size_z} --ligand {ligand_pdbqt_path} --receptor {receptor_pdbqt_path} --log {log_file} --out {output_file}", shell=True)
        best_model, best_score = extract_best_conformation(output_file)
        best_scores.append(best_score)
        best_pdbqt_file = f"{best_pdbqt_dir}/run_{i + 1}_best.pdbqt"
        with open(best_pdbqt_file, 'w') as f:
            f.writelines(best_model)
        best_mol2_file = f"{best_mol2_dir}/run_{i + 1}_best.mol2"
        cmd = [OBABEL, '-ipdbqt', best_pdbqt_file, '-omol2', '-O', best_mol2_file]
        subprocess.run(cmd, check=True)

    with open(f"{task_dir}/best_scores.txt", 'w') as score_file:
        for index, score in enumerate(best_scores, 1):
            score_file.write(f"Run {index}: {score}\n")

if __name__ == '__main__':
    # run_vina_and_extract_best(task_name='dpp4', vina_config_file='dpp4/config.txt', repeat_times=100, )
    def task_batch_dock(task_name):
        df = pd.read_csv(f'task/{task_name}/{task_name}.csv', sep='\t')
        for index,row in df.iterrows():
            code = row['Protein Name']
            os.makedirs(f'task/{task_name}/dock/{code}', exist_ok=True)
            run_vina_and_extract_best(task_dir=f'task/{task_name}/dock/{code}', 
                                    size_x = row['size_x'],
                                    size_y = row['size_y'],
                                    size_z = row['size_z'],
                                    center_x = row['center_x'],
                                    center_y = row['center_y'],
                                    center_z = row['center_z'],
                                    ligand_pdbqt_path=row['Ligand PDBQT Path'], 
                                    receptor_pdbqt_path=row['Protein PDBQT Path'],
                                    repeat_times=100,)
    task_batch_dock('rec-lig')

