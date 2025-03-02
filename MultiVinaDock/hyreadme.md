


## 1. install mgltools
define variable $mgltools_install_dir in pdb2pdbqt.py
```
mgltools_install_dir = '/home/hy/Softwares/MGLTools/'
```
parepare input file:
* have two columns including 'Ligand PDB Path' and 'Protein PDB Path'
* file must be named $task.csv and locate at task/$task/$task.csv

## 2. prepare
```
python prepare.py
```

## 3. dock
```
python dock.py
```


## 4. cluster 
```
python cluster.py
```

