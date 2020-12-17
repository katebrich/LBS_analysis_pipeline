import getopt
import sys
from os.path import isfile
import os

pdb_dir = ""
out_path = ""
list_ligands=False
dataset_path=""

#parse arguments:
try:
    opts, args = getopt.getopt(sys.argv[1:], 'p:o:d:l', ['pdb_dir=', 'output_path=', 'dataset_file=', 'list_ligands'] )
except getopt.GetoptError as err:
    print(err) #unknown option or missing argument
    sys.exit(1)
for opt, arg in opts:
    if opt in ("-p", "--pdb_dir"):
        pdb_dir = arg
    elif opt in ("-o", "--output_path"):
        out_path = arg
    elif opt in ("-d", "--dataset_file"):
        dataset_path = arg
    elif opt in ("-l", "--list_ligands"):
        list_ligands = True

if (pdb_dir == ""):
    print("Argument '-p' ('--pdb_dir') is compulsory.\n")
    sys.exit(1)
if (out_path == ""):
    print("Argument '-o' ('--output_path') is compulsory.\n")
    sys.exit(1)
if (list_ligands and dataset_path == ""):
    print("Argument '-d' ('--dataset_file') is compulsory with option '-l' (--list_ligands).\n")
    sys.exit(1)

out_dir = os.path.dirname(out_path) #go one level up

ligands_dict={}
if list_ligands:
    with open(dataset_path, 'r') as f:
        for line in f:
            columns = line.split()
            ligands_dict[columns[0]+columns[1]] = columns[2]

with open(out_path, 'w') as out_file:
    if list_ligands:
        out_file.write(f"HEADER: protein ligand_codes\n\n")
    for f in os.listdir(pdb_dir):
        f_path = os.path.join(pdb_dir, f)
        if isfile(f_path):
            relpath = os.path.relpath(f_path, out_dir)
            if list_ligands:
                ligands = ligands_dict[f[:5]]
                out_file.write(f"{relpath}\t{ligands}\n")
            else:
                out_file.write(f"{relpath}\n")



