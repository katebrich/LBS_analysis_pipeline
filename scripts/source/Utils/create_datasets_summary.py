import os
import numpy as np
from Bio.PDB import PDBParser

from helper import parse_dataset, get_pdb_path, isPartOfChain, res_mappings_author_to_pdbe, get_mappings_path

P2Rank_path="/home/katebrich/Documents/diplomka/P2Rank"
datasets="mix,mix_filter_p2rank,mix_filter_MOAD"

results = []

for dataset_name in datasets.split(','):
    PDB_dir = f"{P2Rank_path}/datasets/{dataset_name}/PDB"
    LBS_dir = f"{P2Rank_path}/datasets/{dataset_name}/lbs"
    mappings_dir = f"{P2Rank_path}/datasets/{dataset_name}/mappings"
    dataset_file = f"/home/katebrich/Documents/diplomka/LBS_analysis_pipeline/data/datasets/{dataset_name}.txt"

    dataset = parse_dataset(dataset_file)

    #count structures
    str_count = len(dataset)

    lig_count = 0
    binding_count = 0
    nonbinding_count = 0

    #find out if there is ligand filter
    filter = False
    for str in dataset:
        if len(str[2]) != 0:
            filter = True
            break

    #count ligands
    for str in dataset:
        if filter:
            lig_count += len(str[2])
        else:
            pdb_id = str[0]
            chain_id = str[1]
            mappings = dict(
                res_mappings_author_to_pdbe(pdb_id, chain_id, get_mappings_path(mappings_dir, pdb_id, chain_id)))
            parser = PDBParser(PERMISSIVE=0, QUIET=1)
            structure = parser.get_structure(pdb_id + chain_id, get_pdb_path(PDB_dir, pdb_id, chain_id))
            chain = structure[0][chain_id]
            for residue in chain.get_residues():
                if not isPartOfChain(residue, mappings):
                    lig_count += 1

    #count binding and nonbinding residues
    for str in dataset:
        pdb_id = str[0]
        chain_id = str[1]
        file = os.path.join(LBS_dir, f"{pdb_id}{chain_id}.txt")
        lbs = list(np.genfromtxt(file, delimiter=' ', dtype=None)[:,1])
        binding_count += lbs.count(1)
        nonbinding_count += lbs.count(0)

    average_lig = round(lig_count / str_count, 4)
    ratio = round(binding_count / nonbinding_count, 4)

    results.append((dataset_name, str_count, lig_count, average_lig, binding_count, nonbinding_count, ratio))


for r in results:
    print(','. join('{}'.format(x) for x in r))