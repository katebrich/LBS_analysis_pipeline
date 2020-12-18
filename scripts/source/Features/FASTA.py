from Bio import SeqIO
from helper import *
from Features.AA_properties import *

def get_AA_properties(scores_dict, data_dir, pdb_id, chain_id):
    fasta_file = get_fasta_path_long(data_dir, pdb_id, chain_id)
    seq = list(SeqIO.parse(fasta_file, "fasta"))[0]
    return get_AA_scores(scores_dict, seq)

class Hydropathy():
    def get_values(self, data_dir, pdb_id, chain_id):
        return get_AA_properties(hydropathy_kyte_doolitle, data_dir, pdb_id, chain_id)

class Polarity():
    def get_values(self, data_dir, pdb_id, chain_id):
        return get_AA_properties(polarity, data_dir, pdb_id, chain_id)

class Charged():
    def get_values(self, data_dir, pdb_id, chain_id):
        return get_AA_properties(charged, data_dir, pdb_id, chain_id)

class Aromaticity():
    def get_values(self, data_dir, pdb_id, chain_id):
        return get_AA_properties(aromaticity, data_dir, pdb_id, chain_id)

class MolecularWeight():
    def get_values(self, data_dir, pdb_id, chain_id):
        return get_AA_properties(molecular_weight, data_dir, pdb_id, chain_id)

class HBondAtoms():
    def get_values(self, data_dir, pdb_id, chain_id):
        return get_AA_properties(H_bond_atoms, data_dir, pdb_id, chain_id)
