import os
import random

import numpy as np
from Bio import SeqIO

from helper import get_fasta_path_long


class LBS():
    def get_values(self, data_dir, pdb_id, chain_id):
        # get ligand binding sites values
        file = os.path.join(data_dir, "lbs", f"{pdb_id}{chain_id}.txt")
        lbs = tuple(map(tuple, np.genfromtxt(file, delimiter=' ', dtype=None)))
        return lbs

class RandomCont():
    def get_values(self, data_dir, pdb_id, chain_id):
        fasta_file = get_fasta_path_long(data_dir, pdb_id, chain_id) # get fasta sequence to know the length of structure
        fasta = list(SeqIO.parse(fasta_file, "fasta"))[0]
        result = []
        for i in range(1, len(fasta) + 1):
            value = random.uniform(0,10)
            result.append((i, value))
        return result
class RandomBinary():
    def get_values(self, data_dir, pdb_id, chain_id):
        fasta_file = get_fasta_path_long(data_dir, pdb_id, chain_id) # get fasta sequence to know the length of structure
        fasta = list(SeqIO.parse(fasta_file, "fasta"))[0]
        result = []
        for i in range(1, len(fasta) + 1):
            value = random.randint(0,1)
            result.append((i, value))
        return result