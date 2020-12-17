from Bio import SeqIO
from helper import get_fasta_path_long

class AA():
    def get_values(self, data_dir, pdb_id, chain_id):
        fasta_file = get_fasta_path_long(data_dir, pdb_id, chain_id)
        seq = list(SeqIO.parse(fasta_file, "fasta"))[0]
        feature_vals = []
        for i in range(1, len(seq) + 1):
            AA = seq[i - 1]
            if AA in ['B','J','O','U','X','Z']:
                continue
            feature_vals.append((i, AA))
        return feature_vals

class AA_ALA():
    def get_values(self, data_dir, pdb_id, chain_id):
        fasta_file = get_fasta_path_long(data_dir, pdb_id, chain_id)
        seq = list(SeqIO.parse(fasta_file, "fasta"))[0]
        feature_vals = []
        for i in range(1, len(seq) + 1):
            AA = seq[i - 1]
            if AA == 'A':
                feature_vals.append((i, 1))
            else:
                feature_vals.append((i, 0))
        return feature_vals

class AA_CYS():
    def get_values(self, data_dir, pdb_id, chain_id):
        fasta_file = get_fasta_path_long(data_dir, pdb_id, chain_id)
        seq = list(SeqIO.parse(fasta_file, "fasta"))[0]
        feature_vals = []
        for i in range(1, len(seq) + 1):
            AA = seq[i - 1]
            if AA == 'C':
                feature_vals.append((i, 1))
            else:
                feature_vals.append((i, 0))
        return feature_vals

class AA_ASP():
    def get_values(self, data_dir, pdb_id, chain_id):
        fasta_file = get_fasta_path_long(data_dir, pdb_id, chain_id)
        seq = list(SeqIO.parse(fasta_file, "fasta"))[0]
        feature_vals = []
        for i in range(1, len(seq) + 1):
            AA = seq[i - 1]
            if AA == 'D':
                feature_vals.append((i, 1))
            else:
                feature_vals.append((i, 0))
        return feature_vals

class AA_GLU():
    def get_values(self, data_dir, pdb_id, chain_id):
        fasta_file = get_fasta_path_long(data_dir, pdb_id, chain_id)
        seq = list(SeqIO.parse(fasta_file, "fasta"))[0]
        feature_vals = []
        for i in range(1, len(seq) + 1):
            AA = seq[i - 1]
            if AA == 'E':
                feature_vals.append((i, 1))
            else:
                feature_vals.append((i, 0))
        return feature_vals

class AA_PHE():
    def get_values(self, data_dir, pdb_id, chain_id):
        fasta_file = get_fasta_path_long(data_dir, pdb_id, chain_id)
        seq = list(SeqIO.parse(fasta_file, "fasta"))[0]
        feature_vals = []
        for i in range(1, len(seq) + 1):
            AA = seq[i - 1]
            if AA == 'F':
                feature_vals.append((i, 1))
            else:
                feature_vals.append((i, 0))
        return feature_vals

class AA_GLY():
    def get_values(self, data_dir, pdb_id, chain_id):
        fasta_file = get_fasta_path_long(data_dir, pdb_id, chain_id)
        seq = list(SeqIO.parse(fasta_file, "fasta"))[0]
        feature_vals = []
        for i in range(1, len(seq) + 1):
            AA = seq[i - 1]
            if AA == 'G':
                feature_vals.append((i, 1))
            else:
                feature_vals.append((i, 0))
        return feature_vals

class AA_HIS():
    def get_values(self, data_dir, pdb_id, chain_id):
        fasta_file = get_fasta_path_long(data_dir, pdb_id, chain_id)
        seq = list(SeqIO.parse(fasta_file, "fasta"))[0]
        feature_vals = []
        for i in range(1, len(seq) + 1):
            AA = seq[i - 1]
            if AA == 'H':
                feature_vals.append((i, 1))
            else:
                feature_vals.append((i, 0))
        return feature_vals

class AA_ILE():
    def get_values(self, data_dir, pdb_id, chain_id):
        fasta_file = get_fasta_path_long(data_dir, pdb_id, chain_id)
        seq = list(SeqIO.parse(fasta_file, "fasta"))[0]
        feature_vals = []
        for i in range(1, len(seq) + 1):
            AA = seq[i - 1]
            if AA == 'I':
                feature_vals.append((i, 1))
            else:
                feature_vals.append((i, 0))
        return feature_vals

class AA_LYS():
    def get_values(self, data_dir, pdb_id, chain_id):
        fasta_file = get_fasta_path_long(data_dir, pdb_id, chain_id)
        seq = list(SeqIO.parse(fasta_file, "fasta"))[0]
        feature_vals = []
        for i in range(1, len(seq) + 1):
            AA = seq[i - 1]
            if AA == 'K':
                feature_vals.append((i, 1))
            else:
                feature_vals.append((i, 0))
        return feature_vals

class AA_LEU():
    def get_values(self, data_dir, pdb_id, chain_id):
        fasta_file = get_fasta_path_long(data_dir, pdb_id, chain_id)
        seq = list(SeqIO.parse(fasta_file, "fasta"))[0]
        feature_vals = []
        for i in range(1, len(seq) + 1):
            AA = seq[i - 1]
            if AA == 'L':
                feature_vals.append((i, 1))
            else:
                feature_vals.append((i, 0))
        return feature_vals

class AA_MET():
    def get_values(self, data_dir, pdb_id, chain_id):
        fasta_file = get_fasta_path_long(data_dir, pdb_id, chain_id)
        seq = list(SeqIO.parse(fasta_file, "fasta"))[0]
        feature_vals = []
        for i in range(1, len(seq) + 1):
            AA = seq[i - 1]
            if AA == 'M':
                feature_vals.append((i, 1))
            else:
                feature_vals.append((i, 0))
        return feature_vals

class AA_ASN():
    def get_values(self, data_dir, pdb_id, chain_id):
        fasta_file = get_fasta_path_long(data_dir, pdb_id, chain_id)
        seq = list(SeqIO.parse(fasta_file, "fasta"))[0]
        feature_vals = []
        for i in range(1, len(seq) + 1):
            AA = seq[i - 1]
            if AA == 'N':
                feature_vals.append((i, 1))
            else:
                feature_vals.append((i, 0))
        return feature_vals

class AA_PRO():
    def get_values(self, data_dir, pdb_id, chain_id):
        fasta_file = get_fasta_path_long(data_dir, pdb_id, chain_id)
        seq = list(SeqIO.parse(fasta_file, "fasta"))[0]
        feature_vals = []
        for i in range(1, len(seq) + 1):
            AA = seq[i - 1]
            if AA == 'P':
                feature_vals.append((i, 1))
            else:
                feature_vals.append((i, 0))
        return feature_vals

class AA_GLN():
    def get_values(self, data_dir, pdb_id, chain_id):
        fasta_file = get_fasta_path_long(data_dir, pdb_id, chain_id)
        seq = list(SeqIO.parse(fasta_file, "fasta"))[0]
        feature_vals = []
        for i in range(1, len(seq) + 1):
            AA = seq[i - 1]
            if AA == 'Q':
                feature_vals.append((i, 1))
            else:
                feature_vals.append((i, 0))
        return feature_vals

class AA_ARG():
    def get_values(self, data_dir, pdb_id, chain_id):
        fasta_file = get_fasta_path_long(data_dir, pdb_id, chain_id)
        seq = list(SeqIO.parse(fasta_file, "fasta"))[0]
        feature_vals = []
        for i in range(1, len(seq) + 1):
            AA = seq[i - 1]
            if AA == 'R':
                feature_vals.append((i, 1))
            else:
                feature_vals.append((i, 0))
        return feature_vals

class AA_SER():
    def get_values(self, data_dir, pdb_id, chain_id):
        fasta_file = get_fasta_path_long(data_dir, pdb_id, chain_id)
        seq = list(SeqIO.parse(fasta_file, "fasta"))[0]
        feature_vals = []
        for i in range(1, len(seq) + 1):
            AA = seq[i - 1]
            if AA == 'S':
                feature_vals.append((i, 1))
            else:
                feature_vals.append((i, 0))
        return feature_vals

class AA_THR():
    def get_values(self, data_dir, pdb_id, chain_id):
        fasta_file = get_fasta_path_long(data_dir, pdb_id, chain_id)
        seq = list(SeqIO.parse(fasta_file, "fasta"))[0]
        feature_vals = []
        for i in range(1, len(seq) + 1):
            AA = seq[i - 1]
            if AA == 'T':
                feature_vals.append((i, 1))
            else:
                feature_vals.append((i, 0))
        return feature_vals

class AA_VAL():
    def get_values(self, data_dir, pdb_id, chain_id):
        fasta_file = get_fasta_path_long(data_dir, pdb_id, chain_id)
        seq = list(SeqIO.parse(fasta_file, "fasta"))[0]
        feature_vals = []
        for i in range(1, len(seq) + 1):
            AA = seq[i - 1]
            if AA == 'V':
                feature_vals.append((i, 1))
            else:
                feature_vals.append((i, 0))
        return feature_vals

class AA_TRP():
    def get_values(self, data_dir, pdb_id, chain_id):
        fasta_file = get_fasta_path_long(data_dir, pdb_id, chain_id)
        seq = list(SeqIO.parse(fasta_file, "fasta"))[0]
        feature_vals = []
        for i in range(1, len(seq) + 1):
            AA = seq[i - 1]
            if AA == 'W':
                feature_vals.append((i, 1))
            else:
                feature_vals.append((i, 0))
        return feature_vals

class AA_TYR():
    def get_values(self, data_dir, pdb_id, chain_id):
        fasta_file = get_fasta_path_long(data_dir, pdb_id, chain_id)
        seq = list(SeqIO.parse(fasta_file, "fasta"))[0]
        feature_vals = []
        for i in range(1, len(seq) + 1):
            AA = seq[i - 1]
            if AA == 'Y':
                feature_vals.append((i, 1))
            else:
                feature_vals.append((i, 0))
        return feature_vals
