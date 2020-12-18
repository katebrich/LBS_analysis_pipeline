import math
import uuid

from Bio.PDB import *

from AnalysisPipeline import SASA
from helper import *
import os
import Logger
import time
import numpy

logger = Logger.get_logger(os.path.basename(__file__))
counter = None

class LigandBindingSitesComputer():
    dataset_file = ""
    output_dir = ""
    mappings_dir=""
    pdb_dir= ""
    distance_threshold = ""
    filter=""
    moad=""
    total = ""
    def __init__(self, dataset_file, output_dir, mappings_dir, PDB_dir, distance_threshold=4):
        self.output_dir = output_dir
        self.mappings_dir = mappings_dir
        self.pdb_dir = PDB_dir
        self.distance_threshold = distance_threshold
        self.dataset_file = dataset_file

    def run(self, threads):
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

        dataset = parse_dataset(self.dataset_file)

        #find out if there is ligand filter
        self.filter = False
        for str in dataset:
            if len(str[2]) != 0:
                self.filter = True
                break

        start = time.time()
        logger.info(f"Computing ligand binding sites started...")

        self.total = len(dataset)

        from multiprocessing import Pool, Value
        counter = Value('i', 1)
        pool = Pool(int(threads), initializer=self.__pool_init, initargs=(counter,))
        errors = pool.map(self.get_ligand_binding_site, dataset)
        pool.close()
        total_errors = [ent for sublist in errors for ent in sublist]

        if (len(total_errors) == 0):
            logger.info(f"Computing ligand binding sites finished in {math.ceil(time.time() - start)}s. All structures processed successfully.")
        else:
            errors_format = '\n'.join('%s %s' % x for x in total_errors)
            logger.warning(
                f"Computing ligand binding sites finished in {math.ceil(time.time() - start)}s. {len(total_errors)}/{self.total} structures were not processed successfully: \n{errors_format}")

    def get_ligand_binding_site(self, structure):
        pdb_id = structure[0]
        chain_id = structure[1]
        errors = []
        error=False
        pdb_path = get_pdb_path(self.pdb_dir, pdb_id, chain_id)
        try:
            lbs = self.compute_ligand_binding_sites(structure, pdb_path)
            output_file = get_lbs_path(self.output_dir, pdb_id, chain_id)
            with open(output_file, 'w') as f:
                f.write('\n'.join('{} {}'.format(x[0], x[1]) for x in lbs))
        except (KeyboardInterrupt, SystemExit):
            raise
        except Exception as ex:
            error=True
            logger.debug(f"Error while processing {pdb_id} {chain_id}: {ex}", exc_info=True)
        finally:
            global counter
            with counter.get_lock():
                idx = counter.value
                counter.value += 1
            if (error):
                errors.append((structure[0], structure[1]))
                logger.error(f"{idx}/{self.total}: {pdb_id} {chain_id} NOT PROCESSED ! See log for more details.")
            #else:
                #pass
            #    logger.debug(f"{idx}/{self.total}: {pdb_id} {chain_id} processed")
            return errors

    # returns list of tuples:
    # [0] residue number corresponding to PDBe molecule
    # [1] feature value - 0=not binding residue, 1=binding residue
    def compute_ligand_binding_sites(self, structure, pdb_file_path):
        pdb_id = structure[0]
        chain_id = structure[1]

        mappings = dict(res_mappings_author_to_pdbe(pdb_id, chain_id, get_mappings_path(self.mappings_dir, pdb_id, chain_id)))

        AAs = []
        ligands = []
        ligands_valid = None
        if (self.filter): #valid ligands read directly from dataset file
            ligands_valid = structure[2]
        #get all ligands from PDB file
        parser = PDBParser(PERMISSIVE=0, QUIET=1)
        structure = parser.get_structure(pdb_id + chain_id, pdb_file_path)
        chain = structure[0][chain_id]
        for residue in chain.get_residues():
            if not isPartOfChain(residue, mappings):
                #trim spaces at the beginning of ligand resname; bug in PDB parser
                residue.resname = residue.resname.lstrip()
                if (self.filter):
                    if residue.resname in ligands_valid:
                        ligands.append(residue)
                else:
                    ligands.append(residue)

        # get AAs
        # create temporary .pdb file without HETATM lines to compute Solvent Accessible Surface atoms
        temp_file_path = os.path.join(self.output_dir, f"temp_{uuid.uuid1()}")
        try:
            with open(pdb_file_path, 'r') as orig:
                with open(temp_file_path, 'w') as temp:
                    for line in orig:
                        if line.startswith("ATOM"):
                            temp.write(line)
            parser = PDBParser(PERMISSIVE=0, QUIET=1)
            structure = parser.get_structure(pdb_id + chain_id, temp_file_path)
            chain = structure[0][chain_id]
            sr = SASA.ShrakeRupley()
            sr.compute(structure, level="R")
            for residue in chain.get_residues():
                if residue.resname not in sasa_residues: #for example for unknown amino acids (UNK, Xaa)
                    threshold = 5 #default threshold
                else:
                    threshold = 0.05 * sasa_residues[residue.resname]
                if (isPartOfChain(residue, mappings) and residue.sasa > threshold): #leave out residues that are not on the surface
                    AAs.append(residue)
        finally:
            if os.path.exists(temp_file_path):
                os.remove(temp_file_path)

        output = []
        for AA in AAs:
            success = False
            auth_res_num = getFullAuthorResNum(AA.id)
            pdbe_res_num = mappings[auth_res_num]
            for ligand in ligands:
                if (self.isInDistance(AA, ligand)):
                    output.append((pdbe_res_num, 1))
                    success = True
                    break
            if (success == False):
                output.append((pdbe_res_num, 0))
        return output

    def __pool_init(self, args):
        ''' store the counter for later use '''
        global counter
        counter = args

    #True if there exists a heavy (non-hydrogen) atom on SAS of AA residue which is at maximum distance "threshold" from at least one ligand residue heavy atom
    def isInDistance(self, AA_res, ligand_res):
        for atom1 in AA_res.child_list:
            if atom1.element == 'H' or atom1.sasa == 0:  # consider only heavy atoms on Solvent Accessible Surface (SAS)
                continue
            for atom2 in ligand_res.child_list:
                if atom2.element == 'H':
                    continue
                dist = numpy.linalg.norm(atom1.coord - atom2.coord)
                if dist < self.distance_threshold:
                    return True
        return False


sasa_residues = {
    'ALA': 113,   #alanine
    'CYS': 140,   #cysteine
    'ASP': 151,  #aspartic acid
    'GLU': 183,  #glumatic acid
    'PHE': 218,   #phenylalanine
    'GLY': 85,  #glycine
    'HIS': 194,  #histidine
    'ILE': 182,   #isoleucine
    'LYS': 211,  #lysine
    'LEU': 180,   #leucine
    'MET': 204,   #methionine
    'ASN': 158,  #asparagine
    'PRO': 143,  #proline
    'GLN': 189,  #glutamine
    'ARG': 241,  #arginine
    'SER': 122,  #serine
    'THR': 146,  #threonine
    'VAL': 160,   #valine
    'TRP': 259,  #tryptophan
    'TYR': 229,   #tyrosine
}