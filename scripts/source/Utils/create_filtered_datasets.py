import time
import numpy
from Bio.PDB import PDBParser
from scipy.spatial import distance
from Utils.MOAD import MOAD
from helper import parse_dataset, get_pdb_path, res_mappings_author_to_pdbe, isPartOfChain
from multiprocessing import Pool, Value

class DatasetLigandsFilter:
    filter = ""
    pdb_dir = ""

    def __init__(self, filter):
        self.filter = filter

    def run(self, orig_dataset_path, filtered_dataset_path, pdb_dir, threads=1, remove_empty_lines=False):
        self.pdb_dir = pdb_dir
        start = time.time()
        print(f"Filtering ligands in dataset {orig_dataset_path} with filter '{self.filter}' started...")

        if (self.filter == "MOAD"):
            self.moad = MOAD()

        dataset = parse_dataset(orig_dataset_path)

        self.total = len(dataset)

        counter = Value('i', 1)
        pool = Pool(int(threads), initializer=self.__pool_init, initargs=(counter,))
        results = pool.map(self.filter_ligands, dataset)
        pool.close()
        #total_results = [ent for sublist in results for ent in sublist]

        results_filtered = []
        if (remove_empty_lines): #remove structures with no valid ligands
            for res in results:
                if res is None:
                    continue
                if len(res[2]) > 0:
                    results_filtered.append(res)
        else:
            results_filtered = results


        with open(filtered_dataset_path, 'w') as f:
            f.write('\n'.join('{}\t{}\t{}'.format(x[0], x[1], ','.join(x[2])) for x in results_filtered))


        print(f"Finished in {time.time() - start}")

    def filter_ligands(self, structure):
        try:
            pdb_id = structure[0]
            chain_id = structure[1]
            pdb_file = get_pdb_path(self.pdb_dir, pdb_id, chain_id)
            mappings = dict(
                res_mappings_author_to_pdbe(pdb_id, chain_id))

            ligands_all = []
            AAs = []

            parser = PDBParser(PERMISSIVE=0, QUIET=1)
            structure = parser.get_structure(pdb_id + chain_id, pdb_file)
            chain = structure[0][chain_id]
            # get ligands
            for residue in chain.get_residues():
                if not isPartOfChain(residue, mappings):
                    # trim spaces at the beginning of ligand resname; bug in PDB parser
                    residue.resname = residue.resname.lstrip()
                    ligands_all.append(residue)
                else:
                    AAs.append(residue)

            valid_ligands = self.get_valid_ligands(ligands_all, AAs, pdb_id, chain_id)

            global counter
            with counter.get_lock():
                idx = counter.value
                counter.value += 1
            #print(f"{idx}/{self.total}: {pdb_id} {chain_id} processed")

            return (pdb_id, chain_id, valid_ligands)
        except:
            print(f"error {structure[0]} {structure[1]}")
            pass



    def get_valid_ligands(self, ligands, AAs, pdb_id, chain_id):
        result = []
        skipped = []
        small = []
        center_far = []
        far = []
        if (self.filter == ""):
            for ligand in ligands:
                result.append(ligand.resname)
        elif (self.filter == "p2rank"):
            for ligand in ligands:
                #name of the PDB group is not on the list of ignored groups:
                ignored = ["HOH", "DOD", "WAT", "NAG", "MAN", "UNK", "GLC", "ABA", "MPD", "GOL", "SO4", "PO4"]
                if (ligand.resname in ignored):
                    skipped.append(ligand.resname)
                    continue
                #number of ligand atoms is greater than or equal to 5:
                if (len(ligand.child_list) < 5):
                    small.append(ligand.resname)
                    continue
                #distance form the center of the mass of the ligand to the closest protein atom is not greater than 5.5 A
                center = self.__getCenterOfMass(ligand.child_list)
                threshold = 5.5
                success = False
                i = 0
                for AA in AAs:
                    for atom in AA.child_list:
                       if distance.euclidean(atom.coord, center) <= threshold:
                            success = True
                            i += 1
                            break
                if (success == False):
                    center_far.append(ligand.resname)
                    continue
                #distance from any atom of the ligand to the closest protein atom is at least 4 A
                success = False
                for AA in AAs:
                    if (self.__isInDistance(4, AA, ligand)):
                        success = True
                        break
                if (success == False): #no atom in distance 4A was found
                    far.append(ligand.resname)
                    continue

                result.append(ligand.resname)
            #print(
             #   f"{pdb_id} {chain_id}: Remaining {len(result)}/{len(ligands)} - {result}...Center far: {center_far}, Small: {small}, Ignored: {skipped}, Far: {far}")
        elif (self.filter == "MOAD"):
            relevant = self.moad.get_relevant_ligands(pdb_id, chain_id)
            if (relevant == None):
                print(f"{pdb_id} {chain_id} excluded from MOAD - no valid ligands or not an x-ray structure or has resolution higher than 2.5")
                return result #empty list
            filtered = []
            for ligand in ligands:
                if (str(ligand.resname) in relevant): # match only by resname, not id, because of the bug in MOAD (data not updated, old codes remain)
                    result.append(ligand.resname)
                else:
                    filtered.append(ligand.resname)
            #print( f"{pdb_id} {chain_id} RELEVANT : {relevant} ---- FILTERED : {filtered}")
        else:
            raise ValueError(f"Unknown filter level {self.filter}")

        result = list(set(result))  # distinct
        return result

    def __pool_init(self, args):
        ''' store the counter for later use '''
        global counter
        counter = args

    def __getCenterOfMass(self, atoms):
        totalMass = 0.0
        x = 0
        y = 0
        z = 0
        for a in atoms:
            m = a.mass
            totalMass += m
            x += a.coord[0] * m
            y += a.coord[1] * m
            z += a.coord[2] * m
        return (x / totalMass, y / totalMass, z / totalMass)

    def __isInDistance(self, threshold, AA_res, ligand_res):
        for atom1 in AA_res.child_list:
            if atom1.element == 'H': #consider only heavy atoms
                continue
            for atom2 in ligand_res.child_list:
                if atom2.element == 'H':
                    continue
                dist = numpy.linalg.norm(atom1.coord - atom2.coord)
                if dist < threshold:
                    return True
        return False


#######################################################################

filter="MOAD" #p2rank
dataset_file=f"" #TO EDIT
output_file= f"" #TO EDIT
pdb_dir=f"" #TO EDIT
threads=4


df = DatasetLigandsFilter(filter)
df.run(dataset_file, output_file, pdb_dir, threads, remove_empty_lines=True)