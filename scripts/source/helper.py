#from __future__ import print_function
import os
import urllib.request
import json
import numpy as np
from urllib.error import HTTPError
from Bio.PDB import is_aa


#returns uniprotID, entityID, start, end
def get_uniprot_segments(pdb_id, chain_id):
    url = f'https://www.ebi.ac.uk/pdbe/api/mappings/uniprot_segments/{pdb_id}'

    try:
        parsedResponse = restAPI_get_json(url)
    except HTTPError as err:
        if err.code == 404:
            raise ValueError(f"HTTP Error 404 - Not Found: structure {pdb_id} {chain_id} probably does not exist in UniProt or there is a problem with connection to the server.")
        else:
            raise
    uniprotRecords = parsedResponse[f"{pdb_id}"]["UniProt"]

    segments = []
    for id in uniprotRecords:
        uniprot_ID = id
        mappings = uniprotRecords[uniprot_ID]["mappings"]
        for entity in mappings:
            if entity["chain_id"] == chain_id:
                start_res_num = entity["start"]["residue_number"]
                end_res_num = entity["end"]["residue_number"]
                unp_start = entity["unp_start"]
                unp_end = entity["unp_end"]
                if (unp_end - unp_start != end_res_num - start_res_num):
                    raise ValueError(
                        f"PDBe bug: Incorrect segments mapping: {pdb_id} {chain_id}, {uniprot_ID}: unp {unp_start}-{unp_end}, pdbe {start_res_num}-{end_res_num}")
                segments.append((uniprot_ID, unp_start, unp_end, start_res_num, end_res_num))

    if len(segments) < 1:
        raise ValueError(f"Structure {pdb_id} {chain_id} does not exist in UniProt. No uniprot segments were found. ")
    else:
        return segments

#returns parsed .json response (dictionary)
def restAPI_get(url):
    req = urllib.request.Request(url)
    with urllib.request.urlopen(req) as f:
        responseBody = f.read()
    return responseBody

#returns parsed .json response (dictionary)
def restAPI_get_json(url):
    return json.loads(restAPI_get(url))


def get_fasta_path_long(data_dir, pdb_id, chain_id):
    return f"{data_dir}/FASTA/{pdb_id}{chain_id}.fasta"

def get_pdb_path_long(data_dir, pdb_id, chain_id):
    return f"{data_dir}/PDB/{pdb_id}{chain_id}.pdb"

def get_pdb_path(data_dir, pdb_id, chain_id):
    return f"{data_dir}/{pdb_id}{chain_id}.pdb"

def get_lbs_path(data_dir, pdb_id, chain_id):
    return f"{data_dir}/{pdb_id}{chain_id}.txt"

def get_sasa_path_long(data_dir, pdb_id, chain_id):
    return f"{data_dir}/sasa/{pdb_id}{chain_id}.txt"

def get_mappings_path_long(data_dir, pdb_id, chain_id):
    return f"{data_dir}/mappings/{pdb_id}{chain_id}.txt"

def get_mappings_path(data_dir, pdb_id, chain_id):
    return f"{data_dir}/{pdb_id}{chain_id}.txt"

def get_feature_path_long(data_dir, feature, pdb_id, chain_id):
    return f"{data_dir}/features/{feature}/{pdb_id}{chain_id}.txt"


def get_entity_id(pdb_id, chain_id):
    url = f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/{pdb_id}"
    response = restAPI_get_json(url)
    molecules = response[pdb_id]
    entity_ids = []
    for molecule in molecules:
        if "polypeptide" in molecule["molecule_type"]:
            if chain_id in molecule["in_chains"]:
                entity_ids.append(molecule["entity_id"])
    if len(entity_ids) != 1:
        raise ValueError(f"{pdb_id} {chain_id}: Wrong number of entity IDs! {len(entity_ids)}")
    return int(entity_ids[0])

def res_mappings_author_to_pdbe(pdb_id, chain_id, cache_file=""):
    if (cache_file != ""):
        mappings = np.genfromtxt(cache_file, delimiter=' ', dtype='str')
        return list(mappings)
    else:
        mappings = []
        response = restAPI_get_json(f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/residue_listing/{pdb_id}/chain/{chain_id}")
        entity_id = get_entity_id(pdb_id, chain_id)
        molecules = response[pdb_id]["molecules"]
        count = 0
        for molecule in molecules:
            if molecule["entity_id"] == entity_id:
                for residue in molecule["chains"][0]["residues"]:
                    key = str(residue["author_residue_number"]) + residue["author_insertion_code"] #author residue number
                    val = residue["residue_number"] #pdbe molecule residue number
                    mappings.append((key, val))
                count += 1
        if count != 1:
            raise ValueError(f"Error: More or less than one molecule with entity number {entity_id} was found.")
        return list(mappings)


def parse_dataset_not_split_chains(filepath):
    list = []
    with open(filepath) as f:
        for line in f:
            line = line.split()
            pdb_id = line[0]
            chain_ids = line[1]
            chain_ids = chain_ids.upper()
            list.append((pdb_id, chain_ids))
    return list

def parse_dataset(filepath):
    list = []
    with open(filepath) as f:
        for line in f:
            line = line.split()
            pdb_id = line[0]
            chain_ids = line[1]
            ligands = []
            if len(line) == 3 and not line[2].startswith('#'):  # ligands are specified in dataset file
                ligands = line[2].split(',')
            chain_ids = chain_ids.upper()
            for chain in chain_ids.split(','):
                list.append((pdb_id, chain, ligands))
    return list

def parse_prank_dataset(filepath):
    pdb_files = []
    with open(filepath) as f:
        for line in f:
            if (line[0] != '#'): #skip comments
                words = line.split()
                if (len(words) > 0 and words[0] != 'HEADER:'): # skip empty lines and header
                    pdb_files.append(words[0])
    return pdb_files

def getFullAuthorResNum(residue_id):
    auth_res_num = str(residue_id[1])
    if not (residue_id[2].isspace()): # biopython returns a space instead of empty string
        auth_res_num += str(residue_id[2])  # insertion code
    return auth_res_num

def isPartOfChain(residue, mappings):
    if (residue.id[0].isspace()):  # hetero flag is empty
        return True
    elif (is_aa(residue, standard=False)):  # HETATM, but nonstandard AA code (MSE, LYZ etc.)
        # check if mapping exist for this residue number:
        auth_res_num = getFullAuthorResNum(residue.id)
        if (auth_res_num in mappings):
            return True # is part of the chain, not ligand
        else:
            return False  # is truly a ligand
    else:
        return False

def getStructuresFromDirectory(dir):
    list = []
    for filename in os.listdir(dir):
        pdb_id = filename[:4]
        chain_id = filename[4:5]
        list.append((pdb_id, chain_id))
    return list
