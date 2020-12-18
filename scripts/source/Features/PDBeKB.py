from helper import *

def get_funPDBe(resource, label, pdb_id, chain_id):
    url = f"https://www.ebi.ac.uk/pdbe/graph-api/pdb/funpdbe_annotation/{resource}/{pdb_id}"
    feature_vals = {}
    response = restAPI_get_json(url)
    bug = False #todo remove when bug in PDBe-KB is corrected
    for rec in response[pdb_id][0]["annotations"]:
        if rec["label"] == label:
            residues = rec["site_residues"]
            for res in residues:
                if (res["chain_id"] != chain_id):
                    continue
                score = res["raw_score"]
                res_num = res["residue_number"]

                if res_num in feature_vals:
                    bug = True
                feature_vals[res_num] = score

    if bug:
        raise ValueError(f"ERROR: {pdb_id} {chain_id}: {resource}: more values for some residues.")

    if (len(feature_vals) == 0): # no results for this chain
        raise ValueError(f"FunPDBe: {resource} {label}: no values found for {pdb_id} {chain_id}!")

    feature_vals_list = [(k, v) for k, v in feature_vals.items()]
    return feature_vals_list

class Conservation():
    def get_values(self, data_dir, pdb_id, chain_id):
        entity_id = get_entity_id(pdb_id, chain_id)
        url = f"https://www.ebi.ac.uk/pdbe/graph-api/pdb/sequence_conservation/{pdb_id}/{entity_id}"
        response = restAPI_get_json(url)
        feature_vals = []
        for resi in response[pdb_id]["data"]:
            feat_begin = int(resi["start"])
            cons_score = int(resi["conservation_score"])
            if (cons_score >= 4):
                cons_score = "4"
            feature_vals.append((feat_begin, cons_score))

        return feature_vals

class ConservationAllVals():
    def get_values(self, data_dir, pdb_id, chain_id):
        entity_id = get_entity_id(pdb_id, chain_id)
        url = f"https://www.ebi.ac.uk/pdbe/graph-api/pdb/sequence_conservation/{pdb_id}/{entity_id}"
        response = restAPI_get_json(url)

        feature_vals = []
        for resi in response[pdb_id]["data"]:
            feat_begin = int(resi["start"])
            cons_score = int(resi["conservation_score"])
            feature_vals.append((feat_begin, cons_score))

        return feature_vals

class Dynamine():
    def get_values(self, data_dir, pdb_id, chain_id):
        return get_funPDBe("dynamine", "backbone", pdb_id, chain_id)

class Efoldmine():
    def get_values(self, data_dir, pdb_id, chain_id):
        return get_funPDBe("dynamine", "efoldmine", pdb_id, chain_id)

class Depth():
    def get_values(self, data_dir, pdb_id, chain_id):
        return get_funPDBe("depth", "complex_residue_depth", pdb_id, chain_id)