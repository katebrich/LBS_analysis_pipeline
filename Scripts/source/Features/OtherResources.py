from helper import *

class MobiDB():
    def get_values(self, data_dir, pdb_id, chain_id):
        entities = get_uniprot_segments(pdb_id, chain_id)
        feature_vals = []
        for entity in entities:
            uniprot_id = entity[0]
            unp_start = entity[1]
            unp_end = entity[2]
            start_res_num = entity[3]
            end_res_num = entity[4]

            url = f"https://mobidb.bio.unipd.it/ws/{uniprot_id}/consensus"
            response = restAPI_get_json(url)

            pred = response["mobidb_consensus"]["disorder"]["predictors"]
            for p in pred:
                if (p["method"] == "mobidb_lite"):
                    vals = p["scores"]
                    break

            i = 0
            for res_num in range(start_res_num, end_res_num + 1):
                feature_vals.append((res_num, vals[unp_start + i - 1]))
                i += 1

        return feature_vals

class Conservation():
    def get_values(self, data_dir, pdb_id, chain_id):
        filepath = os.path.join(data_dir, "conservation", pdb_id + chain_id + ".json")
        feature_vals = []
        with open(filepath) as json_file:
            data = json.load(json_file)
            scores = data["conservation"]
            for i in range(0, len(scores)):
                score = scores[i]
                if (score < 0):
                    score = 0
                feature_vals.append((i + 1, score))
        return feature_vals

class PhiAngle():
    def get_values(self, data_dir, pdb_id, chain_id):
        entity_id = get_entity_id(pdb_id, chain_id)
        url = f"https://www.ebi.ac.uk/pdbe/api/validation/rama_sidechain_listing/entry/{pdb_id}"
        response = restAPI_get_json(url)
        molecules = response[pdb_id]["molecules"]
        count = 0
        feature_vals = []
        for molecule in molecules:
            if molecule["entity_id"] == entity_id:
                for chain in molecule["chains"]:
                    if (chain["chain_id"] == chain_id):
                        for resi in chain["models"][0]["residues"]:
                            val = resi["phi"]
                            if (val != None):
                                feature_vals.append((resi["residue_number"], val))
                        count += 1
        if count != 1:
            raise ValueError(f"Error: More or less than one molecule with entity number {entity_id} was found.")

        return feature_vals

class PsiAngle():
    def get_values(self, data_dir, pdb_id, chain_id):
        entity_id = get_entity_id(pdb_id, chain_id)
        url = f"https://www.ebi.ac.uk/pdbe/api/validation/rama_sidechain_listing/entry/{pdb_id}"
        response = restAPI_get_json(url)
        molecules = response[pdb_id]["molecules"]
        count = 0
        feature_vals = []
        for molecule in molecules:
            if molecule["entity_id"] == entity_id:
                for chain in molecule["chains"]:
                    if (chain["chain_id"] == chain_id):
                        for resi in chain["models"][0]["residues"]:
                            val = resi["psi"]
                            if (val != None):
                                feature_vals.append((resi["residue_number"], val))
                        count += 1
        if count != 1:
            raise ValueError(f"Error: More or less than one molecule with entity number {entity_id} was found.")

        return feature_vals

class CisPeptide():
    def get_values(self, data_dir, pdb_id, chain_id):
        entity_id = get_entity_id(pdb_id, chain_id)
        url = f"https://www.ebi.ac.uk/pdbe/api/validation/rama_sidechain_listing/entry/{pdb_id}"
        response = restAPI_get_json(url)
        molecules = response[pdb_id]["molecules"]
        feature_vals = []
        for molecule in molecules:
            if molecule["entity_id"] == entity_id:
                for chain in molecule["chains"]:
                    if (chain["chain_id"] == chain_id):
                        for resi in chain["models"][0]["residues"]:
                            val = resi["cis_peptide"]
                            if (val == 'Y'):
                                feature_vals.append((resi["residue_number"], 1))
                            else:
                                feature_vals.append((resi["residue_number"], 0))
        return feature_vals
