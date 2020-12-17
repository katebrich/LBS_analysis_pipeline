
from helper import *

def get_uniprot_binary_type(type, pdb_id, chain_id):
    segments = get_uniprot_segments(pdb_id, chain_id)
    feature_vals = []
    for seg in segments:
        uniprot_id = seg[0]
        segment_begin = seg[1]
        segment_end = seg[2]
        start_res_num = seg[3]
        end_res_num = seg[4]
        url = f"https://www.ebi.ac.uk/proteins/api/features/{uniprot_id}?types={type}"
        response = restAPI_get_json(url)

        feature_vector = [0] * (segment_end - segment_begin + 1)  # including both start and end AAs
        for feature in response["features"]:
            feat_begin = int(feature["begin"])
            feat_end = int(feature["end"])
            for i in range(max(feat_begin, segment_begin), min(feat_end, segment_end) + 1): #take only parts of features that overlap with uniprot segment
                res = i - segment_begin  # mapping pdb residues to uniprot entry, counting from 0
                feature_vector[res] = 1

        i = 0
        for res_num in range(start_res_num, end_res_num + 1):
            feature_vals.append((res_num, feature_vector[i]))
            i += 1

    return feature_vals

class PTM():
    def get_values(self, data_dir, pdb_id, chain_id):
        return get_uniprot_binary_type("CARBOHYD%2CMOD_RES%2CLIPID", pdb_id, chain_id)  # glycosylation+lipidation+mod_res

class Glycosylation():
    def get_values(self, data_dir, pdb_id, chain_id):
        return get_uniprot_binary_type("CARBOHYD", pdb_id, chain_id)

class Lipidation():
    def get_values(self, data_dir, pdb_id, chain_id):
        return get_uniprot_binary_type("LIPID", pdb_id, chain_id)

class ModRes():
    def get_values(self, data_dir, pdb_id, chain_id):
        return get_uniprot_binary_type("MOD_RES", pdb_id, chain_id)

class NonStandard():
    def get_values(self, data_dir, pdb_id, chain_id):
        return get_uniprot_binary_type("NON_STD", pdb_id, chain_id)

class NaturalVariant():
    def get_values(self, data_dir, pdb_id, chain_id):
        return get_uniprot_binary_type("VARIANT", pdb_id, chain_id)

class Compbias():
    def get_values(self, data_dir, pdb_id, chain_id):
        return get_uniprot_binary_type("COMPBIAS", pdb_id, chain_id)

class Disulfid():
    def get_values(self, data_dir, pdb_id, chain_id):
        segments = get_uniprot_segments(pdb_id, chain_id)
        feature_vals = []
        for seg in segments:
            uniprot_id = seg[0]
            segment_begin = seg[1]
            segment_end = seg[2]
            start_res_num = seg[3]
            end_res_num = seg[4]
            url = f"https://www.ebi.ac.uk/proteins/api/features/{uniprot_id}?types=DISULFID"
            response = restAPI_get_json(url)

            feature_vector = [0] * (segment_end - segment_begin + 1)  # including both start and end AAs
            for feature in response["features"]:
                feat_begin = int(feature["begin"])
                feat_end = int(feature["end"])
                for i in (feat_begin, feat_end):  # take only parts of features that overlap with uniprot segment
                    res = i - segment_begin  # mapping pdb residues to uniprot entry, counting from 0
                    if res >= 0 and res < segment_end - segment_begin + 1:
                        feature_vector[res] = 1

            i = 0
            for res_num in range(start_res_num, end_res_num + 1):
                feature_vals.append((res_num, feature_vector[i]))
                i += 1

        return feature_vals

class Variation():
    def get_values(self, data_dir, pdb_id, chain_id):
        segments = get_uniprot_segments(pdb_id, chain_id)
        feature_vals = []
        for seg in segments:
            uniprot_id = seg[0]
            segment_begin = seg[1]
            segment_end = seg[2]
            start_res_num = seg[3]
            end_res_num = seg[4]

            url = f"https://www.ebi.ac.uk/proteins/api/variation/{uniprot_id}"
            response = restAPI_get_json(url)

            feature_vector = [0] * (segment_end - segment_begin + 1)  # including both start and end AAs

            for feature in response["features"]:
                feat_begin = int(feature["begin"])
                feat_end = int(feature["end"])
                rng = range(feat_begin, feat_end + 1)
                for i in rng:
                    res = i - segment_begin  # mapping pdb residues to uniprot entry, counting from 0
                    if res >= 0 and res < segment_end - segment_begin + 1:
                        feature_vector[res] = 1
            i = 0
            for res_num in range(start_res_num, end_res_num + 1):
                feature_vals.append((res_num, feature_vector[i]))
                i += 1

        return feature_vals

class SecondaryStructure():
    def get_values(self, data_dir, pdb_id, chain_id):
        segments = get_uniprot_segments(pdb_id, chain_id)
        feature_vals = []
        for seg in segments:
            uniprot_id = seg[0]
            segment_begin = seg[1]
            segment_end = seg[2]
            start_res_num = seg[3]
            end_res_num = seg[4]
            url = f"https://www.ebi.ac.uk/proteins/api/features/{uniprot_id}?categories=STRUCTURAL"
            response = restAPI_get_json(url)

            feature_vector = ['X'] * (segment_end - segment_begin + 1)  # including both start and end AAs
            for feature in response["features"]:
                feat_begin = int(feature["begin"])
                feat_end = int(feature["end"])
                type = feature["type"]
                if type == "HELIX":
                    val = 'H'
                elif type == "STRAND":
                    val = 'S'
                elif type == "TURN":
                    val = 'T'
                else:
                    raise ValueError(f"Unknown feature type '{type}'")
                for i in range(max(feat_begin, segment_begin), min(feat_end, segment_end) + 1):
                    res = i - segment_begin  # mapping pdb residues to uniprot entry, counting from 0
                    feature_vector[res] = val

            i = 0
            for res_num in range(start_res_num, end_res_num + 1):
                feature_vals.append((res_num, feature_vector[i]))
                i += 1

        return feature_vals

class Helix():
    def get_values(self, data_dir, pdb_id, chain_id):
        segments = get_uniprot_segments(pdb_id, chain_id)
        feature_vals = []
        for seg in segments:
            uniprot_id = seg[0]
            segment_begin = seg[1]
            segment_end = seg[2]
            start_res_num = seg[3]
            end_res_num = seg[4]
            url = f"https://www.ebi.ac.uk/proteins/api/features/{uniprot_id}?categories=STRUCTURAL"
            response = restAPI_get_json(url)

            feature_vector = [0] * (segment_end - segment_begin + 1)  # including both start and end AAs
            for feature in response["features"]:
                feat_begin = int(feature["begin"])
                feat_end = int(feature["end"])
                type = feature["type"]
                if type == "HELIX":
                    val = 1
                else:
                    val = 0

                for i in range(max(feat_begin, segment_begin), min(feat_end, segment_end) + 1):
                    res = i - segment_begin  # mapping pdb residues to uniprot entry, counting from 0
                    feature_vector[res] = val

            i = 0
            for res_num in range(start_res_num, end_res_num + 1):
                feature_vals.append((res_num, feature_vector[i]))
                i += 1

        return feature_vals

class Strand():
    def get_values(self, data_dir, pdb_id, chain_id):
        segments = get_uniprot_segments(pdb_id, chain_id)
        feature_vals = []
        for seg in segments:
            uniprot_id = seg[0]
            segment_begin = seg[1]
            segment_end = seg[2]
            start_res_num = seg[3]
            end_res_num = seg[4]
            url = f"https://www.ebi.ac.uk/proteins/api/features/{uniprot_id}?categories=STRUCTURAL"
            response = restAPI_get_json(url)

            feature_vector = [0] * (segment_end - segment_begin + 1)  # including both start and end AAs
            for feature in response["features"]:
                feat_begin = int(feature["begin"])
                feat_end = int(feature["end"])
                type = feature["type"]
                if type == "STRAND":
                    val = 1
                else:
                    val = 0

                for i in range(max(feat_begin, segment_begin), min(feat_end, segment_end) + 1):
                    res = i - segment_begin  # mapping pdb residues to uniprot entry, counting from 0
                    feature_vector[res] = val

            i = 0
            for res_num in range(start_res_num, end_res_num + 1):
                feature_vals.append((res_num, feature_vector[i]))
                i += 1

        return feature_vals

class Turn():
    def get_values(self, data_dir, pdb_id, chain_id):
        segments = get_uniprot_segments(pdb_id, chain_id)
        feature_vals = []
        for seg in segments:
            uniprot_id = seg[0]
            segment_begin = seg[1]
            segment_end = seg[2]
            start_res_num = seg[3]
            end_res_num = seg[4]
            url = f"https://www.ebi.ac.uk/proteins/api/features/{uniprot_id}?categories=STRUCTURAL"
            response = restAPI_get_json(url)

            feature_vector = [0] * (segment_end - segment_begin + 1)  # including both start and end AAs
            for feature in response["features"]:
                feat_begin = int(feature["begin"])
                feat_end = int(feature["end"])
                type = feature["type"]
                if type == "TURN":
                    val = 1
                else:
                    val = 0

                for i in range(max(feat_begin, segment_begin), min(feat_end, segment_end) + 1):
                    res = i - segment_begin  # mapping pdb residues to uniprot entry, counting from 0
                    feature_vector[res] = val

            i = 0
            for res_num in range(start_res_num, end_res_num + 1):
                feature_vals.append((res_num, feature_vector[i]))
                i += 1

        return feature_vals