'''
Place for defining new features.

The class with new feature must implement
        get_values(self, data_dir, pdb_id, chain_id)
where data_dir is the top folder given by -o argument

Example:
    class NewFeature():
        def get_values(self, data_dir, pdb_id, chain_id):
            *** computation here ***
            return  computed_values

Returned feature values must be tuples where first item is residue number (PDB molecule numbering) and the second item is feature value for the residue.

Edit of config.json is needed to load the class! Add new feature to "features" node.

Example:
    "new_feature": {
			"import_path": "Features.Custom.NewFeature",
			"type": "binary",
			"default": "0"
	}

'''




import os
from helper import res_mappings_author_to_pdbe, get_mappings_path_long
class INTAAConservation():
    def get_values(self, data_dir, pdb_id, chain_id):
        filepath = os.path.join(data_dir, "INTAA_conservation", pdb_id + chain_id + ".pdb.ic")
        feature_vals = []
        mappings = dict(
            res_mappings_author_to_pdbe(pdb_id, chain_id, get_mappings_path_long(data_dir, pdb_id, chain_id)))
        with open(filepath) as file:
            for line in file:
                line = line.split('\t')
                score = line[4]
                auth_res_num = line[2]+line[3].strip() # author res num + insertion code
                pdbe_res_num = mappings[auth_res_num]
                feature_vals.append((pdbe_res_num, score))
        return feature_vals