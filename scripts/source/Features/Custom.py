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