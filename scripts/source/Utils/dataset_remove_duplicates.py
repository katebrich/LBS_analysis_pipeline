
dataset_name="" #TO EDIT
new_name="" #TO EDIT
dataset_file="" #TO EDIT
output_file= "" #TO EDIT

structures = {}
with open(dataset_file, 'r') as input:
    with open(output_file, 'w') as out:
        for line in input:
            s = line.split()
            if (s[0] in structures):
                if (structures[s[0]][0] == s[1]):
                    continue
            else:
                structures[s[0]] = (s[1], s[2])

with open(output_file, 'w') as out:
    for key, value in structures.items():
        out.write(f"{key}\t{value[0]}\t{value[1]}\n")