
dataset_name="" #todo
new_name="" #todo
dataset_file="" #todo
output_file= "" #todo

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