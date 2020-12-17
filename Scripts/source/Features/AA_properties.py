hydropathy_kyte_doolitle = {
    'A': 1.8,   #alanine
    'C': 2.5,   #cysteine
    'D': -3.5,  #aspartic acid
    'E': -3.5,  #glumatic acid
    'F': 2.8,   #phenylalanine
    'G': -0.4,  #glycine
    'H': -3.2,  #histidine
    'I': 4.5,   #isoleucine
    'K': -3.9,  #lysine
    'L': 3.8,   #leucine
    'M': 1.9,   #methionine
    'N': -3.5,  #asparagine
    'P': -1.6,  #proline
    'Q': -3.5,  #glutamine
    'R': -4.5,  #arginine
    'S': -0.8,  #serine
    'T': -0.7,  #threonine
    'V': 4.2,   #valine
    'W': -0.9,  #tryptophan
    'Y': -1.3  #tyrosine
}

molecular_weight = {
    'A': 71.1,   #alanine
    'C': 103.1,   #cysteine
    'D': 115.1,  #aspartic acid
    'E': 129.1,  #glumatic acid
    'F': 147.2,   #phenylalanine
    'G': 57.0,  #glycine
    'H': 137.1,  #histidine
    'I': 113.2,   #isoleucine
    'K': 128.2,  #lysine
    'L': 113.2,   #leucine
    'M': 131.2,   #methionine
    'N': 114.1,  #asparagine
    'P': 97.1,  #proline
    'Q': 128.1,  #glutamine
    'R': 156.2,  #arginine
    'S': 87.1,  #serine
    'T': 101.1,  #threonine
    'V': 99.1,   #valine
    'W': 186.2,  #tryptophan
    'Y': 163.2   #tyrosine
}

polarity = {
    'A': "nonpolar",   #alanine
    'C': "polar_uncharged",   #cysteine
    'D': "polar",  #aspartic acid
    'E': "polar",  #glumatic acid
    'F': "nonpolar",   #phenylalanine
    'G': "nonpolar",  #glycine
    'H': "polar",  #histidine
    'I': "nonpolar",   #isoleucine
    'K': "polar",  #lysine
    'L': "nonpolar",   #leucine
    'M': "nonpolar",   #methionine
    'N': "polar_uncharged",  #asparagine
    'P': "nonpolar",  #proline
    'Q': "polar_uncharged",  #glutamine
    'R': "polar",  #arginine
    'S': "polar_uncharged",  #serine
    'T': "polar_uncharged",  #threonine
    'V': "nonpolar",   #valine
    'W': "nonpolar",  #tryptophan
    'Y': "polar_uncharged"   #tyrosine
}

charged = {
    'A': 0,   #alanine
    'C': 0,   #cysteine
    'D': 1,  #aspartic acid
    'E': 1,  #glumatic acid
    'F': 0,   #phenylalanine
    'G': 0,  #glycine
    'H': 1,  #histidine
    'I': 0,   #isoleucine
    'K': 1,  #lysine
    'L': 0,   #leucine
    'M': 0,   #methionine
    'N': 0,  #asparagine
    'P': 0,  #proline
    'Q': 0,  #glutamine
    'R': 1,  #arginine
    'S': 0,  #serine
    'T': 0,  #threonine
    'V': 0,   #valine
    'W': 0,  #tryptophan
    'Y': 0   #tyrosine
}

aromaticity = {
    'A': 0,   #alanine
    'C': 0,   #cysteine
    'D': 0,  #aspartic acid
    'E': 0,  #glumatic acid
    'F': 1,   #phenylalanine
    'G': 0,  #glycine
    'H': 1,  #histidine
    'I': 0,   #isoleucine
    'K': 0,  #lysine
    'L': 0,   #leucine
    'M': 0,   #methionine
    'N': 0,  #asparagine
    'P': 0,  #proline
    'Q': 0,  #glutamine
    'R': 0,  #arginine
    'S': 0,  #serine
    'T': 0,  #threonine
    'V': 0,   #valine
    'W': 1,  #tryptophan
    'Y': 1   #tyrosine
}

H_bond_atoms = {
    'A': 0,   #alanine
    'C': 0,   #cysteine
    'D': 4,  #aspartic acid
    'E': 4,  #glumatic acid
    'F': 0,   #phenylalanine
    'G': 0,  #glycine
    'H': 4,  #histidine
    'I': 0,   #isoleucine
    'K': 3,  #lysine
    'L': 0,   #leucine
    'M': 0,   #methionine
    'N': 4,  #asparagine
    'P': 0,  #proline
    'Q': 4,  #glutamine
    'R': 5,  #arginine
    'S': 3,  #serine
    'T': 3,  #threonine
    'V': 0,   #valine
    'W': 1,  #tryptophan
    'Y': 2   #tyrosine
}

def get_AA_scores(scores_dict, fasta):
    result = []
    for i in range(1, len(fasta) + 1):
        try:
            AA_score = scores_dict[fasta[i-1]] #unknown AA
        except:
            continue
        result.append((i, AA_score))
    return result
