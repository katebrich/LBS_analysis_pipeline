import math
import os
import random
from matplotlib import pyplot, transforms
import numpy as np
from scipy import stats

iterations = 100

feature = "bfactor" #todo
data_path = "/home/katebrich/Documents/diplomka/P2Rank/datasets_final/mix_filter_p2rank/analysis_0" #todo
pairs_path = os.path.join(data_path, feature, "pairs.txt")
pairs = np.genfromtxt(pairs_path, delimiter=' ', dtype=None, encoding=None)
pairs = list(pairs)

data_binding = [x[1] for x in pairs if x[0] == 1]
data_nonbinding = [x[1] for x in pairs if x[0] == 0]

res_sizes = []
res_p_vals = []

sample_size = 10
step = 10
while sample_size < 1000:#2/3 * len(data_binding):
    p_values = []
    for i in range(1, iterations + 1):
        sample_binding = random.sample(data_binding, sample_size)
        sample_nonbinding = random.sample(data_nonbinding, sample_size)
        out = stats.ttest_ind(sample_binding, sample_nonbinding, axis=0, equal_var=False,
                              nan_policy='raise')  # Welch's test
        p_values.append(float(out[1]))
    mean_p_value = np.mean(p_values)
    res_sizes.append(sample_size)
    res_p_vals.append(mean_p_value)
    print(f"Sample {sample_size} done")
    sample_size = sample_size + step

alpha = 0.05
pyplot.clf()
fig, ax = pyplot.subplots()
pyplot.scatter(res_sizes, res_p_vals)
pyplot.xlabel('Sample size')
pyplot.ylabel('P-value')
ax.axhline(y=alpha, color="red", linestyle='-') #draw horizontal line for alpha value
trans = transforms.blended_transform_factory(
ax.get_yticklabels()[0].get_transform(), ax.transData)
ax.text(0, alpha, "{:.2f}".format(alpha), color="red", transform=trans, ha="right", va="center")
pyplot.savefig(os.path.join(data_path, feature, "pValue_deflation.svg"))
pyplot.clf()
pyplot.close('all')


