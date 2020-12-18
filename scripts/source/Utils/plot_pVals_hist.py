import os
import numpy as np

from AnalysisPipeline import Plots

dataset_folder = "/home/katebrich/Documents/diplomka/P2Rank/datasets_final/mix_filter_MOAD/analysis_1000"

iterations=1000

for feature in os.listdir(dataset_folder):
    try:
        feature_output_dir=os.path.join(dataset_folder, feature)
        p_values=list(np.genfromtxt(os.path.join(feature_output_dir, "p_values.txt"), delimiter=' ', dtype=None, encoding=None))
        Plots.plot_pvalues_histogram(p_values, iterations,
                                         os.path.join(feature_output_dir, f"{feature}_pValues_histogram.png"))
    except:
        print(f"ERROR: {feature}")