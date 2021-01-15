import math
import os
import random
import shutil
import time
from multiprocessing import Pool
import numpy as np
from scipy import stats
from collections import Counter
from AnalysisPipeline import Plots
from helper import getStructuresFromDirectory

import Logger

logger = Logger.get_logger(os.path.basename(__file__))

class AnalysisComputer():
    output_dir = ""
    lbs_dir = ""
    features_dir = ""
    config = ""
    p_values = []
    b_ratios = []
    p_values_means = []
    p_val_perc = []
    cohens_d = []
    cohens_d_means = []
    cohens_w = []
    cohens_w_means = []
    means = []
    errors = []
    features_list = []
    feature_dir = ""
    feature = ""
    lbs_dicts = {}

    def __init__(self, output_dir, lbs_dir, features_dir, features_list, config):
        self.output_dir = output_dir
        self.lbs_dir = lbs_dir
        self.config = config
        self.features_list = features_list
        self.features_dir = features_dir
        random.seed(42)

    def run(self, sample_size, iterations, balance_binding_ratio, draw_plots, alpha, threads):
        start = time.time()
        logger.info(f"Running analysis started...")

        self.prepare_lbs_dicts()

        for feature in self.features_list:
            feature_output_dir = os.path.join(self.output_dir, feature)
            if os.path.exists(feature_output_dir):
                shutil.rmtree(feature_output_dir)
            os.makedirs(feature_output_dir)
            try:
                if self.process_feature(feature, sample_size, iterations, balance_binding_ratio, draw_plots, alpha, threads):
                    logger.info(f"Feature {feature}: done")
            except Exception as ex:
                logger.error(f"Feature {feature}: ERROR: {ex}", exc_info=True)
                self.errors.append(feature)
                continue

        logger.info(f"Running analysis finished in {math.ceil(time.time() - start)}s.")

    def process_feature(self, feature, sample_size, iterations, balance_binding_ratio, draw_plots, alpha, threads):
        feature_dir = os.path.join(self.features_dir, feature)
        dataset = getStructuresFromDirectory(
            feature_dir)  # compute only for structures that have the feature computed, ignore the rest
        feature_output_dir = os.path.join(self.output_dir, feature)
        feature_type = self.config.get_feature_type(feature) # binary/continuous/categorical/ordinal
        statistics = []
        p_values = []
        cohens_d = []
        cohens_w = []

        self.feature_dir = feature_dir
        self.feature = feature
        pool = Pool(int(threads))
        pairs_part = pool.map(self.compute_pairs, dataset)
        pool.close()
        pairs = [ent for sublist in pairs_part for ent in sublist]

        #save pairs
        file = os.path.join(feature_output_dir, f"pairs.txt")
        with open(file, 'w') as f:
            f.write('\n'.join('{} {}'.format(x[0], x[1]) for x in pairs))

        data_binding = [x[1] for x in pairs if x[0] == 1]
        data_nonbinding = [x[1] for x in pairs if x[0] == 0]

        if sample_size == 0:
            iterations = 1

        with open(os.path.join(feature_output_dir, f"iterations.txt"), 'w') as f:
            for i in range(1, iterations + 1):
                f.write(f"********************************\n")
                f.write(f"**        ITERATION {i}       **\n")
                f.write(f"********************************\n")

                if (sample_size == 0): #no sampling
                    if len(data_binding) < 10 or len(data_nonbinding) < 10:
                        logger.error(f"Feature {feature}: Not sufficient number of observations to perform the analysis.")
                        self.errors.append(feature)
                        return False
                    sample_binding = data_binding
                    sample_nonbinding = data_nonbinding
                elif not balance_binding_ratio: #do not balance the ratio, sample from the whole dataset
                    if len(pairs) < sample_size:
                        logger.error(f"Feature {feature}: The sample size is bigger than number of observations.")
                        self.errors.append(feature)
                        return False
                    sample = random.sample(pairs, sample_size)
                    sample_binding = [x[1] for x in sample if x[0] == 1]
                    sample_nonbinding = [x[1] for x in sample if x[0] == 0]
                else: #simple sampling without replacement
                    if len(data_binding) < sample_size or len(data_nonbinding) < sample_size:
                        logger.error(f"Feature {feature}: The sample size is bigger than number of observations.")
                        self.errors.append(feature)
                        return False
                    sample_binding = random.sample(data_binding, sample_size)
                    sample_nonbinding = random.sample(data_nonbinding, sample_size)

                #run analysis and save results to file
                if (feature_type == "continuous"):
                    res = self.welchs_t_test(sample_binding, sample_nonbinding, f)
                    if res == None:
                        self.errors.append(feature)
                        return False
                    cohens_d.append(self.get_cohens_d(sample_binding, sample_nonbinding))
                elif (feature_type == "categorical" or feature_type == "ordinal" or feature_type == "binary"):
                    res = self.chi_squared_test(sample_binding, sample_nonbinding, f, feature)
                    if res == None:
                        self.errors.append(feature)
                        return False
                    cohens_w.append(self.get_cohens_w(res[0], len(sample_binding) + len(sample_nonbinding)))
                else:
                    logger.error(f"Unknown type of feature '{feature_type}'. Please specify the type in config.")
                    self.errors.append(feature)
                    return False

                statistics.append(res[0])
                p_values.append(res[1])

        if draw_plots:
            self.draw_plots(feature_type, data_binding, data_nonbinding, feature_output_dir, feature)

        #save values for summary
        self.p_values.append((feature, p_values))
        self.b_ratios.append((feature, len(data_binding) / (len(data_binding) + len(data_nonbinding))))
        self.p_values_means.append((feature, (np.mean(p_values), np.std(p_values))))
        if (feature_type == "continuous"):
            self.cohens_d.append((feature, cohens_d))
            self.cohens_d_means.append((feature, (np.mean(cohens_d), np.std(cohens_d))))
        else:
            self.cohens_w.append((feature, cohens_w))
            self.cohens_w_means.append((feature, (np.mean(cohens_w), np.std(cohens_w))))

        #todo
        ##calculate effect sizes
        #if (feature_type == "continuous"):
        #    self.cohens_d_vals.append((feature, self.cohens_d(data_binding, data_nonbinding)))
        #else:
        #    self.cohens_w_vals.append((feature, self.cohens_w(data_binding, data_nonbinding)))

        #save difference of means and variance for summary
        if (feature_type == "continuous"):
            mean_binding = np.mean(data_binding)
            mean_nonbinding = np.mean(data_nonbinding)
            means_dif =  abs(mean_binding - mean_nonbinding)
            variance_binding = np.var(data_binding)
            variance_nonbinding = np.var(data_nonbinding)
            self.means.append((feature, (means_dif, mean_binding, mean_nonbinding, variance_binding, variance_nonbinding)))

        #calculate percentage of iterations for which p-value was under significance level alpha
        under_alpha = 0
        for p in p_values:
            if p <= alpha:
                under_alpha = under_alpha + 1
        perc = (under_alpha / iterations) * 100
        self.p_val_perc.append((feature, round(perc,2)))

        #plot p-values
        Plots.plot_pvalues_scatter(p_values, alpha, os.path.join(feature_output_dir, f"{feature}_pValues_scatter.png"))
        Plots.plot_pvalues_histogram(p_values, iterations,
                                     os.path.join(feature_output_dir, f"{feature}_pValues_histogram.png"))

        #save p-values to file
        with open(os.path.join(feature_output_dir, f"p_values.txt"), 'w') as f:
            f.write('\n'.join(f"{x}" for x in p_values))

        return True

    def prepare_lbs_dicts(self):
        # get ligand binding sites values
        dataset = getStructuresFromDirectory(self.lbs_dir)
        for structure in dataset:
            pdb_id = structure[0]
            chain_id = structure[1]
            file = os.path.join(self.lbs_dir, f"{pdb_id}{chain_id}.txt")
            lbs = np.genfromtxt(file, delimiter=' ', dtype=None)
            lbs_dict = dict(lbs)
            self.lbs_dicts[f"{pdb_id}{chain_id}"] = lbs_dict

    def compute_pairs(self, structure):
        pdb_id = structure[0]
        chain_id = structure[1]
        error = False
        pairs = []

        try:
            lbs_dict = self.lbs_dicts[f"{pdb_id}{chain_id}"]
            # get feature values
            file = os.path.join(self.feature_dir, f"{pdb_id}{chain_id}.txt")
            feature = np.genfromtxt(file, delimiter=' ', dtype=None, encoding=None)
            feature_vals = list(feature)

            for val in feature_vals:
                res_num = val[0]
                feature_val = val[1]
                if res_num in lbs_dict:
                    lbs_val = lbs_dict[res_num]
                else:
                    continue
                pairs.append((lbs_val, feature_val))

            return pairs
        except (KeyboardInterrupt, SystemExit):
            raise
        except Exception as ex:
            error = True
            logger.debug(f"{self.feature}: Error while processing {pdb_id} {chain_id}: {ex}")
        finally:
            if (error):
                logger.error(f"{self.feature}: {pdb_id} {chain_id} NOT PROCESSED ! See log for more details.")

    def draw_plots(self, feature_type, data_binding, data_nonbinding, feature_output_dir, feature_name):
        if feature_type == "categorical":
            Plots.plot_binding_nonbinding_ratios(data_binding, data_nonbinding,
                                                 os.path.join(feature_output_dir, f"{feature_name}_ratios"), 1)
            Plots.plot_frequencies(data_binding, data_nonbinding, os.path.join(feature_output_dir, f"{feature_name}_frequencies"))
            Plots.plot_counts(data_binding, data_nonbinding,
                                   os.path.join(feature_output_dir, f"{feature_name}_counts"))
        elif feature_type == "continuous":
            Plots.plot_histogram(data_binding, data_nonbinding, 40,
                                 os.path.join(feature_output_dir, f"{feature_name}_hist_bins_40"))
            Plots.plot_histogram(data_binding, data_nonbinding, 60,
                                 os.path.join(feature_output_dir, f"{feature_name}_hist_bins_60"))
            Plots.plot_histogram(data_binding, data_nonbinding, 100,
                                 os.path.join(feature_output_dir, f"{feature_name}_hist_bins_100"))
            Plots.plot_binding_nonbinding_ratios_continuous(data_binding, data_nonbinding, os.path.join(feature_output_dir, f"{feature_name}_ratios"), 20) #todo smazat

        elif feature_type == "ordinal":
            Plots.plot_binding_nonbinding_ratios(data_binding, data_nonbinding,
                                                 os.path.join(feature_output_dir, f"{feature_name}_ratios"), 0)
            Plots.plot_frequencies(data_binding, data_nonbinding, os.path.join(feature_output_dir, f"{feature_name}_frequencies"))
            Plots.plot_counts(data_binding, data_nonbinding, os.path.join(feature_output_dir, f"{feature_name}_counts"))
        elif feature_type == "binary":
            Plots.plot_binding_nonbinding_ratios(data_binding, data_nonbinding,
                                                 os.path.join(feature_output_dir, f"{feature_name}_ratios"), 0)
            Plots.plot_positives_ratio(data_binding, data_nonbinding, os.path.join(feature_output_dir, f"{feature_name}_positivesRatio"))
            Plots.plot_counts(data_binding, data_nonbinding, os.path.join(feature_output_dir, f"{feature_name}_counts"))


    def write_summary(self):
        self.p_values.sort(key=lambda x: x[0])  # sort by name
        with open(os.path.join(self.output_dir, f"p_values.csv"), 'w') as f:
            for line in self.p_values:
                f.write(f"{str(line[0])}, ")
                f.write(', '.join('{0:.4g}'.format(x) for x in line[1]))
                f.write('\n')

        self.p_values_means.sort(key=lambda x: x[0])  # sort by name
        with open(os.path.join(self.output_dir, f"p_values_means.csv"), 'w') as f:
            f.write(f"\"feature\",\"mean\",\"std_dev\"\n") # header
            for line in self.p_values_means:
                f.write(f"{str(line[0])}, ")
                f.write(', '.join('{0:.4g}'.format(x) for x in line[1]))
                f.write('\n')

        self.cohens_d.sort(key=lambda x: x[0])  # sort by name
        with open(os.path.join(self.output_dir, f"cohens_d.csv"), 'w') as f:
            for line in self.cohens_d:
                f.write(f"{str(line[0])}, ")
                f.write(', '.join('{0:.4g}'.format(x) for x in line[1]))
                f.write('\n')

        self.cohens_d_means.sort(key=lambda x: x[0])  # sort by name
        with open(os.path.join(self.output_dir, f"cohens_d_means.csv"), 'w') as f:
            f.write(f"\"feature\",\"mean\",\"std_dev\"\n")  # header
            for line in self.cohens_d_means:
                f.write(f"{str(line[0])}, ")
                f.write(', '.join('{0:.4g}'.format(x) for x in line[1]))
                f.write('\n')

        self.cohens_w.sort(key=lambda x: x[0])  # sort by name
        with open(os.path.join(self.output_dir, f"cohens_w.csv"), 'w') as f:
            for line in self.cohens_w:
                f.write(f"{str(line[0])}, ")
                f.write(', '.join('{0:.4g}'.format(x) for x in line[1]))
                f.write('\n')

        self.cohens_w_means.sort(key=lambda x: x[0])  # sort by name
        with open(os.path.join(self.output_dir, f"cohens_w_means.csv"), 'w') as f:
            f.write(f"\"feature\",\"mean\",\"std_dev\"\n")  # header
            for line in self.cohens_w_means:
                f.write(f"{str(line[0])}, ")
                f.write(', '.join('{0:.4g}'.format(x) for x in line[1]))
                f.write('\n')

        '''self.cohens_d_vals.sort(key=lambda x: x[1], reverse=True)  # sort by value
        with open(os.path.join(self.output_dir, f"cohens_d.csv"), 'w') as f:
            for line in self.cohens_d_vals:
                f.write(f"{str(line[0])}, ")
                f.write('{0:.4g}'.format(line[1]))
                f.write('\n')

        self.cohens_w_vals.sort(key=lambda x: x[1], reverse=True)  # sort by value
        with open(os.path.join(self.output_dir, f"cohens_w.csv"), 'w') as f:
            for line in self.cohens_w_vals:
                f.write(f"{str(line[0])}, ")
                f.write('{0:.4g}'.format(line[1]))
                f.write('\n')'''

        self.b_ratios.sort(key=lambda x: x[0])  # sort by name
        with open(os.path.join(self.output_dir, f"binding_ratios.csv"), 'w') as f:
            #f.write(f"Feature, # Binding, # Nonbinding, Binding/Total\n") # header
            for line in self.b_ratios:
                f.write(f"{str(line[0])}, ")
                f.write('{0:.4g}'.format(line[1]))
                f.write('\n')

        self.means.sort(key=lambda x: x[1][0])  # sort by means difference
        with open(os.path.join(self.output_dir, f"means_difference.csv"), 'w') as f:
            f.write(f"Feature, Means difference, Mean binding, Mean nonbinding, Variance binding, Variance nonbinding\n") # header
            for line in self.means:
                f.write(f"{str(line[0])}, ")
                f.write(', '.join('{0:.4f}'.format(x) for x in line[1]))
                f.write('\n')

        self.p_val_perc.sort(key=lambda x: x[1], reverse = True)  # sort by percentage
        with open(os.path.join(self.output_dir, f"p_vals_perc.csv"), 'w') as f:
            for line in self.p_val_perc:
                f.write(f"{str(line[0])}, {str(line[1])}%")
                f.write('\n')
        with open(os.path.join(self.output_dir, f"errors.txt"), 'w') as f:
            for line in self.errors:
                f.write(f"{line}\n")

    def welchs_t_test(self, sample_binding, sample_nonbinding, f):
        out = stats.ttest_ind(sample_binding, sample_nonbinding, axis=0, equal_var=False, nan_policy='raise') #equal_var is False to run Welch's test instead of T-test
        t_statistic = out[0]
        p_value = float(out[1])

        f.write(f"T-statistic: {t_statistic}\n")
        f.write(f"P-value: {p_value}\n")
        f.write(f"Sample binding -> mean: {np.mean(sample_binding)}; variance: {np.var(sample_binding)}\n")
        f.write(f"Sample non-binding -> mean: {np.mean(sample_nonbinding)}; variance: {np.var(sample_nonbinding)}\n")
        f.write(f"Difference of means: {np.mean(sample_binding) - np.mean(sample_nonbinding)}\n")

        return (t_statistic, p_value)

    def chi_squared_test(self, sample_binding, sample_nonbinding, f, feature):
        sample_binding = [str(i) for i in sample_binding] #fix for Counter
        sample_nonbinding = [str(i) for i in sample_nonbinding]
        binding_counts = Counter(sample_binding)
        nonbinding_counts = Counter(sample_nonbinding)

        keys1 = binding_counts.keys()
        keys2 = nonbinding_counts.keys()
        categories = list(set().union(keys1, keys2))
        categories = sorted(categories)

        binding = []
        nonbinding = []
        binding_counts_sorted = []
        nonbinding_counts_sorted = []
        for cat in categories:
            binding.append(binding_counts[cat])
            nonbinding.append(nonbinding_counts[cat])
            binding_counts_sorted.append((cat, binding_counts[cat]))
            nonbinding_counts_sorted.append((cat, nonbinding_counts[cat]))

        obs = np.array([np.array(binding), np.array(nonbinding)])
        chi2, p_value, dof, expected = stats.chi2_contingency(obs)

        f.write(f"Binding counts:\n")
        f.write(f"{binding_counts_sorted}\n")
        f.write(f"Non-binding counts:\n")
        f.write(f"{nonbinding_counts_sorted}\n")
        f.write(f"Expected values: \n")
        f.write(f"{expected}\n")

        error = False
        invalid_count = 0
        #check assumptions
        for cell in np.nditer(expected):
            if (cell <= 1): #no cell can have expected value less than 1
                error = True
            if (cell <=5):
                invalid_count = invalid_count + 1
        if (invalid_count >= 0.2 * len(categories) * 2):
            error = True
        if error:
            logger.error(f"Feature {feature}: Chi-squared test assumptions not met. More than 20% of cells have expected value of 5 or less, and/or at least one cell has expected value of 1 or less.")
            return None

        p_value = float(p_value)

        f.write(f"Chi squared: {chi2}\n")
        f.write(f"p-value: {p_value}\n")

        return (chi2, p_value)

    def get_cohens_d(self, sample1, sample2):
        n1, n2 = len(sample1), len(sample2)
        s1, s2 = np.var(sample1, ddof=1), np.var(sample2, ddof=1)
        s = math.sqrt(((n1 - 1) * s1 + (n2 - 1) * s2) / (n1 + n2 - 2))
        u1, u2 = np.mean(sample1), np.mean(sample2)
        return abs((u1 - u2) / s)

    def get_cohens_w(self, chi2, N):
        '''sample_binding = [str(i) for i in sample_binding]  # fix for Counter
        sample_nonbinding = [str(i) for i in sample_nonbinding]
        binding_counts = Counter(sample_binding)
        nonbinding_counts = Counter(sample_nonbinding)
        keys1 = binding_counts.keys()
        keys2 = nonbinding_counts.keys()
        categories = list(set().union(keys1, keys2))
        categories = sorted(categories)
        binding = []
        nonbinding = []
        binding_counts_sorted = []
        nonbinding_counts_sorted = []
        for cat in categories:
            binding.append(binding_counts[cat])
            nonbinding.append(nonbinding_counts[cat])
            binding_counts_sorted.append((cat, binding_counts[cat]))
            nonbinding_counts_sorted.append((cat, nonbinding_counts[cat]))
        obs = np.array([np.array(binding), np.array(nonbinding)])
        chi2, p_value, dof, expected = stats.chi2_contingency(obs)
        N = len(sample_binding) + len(sample_nonbinding)'''  #todo
        return math.sqrt(chi2/N)


    '''def fischers_exact_test(self, data, results_file, feature):
            with open(results_file, 'w') as f:
                counts = Counter(data)

                # todo opravit..kdyz neni ani jedna 1, pada
                binding_positive = counts[(1, 1)]
                binding_negative = counts[(1, 0)]
                non_binding_positive = counts[(0, 1)]
                non_binding_negative = counts[(0, 0)]

                odds_ratio, p_value = stats.fisher_exact([[binding_positive, binding_negative],
                                                          [non_binding_positive, non_binding_negative]])
                p_value = float(p_value)

                f.write(f"Binding positive count: {binding_positive}\n")
                f.write(f"Binding negative count: {binding_negative}\n")
                f.write(f"Nonbinding positive count: {non_binding_positive}\n")
                f.write(f"Nonbinding negative count: {non_binding_negative}\n")

                f.write(f"odds ratio: {odds_ratio}\n")
                f.write(f"p-value: {p_value}\n")
                f.write(f"Binding sites -> positives ratio: {binding_positive / (binding_positive + binding_negative)}\n")
                f.write(
                    f"Non-binding sites -> positives ratio: {non_binding_positive / (non_binding_positive + non_binding_negative)}\n")

            return (feature, p_value, binding_positive + binding_negative, non_binding_positive + non_binding_negative, "Fisher's exact")'''
