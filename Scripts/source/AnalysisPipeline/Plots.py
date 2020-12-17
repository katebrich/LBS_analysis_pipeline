import math
import operator
import numpy as np
from collections import Counter
from matplotlib import pyplot, transforms
from matplotlib.ticker import FormatStrFormatter

def plot_pvalues_scatter(p_values, alpha, output_file):
    pyplot.clf()
    fig, ax = pyplot.subplots()
    pyplot.scatter(range(1, len(p_values) + 1), p_values)
    axes = pyplot.gca()
    axes.set_ylim([0, 1]) #set range of y values
    pyplot.xlabel('Iteration', fontsize=13)
    pyplot.ylabel('P-value', fontsize=13)
    ax.axhline(y=alpha, color="red", linestyle='-') #draw horizontal line for alpha value
    trans = transforms.blended_transform_factory(
        ax.get_yticklabels()[0].get_transform(), ax.transData)
    ax.text(0, alpha, "{:.2f}".format(alpha), color="red", transform=trans,
            ha="right", va="center")
    pyplot.savefig(output_file)
    pyplot.clf()
    pyplot.close('all')

def plot_pvalues_histogram(p_values, iterations, output_file):
    pyplot.clf()
    bins = math.ceil(math.sqrt(iterations))
    pyplot.hist(p_values, bins=bins, color=('cornflowerblue'), alpha=1)
    pyplot.xlabel('P-value', fontsize=13)
    pyplot.ylabel('Frequency', fontsize=13)
    axes = pyplot.gca()
    #axes.set_xlim([0, 1])
    pyplot.savefig(output_file)
    pyplot.clf()
    pyplot.close('all')

def plot_binding_nonbinding_ratios(binding_data, nonbinding_data, output_file, sort_item):
    pyplot.clf()
    binding_data = [str(i) for i in binding_data]
    nonbinding_data = [str(i) for i in nonbinding_data]
    binding_counts = Counter(binding_data)
    nonbinding_counts = Counter(nonbinding_data)
    if (len(binding_counts) < 2 or len(nonbinding_counts) < 2):
        return

    ratios={}
    for k in nonbinding_counts.keys():
        ratio = binding_counts[k] / nonbinding_counts[k]
        ratios[k] = ratio

    if (sort_item == 0): #sort by name
        sorted_ratios = sorted(ratios.items(), key=operator.itemgetter(0))
    elif (sort_item == 1): #sort by ratios values
        sorted_ratios = sorted(ratios.items(), key=operator.itemgetter(1), reverse=True)
    else: #no sorting
        sorted_ratios = ratios.items()

    barWidth = 0.45

    # set height of bar
    bars = [x[1] for x in sorted_ratios]
    # Set position of bar on X axis
    r1 = np.arange(len(bars))
    r2 = [x + barWidth/2 for x in r1]

    fig, ax = pyplot.subplots()
    pyplot.bar(r2, bars, color='cornflowerblue', width=barWidth)

    pyplot.xlabel('Value', fontsize=13)
    pyplot.ylabel('Binding/nonbinding ratio', fontsize=13)
    # Add xticks in the middle of bars
    pyplot.xticks([r + barWidth/2 for r in range(len(bars))], [x[0] for x in sorted_ratios])

    #Draw total binding/nonbinding ratio for comparison
    ratio = len(binding_data) / len(nonbinding_data)
    ax.axhline(y=ratio, color="red", linestyle='--' )

    trans = transforms.blended_transform_factory(
        ax.get_yticklabels()[0].get_transform(), ax.transData)
    ax.text(0, ratio, "{:.3f}".format(ratio), color="red", transform=trans,
            ha="right", va="center")

    #round y axis labels to 3 places
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))

    #save figure and clear
    pyplot.savefig(output_file)
    pyplot.clf()
    pyplot.close('all')


def plot_histogram(binding_data, nonbinding_data, bins, output_file):
    pyplot.clf()
    weights_binding = np.ones_like(binding_data)/len(binding_data)
    weights_nonBinding = np.ones_like(nonbinding_data) / len(nonbinding_data)
    pyplot.hist((binding_data, nonbinding_data), bins=bins, weights=(weights_binding, weights_nonBinding), color=('green', 'red'), alpha=0.5,
                label=('Binding', 'Non-binding'))

    pyplot.legend(loc='upper right')
    pyplot.xlabel('Value', fontsize=13)
    pyplot.ylabel('Density', fontsize=13)
    pyplot.savefig(output_file)
    pyplot.clf()
    pyplot.close('all')

def plot_frequencies(binding_data, nonbinding_data, output_file):
    pyplot.clf()
    binding_data = [str(i) for i in binding_data]
    nonbinding_data = [str(i) for i in nonbinding_data]
    binding_counts = Counter(binding_data)
    nonBinding_counts = Counter(nonbinding_data)

    binding_frequencies = {k: v / len(binding_data) for k, v in binding_counts.items()}
    nonBinding_frequencies = {k: v / len(nonbinding_data) for k, v in nonBinding_counts.items()}

    binding_frequencies_sorted = {}
    nonBinding_frequencies_sorted = {}
    keys = sorted(binding_frequencies.keys()) # categories sorted by name
    for key in keys:
        freq = binding_frequencies.get(key, 0) #default is 0
        binding_frequencies_sorted[key] = freq
        freq = nonBinding_frequencies.get(key, 0)  # default is 0
        nonBinding_frequencies_sorted[key] = freq

    #plot the results
    barWidth = 0.33

    # set height of bar
    bars1 = binding_frequencies_sorted.values()
    bars2 = nonBinding_frequencies_sorted.values()

    # Set position of bar on X axis
    r1 = np.arange(len(bars1))
    r2 = [x + barWidth for x in r1]

    pyplot.bar(r1, bars1, color='g', width=barWidth, edgecolor='white', label='Binding', alpha=0.6)
    pyplot.bar(r2, bars2, color='r', width=barWidth, edgecolor='white', label='Non-binding', alpha=0.6)

    # Add xticks in the middle of bars
    pyplot.xlabel('Value', fontsize=13)
    pyplot.ylabel('Frequency', fontsize=13)
    pyplot.xticks([r + barWidth/2 for r in range(len(bars1))], keys)
    pyplot.legend(fontsize=10, loc='upper right')

    pyplot.savefig(output_file)
    pyplot.clf()
    pyplot.close('all')


def plot_counts(binding_data, nonbinding_data, output_file):
    pyplot.clf()
    binding_data = [str(i) for i in binding_data]
    nonbinding_data = [str(i) for i in nonbinding_data]
    binding_counts = Counter(binding_data)
    nonBinding_counts = Counter(nonbinding_data)

    binding_counts_sorted = {}
    nonBinding_counts_sorted = {}
    keys = sorted(binding_counts.keys()) # categories sorted by name
    for key in keys:
        val = binding_counts.get(key, 0) #default is 0
        binding_counts_sorted[key] = val
        val = nonBinding_counts.get(key, 0)  # default is 0
        nonBinding_counts_sorted[key] = val

    #plot the results
    barWidth = 0.33

    # set height of bar
    bars1 = binding_counts_sorted.values()
    bars2 = nonBinding_counts_sorted.values()

    # Set position of bar on X axis
    r1 = np.arange(len(bars1))
    r2 = [x + barWidth for x in r1]

    pyplot.bar(r1, bars1, color='g', width=barWidth, edgecolor='white', label='Binding sites', alpha=0.6)
    pyplot.bar(r2, bars2, color='r', width=barWidth, edgecolor='white', label='Non-binding sites', alpha=0.6)

    # Add xticks in the middle of bars
    pyplot.xlabel('Value', fontsize=13)
    pyplot.ylabel('Count', fontsize=13)
    pyplot.xticks([r + barWidth/2 for r in range(len(bars1))], keys)
    pyplot.legend()

    pyplot.savefig(output_file)
    pyplot.clf()
    pyplot.close('all')

def plot_positives_ratio(binding_data, nonbinding_data, output_file):
    pyplot.clf()
    binding_data = [str(i) for i in binding_data]
    nonbinding_data = [str(i) for i in nonbinding_data]
    binding_counts = Counter(binding_data)
    nonbinding_counts = Counter(nonbinding_data)

    if (nonbinding_counts['0'] == 0 or nonbinding_counts['1'] == 0):
        return

    #plot the results
    barWidth = 0.33

    # set height of bar as positives ratio
    bars = (binding_counts['1']/(binding_counts['0']+binding_counts['1']),
             nonbinding_counts['1']/(nonbinding_counts['0']+nonbinding_counts['1']))

    # Set position of bar on X axis
    r1 = np.arange(len(bars))
    r2 = [x + barWidth / 2 for x in r1]

    pyplot.bar(r2, bars, width=barWidth, edgecolor='white', color='cornflowerblue')

    # Add xticks in the middle of bars
    pyplot.ylabel('Positives ratio', fontsize=13)
    pyplot.xticks([r + barWidth/2 for r in range(len(bars))], ("Binding","Nonbinding"))

    pyplot.savefig(output_file)
    pyplot.clf()
    pyplot.close('all')