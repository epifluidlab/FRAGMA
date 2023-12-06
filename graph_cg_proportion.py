import ast
import matplotlib.pyplot as plt
import numpy as np

def add_lists(list1, list2):
    if len(list1) == len(list2):
        for i in range(len(list1)):
            list1[i] += list2[i]
    return list1

def graph_cg_proportion(cg_file):
    combined_array_methylated=[0,0,0,0,0,0,0,0,0,0,0]
    total_ads_methylated = 0
    combined_array_unmethylated=[0,0,0,0,0,0,0,0,0,0,0]
    total_ads_unmethylated = 0
    with open(cg_file, "r") as cg_sites:
            for line in cg_sites:
                fields = line.strip().split("\t")
                methylation = fields[2]
                positive_array = fields[6]
                positive_array_frag = int(fields[7])
                negative_array = fields[8]
                negative_array_frag = int(fields[9])
                positive_array = ast.literal_eval(positive_array)
                negative_array = ast.literal_eval(negative_array)
                positive_array = [x / positive_array_frag for x in positive_array]
                negative_array = [x / negative_array_frag for x in negative_array]
                if float(methylation)>=70:
                    combined_array_methylated = add_lists(combined_array_methylated, positive_array)
                    combined_array_methylated = add_lists(combined_array_methylated, negative_array)
                    total_ads_methylated += 2
                else:
                    combined_array_unmethylated = add_lists(combined_array_unmethylated, positive_array)
                    combined_array_unmethylated = add_lists(combined_array_unmethylated, negative_array)
                    total_ads_unmethylated += 2
    combined_array_methylated = [x / total_ads_methylated for x in combined_array_methylated]
    combined_array_unmethylated = [x / total_ads_unmethylated for x in combined_array_unmethylated]
    
    x_positions = np.arange(-5, 6)
    plt.plot(x_positions, combined_array_methylated, marker='o', linestyle='-', color='red', label='Methylated')
    plt.plot(x_positions, combined_array_unmethylated, marker='o', linestyle='-', color='blue', label='Unmethylated')
    plt.xlabel('11-nt window')
    plt.ylabel('Cleavage Proportion')
    plt.legend()
    plt.savefig('Cleavage_Proportion_11_UNMASKED_FILT.png')



graph_cg_proportion("/jet/home/rbandaru/ocean_rbandaru/FRAGMA/UNMASKED.BED")