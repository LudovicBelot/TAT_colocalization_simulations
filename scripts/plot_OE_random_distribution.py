#!/usr/bin/env python

import pandas as pd
import numpy as np
import re
from matplotlib import pyplot as plt
import seaborn as sns

#cl: 
#python script/plot_OE_random_distribution.py test/input/for_bargraph/bargraph_random.tsv test/0-39_genomes_Ph25/TAT_defense_results/simulation/OE_plot_2/1.Random_spot_distribution_all_elements




def main(tsv_input, outfile):
    """
    """
    
    df_input_files = pd.read_csv(tsv_input, sep = "\t")
    list_cats, list_OE, list_colors, list_pvalues = [], [], [], [] #to produce a single bargraph with all results 
    
    for _, row in df_input_files.iterrows():
            
        #Loading data
        df_real = pd.read_csv(row["real_file"], sep = "\t")
        df_sim = pd.read_csv(row["sim_file"], sep = "\t", index_col = 0)
        df_real = df_real[df_real["Element_type"] == row["elem_type"]]
        real_data = df_real[f"N_NR_genes_in_spot"]
        
        # 1. Moyennes des simulations
        simulation_means = df_sim.mean(axis=0)
        simulation_global_mean = simulation_means.mean()
        
        # 2. Moyenne des données réelles
        real_data = real_data.squeeze()
        real_mean = real_data.mean()
        
        p_value_lower = (simulation_means <= real_mean).sum() / len(simulation_means)
        p_value_higher = (simulation_means >= real_mean).sum() / len(simulation_means)
        OE_value = (real_mean - simulation_global_mean) / (real_mean + simulation_global_mean)       
        
        #Storing data for the complete
        if row["elem_type"] == "IS":
            list_cats.append("transposase")
        elif row["elem_type"] == "defense":
            list_cats.append("antiphage")
        else :     
            list_cats.append(row["elem_type"])
        list_OE.append(OE_value)
        list_colors.append(row["color"])
        list_pvalues.append(min(p_value_lower, p_value_higher))
    
    plot_full_OE(list_cats, list_OE, list_pvalues, list_colors, outfile)




def plot_full_OE(list_cats, list_OE_vals, list_pvalues, list_colors, outfile):
    """
    """
    
    #Pvalues_stars 
    list_stars = []
    for p in list_pvalues:
        list_stars.append('***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns')

    if len(list_cats) <= 4:
        fig, ax = plt.subplots(figsize=(3.5, 3.5), dpi=150)
    else :
        fig, ax = plt.subplots(figsize=(6, 3.5), dpi=150)
    
    x = np.arange(len(list_OE_vals)) * 0.5
    
    bars = ax.bar(
        x, list_OE_vals, width=0.25,
        color=list_colors, edgecolor="black", linewidth=0.5
    )
    
    list_labels = [re.sub(r'\s*&\s*', '\n', c) for c in list_cats]
    ax.set_xticks(x)
    ax.tick_params(axis='x', labeltop=True, labelbottom=False)
    ax.tick_params(axis='x', which='both', length=0)
    ax.set_xticklabels(list_labels, fontname = "Arial", fontsize = 8.5)

    #plot stars:
    for i, (bar, star) in enumerate(zip(bars, list_stars)):
        x_i = bar.get_x() + bar.get_width() / 2 
        y_i = list_OE_vals[i] + 0.0075 if list_OE_vals[i] >= 0 else list_OE_vals[i] - 0.05
        ax.text(x_i, y_i, star, ha='center', va='bottom')
    
    # Draw bracket and label
    ax.set_ylim(-0.1, 0.2)
    ax.set_yticks(np.arange(-0.1, 0.21, 0.1))
    #ax.set_ylabel("Average number of a specific element within the spot given another type\n(O-E)/(O+E)")
    ax.set_ylabel("(O-E)/(O+E)\nAverage spot size associated to an element", fontname = "Arial", fontsize = 12)
    ax.axhline(0, color='black', linewidth=0.8)
    sns.despine(ax=ax, top=True, right=True, bottom = True)
    plt.tight_layout()
    plt.savefig(outfile+"_All_distribution_results_from_"+ ("-").join([x for x in set(list_cats)]) + "_perspectives_OE_graph.svg", format = "svg")
    plt.close(fig)














if __name__ == "__main__":
    import sys
    main(sys.argv[1], sys.argv[2])