#!/usr/bin/env python

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import re
from collections import OrderedDict

#Cl:
#python script/plot_OE_coloc_distribution.py test/input/for_bargraph/bargraph_coloc_TA_perspectives.tsv test/0-39_genomes_Ph25/TAT_defense_results/simulation/OE_plot/1.Coloc_simulations

#cl: python script/plot_OE_coloc_distribution.py test/input/for_bargraph/bargraph_coloc_defense_perspectives.tsv test/0-39_genomes_Ph25/TAT_defense_results/simulation/OE_plot/2.Coloc_simulations
#cl: python script/plot_OE_coloc_distribution.py test/input/for_bargraph/bargraph_coloc_IS_perspectives.tsv test/0-39_genomes_Ph25/TAT_defense_results/simulation/OE_plot/3.Coloc_simulations
#cl: python script/plot_OE_coloc_distribution.py test/input/for_bargraph/bargraph_coloc_integrase_perspectives.tsv test/0-39_genomes_Ph25/TAT_defense_results/simulation/OE_plot/4.Coloc_simulations

#cl: python script/plot_OE_coloc_distribution.py test/input/for_bargraph/bargraph_coloc_TA_perspectives_IS_categories.tsv test/0-39_genomes_Ph25/TAT_defense_results/simulation/OE_plot/5.Coloc_TA_with_simulations_IS_categories
#cl: python script/plot_OE_coloc_distribution.py test/input/for_bargraph/bargraph_coloc_defense_perspectives_IS_categories.tsv test/0-39_genomes_Ph25/TAT_defense_results/simulation/OE_plot/6.Coloc_defense_with_simulations_IS_categories

#cl: python script/plot_OE_coloc_distribution.py test/input/for_bargraph/bargraph_coloc_TA_perspectives_integrase_categories.tsv test/0-39_genomes_Ph25/TAT_defense_results/simulation/OE_plot/7.Coloc_TA_with_simulations_integrase_categories
#cl: python script/plot_OE_coloc_distribution.py test/input/for_bargraph/bargraph_coloc_defense_perspectives_integrase_categories.tsv test/0-39_genomes_Ph25/TAT_defense_results/simulation/OE_plot/8.Coloc_defense_with_simulations_integrase_categories

#comparison between TA & defense with MGEs 
#cl: python script/plot_OE_coloc_distribution.py test/input/for_bargraph/bargraph_coloc_comparison_TA_defense_with_MGE_categories.tsv test/0-39_genomes_Ph25/TAT_defense_results/simulation/OE_plot/9.Coloc_with_specific_MGEs_comparison_TA_defense

#Colocalization TA/defense with other elements same graph
#python script/plot_OE_coloc_distribution.py test/input/for_bargraph/bargraph_coloc_TA-defense_perspectives.tsv test/0-39_genomes_Ph25/TAT_defense_results/simulation/OE_plot_2/10.Coloc_TA-defense_perspectives

#IS cointegrate 
#python script/plot_OE_coloc_distribution.py test/input/for_bargraph/bargraph_coloc_TA-defense_perspectives_IS_categories.tsv test/0-39_genomes_Ph25/TAT_defense_results/simulation/OE_plot_2/11.Coloc_TA-defense_IS_all_categories



def main(tsv_input, outfile):
    """
    """
    
    df_input_files = pd.read_csv(tsv_input, sep = "\t")
    list_cats, list_elem1, list_OE_vals, list_colors, list_pvalues = [], [], [], [], [] #to produce a single bargraph with all results 
    
    for _, row in df_input_files.iterrows():
            
        #Loading data
        df_real = pd.read_csv(row["real_file"], sep = "\t")
        df_sim = pd.read_csv(row["sim_file"], sep = "\t", index_col = 0)
        df_real = df_real[df_real["Element_type"] == row["elem1_type"]]
        
        #Added a few lines that if eleme2_type not in ["TA", "defense", "IS", "integrase"] then it will consider the custom group
        if row['elem2_type'] in ["TA", "defense", "IS", "integrase", "Phage_integrase", "Resolvase"]:
            elem2_type = f"N_{row['elem2_type']}_colocalizing"
        else :
            elem2_type = "N_custom_group_colocalization"
        
        real_data = df_real[elem2_type]

        simulation_means = df_sim.mean(axis=0)
        simulation_global_mean = simulation_means.mean()

        real_data = real_data.squeeze()
        real_mean = real_data.mean()
        p_value_lower = (simulation_means <= real_mean).sum() / len(simulation_means)
        p_value_higher = (simulation_means >= real_mean).sum() / len(simulation_means)
        OE_value = (real_mean - simulation_global_mean) / (real_mean + simulation_global_mean)
        #plot_OE_graph(OE_value, p_value_lower, p_value_higher, row["elem1_type"], row["elem2_type"], row["color"], outfile)
                
        list_cats.append(row["elem2_type"])
        list_OE_vals.append(OE_value)
        list_colors.append(row["color"])
        list_pvalues.append(min(p_value_lower, p_value_higher))
        list_elem1.append(row["elem1_type"])
    
    
    plot_full_OE(list_cats, list_OE_vals, list_pvalues, list_colors, list_elem1, outfile)
    



def plot_OE_graph(OE_value, p_value_lower, p_value_higher, elem1_type, elem2_type, color, outfile):


    fig, ax = plt.subplots(figsize=(2, 4), dpi=300)
    heights = [OE_value]

    x = np.array([1])
    bars = ax.bar(
        x, heights, width=0.5,
        color=color, edgecolor="black", linewidth=0.7
    )
    
    ax.set_xticks(x)
    ax.set_xticklabels([elem2_type])

    # Significance bracket using the lower p-value
    p = min(p_value_lower, p_value_higher)
    stars = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'

    y_max = OE_value
    h = 0.05 * y_max
    y = y_max + h

    ax.text((x), y + h, stars, ha='center', va='bottom')
    ax.set_ylim(-1, 1)
    ax.set_ylabel(f"Average number of {elem2_type} within the spot given a {elem1_type}\n(O-E)/(O+E)")
    #ax.set_title(f"{elem_type} – Real vs Simulation")
    ax.axhline(0, color='black', linewidth=0.8)
    sns.despine(ax=ax, top=True, right=True)
    plt.tight_layout()
    plt.savefig(outfile + "_distribution_of_"+ elem2_type + "_given_a_" + elem1_type + "_OE_plot.png", format = "svg")
    plt.close(fig)

def plot_full_OE(list_cats, list_OE_vals, list_pvalues, list_colors, list_elem_perspective_type, outfile):
    """
    """

    #Pvalues_stars 
    list_stars = []
    for p in list_pvalues:
        list_stars.append('***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns')

    if len(list_cats) <= 2:
        fig, ax = plt.subplots(figsize=(1.75, 3.5), dpi=150)
    elif len(list_cats) <= 4:
        fig, ax = plt.subplots(figsize=(5.5, 3.5), dpi=150)
    else :
        fig, ax = plt.subplots(figsize=(7, 3.5), dpi=150)


    bar_spacing = 0.325   # spacing between bars within a same group
    group_gap = 0.6     # spacing between groups

    groups = OrderedDict()
    for i, g in enumerate(list_elem_perspective_type):
        groups.setdefault(g, []).append(i)

    x = np.zeros(len(list_OE_vals))
    current_x = 0
    group_bounds = []  # store (group_label, start_idx, end_idx)

    for group_label, indices in groups.items():
        for j, idx in enumerate(indices):
            x[idx] = current_x + j * bar_spacing

        group_bounds.append((group_label, indices[0], indices[-1]))
        current_x = x[indices[-1]] + group_gap


    separator_x = []
    for i in range(len(group_bounds) - 1):
        _, _, end_idx = group_bounds[i]
        _, start_next, _ = group_bounds[i + 1]
        separator_x.append((x[end_idx] + x[start_next]) / 2)

    bars = ax.bar(
        x, list_OE_vals, width=0.175,
        color=list_colors, edgecolor="black", linewidth=0.5
    )

    #Comment or uncomment here depending if you have subclasses to plot
    #list_cats = [re.sub(r'\s*&\s*', '\n', c) for c in list_cats]
    list_cats = [x.split("(")[1].rsplit(")")[0] for x in list_cats]

    #Added few lines to change the names of specific categories (defense -> antiphages, IS -> transposases)
    list_cats_formated = []
    for i in list_cats:
        if i == "defense":
            list_cats_formated.append("antiphage")
        elif i == "IS":
            list_cats_formated.append("transposase")
        else :
            list_cats_formated.append(i)

    list_cats = list_cats_formated

    ax.set_xticks(x)
    ax.tick_params(axis='x', labeltop=True, labelbottom=False)
    ax.tick_params(axis='x', which='both', length=0)
    ax.set_xticklabels(list_cats, fontsize = 8.5, fontname = "Arial")

    #plot stars:
    for i, (bar, star) in enumerate(zip(bars, list_stars)):
        x_i = bar.get_x() + bar.get_width() / 2 
        y_i = list_OE_vals[i] + 0.025 if list_OE_vals[i] >= 0 else list_OE_vals[i] - 0.125
        ax.text(x_i, y_i, star, ha='center', va='bottom')


    #Comment or uncomment depending on the y limits you choose below
    y_bracket = 0.9
    #y_bracket = 0.45

    for label, start, end in group_bounds:
        x_start = x[start]
        x_end = x[end]
        x_center = (x_start + x_end) / 2
        
        if label == "defense":
            label = "antiphage"


        ax.plot([x_start, x_end], [y_bracket, y_bracket],
                color='black', lw=0.8, clip_on=False)

        ax.plot([x_start, x_start], [y_bracket, y_bracket - 0.02],
                color='black', lw=0.8, clip_on=False)
        ax.plot([x_end, x_end], [y_bracket, y_bracket - 0.02],
                color='black', lw=0.8, clip_on=False)

        ax.text(x_center, y_bracket + 0.02,
                f"From {label} perspective",
                ha='center', va='bottom', fontsize=10, fontname = "Arial")

    #Comment or uncomment the choosen y limites
    ax.set_ylim(-0.35, 0.8)
    #ax.set_ylim(-0.35, 0.35)
    
    for x_sep in separator_x:
        ax.axvline(
            x=x_sep,
            color='grey',
            linestyle='--',
            linewidth=0.6,
            alpha=0.7,
            zorder=0
        )
    
    ax.set_yticks(np.arange(-0.3, 0.81, 0.15))
    #ax.set_yticks(np.arange(-0.3, 0.36, 0.15))

    #ax.set_ylabel("Average number of a specific element within the spot given another type\n(O-E)/(O+E)")
    #Comment or uncomment here depending on what Y axis title you want
    #ax.set_ylabel("(O-E)/(O+E)\nCo-occuring elements", fontname = "Arial", fontsize = 12)
    ax.set_ylabel("(O-E)/(O+E)\nCo-occuring transposase genes", fontname = "Arial", fontsize = 12)

    ax.axhline(0, color='black', linewidth=0.8)
    sns.despine(ax=ax, top=True, right=True, bottom = True)
    plt.tight_layout()

    plt.savefig(
        outfile+"_All_coloc_results_from_"+("-").join(groups.keys())+"_perspectives_OE_graph.svg",
        format="svg"
    )
    plt.close(fig)





if __name__ == "__main__":
    import sys
    main(sys.argv[1], sys.argv[2])