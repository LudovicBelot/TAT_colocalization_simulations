#!/usr/bin/env python

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import re

#Cl:
#python script/plot_violin_coloc_distribution.py test/input/for_bargraph/bargraph_coloc_TA_perspectives.tsv test/0-39_genomes_Ph25/TAT_defense_results/simulation/violin_plot/1.Coloc_simulations

#cl: python script/plot_violin_coloc_distribution.py test/input/for_bargraph/bargraph_coloc_defense_perspectives.tsv test/0-39_genomes_Ph25/TAT_defense_results/simulation/violin_plot/2.Coloc_simulations
#cl: python script/plot_violin_coloc_distribution.py test/input/for_bargraph/bargraph_coloc_IS_perspectives.tsv test/0-39_genomes_Ph25/TAT_defense_results/simulation/violin_plot/3.Coloc_simulations
#cl: python script/plot_violin_coloc_distribution.py test/input/for_bargraph/bargraph_coloc_integrase_perspectives.tsv test/0-39_genomes_Ph25/TAT_defense_results/simulation/violin_plot/4.Coloc_simulations

#cl: python script/plot_violin_coloc_distribution.py test/input/for_bargraph/bargraph_coloc_TA_perspectives_IS_categories.tsv test/0-39_genomes_Ph25/TAT_defense_results/simulation/violin_plot/5.Coloc_TA_with_simulations_IS_categories
#cl: python script/plot_violin_coloc_distribution.py test/input/for_bargraph/bargraph_coloc_defense_perspectives_IS_categories.tsv test/0-39_genomes_Ph25/TAT_defense_results/simulation/violin_plot/6.Coloc_defense_with_simulations_IS_categories

#cl: python script/plot_violin_coloc_distribution.py test/input/for_bargraph/bargraph_coloc_TA_perspectives_integrase_categories.tsv test/0-39_genomes_Ph25/TAT_defense_results/simulation/violin_plot/7.Coloc_TA_with_simulations_integrase_categories
#cl: python script/plot_violin_coloc_distribution.py test/input/for_bargraph/bargraph_coloc_defense_perspectives_integrase_categories.tsv test/0-39_genomes_Ph25/TAT_defense_results/simulation/violin_plot/8.Coloc_defense_with_simulations_integrase_categories

#comparison between TA & defense with MGEs 
#cl: python script/plot_violin_coloc_distribution.py test/input/for_bargraph/bargraph_coloc_comparison_TA_defense_with_MGE_categories.tsv test/0-39_genomes_Ph25/TAT_defense_results/simulation/violin_plot/9.Coloc_with_specific_MGEs_comparison_TA_defense



def main(tsv_input, outfile):
    """
    """
    
    df_input_files = pd.read_csv(tsv_input, sep = "\t")
    list_cats, list_elem1, list_sim_values, list_normalized_sim_values, list_real_values, list_normalized_real_values, list_colors, list_pvalues = [], [], [], [], [], [], [], [] #to produce a single bargraph with all results 
    
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

        # 1. Moyennes des simulations
        simulation_means = df_sim.mean(axis=0)

        # 2. Moyenne des données réelles
        real_data = real_data.squeeze()
        real_mean = real_data.mean()
        p_value_lower = (simulation_means <= real_mean).sum() / len(simulation_means)
        p_value_higher = (simulation_means >= real_mean).sum() / len(simulation_means)
        
        #Storing data for the complete
        list_cats.append(row["elem2_type"])
        list_sim_values.append(simulation_means)
        list_normalized_sim_values.append([x/max(simulation_means) for x in simulation_means])
        list_real_values.append(real_mean)
        list_normalized_real_values.append(real_mean/max(simulation_means))
        list_colors.append(row["color"])
        list_pvalues.append(min(p_value_lower, p_value_higher))
        list_elem1.append(row["elem1_type"])
    
    
    plot_full_violin(list_cats, list_sim_values, list_real_values, list_pvalues, list_colors, list_elem1, outfile)
    plot_normalized_full_violin(list_cats, list_normalized_sim_values, list_normalized_real_values, list_pvalues, list_colors, list_elem1, outfile)


def plot_full_violin(list_cats, list_sim_values, list_real_values, list_pvalues, list_colors, list_elem_perspective_type, outfile):
    """
    """
    
    #Pvalues_stars 
    list_stars = []
    for p in list_pvalues:
        list_stars.append('***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns')

    if len(list_cats) <= 4:
        fig, ax = plt.subplots(figsize=(6.25, 3.5), dpi=150)
    else :
        fig, ax = plt.subplots(figsize=(8, 3.5), dpi=150)
    
    positions = np.arange(1, len(list_sim_values) + 1)
    parts = ax.violinplot(list_sim_values, positions=positions, widths=0.8, showmeans=False, showextrema=False)
    ymax = 0
    for i in list_sim_values:
        if max(i) > ymax:
            ymax = max(i)
    ymax = max(ymax, max(list_real_values)) + 0.05
    ax.set_ylim(bottom=0, top = ymax)
    
    list_labels = [re.sub(r'\s*&\s*', '\n', c) for c in list_cats]
    set_axis_style(ax, list_labels)
    
    #Changing the colors according to the element perspectives:
    for i, pc in enumerate(parts['bodies']):
        pc.set_facecolor(list_colors[i])
        pc.set_edgecolor('black')
        pc.set_alpha(0.65)
    
    #plot stars:
    for i, star in enumerate(list_stars):
        sim_max = float(np.nanmax(np.atleast_1d(list_sim_values[i])))
        real_i = float(list_real_values[i])
        y_i = max(sim_max, real_i)
        ax.text(positions[i], y_i, star, ha='center', va='bottom', clip_on=False, zorder=5)
    
    
    #Adding whiskers representing the quartiles
    quartile1, quartile3 = np.percentile(list_sim_values, [25, 75], axis=1)
    whiskers_min = [min(simulation) for simulation in list_sim_values]
    whiskers_max = [max(simulation) for simulation in list_sim_values]
    
    ax.vlines(positions, quartile1, quartile3, color='k', linestyle='-', lw=3.5)
    ax.vlines(positions, whiskers_min, whiskers_max, color='k', linestyle='-', lw=0.5)
    
    #Adding a red sign representing the real value
    real = np.asarray(list_real_values, dtype=float)
    mask = np.isfinite(real)
    ax.scatter(positions[mask], real[mask], marker='x', color='red', s=40, linewidths=1, zorder=6)
    
    #ax.set_ylabel("Average number of a specific element within the spot given another type\n(O-E)/(O+E)")
    ax.set_ylabel(f"Average number of an element\nbeing retrieved in the same spot")
    sns.despine(ax=ax, top=True, right=True)
    plt.tight_layout()
    plt.savefig(outfile+"_All_coloc_results_from_"+ ("-").join([x for x in set(list_elem_perspective_type)]) + "_perspectives_violin_graph.png")
    plt.close(fig)


def plot_normalized_full_violin(list_cats, list_sim_values, list_real_values, list_pvalues, list_colors, list_elem_perspective_type, outfile):
    """
    """
    
    #Pvalues_stars 
    list_stars = []
    for p in list_pvalues:
        list_stars.append('***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns')

    if len(list_cats) <= 4:
        fig, ax = plt.subplots(figsize=(6.25, 3.5), dpi=150)
    else :
        fig, ax = plt.subplots(figsize=(8, 3.5), dpi=150)
    
    positions = np.arange(1, len(list_sim_values) + 1)
    parts = ax.violinplot(list_sim_values, positions=positions, widths=0.8, showmeans=False, showextrema=False)
    list_labels = [re.sub(r'\s*&\s*', '\n', c) for c in list_cats]
    
    ymax = 0
    for i in list_sim_values:
        if max(i) > ymax:
            ymax = max(i)
    ymax = max(ymax, max(list_real_values))
    ymax = ymax + 0.05
    ax.set_ylim(bottom=0, top = ymax)
    set_axis_style(ax, list_labels)

    #Changing the colors according to the element perspectives:
    for i, pc in enumerate(parts['bodies']):
        pc.set_facecolor(list_colors[i])
        pc.set_edgecolor('black')
        pc.set_alpha(0.65)
    
    #plot stars:
    for i, star in enumerate(list_stars):
        sim_max = float(np.nanmax(np.atleast_1d(list_sim_values[i])))
        real_i = float(list_real_values[i])
        y_i = max(sim_max, real_i)
        ax.text(positions[i], y_i, star, ha='center', va='bottom', clip_on=False, zorder=5)
    
    
    #Adding whiskers representing the quartiles
    quartile1, quartile3 = np.percentile(list_sim_values, [25, 75], axis=1)
    whiskers_min = [min(simulation) for simulation in list_sim_values]
    whiskers_max = [max(simulation) for simulation in list_sim_values]
    
    ax.vlines(positions, quartile1, quartile3, color='k', linestyle='-', lw=3.5)
    ax.vlines(positions, whiskers_min, whiskers_max, color='k', linestyle='-', lw=0.5)
    
    #Adding a red sign representing the real value
    real = np.asarray(list_real_values, dtype=float)
    mask = np.isfinite(real)
    ax.scatter(positions[mask], real[mask], marker='x', color='red', s=40, linewidths=1, zorder=6)
    
    #y axis title
    ax.set_ylabel(f"Average number of an element\nnormalized with the highest value\nobtained through random distribution")
    sns.despine(ax=ax, top=True, right=True)
    plt.tight_layout()
    plt.savefig(outfile+"_All_coloc_results_from_"+ ("-").join([x for x in set(list_elem_perspective_type)]) + "_perspectives_normalized_violin_graph.png")
    plt.close(fig)







def adjacent_values(vals, q1, q3):
    upper_adjacent_value = q3 + (q3 - q1) * 1.5
    upper_adjacent_value = np.clip(upper_adjacent_value, q3, vals.iloc[-1])
    lower_adjacent_value = q1 - (q3 - q1) * 1.5
    lower_adjacent_value = np.clip(lower_adjacent_value, vals.iloc[0], q1)
    return lower_adjacent_value, upper_adjacent_value

def set_axis_style(ax, labels):
    ax.set_xticks(np.arange(1, len(labels) + 1), labels=labels)
    ax.set_xlim(0.25, len(labels) + 0.75)
    ax.set_xlabel('Element association')


if __name__ == "__main__":
    import sys
    main(sys.argv[1], sys.argv[2])