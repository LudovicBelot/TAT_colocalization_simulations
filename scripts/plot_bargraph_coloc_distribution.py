#!/usr/bin/env python

import pandas as pd
import numpy as np
from tqdm import tqdm
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import seaborn as sns
import re


#cl: python script/plot_bargraph_coloc_distribution.py test/input/for_bargraph/bargraph_coloc_TA_perspectives.tsv test/0-39_genomes_Ph25/TAT_defense_results/simulation/bargraph_plot/1.Coloc_simulations
#cl: python script/plot_bargraph_coloc_distribution.py test/input/for_bargraph/bargraph_coloc_defense_perspectives.tsv test/0-39_genomes_Ph25/TAT_defense_results/simulation/bargraph_plot/2.Coloc_simulations
#cl: python script/plot_bargraph_coloc_distribution.py test/input/for_bargraph/bargraph_coloc_IS_perspectives.tsv test/0-39_genomes_Ph25/TAT_defense_results/simulation/bargraph_plot/3.Coloc_simulations
#cl: python script/plot_bargraph_coloc_distribution.py test/input/for_bargraph/bargraph_coloc_integrase_perspectives.tsv test/0-39_genomes_Ph25/TAT_defense_results/simulation/bargraph_plot/4.Coloc_simulations

#cl: python script/plot_bargraph_coloc_distribution.py test/input/for_bargraph/bargraph_coloc_TA_perspectives_IS_categories.tsv test/0-39_genomes_Ph25/TAT_defense_results/simulation/bargraph_plot/5.Coloc_TA_with_simulations_IS_categories
#cl: python script/plot_bargraph_coloc_distribution.py test/input/for_bargraph/bargraph_coloc_defense_perspectives_IS_categories.tsv test/0-39_genomes_Ph25/TAT_defense_results/simulation/bargraph_plot/6.Coloc_defense_with_simulations_IS_categories

#cl: python script/plot_bargraph_coloc_distribution.py test/input/for_bargraph/bargraph_coloc_TA_perspectives_integrase_categories.tsv test/0-39_genomes_Ph25/TAT_defense_results/simulation/bargraph_plot/7.Coloc_TA_with_simulations_integrase_categories
#cl: python script/plot_bargraph_coloc_distribution.py test/input/for_bargraph/bargraph_coloc_defense_perspectives_integrase_categories.tsv test/0-39_genomes_Ph25/TAT_defense_results/simulation/bargraph_plot/8.Coloc_defense_with_simulations_integrase_categories

#comparison between TA & defense with MGEs 
#cl: python script/plot_bargraph_coloc_distribution.py test/input/for_bargraph/bargraph_coloc_comparison_TA_defense_with_MGE_categories.tsv test/0-39_genomes_Ph25/TAT_defense_results/simulation/bargraph_plot/9.Coloc_with_specific_MGEs_comparison_TA_defense


def main(tsv_input, outfile):
    """
    """
    
    df_input_files = pd.read_csv(tsv_input, sep = "\t")
    list_cats, list_elem1, list_real_vals, list_sim_vals, list_sim_errs, list_colors, list_pvalues = [], [], [], [], [], [], [] #to produce a single bargraph with all results 
    
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

        n_above = (simulation_means > real_mean).sum()
        n_below = (simulation_means < real_mean).sum()
        n_equal = (simulation_means == real_mean).sum()
        p_value_lower = (simulation_means <= real_mean).sum() / len(simulation_means)
        p_value_higher = (simulation_means >= real_mean).sum() / len(simulation_means)
        standart_deviation = simulation_means.std()
        z_score = (real_mean - simulation_means.mean()) / standart_deviation  
        plot_bargraph(real_mean, simulation_means, p_value_lower, p_value_higher, standart_deviation, row["elem1_type"], row["elem2_type"], row["color"], outfile)        
        
        #Storing data for the complete
        list_cats.append(row["elem2_type"])
        list_real_vals.append(real_mean)
        list_sim_vals.append(simulation_means.mean())
        list_sim_errs.append(standart_deviation)
        list_colors.append(row["color"])
        list_pvalues.append(min(p_value_lower, p_value_higher))
        list_elem1.append(row["elem1_type"])
    
    
    plot_full_bargraph(list_cats, list_real_vals, list_sim_vals, list_sim_errs, list_pvalues, list_colors, list_elem1, outfile)
        


def plot_bargraph(real_mean, simulation_means, p_value_lower, p_value_higher, standart_deviation, elem1_type, elem2_type, color, outfile):


    # Plot: Real vs Simulation (std as error bar)
    fig, ax = plt.subplots(figsize=(3, 4), dpi=300)

    labels = ["Observed", "Simulations"]
    heights = [real_mean, simulation_means.mean()]
    yerr = [0.0, standart_deviation]

    x = np.array([0.1, 0.6])
    bars = ax.bar(
        x, heights, yerr=yerr, capsize=3, width=0.25,
        color=[color, color], edgecolor="black", linewidth=0.7
    )
    
    bars[1].set_hatch('///')
    bars[1].set_edgecolor('black')
    bars[1].set_linewidth(1.0)
    bars[1].set_alpha(0.65)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)

    # Significance bracket using the lower p-value
    p = min(p_value_lower, p_value_higher)
    stars = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'

    # X positions at bar centers
    x1 = bars[0].get_x() + bars[0].get_width()/2
    x2 = bars[1].get_x() + bars[1].get_width()/2

    # Y level above tallest bar (include sim std)
    y_max = max(real_mean, simulation_means.mean() + standart_deviation)
    h = 0.05 * y_max
    y = y_max + h

    # Draw bracket and label
    ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1, c='black')
    ax.text((x1 + x2) / 2, y + h, stars, ha='center', va='bottom')
    ax.set_ylim(top=max(ax.get_ylim()[1], y + h * 1.5))
    ax.set_ylabel(f"Average number of {elem2_type} within the spot given a {elem1_type}")
    #ax.set_title(f"{elem_type} – Real vs Simulation")
    sns.despine(ax=ax)
    plt.tight_layout()
    plt.savefig(outfile + "_distribution_of_"+ elem2_type + "_given_a_" + elem1_type + ".png")
    plt.close(fig)


def plot_full_bargraph(list_cats, list_real_vals, list_sim_vals, list_sim_errs, list_pvalues, list_colors, list_elem_perspective_type, outfile):
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
        
    x = np.arange(len(list_real_vals))
    width = 0.28
    pair_sep = 0.35
    bars_r = ax.bar(x - pair_sep/2, list_real_vals, width, color=list_colors, edgecolor='black', label='Observed average')
    bars_s = ax.bar(x + pair_sep/2, list_sim_vals,  width, yerr=list_sim_errs, capsize=3,
                    color=list_colors, edgecolor='black', alpha=0.65, label='Random simulations',  hatch='///')  
    
    #adding the correct X axis labelling
    n = len(list_real_vals)
    pos = np.ravel(np.column_stack([x - pair_sep/2, x + pair_sep/2]))
    
    #Splitting labels that are too long in two lines
    list_cats = [re.sub(r'\s*&\s*', '\n', c) for c in list_cats]
    
    ax.set_xticks(x)
    ax.set_xticklabels(list_cats)
    fig.subplots_adjust(bottom=0.25)
    
    for br, bs, star in zip(bars_r, bars_s, list_stars):
        label = star
        x1 = br.get_x() + br.get_width()/2
        x2 = bs.get_x() + bs.get_width()/2
        y  = max(br.get_height(), bs.get_height())
        h  = 0.05 * (ax.get_ylim()[1] or 1)
        ax.plot([x1, x1, x2, x2], [y+h, y+2*h, y+2*h, y+h], c='black', lw=1)
        ax.text((x1+x2)/2, y+2*h, label, ha='center', va='bottom') 

    
    if len(list_elem_perspective_type) == 1:
        ax.set_ylabel("Average number of element within the spot\ngiven a " + list_elem_perspective_type[0])
        legend = ax.legend(fontsize = 7.75)
        #changing the color of the legend in grey
        for h in legend.legend_handles:
            if isinstance(h, mpatches.Patch):
                h.set_facecolor('#D3D3D3')
                h.set_alpha(1.0)
        
        plt.tight_layout()
        plt.savefig(outfile+"_All_coloc_results_from_"+ list_elem_perspective_type[0] + "_perspectives.png")
    
    else :
        legend = ax.legend(fontsize = 7.75)
        ax.set_ylabel("Average number of element within the spot\ngiven a specific element")
        handles = legend.legend_handles 
        current_legend_label = [t.get_text() for t in legend.texts]
        list_extra_labels, list_extra_colors = get_associated_label_and_color(list_elem_perspective_type, list_colors)
        list_extra_handles = [mpatches.Patch(facecolor=list_colors[i], edgecolor='black', label=list_extra_labels[i]) for i in range(len(list_extra_labels))]       
        
        legend = ax.legend(handles + list_extra_handles, current_legend_label + list_extra_labels, fontsize=7.75)
        #changing the color of the legend in grey
        for h in legend.legend_handles :
            if isinstance(h, mpatches.Patch) and h.get_label() not in list_extra_labels:
                h.set_facecolor('#D3D3D3')
                h.set_alpha(1.0)
        

        plt.tight_layout()
        plt.savefig(outfile+"_All_coloc_results_from_"+ ("-").join([x for x in set(list_elem_perspective_type)]) + "_perspectives.png")
    
    plt.close(fig)


def get_associated_label_and_color(list_elem_perspective_type, list_colors):
    """
    """
    list_already_assigned = []
    list_extra_labels = []
    list_extra_colors = []
    
    for i in range (len(list_elem_perspective_type)):
        if list_elem_perspective_type[i] not in list_already_assigned:
            list_extra_labels.append("From " + list_elem_perspective_type[i] + " perspectives")
            list_extra_colors.append(list_colors[i])
            list_already_assigned.append(list_elem_perspective_type[i])
    
    return list_extra_labels, list_extra_colors






if __name__ == "__main__":
    import sys
    main(sys.argv[1], sys.argv[2])