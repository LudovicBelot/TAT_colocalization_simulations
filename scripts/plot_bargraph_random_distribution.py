#!/usr/bin/env python

import pandas as pd
import numpy as np
from tqdm import tqdm
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

#cl: python script/plot_bargraph_random_distribution.py test/input/for_bargraph/bargraph_random.tsv test/0-39_genomes_Ph25/TAT_defense_results/simulation_v2_from_type_perspectives/new_plot/1.1.2-Random_spot_distribution_TA_gene_numbers_within_spot_bargraph


def main(tsv_input, outfile):
    """
    """
    
    df_input_files = pd.read_csv(tsv_input, sep = "\t")
    list_cats, list_OE, list_sim_errs, list_colors, list_pvalues = [], [], [], [], [] #to produce a single bargraph with all results 
    
    for _, row in df_input_files.iterrows():
            
        #Loading data
        df_real = pd.read_csv(row["real_file"], sep = "\t")
        df_sim = pd.read_csv(row["sim_file"], sep = "\t", index_col = 0)
        df_real = df_real[df_real["Element_type"] == row["elem_type"]]
        real_data = df_real[f"N_NR_genes_in_spot"]
        
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
        #plot_bargraph(real_mean, simulation_means, p_value_lower, p_value_higher, standart_deviation, row["elem_type"], row["color"], outfile)        
        
        #Storing data for the complete
        list_cats.append(row["elem_type"])
        list_real_vals.append(real_mean)
        list_sim_vals.append(simulation_means.mean())
        list_sim_errs.append(standart_deviation)
        list_colors.append(row["color"])
        list_pvalues.append(min(p_value_lower, p_value_higher))
    
    plot_full_bargraph(list_cats, list_real_vals, list_sim_vals, list_sim_errs, list_pvalues, list_colors, outfile)
        


def plot_bargraph(real_mean, simulation_means, p_value_lower, p_value_higher, standart_deviation, group_name, color, outfile):


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
    ax.set_ylabel(f"Average number of non-redundant\ngenes within the spot\ngiven a {group_name} system")
    #ax.set_title(f"{elem_type} – Real vs Simulation")
    sns.despine(ax=ax)
    plt.tight_layout()
    plt.savefig(outfile + "_"+ group_name + ".png")
    plt.close(fig)


def plot_full_bargraph(list_cats, list_real_vals, list_sim_vals, list_sim_errs, list_pvalues, list_colors, outfile):
    """
    """
    #Pvalues_stars 
    list_stars = []
    for p in list_pvalues:
        list_stars.append('***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns')

    fig, ax = plt.subplots(figsize=(6.25, 3.5), dpi=150)
    x = np.arange(len(list_real_vals))
    width = 0.28
    pair_sep = 0.35
    bars_r = ax.bar(x - pair_sep/2, list_real_vals, width, color=list_colors, edgecolor='black', label='Observed average')
    bars_s = ax.bar(x + pair_sep/2, list_sim_vals,  width, yerr=list_sim_errs, capsize=3,
                    color=list_colors, edgecolor='black', alpha=0.65, label='Random simulations',  hatch='///')  
    
    #adding the correct X axis labelling
    n = len(list_real_vals)
    pos = np.ravel(np.column_stack([x - pair_sep/2, x + pair_sep/2]))
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

    ax.set_ylabel("Average number of non-redundant genes\nwithin the spot given an element")
    legend = ax.legend(fontsize = 7.75)
    #changing the color of the legend in grey
    for h in legend.legend_handles:
        if isinstance(h, mpatches.Patch):
            h.set_facecolor('#D3D3D3')
            h.set_alpha(1.0)
            
    plt.tight_layout()
    plt.savefig(outfile+"_All_results.png")
    plt.close(fig)





if __name__ == "__main__":
    import sys
    main(sys.argv[1], sys.argv[2])