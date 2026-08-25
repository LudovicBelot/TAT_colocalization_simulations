#!/usr/bin/env python

import pandas as pd
import numpy as np
from tqdm import tqdm
from matplotlib import pyplot as plt
import seaborn as sns


#cl: python script/plot_colocalization.py test/0-39_genomes_Ph25/TAT_defense_results/simulation/0-Elements_genomic_features.tsv test/0-39_genomes_Ph25/TAT_defense_results/simulation/data/1.1.3_Random_distribution_of_TA_from_TA_perspective.tsv TA TA test/0-39_genomes_Ph25/TAT_defense_results/simulation/new_plot/1.1.3_Random_distribution_of_TA_from_TA_perspective.png
#cl: python script/plot_colocalization.py test/0-39_genomes_Ph25/TAT_defense_results/simulation/0-Elements_genomic_features.tsv test/0-39_genomes_Ph25/TAT_defense_results/simulation/data/1.1.4_Random_distribution_of_defense_from_TA_perspective.tsv TA defense test/0-39_genomes_Ph25/TAT_defense_results/simulation/new_plot/1.1.4_Random_distribution_of_defense_from_TA_perspective.png
#cl: python script/plot_colocalization.py test/0-39_genomes_Ph25/TAT_defense_results/simulation/0-Elements_genomic_features.tsv test/0-39_genomes_Ph25/TAT_defense_results/simulation/data/1.1.6_Random_distribution_of_integrase_from_TA_perspective.tsv TA integrase test/0-39_genomes_Ph25/TAT_defense_results/simulation/new_plot/1.1.6_Random_distribution_of_integrase_from_TA_perspective.png
#cl: python script/plot_colocalization.py test/0-39_genomes_Ph25/TAT_defense_results/simulation/0-Elements_genomic_features.tsv test/0-39_genomes_Ph25/TAT_defense_results/simulation/data/1.1.5_Random_distribution_of_IS_from_TA_perspective.tsv TA IS test/0-39_genomes_Ph25/TAT_defense_results/simulation/new_plot/1.1.5_Random_distribution_of_IS_from_TA_perspective.png

#cl: python script/plot_colocalization.py test/0-39_genomes_Ph25/TAT_defense_results/simulation/0-Elements_genomic_features.tsv test/0-39_genomes_Ph25/TAT_defense_results/simulation/data/2.1.3_Random_distribution_of_TA_from_defense_perspective.tsv defense TA test/0-39_genomes_Ph25/TAT_defense_results/simulation/new_plot/2.1.3_Random_distribution_of_TA_from_defense_perspective.png
#cl: python script/plot_colocalization.py test/0-39_genomes_Ph25/TAT_defense_results/simulation/0-Elements_genomic_features.tsv test/0-39_genomes_Ph25/TAT_defense_results/simulation/data/2.1.4_Random_distribution_of_defense_from_defense_perspective.tsv defense defense test/0-39_genomes_Ph25/TAT_defense_results/simulation/new_plot/2.1.4_Random_distribution_of_defense_from_defense_perspective.png
#cl: python script/plot_colocalization.py test/0-39_genomes_Ph25/TAT_defense_results/simulation/0-Elements_genomic_features.tsv test/0-39_genomes_Ph25/TAT_defense_results/simulation/data/2.1.5_Random_distribution_of_IS_from_defense_perspective.tsv defense IS test/0-39_genomes_Ph25/TAT_defense_results/simulation/new_plot/2.1.5_Random_distribution_of_IS_from_defense_perspective.png
#cl: python script/plot_colocalization.py test/0-39_genomes_Ph25/TAT_defense_results/simulation/0-Elements_genomic_features.tsv test/0-39_genomes_Ph25/TAT_defense_results/simulation/data/2.1.6_Random_distribution_of_integrase_from_defense_perspective.tsv defense integrase test/0-39_genomes_Ph25/TAT_defense_results/simulation/new_plot/2.1.6_Random_distribution_of_integrase_from_defense_perspective.png

#cl: python script/plot_colocalization.py test/0-39_genomes_Ph25/TAT_defense_results/simulation/0-Elements_genomic_features.tsv test/0-39_genomes_Ph25/TAT_defense_results/simulation/data/3.1.3_Random_distribution_of_TA_from_IS_perspective.tsv IS TA test/0-39_genomes_Ph25/TAT_defense_results/simulation/new_plot/3.1.3_Random_distribution_of_TA_from_IS_perspective.png
#cl: python script/plot_colocalization.py test/0-39_genomes_Ph25/TAT_defense_results/simulation/0-Elements_genomic_features.tsv test/0-39_genomes_Ph25/TAT_defense_results/simulation/data/3.1.4_Random_distribution_of_defense_from_IS_perspective.tsv IS defense test/0-39_genomes_Ph25/TAT_defense_results/simulation/new_plot/3.1.4_Random_distribution_of_defense_from_IS_perspective.png
#cl: python script/plot_colocalization.py test/0-39_genomes_Ph25/TAT_defense_results/simulation/0-Elements_genomic_features.tsv test/0-39_genomes_Ph25/TAT_defense_results/simulation/data/3.1.5_Random_distribution_of_IS_from_IS_perspective.tsv IS IS test/0-39_genomes_Ph25/TAT_defense_results/simulation/new_plot/3.1.5_Random_distribution_of_IS_from_IS_perspective.png
#cl: python script/plot_colocalization.py test/0-39_genomes_Ph25/TAT_defense_results/simulation/0-Elements_genomic_features.tsv test/0-39_genomes_Ph25/TAT_defense_results/simulation/data/3.1.6_Random_distribution_of_integrase_from_IS_perspective.tsv IS integrase test/0-39_genomes_Ph25/TAT_defense_results/simulation/new_plot/3.1.6_Random_distribution_of_integrase_from_IS_perspective.png

#cl: python script/plot_colocalization.py test/0-39_genomes_Ph25/TAT_defense_results/simulation/0-Elements_genomic_features.tsv test/0-39_genomes_Ph25/TAT_defense_results/simulation/data/4.1.3_Random_distribution_of_TA_from_integrase_perspective.tsv integrase TA test/0-39_genomes_Ph25/TAT_defense_results/simulation/new_plot/14.1.3_Random_distribution_of_TA_from_integrase_perspective.png
#cl: python script/plot_colocalization.py test/0-39_genomes_Ph25/TAT_defense_results/simulation/0-Elements_genomic_features.tsv test/0-39_genomes_Ph25/TAT_defense_results/simulation/data/4.1.4_Random_distribution_of_defense_from_integrase_perspective.tsv integrase defense test/0-39_genomes_Ph25/TAT_defense_results/simulation/new_plot/4.1.4_Random_distribution_of_defense_from_integrase_perspective.png
#cl: python script/plot_colocalization.py test/0-39_genomes_Ph25/TAT_defense_results/simulation/0-Elements_genomic_features.tsv test/0-39_genomes_Ph25/TAT_defense_results/simulation/data/4.1.5_Random_distribution_of_IS_from_integrase_perspective.tsv integrase IS test/0-39_genomes_Ph25/TAT_defense_results/simulation/new_plot/4.1.5_Random_distribution_of_IS_from_integrase_perspective.png
#cl: python script/plot_colocalization.py test/0-39_genomes_Ph25/TAT_defense_results/simulation/0-Elements_genomic_features.tsv test/0-39_genomes_Ph25/TAT_defense_results/simulation/data/4.1.6_Random_distribution_of_integrase_from_integrase_perspective.tsv integrase integrase test/0-39_genomes_Ph25/TAT_defense_results/simulation/new_plot/4.1.6_Random_distribution_of_integrase_from_integrase_perspective.png



def main(tsv_real, tsv_sim, elem1_type, elem2_type, outfile):
    """
    """

    df_real = pd.read_csv(tsv_real, sep = "\t")
    df_sim = pd.read_csv(tsv_sim, sep = "\t", index_col = 0)
    df_real = df_real[df_real["Element_type"] == elem1_type]
    real_data = df_real[f"N_{elem2_type}_colocalizing"]
    
    # 1. Moyennes des simulations
    simulation_means = df_sim.mean(axis=0)

    # 2. Moyenne des données réelles
    # Adapter cette ligne si ta colonne a un nom spécifique
    real_data = real_data.squeeze()
    real_mean = real_data.mean()

    n_above = (simulation_means > real_mean).sum()
    n_below = (simulation_means < real_mean).sum()
    n_equal = (simulation_means == real_mean).sum()
    p_value_lower = (simulation_means <= real_mean).sum() / len(simulation_means)
    p_value_higher = (simulation_means >= real_mean).sum() / len(simulation_means)
    z_score = (real_mean - simulation_means.mean()) / simulation_means.std()
    
    # 3. Plot
    plt.figure(figsize=(10, 6))

    # Histogramme ou KDE des moyennes simulées
    if elem1_type == "TA":
        color_sim = "blue"
    elif elem1_type == "defense":
        color_sim = "#a00000"
    elif elem1_type == "IS":
        color_sim = "#00c000"
    elif elem1_type == "integrase":
        color_sim = "#ad00e3"
    
    #Uncomment here for bargraph + kde
    #sns.histplot(simulation_means, kde=True, bins=30, color=color_sim, alpha=0.6, label="Mean simulation distribution")

    #Kde + area filled
    #sns.kdeplot(simulation_means, fill=True, color=color_sim, alpha=0.6, label="Mean simulation distribution")

    #Kde + area filled with the N_simulations axis kept instead of the density probability
    # First plot histogram to get the count scale
    counts, bins, patches = plt.hist(simulation_means, bins=30, color=None, alpha=0, label="Mean simulation distribution")
    ax1 = plt.gca()
    ax1.set_ylabel("N simulation", fontsize=14)
    ax1.set_xlabel(f"Mean number of {elem2_type} associated with a given {elem1_type}", fontsize=14)
    ax1.tick_params(axis='both', which='major', labelsize=14)
    
    # Then overlay a KDE if you want the smooth curve
    ax2 = plt.gca().twinx()
    sns.kdeplot(simulation_means, ax=ax2, color=color_sim, alpha=0.3, fill=True)
    ax2.set_ylabel('')
    ax2.tick_params(right=False, labelright=False)
    ax2.spines['right'].set_visible(False)
    
    # Ligne verticale pour la moyenne réelle
    #plt.axvline(real_mean, color="red", linestyle="-", linewidth=4, label=f"Real mean = {real_mean:.4f}")
    plt.axvline(real_mean, color="red", linestyle="--", linewidth=3, label=f"Real mean = {real_mean:.4f}")
    
    # Texte indicatif
    """
    plt.text(0.95, 0.8,
             f"N simulation ≥ real : {n_above+n_equal} ({round(n_above/len(simulation_means)*100, 4)}%)\n"
             f"Unilateral p-value ≥ real : {p_value_higher}\n"
             f"N simulations ≤ real : {n_below + n_equal} ({round((n_below + n_equal)/len(simulation_means)*100, 4)}%)\n"
             f"Unilateral p-value ≤ real : {p_value_lower}\n"
             f"Z-score : {z_score:.4f}",
             transform=plt.gca().transAxes,
             fontsize=10,
             verticalalignment='top',
             horizontalalignment='right',
             bbox=dict(boxstyle="round", facecolor="white", edgecolor="gray", alpha=0.8))
    """

    # Remove the box around the plot
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Mise en forme
    #plt.xlabel(f"Mean number of {elem2_type} associated with a given {elem1_type}", fontsize=14)
    #plt.ylabel("N simulation", fontsize=14)
    plt.tick_params(axis='both', which='major', labelsize=14)
    #plt.title(f"{elem1_type} colocalization perspectives with (or not) randomly distributed {elem2_type}")
    #plt.legend()
    plt.grid(False)
    plt.tight_layout()
    plt.savefig(outfile)
    plt.close()





if __name__ == "__main__":
    import sys
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])