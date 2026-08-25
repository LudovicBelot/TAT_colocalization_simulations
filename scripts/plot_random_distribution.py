#!/usr/bin/env python

import pandas as pd
import numpy as np
from tqdm import tqdm
from matplotlib import pyplot as plt
import seaborn as sns

#cl: python script/plot_random_distribution.py test/0-39_genomes_Ph25/TAT_defense_results/simulation_v2_from_type_perspectives/0-Elements_genomic_features.tsv test/0-39_genomes_Ph25/TAT_defense_results/simulation_v2_from_type_perspectives/data/1.1.2-Random_spot_distribution_TA_gene_numbers_within_spot.tsv TA test/0-39_genomes_Ph25/TAT_defense_results/simulation_v2_from_type_perspectives/new_plot/1.1.2-Random_spot_distribution_TA_gene_numbers_within_spot.png
#cl: python script/plot_random_distribution.py test/0-39_genomes_Ph25/TAT_defense_results/simulation_v2_from_type_perspectives/0-Elements_genomic_features.tsv test/0-39_genomes_Ph25/TAT_defense_results/simulation_v2_from_type_perspectives/data/2.1.2-Random_spot_distribution_defense_gene_numbers_within_spot.tsv defense test/0-39_genomes_Ph25/TAT_defense_results/simulation_v2_from_type_perspectives/new_plot/1.0.2-Random_spot_distribution_defense_gene_numbers_within_spot.png
#cl: python script/plot_random_distribution.py test/0-39_genomes_Ph25/TAT_defense_results/simulation_v2_from_type_perspectives/0-Elements_genomic_features.tsv test/0-39_genomes_Ph25/TAT_defense_results/simulation_v2_from_type_perspectives/data/3.1.2-Random_spot_distribution_IS_gene_numbers_within_spot.tsv IS test/0-39_genomes_Ph25/TAT_defense_results/simulation_v2_from_type_perspectives/new_plot/1.0.3-Random_spot_distribution_IS_gene_numbers_within_spot.png
#cl: python script/plot_random_distribution.py test/0-39_genomes_Ph25/TAT_defense_results/simulation_v2_from_type_perspectives/0-Elements_genomic_features.tsv test/0-39_genomes_Ph25/TAT_defense_results/simulation_v2_from_type_perspectives/data/4.1.2-Random_spot_distribution_integrase_gene_numbers_within_spot.tsv integrase test/0-39_genomes_Ph25/TAT_defense_results/simulation_v2_from_type_perspectives/new_plot/1.0.4-Random_spot_distribution_integrase_gene_numbers_within_spot.png

def main(tsv_real, tsv_sim, elem_type, outfile):
    """
    """
       
    df_real = pd.read_csv(tsv_real, sep = "\t")
    df_sim = pd.read_csv(tsv_sim, sep = "\t", index_col = 0)
    df_real = df_real[df_real["Element_type"] == elem_type]
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
    z_score = (real_mean - simulation_means.mean()) / simulation_means.std()    
    
    # 3. Plot
    plt.figure(figsize=(10, 6))

    # Histogramme ou KDE des moyennes simulées
    if elem_type == "TA":
        color_sim = "blue"
    elif elem_type == "defense":
        color_sim = "#a00000"
    elif elem_type == "IS":
        color_sim = "#00c000"
    elif elem_type == "integrase":
        color_sim = "#ad00e3"
  
    sns.histplot(simulation_means, kde=True, bins=30, color=color_sim, alpha=0.6, label="Moyennes (simulations)")

    # Ligne verticale pour la moyenne réelle
    plt.axvline(real_mean, color="red", linestyle="--", linewidth=2, label=f"Real mean = {real_mean:.4f}")

    # Texte indicatif
    plt.text(0.95, 0.8,
             f"Simulations ≥ real : {n_above+n_equal} ({round(n_above/len(simulation_means)*100, 4)}%)\n"
             f"Unilateral p-value ≥ real : {p_value_higher}\n"
             f"Simulations ≤ real : {n_below + n_equal} ({round((n_below + n_equal)/len(simulation_means)*100, 4)}%)\n"
             f"Unilateral p-value ≤ real : {p_value_lower}\n"
             f"Z-score : {z_score:.4f}",
             transform=plt.gca().transAxes,
             fontsize=10,
             verticalalignment='top',
             horizontalalignment='right',
             bbox=dict(boxstyle="round", facecolor="white", edgecolor="gray", alpha=0.8))
    
    
    # Mise en forme
    plt.xlabel("Length of the spots (in number of non-redundant genes)")
    plt.ylabel("N simulations")
    plt.title(f"{elem_type} random distribution within Photorhabdus Spot-pangenome")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(outfile)
    plt.close()



if __name__ == "__main__":
    import sys
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])