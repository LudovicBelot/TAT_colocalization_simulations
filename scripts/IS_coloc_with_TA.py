#!/usr/bin/env python

import pandas as pd
from tqdm import tqdm
import random
import copy
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline
import seaborn as sns

"""
Standalone script which aims to plot the colocalization of TA systems with each IS family without rerunning the whole TAT_defense.py script
"""



def main(outdir):
    """
    """
    
    df_spot_info = pd.read_csv(f"{outdir}/TAT_defense_results/core_mapping/3-all_spots_features.tsv", sep = "\t")
    df_all_elements = pd.read_csv(f"{outdir}/TAT_defense_results/simulation_TAs_perspectives_all_IS_families/0-Elements_genomic_features.tsv", sep = "\t")
    list_IS_families_to_simulate = df_all_elements[df_all_elements["Element_type"] == "IS"]["Element_family"].unique().tolist()
    #df_all_elements = fetch_coloc_features(df_all_elements, list_IS_families_to_simulate)
    #df_all_elements.to_csv(f"{outdir}/TAT_defense_results/simulation_TAs_perspectives_all_IS_families/0.1-Elements_genomic_features_all_IS_families.tsv", sep = "\t", index = False)
    df_all_elements = pd.read_csv(f"{outdir}/TAT_defense_results/simulation_TAs_perspectives_all_IS_families/0.1-Elements_genomic_features_all_IS_families.tsv", sep = "\t")
    
    init_global_res_file(outdir)
    n=0
    for elem_family in tqdm(list_IS_families_to_simulate, desc = "Running the simulations for the IS genes (either all, or based on their family)", total = len(list_IS_families_to_simulate)):
        n += 1
        """
        df_simulations_IS = run_simulation(df_spot_info, df_all_elements[df_all_elements["Element_family"] == elem_family], n_simulations = 10000)
        df_simulations_IS.to_csv(f"{outdir}/TAT_defense_results/simulation/data/3.{n}.1-Random_spot_distribution_{elem_family}.tsv", sep = "\t")
        df_elem_family_with_TA = convert_spot_to_number_of_elements_simulated(df_simulations_IS, df_all_elements[df_all_elements["Element_type"] == "TA"], elem_family)
        df_elem_family_with_TA.to_csv(f"{outdir}/TAT_defense_results/simulation/data/3.{n}.2-Colocalization_perspectives_from_TA_point_of_view_with_{elem_family}.tsv", sep = "\t")
        """
        #Remove the next line and uncomment before
        df_elem_family_with_TA = pd.read_csv(f"{outdir}/TAT_defense_results/simulation_TAs_perspectives_all_IS_families/data/3.{n}.2-Colocalization_perspectives_from_TA_point_of_view_with_{elem_family}.tsv", sep = "\t", index_col= 0)
        tmp_global_res = plot_mean_distribution_from_perspective(df_elem_family_with_TA,
                                                                df_all_elements[df_all_elements["Element_type"] == "TA"],
                                                                elem_family,
                                                                "TA",
                                                                f"{outdir}/TAT_defense_results/simulation_TAs_perspectives_all_IS_families/plot/3.{n}.3_Random_distribution_of_{elem_family}_from_TA_perspective.png")
        
        tmp_global_res.insert(1, str(len(df_all_elements[df_all_elements["Element_family"] == elem_family])))
        write_global_res(tmp_global_res, outdir)



def fetch_coloc_features(df_element, list_elem_families):
    """
    Function which add a column representing the number of element retrieved within the spot of the studied element.
    One column per element tested.
    Note that, when establishing the number of elements of the same type than the one studied (e.g. looking at the number of TA colocalizing with a TA), we substract one.
    So we do not count the element itself.
    """
    
    for family in tqdm(list_elem_families, desc = "Fetching coloc features for each element", position= 1, leave = False):
        tmp_list_res = []
        for row in df_element.iterrows():
            if row[1]["Element_family"] == family: #then we need to remove 1 to avoid counting the studied element itself
                tmp_list_res.append(len(df_element[(df_element["Element_family"] == family) & (df_element["Spot_number"] == row[1]["Spot_number"])]) - 1)
            else : 
                tmp_list_res.append(len(df_element[(df_element["Element_family"] == family) & (df_element["Spot_number"] == row[1]["Spot_number"])]))
        
        df_element[f"N_{family}_colocalizing"] = tmp_list_res
                    
    return df_element



def run_simulation(df_spot_info, df_element_to_test, n_simulations = 10000):
    """
    Function that will perform the simulation
    Parameters
    ----------
    df_spot_info : pandas dataframe, dataframe with the spot information
    df_element_to_test : pandas dataframe, dataframe with the search results for one of the element to randomly distribute
    n_simulations : int, number of simulations to perform
    """

    # create a dict with the spot informations so we can perform the distribution simulation in the given spot-pangenome
    d_spot_info = {} #k: lispot_number, v: number of non redundant genes in the spot
    list_spot_number = [] #to ease operations with random distribution
    
    for row in df_spot_info.iterrows():
        d_spot_info[row[1]["interval_number"]] = int(row[1]["N_NR_genes"])
        
        for _ in range(int(row[1]["N_NR_genes"])): # A large spot (in number of genes), should be more likely to have a specific genes if randomly distributed
            list_spot_number.append(row[1]["interval_number"])
    
    # Getting the informations about the elements to distribute:
    # i.e. how many elements to distribute, and the length of each element (in term of number of genes)
    # e.g. a TA is composed most of the time of 5 genes
    
    d_elem_info = {} #k : element_name, v : length of the element
    for row in df_element_to_test.iterrows():
        d_elem_info[row[1]["Element_name"]] = int(row[1]["Element_length"])
    
    d_res_simulation = {} #k : element randomly distributed, v d{simulation number : spot number in which the element was distributed}
    for elem in d_elem_info.keys():
        d_res_simulation[elem] = {}    
    
    for i in tqdm(range(n_simulations), desc = "Performing the random distribution within the spot-pangenome", total = n_simulations, position = 1, leave = False):
        
        # initialing the dictionnary
        d_current_simulation = copy.deepcopy(d_spot_info)
        
        for element_name, element_length in d_elem_info.items():
            
            # check if the spot is large enough to host the element
            is_spot_large_enough = False
            while not is_spot_large_enough:
                
                # randomly select a spot
                random_spot = random.choice(list_spot_number)
                
                # if the spot is large enough, we assign the element to it and soustract the length of the element from the spot
                if d_current_simulation[random_spot] >= element_length:    
                    d_res_simulation[element_name][i] = int(random_spot)
                    d_current_simulation[random_spot] -= element_length
                    is_spot_large_enough = True

    
    df_res_simulation = pd.DataFrame.from_dict(d_res_simulation, orient = "index")
    
    return df_res_simulation




def convert_spot_to_number_of_elements_simulated(df_sim, df_all_elements, elem_of_interest):
    """
    Estimating the number of randomly distributed elem_of_interest that would colocalize with each TA systems so we can plot it later
    """
    
    d_res = {}
    list_simulation = df_sim.columns.tolist()

    list_index_real = df_all_elements.index.tolist()
    for index in tqdm(list_index_real, total = len(list_index_real), desc = f"Preparing the colocalization data for plotting the distribution between TAs and {elem_of_interest}", position = 1, leave = False):
        d_res[index] = {}
        for sim in tqdm(list_simulation, total = len(list_simulation), desc = "Doing it for each simulation", position = 2, leave = False):
            d_res[index][sim] = int(df_sim[sim].value_counts().get(df_all_elements.loc[index, "Spot_number"], 0))
    
    return pd.DataFrame.from_dict(d_res, orient = "index")




def plot_mean_distribution_from_perspective(df_sim, df_real, family_colocalizing, elem_perspective_name, outname):
    
    """
    Plot the mean distribution of the real data and simulations from a same element perspectives.
    """
    
    real_data = df_real[f"N_{family_colocalizing}_colocalizing"]
    
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
    # Calculate Z-score
    sim_std = simulation_means.std()
    z_score = (real_mean - simulation_means.mean()) / sim_std
    
    # 3. Plot
    plt.figure(figsize=(10, 6))

    # Histogramme ou KDE des moyennes simulées
    sns.histplot(simulation_means, kde=True, bins="auto", color="blue", alpha=0.6, label="Mean distribution (simulations)")

    # Ligne verticale pour la moyenne réelle
    plt.axvline(real_mean, color="red", linestyle="--", linewidth=2, label=f"Real mean = {real_mean:.4f}")
    
    #Setting the limits of the x axis
    x_min = min(simulation_means.min(), real_mean)
    x_max = max(simulation_means.max(), real_mean)
    padding = (x_max - x_min) * 0.05  # 5% padding
    plt.xlim(0, 2) #Fixed limits so it is easier to compare them
    
    # Texte indicatif
    plt.text(0.95, 0.8,
             f"Simulations ≥  real : {n_above} ({round((n_above + n_equal)/len(simulation_means)*100, 4)}%)\n"
             f"Unilateral p-value ≥  real : {p_value_higher}\n"
             f"Simulations ≤ real : {n_below + n_equal} ({round((n_below + n_equal)/len(simulation_means)*100, 4)}%)\n"
             f"Unilateral p-value ≤ real : {p_value_lower}\n"
             "Z-score : {:.4f}".format(z_score),
             transform=plt.gca().transAxes,
             fontsize=10,
             verticalalignment='top',
             horizontalalignment='right',
             bbox=dict(boxstyle="round", facecolor="white", edgecolor="gray", alpha=0.8))
    
    
    # Mise en forme
    plt.xlabel("Mean element number from TA perspectives")
    plt.ylabel("N simulations")
    plt.title(f"{elem_perspective_name} colocalization perspectives with (or not) randomly distributed {family_colocalizing}")
    plt.legend(loc='upper right')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(outname)
    plt.close()
    
    return [family_colocalizing, str(real_mean), str(n_below), str(n_equal), str(n_above), str(round(p_value_lower,4)), str(round(p_value_higher,4)), str(round(z_score,4))]



def init_global_res_file(outdir):
    """
    Writing the head of the final res file (tsv file that will be updated for each IS family)
    """

    with open(outdir + "/TAT_defense_results/simulation_TAs_perspectives_all_IS_families/data/3.0-All_IS_families_colocalization_from_TA_perspectives.tsv", "w") as f:
        f.write("IS_family\tN_members_in_family\treal_mean_value\tn_below_real_mean_value\tn_equal_real_mean_value\tn_over_real_mean_value\tp-value_below_equal_real\tp-value_over_equal_real\tZ-score\n")
    

def write_global_res(res2add, outdir):
    """
    Adding one result line to the global file
    """
    
    with open(outdir + "/TAT_defense_results/simulation_TAs_perspectives_all_IS_families/data/3.0-All_IS_families_colocalization_from_TA_perspectives.tsv", "a") as f:
        f.write("\t".join(res2add) + "\n")






if __name__ == "__main__":
    import sys
    main(sys.argv[1])