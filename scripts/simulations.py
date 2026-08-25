#!/usr/bin/env python

import pandas as pd
import copy
import random
from tqdm import tqdm
import os
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import make_interp_spline
import seaborn as sns


def main(outdir, to_test = ["TA", "defense", "IS", "integrase"], to_consider = "all", n_cpu = 1, n_simulations = 10000):
    """
    Mais function that will performs random distribution of TAs, defense systems, integrase genes and IS transposase genes.
    Parameters
    ----------
    indir : str, path to the directory containing the genomes
    to_test: list of str, list of elements to test (TA, defense, integrase, IS transposase)
    outdir : str, path to the directory where the results will be saved
    to_consider : str, either "all" or "per_family". "All" will consider all "TAs" as a whole while "per_family" will separate them by their family before proceedding to the simulations (same operations with the other types) 
    n_cpu : int, number of cpu to use
    
    """
    
    # load the tsv file with the spot informations (number of non redundant genes per spot)
    df_spot_info = pd.read_csv(f"{outdir}/TAT_defense_results/core_mapping/3-all_spots_features.tsv", sep = "\t")
    
    if os.path.exists(f"{outdir}/TAT_defense_results/simulation/0-Elements_genomic_features.tsv"):
        df_all_elements = pd.read_csv(f"{outdir}/TAT_defense_results/simulation/0-Elements_genomic_features.tsv", sep = "\t")
    
    else :
        df_all_elements = fetch_elements_features(outdir, to_test, df_spot_info)
        df_all_elements = fetch_coloc_features(df_all_elements, outdir, to_test)

    
    #Now running the simulations
    d_df_res_simulations = {} #k : string "TA"/"defense"/"IS"/"integrase", v: dataframe of the associated element randomly distributed
    for simulation_to_run in to_test:
        
        if simulation_to_run == "TA":
            
            list_element_type_to_simulate = []
            if to_consider == "all":
                list_element_type_to_simulate = ["TA"]
                family_or_type = "Element_type"
            elif to_consider == "per_family":
                list_element_type_to_simulate = df_all_elements[df_all_elements["Element_type"] == "TA"]["Element_family"].drop_duplicates().tolist()
                family_or_type = "Element_family"
            
            print("Initializing the random distribution of the TA systems")
            n=0
            for elem in tqdm(list_element_type_to_simulate, desc = "Running the simulations for the TAs (either all, or based on their family)", total = len(list_element_type_to_simulate)):
                n += 1
                df_simulations_TAs = run_simulation(df_spot_info, df_all_elements[df_all_elements[family_or_type] == elem], n_simulations = n_simulations)
                d_df_res_simulations["TA"] = df_simulations_TAs
                df_simulations_TAs.to_csv(f"{outdir}/TAT_defense_results/simulation/data/1.{n}.1-Random_spot_distribution_{elem}.tsv", sep = "\t")
                df_simulations_TAs_gene_numbers = convert_spot_to_number_of_genes(df_simulations_TAs, df_spot_info)
                df_simulations_TAs_gene_numbers.to_csv(f"{outdir}/TAT_defense_results/simulation/data/1.{n}.2-Random_spot_distribution_{elem}_gene_numbers_within_spot.tsv", sep = "\t")
                plot_CDF_distribution(f"{outdir}/TAT_defense_results/simulation/plot/1.{n}.1-Random_spot_distribution_{elem}_plot.png", df_simulations_TAs_gene_numbers, df_all_elements[df_all_elements[family_or_type] == elem])
                plot_mean_distribution(f"{outdir}/TAT_defense_results/simulation/plot/1.{n}.2-Random_spot_mean_distribution_{elem}_plot.png", df_simulations_TAs_gene_numbers, df_all_elements[df_all_elements[family_or_type] == elem])
                
                """
                #Estimating the % of TA colocalizing with each element and comparing them with 
                d_simulations_colocalization_with_TA = convert_spot_to_number_of_elements(df_simulations_TAs, outdir, "TAs")
                n_2 = 2
                for elem2, associated_df in  d_simulations_colocalization_with_TA.items():
                    n_2 += 1
                    associated_df.to_csv(f"{outdir}/TAT_defense_results/simulation/data/1.{n}.{n_2}-Random_distribution_of_TA_coloc_with_{elem2}_number_per_spot.tsv", sep = "\t")
                    plot_coloc_mean_distribution(f"{outdir}/TAT_defense_results/simulation/plot/1.{n}.{n_2}-Random_spot_mean_distribution_TA_coloc_with_{elem2}_plot.png", associated_df, df_all_elements[df_all_elements[family_or_type] == elem], elem2)   
                """ 
                    
        elif simulation_to_run == "defense":
            
            list_element_type_to_simulate = []
            if to_consider == "all":
                list_element_type_to_simulate = ["defense"]
                family_or_type = "Element_type"
            elif to_consider == "per_family":
                list_element_type_to_simulate = df_all_elements[df_all_elements["Element_type"] == "defense"]["Element_family"].drop_duplicates().tolist()
                family_or_type = "Element_family"
            
            print("Initializing the random distribution of the defense systems")
            n=0
            for elem in tqdm(list_element_type_to_simulate, desc = "Running the simulations for the defense systems (either all, or based on their family)", total = len(list_element_type_to_simulate)):
                n += 1
                df_simulations_defense = run_simulation(df_spot_info, df_all_elements[df_all_elements[family_or_type] == elem], n_simulations = n_simulations)
                d_df_res_simulations["defense"] = df_simulations_defense
                df_simulations_defense.to_csv(f"{outdir}/TAT_defense_results/simulation/data/2.{n}.1-Random_spot_distribution_{elem}.tsv", sep = "\t")
                df_simulations_defenses_gene_numbers = convert_spot_to_number_of_genes(df_simulations_defense, df_spot_info)
                df_simulations_defenses_gene_numbers.to_csv(f"{outdir}/TAT_defense_results/simulation/data/2.{n}.2-Random_spot_distribution_{elem}_gene_numbers_within_spot.tsv", sep = "\t")
                plot_CDF_distribution(f"{outdir}/TAT_defense_results/simulation/plot/2.{n}.1-Random_spot_distribution_{elem}_plot.png", df_simulations_defenses_gene_numbers, df_all_elements[df_all_elements[family_or_type] == elem])
                plot_mean_distribution(f"{outdir}/TAT_defense_results/simulation/plot/2.{n}.2-Random_spot_mean_distribution_{elem}_plot.png", df_simulations_defenses_gene_numbers, df_all_elements[df_all_elements[family_or_type] == elem])

                """
                #Estimating the % of defense systems colocalizing with each element and comparing them with 
                d_simulations_colocalization_with_defense = convert_spot_to_number_of_elements(df_simulations_defense, outdir, "defense")
                n_2 = 2
                for elem2, associated_df in  d_simulations_colocalization_with_defense.items():
                    n_2 += 1
                    associated_df.to_csv(f"{outdir}/TAT_defense_results/simulation/data/2.{n}.{n_2}-Random_distribution_of_defense_coloc_with_{elem2}_number_per_spot.tsv", sep = "\t")
                    plot_coloc_mean_distribution(f"{outdir}/TAT_defense_results/simulation/plot/2.{n}.{n_2}-Random_spot_mean_distribution_defense_coloc_with_{elem2}_plot.png", associated_df, df_all_elements[df_all_elements[family_or_type] == elem], elem2)   
                """


        elif simulation_to_run == "IS":
            
            list_element_type_to_simulate = []
            if to_consider == "all":
                list_element_type_to_simulate = ["IS"]
                family_or_type = "Element_type"
            elif to_consider == "per_family":
                list_element_type_to_simulate = df_all_elements[df_all_elements["Element_type"] == "IS"]["Element_family"].drop_duplicates().tolist()
                family_or_type = "Element_family"
            
            print("Initializing the random distribution of the IS transposase genes")
            n=0
            for elem in tqdm(list_element_type_to_simulate, desc = "Running the simulations for the IS genes (either all, or based on their family)", total = len(list_element_type_to_simulate)):
                n += 1
                df_simulations_IS = run_simulation(df_spot_info, df_all_elements[df_all_elements[family_or_type] == elem], n_simulations = n_simulations)
                d_df_res_simulations["IS"] = df_simulations_IS
                df_simulations_IS.to_csv(f"{outdir}/TAT_defense_results/simulation/data/3.{n}.1-Random_spot_distribution_{elem}.tsv", sep = "\t")
                df_simulations_IS_gene_numbers = convert_spot_to_number_of_genes(df_simulations_IS, df_spot_info)
                df_simulations_IS_gene_numbers.to_csv(f"{outdir}/TAT_defense_results/simulation/data/3.{n}.2-Random_spot_distribution_{elem}_gene_numbers_within_spot.tsv", sep = "\t")
                plot_CDF_distribution(f"{outdir}/TAT_defense_results/simulation/plot/3.{n}.1-Random_spot_distribution_{elem}_plot.png", df_simulations_IS_gene_numbers, df_all_elements[df_all_elements[family_or_type] == elem])
                plot_mean_distribution(f"{outdir}/TAT_defense_results/simulation/plot/3.{n}.2-Random_spot_mean_distribution_{elem}_plot.png", df_simulations_IS_gene_numbers, df_all_elements[df_all_elements[family_or_type] == elem])

                """
                #Estimating the % of IS genes colocalizing with each element and comparing them with 
                d_simulations_colocalization_with_IS = convert_spot_to_number_of_elements(df_simulations_IS, outdir, "IS")
                n_2 = 2
                for elem2, associated_df in  d_simulations_colocalization_with_IS.items():
                    n_2 += 1
                    associated_df.to_csv(f"{outdir}/TAT_defense_results/simulation/data/3.{n}.{n_2}-Random_distribution_of_IS_coloc_with_{elem2}_number_per_spot.tsv", sep = "\t")
                    plot_coloc_mean_distribution(f"{outdir}/TAT_defense_results/simulation/plot/3.{n}.{n_2}-Random_spot_mean_distribution_IS_coloc_with_{elem2}_plot.png", associated_df, df_all_elements[df_all_elements[family_or_type] == elem], elem2)   
                """

    
        elif simulation_to_run == "integrase":
            
            list_element_type_to_simulate = []
            if to_consider == "all":
                list_element_type_to_simulate = ["integrase"]
                family_or_type = "Element_type"
            elif to_consider == "per_family":
                list_element_type_to_simulate = df_all_elements[df_all_elements["Element_type"] == "integrase"]["Element_family"].drop_duplicates().tolist()
                family_or_type = "Element_family"
                
            print("Initializing the random distribution of the integrase transposase genes")
            n=0
            for elem in tqdm(list_element_type_to_simulate, desc = "Running the simulations for the integrase genes (either all, or based on their family)", total = len(list_element_type_to_simulate)):
                n += 1
                df_simulations_integrase = run_simulation(df_spot_info, df_all_elements[df_all_elements[family_or_type] == elem], n_simulations = n_simulations)
                d_df_res_simulations["integrase"] = df_simulations_integrase
                df_simulations_integrase.to_csv(f"{outdir}/TAT_defense_results/simulation/data/4.{n}.1-Random_spot_distribution_{elem}.tsv", sep = "\t")
                df_simulations_integrase_gene_numbers = convert_spot_to_number_of_genes(df_simulations_integrase, df_spot_info)
                df_simulations_integrase_gene_numbers.to_csv(f"{outdir}/TAT_defense_results/simulation/data/4.{n}.2-Random_spot_distribution_{elem}_gene_numbers_within_spot.tsv", sep = "\t")
                plot_CDF_distribution(f"{outdir}/TAT_defense_results/simulation/plot/4.{n}.1-Random_spot_distribution_{elem}_plot.png", df_simulations_integrase_gene_numbers, df_all_elements[df_all_elements[family_or_type] == elem])
                plot_mean_distribution(f"{outdir}/TAT_defense_results/simulation/plot/4.{n}.2-Random_spot_mean_distribution_{elem}_plot.png", df_simulations_integrase_gene_numbers, df_all_elements[df_all_elements[family_or_type] == elem])

                """
                #Estimating the % of integrase genes colocalizing with each element and comparing them with 
                d_simulations_colocalization_with_integrase = convert_spot_to_number_of_elements(df_simulations_integrase, outdir, "integrases")
                n_2 = 2
                for elem2, associated_df in  d_simulations_colocalization_with_integrase.items():
                    n_2 += 1
                    associated_df.to_csv(f"{outdir}/TAT_defense_results/simulation/data/4.{n}.{n_2}-Random_distribution_of_integrase_coloc_with_{elem2}_number_per_spot.tsv", sep = "\t")
                    plot_coloc_mean_distribution(f"{outdir}/TAT_defense_results/simulation/plot/4.{n}.{n_2}-Random_spot_mean_distribution_integrase_coloc_with_{elem2}_plot.png", associated_df, df_all_elements[df_all_elements[family_or_type] == elem], elem2)   
                """


    #Now estimating the % of systems colocaling with others, and comparing with the % of systems colocalizing with the randomly distributed system
    for element_type in to_test:
        d_simulated_spot_from_element_perspective = convert_spot_to_number_of_elements_simulated(d_df_res_simulations, df_all_elements, element_type)
        n_1 = initialize_n(element_type)
        n_2 = 2
        for k, v in d_simulated_spot_from_element_perspective.items():
            n_2 += 1
            v.to_csv(f"{outdir}/TAT_defense_results/simulation/data/{n_1}.1.{n_2}_Random_distribution_of_{k}_from_{element_type}_perspective.tsv", sep = "\t")
            plot_mean_distribution_from_perspective(v,
                                                    df_all_elements[df_all_elements[family_or_type] == element_type],
                                                    k,
                                                    element_type,
                                                    f"{outdir}/TAT_defense_results/simulation/data/{n_1}.1.{n_2}_Random_distribution_of_{k}_from_{element_type}_perspective.png")






def initialize_n(elem):
    """
    """
    
    if elem == "TA":
        return 1
    elif elem == "defense":
        return 2
    elif elem == "IS":
        return 3
    elif elem == "integrase":
        return 4


def fetch_elements_features(outdir, list_to_test, df_spot_info):
    """
    Small function which aims to fetch all necessary informations from the TAT_defense_search to perform the simulations.
    Will return a pandas dataframe (and save it in the simulation folder) with five columns:
    First : Element_name, Second : Element_type (either IS/Integrase/TA/defense), Third : Element_length (in number of genes), Fourth: Spot_number, Fith: Element_genes
    """
    
    df_spot_info = df_spot_info.set_index("interval_number")    
    d_res = {}
    
    for element_to_fetch in list_to_test:
        
        d_n_elements = {} # dict to store and increase the number of each given element
        
        if element_to_fetch == "TA":
            print("Fetching the necessary genomic informations of each TA for simulations")
            
            df_TA = pd.read_csv(f"{outdir}/TAT_defense_results/TA_defense_search/1.1-TA_spot_pangenome.tsv", sep = "\t")
            for row in df_TA[df_TA["Representative"] == True].iterrows():
                if row[1]["TA_family"] not in d_n_elements.keys():
                    d_n_elements[row[1]["TA_family"]] = 1
                else :
                    d_n_elements[row[1]["TA_family"]] += 1

                d_res[str(row[1]["TA_family"]) + "_" + str(d_n_elements[row[1]["TA_family"]])] = {"Element_name" : str(row[1]["TA_family"]) + "_" + str(d_n_elements[row[1]["TA_family"]]),
                                                                                                  "Element_type": "TA",
                                                                                                  "Element_family": str(row[1]["TA_family"]).replace("/", "_"),
                                                                                                  "Element_length": len(row[1]["Ref_TA"].split("-")),
                                                                                                  "Spot_number": str(row[1]["SPOT"]),
                                                                                                  "Element_genes": row[1]["Ref_TA"],
                                                                                                  "N_NR_genes_in_spot": df_spot_info.loc[int(row[1]["SPOT"]), "N_NR_genes"]
                                                                                                  }

        elif element_to_fetch == "defense":
            print("Fetching the necessary genomic informations of each defense system for simulations")
            
            df_defense = pd.read_csv(f"{outdir}/TAT_defense_results/TA_defense_search/2.2-Spot_defense_finder_systems.tsv", sep = "\t")
            df_defense_finder_search = pd.read_csv(f"{outdir}/TAT_defense_results/TA_defense_search/all_prt_formated_defense_finder_systems.tsv", sep = "\t")
            df_names = pd.read_csv(f"{outdir}/TAT_defense_results/tmp/gene_names_formated.tsv", sep = "\t", names = ["formated_name","real_name"])
            df_real2formated = df_names.set_index("real_name")
            df_formated2real = df_names.set_index("formated_name")
            
            for defense_type in [x for x in df_defense.columns if x not in ["Spot_number", "MazEF"]]:
                d_n_elements[defense_type] = 0
                for row in df_defense.iterrows():
                    if isinstance(row[1][defense_type], float):
                        continue
                    else : #There is this type of defense in the spot, so we count their number within the spot
                        tmp_list_gene_already_attributed = []
                        
                        for gene in row[1][defense_type].split(","):
                            if gene not in tmp_list_gene_already_attributed:
                                
                                index_gene_search = df_defense_finder_search[df_defense_finder_search["protein_in_syst"].str.contains(str(df_real2formated.loc[gene, "formated_name"]), na=False)].index[0]
                                tmp_defense_system_gene_list = [df_formated2real.loc[x, "real_name"] for x in df_defense_finder_search.loc[index_gene_search, "protein_in_syst"].split(",")]
                                tmp_list_gene_already_attributed += tmp_defense_system_gene_list
                            
                                d_n_elements[defense_type] += 1
                                d_res[str(defense_type) + "_" + str(d_n_elements[defense_type])] = {"Element_name" : str(defense_type) + "_" + str(d_n_elements[defense_type]),
                                                                                                "Element_type": "defense",
                                                                                                "Element_family": str(defense_type).replace("/", "_"),
                                                                                                "Element_length": len(tmp_defense_system_gene_list),
                                                                                                "Spot_number": row[1]["Spot_number"],
                                                                                                "Element_genes": ("-").join(tmp_defense_system_gene_list),
                                                                                                "N_NR_genes_in_spot": df_spot_info.loc[int(row[1]["Spot_number"]), "N_NR_genes"]
                                                                                                }
        
        elif element_to_fetch == "IS": 
            print("Fetching the necessary genomic informations of each IS for simulations")
            
            df_IS = pd.read_csv(f"{outdir}/TAT_defense_results/TA_defense_search/3.1-IS_spot_pangenome.tsv", sep = "\t")
            for row in df_IS.iterrows():
                if row[1]["mge_hit"] not in d_n_elements.keys():
                    d_n_elements[row[1]["mge_hit"]] = 1
                else :
                    d_n_elements[row[1]["mge_hit"]] += 1

                d_res[str(row[1]["mge_hit"]) + "_" + str(d_n_elements[row[1]["mge_hit"]])] = {"Element_name" : str(row[1]["mge_hit"]) + "_" + str(d_n_elements[row[1]["mge_hit"]]),
                                                                                                  "Element_type": "IS",
                                                                                                  "Element_family": str(row[1]["mge_hit"]).replace("/", "_"),
                                                                                                  "Element_length": 1,
                                                                                                  "Spot_number": str(row[1]["spot_number"]),
                                                                                                  "Element_genes": row[1]["gene_name"],
                                                                                                  "N_NR_genes_in_spot": df_spot_info.loc[int(row[1]["spot_number"]), "N_NR_genes"]
                                                                                                  }
        elif element_to_fetch == "integrase":
            print("Fetching the necessary genomic informations of each integrase for simulations")
            
            df_integrase = pd.read_csv(f"{outdir}/TAT_defense_results/TA_defense_search/4.1-Integrase_spot_pangenome.tsv", sep = "\t")
            for row in df_integrase.iterrows():
                if row[1]["mge_hit"] not in d_n_elements.keys():
                    d_n_elements[row[1]["mge_hit"]] = 1
                else :
                    d_n_elements[row[1]["mge_hit"]] += 1

                d_res[str(row[1]["mge_hit"]) + "_" + str(d_n_elements[row[1]["mge_hit"]])] = {"Element_name" : str(row[1]["mge_hit"]) + "_" + str(d_n_elements[row[1]["mge_hit"]]),
                                                                                                  "Element_type": "integrase",
                                                                                                  "Element_family": str(row[1]["mge_hit"]).replace("/", "_"),
                                                                                                  "Element_length": 1,
                                                                                                  "Spot_number": str(row[1]["spot_number"]),
                                                                                                  "Element_genes": row[1]["gene_name"],
                                                                                                  "N_NR_genes_in_spot": df_spot_info.loc[int(row[1]["spot_number"]), "N_NR_genes"]
                                                                                                  }
    
    df_res = pd.DataFrame.from_dict(d_res, orient = "index")
    #df_res.to_csv(f"{outdir}/TAT_defense_results/simulation/0-Elements_genomic_features.tsv", sep = "\t", index = False)
    
    return df_res



def fetch_coloc_features(df_in, outdir, to_test):
    """
    Function which add a column representing the number of element retrieved within the spot of the studied element.
    One column per element tested.
    Note that, when establishing the number of elements of the same type than the one studied (e.g. looking at the number of TA colocalizing with a TA), we substract one.
    So we do not count the element itself.
    """
    
    df_search_results = pd.read_csv(f"{outdir}/TAT_defense_results/TA_defense_search/5.0-Spot_merged_results_only_total.tsv", sep = "\t")
    df_search_results.rename(columns={"n_TAs":"n_TA", "n_integrases": "n_integrase"}, inplace= True) # Messy but quick fix to avoid problem later
    df_search_results = df_search_results.set_index("Spot_number")
    df_search_results = df_search_results.fillna(0)
    df_out = df_in.copy(deep = True)
    
    for elem in to_test:
        list_column2add = []
        for row in df_in.iterrows():
            if row[1]["Element_type"] == elem:
                list_column2add.append(int(df_search_results.loc[int(row[1]["Spot_number"]), f"n_{elem}"]) - 1)
            else :
                list_column2add.append(int(df_search_results.loc[int(row[1]["Spot_number"]), f"n_{elem}"]))
        
        df_out[f"N_{elem}_colocalizing"] = list_column2add
    
    df_out.to_csv(f"{outdir}/TAT_defense_results/simulation/0-Elements_genomic_features.tsv", sep = "\t", index = False)
    return df_out



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


def convert_spot_to_number_of_genes(df_in, df_spot_info):
    """
    """
    df_spot_info = df_spot_info.set_index("interval_number")
    list_simulations = df_in.columns.tolist()
    d_out = {}
    
    for row in df_in.iterrows():
        d_out[row[0]] = {}
        for i in list_simulations:
            d_out[row[0]][i] = df_spot_info.loc[row[1][i], "N_NR_genes"]
    
    df_out = pd.DataFrame.from_dict(d_out, orient = "index")
    
    return df_out

def convert_spot_to_number_of_elements(df_in, outdir, studied_element):
    """
    return a Dict containing all df_in converted with the number of elements retrieved in the spot.
    Dict: k = TA/defense/IS/integrase, v = dataframe converted for the corresponding element
    """
    
    df_search_results = pd.read_csv(f"{outdir}/TAT_defense_results/TA_defense_search/5.0-Spot_merged_results_only_total.tsv", sep = "\t")
    df_search_results["Spot_number"] = df_search_results["Spot_number"].astype('int64')
    df_search_results = df_search_results.set_index("Spot_number")
    df_search_results = df_search_results.fillna(0)
    list_simulations = df_in.columns.tolist()
    d_return = {}
    list_element_to_process = [x for x in df_search_results.columns.tolist() if x != "n_" + studied_element]
    
    # First we determine the number of the studied elements that colocalize with other(s)
    if studied_element == "TAs":
        studied_element = "TA" #Messy & quick fix to avoid problem later
    if studied_element == "integrases":
        studied_element = "integrase" #Messy & quick fix to avoid problem later
    
    d_out = {}
    for row in df_in.iterrows():
        d_out[row[0]] = {}
        for i in list_simulations:
            #we count the number of occurences that this elements was randomly distributed in this spot, and substract one to not count the element itself
            d_out[row[0]][i] = int(df_in[i].value_counts().get(row[1][i], 0)) - 1
    
    d_return[studied_element] = pd.DataFrame.from_dict(d_out, orient = "index")
    
    
    for element in list_element_to_process:
        d_out = {}
        
        if element == "n_TAs":
            for row in df_in.iterrows():
                d_out[row[0]] = {}
                for i in list_simulations:
                    if row[1][i] in df_search_results.index.tolist():
                        d_out[row[0]][i] = df_search_results.loc[row[1][i], "n_TAs"]
                    else :
                        d_out[row[0]][i] = 0

            d_return["TA"] =  pd.DataFrame.from_dict(d_out, orient = "index")
            
        elif element == "n_defense":
            for row in df_in.iterrows():
                d_out[row[0]] = {}
                for i in list_simulations:
                    if row[1][i] in df_search_results.index.tolist():
                        d_out[row[0]][i] = df_search_results.loc[row[1][i], "n_defense"]
                    else :
                        d_out[row[0]][i] = 0
                        
            d_return["defense"] =  pd.DataFrame.from_dict(d_out, orient = "index")
        
        elif element == "n_IS":
            for row in df_in.iterrows():
                d_out[row[0]] = {}
                for i in list_simulations:
                    if row[1][i] in df_search_results.index.tolist():
                        d_out[row[0]][i] = df_search_results.loc[row[1][i], "n_IS"]
                    else :
                        d_out[row[0]][i] = 0
                        
            d_return["IS"] =  pd.DataFrame.from_dict(d_out, orient = "index")
            
        elif element == "n_integrases":
            for row in df_in.iterrows():
                d_out[row[0]] = {}
                for i in list_simulations:
                    if row[1][i] in df_search_results.index.tolist():
                        d_out[row[0]][i] = df_search_results.loc[row[1][i], "n_integrases"]
                    else :
                        d_out[row[0]][i] = 0

            d_return["integrase"] =  pd.DataFrame.from_dict(d_out, orient = "index")
        
    
    return d_return



def convert_spot_to_number_of_elements_simulated(d_simulations, df_all_elements, elem_of_interest):
    """
    
    """
    
    d_out = {}
    
    for k, df_sim in d_simulations.items():
        d_res = {}
        list_simulation = df_sim.columns.tolist()
        if k == elem_of_interest: #we count the number of occurences that this elements was randomly distributed in this spot, and substract one to not count the element itself       
            list_index_sim = df_sim.index.tolist()
            for index in tqdm(list_index_sim, total = len(list_index_sim), desc = f"Fetching the number of {k} elements being randomly distributed in the spots containing a real {elem_of_interest}"):
                d_res[index] = {}            
                for sim in list_simulation:
                    d_res[index][sim] = int(df_sim[sim].value_counts().get(df_sim.loc[index, sim], 0)) - 1
        
        else :
            list_index_real = df_all_elements[df_all_elements["Element_type"] == elem_of_interest].index.tolist()
            for index in tqdm(list_index_real, desc = f"Fetching the number of {k} elements being randomly distributed in the spots containing {elem_of_interest}", total = len(list_index_real)):
                d_res[index] = {}
                for sim in list_simulation:
                    d_res[index][sim] = int(df_sim[sim].value_counts().get(df_all_elements.loc[index, "Spot_number"], 0))
        
        d_out[k] = pd.DataFrame.from_dict(d_res, orient = "index")
    
    
    return d_out
        






def plot_CDF_distribution(outname, df_sim, real_data):
    """
    Function which will plot the distribution of random simulations
    """
    
    real_data = real_data["N_NR_genes_in_spot"]
    
    all_region_sizes = sorted(set(df_sim.values.flatten()).union(set(real_data.values)))

    cdf_list = []
    for col in df_sim.columns:
        counts = df_sim[col].value_counts().sort_index()
        cum_values = counts.cumsum()
        cum_percent = cum_values / cum_values.max() * 100
        cum_percent_full = cum_percent.reindex(all_region_sizes, method='ffill').fillna(0)
        cdf_list.append(cum_percent_full)

    cdf_df = pd.DataFrame(cdf_list)
    cdf_mean = cdf_df.mean()
    cdf_std = cdf_df.std()

    real_counts = real_data.value_counts().sort_index()
    real_cum = real_counts.cumsum()
    real_cdf = real_cum / real_cum.max() * 100
    real_cdf_full = real_cdf.reindex(all_region_sizes, method='ffill').fillna(0)

    plt.figure(figsize=(10, 6))

    plt.plot(all_region_sizes, cdf_mean, label="CDF moyenne (simulations)", color="blue")
    plt.fill_between(all_region_sizes, cdf_mean - cdf_std, cdf_mean + cdf_std,
                    color="blue", alpha=0.3, label="±1 écart-type")


    plt.plot(all_region_sizes, real_cdf_full, label="Données réelles", color="red", linewidth=2)

    plt.xlabel("Taille de la région (nb de gènes accessoires)")
    plt.ylabel("Densité cumulative (%)")
    plt.title("CDF cumulative : simulations vs données réelles")
    plt.ylim(0, 100)
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(outname)
    plt.close()

def plot_mean_distribution(outname, df_sim, real_data):
    """
    Function which will plot the distribution of the mean of the random simulations
    """

    real_data = real_data["N_NR_genes_in_spot"]
    
    simulation_means = df_sim.mean(axis=0)

    real_data = real_data.squeeze()
    real_mean = real_data.mean()

    n_above = (simulation_means > real_mean).sum()
    n_below = (simulation_means < real_mean).sum()
    n_equal = (simulation_means == real_mean).sum()
    p_value_lower = (simulation_means <= real_mean).sum() / len(simulation_means)
    p_value_higher = (simulation_means >= real_mean).sum() / len(simulation_means)
    
    plt.figure(figsize=(10, 6))
    sns.histplot(simulation_means, kde=True, bins=30, color="blue", alpha=0.6, label="Moyennes (simulations)")
    plt.axvline(real_mean, color="red", linestyle="--", linewidth=2, label=f"Moyenne réelle = {real_mean:.2f}")

    plt.text(0.95, 0.8,
             f"Simulations > real : {n_above} ({round(n_above/len(simulation_means)*100, 2)}%)\n"
             f"Unilateral p-value > real : {p_value_higher}\n"
             f"Simulations ≤ réelle : {n_below + n_equal} ({round((n_below + n_equal)/len(simulation_means)*100, 2)}%)\n"
             f"Unilateral p-value ≤ real : {p_value_lower}",
             transform=plt.gca().transAxes,
             fontsize=10,
             verticalalignment='top',
             horizontalalignment='right',
             bbox=dict(boxstyle="round", facecolor="white", edgecolor="gray", alpha=0.8))
    
    
    plt.xlabel("Moyenne de taille de région (nb de gènes accessoires)")
    plt.ylabel("Fréquence")
    plt.title("Distribution des moyennes par simulation vs moyenne réelle")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(outname)
    plt.close()


def plot_coloc_mean_distribution(outname, df_sim, real_data, elem_colocalizing):
    """
    """
    
    real_data = real_data[f"N_{elem_colocalizing}_colocalizing"]
    
    simulation_means = df_sim.mean(axis=0)

    real_data = real_data.squeeze()
    real_mean = real_data.mean()

    n_above = (simulation_means > real_mean).sum()
    n_below = (simulation_means < real_mean).sum()
    n_equal = (simulation_means == real_mean).sum()
    p_value_lower = (simulation_means <= real_mean).sum() / len(simulation_means)
    p_value_higher = (simulation_means >= real_mean).sum() / len(simulation_means)

    plt.figure(figsize=(10, 6))

    sns.histplot(simulation_means, kde=True, bins=30, color="blue", alpha=0.6, label="Moyennes (simulations)")

    plt.axvline(real_mean, color="red", linestyle="--", linewidth=2, label=f"Moyenne réelle = {real_mean:.2f}")

    plt.text(0.95, 0.8,
             f"Simulations > real : {n_above} ({round(n_above/len(simulation_means)*100, 2)}%)\n"
             f"Unilateral p-value > real : {p_value_higher}\n"
             f"Simulations ≤ réelle : {n_below + n_equal} ({round((n_below + n_equal)/len(simulation_means)*100, 2)}%)\n"
             f"Unilateral p-value ≤ real : {p_value_lower}",
             transform=plt.gca().transAxes,
             fontsize=10,
             verticalalignment='top',
             horizontalalignment='right',
             bbox=dict(boxstyle="round", facecolor="white", edgecolor="gray", alpha=0.8))
    
    plt.xlabel("Moyenne de systems qui colocalize avec")
    plt.ylabel("Fréquence")
    plt.title(f"Distribution des moyennes par simulation vs moyenne réelle de colocalization avec {elem_colocalizing}")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(outname)
    plt.close()
    

def plot_mean_distribution_from_perspective(df_sim, df_real, elem_colocalizing, elem_perspective_name, outname):
    
    """
    Plot the mean distribution of the real data and simulations from a same element perspectives.
    """
    
    real_data = df_real[f"N_{elem_colocalizing}_colocalizing"]
    
    simulation_means = df_sim.mean(axis=0)

    real_data = real_data.squeeze()
    real_mean = real_data.mean()

    n_above = (simulation_means > real_mean).sum()
    n_below = (simulation_means < real_mean).sum()
    n_equal = (simulation_means == real_mean).sum()
    p_value_lower = (simulation_means <= real_mean).sum() / len(simulation_means)
    p_value_higher = (simulation_means >= real_mean).sum() / len(simulation_means)
    
    plt.figure(figsize=(10, 6))
    sns.histplot(simulation_means, kde=True, bins=30, color="blue", alpha=0.6, label="Moyennes (simulations)")
    plt.axvline(real_mean, color="red", linestyle="--", linewidth=2, label=f"Moyenne réelle = {real_mean:.2f}")

    plt.text(0.95, 0.8,
             f"Simulations > real : {n_above} ({round(n_above/len(simulation_means)*100, 2)}%)\n"
             f"Unilateral p-value > real : {p_value_higher}\n"
             f"Simulations ≤ réelle : {n_below + n_equal} ({round((n_below + n_equal)/len(simulation_means)*100, 2)}%)\n"
             f"Unilateral p-value ≤ real : {p_value_lower}",
             transform=plt.gca().transAxes,
             fontsize=10,
             verticalalignment='top',
             horizontalalignment='right',
             bbox=dict(boxstyle="round", facecolor="white", edgecolor="gray", alpha=0.8))
    
    
    plt.xlabel("Moyenne de systems qui colocalize avec")
    plt.ylabel("Fréquence")
    plt.title(f"{elem_perspective_name} colocalization perspectives with (or not) randomly distributed {elem_colocalizing}")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(outname)
    plt.close()
    
    
    
    
    


if __name__ == "__main__":
    import sys
    main(sys.argv[1], sys.argv[2], sys.argv[3])