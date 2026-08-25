#!/usr/bin/env python

import pandas as pd
import numpy as np
from tqdm import tqdm
import os
import random
import copy
import matplotlib.pyplot as plt
import seaborn as sns

"""
Standalone script which aims to use the data generated from TAT_defense.py core mapping and search modules to perform random simulations of a given element group and plot the spot-colocalization with TAs.

cl : 

python script/Custom_sim.py test/0-39_genomes_Ph25 test/input/Custom_group/Defense_all.lst test/0-39_genomes_Ph25/TAT_defense_results/1.0-Coloc_simulation_TA_all_defense Defense[All] TA
python script/Custom_sim.py test/0-39_genomes_Ph25 test/input/Custom_group/IS_all.lst test/0-39_genomes_Ph25/TAT_defense_results/1.2-Coloc_simulation_TA_IS_all IS[all] TA
python script/Custom_sim.py test/0-39_genomes_Ph25 test/input/Custom_group/All_TAs.lst test/0-39_genomes_Ph25/TAT_defense_results/1.1-Coloc_simulation_TA_TA TA[all] TA
python script/Custom_sim.py test/0-39_genomes_Ph25 test/input/Custom_group/Integrase_all.lst test/0-39_genomes_Ph25/TAT_defense_results/1.4-Coloc_simulation_TA_integrase Integrase[all] TA
python script/Custom_sim.py test/0-39_genomes_Ph25 test/input/Custom_group/Integrase_resolvase.lst test/0-39_genomes_Ph25/TAT_defense_results/1.4-Coloc_simulation_TA_integrase_resolvase Integrase[resolvase] TA
python script/Custom_sim.py test/0-39_genomes_Ph25 test/input/Custom_group/Integrase_phage.lst test/0-39_genomes_Ph25/TAT_defense_results/1.5-Coloc_simulation_TA_integrase_phage Integrase[phage] TA
python script/Custom_sim.py test/0-39_genomes_Ph25 test/input/Custom_group/IS_all.lst test/0-39_genomes_Ph25/TAT_defense_results/1.6-Coloc_simulation_TA_IS IS[all] TA
python script/Custom_sim.py test/0-39_genomes_Ph25 test/input/Custom_group/IS_cointegrate.lst test/0-39_genomes_Ph25/TAT_defense_results/1.7-Coloc_simulation_TA_IS_cointegrate IS[cointegrate] TA
python script/Custom_sim.py test/0-39_genomes_Ph25 test/input/Custom_group/IS_copy.lst test/0-39_genomes_Ph25/TAT_defense_results/1.8-Coloc_simulation_TA_IS_copy IS[copy] TA
python script/Custom_sim.py test/0-39_genomes_Ph25 test/input/Custom_group/IS_cut_and_peel.lst test/0-39_genomes_Ph25/TAT_defense_results/1.9-Coloc_simulation_TA_IS_cut_peel IS[Cut-Peel] TA
python script/Custom_sim.py test/0-39_genomes_Ph25 test/input/Custom_group/IS_unknown_mechanism.lst test/0-39_genomes_Ph25/TAT_defense_results/1.10-Coloc_simulation_TA_IS_unknown IS[Unkown] TA


python script/Custom_sim.py test/0-39_genomes_Ph25 test/input/Custom_group/Defense_all.lst test/0-39_genomes_Ph25/TAT_defense_results/2.0-Coloc_simulation_defense_defense defense[all] defense
python script/Custom_sim.py test/0-39_genomes_Ph25 test/input/Custom_group/All_TAs.lst test/0-39_genomes_Ph25/TAT_defense_results/2.1-Coloc_simulation_defense_TA TA[all] defense
python script/Custom_sim.py test/0-39_genomes_Ph25 test/input/Custom_group/IS_all.lst test/0-39_genomes_Ph25/TAT_defense_results/2.2-Coloc_simulation_defense_IS IS[all] defense
python script/Custom_sim.py test/0-39_genomes_Ph25 test/input/Custom_group/Integrase_all.lst test/0-39_genomes_Ph25/TAT_defense_results/2.4-Coloc_simulation_defense_integrase Integrase[all] defense
python script/Custom_sim.py test/0-39_genomes_Ph25 test/input/Custom_group/IS_cointegrate.lst test/0-39_genomes_Ph25/TAT_defense_results/2.5-Coloc_simulation_defense_IS_cointegrate IS[cointegrate] defense
python script/Custom_sim.py test/0-39_genomes_Ph25 test/input/Custom_group/IS_copy.lst test/0-39_genomes_Ph25/TAT_defense_results/2.6-Coloc_simulation_defense_IS_copy IS[copy] defense
python script/Custom_sim.py test/0-39_genomes_Ph25 test/input/Custom_group/IS_cut_and_peel.lst test/0-39_genomes_Ph25/TAT_defense_results/2.7-Coloc_simulation_defense_IS_cut_peel IS[Cut-Peel] defense
python script/Custom_sim.py test/0-39_genomes_Ph25 test/input/Custom_group/IS_unknown_mechanism.lst test/0-39_genomes_Ph25/TAT_defense_results/2.8-Coloc_simulation_defense_IS_unknown IS[unknown] defense
python script/Custom_sim.py test/0-39_genomes_Ph25 test/input/Custom_group/Integrase_phage.lst test/0-39_genomes_Ph25/TAT_defense_results/2.9-Coloc_simulation_defense_integrase_phage Integrase[integrase-recombinase] defense
python script/Custom_sim.py test/0-39_genomes_Ph25 test/input/Custom_group/Integrase_resolvase.lst test/0-39_genomes_Ph25/TAT_defense_results/2.10-Coloc_simulation_defense_integrase_resolvase Integrase[resolvase] defense

python script/Custom_sim.py test/0-39_genomes_Ph25 test/input/Custom_group/IS_cut_and_peel.lst test/0-39_genomes_Ph25/TAT_defense_results/13.3-Coloc_simulation_defense_IS_cut_peel IS[cut/peel] defense
python script/Custom_sim.py test/0-39_genomes_Ph25 test/input/Custom_group/IS_cointegrate_copy.lst test/0-39_genomes_Ph25/TAT_defense_results/13.4-Coloc_simulation_defense_IS_cointegrate_copy IS[cointegrate/copy] defense
python script/Custom_sim.py test/0-39_genomes_Ph25 test/input/Custom_group/IS_cointegrate.lst test/0-39_genomes_Ph25/TAT_defense_results/13.4-Coloc_simulation_defense_IS_cointegrate IS[cointegrate] defense
python script/Custom_sim.py test/0-39_genomes_Ph25 test/input/Custom_group/IS_unknown_mechanism.lst test/0-39_genomes_Ph25/TAT_defense_results/13.6-Coloc_simulation_defense_IS_unknown IS[unknown] defense
python script/Custom_sim.py test/0-39_genomes_Ph25 test/input/Custom_group/IS_copy.lst test/0-39_genomes_Ph25/TAT_defense_results/13.7-Coloc_simulation_defense_IS_copy IS[copy] defense
python script/Custom_sim.py test/0-39_genomes_Ph25 test/input/Custom_group/Integrase_resolvase.lst test/0-39_genomes_Ph25/TAT_defense_results/13.8-Coloc_simulation_defense_integrase_resolvase Integrase[resolvase] defense
python script/Custom_sim.py test/0-39_genomes_Ph25 test/input/Custom_group/Integrase_phage.lst test/0-39_genomes_Ph25/TAT_defense_results/13.9-Coloc_simulation_defense_integrase_phage Integrase[phage] defense
python script/Custom_sim.py test/0-39_genomes_Ph25 test/input/Custom_group/Integrase_all.lst test/0-39_genomes_Ph25/TAT_defense_results/2.4-Coloc_simulation_defense_integrase_integrase Integrase[all] defense

python script/Custom_sim.py test/0-39_genomes_Ph25 test/input/Custom_group/All_TAs.lst test/0-39_genomes_Ph25/TAT_defense_results/3.1-Coloc_simulation_IS_all_TA TA[all] IS
python script/Custom_sim.py test/0-39_genomes_Ph25 test/input/Custom_group/Defense_all.lst test/0-39_genomes_Ph25/TAT_defense_results/3.2-Coloc_simulation_IS_all_defense Defense[all] IS
python script/Custom_sim.py test/0-39_genomes_Ph25 test/input/Custom_group/IS_all.lst test/0-39_genomes_Ph25/TAT_defense_results/3.3-Coloc_simulation_IS_IS IS[all] IS
python script/Custom_sim.py test/0-39_genomes_Ph25 test/input/Custom_group/Integrase_all.lst test/0-39_genomes_Ph25/TAT_defense_results/3.4-Coloc_simulation_IS_integrase Integrase[all] IS

python script/Custom_sim.py test/0-39_genomes_Ph25 test/input/Custom_group/Defense_all.lst test/0-39_genomes_Ph25/TAT_defense_results/4.2-Coloc_simulation_integrase_all_defense Defense[all] integrase
python script/Custom_sim.py test/0-39_genomes_Ph25 test/input/Custom_group/All_TAs.lst test/0-39_genomes_Ph25/TAT_defense_results/4.3-Coloc_simulation_integrase_all_TA TA[all] integrase
python script/Custom_sim.py test/0-39_genomes_Ph25 test/input/Custom_group/Integrase_all.lst test/0-39_genomes_Ph25/TAT_defense_results/4.0-Coloc_simulation_integrase_integrase Integrase[all] integrase
python script/Custom_sim.py test/0-39_genomes_Ph25 test/input/Custom_group/IS_all.lst test/0-39_genomes_Ph25/TAT_defense_results/4.1-Coloc_simulation_integrase_IS IS[all] integrase


"""


def main(TAT_directory, file_custom_group, outdir, group_name, elem_of_interest, to_test = ["TA", "defense", "IS", "integrase"], plot_for_article = True, plot_legend_article_figure = True):
    """
    TAT_directory : str, path to the directory containing the results from the TAT_defense.py script (core mapping and search modules)
    file_custom_group : path to a text file containing the list of elements to put in a same group to run the simulation
    outdir : str, path to the directory where the simulation will be saved
    """
    
    #first, creating the outdir
    create_dir(outdir)
    
    #get the list of elements to group for the simulation
    list_elem2 = fetch_elem_list(file_custom_group)
    
    # load the tsv file with the spot informations (number of non redundant genes per spot)
    df_spot_info = pd.read_csv(f"{TAT_directory}/TAT_defense_results/core_mapping/3-all_spots_features.tsv", sep = "\t")
    
    if os.path.exists(f"{outdir}/0-Elements_genomic_features.tsv"):
        df_all_elements = pd.read_csv(f"{outdir}/0-Elements_genomic_features.tsv", sep = "\t")
    
    else :
        df_all_elements = fetch_elements_features(TAT_directory, to_test, df_spot_info)
        df_all_elements = fetch_coloc_features(df_all_elements, outdir, to_test, TAT_directory)
        df_all_elements = compute_sum_custom_group_coloc(df_all_elements, list_elem2)
        df_all_elements.to_csv(f"{outdir}/0-Elements_genomic_features.tsv", sep = "\t", index = False)
        
    #run the simulation
    if os.path.exists(f"{outdir}/1-Custom_group_random_distribution.tsv") == False :
        df_simulations = run_simulation(df_spot_info, df_all_elements[df_all_elements["Element_family"].isin(list_elem2)], n_simulations = 10000)
        df_simulations.to_csv(f"{outdir}/1-Custom_group_random_distribution.tsv", sep = "\t")
    else :
        df_simulations = pd.read_csv(f"{outdir}/1-Custom_group_random_distribution.tsv", sep = "\t", index_col=0)
    
    #Convert the df_simulations spot_numbers to the numbers NR of genes within the associated spot
    if os.path.exists(f"{outdir}/1.2-Custom_group_random_distribution_with_spot_numbers.tsv") == False:
        df_simulations_with_spot_numbers = convert_spot_to_number_NRg(df_simulations, TAT_directory)
        df_simulations_with_spot_numbers.to_csv(f"{outdir}/1.2-Custom_group_random_distribution_with_spot_numbers.tsv", sep = "\t")
    else :
        df_simulations_with_spot_numbers = pd.read_csv(f"{outdir}/1.2-Custom_group_random_distribution_with_spot_numbers.tsv", sep = "\t", index_col=0)
    
    df_all_elements = pd.read_csv(f"{outdir}/0-Elements_genomic_features.tsv", sep = "\t") #Messy fix to avoid problems with the dataframe cell type later
    
    if os.path.exists(f"{outdir}/2-Colocalization_perspectives_from_{elem_of_interest}_point_of_view_with_custom_group.tsv") == False:
        df_elem_family_with_TA = convert_spot_to_number_of_elements_simulated(df_simulations, df_all_elements[df_all_elements["Element_type"] == elem_of_interest])
        df_elem_family_with_TA.to_csv(f"{outdir}/2-Colocalization_perspectives_from_{elem_of_interest}_point_of_view_with_custom_group.tsv", sep = "\t")
    
    else :
        df_elem_family_with_TA = pd.read_csv(f"{outdir}/2-Colocalization_perspectives_from_{elem_of_interest}_point_of_view_with_custom_group.tsv", index_col= 0, sep = "\t")
    
    
    if plot_for_article == False :
        plot_spot_length_distribution_random_group(df_simulations_with_spot_numbers, 
                                       df_all_elements[df_all_elements["Element_family"].isin(list_elem2)],
                                       outdir,
                                       group_name)

        plot_CDF_distribution(f"{outdir}/2.3-plot_length_{group_name}_distribution_simulation_for_article.png", 
                              df_simulations_with_spot_numbers, 
                              df_all_elements[df_all_elements["Element_family"].isin(list_elem2)], 
                              group_name)

        global_res = plot_mean_distribution_from_perspective(df_elem_family_with_TA,
                                                            df_all_elements[df_all_elements["Element_type"] == elem_of_interest],
                                                            elem_of_interest,
                                                            f"{outdir}/3-Colocalization_perspectives_from_{elem_of_interest}_point_of_view_with_custom_group.png",
                                                            group_name
                                                            )
    else :
        plot_spot_length_distribution_random_group_for_article(df_simulations_with_spot_numbers, 
                                                                df_all_elements[df_all_elements["Element_family"].isin(list_elem2)],
                                                                outdir,
                                                                group_name,
                                                                plot_legend_article_figure = plot_legend_article_figure)
        
        plot_CDF_distribution(f"{outdir}/2.3-plot_length_{group_name}_distribution_simulation_for_article.png", 
                              df_simulations_with_spot_numbers, 
                              df_all_elements[df_all_elements["Element_family"].isin(list_elem2)], 
                              group_name)
        
        
        global_res = plot_mean_distribution_from_perspective_for_article(df_elem_family_with_TA,
                                                                        df_all_elements[df_all_elements["Element_type"] == elem_of_interest],
                                                                        elem_of_interest,
                                                                        f"{outdir}/3-Colocalization_perspectives_from_{elem_of_interest}_point_of_view_with_custom_group.png",
                                                                        group_name,
                                                                        plot_legend_article_figure = plot_legend_article_figure
                                                                        )
    
    global_res.insert(1, str(len(df_all_elements[df_all_elements["Element_family"].isin(list_elem2)])))
    write_global_res(global_res, outdir, elem_of_interest, )





def fetch_elements_features(TA_defense_folder, list_to_test, df_spot_info):
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
            
            df_TA = pd.read_csv(f"{TA_defense_folder}/TAT_defense_results/TA_defense_search/1.1-TA_spot_pangenome.tsv", sep = "\t")
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
            
            df_defense = pd.read_csv(f"{TA_defense_folder}/TAT_defense_results/TA_defense_search/2.2-Spot_defense_finder_systems.tsv", sep = "\t")
            df_defense_finder_search = pd.read_csv(f"{TA_defense_folder}/TAT_defense_results/TA_defense_search/all_prt_formated_defense_finder_systems.tsv", sep = "\t")
            df_names = pd.read_csv(f"{TA_defense_folder}/TAT_defense_results/tmp/gene_names_formated.tsv", sep = "\t", names = ["formated_name","real_name"])
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
            
            df_IS = pd.read_csv(f"{TA_defense_folder}/TAT_defense_results/TA_defense_search/3.1-IS_spot_pangenome.tsv", sep = "\t")
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
            
            df_integrase = pd.read_csv(f"{TA_defense_folder}/TAT_defense_results/TA_defense_search/4.1-Integrase_spot_pangenome.tsv", sep = "\t")
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



def fetch_coloc_features(df_in, outdir, to_test, TA_defense_folder):
    """
    Function which add a column representing the number of element retrieved within the spot of the studied element.
    One column per element tested.
    Note that, when establishing the number of elements of the same type than the one studied (e.g. looking at the number of TA colocalizing with a TA), we substract one.
    So we do not count the element itself.
    """
    
    df_search_results = pd.read_csv(f"{TA_defense_folder}/TAT_defense_results/TA_defense_search/5.1-Spot_merged_results_detailled.tsv", sep = "\t")   
    df_search_results.rename(columns={"n_TAs":"n_TA", "n_integrases": "n_integrase"}, inplace= True) # Messy but quick fix to avoid problem later
    df_search_results = df_search_results.set_index("Spot_number")
    df_search_results = df_search_results.fillna(0)
    df_IS_search = pd.read_csv(f"{TA_defense_folder}/TAT_defense_results/TA_defense_search/3.1-IS_spot_pangenome.tsv", sep = "\t")
    df_integrase_search = pd.read_csv(f"{TA_defense_folder}/TAT_defense_results/TA_defense_search/4.1-Integrase_spot_pangenome.tsv", sep = "\t")           
    
    df_out = df_in.copy(deep = True)
    
    for elem in to_test:
        list_column2add = []
        for row in df_in.iterrows():
            if row[1]["Element_type"] == elem:
                list_column2add.append(int(df_search_results.loc[int(row[1]["Spot_number"]), f"n_{elem}"]) - 1)
            else :
                list_column2add.append(int(df_search_results.loc[int(row[1]["Spot_number"]), f"n_{elem}"]))
        
        df_out[f"N_{elem}_colocalizing"] = list_column2add
    
    #Adding each TA family
    for TA_family in ["HipA", "ParE_RelE-like", "VapC", "ataT", "FIC", "Other", "CcdB_MazF", "RES", "HicA", "YafO"]:
        list_column2add = []
        
        if TA_family == "ParE_RelE-like":
            TA_family_name_formated = "ParE/RelE-like"
        elif TA_family == "CcdB_MazF":
            TA_family_name_formated = "CcdB/MazF"
        else:
            TA_family_name_formated = TA_family

        for row in df_in.iterrows():
            if row[1]["Element_family"] == TA_family:
                list_column2add.append(int(df_search_results.loc[int(row[1]["Spot_number"]), f"{TA_family_name_formated}"]) - 1)
            else :
                list_column2add.append(int(df_search_results.loc[int(row[1]["Spot_number"]), f"{TA_family_name_formated}"]))
        
        df_out[f"N_{TA_family}_colocalizing"] = list_column2add
        df_out = df_out.copy()

    #Adding the different defense families
    for defense_family in df_in[df_in["Element_type"] == "defense"]["Element_family"].unique():
        list_column2add = []
        for row in df_in.iterrows():
            if row[1]["Element_family"] == defense_family:
                list_column2add.append(int(df_search_results.loc[int(row[1]["Spot_number"]), defense_family]) - 1)
            else :
                list_column2add.append(int(df_search_results.loc[int(row[1]["Spot_number"]), defense_family]))
        
        df_out[f"N_{defense_family}_colocalizing"] = list_column2add
        df_out = df_out.copy()
    
    #Adding each IS family
    for IS_family in df_IS_search["mge_hit"].unique():
        list_column2add = []
        for row in df_in.iterrows():
            if row[1]["Element_family"] == IS_family:
                list_column2add.append(int(len(df_IS_search[(df_IS_search["mge_hit"] == IS_family) & (df_IS_search["spot_number"] == int(row[1]["Spot_number"]))]) - 1))   
            else :
                list_column2add.append(int(len(df_IS_search[(df_IS_search["mge_hit"] == IS_family) & (df_IS_search["spot_number"] == int(row[1]["Spot_number"]))])))
    
        df_out[f"N_{IS_family}_colocalizing"] = list_column2add

    #Adding each integrase family
    for integrase_family in df_integrase_search["mge_hit"].unique():
        list_column2add = []
        for row in df_in.iterrows():
            if row[1]["Element_family"] == integrase_family:
                list_column2add.append(int(len(df_integrase_search[(df_integrase_search["mge_hit"] == integrase_family) & (df_integrase_search["spot_number"] == int(row[1]["Spot_number"]))]) - 1))   
            else :
                list_column2add.append(int(len(df_integrase_search[(df_integrase_search["mge_hit"] == integrase_family) & (df_integrase_search["spot_number"] == int(row[1]["Spot_number"]))])))
    
        df_out[f"N_{integrase_family}_colocalizing"] = list_column2add

    return df_out


def compute_sum_custom_group_coloc(df_all_elements, list_elem2):
    """
    """
    
    list_custom_res = []
    for _, row in df_all_elements.iterrows():
        n_custom = 0
        for elem in list_elem2:
            n_custom += int(row[f"N_{elem}_colocalizing"])
        list_custom_res.append(n_custom)
    
    df_all_elements["N_custom_group_colocalization"] = list_custom_res
    
    return df_all_elements




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
    d_spot_info = {} #k: listpot_number, v: number of non redundant genes in the spot
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


def convert_spot_to_number_NRg(df_simulations, TAT_defense_dir):
    """
    """
    
    df_spot_tmp = pd.read_csv(f"{TAT_defense_dir}/TAT_defense_results/core_mapping/3-all_spots_features.tsv", sep = "\t", index_col= 0)
    list_sim = df_simulations.columns.tolist()
    d_spot_NRg = {} #to store the results
    
    for elem, row in tqdm(df_simulations.iterrows(), "Converting the spot number to the number of NR genes for each element randomly distributed", position=0, leave = False, total = len(df_simulations)):
        if elem not in d_spot_NRg.keys():
            d_spot_NRg[elem] = {}
        for n_sim in tqdm(list_sim, desc = "Doing so for each simulation", position = 1, leave = False):
            d_spot_NRg[elem][n_sim] = df_spot_tmp.loc[row[n_sim], "N_NR_genes"]
        
    df_spot_NRg = pd.DataFrame.from_dict(d_spot_NRg, orient = "index")
    
    return df_spot_NRg


def convert_spot_to_number_of_elements_simulated(df_sim, df_all_elements, is_same_elem_type = False): #Change is_same_elem_type to True if your custom group includes TAs
    """
    Estimating the number of randomly distributed elem_of_interest that would colocalize with each TA systems so we can plot it later
    """
    
    list_simulation = df_sim.columns.tolist()
    df_all_elements = df_all_elements.set_index("Element_name")
    list_index_real = df_all_elements.index.tolist()

    if is_same_elem_type == False:
        
        # Get spot numbers for the indices we need
        spot_numbers = df_all_elements.loc[list_index_real, "Spot_number"]
        
        # Use pandas apply with pre-computed spot numbers
        def map_counts(col):
            """
            Small helper function to speed up the conversion of the dataframe
            """
            value_counts = col.value_counts()
            return spot_numbers.map(value_counts).fillna(0).astype(int)
                
        # Apply to all simulation columns at once
        result_df = df_sim[list_simulation].apply(map_counts)
        result_df.index = list_index_real
    
        return result_df


    else :
        #In this case we wants for each element computes the number of same family elements within the same spot and substract 1 to avoid counting it
        #We do not refer to the real data here as we study self spot colocalization
        
        df_sim_res = df_sim.copy(deep = True)
        for col in df_sim_res.columns:
            value_counts = df_sim_res[col].value_counts()
            df_sim_res[col] = df_sim_res[col].map(value_counts) - 1
                    
        return df_sim_res


def map_counts(col, spot_numbers):
    """
    small helper function
    """
    value_counts = col.value_counts()
    return spot_numbers.map(value_counts).fillna(0).astype(int)


def plot_spot_length_distribution_random_group(df_sim_NRg, df_real, outname, group_name):
    """
    """

    real_data = df_real[f"N_NR_genes_in_spot"]
    real_data = real_data.squeeze()
    real_mean = real_data.mean()

    simulation_means = df_sim_NRg.mean(axis=0)
    
    n_above = (simulation_means > real_mean).sum()
    n_below = (simulation_means < real_mean).sum()
    n_equal = (simulation_means == real_mean).sum()
    p_value_lower = (simulation_means <= real_mean).sum() / len(simulation_means)
    p_value_higher = (simulation_means >= real_mean).sum() / len(simulation_means)
    sim_std = simulation_means.std()
    z_score = (real_mean - simulation_means.mean()) / sim_std
    
    plt.figure(figsize=(10, 6))
    sns.histplot(simulation_means, kde=True, bins="auto", color="blue", alpha=0.6, label="Mean distribution (simulations)")
    plt.axvline(real_mean, color="red", linestyle="--", linewidth=2, label=f"Real mean = {real_mean:.4f}")
    
    x_min = min(simulation_means.min(), real_mean)
    x_max = max(simulation_means.max(), real_mean)
    padding = (x_max - x_min) * 0.05  # 5% padding
    plt.xlim(x_min - padding, x_max + padding) #Fixed limits so it is easier to compare them
    
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
    
    
    plt.xlabel(f"Mean number of NRgenes within a given {group_name} associated spot")
    plt.ylabel("N simulations")
    plt.title(f"Mean number of NRgenes within a given {group_name} associated spot")
    plt.legend(loc='upper right')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(outname + "/2-splot_length_custom_group_distribution_simulation.png", dpi = 600)
    plt.close()


def plot_spot_length_distribution_random_group_for_article(df_sim_NRg, df_real, outname, group_name, plot_legend_article_figure = False):
    """
    """
    
    #setting up the color to use for the plot
    group_name_formated = group_name.split('[')[0]        
    if group_name_formated == "TA":
        color_plot = "blue"
    elif group_name_formated == "defense" or group_name_formated == "Defense":
        color_plot = "#a00000"
    elif group_name_formated == "IS":
        color_plot = "#00c000"
    elif group_name_formated == "integrase" or group_name_formated == "Integrase":
        color_plot = "#ad00e3"
    

    real_data = df_real[f"N_NR_genes_in_spot"]
    real_data = real_data.squeeze()
    real_mean = real_data.mean()

    simulation_means = df_sim_NRg.mean(axis=0)
    plt.figure(figsize=(10, 6))

    counts, bins, patches = plt.hist(simulation_means, bins="auto", color=None, alpha=0)
    ax1 = plt.gca()
    ax1.set_ylabel("N simulation", fontsize=14)
    ax1.set_xlabel(f"Mean number of NRgenes associated to a given {group_name}", fontsize=14) 
    ax1.tick_params(axis='both', which='major', labelsize=14)

    ax2 = plt.gca().twinx()
    sns.kdeplot(simulation_means, ax=ax2, alpha=0.3, fill=True, color = color_plot, label=f"Distribution mean simulations")
    ax2.set_ylabel('')
    ax2.tick_params(right=False, labelright=False)
    ax2.spines['right'].set_visible(False)
    
    x_min = min(simulation_means.min(), real_mean)
    x_max = max(simulation_means.max(), real_mean)
    padding = (x_max - x_min) * 0.05  # 5% padding
    plt.xlim(x_min - padding, x_max + padding) #Fixed limits so it is easier to compare them
    
    # Vertical line with the real mean value
    plt.axvline(real_mean, color="red", linestyle="--", linewidth=3, label=f"Real mean = {real_mean:.4f}")    
    
    if plot_legend_article_figure == True:
        
        n_above = (simulation_means > real_mean).sum()
        n_below = (simulation_means < real_mean).sum()
        n_equal = (simulation_means == real_mean).sum()
        p_value_lower = (simulation_means <= real_mean).sum() / len(simulation_means)
        p_value_higher = (simulation_means >= real_mean).sum() / len(simulation_means)
        sim_std = simulation_means.std()
        z_score = (real_mean - simulation_means.mean()) / sim_std

        #Plotting stats values + legend
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
        
        plt.legend(loc='upper right')
        
    plt.tick_params(axis='both', which='major', labelsize=14)
    plt.grid(False)
    plt.tight_layout()
    plt.savefig(outname + "/2.2-plot_length_custom_group_distribution_simulation_for_article.png", dpi = 600)
    plt.close()



def plot_mean_distribution_from_perspective(df_sim, df_real, elem_perspective_name, outname, group_name):
    
    """
    Plot the mean distribution of the real data and simulations from a same element perspectives.
    """

    real_data = df_real[f"N_custom_group_colocalization"]
    

    simulation_means = df_sim.mean(axis=0)
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
    

    plt.figure(figsize=(10, 6))


    sns.histplot(simulation_means, kde=True, bins="auto", color="blue", alpha=0.6, label="Mean distribution (simulations)")
    plt.axvline(real_mean, color="red", linestyle="--", linewidth=2, label=f"Real mean = {real_mean:.4f}")
    

    x_min = min(simulation_means.min(), real_mean)
    x_max = max(simulation_means.max(), real_mean)
    padding = (x_max - x_min) * 0.05  # 5% padding
    plt.xlim(x_min - padding, x_max + padding)
    
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
    
    
    plt.xlabel(f"Mean number of {group_name} associated with a given {elem_perspective_name}")
    plt.ylabel("N simulations")
    plt.title(f"Mean number of {group_name} associated with a given {elem_perspective_name}")
    plt.legend(loc='upper right')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(outname)
    plt.close()
    
    return [group_name, str(real_mean), str(n_below), str(n_equal), str(n_above), str(round(p_value_lower,4)), str(round(p_value_higher,4)), str(round(z_score,4))]



def plot_mean_distribution_from_perspective_for_article(df_sim, df_real, elem_perspective, outname, group_name, plot_legend_article_figure = False):
    
    """
    Plot the mean distribution of the real data and simulations from a same element perspectives.
    """
   
    if elem_perspective == "TA":
        color_plot = "blue"
    elif elem_perspective == "defense" or elem_perspective == "Defense":
        color_plot = "#a00000"
    elif elem_perspective == "IS":
        color_plot = "#00c000"
    elif elem_perspective == "integrase" or elem_perspective == "Integrase":
        color_plot = "#ad00e3"


    real_data = df_real[f"N_custom_group_colocalization"]
    
    simulation_means = df_sim.mean(axis=0)

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
    
    plt.figure(figsize=(10, 6))
    
    counts, bins, patches = plt.hist(simulation_means, bins="auto", color=None, alpha=0)
    ax1 = plt.gca()
    ax1.set_ylabel("N simulation", fontsize=14)
    ax1.set_xlabel(f"Mean number of {group_name} associated with a given {elem_perspective}", fontsize=14) #change here for the name of the custom group
    ax1.tick_params(axis='both', which='major', labelsize=14)

    ax2 = plt.gca().twinx()
    sns.kdeplot(simulation_means, ax=ax2, alpha=0.3, fill=True, color = color_plot)
    ax2.set_ylabel('')
    ax2.tick_params(right=False, labelright=False)
    ax2.spines['right'].set_visible(False)

    plt.axvline(real_mean, color="red", linestyle="--", linewidth=3, label=f"Real mean = {real_mean:.4f}")
    
    x_min = min(simulation_means.min(), real_mean)
    x_max = max(simulation_means.max(), real_mean)
    padding = (x_max - x_min) * 0.05  # 5% padding
    plt.xlim(x_min - padding, x_max + padding)   

    if plot_legend_article_figure == True:
        
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
        
        plt.legend(loc='upper right')


    plt.tick_params(axis='both', which='major', labelsize=14)
    plt.grid(False)
    plt.tight_layout()
    plt.savefig(outname.rsplit(".", 1)[0] + "_for_article.png", dpi = 300)
    plt.close()
    
    return [group_name, str(real_mean), str(n_below), str(n_equal), str(n_above), str(round(p_value_lower,4)), str(round(p_value_higher,4)), str(round(z_score,4))]



def plot_CDF_distribution(outname, df_sim, real_data, group_name):
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
    plt.title(f"CDF cumulative : simulations vs données réelles : {group_name}")
    plt.ylim(0, 100)
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(outname)
    plt.close()



def write_global_res(res2add, outdir, elem_perspective_name):
    """
    """
    res2add.insert(0, elem_perspective_name)
    with open(outdir + f"/4-Global_results_colocalization_{elem_perspective_name}_perspectives_with_custom_group.tsv", "w") as f:
        f.write("Element_perspectives\tCustom_group\tN_members_in_family\treal_mean_value\tn_below_real_mean_value\tn_equal_real_mean_value\tn_over_real_mean_value\tp-value_below_equal_real\tp-value_over_equal_real\tZ-score\n")
        f.write("\t".join(res2add) + "\n")



def fetch_elem_list(file_custom_group):
    """
    """
    
    list_res = []
    with open(file_custom_group, "r") as f:
        for line in f:
            if line.strip() not in list_res: #safeguard to ensure there are no duplicate
                list_res.append(line.strip())
    
    return list_res


def create_dir(outdir):
    """
    """
    
    
    new_list_dir = [outdir]
    for dir in new_list_dir:
        if not os.path.exists(dir):
            os.mkdir(dir)




if __name__ == "__main__":
    import sys
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])