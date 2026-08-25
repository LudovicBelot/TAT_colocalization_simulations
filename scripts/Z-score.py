#!/usr/bin/env python

import pandas as pd
import numpy as np
from tqdm import tqdm

"""
Small script which compute the Z score between each simulation and the real data.
cl example : python script/Z-score.py test/0-39_genomes_Ph25/TAT_defense_results/simulation/data/3.0-All_IS_families_colocalization_from_TA_perspectives.tsv test/0-39_genomes_Ph25/TAT_defense_results/simulation/data
"""

def main(res_file, res_dir):
    """
    """

    df_res = pd.read_csv(res_file, sep = "\t")
    list_elements_tested = df_res["IS_family"].tolist()
    
    n = 0
    list_mean_sim = []
    list_standard_dev = []
    list_Z_score = []
    for elem in tqdm(list_elements_tested, total = len(list_elements_tested)):
        n += 1
        df_sim = pd.read_csv(res_dir + "/" + f"3.{n}.2-Colocalization_perspectives_from_TA_point_of_view_with_{elem}.tsv", sep = "\t", index_col = 0)
        mean_each_sim = df_sim.mean()
        overall_mean_sim = np.mean(mean_each_sim)
        standard_dev = np.std(mean_each_sim, ddof=1)
        real_mean = df_res[df_res["IS_family"] == elem]["real_mean_value"].tolist()[0]
        Z_score = (real_mean - overall_mean_sim)/standard_dev
        
        list_mean_sim.append(overall_mean_sim)
        list_standard_dev.append(standard_dev)
        list_Z_score.append(Z_score)
    
    
    df_res["mean_simulations"] = list_mean_sim
    df_res["standard_deviation_simulations"] = list_standard_dev
    df_res["Z-score"] = list_Z_score
    
    df_res.to_csv(res_file, sep = "\t")



if __name__ == "__main__":
    import sys
    main(sys.argv[1], sys.argv[2])