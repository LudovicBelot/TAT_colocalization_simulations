#!/usr/bin/env python

import pandas as pd
import sys
from tqdm import tqdm


def main(infile, outfile):
    """
    Small script that aims to estimates the percentage of each element type (TA/defense/IS/integrase) to colocalize with each element families (e.g. ParE/RelE for TAs)
    """
    
    df_in = pd.read_csv(infile, sep = "\t")
    
    list_element_type = df_in["Element_type"].unique().tolist()
    list_element_families = df_in["Element_family"].unique().tolist()
    list_spot_TAs_no_defense = [x for x in df_in[df_in["Element_type"] == "TA"]["Spot_number"].unique().tolist() if x not in df_in[df_in["Element_type"] == "defense"]["Spot_number"].unique().tolist()]
    list_spot_defense_no_TAs = [x for x in df_in[df_in["Element_type"] == "defense"]["Spot_number"].unique().tolist() if x not in df_in[df_in["Element_type"] == "TA"]["Spot_number"].unique().tolist()]

    df_in[df_in["Element_type"] == "defense"]["Spot_number"].unique().tolist()
    d_out = {}
    
    for family in tqdm(list_element_families, desc = "Determining the % of each elements families to colocalize with a each element type", total = len(list_element_families)):
        family_type = df_in[df_in["Element_family"] == family]["Element_type"].tolist()[0]
        d_out[family] = {"Element_type": family_type}
        n_members_in_family = len(df_in[df_in["Element_family"] == family])
        for elem_type in tqdm(list_element_type, desc = "Doing so for each element type", position = 1, leave = False):
            list_spot_with_elem_type = df_in[df_in["Element_type"] == elem_type]["Spot_number"].unique().tolist()
            n_members_also_in_type_spot = len(df_in[(df_in["Element_family"] == family) & (df_in["Spot_number"].isin(list_spot_with_elem_type))])
            d_out[family][elem_type] = f"{round(n_members_also_in_type_spot/n_members_in_family*100, 2)} % ({n_members_also_in_type_spot}/{n_members_in_family})"

        n_members_with_TA_but_not_defense = len(df_in[(df_in["Element_family"] == family) & (df_in["Spot_number"].isin(list_spot_TAs_no_defense))])
        n_members_with_defense_but_not_TA= len(df_in[(df_in["Element_family"] == family) & (df_in["Spot_number"].isin(list_spot_defense_no_TAs))])
        d_out[family]["Coloc_with_TA_but_not_defense"] = f"{round(n_members_with_TA_but_not_defense/n_members_in_family*100, 2)} % ({n_members_with_TA_but_not_defense}/{n_members_in_family})"
        d_out[family]["Coloc_with_defense_but_not_TA"] = f"{round(n_members_with_defense_but_not_TA/n_members_in_family*100, 2)} % ({n_members_with_defense_but_not_TA}/{n_members_in_family})"
        
    df_out = pd.DataFrame.from_dict(d_out, orient = "index")
    df_out.to_csv(outfile, sep = "\t")
        



if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])