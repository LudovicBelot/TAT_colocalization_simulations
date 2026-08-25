import subprocess
import pandas as pd
import numpy as np
from tqdm import tqdm
import os
import multiprocessing
from Bio import SeqIO
import copy

def main(indir, outdir, list_to_search = ["TA", "defense", "IS", "integrase"], custom_TA = False, n_cpu = 1):
    """
    Main function that will searches for canonical type II TAs, defense systems, IS and integrases.
    Will also determine how many systems are redundants and only keep representatives ones.
    Using custom scripts/HMM profiles for TA search.
    Defense finder for defense systems (https://github.com/mdmparis/defense-finder)
    HMMsearch for IS and integrases
    
    Note that I re-used these functions from another project, therefore they are not optimized in this framework
    """
    
    for system in list_to_search:
        if system == "TA":
            df_TAsearch = search_TA(indir, outdir, n_cpu = n_cpu)
            df_TAsearch.to_csv(f"{outdir}/TAT_defense_results/TA_defense_search/0-TAsearch_all.tsv", sep = "\t", index = False)
            screened_TA2spot(df_TAsearch, outdir)

        elif system == "defense":
            search_defense(outdir, n_cpu = n_cpu)

        elif system == "IS":
            search_IS(outdir, n_cpu = n_cpu)

        elif system == "integrase":
            search_integrase(outdir, n_cpu = n_cpu)

    if custom_TA != False:
        df_custom_TA = pd.read_csv(custom_TA, sep = "\t")
        custom_TA2spot(df_custom_TA, outdir)
        list_to_search.insert(0,"TA")
    
    #Combine all results together
    df_final_results = parse_spot_results(outdir, list_to_search)
    #Compute the % of each element type to spot-colocalize with the others
    compute_colocalization(outdir, df_final_results)





def custom_TA2spot(df_custom_TA, outdir):
    """
    Function which based on the coordinates of the provided custom TAs, attributes them to a spot.
    Note that it works only with the reference genome.
    """
    
    df_ref_intervals = pd.read_csv(f"{outdir}/TAT_defense_results/core_mapping/0-Ref_genome_intervals.tsv", sep = "\t")
    
    d_res_custom = {}
    d_res_to_continue_with = {}
    
    for TA_row in tqdm(df_custom_TA.iterrows(), desc = "Associating each custom TA to a spot"):
        df_tmp_interval = df_ref_intervals[(df_ref_intervals["left_c"] <= TA_row[1]["left_c"]) & (df_ref_intervals["right_c"] >= TA_row[1]["right_c"])]
        
        if len(df_tmp_interval) == 0: #means coordinates at least partially overlaps with core genes
            d_res_custom[TA_row[1]["Ref_TA"]] = {"Custom_TA": TA_row[1]["Ref_TA"],
                                          "Interval": "Core"}
        
        elif len(df_tmp_interval) == 1: #Then fetch the correct interval associated
            d_res_custom[TA_row[1]["Ref_TA"]] = {"Custom_TA": TA_row[1]["Ref_TA"],
                                    "Interval": df_tmp_interval["interval_number"].tolist()[0]}
            
            d_res_to_continue_with[TA_row[1]["Ref_TA"]] = {"Ref_TA": TA_row[1]["Ref_TA"],
                                                           "TA_family": TA_row[1]["Family"] if "Family" in df_custom_TA.columns else "Custom",
                                                           "SPOT": df_tmp_interval["interval_number"].tolist()[0],
                                                           "Representative": "True" 
                                                            }
        
        elif len(df_tmp_interval) > 1: #means the given custom TA is spread on more than one interval with the given coordinates. It should not happen normally
            d_res_custom[TA_row[1]["Ref_TA"]] = {"Custom_TA": TA_row[1]["Ref_TA"],
                                    "Interval": "Multiple_interval"}

    df_res_custom = pd.DataFrame.from_dict(d_res_custom, orient = "index")
    df_res_to_continue_with = pd.DataFrame.from_dict(d_res_to_continue_with, orient = "index")
    
    d_res_per_spot = {}
    for spot in df_res_to_continue_with["SPOT"].unique():
        df_tmp = df_res_to_continue_with[df_res_to_continue_with["SPOT"] == spot]
        d_res_per_spot[spot] = {"Spot_number": spot,
                                "TAlist":(",").join(df_tmp["Ref_TA"].tolist()),
                                "TA_families":(",").join(df_tmp["TA_family"].tolist())
                                }
        
        for family in df_res_to_continue_with["TA_family"]:
            d_res_per_spot[spot][family] = len(df_tmp[df_tmp["TA_family"] == family])

        d_res_per_spot[spot]["n_TAs"] = len(df_tmp)
    
    df_res_per_spot =  pd.DataFrame.from_dict(d_res_per_spot, orient = "index")
    
    df_res_custom.to_csv(f"{outdir}/TAT_defense_results/TA_defense_search/1.0-Custom_TA_to_spots.tsv", sep = "\t",index = False)
    df_res_to_continue_with.to_csv(f"{outdir}/TAT_defense_results/TA_defense_search/1.1-TA_spot_pangenome.tsv", sep = "\t", index = False)
    df_res_per_spot.to_csv(f"{outdir}/TAT_defense_results/TA_defense_search/1.2-TAlist_per_spot_NR.tsv", sep = "\t", index = False)




def search_TA(indir, outdir, n_cpu = 1):
    """
    Search for canonical type II TAs using custom scripts/HMM profiles.
    """

    list_file_processed = []

    for file in tqdm([x for x in os.listdir(f"{indir}/Proteins") if x.endswith(".prt")], desc = "Running Hmmsearch for TAs detection"):
        #Hmmsearch toxins
        subprocess.call(f"hmmsearch --cpu {n_cpu} --tblout {outdir}/TAT_defense_results/tmp/TA_defense_search/{file}_toxins.tsv -o {outdir}/TAT_defense_results/tmp/TA_defense_search/trash.tsv {os.path.expanduser('~')}/2025-01-08_all_redo/4-TAT_defense/script/db/TA/TypeII_toxins.hmm {indir}/Proteins/{file}",
                        shell = True)

        #Hmmsearch Antitoxins
        subprocess.call(f"hmmsearch --cpu {n_cpu} --tblout {outdir}/TAT_defense_results/tmp/TA_defense_search/{file}_antitoxins.tsv -o {outdir}/TAT_defense_results/tmp/TA_defense_search/trash.tsv {os.path.expanduser('~')}/2025-01-08_all_redo/4-TAT_defense/script/db/TA/TypeII_antitoxins.hmm {indir}/Proteins/{file}",
                        shell = True)

        list_file_processed.append(file)


    list_multi = [[x, indir, outdir] for x in np.array_split(list_file_processed, n_cpu) if x.size != 0]
    
    with multiprocessing.Pool(n_cpu) as pool:
        list_df_TA = pool.starmap(func= TAsearch_associate_TandA, iterable= list_multi)
    
    df_all_TAs = pd.concat(list_df_TA)
    
    return df_all_TAs



def TAsearch_associate_TandA(list_files, indir, outdir, space_interTA = 350):
    """
    Function that determines which toxins and antitoxins are associated (within 350bp by default).
    Note that I re-used these functions from another project, therefore they are not optimized in this framework
    """
    
    hmmsearch_columns = ["target name", "accession", "query_name", "accession_2", "E-value", "score", "biais"]
    d_res = {}
    
    for file in tqdm(list_files, desc= f"Associating each toxins to their potential antitoxins: ", leave = False):
    
        df_tox = pd.read_csv(f"{outdir}/TAT_defense_results/tmp/TA_defense_search/{file}_toxins.tsv", sep ="\s+", comment = "#", usecols = [0,1,2,3,4,5,6], names = hmmsearch_columns)
        df_antitox = pd.read_csv(f"{outdir}/TAT_defense_results/tmp/TA_defense_search/{file}_antitoxins.tsv", sep ="\s+", comment = "#", usecols = [0,1,2,3,4,5,6], names = hmmsearch_columns)

        #Keeping only the best hits
        df_tox = df_tox.sort_values(by = "E-value").drop_duplicates(subset = ["target name"])
        df_antitox = df_antitox.sort_values(by = "E-value").drop_duplicates(subset = ["target name"])

        df_gff = pd.read_csv(f"{indir}/gff3/{file.rsplit('.',1)[0]}.gff", sep = "\t", comment = "#", names= ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"])
        df_gff["attributes"] = df_gff["attributes"].map(lambda x: x.split("ID=",1)[1].split(";")[0])
        df_gff["short_id"] = df_gff["attributes"].map(lambda x: x.rsplit(".",1)[-1])

        #adding contig/left_c/right_c for every hit
        df_tox = hmmsearch_add_cds_features(df_tox, df_gff)
        df_antitox = hmmsearch_add_cds_features(df_antitox, df_gff)
        
        n_index = 0
        list_antitoxins_already_attributed = []
        
        if df_tox.empty == True or df_antitox.empty == True :
            continue
        
        for row in df_tox.iterrows(): #checking for each toxin if we can associate them with an antitoxin (from the best toxin evalue to the worst)
            
            if df_antitox[(df_antitox["seqid"] == row[1]["seqid"]) & 
                            (df_antitox["strand"] == row[1]["strand"]) &
                            (df_antitox["end"] >= int(row[1]["start"])-space_interTA) &
                            (df_antitox["start"] <= int(row[1]["end"])+space_interTA) &
                            (~df_antitox["target name"].isin(list_antitoxins_already_attributed))
                            ].empty == False:
                
                
                df_tmp_antitox = df_antitox[(df_antitox["seqid"] == row[1]["seqid"]) & 
                            (df_antitox["strand"] == row[1]["strand"]) &
                            (df_antitox["end"] >= int(row[1]["start"])-space_interTA) &
                            (df_antitox["start"] <= int(row[1]["end"])+space_interTA) &
                            (~df_antitox["target name"].isin(list_antitoxins_already_attributed))
                            ]

                d_res[f"{file}_{n_index}"] = {"genome": f"{file.rsplit('.',1)[0]}",
                                    "contig":row[1]["seqid"],
                                    "toxin_id": row[1]["target name"],
                                    "toxin_hit": row[1]["query_name"],
                                    "E-value_tox":row[1]["E-value"],
                                    "antitoxin_id": df_tmp_antitox["target name"].tolist()[0],
                                    "antitoxin_hit": df_tmp_antitox["query_name"].tolist()[0],
                                    "E-value_antitox": df_tmp_antitox["E-value"].tolist()[0],
                                    "start_tox":row[1]["start"],
                                    "end_tox":row[1]["end"],
                                    "start_antitox":df_tmp_antitox["start"].tolist()[0],
                                    "end_antitox":df_tmp_antitox["end"].tolist()[0]
                                    }

                list_antitoxins_already_attributed.append(df_tmp_antitox["target name"].tolist()[0])
                n_index +=1
    
    
    if d_res != {}:
        return pd.DataFrame.from_dict(d_res, orient = "index")
    else:
        return pd.DataFrame()


def hmmsearch_add_cds_features(df2change, df_gff):

    df2change["short_id"] = df2change.apply(lambda x: f"{x['target name'].rsplit('.',1)[-1]}", axis = 1)
    df2return = pd.merge(df2change, df_gff, on = "short_id")
    return df2return



def screened_TA2spot(df_TAsearch, outdir):
    """
    Function which determine whether each found TAs are part of the spot pangenome, the core genome or within a breakpoint interval.
    Then, for each spot it will run mmseqs and count/keep only NR_TAs in each spot (0.8 cov/0.8 id, one gene in common between two system is enough to consider them redundant).
    """
    
    df_intervals = pd.read_csv(f"{outdir}/TAT_defense_results/core_mapping/2-all_genomes_intervals.tsv", sep = "\t", comment = "#")    
    df_family_index = pd.read_csv(f"{os.path.expanduser('~')}/2025-01-08_all_redo/4-TAT_defense/script/db/TA/TA_family.idx", sep = ",", comment = "#")
    
    d_res_perTA = {}
    d_res_perSPOT = {}
    with open(f"{outdir}/TAT_defense_results/tmp/list_allcore_genes.lst", "r") as f:
        for line in f:
            list_core_genes = line.split(",")



    for row in df_TAsearch.iterrows():
        if row[1]["toxin_id"] in list_core_genes or row[1]["antitoxin_id"] in list_core_genes:
            d_res_perTA[f"{row[1]['toxin_id']}-{row[1]['antitoxin_id']}"] = {"TA": f"{row[1]['toxin_id']}-{row[1]['antitoxin_id']}",
                                    "family": df_family_index[df_family_index["hmm_name"] == row[1]['toxin_hit']]["TA_family"].tolist()[0] if row[1]["toxin_hit"] in df_family_index["hmm_name"].tolist() else "Other",
                                    "SPOT": "CORE"}
        
        elif df_intervals[df_intervals[row[1]["genome"]].str.contains(str(row[1]["toxin_id"]), na=False)].empty == False:
            tmp_index = df_intervals[df_intervals[row[1]["genome"]].str.contains(str(row[1]["toxin_id"]), na=False)].index.tolist()[0]
            
            if str(row[1]["antitoxin_id"]) in df_intervals.loc[tmp_index, row[1]["genome"]].split(","):
                d_res_perTA[f"{row[1]['toxin_id']}-{row[1]['antitoxin_id']}"] = {"TA": f"{row[1]['toxin_id']}-{row[1]['antitoxin_id']}",
                                        "family": df_family_index[df_family_index["hmm_name"] == row[1]['toxin_hit']]["TA_family"].tolist()[0] if row[1]["toxin_hit"] in df_family_index["hmm_name"].tolist() else "Other",
                                        "SPOT": df_intervals.loc[tmp_index,"interval_number"]}
            
            else :
                d_res_perTA[f"{row[1]['toxin_id']}-{row[1]['antitoxin_id']}"] = {"TA": f"{row[1]['toxin_id']}-{row[1]['antitoxin_id']}",
                                        "family": df_family_index[df_family_index["hmm_name"] == row[1]['toxin_hit']]["TA_family"].tolist()[0] if row[1]["toxin_hit"] in df_family_index["hmm_name"].tolist() else "Other",
                                        "SPOT": "BI"}
        
        else:
            d_res_perTA[f"{row[1]['toxin_id']}-{row[1]['antitoxin_id']}"] = {"TA": f"{row[1]['toxin_id']}-{row[1]['antitoxin_id']}",
                                    "family": df_family_index[df_family_index["hmm_name"] == row[1]['toxin_hit']]["TA_family"].tolist()[0] if row[1]["toxin_hit"] in df_family_index["hmm_name"].tolist() else "Other",
                                    "SPOT": "BI"}

    # Now getting the list of all TAs in each SPOT to determine which are redundants or not
    for v in d_res_perTA.values():
        if v["SPOT"] not in ["BI","CORE"]:
            if v["SPOT"] not in d_res_perSPOT.keys():
                d_res_perSPOT[v["SPOT"]] = {"SPOT":v["SPOT"], "TAlist":v["TA"]}
            else:
                d_res_perSPOT[v["SPOT"]]["TAlist"] += "," + str(v["TA"])
    
      
    all_prt_seqIO = SeqIO.index(f"{outdir}/TAT_defense_results/tmp/all_prt_cat.prt", "fasta")
    final_list2keep = []
    
    for spot, v in tqdm(d_res_perSPOT.items(), desc = "Running MMseqs to remove spot redundants TAs"):
        d_tmp_family = {}
        for TA in v["TAlist"].split(","):
            if d_res_perTA[TA]["family"] not in d_tmp_family.keys():
                d_tmp_family[d_res_perTA[TA]["family"]] = [TA]
            else:
                d_tmp_family[d_res_perTA[TA]["family"]] += [TA]
        
        #running mmseqs to determine if there are redundants systems
        for TA_type in d_tmp_family.keys():
            if len(d_tmp_family[TA_type]) <= 1:
                for x in d_tmp_family[TA_type]:
                    final_list2keep.append(x)
                continue
            
            tmp_str_tox = ""
            tmp_str_antitox = ""
            for TA2test in d_tmp_family[TA_type]:
                tmp_str_tox += f">{TA2test.split('-')[0]}\n{all_prt_seqIO[TA2test.split('-')[0]].seq}\n"
                tmp_str_antitox += f">{TA2test.split('-')[1]}\n{all_prt_seqIO[TA2test.split('-')[1]].seq}\n"
                
    
            if  os.path.exists(f"{outdir}/TAT_defense_results/tmp/TA_defense_search/toxins_seq.tmp_TA_cluster4mmseqs.faa") == True:
                os.remove(f"{outdir}/TAT_defense_results/tmp/TA_defense_search/toxins_seq.tmp_TA_cluster4mmseqs.faa") 
            with open(f"{outdir}/TAT_defense_results/tmp/TA_defense_search/toxins_seq.tmp_TA_cluster4mmseqs.faa", "w") as f:
                f.write(tmp_str_tox)
            
            if  os.path.exists(f"{outdir}/TAT_defense_results/tmp/TA_defense_search/antitoxins_seq.tmp_TA_cluster4mmseqs.faa") == True:
                os.remove(f"{outdir}/TAT_defense_results/tmp/TA_defense_search/antitoxins_seq.tmp_TA_cluster4mmseqs.faa") 
            with open(f"{outdir}/TAT_defense_results/tmp/TA_defense_search/antitoxins_seq.tmp_TA_cluster4mmseqs.faa", "w") as f:
                f.write(tmp_str_antitox)
        
            subprocess.call(f"mmseqs easy-cluster {outdir}/TAT_defense_results/tmp/TA_defense_search/toxins_seq.tmp_TA_cluster4mmseqs.faa {outdir}/TAT_defense_results/tmp/TA_defense_search/spot_toxins {outdir}/TAT_defense_results/tmp/TA_defense_search --min-seq-id 0.8 -v 0 --remove-tmp-files", shell = True)
            df_mmseqs_toxin = pd.read_csv(f"{outdir}/TAT_defense_results/tmp/TA_defense_search/spot_toxins_cluster.tsv", sep = "\t", names = ["query","subject"])
            list_rep_toxin = df_mmseqs_toxin["query"].drop_duplicates().tolist()
            
            subprocess.call(f"mmseqs easy-cluster {outdir}/TAT_defense_results/tmp/TA_defense_search/antitoxins_seq.tmp_TA_cluster4mmseqs.faa {outdir}/TAT_defense_results/tmp/TA_defense_search/spot_antitoxins {outdir}/TAT_defense_results/tmp/TA_defense_search --min-seq-id 0.8 -v 0 --remove-tmp-files", shell = True)
            df_mmseqs_antitoxin = pd.read_csv(f"{outdir}/TAT_defense_results/tmp/TA_defense_search/spot_antitoxins_cluster.tsv", sep = "\t", names = ["query","subject"])
            list_rep_antitoxin = df_mmseqs_antitoxin["query"].drop_duplicates().tolist()

            for TA2test in d_tmp_family[TA_type]:
                if TA2test.split('-')[0] in list_rep_toxin and TA2test.split('-')[1] in list_rep_antitoxin:
                    final_list2keep.append(TA2test)


    for TA_key in d_res_perTA.keys():
        if d_res_perTA[TA_key]["TA"] in final_list2keep:
            d_res_perTA[TA_key]["Representative"] = "True"
        else :
            d_res_perTA[TA_key]["Representative"] = "False"

    d_res_perSPOT = {} #redoing it so we keep only representatives TAs
    for v in d_res_perTA.values():
        if v["Representative"] == "True":
            if v["SPOT"] not in d_res_perSPOT.keys():
                d_res_perSPOT[v["SPOT"]] = {"Spot_number":v["SPOT"], 
                                            "TAlist":v["TA"], 
                                            "TA_families": v["family"]
                                            }
            else:
                d_res_perSPOT[v["SPOT"]]["TAlist"] += ","+str(v["TA"])
                d_res_perSPOT[v["SPOT"]]["TA_families"] += "," + str(v["family"])
                

            #added these lines to also count the number of each TA family within the spot
            if v["family"] not in d_res_perSPOT[v["SPOT"]].keys():
                d_res_perSPOT[v["SPOT"]][v["family"]] = 1
            
            else :
                d_res_perSPOT[v["SPOT"]][v["family"]] += 1

    
    df_res_perTA = pd.DataFrame.from_dict(d_res_perTA, orient = "index")
    df_res_perSPOT= pd.DataFrame.from_dict(d_res_perSPOT, orient = "index")
    df_res_perSPOT["n_TAs"] = df_res_perSPOT.apply(lambda x: len(x["TAlist"].split(",")), axis = 1)
    df_res_perTA.to_csv(f"{outdir}/TAT_defense_results/TA_defense_search/1.0-TA_spot_screening_redundants.tsv", index = False, sep = "\t")
    
    #formating some columns names for further operations
    df_res_perTA.rename(columns = {"TA":"Ref_TA", "family":"TA_family"}, inplace = True)
    df_res_perTA.to_csv(f"{outdir}/TAT_defense_results/TA_defense_search/1.1-TA_spot_pangenome.tsv", sep = "\t", index = False)
    df_res_perSPOT.to_csv(f"{outdir}/TAT_defense_results/TA_defense_search/1.2-TAlist_per_spot_NR.tsv", index = False, sep = "\t")

    return df_res_perSPOT


""""""


def search_IS(outdir, n_cpu = 1):
    """
    Function which searches for IS transposase genes within the complete proteome of the studied phylogeny.
    Then runs mmseqs to keep only non-redundant IS in each spot.
    
    """

    df_SP = pd.read_csv(f"{outdir}/TAT_defense_results/core_mapping/2-all_genomes_intervals.tsv", sep = "\t")
    #we concatenate all interval into one column to ease latter operations
    df_SP["all"] = df_SP[df_SP.columns[1:]].apply(lambda x: ",".join(x.dropna().astype(str)), axis = 1)

    with open(f"{outdir}/TAT_defense_results/tmp/list_allcore_genes.lst", "r") as f:
        for line in f:
            list_core_genes = line.split(",")
            

    #Hmmsearch for IS
    print("Using Hmmsearch to search for IS transposase genes")
    subprocess.call(f"hmmsearch --cpu {n_cpu} --tblout {outdir}/TAT_defense_results/tmp/TA_defense_search/hmmsearch_IS.tsv -o {outdir}/TAT_defense_results/tmp/TA_defense_search/trash.tsv {os.path.expanduser('~')}/2025-01-08_all_redo/4-TAT_defense/script/db/mge/23-08-23_ISescan_clusters_converted.hmm {outdir}/TAT_defense_results/tmp/all_prt_cat.prt", 
                    shell = True)

    hmmsearch_columns = ["target name", "accession", "query_name", "accession_2", "E-value", "score", "biais"]
    df_IS_hmmsearch = pd.read_csv(f"{outdir}/TAT_defense_results/tmp/TA_defense_search/hmmsearch_IS.tsv", sep ="\s+", comment = "#", usecols = [0,1,2,3,4,5,6], names = hmmsearch_columns)
    df_IS_hmmsearch = df_IS_hmmsearch.sort_values(by = "E-value")
    df_IS_hmmsearch = df_IS_hmmsearch[df_IS_hmmsearch["E-value"] <= 0.001].reset_index(drop=True)
    df_IS_hmmsearch = df_IS_hmmsearch.drop_duplicates(subset= "target name")

    #create two tsv files, first with each mge hit associated to a spot
    #and the second one with all hits in BI or in core
    list_IS_res = []
    #dividing the dataframe to allow multi cpu parsing of the results
    list_iterables4multi_cpu = []
    df_IS_hmmsearch_splitted = [x for x in np.array_split(df_IS_hmmsearch, n_cpu) if x.empty == False]
    for i in range(0,len(df_IS_hmmsearch_splitted)):
        list_iterables4multi_cpu.append([df_IS_hmmsearch_splitted[i], df_SP, list_core_genes, i])

    with multiprocessing.Pool(n_cpu) as pool:
        list_IS_res = [[x[0],x[1]] for x in pool.starmap(func = mge2SP, iterable= list_iterables4multi_cpu)]
    
    d_IS_res = {}
    d_IS_outofspot = {}

    for i in list_IS_res:
        d_IS_res = d_IS_res | i[0]
        d_IS_outofspot = d_IS_outofspot | i[1]
    
    df_IS_res = pd.DataFrame.from_dict(d_IS_res, orient = "index")
    df_IS_outofspot = pd.DataFrame.from_dict(d_IS_outofspot, orient = "index")
    
    #using mmseqs2 to remove redundants IS hits (in case they are within the spot genome)
    seqIO_all_proteome = SeqIO.index(f"{outdir}/TAT_defense_results/tmp/all_prt_cat.prt", "fasta")

    list_genes2keep = []
    for spot in tqdm(df_IS_res["spot_number"].drop_duplicates().tolist(), desc = "Removing spot redundant IS hit"):
        if len(df_IS_res[df_IS_res["spot_number"] == spot]["gene_name"]) == 1:
            list_genes2keep.append(df_IS_res[df_IS_res["spot_number"] == spot]["gene_name"].tolist()[0])
            continue

        #creating a tmp file to mmseqs2
        tmp_fasta4mmseq = ""
        for gene in df_IS_res[df_IS_res["spot_number"] == spot]["gene_name"].tolist():
            tmp_fasta4mmseq += f">{gene}\n{seqIO_all_proteome[gene].seq}\n"

        if os.path.exists(f"{outdir}/TAT_defense_results/tmp/TA_defense_search/tmp_IS_SPcluster4mmseqs.faa") == True:
            os.remove(f"{outdir}/TAT_defense_results/tmp/TA_defense_search/tmp_IS_SPcluster4mmseqs.faa")  
     
        with open (f"{outdir}/TAT_defense_results/tmp/TA_defense_search/tmp_IS_SPcluster4mmseqs.faa", "w") as f:
            f.write(tmp_fasta4mmseq)

        #now doing the mmseqs
        subprocess.call(f"mmseqs easy-cluster {outdir}/TAT_defense_results/tmp/TA_defense_search/tmp_IS_SPcluster4mmseqs.faa {outdir}/TAT_defense_results/tmp/TA_defense_search/IS_in_SP {outdir}/TAT_defense_results/tmp/TA_defense_search --min-seq-id 0.8 -v 0 --remove-tmp-files", shell = True)
        df_tmp = pd.read_csv(f"{outdir}/TAT_defense_results/tmp/TA_defense_search/IS_in_SP_cluster.tsv", sep = "\t", names=["query","target"])
        
        if df_tmp[df_tmp["query"] != df_tmp["target"]].empty:
            for i in df_tmp["query"].tolist():
                list_genes2keep.append(i)
        else :
            tmp_list_query = df_tmp["query"].drop_duplicates().tolist()            
            #Adding all genes but checking at the same time if multiples copies exists within a single genome, will consider the max number of copies in a single genome.
            for i in tmp_list_query:
                list_copies = max_genes_copy_in_a_single_genome(df_tmp[df_tmp["query"] == i]["target"].tolist())
                list_genes2keep += list_copies


    df_IS_final = df_IS_res[df_IS_res["gene_name"].isin(list_genes2keep)]
    df_IS_final.to_csv(f"{outdir}/TAT_defense_results/TA_defense_search/3.1-IS_spot_pangenome.tsv", sep= "\t", index = False)
    df_IS_outofspot.to_csv(f"{outdir}/TAT_defense_results/TA_defense_search/3.2-IS_outofspotpangenome.tsv", sep= "\t", index = False)

    d_IS_per_spot = {}
    for row in df_IS_final.iterrows():
        if row[1]["spot_number"] not in d_IS_per_spot.keys():
            d_IS_per_spot[row[1]["spot_number"]] = {"Spot_number": row[1]["spot_number"], "n_IS": 1}
        else :
            d_IS_per_spot[row[1]["spot_number"]]["n_IS"] += 1
    
    df_IS_per_spot = pd.DataFrame.from_dict(d_IS_per_spot, orient = "index")
    df_IS_per_spot.to_csv(f"{outdir}/TAT_defense_results/TA_defense_search/3.3-IS_per_spot_pangenome.tsv", sep= "\t", index = False)
    
    return df_IS_per_spot


def search_integrase(outdir, n_cpu = 1):
    """
    Function which searches for integrase genes within the complete proteome of the studied phylogeny.
    Then runs mmseqs to keep only non-redundant IS in each spot.
    
    """

    df_SP = pd.read_csv(f"{outdir}/TAT_defense_results/core_mapping/2-all_genomes_intervals.tsv", sep = "\t")
    #we concatenate all interval into one column to ease latter operations
    df_SP["all"] = df_SP[df_SP.columns[1:]].apply(lambda x: ",".join(x.dropna().astype(str)), axis = 1)

    with open(f"{outdir}/TAT_defense_results/tmp/list_allcore_genes.lst", "r") as f:
        for line in f:
            list_core_genes = line.split(",")
            

    #Hmmsearch for integrases (Tyr recombinase pfam : PF00589, and Ser recombinases pfam: PF00239 / PF07508)
    print("Using Hmmsearch to search for integrase genes")
    subprocess.call(f"hmmsearch --cpu {n_cpu} --tblout {outdir}/TAT_defense_results/tmp/TA_defense_search/hmmsearch_integrase.tsv -o {outdir}/TAT_defense_results/tmp/TA_defense_search/trash.tsv {os.path.expanduser('~')}/2025-01-08_all_redo/4-TAT_defense/script/db/mge/23-05-25_integrases_recombinases.hmm {outdir}/TAT_defense_results/tmp/all_prt_cat.prt", 
                    shell = True)

    hmmsearch_columns = ["target name", "accession", "query_name", "accession_2", "E-value", "score", "biais"]
    df_integrase_hmmsearch = pd.read_csv(f"{outdir}/TAT_defense_results/tmp/TA_defense_search/hmmsearch_integrase.tsv", sep ="\s+", comment = "#", usecols = [0,1,2,3,4,5,6], names = hmmsearch_columns)
    df_integrase_hmmsearch = df_integrase_hmmsearch.sort_values(by = "E-value")
    df_integrase_hmmsearch = df_integrase_hmmsearch[df_integrase_hmmsearch["E-value"] <= 0.001].reset_index(drop=True)
    df_integrase_hmmsearch = df_integrase_hmmsearch.drop_duplicates(subset= "target name")

    #create two tsv files, first with each mge hit associated to a spot
    #and the second one with all hits in BI or in core
    list_integrase_res = []
    #dividing the dataframe to allow multi cpu parsing of the results
    list_iterables4multi_cpu = []
    df_integrase_hmmsearch_splitted = [x for x in np.array_split(df_integrase_hmmsearch, n_cpu) if x.empty == False]
    for i in range(0,len(df_integrase_hmmsearch_splitted)):
        list_iterables4multi_cpu.append([df_integrase_hmmsearch_splitted[i], df_SP, list_core_genes, i])

    with multiprocessing.Pool(n_cpu) as pool:
        list_integrase_res = [[x[0],x[1]] for x in pool.starmap(func = mge2SP, iterable= list_iterables4multi_cpu)]
    
    d_integrase_res = {}
    d_integrase_outofspot = {}

    for i in list_integrase_res:
        d_integrase_res = d_integrase_res | i[0]
        d_integrase_outofspot = d_integrase_outofspot | i[1]
    
    df_integrase_res = pd.DataFrame.from_dict(d_integrase_res, orient = "index")
    df_integrase_outofspot = pd.DataFrame.from_dict(d_integrase_outofspot, orient = "index")
    
    #using mmseqs2 to remove redundants IS hits (in case they are within the spot genome)
    seqIO_all_proteome = SeqIO.index(f"{outdir}/TAT_defense_results/tmp/all_prt_cat.prt", "fasta")

    list_genes2keep = []
    for spot in tqdm(df_integrase_res["spot_number"].drop_duplicates().tolist(), desc = "Removing spot redundant integrase gene hits"):
        if len(df_integrase_res[df_integrase_res["spot_number"] == spot]["gene_name"]) == 1:
            list_genes2keep.append(df_integrase_res[df_integrase_res["spot_number"] == spot]["gene_name"].tolist()[0])
            continue

        #creating a tmp file to mmseqs2
        tmp_fasta4mmseq = ""
        for gene in df_integrase_res[df_integrase_res["spot_number"] == spot]["gene_name"].tolist():
            tmp_fasta4mmseq += f">{gene}\n{seqIO_all_proteome[gene].seq}\n"

        if os.path.exists(f"{outdir}/TAT_defense_results/tmp/TA_defense_search/tmp_integrase_SPcluster4mmseqs.faa") == True:
            os.remove(f"{outdir}/TAT_defense_results/tmp/TA_defense_search/tmp_integrase_SPcluster4mmseqs.faa")  
     
        with open (f"{outdir}/TAT_defense_results/tmp/TA_defense_search/tmp_integrase_SPcluster4mmseqs.faa", "w") as f:
            f.write(tmp_fasta4mmseq)

        #now doing the mmseqs
        subprocess.call(f"mmseqs easy-cluster {outdir}/TAT_defense_results/tmp/TA_defense_search/tmp_integrase_SPcluster4mmseqs.faa {outdir}/TAT_defense_results/tmp/TA_defense_search/integrase_in_SP {outdir}/TAT_defense_results/tmp/TA_defense_search --min-seq-id 0.8 -v 0 --remove-tmp-files", shell = True)
        df_tmp = pd.read_csv(f"{outdir}/TAT_defense_results/tmp/TA_defense_search/integrase_in_SP_cluster.tsv", sep = "\t", names=["query","target"])
        
        if df_tmp[df_tmp["query"] != df_tmp["target"]].empty:
            for i in df_tmp["query"].tolist():
                list_genes2keep.append(i)
        else :
            tmp_list_query = df_tmp["query"].drop_duplicates().tolist()            
            #Adding all genes but checking at the same time if multiples copies exists within a single genome, will consider the max number of copies in a single genome.
            for i in tmp_list_query:
                list_copies = max_genes_copy_in_a_single_genome(df_tmp[df_tmp["query"] == i]["target"].tolist())
                list_genes2keep += list_copies


    df_integrase_final = df_integrase_res[df_integrase_res["gene_name"].isin(list_genes2keep)]
    df_integrase_final.to_csv(f"{outdir}/TAT_defense_results/TA_defense_search/4.1-Integrase_spot_pangenome.tsv", sep= "\t", index = False)
    df_integrase_outofspot.to_csv(f"{outdir}/TAT_defense_results/TA_defense_search/4.2-Integrase_outofspotpangenome.tsv", sep= "\t", index = False)

    d_integrase_per_spot = {}
    for row in df_integrase_final.iterrows():
        if row[1]["spot_number"] not in d_integrase_per_spot.keys():
            d_integrase_per_spot[row[1]["spot_number"]] = {"Spot_number": row[1]["spot_number"], "n_integrases": 1}
        else :
            d_integrase_per_spot[row[1]["spot_number"]]["n_integrases"] += 1
    
    df_integrase_per_spot = pd.DataFrame.from_dict(d_integrase_per_spot, orient = "index")
    df_integrase_per_spot.to_csv(f"{outdir}/TAT_defense_results/TA_defense_search/4.3-Integrase_per_spot_pangenome.tsv", sep= "\t", index = False)
    
    return df_integrase_per_spot






def mge2SP(df_mge, df_SP, list_core_genes, index_cpu):

    d_tmp_res = {}
    d_tmp_outofspot = {}

    for row in tqdm(df_mge.iterrows(), desc = f"Associating each mge hit to a spot (cpu:{index_cpu})", leave=False):
        if df_SP[df_SP["all"].str.contains(row[1]["target name"])].empty == False:
            d_tmp_res[row[1]["target name"]] = {"gene_name": row[1]["target name"],
                                                "spot_number": df_SP.loc[df_SP[df_SP["all"].str.contains(row[1]["target name"])].index.tolist()[0],"interval_number"],
                                                "mge_hit":row[1]["query_name"],
                                                "E-value": row[1]["E-value"]
                                                }

        else : #need to determine whether the mge hit is part of BI or core genome
            if row[1]["target name"] in list_core_genes:
                d_tmp_outofspot[row[1]["target name"]] = {"gene_name": row[1]["target name"], "core_or_BI": "core", "mge_hit": row[1]["query_name"], "E-value": row[1]["E-value"]}
            else :
                d_tmp_outofspot[row[1]["target name"]] = {"gene_name": row[1]["target name"], "core_or_BI": "BI", "mge_hit": row[1]["query_name"], "E-value": row[1]["E-value"]}
    
    return [d_tmp_res,d_tmp_outofspot]


def max_genes_copy_in_a_single_genome(list2test):
    """
    Function which determine the max number of gene copy it exists in a single genome. 
    Will return this max number will return the list of genes considered as different copies 
    """
    max_copies = 0
    list_trimmed = [x.rsplit(".",1)[0] for x in list2test]
    for item in list_trimmed:
        if list_trimmed.count(item) > max_copies:
            max_copies = list_trimmed.count(item)
            genome_with_max = item
    
    list_index = []
    for index in range(len(list_trimmed)):
        if list_trimmed[index] == genome_with_max:
            list_index.append(index)
    
    list_max_copies = [list2test[i] for i in list_index]    
        
    return list_max_copies



def search_defense(outdir, n_cpu = 1):
    """
    Function which uses defense finder to search for defense systems in the given pangenome.
    Then uses MMseqs2 to determine the non redundant set of defense systems per spot.
    """
    
    #Running defense finder
    #uncomment here when testing is done
    subprocess.call(f"defense-finder run -o {outdir}/TAT_defense_results/TA_defense_search --db-type gembase {outdir}/TAT_defense_results/tmp/all_prt_formated.prt", shell = True)

    #Parsing the results + removing the redundant systems
    df_spot_pangenome = pd.read_csv(f"{outdir}/TAT_defense_results/core_mapping/2-all_genomes_intervals.tsv", sep = "\t", comment = "#")
    df_spot_pangenome.set_index("interval_number", inplace = True)
    df_spot_pangenome.dropna(axis = 0, how = "all", inplace = True)
    
    df_summary, df_spot_defense_without_duplicate, df_spot_defense_without_duplicate_numbers, df_systems_out_of_spotpangenome = defense_finder_parser(f"{outdir}/TAT_defense_results/TA_defense_search/all_prt_formated_defense_finder_systems.tsv", 
                                                                                                                                                      df_spot_pangenome, 
                                                                                                                                                      f"{outdir}/TAT_defense_results/tmp/all_prt_cat.prt", 
                                                                                                                                                      outdir)
    
    df_summary.to_csv(f"{outdir}/TAT_defense_results/TA_defense_search/2.1-Summary_defense_finder_systems.tsv", sep = "\t", index = False)
    df_spot_defense_without_duplicate.to_csv(f"{outdir}/TAT_defense_results/TA_defense_search/2.2-Spot_defense_finder_systems.tsv", sep = "\t", index = False)
    df_spot_defense_without_duplicate_numbers.to_csv(f"{outdir}/TAT_defense_results/TA_defense_search/2.3-Spot_defense_finder_systems_numbers.tsv", sep = "\t", index = False)
    df_systems_out_of_spotpangenome.to_csv(f"{outdir}/TAT_defense_results/TA_defense_search/2.4-Systems_out_of_spot_pangenome.tsv", sep = "\t", index = False)






def defense_finder_parser(defense_finder_systems_file, df_spot_pangenome, cat_prt_file, outdir):

    #will create two files => one summary with the numbers in each genome of each type of Defense system
    # second file will be the corresponding spot to each system.
    df_defense = pd.read_csv(defense_finder_systems_file, sep = "\t")
    
    #Note added a slight step to convert the names of the initial genes back to their real names , otherwise when assigning the systems to their spot, we might miss those that were initially in borders of contigs
    df_real_gene_names = pd.read_csv(f"{outdir}/TAT_defense_results/tmp/gene_names_formated.tsv", sep = "\t", comment = "#", names = ["gene_name_formated", "gene_name_real"])
    d_real_gene_names = dict(zip(df_real_gene_names["gene_name_formated"], df_real_gene_names["gene_name_real"]))
    list_sys_beg = [d_real_gene_names[x] for x in df_defense["sys_beg"].tolist()]
    list_sys_end = [d_real_gene_names[x] for x in df_defense["sys_end"].tolist()]
    list_protein_in_syst = [(",").join([d_real_gene_names[y] for y in x.split(",")]) for x in df_defense["protein_in_syst"].tolist()]
    
    df_defense["sys_beg"] = list_sys_beg
    df_defense["sys_end"] = list_sys_end
    df_defense["protein_in_syst"] = list_protein_in_syst

    # creating first file
    df_defense["genome"] = df_defense.apply(lambda x: x["sys_id"].rsplit(".", 1)[0], axis = 1)

    d_summary = {}
    for genome in df_spot_pangenome.columns.tolist():
        d_summary[genome] = {}
        d_summary[genome]["Total"] = len(df_defense[df_defense["genome"] == genome])
        for defense_type in df_defense["type"].drop_duplicates().tolist():
            d_summary[genome][defense_type] = len(df_defense[(df_defense["genome"] == genome) & (df_defense["type"] == defense_type)])

    df_summary = pd.DataFrame.from_dict(d_summary, orient = "index")
    df_summary.reset_index(inplace = True)
    df_summary.rename(columns={"index":"Genome"}, inplace= True)

    # now associating each defense system identified to a spot.
    # for that we detect which interval contains both first gene and last gene of the system
    # note will create aswell a txt file with the systems that are not totally part of the spot-pangenome

    d_systems_out_of_spotpangenome = {}
    d_spot_defense = {}
    for i in range(df_spot_pangenome.index[-1]+1):
        d_spot_defense[i]= {}
        for defense_type in df_defense["type"].drop_duplicates().tolist():
            d_spot_defense[i][defense_type] = []

    for row in tqdm(df_defense.iterrows(), desc = "Checking whether each out of spot pangenome defense hit are part of core genome or Breakpoint interval"):
        if df_spot_pangenome[(df_spot_pangenome[row[1]["genome"]].str.contains(f"{row[1]['sys_beg']}", na = False)) & (df_spot_pangenome[row[1]["genome"]].str.contains(f"{row[1]['sys_end']}", na = False))].empty == False:
            d_spot_defense[df_spot_pangenome[(df_spot_pangenome[row[1]["genome"]].str.contains(f"{row[1]['sys_beg']}", na = False)) & (df_spot_pangenome[row[1]["genome"]].str.contains(f"{row[1]['sys_end']}", na = False))].index[0]][row[1]["type"]].append(row[1]["protein_in_syst"])
        else :
            d_systems_out_of_spotpangenome[row[1]['sys_id']] = {"type": row[1]['type'], "protein_in_syst": row[1]['protein_in_syst']}

    df_spot_defense = pd.DataFrame.from_dict(d_spot_defense, orient = "index")
    df_spot_defense = df_spot_defense.loc[~(df_spot_defense == 0).all(axis=1)]
    df_spot_defense.reset_index(inplace = True)
    df_spot_defense.rename(columns={"index":"Spot_number"}, inplace= True)
    df_spot_defense_without_duplicate, df_spot_defense_without_duplicate_numbers = remove_redundants(df_spot_defense, cat_prt_file, outdir)

    list_total_number_per_spot =[]
    for spot in df_spot_defense_without_duplicate_numbers.Spot_number.tolist():
        list_total_number_per_spot.append(len(df_spot_defense_without_duplicate_numbers[df_spot_defense_without_duplicate_numbers["Spot_number"] == spot]))
    
    
    #Determine whether each out of spotpangenome systems are part of the core genome or in breakpoint interval or in a contig without core genes
    df_systems_out_of_spotpangenome = analyze_outofspot(d_systems_out_of_spotpangenome, outdir)


    return df_summary, df_spot_defense_without_duplicate, df_spot_defense_without_duplicate_numbers , df_systems_out_of_spotpangenome




def remove_redundants(df_spot_defense, cat_prt_file, outdir):

    """
    Function which remove redundants systems in each spot.
    Note that we consider two systems as redundant if at least one of their proteins is clustered through mmseqs
    """
    
    #first we create a tmp file for the mmseq alignment with the prt sequences of each systems in each spot
    prt_SeqIO = SeqIO.index(cat_prt_file, "fasta")
    d_prt2keep = {}
    for defense_type in df_spot_defense.loc[:,df_spot_defense.columns != "Spot_number"].columns.tolist():
        d_prt2keep[defense_type] = [x for x in df_spot_defense[defense_type].tolist() if x != []]
    
    #Now using mmseqs for each type of systems to determine if we have redundants systems in our dataset
    d_defense_representatives = {} #k = defense type, v = list of representatives systems
    for defense_type, list_of_prt_lists in d_prt2keep.items():
        d_defense_representatives[defense_type] = []
        for list_prt in list_of_prt_lists:
            str_tmp_prt = ""
            for i in list_prt:
                for x in i.strip().split(","):
                    str_tmp_prt += f">{x}\n{prt_SeqIO[x].seq}\n"
        
            if str_tmp_prt == "":
                continue
            else :
                if  os.path.exists(cat_prt_file+".tmp_spot_defense_cluster4mmseqs.faa") == True:
                    os.remove(cat_prt_file+".tmp_spot_defense_cluster4mmseqs.faa") # Idk why I have some problems sometimes here with the tmp file overwriting so added this line
                with open(cat_prt_file+".tmp_spot_defense_cluster4mmseqs.faa", "w") as f:
                    f.write(str_tmp_prt)

                subprocess.call(f"mmseqs easy-cluster {cat_prt_file+'.tmp_spot_defense_cluster4mmseqs.faa'} {outdir}/TAT_defense_results/tmp/TA_defense_search/spot_defense {outdir}/TAT_defense_results/tmp/TA_defense_search --min-seq-id 0.8 -v 0 --remove-tmp-files", shell = True)
                d_defense_representatives[defense_type] += parse_representatives_mmseqs(f"{outdir}/TAT_defense_results/tmp/TA_defense_search/spot_defense_cluster.tsv", list_prt)

    #now we localize again each representative system within the spot pangenome
    d_representatives_res = {}
    for i in range(len(df_spot_defense)):
        d_representatives_res[i] = {}
        for defense_type in d_defense_representatives.keys():
            d_representatives_res[i][defense_type] = []
        
    df_spot_defense.loc[:,df_spot_defense.columns != "Spot_number"] = df_spot_defense.loc[:,df_spot_defense.columns != "Spot_number"].applymap(lambda x: ",".join(x))
    for defense_type, list_str_representatives in d_defense_representatives.items():
        for str_representatives in list_str_representatives:
            d_representatives_res[df_spot_defense[df_spot_defense[defense_type].str.contains(str_representatives)]["Spot_number"].values[0]][defense_type] += [str_representatives]

    df_representatives_res = pd.DataFrame.from_dict(d_representatives_res, orient = "index")
    df_representatives_res = df_representatives_res.loc[~(df_representatives_res.astype(str) == "[]").all(axis=1)]
    df_representatives_res.reset_index(inplace = True)
    df_representatives_res.rename(columns={"index":"Spot_number"}, inplace= True)
    df_representatives_res_numbers = df_representatives_res.copy(deep = True)
    df_representatives_res_numbers.loc[:,df_representatives_res_numbers.columns != "Spot_number"] = df_representatives_res_numbers.loc[:,df_representatives_res_numbers.columns != "Spot_number"].applymap(lambda x: len(x))
    df_representatives_res.loc[:,df_representatives_res.columns != "Spot_number"] = df_representatives_res.loc[:,df_representatives_res.columns != "Spot_number"].applymap(lambda x: ",".join(x))

    return df_representatives_res, df_representatives_res_numbers



def parse_representatives_mmseqs(mmseqs_tsv, list_defense_this_spot):
    
    df_mmseqs = pd.read_csv(mmseqs_tsv, sep = "\t", names = ["query","subject"])
    list_defense_representatives = copy.deepcopy(list_defense_this_spot) #we will remove systems from this list in case there are redundance with others systems
    #here we considered as redundant, systems with at least one gene matching during the mmseqs analysis

    for defense_system in list_defense_this_spot:
        if defense_system in list_defense_representatives:
            for gene in defense_system.strip().split(","):
                if len(df_mmseqs[df_mmseqs["query"] == gene]) == 1:
                    continue
                else:
                    for row in df_mmseqs[(df_mmseqs["query"] == gene) & (df_mmseqs["subject"] != gene)].iterrows():
                        list_defense_representatives = remove_redundants_systems(row[1]["query"], row[1]["subject"], list_defense_representatives)

    return list_defense_representatives


def remove_redundants_systems(query, subject, list_defense_representatives):
    for i in list_defense_representatives:
        if subject in i and query not in i:
            list_defense_representatives.remove(i)

    return list_defense_representatives




def analyze_outofspot(d_outofspot, outdir):

    list_all_core_genes = []
    with open(f"{outdir}/TAT_defense_results/tmp/list_allcore_genes.lst", "r") as f:
        for line in f:
            list_all_core_genes += line.split(",")
    
    seqIO_allproteome = SeqIO.index(f"{outdir}/TAT_defense_results/tmp/all_prt_cat.prt", "fasta")

    d_sys_in_core = {}
    str_core_sys_proteome = ""
    d_sys_in_breakpoint_interval = {}
    str_BI_sys_proteome = ""
    for system_id, system_values in d_outofspot.items():
        is_in_core = False
        for gene_component in system_values["protein_in_syst"].split(","):
            if gene_component in list_all_core_genes:
                is_in_core = True
                break
        
        if is_in_core == True:
            d_sys_in_core[system_id] = {"type": system_values["type"], "core_or_BI": "core", "protein_in_syst":system_values["protein_in_syst"]}
            for gene in system_values["protein_in_syst"].split(","):
                str_core_sys_proteome += f">{gene}\n{seqIO_allproteome[gene].seq}\n"

        else:
            d_sys_in_breakpoint_interval[system_id] = {"type": system_values["type"], "core_or_BI": "BI", "protein_in_syst":system_values["protein_in_syst"]}
            for gene in system_values["protein_in_syst"].split(","):
                str_BI_sys_proteome += f">{gene}\n{seqIO_allproteome[gene].seq}\n"

    #running mmseqs for each system either in core or in BI to eliminate redundants systems
    if os.path.exists(f"{outdir}/TAT_defense_results/tmp/TA_defense_search/sys_in_core.tmp_cluster4mmseqs.faa") == True:
        os.remove(f"{outdir}/TAT_defense_results/tmp/TA_defense_search/sys_in_core.tmp_cluster4mmseqs.faa")
    if os.path.exists(f"{outdir}/TAT_defense_results/tmp/TA_defense_search/sys_in_BI.tmp_cluster4mmseqs.faa") == True:
        os.remove(f"{outdir}/TAT_defense_results/tmp/TA_defense_search/sys_in_BI.tmp_cluster4mmseqs.faa")
    
    with open(f"{outdir}/TAT_defense_results/tmp/TA_defense_search/sys_in_core.tmp_cluster4mmseqs.faa", "w") as f:
        f.write(str_core_sys_proteome)
    with open(f"{outdir}/TAT_defense_results/tmp/TA_defense_search/sys_in_BI.tmp_cluster4mmseqs.faa", "w") as f:
        f.write(str_BI_sys_proteome)
    
    
    if str_core_sys_proteome != "":    
        subprocess.call(f"mmseqs easy-cluster {outdir}/TAT_defense_results/tmp/TA_defense_search/sys_in_core.tmp_cluster4mmseqs.faa {outdir}/TAT_defense_results/tmp/TA_defense_search/sys_in_core {outdir}/TAT_defense_results/tmp/TA_defense_search --min-seq-id 0.8 -v 0 --remove-tmp-files", shell = True)
        df_mmseqs_sys_in_core = pd.read_csv(f"{outdir}/TAT_defense_results/tmp/TA_defense_search/sys_in_core_cluster.tsv", sep = "\t", names = ["query","subject"])
        
    if str_BI_sys_proteome != "":
        subprocess.call(f"mmseqs easy-cluster {outdir}/TAT_defense_results/tmp/TA_defense_search/sys_in_BI.tmp_cluster4mmseqs.faa {outdir}/TAT_defense_results/tmp/TA_defense_search/sys_in_BI {outdir}/TAT_defense_results/tmp/TA_defense_search --min-seq-id 0.8 -v 0 --remove-tmp-files", shell = True)
        df_mmseqs_sys_in_BI = pd.read_csv(f"{outdir}/TAT_defense_results/tmp/TA_defense_search/sys_in_BI_cluster.tsv", sep = "\t", names = ["query","subject"])
    

    #changing dict into df so it will ease next operations
    df_in_core = pd.DataFrame.from_dict(d_sys_in_core, orient = "index")
    df_in_BI = pd.DataFrame.from_dict(d_sys_in_breakpoint_interval, orient = "index")

    d_final_core_BI_sys = {}
    list_redundant_sys = []
    n_index_final = 0
    if df_in_core.empty == False:
        for row in df_in_core.iterrows():
            tmp_list_redundant= []
            if row[0] in list_redundant_sys:
                continue

            for gene_in_sys in row[1]["protein_in_syst"].split(","):
                if len(df_mmseqs_sys_in_core[df_mmseqs_sys_in_core["query"] == gene_in_sys]) == 1:
                    continue
                else :
                    for row_mmseqs in df_mmseqs_sys_in_core[(df_mmseqs_sys_in_core["query"] == gene_in_sys) & (df_mmseqs_sys_in_core["subject"] != gene_in_sys)].iterrows():
                        if df_in_core[df_in_core["protein_in_syst"].str.contains(row_mmseqs[1]["subject"])].index.tolist()[0] not in list_redundant_sys:
                            list_redundant_sys.append(df_in_core[df_in_core["protein_in_syst"].str.contains(row_mmseqs[1]["subject"])].index.tolist()[0])
                            tmp_list_redundant.append(df_in_core[df_in_core["protein_in_syst"].str.contains(row_mmseqs[1]["subject"])].index.tolist()[0])
                
            d_final_core_BI_sys[n_index_final] = {"system_id":row[0], "BI_or_core": "core", "type":row[1]["type"], "number_of_redundant": len(tmp_list_redundant)+1, "list_redundants":tmp_list_redundant}
            n_index_final += 1

    #same but for systems in breakpoint interval
    if df_in_BI.empty == False:
        for row in df_in_BI.iterrows():
            tmp_list_redundant= []
            if row[0] in list_redundant_sys:
                continue

            for gene_in_sys in row[1]["protein_in_syst"].split(","):
                if len(df_mmseqs_sys_in_BI[df_mmseqs_sys_in_BI["query"] == gene_in_sys]) == 1:
                    continue
                else :
                    for row_mmseqs in df_mmseqs_sys_in_BI[(df_mmseqs_sys_in_BI["query"] == gene_in_sys) & (df_mmseqs_sys_in_BI["subject"] != gene_in_sys)].iterrows():
                        if df_in_BI[df_in_BI["protein_in_syst"].str.contains(row_mmseqs[1]["subject"])].index.tolist()[0] not in list_redundant_sys:
                            list_redundant_sys.append(df_in_BI[df_in_BI["protein_in_syst"].str.contains(row_mmseqs[1]["subject"])].index.tolist()[0])
                            tmp_list_redundant.append(df_in_BI[df_in_BI["protein_in_syst"].str.contains(row_mmseqs[1]["subject"])].index.tolist()[0])
                
            d_final_core_BI_sys[n_index_final] = {"system_id":row[0], "BI_or_core": "BI", "type":row[1]["type"], "number_of_redundant": len(tmp_list_redundant)+1, "list_redundants":tmp_list_redundant}
            n_index_final += 1

    if df_in_core.empty == False or df_in_BI.empty == False:
        df_final_core_BI_sys = pd.DataFrame.from_dict(d_final_core_BI_sys, orient = "index")
    else :
        df_final_core_BI_sys = pd.DataFrame()
    
    return df_final_core_BI_sys




def parse_spot_results(outdir, list_to_search):
    
    """
    Function which for each spot, determine the number of each elements in it.
    Note that MazEF systems are now also considered as anti-phage systems. Due to their involvement in the regulation of the SP-beta phage cycle in bacillus.
    However, here they are classified as TAs, and those that were found by defense finder are not considered.
    
    """
    
    list_df_to_merge_detailled = [] #to store the results with all families for TAs and defense
    list_df_to_merge_only_total = [] #to store only the results of the total number of TAs and defense in each spot
    
    for searched in list_to_search :
        
        if searched == "TA":
            df_TA = pd.read_csv(outdir + "/TAT_defense_results/TA_defense_search/1.2-TAlist_per_spot_NR.tsv",
                                sep = "\t")
            df_TA = df_TA.drop(["TAlist", "TA_families"], axis = 1)
            list_df_to_merge_detailled.append(df_TA)
            list_df_to_merge_only_total.append(df_TA.iloc[:, [0, -1]])
            
        
        elif searched == "defense":
            df_defense = pd.read_csv(outdir + "/TAT_defense_results/TA_defense_search/2.3-Spot_defense_finder_systems_numbers.tsv",
                                    sep = "\t")

            #removing the MazEF results from the defense search
            if "MazEF" in df_defense.columns:
                df_defense = df_defense.drop(["MazEF"], axis = 1)
            
            cols_to_sum = [col for col in df_defense.columns if col != "Spot_number"]
            df_defense["n_defense"] = df_defense[cols_to_sum].sum(axis=1)
            list_df_to_merge_detailled.append(df_defense)
            list_df_to_merge_only_total.append(df_defense.iloc[:, [0, -1]])


        elif searched == "IS":
            df_IS = pd.read_csv(outdir + "/TAT_defense_results/TA_defense_search/3.3-IS_per_spot_pangenome.tsv",
                                sep = "\t")
            list_df_to_merge_detailled.append(df_IS)
            #Here the next line should not change much the df_IS, but might be useful if we want to separate each IS families later on.
            list_df_to_merge_only_total.append(df_IS.iloc[:, [0, -1]])
            
            
        elif searched == "integrase":
            df_integrase = pd.read_csv(outdir + "/TAT_defense_results/TA_defense_search/4.3-Integrase_per_spot_pangenome.tsv",
                                sep = "\t")
            list_df_to_merge_detailled.append(df_integrase)
            #Same remark than for IS
            list_df_to_merge_only_total.append(df_integrase.iloc[:, [0, -1]])
    
    
    df_all_results_detailled = list_df_to_merge_detailled[0]
    if len(list_df_to_merge_detailled) > 1 :
        for df2merge in list_df_to_merge_detailled[1:]:
            df_all_results_detailled = df_all_results_detailled.merge(df2merge,
                                 on = "Spot_number",
                                 how = "outer")
    
    
    df_all_results_only_total = list_df_to_merge_only_total[0]
    if len(list_df_to_merge_only_total) > 1 :
        for df2merge in list_df_to_merge_only_total[1:]:
            df_all_results_only_total = df_all_results_only_total.merge(df2merge,
                                 on = "Spot_number",
                                 how = "outer")
    
    df_all_results_only_total.to_csv(outdir + "/TAT_defense_results/TA_defense_search/5.0-Spot_merged_results_only_total.tsv", sep = "\t", index = False)
    df_all_results_detailled.to_csv(outdir + "/TAT_defense_results/TA_defense_search/5.1-Spot_merged_results_detailled.tsv", sep = "\t", index = False)

    return df_all_results_only_total



def compute_colocalization(outdir, df_final_results):
    """
    Function which determines the tendency of each element type to colocalize with others (i.e. TAs, defenses, IS, integrases).
    Note that here we only considers
    """
    
    list_columns = [x for x in df_final_results.columns if x != "Spot_number"]
    d_coloc = {}
    
    for col in list_columns:
        
        sum_total_col1 = df_final_results[col].sum()
        sum_col1_self_coloc = df_final_results[df_final_results[col] >= 2][col].sum()
        d_coloc[col.split("_", 1)[-1]] = {"N_total": sum_total_col1,
                                          "N_self_colocalization": sum_col1_self_coloc,
                                          "%_self_colocalization": sum_col1_self_coloc/sum_total_col1*100,
                                          "list_other_elements": []
                                          }
        
        #Now determining the colocalization with the others type of elements
        for col2 in list_columns:
            if col2 != col:
                sum_coloc_with_col2_element = df_final_results[df_final_results[col2] > 0][col].sum()
                d_coloc[col.split("_", 1)[-1]]["N_colocalization_with_" + col2.split("_", 1)[-1]] = sum_coloc_with_col2_element
                d_coloc[col.split("_", 1)[-1]]["%_colocalization_with_" + col2.split("_", 1)[-1]] = sum_coloc_with_col2_element/sum_total_col1*100
                d_coloc[col.split("_", 1)[-1]]["list_other_elements"].append(col2.split("_", 1)[-1])
    
    # Now writing the results in the outfile
    Str2write_res = ("Results of the colocalization analysis considering the following element types:\n" 
                    f"{('/').join([x.split('_', 1)[-1] for x in list_columns])}"
                    "\n\n"
                    )

    for k, v in d_coloc.items():
        Str2write_res += (f"####    Colocalization results from {k} perspectives :  ####\n"
                          f"N_total of {k} elements :    {int(v['N_total'])}\n"
                          f"N_self_colocalization :    {int(v['N_self_colocalization'])}\n"
                          f"%_self_colocalization :    {float(v['%_self_colocalization'])}\n"
                          )
        
        for other_element_type in v["list_other_elements"]:
            Str2write_res += (f"--> Colocalization of {k} with {other_element_type} :\n"
                              f"N_colocalization :    {int(v['N_colocalization_with_' + other_element_type])}\n"
                              f"%_colocalization :    {int(v['%_colocalization_with_' + other_element_type])}\n"
                              )
        
        Str2write_res += "\n\n"
    
    
    with open(outdir + "/TAT_defense_results/TA_defense_search/5.2-Colocalization_results.txt", "w") as f:
        f.write(Str2write_res)