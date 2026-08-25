import pandas as pd
import numpy as np
import os
from tqdm import tqdm
import multiprocessing as mp
from Bio import SeqIO
import subprocess

def main(indir, core_file, outdir, ref_genome = "First",  n_cpu = 1):
    """
    main function to establish the core genome mapping and the spot-pangenome within the provided files.
    Note that it requires as the indir the directoy containing the two folders Proteins/gff3 from PanACoTA annotate + the core genome file from PanACoTA corepers
    """
    
    #Checking whether the provided ref_genome to establish the mapping is within the gff3 folder, if none are provided use the first (alphabetically) file
    gff_indir = f"{indir}/gff3"
    list_genomes = [x.rsplit(".",1)[0] for x in os.listdir(gff_indir)]
    
    if ref_genome == "First":
        ref_genome = list_genomes[0]
        
    else:
        if ref_genome not in list_genomes:
            print(f"Could not find the {ref_genome} gff file in the gff3 folder, please make you sure the associated file is within the gff3 folder")
            exit()
            
        else :
            ref_genome = ref_genome

    #Then using the PanACoTA corepers file to create a tsv file with all core gene coordinates
    print("Fetching and associating the features to each core gene") #note: might do a logger later
    if os.path.exists(f"{outdir}/TAT_defense_results/core_mapping/core_genome_features.tsv"):
        df_core = pd.read_csv(f"{outdir}/TAT_defense_results/core_mapping/core_genome_features.tsv", sep = "\t")
    else :    
        df_core = corefile2tsv(core_file, gff_indir, outdir)

    #Using the newly created tsv file to create the reference interval file
    print(f"Establishing the reference intervals based on the reference genome : {ref_genome}")
    if os.path.exists(f"{outdir}/TAT_defense_results/core_mapping/0-Ref_genome_intervals.tsv"):
        df_ref_intervals = pd.read_csv(f"{outdir}/TAT_defense_results/core_mapping/0-Ref_genome_intervals.tsv", sep = "\t")
    else :
        df_ref_intervals = ref_interval(ref_genome, df_core, outdir)

    #Now based on the reference intervals, check whether or not we can recreate them in the other genomes
    df_gff = pd.read_csv(outdir + "/TAT_defense_results/tmp/all_gff_cat.gff", comment = "#", sep = "\t")
    
    print(f"Checking for each interval and in each genome whether they can be recreated or should be considered as Breakpoint Interval")
    
    if os.path.exists(f"{outdir}/TAT_defense_results/core_mapping/1-all_genomes_check_corespot.tsv"):
        df_check = pd.read_csv(f"{outdir}/TAT_defense_results/core_mapping/1-all_genomes_check_corespot.tsv", sep = "\t")
    else :
        df_check = check_interval(ref_genome, df_core, df_ref_intervals, df_gff, list_genomes, outdir)
    
    #For those accepted, associate each intervals to their accessory genes
    print(f"Fetching and associating the accessory genes to their intervals (except for breakpoint interval)")    
    if os.path.exists(f"{outdir}/TAT_defense_results/core_mapping/2-all_genomes_intervals.tsv"):
        df_all_intervals = pd.read_csv(f"{outdir}/TAT_defense_results/core_mapping/2-all_genomes_intervals.tsv", sep = "\t")
    else :
        df_all_intervals = recreate_intervals(df_core, df_gff, df_check, df_ref_intervals, list_genomes, outdir)
    
    #Finally, we use mmseqs2 to obtain all necessary informations for each genomic spot
    print("Using mmseqs2 to determine how many accessory gene families are present in each genomic spot")
    get_spot_information(df_all_intervals, outdir)
        
        
        
        

def corefile2tsv(core_file, gff_indir, outdir):
    
    """
    Function that create a tsv file with all informations associated to each genes (e.g. in which genome, which core family gene it is, ...)
    """
    
    #first we create a dictionnary with store for each core gene (from every genome)
    d_core = corepers2df(core_file)

    #now for every gene in this dict we open the genome file and get the features of the core gene
    list_core = []
    list_genome = []
    for k,v in d_core.items():
        list_core.append(k)
        if v["genome_name"] not in list_genome:
            list_genome.append(v["genome_name"])
    
    gff_columns = ["contig","source","type", "left_coordinate", "right_coordinate", "dot", "strand", "0" ,"description"]

    for g in list_genome:
        df_tmp = pd.read_csv(f"{gff_indir}/{g}.gff", sep = "\t", comment = "#", names = gff_columns)
        df_tmp["description"] = df_tmp["description"].apply(lambda x: x.split(';')[0].split('=')[-1])

        for gene, d_gene in d_core.items():       
            if d_gene["genome_name"] == g:
                d_core[gene]["left_coordinate"] = df_tmp[df_tmp["description"] == gene]["left_coordinate"].item()
                d_core[gene]["right_coordinate"] = df_tmp[df_tmp["description"] == gene]["right_coordinate"].item()
                d_core[gene]["strand"] = df_tmp[df_tmp["description"] == gene]["strand"].item()

    df_core = pd.DataFrame.from_dict(d_core, orient = "index")
    df_core.to_csv(f"{outdir}/TAT_defense_results/core_mapping/core_genome_features.tsv", sep = "\t", index = False)
    
    return df_core



def corepers2df(file2df, exclude = None):

    #return a dataframe with index = gene name, columns: gene name/ comparative genomic family / number of genes in this family /Genome name / number of the ORF gene
    d={}
    with open(file2df, "r") as f:
        for line in f:
            for i in range(1, len(line.split(" "))):
                if exclude == None :
                    d[line.split(" ")[-i].strip()] = {"gene_name": line.split(" ")[-i].strip(), "core_family": line.split(" ")[0], "genome_name": line.split(" ")[-i].rsplit("_",1)[0][:-6],
                                                    "contig": line.split(" ")[-i].rsplit("_",1)[0][:-1], "gene_number": int(line.split(" ")[-i].rsplit("_",1)[1]) }

                else :
                    if line.split(" ")[-i].rsplit("_",1)[0][:-6] != exclude:
                        d[line.split(" ")[-i].strip()] = {"gene_name": line.split(" ")[-i].strip(), "core_family": line.split(" ")[0], "genome_name": line.split(" ")[-i].rsplit("_",1)[0][:-6],
                                                        "contig": line.split(" ")[-i].rsplit("_",1)[0][:-1], "gene_number": int(line.split(" ")[-i].rsplit("_",1)[1]) }

    return d



def ref_interval(ref_genome, df_core, outdir):
    
    #First we divide the reference genome in intervals based on their core genome 
    df_core_ref = df_core[df_core["genome_name"] == ref_genome].sort_values(by = "left_coordinate").reset_index(drop= True)
    #creating a file with all core genes separated by a comma useful for control later
    str_all_core_genes = ""
    tmp_index = 0
    for i in df_core["gene_name"].tolist():
        if tmp_index > 0:
            str_all_core_genes += ","
        str_all_core_genes += i
        tmp_index += 1
    with open(f"{outdir}/TAT_defense_results/tmp/list_allcore_genes.lst", "w") as f:
        f.write(str_all_core_genes)


    d_intervals = {}
    for row_index in df_core_ref.index.tolist():
        if row_index + 1 < len(df_core_ref):
            d_intervals[row_index] = {  "interval_number": row_index,
                                        "left_core": df_core_ref.loc[row_index, "gene_name"],
                                        "left_core_family": df_core_ref.loc[row_index, "core_family"],
                                        "right_core": df_core_ref.loc[row_index+1, "gene_name"],
                                        "right_core_family": df_core_ref.loc[row_index+1, "core_family"],
                                        "left_c": df_core_ref.loc[row_index, "right_coordinate"] ,
                                        "right_c": df_core_ref.loc[row_index+1, "left_coordinate"]
                                        }

        elif row_index + 1 == len(df_core_ref):
            d_intervals[row_index] = {  "interval_number": row_index,
                                        "left_core": df_core_ref.loc[row_index, "gene_name"],
                                        "left_core_family": df_core_ref.loc[row_index, "core_family"],
                                        "right_core": df_core_ref.loc[0, "gene_name"],
                                        "right_core_family": df_core_ref.loc[0, "core_family"],
                                        "left_c": df_core_ref.loc[row_index, "right_coordinate"] ,
                                        "right_c": df_core_ref.loc[0, "left_coordinate"]
                                        }
    
    df_intervals = pd.DataFrame.from_dict(d_intervals, orient = "index")

    df_intervals.to_csv(f"{outdir}/TAT_defense_results/core_mapping/0-Ref_genome_intervals.tsv", sep = "\t", index = False)
    
    return df_intervals


def check_interval(ref_genome, df_core, df_ref, df_gff, list_genomes, outdir):
    # Now based on the intervals generated with the reference genome, we try to rebuild the same interval in the others genomes
    # For that we need to check two things to accept the interval
    # First: in case the genome isn't fully assembled, both core genes used to delimit the interval need to be located on the same contig
    # Second : To make sure there are no important genomic rearrangement, for that we check if there are no others core gene within the same interval in this genome

    #it will return a dataframe with row = intervals , columns = genomes ====> each cell with either : accepted (we're able to rebuild the interval in this genome),
    # rejected (differents contigs) or rejected (rearrangments)
    # last columns count the number of accepted genomes in which we can rebuild the interval

    d_accepted = {}

    for interval_number in tqdm(df_ref.index.tolist(), desc = "Checking each interval in each genome"):
        d_accepted[interval_number] = {}
        for genome in list_genomes:
            if genome == ref_genome:
                d_accepted[interval_number][genome] = "Accepted"
            else :
                core1, core2 = get_core(df_ref.loc[interval_number, "left_core_family"], df_ref.loc[interval_number, "right_core_family"], genome, df_core)
                check_contig_results = check_contig(core1, core2, df_core, genome)
                
                if check_contig_results[0] == False:
                    if check_contig_results[1] == False:
                        d_accepted[interval_number][genome] = "rejected (differents contigs)"

                    elif check_contig_results[1] == True:
                        d_accepted[interval_number][genome] = "Accepted (differents contigs)"

                else :
                    if check_rearrangement(core1, core2, df_gff, df_core, genome) == True:
                        d_accepted[interval_number][genome] = "Accepted"
                    else :
                        d_accepted[interval_number][genome] = "rejected (rearrangments)"


    df_accepted = pd.DataFrame.from_dict(d_accepted, orient = "index")
    list_columns = df_accepted.columns.tolist()
    df_accepted["n_accepted"] = df_accepted.apply(lambda x: count_accepted(x, list_columns), axis = 1)
    df_accepted = df_accepted.reset_index().rename(columns = {"index": "interval_number"})
    df_accepted.to_csv(f"{outdir}/TAT_defense_results/core_mapping/1-all_genomes_check_corespot.tsv", sep = "\t", index = False)
    
    return df_accepted


def get_core(core_family1, core_family2, genome, df_core):

    df_core2search = df_core[df_core["genome_name"] == genome]
    new_core1 = df_core2search[df_core2search["core_family"] == core_family1]["gene_name"].values[0]
    new_core2 = df_core2search[df_core2search["core_family"] == core_family2]["gene_name"].values[0]

    return new_core1, new_core2


def check_contig(core1, core2, df_core, genome) :
    #return a tuple
    #First element is True of False depending on whether or not both core are located within the same contig
    #If both core genes are located on differents contigs, will also check whether the two cores are located at the extremities
    #We're making the assumption here that if both core are located on the extremities it is likely that there was no rearrangement 

    if core1.rsplit("_", 1)[0][:-1] == core2.rsplit("_", 1)[0][:-1]:
        return (True, True)
    else :
        df_tmp_core = df_core[df_core["genome_name"] == genome].sort_values(by = "left_coordinate").reset_index(drop = True)
        
        if is_border(core1, df_tmp_core) == True and is_border(core2, df_tmp_core) == True:
            return (False, True)
        else :
            return (False, False)


def is_border(core_gene, df_tmp_core):

    tmp_contig = df_tmp_core[df_tmp_core["gene_name"] == core_gene]["contig"].tolist()[0]
    df_tmp_core = df_tmp_core[df_tmp_core["contig"] == tmp_contig].sort_values(by = "left_coordinate").reset_index(drop = True)
    
    #In case the core genes is the only one on the contig, we cannot determine whether we should keep the upstream or dw genes for the interval
    #So we exclude this case
    if len(df_tmp_core) == 1 :
        return False

    if df_tmp_core.loc[0, "gene_name"] == core_gene or df_tmp_core.loc[len(df_tmp_core)-1, "gene_name"] == core_gene:
        return True 
    else:
        return False


def check_rearrangement(core1, core2, df_gff, df_core, genome):

    # get the list of genes within the interval (core excluded)
    df_tmp = df_gff[df_gff["genome"] == genome].reset_index(drop = True)
    list_core = df_core[df_core["genome_name"] == genome]["gene_name"].tolist()
    df_core = df_core[df_core["genome_name"] == genome].sort_values(by = "left_coordinate").reset_index(drop = True)
    list_core.remove(core1)
    list_core.remove(core2)

    index_core1 = df_tmp[df_tmp["id"] == core1].index[0]
    index_core2 = df_tmp[df_tmp["id"] == core2].index[0]
    
    #In case the genome is fully assembled and circular => check if both core are first/last core gene in the genome 
    if len(df_tmp["contig"].drop_duplicates().tolist()) == 1 and ((df_core[df_core["gene_name"] == core1].index[0] == 0 and df_core[df_core["gene_name"] == core2].index[0] == len(df_core)-1) or (df_core[df_core["gene_name"] == core2].index[0] == 0 and df_core[df_core["gene_name"] == core1].index[0] == len(df_core)-1)):
        list_genes_interval = []
        if min(index_core1,index_core2) > 0:
            list_genes_interval += df_tmp.iloc[:min(index_core1,index_core2)]["id"].tolist()
        if max(index_core1,index_core2)+1 < len(df_tmp):
            list_genes_interval += df_tmp.iloc[max(index_core1, index_core2)+1:]["id"].tolist()

    else:
        list_genes_interval = df_tmp.iloc[min(index_core1, index_core2)+1:max(index_core1, index_core2)]["id"].tolist()

    #now checking if there are others core genes within the list_genes_interval
    for gene in list_genes_interval:
        if gene in list_core:
            return False
    
    return True

def count_accepted(row, list_columns):
    count = 0
    for columns in list_columns:
        if row[columns] == "Accepted" or row[columns] == "Accepted (differents contigs)":
            count += 1

    return count





def recreate_intervals(df_core, df_gff, df_check, df_ref_intervals, list_genomes, outdir):
    """
    Function which fetch all accessories genes between each core genes pairs. 
    Note that core gene loners (only core gene within a single contig) are considered as Breakpoint intervals and excluded
    """

    df_gff = df_gff[df_gff["type"] == "CDS"].reset_index(drop = True)

    d_intervals = {}
    for interval_number in tqdm(df_ref_intervals["interval_number"].tolist(), desc = "Recreating all intervals in all genomes (if validated)"):
        d_intervals[interval_number] = {}
        for genome in list_genomes:
            if df_check.loc[interval_number,genome] == "Accepted" :
                core1, core2 = get_core(df_ref_intervals.loc[interval_number, "left_core_family"], df_ref_intervals.loc[interval_number, "right_core_family"], genome, df_core)
                d_intervals[interval_number][genome] = get_interval_genes(core1, core2, df_gff, df_core, genome)
            
            elif df_check.loc[interval_number,genome] == "Accepted (differents contigs)":
                core1, core2 = get_core(df_ref_intervals.loc[interval_number, "left_core_family"], df_ref_intervals.loc[interval_number, "right_core_family"], genome, df_core)
                d_intervals[interval_number][genome] = get_interval_genes_diff_contigs(core1, core2, df_gff, df_core, genome)

            else :
                d_intervals[interval_number][genome] = "Rejected"
    
    df_intervals = pd.DataFrame.from_dict(d_intervals, orient = "index")
    df_intervals = df_intervals.reset_index().rename(columns = {"index": "interval_number"})
    df_intervals.to_csv(f"{outdir}/TAT_defense_results/core_mapping/2-all_genomes_intervals.tsv", sep = "\t", index = False)
    
    return df_intervals


def get_interval_genes(core1, core2, df_gff, df_core, genome):

    df_tmp = df_gff[df_gff["genome"] == genome].reset_index(drop = True)
    index_core1 = df_tmp[df_tmp["id"] == core1].index[0]
    index_core2 = df_tmp[df_tmp["id"] == core2].index[0]
    df_core = df_core[df_core["genome_name"] == genome].sort_values(by = "left_coordinate").reset_index(drop = True)

    if len(df_tmp["contig"].drop_duplicates().tolist()) == 1 and ((df_core[df_core["gene_name"] == core1].index[0] == 0 and df_core[df_core["gene_name"] == core2].index[0] == len(df_core)-1) or (df_core[df_core["gene_name"] == core2].index[0] == 0 and df_core[df_core["gene_name"] == core1].index[0] == len(df_core)-1)):
        list_genes_interval = []
        if min(index_core1,index_core2) > 0:
            list_genes_interval += df_tmp.iloc[:min(index_core1,index_core2)]["id"].tolist()
        if max(index_core1,index_core2)+1 < len(df_tmp):
            list_genes_interval += df_tmp.iloc[max(index_core1, index_core2)+1:]["id"].tolist()

        return (",").join(list_genes_interval)

    else :
        return (",").join(df_tmp.iloc[min(index_core1, index_core2)+1:max(index_core1, index_core2)]["id"].tolist())
    


def get_interval_genes_diff_contigs(core1, core2, df_gff, df_core, genome):

    df_tmp = df_gff[df_gff["genome"] == genome].reset_index(drop = True)
    list_genes_interval = []
    df_core = df_core[df_core["genome_name"] == genome].sort_values(by = "left_coordinate").reset_index(drop = True)

    for core_gene in [core1,core2]:
        tmp_contig = df_core[df_core["gene_name"] == core_gene]["contig"].tolist()[0]
        df_tmp_core = df_core[df_core["contig"] == tmp_contig].reset_index(drop = True)
        df_tmp = df_gff[df_gff["contig"] == tmp_contig].sort_values(by = "left_c").reset_index(drop = True)
        tmp_index_core = df_tmp[df_tmp["id"] == core_gene].index[0]

        if df_tmp_core.loc[0,"gene_name"] == core_gene:
            list_genes_interval += df_tmp.loc[:tmp_index_core,"id"].tolist()

        elif df_tmp_core.loc[len(df_tmp_core)-1,"gene_name"] == core_gene:
            list_genes_interval += df_tmp.loc[tmp_index_core:,"id"].tolist()

        list_genes_interval.remove(core_gene)
    
    if len(list_genes_interval) == 0:
        return ""
    else:
        return (",").join(list_genes_interval)



def get_spot_information(df_all_intervals, outdir, n_cpu = 1):
    """
    Function that will use mmseqs2 to determine how many different gene families are present in each spot.
    Also record other information like how many non-empty intervals there are for each spot, the maximal number of genes within a single interval of the spot
    Nb: thought I could use multiprocessing and mmseqs at the same time to speed up the operation, but it leads to segmentation fault errors.
    Therefore we should only use n_cpu = 1
    """
    
    #preparing the necessary variables for parallelization
    list_df_all_intervals_chunks = np.array_split(df_all_intervals, n_cpu)
    
    cpu_number_to_display = 0
    list_args_chunk = []
    for x in list_df_all_intervals_chunks:
        list_args_chunk.append([x, outdir, cpu_number_to_display])
        cpu_number_to_display += 1
    
    with mp.Pool(processes = n_cpu) as pool:
        results = pool.map(analyze_spot_chunk, list_args_chunk)

    d_res_spot = {}
    for result_chunk in results:
        d_res_spot.update(result_chunk)
    
    df_res_spot = pd.DataFrame.from_dict(d_res_spot, orient="index")
    df_res_spot.to_csv(f"{outdir}/TAT_defense_results/core_mapping/3-all_spots_features.tsv", sep = "\t", index = False)




def analyze_spot_chunk(args):
    """
    Parallelized function that will be used to analyze each spot.
    Args should be a list of arguments that will be passed to the function.
    As follow : [df_all_intervals_chunk, outdir, cpu_number]
    cpu_number, just so we can use tqdm in parallel
    """
    
    df_to_analyze = args[0]
    outdir = args[1]
    cpu_number = args[2]
    list_genomes = df_to_analyze.columns.tolist()[1:]
    list_intervals = df_to_analyze["interval_number"].tolist()
    df_to_analyze = df_to_analyze.set_index("interval_number")
    
    with open(f"{outdir}/TAT_defense_results/tmp/all_prt_cat.prt", "r") as handle:
        All_seqs_SeqIO = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))

    d_chunk_res = {}
    
    #Now for every spot, we fetch all accessory gene sequences and paste them in a tmp files
    
    for interval in tqdm(list_intervals, desc = "Clustering the accessory genes in each spot based on their sequence identities", position = cpu_number, total = len(list_intervals)):
        
        # Initializing the necessary variables
        n_genomes_non_empty_intervals = 0
        n_genomes_rejected = 0
        list_n_genes_per_genome = []
        tmp_genes_list = []
        
        for genome in list_genomes:
            if isinstance(df_to_analyze.loc[interval,genome], str) and df_to_analyze.loc[interval,genome] != "Rejected":
                if df_to_analyze.loc[interval,genome].split(",") != [" "]:
                    tmp_genes_list += df_to_analyze.loc[interval,genome].split(",")
                    list_n_genes_per_genome.append(len(df_to_analyze.loc[interval,genome].split(",")))
                    n_genomes_non_empty_intervals += 1

            elif isinstance(df_to_analyze.loc[interval,genome], str) and df_to_analyze.loc[interval,genome] == "Rejected":
                n_genomes_rejected += 1

        #Now loading the SeqIO dict to write only these genes within a tmp file
        
        if len(tmp_genes_list) == 0:
            continue
        
        elif len(tmp_genes_list) == 1 :
            d_chunk_res[interval] = {"interval_number": interval,
                                     "N_accessory_genes": 1,
                                     "N_families": 1,
                                     "N_NR_genes": 1,
                                     "Mean": 1,
                                     "N_genomes_non_empty_intervals": n_genomes_non_empty_intervals,
                                     "N_genomes_rejected": n_genomes_rejected
                                     }
        
        else :
            #If there are more than 1 genes, we write their sequences in a tmp file and run mmseqs2 to cluster them
            tmp_str_seq = ""
            for gene in tmp_genes_list:
                try:
                    tmp_str_seq += f">{gene}\n{All_seqs_SeqIO[gene].seq}\n"
                except KeyError:
                    print(f"Gene {tmp_str_seq} not found in the fasta file")
                    continue
            
            with open(f"{outdir}/TAT_defense_results/tmp/core_mapping/{interval}_seqs_to_cluster.fasta", "w") as f:
                f.write(tmp_str_seq)

            #Now we run mmseqs2 to cluster these sequences
            subprocess.run([
                            "mmseqs", "easy-cluster",
                            f"{outdir}/TAT_defense_results/tmp/core_mapping/{interval}_seqs_to_cluster.fasta",
                            f"{outdir}/TAT_defense_results/tmp/core_mapping/{interval}",
                            f"{outdir}/TAT_defense_results/tmp/core_mapping/trash",
                            "--min-seq-id", "0.8",
                            "-c", "0.8",
                            "--cov-mode", "0",
                            "--threads", "1",
                            "-v", "0"
                            ], check = True)
                
            #collecting data from the mmseqs2 output
            df_mmseqs2_tmp = pd.read_csv(f"{outdir}/TAT_defense_results/tmp/core_mapping/{interval}_cluster.tsv", sep = "\t", names = ["representative","members"])
            tmp_N_NR_genes = count_NR_genes(df_mmseqs2_tmp)
            
            d_chunk_res[interval] = {"interval_number": interval,
                                     "N_accessory_genes": len(tmp_genes_list),
                                     "N_families": len(df_mmseqs2_tmp["representative"].unique()),
                                     "N_NR_genes": tmp_N_NR_genes,
                                     "Mean": len(tmp_genes_list)/n_genomes_non_empty_intervals,
                                     "N_genomes_non_empty_intervals": n_genomes_non_empty_intervals,
                                     "N_genomes_rejected": n_genomes_rejected
                                     }


    return d_chunk_res



def count_NR_genes(df_mmseqs2_tmp):
    """
    Small function that counts non-redundant genes from the mmseqs2 output.
    Non-redundant genes consider the number of families but also whether there are multiples copies of a same families within a same genome.
    For example, if there are 7 families within a spot but that there are 2 copies of the same family within a same genome, then the number of non-redundant genes is 8.
    Note that in case a same family is present in different numbers when comparing the different genomes, only the highest number is considered.
    """
    
    d_families_NR = {}
    for family_representative in df_mmseqs2_tmp["representative"].drop_duplicates().tolist():
        d_tmp = {} #k = genome, v = number of occurence of the gene family
        for member in df_mmseqs2_tmp[df_mmseqs2_tmp["representative"] == family_representative]["members"].tolist():
            if member.rsplit(".",1)[0] not in d_tmp.keys():
                d_tmp[member.rsplit(".",1)[0]] = 1
            else:
                d_tmp[member.rsplit(".",1)[0]] += 1
    
        tmp_max_NR = 0
        for genome in d_tmp.keys():
            if d_tmp[genome] > tmp_max_NR:
                tmp_max_NR = d_tmp[genome]
                
        #Getting the max number of occurence of this family in a single genome
        d_families_NR[family_representative] = tmp_max_NR
                
    #Now counting the total number of non-redundant genes
    total_NR_genes = 0
    for representative_family in d_families_NR.keys():
        total_NR_genes += d_families_NR[representative_family]
    
    return total_NR_genes