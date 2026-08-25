import shutil
import os
import subprocess

def create_dir(outdir,**kwargs):
    
    to_run = kwargs.get("to_run", "all")
    remove_tmp = kwargs.get("remove_tmp", True)
    
    #in case the tmp folder already exist, remove it to avoid analysis problems
    if remove_tmp == True:
        shutil.rmtree(f"{outdir}/tmp", ignore_errors= True)

    list_dir = ["TAT_defense_results","TAT_defense_results/tmp"]
    
    if to_run == "all":
        list_dir.append("TAT_defense_results/core_mapping")
        list_dir.append("TAT_defense_results/tmp/core_mapping")
        list_dir.append("TAT_defense_results/TA_defense_search")
        list_dir.append("TAT_defense_results/tmp/TA_defense_search")
        list_dir.append("TAT_defense_results/simulation")
        list_dir.append("TAT_defense_results/simulation/per_family")
        list_dir.append("TAT_defense_results/simulation/plot")
        list_dir.append("TAT_defense_results/simulation/data")
        list_dir.append("TAT_defense_results/tmp/simulation")
        
    
    elif to_run == "core_mapping":
        list_dir.append(f"TAT_defense_results/{to_run}")
        list_dir.append(f"TAT_defense_results/tmp/{to_run}")

    elif to_run == "search":
        list_dir.append(f"TAT_defense_results/TA_defense_search")
        list_dir.append(f"TAT_defense_results/tmp/TA_defense_search")
        
    elif to_run == "simulation":
        list_dir.append(f"TAT_defense_results/{to_run}")
        list_dir.append(f"TAT_defense_results/{to_run}/plot")
        list_dir.append(f"TAT_defense_results/{to_run}/data")
        list_dir.append(f"TAT_defense_results/tmp/{to_run}")
        

    list_dir.insert(0,outdir)


    for i in list_dir:
        if i == list_dir[0]:
            try :
                os.mkdir(i)
            except FileExistsError:
                continue
        else:
            try :
                os.mkdir(f"{outdir}/{i}")
            except FileExistsError:
                continue


def threads2use(n):
    
    from concurrent.futures import ThreadPoolExecutor
    executor = ThreadPoolExecutor()

    if int(n) != 0 and int(n) <= executor._max_workers:
        return int(n)
    else:
        return int(executor._max_workers)

def cpu2use(n):

    import multiprocessing
    max_cpus = multiprocessing.cpu_count()
    
    if int(n) != 0 and int(n) <= max_cpus:
        return int(n)
    else:
        return max_cpus

    
def prepare_tmp_files(indir, outdir):

    import pandas as pd
    from Bio import SeqIO
    
    #Check whether it exists a prt and a gff file for each genome
    list_prt_genomes = [x.rsplit(".",1)[0] for x in os.listdir(indir + "/Proteins")]
    list_gff_genomes = [x.rsplit(".",1)[0] for x in os.listdir(indir + "/gff3")]
    
    if len(list_prt_genomes) == len(list_gff_genomes):
        for g in list_prt_genomes:
            if g not in list_gff_genomes:
                print(f"Could not find the {g} gff file in the gff3 folder, please there are the same files in the Proteins and gff3 folder")
                exit()
    else : 
        print("There are a different number of files provided in the Proteins and gff3 folder, please make sure there are the same number of files")
        exit()
    
    
    #first creating the tmp gff file containing the data from all gff files
    list_df2concatenate = []
    
    for i in list_gff_genomes :
        df_tmp = pd.read_csv(f"{indir}/gff3/{i}.gff", sep = "\t", comment = "#", names = ["contig", "source", "type", "left_c", "right_c", ".", "strand", "0", "id"])
        df_tmp["genome"] = i
        list_df2concatenate.append(df_tmp)
    
    df_all_gff = pd.concat(list_df2concatenate)
    df_all_gff["id"] = df_all_gff["id"].apply(lambda x: x.split(";",1)[0].split("=")[-1])
    df_all_gff.to_csv(outdir + "/TAT_defense_results/tmp/all_gff_cat.gff", sep = "\t", index = False)

    #Doing the same for the prt files
    subprocess.run(f"cat {indir}/Proteins/*.prt > {outdir}/TAT_defense_results/tmp/all_prt_cat.prt", shell = True, check = True)
    prt_seqIO = SeqIO.parse(f"{outdir}/TAT_defense_results/tmp/all_prt_cat.prt", "fasta")
    str_all_prt = ""
    d_gene_names_formated = {}
    
    #Slight change to the id of the genes to allows the gembase search with defense finder (otherwise genes in the border of contig would be not correctly considered)
    for record in prt_seqIO:
        id_formated = record.id.split(" ")[0].split("_")[0][:-1]+"i_"+record.id.split("_")[1] #PanACoTA considers the genes which are located in the borders of the contig by adding a "b" but it may causes somes troubles later
        str_all_prt += f">{id_formated}\n{record.seq}\n"
        d_gene_names_formated[id_formated] = record.id.split(" ")[0]

    with open(f"{outdir}/TAT_defense_results/tmp/all_prt_formated.prt", "w") as f:
        f.write(str_all_prt)
    
    #saving the d_genes_names_formated in a file to use it later and quickly find the real gene names
    df_gene_names_formated = pd.DataFrame.from_dict(d_gene_names_formated, orient = "index")
    df_gene_names_formated.to_csv(f"{outdir}/TAT_defense_results/tmp/gene_names_formated.tsv", header = False, sep = "\t")