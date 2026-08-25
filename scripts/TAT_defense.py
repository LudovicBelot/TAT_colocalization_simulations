#!/usr/bin/env python

import argparse
import utilities


def main():
    """
    main function that use the other scripts to establish the spot-pangenome, search the defense, TAs, and MGEs.
    Then performs the simulation to determine whether its random or not
    """

    args = get_args()
    utilities.create_dir(args.outdir, to_run = args.script)
    utilities.prepare_tmp_files(args.indir, args.outdir)
    n_cpu = utilities.cpu2use(args.cpu)
    
    if args.script == "all":
        
        import core_mapping
        import TA_defense_search
        import simulations

        core_mapping.main(args.indir, args.corefile, args.outdir, ref_genome = args.ref, n_cpu = n_cpu)
        
        if args.custom_TA == False :
            list_to_search = ["TA", "defense", "IS", "integrase"]
        else:
            list_to_search = ["defense", "IS", "integrase"]
            
        TA_defense_search.main(args.indir, args.outdir, n_cpu = n_cpu, list_to_search = list_to_search, custom_TA = args.custom_TA)
        simulations.main(args.outdir, n_cpu = n_cpu)
    
    elif args.script == "core_mapping":
        import core_mapping
        core_mapping.main(args.indir, args.corefile, args.outdir, ref_genome = args.ref, n_cpu = n_cpu)
    
    elif args.script == "search":
        import TA_defense_search
        
        if args.custom_TA == False :
            list_to_search = ["TA", "defense", "IS", "integrase"]
        else:
            list_to_search = ["defense", "IS", "integrase"]
            
        TA_defense_search.main(args.indir, args.outdir, n_cpu = n_cpu, list_to_search = list_to_search, custom_TA = args.custom_TA)
    
    elif args.script == "simulation":
        import simulations
        simulations.main(args.outdir, n_cpu = n_cpu)




def get_args():
    
    parser = argparse.ArgumentParser()

    parser.add_argument("--script", "-s",
                        help = " (REQUIRED) Which scripts to run (either 'all', 'core_mapping', 'search', 'simulation')",
                        choices=['all', 'core_mapping', 'search', "simulation"], 
                        required = True)
    
    parser.add_argument("--indir", "-i",
                        help = " (REQUIRED) Directory with your genomes sequences (except if you are running only the simulation only in which case, it requires the result folder from the previous scripts)", required = True)
    
    parser.add_argument("--outdir", "-o",
                        help = " (REQUIRED) Out directory to save the results", required = True)
    
    parser.add_argument("--corefile", "-c",
                        help = " (REQUIRED) Corepers file from PanACoTA if 'all' or 'core_mapping', else the core feature file generated from the core_mapping module", required = True)    
    
    parser.add_argument("--ref", "-r",
                        help = " Reference genome to use for the spot pangenome mapping (default = 'First' genome in the gff3 folder)", default = "First")
    parser.add_argument("--custom_TA", "-T",
                        help = "Custom tsv file containing the locations of your custom TA systems in the REFERENCE genome only. Requires at least 3 columns named : 'Ref_TA', 'left_c', 'right_c'. Optionnally, there can be a 4th column with the category/family of the TAs (column name : 'Family').", default = False)
    parser.add_argument("--cpu",
                        help = " Number of cpu to use (default = 1, 0 for all cpus available)", default = 1)

    args = parser.parse_args()
    
    return args

if __name__ == "__main__":
    main()