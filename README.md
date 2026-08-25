These scripts were written with the aim to :

1) Identify TA/Antiphage/IS/Integrases within the Photorhabdus genomes
2) Map them to a spot based on their genomic localization compared to the core genome.
3) Compare their observed distribution to the one expected through random spot assignation.

This allow to both determine whether every of these genetic elements are randomly distributed in the Photorhabdus genomic spots or not and whether they tend to associate with another type of element.

Main script :

=> TAT_defense.py

Arguments :

  "--script", "-s"
      help = " (REQUIRED) Which scripts to run (either 'all', 'core_mapping', 'search', 'simulation')"
      choices=['all', 'core_mapping', 'search', "simulation"]
      required = True

  "--indir", "-i"
      help = " (REQUIRED) Directory with your genomes sequences (except if you are running only the simulation only in which case, it requires the result folder from the previous scripts)"
      required = True

  "--outdir", "-o"
      help = " (REQUIRED) Out directory to save the results"
      required = True

  "--corefile", "-c"
      help = " (REQUIRED) Corepers file from PanACoTA if 'all' or 'core_mapping', else the core feature file generated from the core_mapping module"
      required = True

  "--ref", "-r"
      help = " Reference genome to use for the spot pangenome mapping (default = 'First' genome in the gff3 folder)"
      default = "First"

  "--custom_TA", "-T"
      help = "Custom tsv file containing the locations of your custom TA systems in the REFERENCE genome only. Requires at least 3 columns named : 'Ref_TA', 'left_c', 'right_c'. Optionnally, there can be a 4th column with the category/family of the TAs (column name : 'Family')."
      default = False

  "--cpu"
      help = " Number of cpu to use (default = 1, 0 for all cpus available)"
      default = 1


In addition, the script Custom_sim.py allows for more options than the simulation module if needed.

