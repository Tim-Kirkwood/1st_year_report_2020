# 1st_year_report_2020

## “Cluster Architecture” section
Functional annotations for all DEM32671 cluster proteins were extracted using the “ClusterBlast function annotation”.  Run antiSMASH, download the antiSMASH .zip file and extract it. Then move the single .JSON file to the directory containing your script.  Files not included in this repository as they contain confidential information.  Script was developed for the antiSMASH output of the highly annotated genome used in this report, and may not work for a generic NCBI/RAST annotated genome.

## “Bioinformatic analysis of pbtX” section
The bioinformatic analysis described in the indicated section used the script “pbt ClusterBlast analyser”. To run this code, the antiSMASH file generated after analysis of the NCBI geneome discussed in the report must be extracted, and the .json file placed in the same directory as this script.

## “Mining MIBIG for characterised cluster trends” section
The bioinformatic analysis described in the indicated section used the script “MIBIG miner”. The JSON-format MIBIG database should be downloaded from https://mibig.secondarymetabolites.org/download. The code below should then be executed, ensuring the script is in the same directory as the JSON files. The file path for the .JSON file will need to be entered in line 12 (‘path’ variable).

## “Multiplex PCR primer design” section
The workflow described in this section was run using the python script “PrimerBLAST multiplexer”.  URLs linking to each primerBLAST run should be entered into the script in lines 47 to 58. The URLs will expire, so this program should be used within 24 hours of the original primerBLAST run. Expired URLS are included here as an example, and were generated from primerBLAST runs generating primers of 20, 25 and 30bp for each cluster. 


# REFERENCES

## AntiSMASH reference: 
antiSMASH 5.0: updates to the secondary metabolite genome mining pipeline
Kai Blin, Simon Shaw, Katharina Steinke, Rasmus Villebro, Nadine Ziemert, Sang Yup Lee, Marnix H Medema, & Tilmann Weber
Nucleic Acids Research (2019) doi: 10.1093/nar/gkz310.

## MIBIG refernce:
Satria A Kautsar, Kai Blin, Simon Shaw, Jorge C Navarro-Muñoz, Barbara R Terlouw, Justin J J van der Hooft, Jeffrey A van Santen, Vittorio Tracanna, Hernando G Suarez Duran, Victòria Pascal Andreu, Nelly Selem-Mojica, Mohammad Alanjary, Serina L Robinson, George Lund, Samuel C Epstein, Ashley C Sisto, Louise K Charkoudian, Jérôme Collemare, Roger G Linington, Tilmann Weber, Marnix H Medema, MIBiG 2.0: a repository for biosynthetic gene clusters of known function, Nucleic Acids Research, Volume 48, Issue D1, 08 January 2020, Pages D454–D458, https://doi.org/10.1093/nar/gkz882

## PrimerBLAST reference:
Ye J, Coulouris G, Zaretskaya I, Cutcutache I, Rozen S, Madden T (2012).
Primer-BLAST: A tool to design target-specific primers for polymerase chain reaction.
BMC Bioinformatics. 13:134.
