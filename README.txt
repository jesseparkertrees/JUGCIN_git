# JUGCIN_git
This is the readme for the "JUGCIN_git" repo.
Investigating range-wide patterns of hybridization in Juglans cinerea.

Included are files for data processing, analysis, and visualization.
Data includes range wide butternut genomic data (GBS and microsatellites), butternut health data from select populations, and geospatial datasets.

To view summary tables and figures, visit "./summaries/".

To reproduce analyses, visit "./misc/DIRECTIONS.txt" for instructions on how to proceed.



SUMMARY OF DIRECTORIES AND FILES (not currently up to date 7/16):

"./archive/"	#includes archived files that are not currently relevant

".data/" 	#includes raw data sets and inputs for other analyses
	1. "ButternutLandscapeGPS_Aziz07022026.csv" 
		#table of Aziz GBS samples with coordinates and sample information
	2. "genotype_cats_31.txt"
		#table of genotype frequency classes for running NewHybrids with 31 categories
	3. "genotype_cats_37.txt"
		#table of genotype frequency classes for running NewHybrids with 37 categories
	4. "JUGCIN_Emma_Collected_depDel.csv"
		#table of field collected health data compiled from Hoban lab data (2024-2025) 
	5. "JUGCIN_msatMastSampleSheet.csv"
		#table of Hoban microsatellite samples (circa 2010) with coordinates and sample information
	6. "msat_popfile"
		#population IDs for Hoban microsatellite data
	7. "newhybrids_input_withZ.txt"
		#NewHybrids formatted genotype table for Hoban microsatellites
	8. "Nuclear_Data_F_5_1_10good.csv"
		#genotype data for Hoban microsatellites 
	9. "project_data.str" 
		#structure formatted data for Hoban microsatellites

"./misc/"
	1. "DIRECTIONS.txt"
		#instructions for how to run analyses, order for executing script files, etc.
	2. "fastSTRUCTUREparams.txt"
		#code to be executed in command line for running fastSTRUCTURE on Aziz GBS samples
	3. "msatNewhybridsID_key.csv" 
		#key to Hoban microsatellite NewHybrids IDs and -z parameters
	4. "msatNEWHYBRIDSparams.txt"
		#code to executed in command line to run NewHybrids on Hoban microsatellite samples
	5. "msatSTRUCTUREparams.txt"
		#parameters used for running STRUCTURE in the GUI for Hoban microsatellite samples
 	6. "NEWHYBRIDSparams.txt"
		#code to run NewHybrids in the command line for Aziz GBS samples

"./outputs/"	#includes intermediate and final outputs from STRUCTURE and NewHybrids analyses
	1. "structure/"
		1. ...
	2. "msatNEWHYBRIDS/"
		#contains outputs from running NewHybrids on Hoban microsatellite samples
		1. "31cats/"
			#NewHybrids output using "genotype_cats_31.txt" for the genotype category table
		2. "aa-*" and "EchoedGtypData.txt"
			#NewHybrids outputs using the default 6 hybrid classes
	3. "gbsNEWHYBRIDS/"
		#contain outputs from running NewHybrids on Aziz GBS samples
		1. "37cats_run*/"
			#directories containing the NewHybrids outputs (using "genotype_cats_37.txt" for the genotype categories) 
		2. "inputs/" 
			#includes 5 sets of 200 randomly sampled SNPs for NewHybrids (code for generating these is in "./scripts/JUGCIN_genetic.R"
		3. "run*/"
			#directories containing the NewHybrids outputs (using "genotype_cats_31.txt" for the genotype categories)

"./scripts/"
	1. "JUGCIN_genetic.R"
		#R script for loading GBS genotype data, generating files for STRUCTURE and NewHybrids, and analyzing and summarizing results from STRUCTURE and NewHybrids
	2. "JUGCIN_health.R"
		#R script for analyzing site level health data 
	3. "JUGCIN_prepareData.R"
		#R script for compiling data from raw data files, extracting environmental variables, filtering data, etc.
	4. "JUGCIN_processMsatData.R"
		#R script for processing Hoban microsatellite genotype table, generating files for STRUCTURE and NewHybrids, and summarizing and analyzing their results
	5. "JUGCIN_relatedness.R"
		#R script for CERVUS relatedness analysis; generating CERVUS files, importing CERVUS results and analyzing, etc.
	6. "JUGCIN_spatial.R"	
		#R script for spatial analyses of hybridization

"./summaries/"
	1. "NewHybrids_categories_generations.xlsx"
		#summary of 37 NewHybrids categories with minimum number of generations to produce each
	2. "NewHybridsCategoriesSummaryTableALL.csv"
		#summary of NewHybrids category counts
	3. "NewHybridsCategoriesSummaryTableCounties.csv"
		#summary of NewHybrids category counts by county
	4. "NewHybridsCategoriesSummaryTableStates.csv"
		#summary of NewHybrids category counts by state
	5. ... other table and figures
