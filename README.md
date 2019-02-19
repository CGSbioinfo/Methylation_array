Pipeline for the analysis of the Epic array

The pipeline works using 2 scripts, the function scripts is loaded first and then can be called from the command line or another scripts.
See the description below for the arguments available for each function.

#Description of the functions
##load_data
* working directory: directory where the analysis will take place, should contains a folder with idats file and samplesheet
* sampleSheetPatter: pattern to find the samplesheet, by default the function will look for any ".csv" file in the folder, if several are present it will crash. Use the pattern to identify the one you want to use

The function return an rgSet object (Red and Green set)
