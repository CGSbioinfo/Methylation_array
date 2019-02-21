# Pipeline for the analysis of the Epic array

The pipeline works using 2 scripts, the function scripts is loaded first and then can be called from the command line or another scripts.
See the description below for the arguments available for each function.

# Description of the functions  

## load_data
This function is to load the data, the requirements are:
* Folder containing the **.idats** files
* Sample sheet (GenomeStudio format), matching the arrays to sample names, the sample sheet should contain in addition to the required field:
  * **Gender** information if available
  * **Sample_Group** field for the comparisons
 
**Arguments**:
* *working directory*: directory where the analysis will take place, should contains a folder with idats file and samplesheet
* *sampleSheetPattern*: pattern to find the samplesheet, by default the function will look for any ".csv" file in the folder, if several are present it will crash. Use the pattern to identify the one you want to use.
  
The function return an rgSet object (Red and Green set)

## preprocess_data
This function takes an rgSet (typically the output from load_data) and runs several steps to generated a gmSet:
* *QCinfo*: generate QC information for the data
* *preprocessENmix*: run the preprocessing step of the ENmix package, by default it using the background estimation from the *oob* method, other method available are *est* or *neg*, more details below. The dye biass correction is performed using the *RELIC* method, other methods available are *mean* or *none*. parameter that can be set is the number of cores to use, by default it is set to 10, use *nCores* to adjust it.
* Outlier CpG measurements are remove using the rm.outlier method
* Normalisation of the set is performed using the *preProcessQuantile* method from the minfi package
* Probes with SNPs in their sequence are removed from the data set.

This function return an gmSet. Note on the choice of background correction method: *oob* stands for out-of-band and ENmix will use out-of-band Infinium I intensities ("oob") to estimate normal distribution parameters to model background noise. Option "est" will use combined methylated and unmethylated intensities to estimate background distribution parameters separately for each color channel and each probe type. Option "neg" will use 600 chip internal controls probes to estimate background distribution parameters.
For the dye bias method, the *mean* method correction is based on average of red/green ration. The *RELIC* method is based on the method described in: *Zongli Xu, Sabine A. S. Langie, Patrick De Boever, Jack A. Taylor1 and Liang Niu, RELIC: a novel dye-bias correction method for Illumina Methylation BeadChip, BMC Genomics. 2017; 18: 4.*. 







