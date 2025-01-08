The genetic causality inference model(GCIM) is a statistical method for detecting the direction of causation in GxE interaction studies. 

Author list: Zinabu Fentaw, S.Hong Lee

Package installation
~~~
library(devtools)
install_github("zinabuf/GCIM")
~~~

Load the library
~~~
library(GCIM)
~~~
Data preparation for input data files

The data preparation follows: All data files should be split into two files for discovery and target data for Genetic data, Outcome(phenotype data), exposure(environmental data), and confounder variables. 

Genetic data 

The genetic data should be in Plink binary format(.bed, .bim, and .fam), and then it should be split into the discovery dataset(ideally 80% of the data ) and the target dataset(the remaining 20% of the data)

Outcome data, Exposure or Environmental variables, and confounder variables

The outcome, exposure, and confounder variable or other covariate data should be split into the discovery dataset(ideally 80% of the data ) and the target dataset(the remaining 20% of the data) depending on the type of outcome variables. 


Quick start

GCIM analysis uses PLink2 for the analysis of discovery data and the package is compatible with the Linux operating system. 
1. download the plink2 from the Plink website and then specify the executable Plink file path.
   
~~~
plink_path <- "<plink_path>/plink2"
~~~
2. Set the working directory and run the following R functions
3. 






