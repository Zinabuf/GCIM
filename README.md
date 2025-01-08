The genetic causality inference model(GCIM) is a statistical method for detecting the direction of causation in GxE interaction studies. 

Author list: Zinabu Fentaw S.Hong Lee

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
The genetic data should be in Plink binary format(.bed, .bim, and .fam) then it should be split into the discovery dataset(ideally 80% of the data ) and the target dataset(the remaining target dataset)

