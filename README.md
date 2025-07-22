---
# <h1 align="center">GCIM</h1>
---

The genetic causality inference model(GCIM) is a statistical method for detecting the direction of causation in GxE interaction studies. 
- 
 Authors: Zinabu Fentaw, Dovini Jayasinghe, S.Hong Lee
-

NB: The proposed direction of causation refers to the causal directions of GxE interactions that are the primary focus of the researcher's interest, while the reverse direction of causation examines the opposite directions(changing the role of outcome and exposure variables) of GxE interactions to test their proper causal directions.
   
Package installation 
From GitHub 

~~~
library(devtools)
install_github("DoviniJ/GxEprs")
install_github("Zinabuf/GCIM")
~~~
Or the CRAN version via
~~~
install.packages("GxEprs") 
~~~

Load the library

~~~
library(GxEprs)
~~~

**Data preparations**
The dataset is divided into **discovery (80%)** and **target (20%)** datasets, ensuring consistency across genetic, outcome, exposure, and confounder data. Genetic data, stored in **PLINK binary format** (`.bed`, `.bim`, `.fam`). The outcome variable should include **FID, IID, and phenotype values**, where case-control labels follow PLINK conventions: **1 = Control, 2 = Case** in the discovery dataset, and **0 = Control, 1 = Case** in the target dataset. Exposure and confounder variables are formatted into at least **m columns** (**FID, IID, exposure, exposure_sqaure, confounder1-m **) and partitioned in the same proportions. For a binary outcome, the fourth column (exposure_sqaure) should be constant values(zero variance) such as 2,2,2, ... 2  for all observations. This structured approach ensures compatibility across all data types, thereby maintaining accurate alignment to estimate GxE interactions.


A Guide for GCIM analyses

GCIM analyses use PLink2 to analyze discovery data, and the package is compatible with the Linux operating system. 
1. Download plink2 software from the [Plink](https://www.cog-genomics.org/plink/2.0/) website and then specify the executable Plink file path.
   
~~~
plink_path <- "<plink_path>/plink2"
~~~
**Discovery dataset**: In the discovery dataset, please use the data input from the GxEprs data input format, then construct the PRS of the exposure using [Plink](https://www.cog-genomics.org/plink/2.0/) based on the GWAS summary statistics for the PRS of the target samples.
**Target dataset**: The target dataset for the model is also similar in data format, except, the use of PRS of the exposure variable(PRS of E) rather than the use of entire exposure values, and also includes the constant values for the fourth column, which is the square of the third column in the GxEprs model. 

~~~
