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
The dataset is divided into **discovery (80%)** and **target (20%)** subsets, ensuring consistency across genetic, outcome, exposure, and confounder data. Genetic data, stored in **PLINK binary format** (`.bed`, `.bim`, `.fam`). The outcome variable should include **FID, IID, and phenotype values**, where case-control labels follow PLINK conventions: **1 = Control, 2 = Case** in the discovery dataset, and **0 = Control, 1 = Case** in the target dataset. Exposure and confounder variables are formatted into at least **19 columns** (**FID, IID, exposure, confounder1-n **) and partitioned in the same proportions. This structured approach ensures compatibility across all data types, thereby maintaining accurate alignment to estimate GxE interactions.


A Guide for GCIM analyses

GCIM analyses use PLink2 to analyze discovery data, and the package is compatible with the Linux operating system. 
1. Download plink2 software from the [Plink](https://www.cog-genomics.org/plink/2.0/) website and then specify the executable Plink file path.
   
~~~
plink_path <- "<plink_path>/plink2"
~~~
**Discovery dataset**: In the discovery dataset, please use the data input from the GxEprs data input format, except the value of e(exposure variable)-squared, which should be prepared as constant values to remove its effect from the model. Compute the PRS of the exposure using [Plink](https://www.cog-genomics.org/plink/2.0/) based on the GWAS from the discovery dataset for the PRS of the target samples. Use input and data preparation from GxEprs.
**Target dataset**: The target dataset for the model is also similar in data format, except, the use of PRS of the exposure variable(PRS of E) rather than the use of entire exposure values, and also includes the constant values for the fourth column, which is the square of the third column in the GxEprs model. The same approach will be applied for the reverse direction test.
Example data:
Quantitative outcome with Quantitative exposure:
discovery dataset: 
1. Conduct GWEIS for the outcome
2. Conduct GWAS
3. Estimate the PRS of the exposure
4. Conduct linear regression

# quantitative outcome 
~~~
# Compute PRS of exposure variables
Conduct a GWAS for the exposure adjusted by covariates, and then compute the PRS for the exposure. 
~~~

~~~
# Conduct GWEIS based on the assigned outcome variables, exposure with confounders.
b <- GWEIS_binary(plink_path, "mydata", "Bpd.txt", "PRS_exp")
add <- b[c("ID", "A1", "ADD_BETA")]
gxe <- b[c("ID", "A1", "INTERACTION_BETA")]
~~~

~~~
# Compute the PRS of each summary 
q <- PRS_quantitative(plink_path, "mydata", summary_input = add)
r <- PRS_quantitative(plink_path, "mydata", summary_input = gxe)
~~~

~~~
# Conduct regression analyses for GCIM
#PRS of the exposure(PRS_exp) variables should be the third column, and include the constant values for the fourth column, and then add the  number of covariates for the adjustment. 
y <- summary_regular_quantitative("Qpt.txt", "Qct.txt", add_score = q, gxe_score = r, Model = 4)
~~~

~~~
y$summary
~~~
Perform the reverse analyses by changing the roles of exposure and outcome variables. 

# Binary outcome

**Compute PRS of exposure variables**
Conduct a GWAS for the exposure adjusted by covariates, and then compute the PRS for the exposure.
Conduct a log transformation of the OR if the exposure variable is binary.  Then, compute the PRS of the exposure variables using PRS score values. 

~~~
b <- GWEIS_binary(plink_path, "mydata", "Bpd.txt", "Bcd.txt")
add <- b[c("ID", "A1", "ADD_BETA")]
gxe <- b[c("ID", "A1", "INTERACTION_BETA")]
~~~

~~~
q <- PRS_binary(plink_path, "mydata", summary_input = add)
r <- PRS_binary(plink_path, "mydata", summary_input = gxe)
~~~

~~~
z <- summary_regular_binary("Bpt.txt", "add_score = p", add_score = q, gxe_score = r, Model = 5)
~~~

~~~
z$summary
~~~
