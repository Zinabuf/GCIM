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
3. Check the combination of the outcome and exposure variable types
   3.1. Binary outcome with binary exposure variable any type of confounder for the proposed direction of causation
This function performs genome-wide interaction studies (GWEIS), genome-wide association studies
(GWAS), polygenic risk score (PRS) computation, and regression analysis to determine causal
directions for binary outcomes and binary exposure variables in the reverse direction.

   3.1.1. Performing GWEIS

   ~~~
   bbp_gweis(plink_path, dis_snp, bp_dis_phen, bp_dis_cov, output_dir, confounders)
   ~~~
   3.1.2. Performing GWAS
   ~~~
   bbp_gweis(plink_path, dis_snp, bp_dis_cov, output_dir, confounders)
   ~~~
   3.1.4. GCIM computing estimation for the proposed direction.

   ~~~
     gcim_bbp <- glm(Outcome ~ Additive_PRS + Interaction_PRS + Cov_PRS +
                    Interaction_PRS:Cov_PRS + Confounders,
                  family = "binomial", data = regression_data)

  return(summary(gcim_bbp))
}
  ~~~ 
 3.2. Binary outcome with binary exposure variable any type of confounder for the reverse direction of causation.
 
 
3.3. Binary outcome with quantitative exposure variable any type of confounder for the proposed direction of causation. 
   
   
3.4. Binary outcome with quantitative exposure variable any type of confounder for the reverse direction of causation.


3.5. quantitative outcome with quantitative exposure variable any type of confounder for the proposed direction of causation.


3.6. quantitative outcome with quantitative exposure variable any type of confounder for the proposed direction of causation.


3.7. quantitative outcome with quantitative exposure variable any type of confounder for the proposed direction of causation.


3.8. quantitative outcome with quantitative exposure variable any type of confounder for the proposed direction of causation.







