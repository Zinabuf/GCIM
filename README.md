The genetic causality inference model(GCIM) is a statistical method for detecting the direction of causation in GxE interaction studies. 

- 
 Author list: Zinabu Fentaw, S.Hong Lee
- 
+  
Package installation
+

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

Depending on the type of outcome variables, the outcome, exposure, and confounder variable or other covariate data should be split into the discovery dataset(ideally 80% of the data ) and the target dataset(the remaining 20% of the data). 


Quick start

GCIM analysis uses PLink2 to analyze discovery data, and the package is compatible with the Linux operating system. 
1. download the plink2 from the Plink website and specify the executable Plink file path.
   
~~~
plink_path <- "<plink_path>/plink2"
~~~
2. Set the working directory and run the following R functions 
3. Check the combination of the outcome and exposure variable types
   
   3.1. Binary outcome with binary exposure variable any type of confounder for the proposed direction of causation.
   
This function performs genome-wide interaction studies (GWEIS), genome-wide association studies
(GWAS), polygenic risk score (PRS) computation, and regression analysis to determine causal
directions for binary outcomes and binary exposure variables in the proposed direction.

   3.1.1. Performing GWEIS

   ~~~
   bbp_gweis(plink_path, dis_snp, bp_dis_phen, bp_dis_cov, output_dir, confounders)
   ~~~

   3.1.2. Performing GWAS
   
   ~~~
   bp_gwas(plink_path, dis_snp, bp_dis_cov, output_dir, confounders)
   ~~~

   3.1.3. GCIM analysis for the proposed direction with binary outcomes and binary exposures.

   ~~~
     gcim_bbp <- glm(Outcome ~ Additive_PRS + Interaction_PRS + Cov_PRS +
                    Interaction_PRS:Cov_PRS + Confounders,
                  family = "binomial", data = regression_data)

  summary(gcim_bbp)
  ~~~
 
 3.2. Binary outcome with binary exposure variable any type of confounder for the reverse direction of causation.


This function performs genome-wide interaction studies (GWEIS), genome-wide association studies
(GWAS), polygenic risk score (PRS) computation, and regression analysis to determine causal
directions for binary outcomes and binary exposure variables in the reverse direction.

  3.2.1. Performing GWEIS

   ~~~
   bbr_gweis(plink_path, dis_snp, br_dis_cov, br_dis_phen, output_dir, confounders)
   ~~~

   3.2.2. Performing GWAS
   
   ~~~
   bp_gwas(plink_path, dis_snp, br_dis_phen, output_dir, confounders)
   ~~~

   3.2.3. GCIM analysis for the reverse direction with binary outcomes and binary exposures.

   ~~~
     gcim_bbr <- glm(Outcome ~ Additive_PRS + Interaction_PRS + Cov_PRS +
                    Interaction_PRS:Cov_PRS + Confounders,
                  family = "binomial", data = regression_data)

  summary(gcim_bbr)
   ~~~

3.3. Binary outcome with quantitative exposure variable any type of confounder for the proposed direction of causation.

 This function performs genome-wide interaction studies (GWEIS), genome-wide association studies
(GWAS), polygenic risk score (PRS) computation, and regression analysis to determine causal
directions for binary outcomes and quantitative exposure variables in the proposed direction.

    3.3.1. Performing GWEIS

   ~~~
   bqp_gweis(plink_path, dis_snp, bp_dis_phen, qp_dis_cov, output_dir, confounders)
   ~~~

   3.3.2. Performing GWAS
   
   ~~~
   qp_gwas(plink_path, dis_snp, qp_dis_cov, output_dir, confounders)
   ~~~

   3.3.3. GCIM analysis for the proposed direction with binary outcomes and quantitative exposures.

   ~~~
     gcim_bqp <- glm(Outcome ~ Additive_PRS + Interaction_PRS + Cov_PRS +
                    Interaction_PRS:Cov_PRS + Confounders,
                  family = "binomial", data = regression_data)

  summary(gcim_bqp)
   ~~~

3.4. Binary outcome with quantitative exposure variable any type of confounder for the reverse direction of causation.

 This function performs genome-wide interaction studies (GWEIS), genome-wide association studies
(GWAS), polygenic risk score (PRS) computation, and regression analysis to determine causal
directions for binary outcomes and quantitative exposure variables in the reverse direction.

 3.4.1. Performing GWEIS

   ~~~
   bqr_gweis(plink_path, dis_snp, qr_dis_cov, br_dis_phen, output_dir, confounders)
   ~~~

   3.4.2. Performing GWAS
    
   ~~~
   br_gwas(plink_path, dis_snp, br_dis_cov, output_dir, confounders)
   ~~~

   3.4.3. GCIM analysis for the reverse direction with binary outcomes and quantitative exposures.

   ~~~
     gcim_bqr <- lm(Outcome ~ Additive_PRS + Interaction_PRS + Cov_PRS +
                    Interaction_PRS:Cov_PRS + Confounders,
                   data = regression_data)

  summary(gcim_bqr)
  ~~~

3.5. Quantitative outcome with quantitative exposure variable any type of confounder for the proposed direction of causation.

 This function performs genome-wide interaction studies (GWEIS), genome-wide association studies
(GWAS), polygenic risk score (PRS) computation, and regression analysis to determine causal
directions for quantitative outcomes and quantitative exposure variables in the proposed direction.

   3.5.1. Performing GWEIS

   ~~~
   qqp_gweis(plink_path, dis_snp, qp_dis_phen, qp_dis_cov, output_dir, confounders)
   ~~~

   3.5.2. Performing GWAS
   
   ~~~
   qp_gwas(plink_path, dis_snp, qp_dis_cov, output_dir, confounders)
   ~~~

   3.5.3. GCIM analysis for the proposed direction with quantitative outcomes and quantitative exposures.

   ~~~
     gcim_qqp <- lm(Outcome ~ Additive_PRS + Interaction_PRS + Cov_PRS +
                    Interaction_PRS:Cov_PRS + Confounders,
                   data = regression_data)

  summary(gcim_qqp)
  ~~~

3.6. Quantitative outcome with quantitative exposure variable any type of confounder for the reverse direction of causation.

 This function performs genome-wide interaction studies (GWEIS), genome-wide association studies
(GWAS), polygenic risk score (PRS) computation, and regression analysis to determine causal
directions for quantitative outcomes and quantitative exposure variables in the proposed direction.

 3.6.1. Performing GWEIS

   ~~~
   qqr_gweis(plink_path, dis_snp, qr_dis_cov, qr_dis_phen, output_dir, confounders)
   ~~~

   3.6.2. Performing GWAS
   
   ~~~
   qr_gwas(plink_path, dis_snp, qr_dis_phen, output_dir, confounders)
   ~~~

   3.6.3. GCIM analysis for the reverse direction with quantitative outcomes and quantitative exposures.

   ~~~
     gcim_qqr <- lm(Outcome ~ Additive_PRS + Interaction_PRS + Cov_PRS +
                    Interaction_PRS:Cov_PRS + Confounders,
                  data = regression_data)

  summary(gcim_qqr)
  ~~~

3.7. Quantitative outcome with binary exposure variable any type of confounder for the proposed direction of causation.

 This function performs genome-wide interaction studies (GWEIS), genome-wide association studies
(GWAS), polygenic risk score (PRS) computation, and regression analysis to determine causal
directions for quantitative outcomes and Binary exposure variables in the proposed direction.
 3.7.1. Performing GWEIS

   ~~~
   qbp_gweis(plink_path, dis_snp, bp_dis_phen, bp_dis_cov, output_dir, confounders)
   ~~~

   3.7.2. Performing GWAS
   
   ~~~
   qp_gwas(plink_path, dis_snp, bp_dis_cov, output_dir, confounders)
   ~~~

   3.7.3. GCIM analysis for the proposed direction with quantitative outcomes and binary exposures.

   ~~~
     gcim_qbp <- lm(Outcome ~ Additive_PRS + Interaction_PRS + Cov_PRS +
                    Interaction_PRS:Cov_PRS + Confounders,
                    data = regression_data)

  summary(gcim_qbp)
   ~~~

3.8. Quantitative outcome with binary exposure variable any type of confounder for the proposed direction of causation.

 3.8.1. Performing GWEIS

   ~~~
   qbr_gweis(plink_path, dis_snp, br_dis_cov, qr_dis_phen, output_dir, confounders)
   ~~~

   3.8.2. Performing GWAS
   
   ~~~
   qr_gwas(plink_path, dis_snp, bp_dis_cov, output_dir, confounders)
   ~~~

   3.8.3. GCIM analysis for the proposed direction with binary outcomes and binary exposures.

   ~~~
     gcim_qbr <- glm(Outcome ~ Additive_PRS + Interaction_PRS + Cov_PRS +
                    Interaction_PRS:Cov_PRS + Confounders,
                  family = "binomial", data = regression_data)

  summary(gcim_qbr)
   ~~~




