The genetic causality inference model(GCIM) is a statistical method for detecting the direction of causation in GxE interaction studies. 

- 
 Authors: Zinabu Fentaw, Dovini Jayasinghe, S.Hong Lee
-

NB: The proposed direction of causation refers to the causal directions of GxE interactions that is the primary focus of the researcher's interest, while the reverse direction of causation examines the opposite directions of GxE interactions to test its validity.
   
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

The data preparation follows: All data files should be split into two files for discovery and target data including Genetic data, Outcome(phenotype data), exposure(environmental data), and confounder variables. 

Genetic data 

The genetic data should be in Plink binary format(.bed, .bim, and .fam). Then it should be split into the discovery dataset(ideally 80% of the data ) and the target dataset(the remaining 20% of the data). 

Outcome data, exposure or environmental variables, and confounder variables are also split in proportions similar to those above and are compatible with the Plink data format. 
    #outcome should contain three columns ( FID, IID, and phenotype Value) and the phenotype value for case-control data should be 1 for control and 2 for case and in the target dataset, it should be 0 for controls and 1 for cases.
    #the exposure should contain at least 19 columns (FID, IID, and exposure values, confounder 1 confounder1, ...confounder16)
    

Depending on the type of outcome variables, the outcome, exposure, and confounder variable or other covariate data should be split into the discovery dataset(ideally 80% of the data ) and the target dataset(the remaining 20% of the data). 


Quick start

GCIM analysis uses PLink2 to analyze discovery data, and the package is compatible with the Linux operating system. 
1. download the plink2 from the Plink website and specify the executable Plink file path.
   
~~~
plink_path <- "<plink_path>/plink2"
~~~
2. Set the working directory and run the following R functions

  ~~~
output_dir <- "<output_path>/output_dir"
 ~~~

4. Check the combination of the outcome and exposure variable types
   
   3.1. Binary outcome with binary exposure variable any confounder for the proposed direction of causation.
   
This function performs genome-wide interaction studies (GWEIS), genome-wide association studies
(GWAS), polygenic risk score (PRS) computation, and regression analysis to determine causal
directions for binary outcomes and binary exposure variables in the proposed direction.


   3.1.1. Performing GWEIS

   ~~~
  a <- bbp_gweis(plink_path, dis_snp, bp_dis_phen, bp_dis_cov, output_dir)
   ~~~


   3.1.2. Performing GWAS

   ~~~
  b <- bp_gwas(plink_path, dis_snp, bp_dis_cov, output_dir)
   ~~~

  3.1.3. Compute Polygenic Risk Scores (PRS) for bbp

   ~~~
  c <- bbp_prs(plink_path, tar_snp, output_dir)
   ~~~

   3.1.4. GCIM analysis for the proposed direction with binary outcomes and exposures.

   ~~~
  d <- gcim_bbp(bp_tar_phen, bp_tar_cov, Add_PRS, Int_PRS, Cov_PRS, confounders)
  ~~~

 Open model output as 

 ~~~
  print(d)
 ~~~

 3.2. Binary outcome with binary exposure variable any confounder for the reverse direction of causation.


This function performs genome-wide interaction studies (GWEIS), genome-wide association studies
(GWAS), polygenic risk score (PRS) computation, and regression analysis to determine causal
directions for binary outcomes and binary exposure variables in the reverse direction.

  3.2.1. Performing GWEIS

   ~~~
   e <-  bbr_gweis(plink_path, dis_snp, br_dis_cov, br_dis_phen, output_dir)
   ~~~

   3.2.2. Performing GWAS
   
   ~~~
   f <- br_gwas(plink_path, dis_snp, br_dis_phen, output_dir)
   ~~~
  3.2.3. Compute Polygenic Risk Scores (PRS) for bbr

  ~~~
   g <- bbr_prs(plink_path, tar_snp, output_dir)
  ~~~~

   3.2.4. GCIM analysis for the reverse direction with binary outcomes and exposures.

   ~~~
    h <- gcim_bbr(br_tar_cov, br_tar_phen, Add_PRS, Int_PRS, Cov_PRS, confounders)
   ~~~

3.3. Binary outcome with quantitative exposure variable any type of confounder for the proposed direction of causation.

 This function performs genome-wide interaction studies (GWEIS), genome-wide association studies
(GWAS), polygenic risk score (PRS) computation, and regression analysis to determine causal
directions for binary outcomes and quantitative exposure variables in the proposed direction.

    3.3.1. Performing GWEIS

   ~~~
   bqp_gweis(plink_path, dis_snp, bp_dis_phen, qp_dis_cov, output_dir)
   ~~~

   3.3.2. Performing GWAS
   
   ~~~
   qp_gwas(plink_path, dis_snp, qp_dis_cov, output_dir)
   ~~~

  3.3.3. Compute Polygenic Risk Scores (PRS) for bqp

  ~~~
  bqp_prs(plink_path, tar_snp, output_dir)
  ~~~~


   3.3.4. GCIM analysis for the proposed direction with binary outcomes and quantitative exposures.

   ~~~
gcim_bqp(bp_tar_phen, qp_tar_cov, prs_add_scaled, prs_int_scaled, prs_cov_scaled, confounders)
   ~~~

3.4. Binary outcome with quantitative exposure variable any confounder for the reverse direction of causation.

 This function performs genome-wide interaction studies (GWEIS), genome-wide association studies
(GWAS), polygenic risk score (PRS) computation, and regression analysis to determine causal
directions for binary outcomes and quantitative exposure variables in the reverse direction.

 3.4.1. Performing GWEIS

   ~~~
   bqr_gweis(plink_path, dis_snp, qr_dis_cov, br_dis_phen, output_dir)
   ~~~

   3.4.2. Performing GWAS
    
   ~~~
   br_gwas(plink_path, dis_snp, br_dis_phen, output_dir)
   ~~~

  3.4.3. Compute Polygenic Risk Scores (PRS) for bqr

~~~
bqr_prs(plink_path, tar_snp, output_dir)
~~~~


   3.4.4. GCIM analysis for the reverse direction with binary outcomes and quantitative exposures.

   ~~~
gcim_bqr(qr_tar_cov, br_tar_phen, prs_add_scaled, prs_int_scaled, prs_cov_scaled, confounders)
  ~~~

3.5. Quantitative outcome with quantitative exposure variable any confounder for the proposed direction of causation.

 This function performs genome-wide interaction studies (GWEIS), genome-wide association studies
(GWAS), polygenic risk score (PRS) computation, and regression analysis to determine causal
directions for quantitative outcomes and quantitative exposure variables in the proposed direction.

   3.5.1. Performing GWEIS

   ~~~
 qqp_gweis(plink_path, dis_snp, qp_dis_phen, qp_dis_cov, output_dir)
   ~~~

   3.5.2. Performing GWAS
   
   ~~~
   qp_gwas(plink_path, dis_snp, qp_dis_cov, output_dir, confounders)
   ~~~

  3.5.3. Compute Polygenic Risk Scores (PRS) for qqp

~~~
qqp_prs(plink_path, tar_snp, output_dir)
~~~~

   3.5.4. GCIM analysis for the proposed direction with quantitative outcomes and exposures.

   ~~~
gcim_qqp(qp_tar_phen, qp_tar_cov, prs_add_scaled, prs_int_scaled, prs_cov_scaled, confounders)
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
   qr_gwas(plink_path, dis_snp, qr_dis_phen, output_dir)
   ~~~

  3.6.3. Compute Polygenic Risk Scores (PRS) for qqr

~~~
qqr_prs(plink_path, tar_snp, output_dir)
~~~~

   3.6.4. GCIM analysis for the reverse direction with quantitative outcomes and quantitative exposures.

   ~~~
gcim_qqr(qr_tar_cov, qr_tar_phen, prs_add_scaled, prs_int_scaled, prs_cov_scaled, confounders)
  ~~~

3.7. Quantitative outcome with binary exposure variable any type of confounder for the proposed direction of causation.

 This function performs genome-wide interaction studies (GWEIS), genome-wide association studies
(GWAS), polygenic risk score (PRS) computation, and regression analysis to determine causal
directions for quantitative outcomes and binary exposure variables in the proposed direction.

 3.7.1. Performing GWEIS

   ~~~
   qbp_gweis(plink_path, dis_snp, qp_dis_phen, bp_dis_cov, output_dir)
   ~~~

   3.7.2. Performing GWAS
   
   ~~~
   bp_gwas(plink_path, dis_snp, bp_dis_cov, output_dir)
   ~~~

  3.7.3. Compute Polygenic Risk Scores (PRS) for qbp

~~~
qbp_prs(plink_path, tar_snp, output_dir)
~~~~

   3.7.4. GCIM analysis for the proposed direction with quantitative outcomes and binary exposures.

   ~~~
gcim_qbp(qp_tar_phen, bp_tar_cov, prs_add_scaled, prs_int_scaled, prs_cov_scaled, confounders)
   ~~~

3.8. Quantitative outcome with binary exposure variable any confounder for the reverse direction of causation.

 This function performs genome-wide interaction studies (GWEIS), genome-wide association studies
(GWAS), polygenic risk score (PRS) computation, and regression analysis to determine causal
directions for quantitative outcomes and binary exposure variables in the reverse direction.

 3.8.1. Performing GWEIS

   ~~~
  qbr_gweis(plink_path, dis_snp, br_dis_cov, qr_dis_phen, output_dir)
   ~~~

   3.8.2. Performing GWAS
   
   ~~~
   qr_gwas(plink_path, dis_snp, qr_dis_phen, output_dir)
   ~~~

  3.8.3. Compute Polygenic Risk Scores (PRS) for qbr

~~~
qbr_prs(plink_path, tar_snp, output_dir)
~~~~


   3.8.4. GCIM analysis for the reverse direction with quantitative outcomes and binary exposures.

   ~~~
gcim_qbr(br_tar_cov, qr_tar_phen, prs_add_scaled, prs_int_scaled, prs_cov_scaled, confounders)
   ~~~




