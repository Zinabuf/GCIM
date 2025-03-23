# GCIM
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

**Data preparations**
The dataset is divided into **discovery (80%)** and **target (20%)** subsets, ensuring consistency across genetic, outcome, exposure, and confounder data. Genetic data, stored in **PLINK binary format** (`.bed`, `.bim`, `.fam`). The outcome variable should include **FID, IID, and phenotype values**, where case-control labels follow PLINK conventions: **1 = Control, 2 = Case** in the discovery dataset, and **0 = Control, 1 = Case** in the target dataset. Exposure and confounder variables are formatted into at least **19 columns** (**FID, IID, exposure, confounder1–16**) and partitioned in the same proportions. This structured approach ensures compatibility across all data types, maintaining alignment to estimate GxE interactions accurately.


A quick guide

GCIM analysis uses PLink2 to analyze discovery data, and the package is compatible with the Linux operating system. 
1. download the plink2 from the [Plink](https://www.cog-genomics.org/plink/2.0/) website and specify the executable Plink file path.
   
~~~
plink_path <- "<plink_path>/plink2"
~~~

Specify the plink file prefix, which reduces error in R

~~~
dis_snp <- system.file("data", "dis_snp", package = "GCIM")
tar_snp <- system.file("data", "tar_snp", package = "GCIM")
~~~
 Specifying **inst** directories

 ~~~
 inst_path <- system.file(package = "GCIM") 
~~~

Inst files are available in the directories

~~~
inst_path <- system.file(package = "GCIM") 
bp_dis_cov <- paste0(inst_path, "bp_dis_cov.txt")
bp_dis_phen <- paste0(inst_path, "/bp_dis_phen.txt")
bp_tar_cov <- paste0(inst_path, "/bp_tar_cov.txt")
bp_tar_phen <- paste0(inst_path, "/bp_tar_phen.txt")
br_dis_cov <- paste0(inst_path, "/br_dis_cov.txt")
br_dis_phen <- paste0(inst_path, "/br_dis_phen.txt")
br_tar_cov <- paste0(inst_path, "/br_tar_cov.txt")
br_tar_phen <- paste0(inst_path, "/br_tar_phen.txt")
qp_dis_cov <- paste0(inst_path, "/qp_dis_cov.txt")
qp_dis_phen <- paste0(inst_path, "/qp_dis_phen.txt")
qp_tar_cov <- paste0(inst_path, "/qp_tar_cov.txt")
qp_tar_phen <- paste0(inst_path, "/qp_tar_phen.txt")
qr_dis_cov <- paste0(inst_path, "/qr_dis_cov.txt")
qr_dis_phen <- paste0(inst_path, "/qr_dis_phen.txt")
qr_tar_cov <- paste0(inst_path, "/qr_tar_cov.txt")
qr_tar_phen <- paste0(inst_path, "/qr_tar_phen.txt")
dis_snp <- paste0(inst_path, "/dis_snp")
tar_snp <- paste0(inst_path, "/dis_snp")
~~~
Set the proposed directions and the reverse direction based on the type of outcome variables. Depending on the type of outcome variable, whether binary or quantitative, step-by-step tests can be conducted as a one-time process for one-direction tests.
For instance, if you have a binary outcome with a quantitative exposure variable, follow these steps.

**A. Proposed causal directions**

- Step 1: Conduct b_gweis
- Step 2: Conduct q_gwas
- Step 3: Compute prs_score
- Step 4: Compute a regression using gcim_b
  
**B. The reverse causal direction**

- Step 1: Conduct q_gweis
- Step 2: Conduct b_gwas
- Step 3: Compute prs_score
- Step 4: Compute a regression using gcim_q
Compare the statistical test results between the proposed and reverse causal directions and declare the correct directions.

**Input data**

All the input data should be prepared with the following data dimensions
**Discovery dataset**
The outcome variable in the discovery dataset should be structured with at least three columns: **FID**, **IID**, and **Phenotype**. For binary traits, controls and cases must be coded as **1** and **2**, respectively, while exposure and covariate values should be assigned according to their respective distributions. The corresponding genotype data for the discovery dataset should be formatted in **PLINK binary format** (.bed, .bim, .fam) to ensure compatibility with genome-wide-by-environment interaction analyses.
**Target dataset** 
The outcome variable in the target dataset should be structured with at least three columns: **FID**, **IID**, and **Phenotype**. For binary traits, controls and cases must be coded as **0** and **1**, respectively, while exposure and covariate values should be assigned according to their respective distributions. The corresponding genotype data for the target dataset should be formatted in **PLINK binary format** (.bed, .bim, .fam) to ensure compatibility with Polygenic risk score estimation.

**Start analysis**

**Performing GWEIS**

   ~~~
 a <- b_gweis(plink_path, dis_snp, bp_dis_phen, bp_dis_cov)
   ~~~

   ~~~
 b <- q_gweis(plink_path, dis_snp, qp_dis_phen, qp_dis_cov)
   ~~~

**Performing GWAS**

   ~~~
 c <- b_gwas(plink_path, dis_snp, bp_dis_cov)
   ~~~

   ~~~
 d <- q_gwas(plink_path, dis_snp, qp_dis_cov)
   ~~~

**Compute PRS**

   ~~~
 e <- prs_scores(plink_path, tar_snp)
   ~~~

**Compute GCIM for the required directions**

   ~~~
f <- gcim_b(bp_tar_phen, bp_tar_cov, Add_PRS, Int_PRS, Cov_PRS, confounders)
   ~~~

   ~~~
g <- gcim_q(qp_tar_phen, qp_tar_cov, Add_PRS, Int_PRS, Cov_PRS, confounders)
   ~~~

**Tsee regression outputs**

~~~
print(f)
~~~

~~~
print(g)
~~~

**Example of the analysis from the data in the package**
Quantitative outcome and quantitative exposure in testing in two different causal directions. 
**Proposed causal directions**
Quantitative outcome ~~~qp_dis_phen~~~, quantitative exposure ~~~qp_dis_cov~~~ and genotype data for the discovery dataset, while in the target dataset,
Quantitative outcome ~~~qp_tar_phen~~~, quantitative exposure ~~~qp_tar_cov~~~. 


