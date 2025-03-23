# GCIM
The genetic causality inference model(GCIM) is a statistical method for detecting the direction of causation in GxE interaction studies. 

- 
 Authors: Zinabu Fentaw, Dovini Jayasinghe, S.Hong Lee
-

NB: The proposed direction of causation refers to the causal directions of GxE interactions that are the primary focus of the researcher's interest, while the reverse direction of causation examines the opposite directions of GxE interactions to test its validity.
   
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
The quantitative outcome "qp_dis_phen", quantitative exposure "qp_dis_cov" and genotype data for the discovery dataset, while in the target dataset,
The quantitative outcome is "qp_tar_phen", quantitative exposure is "qp_tar_cov". 
In the discovery dataset used to compute GWEIS using b and GWAS for Exposure using d then use e to compute PRS then g in the proposed causal directions. 
finally, the result from this direction is displayed as: 

~~~
 print(g)

Call:
lm(formula = model_formula, data = regression_data)

Residuals:
    Min      1Q  Median      3Q     Max
-9.3085 -3.1413 -0.7675  2.5247 18.5288

Coefficients:
                  Estimate Std. Error t value Pr(>|t|)
(Intercept)     -4.344e+02  3.452e+02  -1.258 0.208690
Add_PRS         -7.675e-02  1.754e-01  -0.438 0.661813
Int_PRS         -1.037e-01  1.658e-01  -0.625 0.532041
Covariate_Pheno  4.423e-02  1.662e-01   0.266 0.790249
Conf_1           1.032e-02  3.545e-02   0.291 0.770957
Conf_2          -9.487e-02  5.904e-02  -1.607 0.108460
Conf_3          -2.649e-02  2.113e-02  -1.254 0.210258
Conf_4          -1.336e-01  3.441e-01  -0.388 0.697807
Conf_5          -1.779e-01  1.029e-01  -1.729 0.084268 .
Conf_6           7.062e-02  1.258e-01   0.561 0.574716
Conf_7           2.706e-01  1.181e-01   2.292 0.022183 *
Conf_8           1.028e-01  8.872e-02   1.159 0.246768
Conf_9          -4.297e-02  3.793e-02  -1.133 0.257573
Conf_10         -4.065e-01  1.133e-01  -3.587 0.000355 ***
Conf_11          6.330e-02  1.042e-01   0.607 0.543856
Conf_12          2.339e-01  1.046e-01   2.237 0.025583 *
Conf_13         -5.914e-02  3.751e-02  -1.576 0.115340
Conf_14          1.605e-01  8.152e-02   1.969 0.049260 *
Conf_15          4.179e-02  3.136e-02   1.333 0.183079
Conf_16          5.888e-04  1.796e-04   3.278 0.001092 **
Int_PRS:Cov_PRS  6.021e-02  1.539e-01   0.391 0.695799
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 4.646 on 779 degrees of freedom
Multiple R-squared:  0.06594,   Adjusted R-squared:  0.04195
F-statistic:  2.75 on 20 and 779 DF,  p-value: 6.311e-05
~~~

**reverse causal directions**
The quantitative outcome "qr_dis_phen", quantitative exposure "qr_dis_cov" and genotype data for the target dataset, while in the target dataset,
The quantitative outcome is "qr_tar_phen", quantitative exposure is "qr_tar_cov". 
In the target dataset used to compute GWEIS using b and GWAS for Exposure using d then use e to compute PRS then g in the proposed causal directions. 
finally, the result from this direction is displayed as: 

GWEIS
~~~
 b <- q_gweis(plink_path, dis_snp, qr_dis_cov, qr_dis_phen)
~~~
GWAS
~~~
 d <- q_gwas(plink_path, dis_snp, qr_dis_phen)
~~~
PRS values
 ~~~
 e <- prs_scores(plink_path, tar_snp)
 ~~~
Compute GWEIS
   ~~~
g <- gcim_q(qr_tar_cov, qr_tar_phen, Add_PRS, Int_PRS, Cov_PRS, confounders)
   ~~~

~~~
 print(g)

Call:
lm(formula = model_formula, data = regression_data)

Residuals:
      Min        1Q    Median        3Q       Max
-0.192639 -0.065748 -0.000373  0.068433  0.236352

Coefficients:
                  Estimate Std. Error t value Pr(>|t|)
(Intercept)     -3.302e+00  6.910e+00  -0.478   0.6329
Add_PRS          1.076e-03  3.430e-03   0.314   0.7537
Int_PRS          5.608e-04  3.482e-03   0.161   0.8721
Covariate_Pheno -6.186e-04  6.686e-04  -0.925   0.3551
Conf_1          -6.614e-06  7.073e-04  -0.009   0.9925
Conf_2          -7.223e-04  1.184e-03  -0.610   0.5419
Conf_3           5.978e-04  4.231e-04   1.413   0.1581
Conf_4           2.199e-03  6.884e-03   0.319   0.7495
Conf_5           2.123e-03  2.070e-03   1.026   0.3054
Conf_6           1.981e-03  2.519e-03   0.787   0.4318
Conf_7           1.945e-03  2.310e-03   0.842   0.4002
Conf_8          -1.726e-03  1.778e-03  -0.971   0.3320
Conf_9           4.812e-04  7.551e-04   0.637   0.5241
Conf_10         -1.288e-02  2.271e-03  -5.671 2.01e-08 ***
Conf_11         -1.108e-03  2.110e-03  -0.525   0.5998
Conf_12          2.029e-03  2.089e-03   0.971   0.3319
Conf_13         -7.882e-04  7.500e-04  -1.051   0.2936
Conf_14          3.065e-04  1.636e-03   0.187   0.8515
Conf_15          3.790e-04  6.277e-04   0.604   0.5461
Conf_16          6.673e-06  3.603e-06   1.852   0.0644 .
Int_PRS:Cov_PRS -3.898e-04  2.507e-03  -0.156   0.8764
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.09297 on 779 degrees of freedom
Multiple R-squared:  0.0498,    Adjusted R-squared:  0.0254
F-statistic: 2.041 on 20 and 779 DF,  p-value: 0.004647
~~~~



