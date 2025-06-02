---
# <h1 align="center">GCIM</h1>
---

The genetic causality inference model(GCIM) is a statistical method for detecting the direction of causation in GxE interaction studies. 

- 
 Authors: Zinabu Fentaw, Dovini Jayasinghe, S.Hong Lee
-

NB: The proposed direction of causation refers to the causal directions of GxE interactions that are the primary focus of the researcher's interest, while the reverse direction of causation examines the opposite directions of GxE interactions to test its proper causal directions.
   
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
1. Download plink2 software from the [Plink](https://www.cog-genomics.org/plink/2.0/) website and then specify the executable Plink file path.
   
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
bp_dis_cov <- paste0(inst_path, "/bp_dis_cov.txt")
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

**1**  **Quantitative outcome and quantitative exposure** in testing in two different causal directions.

**1.1.** **Proposed causal directions**
The quantitative outcome "qp_dis_phen", quantitative exposure "qp_dis_cov" and genotype data for the discovery dataset, while in the target dataset,
The quantitative outcome is "qp_tar_phen", quantitative exposure is "qp_tar_cov". 
In the discovery dataset used to compute GWEIS using *b* and GWAS for exposure using "d" then use "e" to compute PRS then *g* in the proposed causal directions. 
Finally, the result from this direction is displayed as: 

~~~
 print(g)
~~~

~~~
Call:
lm(formula = model_formula, data = regression_data)

Residuals:
    Min      1Q  Median      3Q     Max
-9.1782 -3.1749 -0.7568  2.3319 18.3563

Coefficients:
                  Estimate Std. Error t value Pr(>|t|)
(Intercept)     -4.274e+02  3.454e+02  -1.237 0.216358
Add_PRS         -3.082e-02  1.690e-01  -0.182 0.855398
Int_PRS          9.378e-02  1.666e-01   0.563 0.573643
Covariate_Pheno  4.888e-02  1.661e-01   0.294 0.768668
Conf_1           1.157e-02  3.545e-02   0.326 0.744260
Conf_2          -9.765e-02  5.924e-02  -1.648 0.099662 .
Conf_3          -2.658e-02  2.115e-02  -1.257 0.209209
Conf_4          -1.277e-01  3.438e-01  -0.371 0.710402
Conf_5          -1.729e-01  1.029e-01  -1.680 0.093297 .
Conf_6           6.221e-02  1.261e-01   0.493 0.621800
Conf_7           2.755e-01  1.177e-01   2.340 0.019517 *
Conf_8           1.086e-01  8.841e-02   1.228 0.219734
Conf_9          -4.362e-02  3.785e-02  -1.152 0.249510
Conf_10         -4.045e-01  1.134e-01  -3.566 0.000385 ***
Conf_11          6.122e-02  1.042e-01   0.588 0.556931
Conf_12          2.353e-01  1.046e-01   2.249 0.024769 *
Conf_13         -6.056e-02  3.752e-02  -1.614 0.106930
Conf_14          1.642e-01  8.164e-02   2.012 0.044607 *
Conf_15          4.116e-02  3.138e-02   1.312 0.190005
Conf_16          5.905e-04  1.807e-04   3.269 0.001128 **
Int_PRS:Cov_PRS -6.615e-02  1.561e-01  -0.424 0.671790
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 4.646 on 779 degrees of freedom
Multiple R-squared:  0.06582,   Adjusted R-squared:  0.04184
F-statistic: 2.744 on 20 and 779 DF,  p-value: 6.52e-05
~~~

**1.2.** **reverse causal directions**
The quantitative outcome "qr_dis_cov", quantitative exposure "qr_dis_phen" and genotype data for the target dataset, while in the target dataset,
The quantitative outcome is "qr_tar_cov", quantitative exposure is "qr_tar_phen". 
In the target dataset used to compute GWEIS using *b* and GWAS for exposure using *d* then use *e* to compute PRS then *g* in the reverse causal directions. 
Finally, the result from this direction is displayed as: 

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
~~~

~~~
Call:
lm(formula = model_formula, data = regression_data)

Residuals:
      Min        1Q    Median        3Q       Max
-0.192060 -0.065682 -0.001634  0.068285  0.238632

Coefficients:
                  Estimate Std. Error t value Pr(>|t|)
(Intercept)     -3.334e+00  6.890e+00  -0.484   0.6286
Add_PRS         -1.610e-03  3.358e-03  -0.479   0.6318
Int_PRS          9.244e-04  3.397e-03   0.272   0.7856
Covariate_Pheno -5.807e-04  6.688e-04  -0.868   0.3856
Conf_1          -5.734e-06  7.081e-04  -0.008   0.9935
Conf_2          -8.280e-04  1.189e-03  -0.697   0.4863
Conf_3           5.926e-04  4.229e-04   1.401   0.1615
Conf_4           2.038e-03  6.884e-03   0.296   0.7673
Conf_5           2.113e-03  2.069e-03   1.021   0.3075
Conf_6           1.902e-03  2.522e-03   0.754   0.4509
Conf_7           1.984e-03  2.310e-03   0.859   0.3906
Conf_8          -1.811e-03  1.771e-03  -1.023   0.3068
Conf_9           5.146e-04  7.539e-04   0.683   0.4950
Conf_10         -1.280e-02  2.266e-03  -5.650 2.24e-08 ***
Conf_11         -1.159e-03  2.109e-03  -0.549   0.5829
Conf_12          1.927e-03  2.089e-03   0.922   0.3566
Conf_13         -7.947e-04  7.498e-04  -1.060   0.2895
Conf_14          2.038e-04  1.637e-03   0.125   0.9009
Conf_15          3.818e-04  6.258e-04   0.610   0.5420
Conf_16          6.807e-06  3.612e-06   1.885   0.0598 .
Int_PRS:Cov_PRS -2.012e-03  2.825e-03  -0.712   0.4766
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.09292 on 779 degrees of freedom
Multiple R-squared:  0.05064,   Adjusted R-squared:  0.02627
F-statistic: 2.078 on 20 and 779 DF,  p-value: 0.003784
~~~~

**2**. **Quantitative outcome and Binary exposure** in testing in two different causal directions.

**2.1.** **Proposed causal directions**
The quantitative outcome "qp_dis_phen", binary exposure "bp_dis_cov" and genotype data for the discovery dataset, while in the target dataset,
The quantitative outcome is "qp_tar_phen", and the binary exposure is "bp_tar_cov". 
In the discovery dataset used to compute GWEIS using *"b"* and GWAS for exposure using *"c"* then use *"e"* to compute PRS then *"g"* in the proposed causal directions. 
finally, the result from this direction is displayed as: 
The analysis will be in the following steps. 

GWEIS
~~~
 b <- q_gweis(plink_path, dis_snp, qp_dis_phen, bp_dis_cov)
~~~
GWAS
~~~
 c <- b_gwas(plink_path, dis_snp, bp_dis_cov)
~~~
PRS values
 ~~~
 e <- prs_scores(plink_path, tar_snp)
 ~~~
Compute GWEIS
~~~
g <- gcim_q(qp_tar_phen, bp_tar_cov, Add_PRS, Int_PRS, Cov_PRS, confounders)
~~~

~~~
print(g)
~~~

~~~
Call:
lm(formula = model_formula, data = regression_data)

Residuals:
    Min      1Q  Median      3Q     Max
-9.2207 -3.1875 -0.7733  2.4365 18.4096

Coefficients:
                  Estimate Std. Error t value Pr(>|t|)
(Intercept)     -4.259e+02  3.449e+02  -1.235 0.217304
Add_PRS         -3.246e-02  1.683e-01  -0.193 0.847155
Int_PRS         -6.700e-02  1.744e-01  -0.384 0.700945
Covariate_Pheno  2.240e-01  5.919e-01   0.378 0.705182
Conf_1           1.175e-02  3.544e-02   0.331 0.740393
Conf_2          -9.467e-02  5.978e-02  -1.584 0.113656
Conf_3          -2.559e-02  2.116e-02  -1.209 0.226959
Conf_4          -1.158e-01  3.429e-01  -0.338 0.735749
Conf_5          -1.704e-01  1.035e-01  -1.646 0.100212
Conf_6           6.825e-02  1.260e-01   0.542 0.588211
Conf_7           2.782e-01  1.156e-01   2.407 0.016315 *
Conf_8           1.113e-01  8.879e-02   1.254 0.210388
Conf_9          -4.563e-02  3.830e-02  -1.191 0.233820
Conf_10         -4.021e-01  1.141e-01  -3.525 0.000448 ***
Conf_11          5.935e-02  1.041e-01   0.570 0.568843
Conf_12          2.323e-01  1.047e-01   2.219 0.026747 *
Conf_13         -6.101e-02  3.790e-02  -1.610 0.107821
Conf_14          1.595e-01  8.141e-02   1.959 0.050445.
Conf_15          4.102e-02  3.133e-02   1.309 0.190875
Conf_16          5.956e-04  1.808e-04   3.295 0.001029 **
Int_PRS:Cov_PRS  2.720e-02  1.312e-01   0.207 0.835835
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 4.647 on 779 degrees of freedom
Multiple R-squared:  0.06546,   Adjusted R-squared:  0.04147
F-statistic: 2.728 on 20 and 779 DF,  p-value: 7.232e-05
~~~

**2.2.** **reverse causal directions**
The binary outcome "br_dis_cov", quantitative exposure "qr_dis_phen" and genotype data for the discovery dataset, while in the target dataset,
the binary outcome is "br_tar_cov", and quantitative exposure is "qr_tar_phen". 
In the discovery dataset used to compute GWEIS using *a* and GWAS for exposure using *d* then use *e* to compute PRS then *f* in the reverse causal directions. 
Finally, the result from this direction is displayed as: 

GWEIS

~~~
 a <- b_gweis(plink_path, dis_snp, br_dis_cov, qr_dis_phen)
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
f <- gcim_b(br_tar_cov, qr_tar_phen, Add_PRS, Int_PRS, Cov_PRS, confounders)
~~~

~~~
print(f)
~~~

~~~
Call:
glm(formula = model_formula, family = binomial(), data = regression_data)

Deviance Residuals:
    Min       1Q   Median       3Q      Max
-1.3233  -0.4885  -0.3860  -0.2777   2.7542

Coefficients:
                  Estimate Std. Error z value Pr(>|z|)
(Intercept)     -5.867e+02  2.563e+02  -2.290   0.0220 *
Add_PRS         -7.340e-02  1.312e-01  -0.559   0.5759
Int_PRS          2.461e-01  1.305e-01   1.885   0.0594 .
Covariate_Pheno  3.773e-03  2.254e-02   0.167   0.8671
Conf_1          -3.237e-02  2.634e-02  -1.229   0.2191
Conf_2           3.452e-02  4.465e-02   0.773   0.4394
Conf_3           9.559e-03  1.604e-02   0.596   0.5511
Conf_4           1.999e-01  2.553e-01   0.783   0.4336
Conf_5          -1.693e-01  7.943e-02  -2.131   0.0331 *
Conf_6          -1.470e-02  9.835e-02  -0.149   0.8812
Conf_7          -1.050e-01  8.957e-02  -1.172   0.2412
Conf_8          -1.739e-02  7.243e-02  -0.240   0.8102
Conf_9          -3.244e-02  2.986e-02  -1.086   0.2773
Conf_10         -2.584e-02  8.463e-02  -0.305   0.7601
Conf_11          2.534e-02  8.102e-02   0.313   0.7544
Conf_12         -1.081e-01  7.822e-02  -1.383   0.1668
Conf_13         -3.873e-02  2.851e-02  -1.359   0.1743
Conf_14         -2.551e-01  6.136e-02  -4.158 3.22e-05 ***
Conf_15          5.283e-02  2.326e-02   2.271   0.0231 *
Conf_16          9.654e-05  1.399e-04   0.690   0.4901
Int_PRS:Cov_PRS -1.035e-02  1.106e-01  -0.094   0.9254
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

(Dispersion parameter for binomial family taken to be 1)

    Null deviance: 520.13  on 799  degrees of freedom
Residual deviance: 477.60  on 779  degrees of freedom
AIC: 519.6

Number of Fisher Scoring iterations: 5
~~~

**3**  **Binary outcome and quantitative exposure** in testing in two different causal directions. 

**3.1.** **Proposed causal directions**
The Binary outcome "bp_dis_phen", quantitative exposure "qp_dis_cov" and genotype data for the discovery dataset, while in the target dataset,
the binary outcome is "bp_tar_phen", and quantitative exposure is "qp_tar_cov". 
In the discovery dataset used to compute GWEIS using *a* and GWAS for exposure using *d* then use *e* to compute PRS then *f* in the proposed causal directions. 
Finally, the result from this direction is displayed as: 

GWEIS

~~~
 b <- b_gweis(plink_path, dis_snp, bp_dis_phen, qp_dis_cov)
~~~

GWAS

~~~
 d <- q_gwas(plink_path, dis_snp, qp_dis_phen)
~~~

PRS values

 ~~~
 e <- prs_scores(plink_path, tar_snp)
 ~~~

Compute regression

 ~~~
f <- gcim_b(bp_tar_phen, qp_tar_cov, Add_PRS, Int_PRS, Cov_PRS, confounders)
~~~

~~~
print(f)
~~~

~~~
Call:
glm(formula = model_formula, family = binomial(), data = regression_data)

Deviance Residuals:
    Min       1Q   Median       3Q      Max
-0.9335  -0.4161  -0.2977  -0.1791   2.6932
Coefficients:
                  Estimate Std. Error z value Pr(>|z|)
(Intercept)      1.407e+02  3.079e+02   0.457  0.64762
Add_PRS         -8.051e-03  1.140e+00  -0.007  0.99436
Int_PRS         -1.336e-01  1.134e+00  -0.118  0.90617
Covariate_Pheno  4.390e-01  1.484e-01   2.959  0.00309 **
Conf_1          -2.408e-02  3.152e-02  -0.764  0.44499
Conf_2           7.504e-02  5.035e-02   1.490  0.13611
Conf_3           7.859e-03  1.930e-02   0.407  0.68393
Conf_4          -1.259e-01  3.081e-01  -0.409  0.68279
Conf_5           3.426e-02  9.506e-02   0.360  0.71853
Conf_6           2.083e-01  1.068e-01   1.951  0.05109 .
Conf_7          -3.741e-01  1.172e-01  -3.192  0.00141 **
Conf_8          -1.351e-01  8.041e-02  -1.680  0.09292 .
Conf_9           4.360e-02  3.341e-02   1.305  0.19198
Conf_10          1.427e-01  9.655e-02   1.478  0.13954
Conf_11          2.227e-01  1.016e-01   2.191  0.02847 *
Conf_12          2.638e-01  9.851e-02   2.678  0.00741 **
Conf_13         -2.797e-02  2.982e-02  -0.938  0.34821
Conf_14          2.101e-01  7.784e-02   2.699  0.00695 **
Conf_15         -1.312e-02  2.798e-02  -0.469  0.63925
Conf_16          4.972e-05  1.652e-04   0.301  0.76349
Int_PRS:Cov_PRS  1.604e-01  2.324e-01   0.690  0.49012

Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

(Dispersion parameter for binomial family taken to be 1)

    Null deviance: 405.82  on 799  degrees of freedom
Residual deviance: 359.71  on 779  degrees of freedom
AIC: 401.71
Number of Fisher Scoring iterations: 6
~~~

**3.2.** **Reverse causal directions**
The quantitative outcome "qr_dis_cov", binary exposure "br_dis_phen" and genotype data for the discovery dataset, while in the target dataset,
the quantitative outcome is "qr_tar_cov", and the binary exposure is "br_tar_phen". 
In the discovery dataset used to compute GWEIS using *b* and GWAS for exposure using *c* then use *e* to compute PRS then *g* in the reverse causal directions. 
finally, the result from this direction is displayed as: 
GWEIS

~~~
 b <- q_gweis(plink_path, dis_snp, qr_dis_cov, br_dis_phen)
~~~

GWAS

~~~
 c <- b_gwas(plink_path, dis_snp, br_dis_phen)
~~~

PRS values

 ~~~
 e <- prs_scores(plink_path, tar_snp)
 ~~~

Compute regression

 ~~~
g <- gcim_q(qr_tar_cov, br_tar_phen, Add_PRS, Int_PRS, Cov_PRS, confounders)
~~~

~~~
print(g)
~~~

~~~
Call:
lm(formula = model_formula, data = regression_data)

Residuals:
      Min        1Q    Median        3Q       Max
-0.193632 -0.067855  0.000162  0.069502  0.241145

Coefficients:
                  Estimate Std. Error t value Pr(>|t|)
(Intercept)     -4.975e+00  6.922e+00  -0.719   0.4725
Add_PRS         -1.585e-03  3.362e-03  -0.471   0.6375
Int_PRS         -1.576e-03  3.579e-03  -0.440   0.6598
Covariate_Pheno  1.943e-02  1.007e-02   1.929   0.0541 .
Conf_1           2.022e-05  7.065e-04   0.029   0.9772
Conf_2          -3.293e-04  1.186e-03  -0.278   0.7814
Conf_3           5.163e-04  4.222e-04   1.223   0.2217
Conf_4           2.585e-03  6.856e-03   0.377   0.7062
Conf_5           1.671e-03  2.059e-03   0.812   0.4173
Conf_6           1.952e-03  2.510e-03   0.778   0.4370
Conf_7           2.388e-03  2.316e-03   1.031   0.3028
Conf_8          -1.268e-03  1.781e-03  -0.712   0.4769
Conf_9           3.559e-04  7.540e-04   0.472   0.6371
Conf_10         -1.233e-02  2.274e-03  -5.420 7.96e-08 ***
Conf_11         -5.593e-04  2.086e-03  -0.268   0.7887
Conf_12          1.551e-03  2.093e-03   0.741   0.4591
Conf_13         -5.983e-04  7.514e-04  -0.796   0.4261
Conf_14          2.250e-05  1.627e-03   0.014   0.9890
Conf_15          5.293e-04  6.288e-04   0.842   0.4002
Conf_16          5.630e-06  3.623e-06   1.554   0.1206
Int_PRS:Cov_PRS  2.516e-05  2.707e-03   0.009   0.9926
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.09277 on 779 degrees of freedom
Multiple R-squared:  0.05386,   Adjusted R-squared:  0.02957
F-statistic: 2.217 on 20 and 779 DF,  p-value: 0.001689
~~~

**4**  **Binary outcome and Binary exposure** in testing in two different causal directions. 

**4.1.** **Proposed causal directions**
The Binary outcome "bp_dis_phen", binary exposure "bp_dis_cov" and genotype data for the discovery dataset, while in the target dataset,
the binary outcome is "bp_tar_phen", and binary exposure is "bp_tar_cov". 
In the discovery dataset used to compute GWEIS using *a* and GWAS for exposure using *c*, *e* is used to calculate PRS, and then *f* is used in the proposed causal directions. 
Finally, the result from this direction is displayed as: 

GWEIS

~~~
 b <- b_gweis(plink_path, dis_snp, bp_dis_phen, bp_dis_cov)
~~~

GWAS

~~~
 c <- b_gwas(plink_path, dis_snp, bp_dis_cov)
~~~

PRS values

 ~~~
 e <- prs_scores(plink_path, tar_snp)
 ~~~

Compute regression

 ~~~
f <- gcim_b(bp_tar_phen, bp_tar_cov, Add_PRS, Int_PRS, Cov_PRS, confounders)
~~~

~~~
print(f)
~~~

**4.2.** **reverse causal directions**
The Binary outcome "br_dis_cov", binary exposure "br_dis_phen" and genotype data for the discovery dataset, while in the target dataset,
the binary outcome is "br_tar_cov", and binary exposure is "br_tar_phen". 
In the discovery dataset used to compute GWEIS using *a* and GWAS for exposure using *c*, *e* is used to calculate PRS, and then *f* is used in the proposed causal directions. 

GWEIS

~~~
 b <- b_gweis(plink_path, dis_snp, br_dis_cov, br_dis_phen)
~~~

GWAS

~~~
 c <- b_gwas(plink_path, dis_snp, br_dis_phen)
~~~

PRS values

 ~~~
 e <- prs_scores(plink_path, tar_snp)
 ~~~

Compute regression

 ~~~
f <- gcim_b(br_tar_cov, br_tar_phen, Add_PRS, Int_PRS, Cov_PRS, confounders)
~~~

~~~
print(f)
~~~
Finally, the result from this direction is displayed as: 
~~~
Call:
glm(formula = model_formula, family = binomial(), data = regression_data)

Deviance Residuals:
    Min       1Q   Median       3Q      Max
-1.3883  -0.4707  -0.3543  -0.2403   2.4512

Coefficients:
                  Estimate Std. Error z value Pr(>|z|)
(Intercept)     -7.458e+02  2.672e+02  -2.791  0.00526 **
Add_PRS         -8.492e-02  1.334e-01  -0.637  0.52438
Int_PRS          3.214e-02  1.358e-01   0.237  0.81299
Covariate_Pheno  1.418e+00  2.990e-01   4.743 2.10e-06 ***
Conf_1          -3.053e-02  2.697e-02  -1.132  0.25767
Conf_2           5.925e-02  4.572e-02   1.296  0.19503
Conf_3           8.024e-03  1.616e-02   0.497  0.61946
Conf_4           1.779e-01  2.570e-01   0.692  0.48868
Conf_5          -2.054e-01  8.117e-02  -2.530  0.01141 *
Conf_6          -1.404e-02  1.030e-01  -0.136  0.89166
Conf_7          -4.385e-02  8.766e-02  -0.500  0.61690
Conf_8           3.308e-02  7.581e-02   0.436  0.66262
Conf_9          -4.693e-02  3.036e-02  -1.546  0.12214
Conf_10          3.002e-02  8.769e-02   0.342  0.73210
Conf_11          5.073e-02  7.851e-02   0.646  0.51819
Conf_12         -1.447e-01  7.970e-02  -1.816  0.06935 .
Conf_13         -2.721e-02  2.995e-02  -0.908  0.36365
Conf_14         -2.839e-01  6.298e-02  -4.507 6.57e-06 ***
Conf_15          6.724e-02  2.427e-02   2.771  0.00559 **
Conf_16          3.650e-05  1.422e-04   0.257  0.79739
Int_PRS:Cov_PRS  1.219e-01  9.980e-02   1.221  0.22194
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

(Dispersion parameter for binomial family taken to be 1)

    Null deviance: 520.13  on 799  degrees of freedom
Residual deviance: 460.29  on 779  degrees of freedom
AIC: 502.29

Number of Fisher Scoring iterations: 6
~~~
