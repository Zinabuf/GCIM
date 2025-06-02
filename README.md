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
inst_path <- system.file(package = "GCIM") 
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
-0.192525 -0.063515  0.000921  0.067702  0.237396

Coefficients:
                  Estimate Std. Error t value Pr(>|t|)
(Intercept)     -3.531e+00  6.890e+00  -0.513   0.6084
Add_PRS         -1.811e-03  2.650e-02  -0.068   0.9455
Int_PRS         -3.138e-03  2.636e-02  -0.119   0.9052
Covariate_Pheno  1.145e-02  1.315e-02   0.870   0.3845
Conf_1           5.196e-05  7.076e-04   0.073   0.9415
Conf_2          -4.941e-04  1.178e-03  -0.419   0.6751
Conf_3           5.880e-04  4.220e-04   1.393   0.1639
Conf_4           2.350e-03  6.851e-03   0.343   0.7317
Conf_5           1.940e-03  2.056e-03   0.944   0.3457
Conf_6           2.055e-03  2.521e-03   0.815   0.4151
Conf_7           2.296e-03  2.316e-03   0.992   0.3217
Conf_8          -1.437e-03  1.775e-03  -0.809   0.4187
Conf_9           3.786e-04  7.541e-04   0.502   0.6158
Conf_10         -1.300e-02  2.266e-03  -5.738 1.37e-08 ***
Conf_11         -1.065e-03  2.086e-03  -0.510   0.6100
Conf_12          1.744e-03  2.093e-03   0.833   0.4050
Conf_13         -7.306e-04  7.489e-04  -0.976   0.3296
Conf_14         -2.540e-05  1.636e-03  -0.016   0.9876
Conf_15          3.981e-04  6.259e-04   0.636   0.5249
Conf_16          6.099e-06  3.593e-06   1.697   0.0901 .
Int_PRS:Cov_PRS  4.631e-03  3.574e-03   1.296   0.1955
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.09284 on 779 degrees of freedom
Multiple R-squared:  0.05242,   Adjusted R-squared:  0.02809
F-statistic: 2.155 on 20 and 779 DF,  p-value: 0.002433

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
