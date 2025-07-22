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
install.packages(GxEprs)
install.packages(GCIM) 
~~~

Load the library

~~~
library(GxEprs)
~~~

**Data preparations**
The dataset is divided into **discovery (80%)** and **target (20%)** datasets, ensuring consistency across genetic, outcome, exposure, and confounder data. Genetic data, stored in **PLINK binary format** (`.bed`, `.bim`, `.fam`). The outcome variable should include **FID, IID, and phenotype values**, where case-control labels follow PLINK conventions: **1 = Control, 2 = Case** in the discovery dataset, and **0 = Control, 1 = Case** in the target dataset. Exposure and confounder variables are formatted into at least **m columns** (**FID, IID, exposure, exposure_sqaure, confounder1-m **) and partitioned in the same proportions. For a binary outcome, the fourth column (exposure_square) should contain constant values (zero variance), such as 2, 2, 2, ..., 2,  for all observations. This structured approach ensures compatibility across all data types, thereby maintaining accurate alignment to estimate GxE interactions.


A Guide for GCIM analyses

GCIM analyses use PLink2 to analyze discovery data, and the package is compatible with the Linux operating system. 
1. Download plink2 software from the [Plink](https://www.cog-genomics.org/plink/2.0/) website and then specify the executable Plink file path.
   
~~~
plink_path <- "<plink_path>/plink2"
~~~
**Discovery dataset**: In the discovery dataset, please use the data input from the GxEprs data input format, then construct the PRS of the exposure using [Plink](https://www.cog-genomics.org/plink/2.0/) based on the GWAS summary statistics for the PRS of the target samples.
**Target dataset**: The target dataset for the model is also similar in data format, except, the use of PRS of the exposure variable(PRS of E) rather than the use of entire exposure values, and also includes the constant values for the fourth column, which is the square of the third column in the GxEprs model. 

**Example data**
To conduct a GCIM analysis, we must use the same data format as GxEprs, especially in the discovery dataset, but there is a slight difference in the target dataset; no need for the square of the exposure variables. Here is an example of analysis using the R script in the exaple directories. 
1. Quantitative outcome
1.1. quantitative exposure
   
 ~~~
# Load required libraries
library(GxEprs)
library(GCIM)
# Set plink path
plink_path <- "/data/alh-admzw/plink2"

# For quantitative traits, use corresponding functions
a <- GWAS_quantitative(plink_path, "DummyData", "Qexp_disc.txt", "Qcov_disc.txt")
b <- GWEIS_quantitative(plink_path, "DummyData", "Qphe_discovery.txt", "Qcov_discovery.txt")

# Extract and compute PRS
trd <- a[c("ID", "A1", "BETA")]
add <- b[c("ID", "A1", "ADD_BETA")]
gxe <- b[c("ID", "A1", "INTERACTION_BETA")]

q <- PRS_quantitative(plink_path, "DummyData", summary_input = add)
r <- PRS_quantitative(plink_path, "DummyData", summary_input = gxe) 
p <- PRS_quantitative(plink_path, "DummyData", summary_input = trd)

# Run GCIM for continuous outcome
result_qnt <- gcim_q("Qphe_target.txt", "Bexp_target.txt",
  Add_PRS,
  Int_PRS,
  Cov_PRS
)
print(result_qnt$model_summary)
~~~

 Result displayed using this approach

~~~
Coefficients:
                  Estimate Std. Error t value Pr(>|t|)
(Intercept)      2.777e+01  3.599e+00   7.715 7.54e-13 ***
Add_PRS         -4.979e+01  8.691e+01  -0.573   0.5674
Int_PRS         -3.565e+01  7.775e+01  -0.459   0.6471
Covariate_Pheno  3.492e-01  6.994e-01   0.499   0.6182
Conf_1          -7.260e-02  1.286e-01  -0.565   0.5730
Conf_2          -1.286e-02  4.087e-02  -0.315   0.7534
Conf_3          -1.298e-02  2.344e-01  -0.055   0.9559
Conf_4          -1.162e-02  2.569e-01  -0.045   0.9640
Conf_5           1.188e-01  2.266e-01   0.524   0.6009
Conf_6          -3.234e-02  1.666e-01  -0.194   0.8463
Conf_7          -1.216e-01  8.028e-02  -1.515   0.1315
Conf_8           5.611e-01  2.174e-01   2.580   0.0107 *
Conf_9          -6.929e-02  1.884e-01  -0.368   0.7135
Conf_10         -2.249e-01  1.933e-01  -1.163   0.2463
Conf_11         -8.115e-02  8.315e-02  -0.976   0.3304
Conf_12         -1.349e-01  1.760e-01  -0.766   0.4444
Int_PRS:Cov_PRS -1.125e+04  9.054e+04  -0.124   **0.9012**
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 4.735 on 183 degrees of freedom
Multiple R-squared:  0.08653,   Adjusted R-squared:  0.006659
F-statistic: 1.083 on 16 and 183 DF,  p-value: 0.3737
~~~
1.2. Binary exposure variables

~~~
Coefficients:
                  Estimate Std. Error t value Pr(>|t|)
(Intercept)      2.777e+01  3.599e+00   7.715 7.54e-13 ***
Add_PRS         -4.979e+01  8.691e+01  -0.573   0.5674
Int_PRS         -3.565e+01  7.775e+01  -0.459   0.6471
Covariate_Pheno  3.492e-01  6.994e-01   0.499   0.6182
Conf_1          -7.260e-02  1.286e-01  -0.565   0.5730
Conf_2          -1.286e-02  4.087e-02  -0.315   0.7534
Conf_3          -1.298e-02  2.344e-01  -0.055   0.9559
Conf_4          -1.162e-02  2.569e-01  -0.045   0.9640
Conf_5           1.188e-01  2.266e-01   0.524   0.6009
Conf_6          -3.234e-02  1.666e-01  -0.194   0.8463
Conf_7          -1.216e-01  8.028e-02  -1.515   0.1315
Conf_8           5.611e-01  2.174e-01   2.580   0.0107 *
Conf_9          -6.929e-02  1.884e-01  -0.368   0.7135
Conf_10         -2.249e-01  1.933e-01  -1.163   0.2463
Conf_11         -8.115e-02  8.315e-02  -0.976   0.3304
Conf_12         -1.349e-01  1.760e-01  -0.766   0.4444
Int_PRS:Cov_PRS -1.125e+04  9.054e+04  -0.124   **0.9012**
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 4.735 on 183 degrees of freedom
Multiple R-squared:  0.08653,   Adjusted R-squared:  0.006659
F-statistic: 1.083 on 16 and 183 DF,  p-value: 0.3737
~~~

 2. Binary outcome
    2.1. Quantitative exposure

    ~~~
    Coefficients:
                  Estimate Std. Error z value Pr(>|z|)
(Intercept)      3.499e+00  3.429e+00   1.020   0.3075
Add_PRS         -6.144e+01  4.595e+01  -1.337   0.1811
Int_PRS         -9.495e+01  6.858e+01  -1.384   0.1662
Covariate_Pheno  5.700e-02  4.689e-01   0.122   0.9032
Conf_1          -9.505e-03  2.292e-01  -0.041   0.9669
Conf_2           7.363e-02  1.018e-01   0.723   0.4696
Conf_3          -3.781e-02  3.926e-02  -0.963   0.3355
Conf_4           2.544e-01  2.316e-01   1.099   0.2720
Conf_5           1.308e-01  2.366e-01   0.553   0.5803
Conf_6           1.027e-01  1.983e-01   0.518   0.6044
Conf_7          -7.876e-02  1.587e-01  -0.496   0.6196
Conf_8           5.808e-02  8.538e-02   0.680   0.4963
Conf_9           9.249e-02  2.134e-01   0.433   0.6648
Conf_10          2.398e-01  1.701e-01   1.410   0.1587
Conf_11          5.093e-02  1.688e-01   0.302   0.7628
Conf_12          2.992e-03  8.217e-02   0.036   0.9710
Conf_13          6.016e-02  1.583e-01   0.380   0.7039
Conf_14          1.159e+00  6.839e-01   1.694   0.0902 .
Conf_15         -1.936e-01  7.891e-02  -2.453   0.0142 *
Int_PRS:Cov_PRS -2.775e+02  6.966e+03  -0.040   **0.9682**
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

(Dispersion parameter for binomial family taken to be 1)

    Null deviance: 106.554  on 199  degrees of freedom
Residual deviance:  85.527  on 180  degrees of freedom
AIC: 125.53
~~~

2.2. Binary exposure variables

~~~
Coefficients:
                  Estimate Std. Error z value Pr(>|z|)
(Intercept)     -6.423e-01  2.832e+00  -0.227    0.821
Add_PRS         -6.132e+01  4.482e+01  -1.368    0.171
Int_PRS         -7.663e+01  8.895e+01  -0.862    0.389
Covariate_Pheno  9.287e-01  6.087e-01   1.526    0.127
Conf_1           8.027e-02  9.956e-02   0.806    0.420
Conf_2          -1.717e-02  3.633e-02  -0.473    0.636
Conf_3           2.175e-01  2.060e-01   1.056    0.291
Conf_4           1.455e-01  2.213e-01   0.658    0.511
Conf_5           7.353e-02  1.952e-01   0.377    0.706
Conf_6          -5.747e-02  1.476e-01  -0.389    0.697
Conf_7           3.685e-02  7.686e-02   0.479    0.632
Conf_8          -1.000e-02  1.863e-01  -0.054    0.957
Conf_9           2.306e-01  1.547e-01   1.490    0.136
Conf_10          1.131e-01  1.476e-01   0.766    0.444
Conf_11          1.993e-03  7.720e-02   0.026    0.979
Conf_12         -9.569e-03  1.475e-01  -0.065    0.948
Int_PRS:Cov_PRS  1.830e+03  4.816e+03   0.380    **0.704**

(Dispersion parameter for binomial family taken to be 1)

    Null deviance: 106.554  on 199  degrees of freedom
Residual deviance:  93.325  on 183  degrees of freedom
AIC: 127.32
~~~    

