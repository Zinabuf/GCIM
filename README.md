---
# <h1 align="center">GCIM</h1>
---

The **genetic causality inference model (GCIM)** is a statistical method for detecting the causal direction of Genotype-by-environment(GxE) interaction studies 



 ## Authors: Zinabu Fentaw, Dovini Jayasinghe, S.Hong Lee


_<div align="justify">GCIM is a novel statistical method that extends beyond traditional PRS×E approaches.
    It systematically evaluates both the proposed and reverse causal directions, rather than relying
    on the assumptions of prior causal directions. By explicitly testing both directions, GCIM offers researchers
    data-driven insight into the likely causal directions of GxE interactions.</div>_

_<div align="justify">NB: The proposed direction of causation refers to the hypothesized G×E interaction in which the exposure affects the outcome, aligning with the researcher’s primary interest. In contrast, the reverse direction test evaluates the opposite relationship, switching the roles of exposure and outcome to assess whether the assumed causal direction is valid. Once the data has been properly prepared, begin by testing the causal direction based on the proposed direction of interest, that is the direction specified by the researcher. After this primary analysis, assess the reverse causal direction by switching the roles of the exposure and outcome variables in both the discovery and target datasets. This involves treating the original outcome as the exposure and the original exposure as the outcome, ensuring consistency in data structure and formatting across both analyses.</div>_
   
## I. Package installation 
From GitHub 

~~~
library(devtools)
install_github("DoviniJ/GxEprs")
install_github("Zinabuf/GCIM")
~~~

Or the CRAN version via

~~~
install.packages(GxEprs)
~~~

## II. Load the library

~~~
library(GxEprs)
library(GCIM)
~~~

## III. Data Preparation

To ensure consistent and reliable estimation of G×E, the dataset should be split into two subsets: **a discovery dataset (80%)** and **a target dataset (20%)**. This split should maintain consistency across genetic data, outcome variables, exposure variables, and potential confounders. The genetic data must be in **PLINK binary format**, comprising three files: `.bed`, `.bim`, and `.fam`. The **outcome file** should include `FID`, `IID`, and the outcome variable. For binary outcomes, follow standard coding conventions: use **PLINK’s default coding (1 = Control, 2 = Case)** for the **discovery dataset**, and use **binary coding (0 = Control, 1 = Case)** for the **target dataset**. The **exposure and confounder file for the discovery dataset** should contain at least `m` columns (Minimum of 3) with the following format: `FID`, `IID`, `exposure`, `exposure_squared`, `confounder_1`, ..., `confounder_m`. Note that **exposure-squared is not required** in the **target dataset**. Additionally, for conducting a genome-wide association study (GWAS) of the exposure variable in the discovery dataset (for use in computing the Polygenic Risk Score of the exposure), the exposure data should also be formatted separately as:

* `FID`, `IID`, and **exposure values**
* A covariate file in the format: `FID`, `IID`, `confounder_1`, `exposure_square`, `confounder_2`, ..., `confounder_m`

This standardized format ensures that all variables, such as genetic, exposure, outcome, and confounders are properly aligned across discovery and target datasets, facilitating valid and replicable G×E interaction analysis.
All GWAS, GWEIS, and polygenic risk score (PRS) construction steps are performed using the [GxEprs](https://github.com/DoviniJ/GxEprs) R package 
, while the regression analyses for both binary and quantitative outcomes are conducted using the GCIM R package. 


**Example data**
<div align="justify">To conduct a GCIM analysis, the input data must follow the same format as required for GxEprs, particularly in the discovery dataset. The only distinction arises in the target dataset, where the squared term of the exposure variable is not necessary. I've included an example analysis using the accompanying R script below to show the implementation.
 ##Note: The squared term of the exposure variable has no role in modeling only for the dataframe, so it can be any value.</div>  

Data structure
<div align="justify">After splitting the data into two independent subsets, we designated one as the discovery dataset and the other as the target dataset. The discovery dataset contains all necessary inputs, including genotype data, the outcome and exposure phenotypes, as well as covariates used for adjustment. This dataset is used to conduct both GWAS (for the exposure variable) and GWEIS (for the outcome variable). The target dataset, in contrast, is reserved for detecting the direction of GxE interactions using the PRS derived from the analyses in the discovery dataset. 
 To construct PRS for the exposure variable, we first performed a GWAS on the quantitative exposure phenotype, adopting the same input data format required by the GxEprs framework. In this procedure, the exposure is treated as the outcome variable in the GWAS to obtain SNP effect estimates. For reproducibility, the exposure and covariates should be stored in a separate file. Example:</div> 'Qcov_discovery_phen.txt'.
 
~~~
  FID  IID   Exposure
1 ID_1 ID_1 -0.64402046
2 ID_2 ID_2 -0.02786981
3 ID_3 ID_3  2.12865748
4 ID_4 ID_4  2.12865748
5 ID_5 ID_5 -0.95209579
6 ID_6 ID_6 -0.02786981
~~~

The covariates used for adjustment should be provided in a separate file. Example `Qcov_discovery_cov.txt`

~~~
   FID   IID   Conf_1 Exposure_square Conf_2 Conf_3 Conf_4 Conf_5    Conf_6
1 ID_1 ID_1 -3.831420 0.4147623548 64 -14.03640 5.517420  0.0714337  5.662630
2 ID_2 ID_2  0.614044 0.0007767262 66 -10.85050 2.119980 -0.8828830 -0.441662
3 ID_3 ID_3 -0.237792 4.5311826719 55  -9.75369 3.183430 -2.0979300  6.873450
4 ID_4 ID_4  6.698660 4.5311826719 47  -9.07045 0.956878 -2.4840700  1.063590
5 ID_5 ID_5 -1.614230 0.9064863904 59 -12.93790 1.294610 -1.7997300  1.444040
6 ID_6 ID_6 -4.389270 0.0007767262 52 -11.85160 0.888978 -2.7231000  1.116810
   Conf_7    Conf_8     Conf_9     Conf_10    Conf_11   Conf_12  Conf_13 Conf_14
1  0.865562 -2.269570 -0.09658590 -2.354970  1.0588900  0.195302   0   7
2 -2.641770  2.789440  0.52458600  2.671340 -2.6372400 -0.998764   1  20
3 11.377700  2.969610 -1.11879000  0.873649  3.3552300 -4.578310   1  10
4 -3.132470  2.123200 -0.00976751  0.820582  0.0305345  1.630300   1  20
5 -6.828980 -2.967950 -2.91577000 -1.828810  7.1589200  2.109160   1  20
6 -3.646760 -0.594538 -1.75430000 -0.716014 -2.3906700  1.312950   1  10
~~~

 For GWEIS analyses, the outcome variables 
<div align="justify">To generate both the additive and interaction polygenic risk scores (PRS), we performed a genome-wide environment interaction study (GWEIS) using the GxEprs data framework. When conducting a GWEIS with a quantitative outcome, the input data must follow the same format as required for the GxEprs framework. For reproducibility, the outcome data should be organized in a dedicated file, for example:</div> 'Qphe_discovery.txt'.

 ~~~
FID   IID    Outcome
1 ID_1 ID_1 31.6534
2 ID_2 ID_2 25.5035
3 ID_3 ID_3 26.7391
4 ID_4 ID_4 25.5271
5 ID_5 ID_5 26.7165
6 ID_6 ID_6 38.8272
~~~

<div align="justify">The exposure variable and the covariate that used to addust also should look like the following data format as expressed in GxEprs 
For reproducibility, the exposure and covariate data should be organized in the following file format, for example:</div> 'Qcov_discovery.txt'.

~~~
   FID IID   Exposure    Exposure_square Conf_1 Conf_2 Conf_3 Conf_4  Conf_5 
1 ID_1 ID_1 -0.64402046 0.4147623548 -3.831420 64 -14.03640 5.517420  0.0714337
2 ID_2 ID_2 -0.02786981 0.0007767262  0.614044 66 -10.85050 2.119980 -0.8828830
3 ID_3 ID_3  2.12865748 4.5311826719 -0.237792 55  -9.75369 3.183430 -2.0979300
4 ID_4 ID_4  2.12865748 4.5311826719  6.698660 47  -9.07045 0.956878 -2.4840700
5 ID_5 ID_5 -0.95209579 0.9064863904 -1.614230 59 -12.93790 1.294610 -1.7997300
6 ID_6 ID_6 -0.02786981 0.0007767262 -4.389270 52 -11.85160 0.888978 -2.7231000
      Conf_6  Conf_7    Conf_8    Conf_9     Conf_10   Conf_11      Conf_12  Conf_13
1  5.662630  0.865562 -2.269570 -0.09658590 -2.354970  1.0588900  0.195302   0
2 -0.441662 -2.641770  2.789440  0.52458600  2.671340 -2.6372400 -0.998764   1
3  6.873450 11.377700  2.969610 -1.11879000  0.873649  3.3552300 -4.578310   1
4  1.063590 -3.132470  2.123200 -0.00976751  0.820582  0.0305345  1.630300   1
5  1.444040 -6.828980 -2.967950 -2.91577000 -1.828810  7.1589200  2.109160   1
6  1.116810 -3.646760 -0.594538 -1.75430000 -0.716014 -2.3906700  1.312950   1
  Conf_14
1   7
2  20
3  10
4  20
5  20
6  10
~~~

Target data set 
The quantitative outcome data should be organized in a separate file, for example: `Qphe_target.txt`

~~~
   FID     IID  Outcome
1 ID_801 ID_801 26.5723
2 ID_802 ID_802 20.2632
3 ID_803 ID_803 27.7365
4 ID_804 ID_804 18.7500
5 ID_805 ID_805 23.3025
6 ID_806 ID_806 24.9871
~~~

 The exposure variable and other covariates for the adjustments are for the target dataset and should be provided in a separate file, for example: `Qexp_target.tx`

 ~~~
   FID   IID      Exposure    Conf_1  Conf_2 Conf_3 Conf_4  Conf_5   Conf_6
1 ID_801 ID_801 -0.64402046 -3.826590 69 -13.8514 3.96080 -1.788050 0.0692473
2 ID_802 ID_802 -0.95209579  2.065150 60 -12.2438 4.04169 -0.905739 5.9656000
3 ID_803 ID_803 -0.02786981 -0.795863 62 -10.9195 6.91985 -2.920880 1.2601900
4 ID_804 ID_804 -0.64402046 -2.620880 67  -9.9271 4.10960 -2.354540 0.7190210
5 ID_805 ID_805  0.28020552 -3.331640 67 -11.8637 5.88272  1.072880 2.7448800
6 ID_806 ID_806  0.89635617  3.252030 61 -11.5364 5.79318 -2.311240 4.5023300
   Conf_7    Conf_8    Conf_9   Conf_10    Conf_11  Conf_12   Conf_13 Conf_14
1 -6.32556  2.853590  1.0851600 -1.303040  3.41659  1.415770   0   7
2  8.35545 -1.435760 -0.6181530  0.746918  5.11019 -0.207188   1  19
3 -5.56624 -0.552624 -0.0756095 -0.910047 -1.33896  1.726360   0   7
4 -1.82806 -1.821070  1.2157400 -3.566930 -7.91232  2.710110   0  10
5 -7.32776 -2.394770 -3.0798300 -1.436250  2.08822  1.429390   1  15
6  1.53227 -1.898840 -0.6726290  0.826352  4.01520  0.972757   1   7
~~~


## IV. Analysis workflow

GCIM analyses use PLink2 to analyze discovery data. 
1. Download plink2 software from the [Plink](https://www.cog-genomics.org/plink/2.0/) website and then specify the executable Plink software path.
   
~~~
plink_path <- "<plink_path>/plink2"
~~~ 
### 1. proposed causal direction 
#### 1.1. Quantitative outcome 
   
 ~~~
# Load required libraries
library(GxEprs)
library(GCIM)
# Set plink path
plink_path <- "<plink_path>/plink2"

# For quantitative traits, use corresponding functions
a <- GWAS_quantitative(plink_path, "DummyData", "Qcov_discovery_phen.txt", "Qcov_discovery_cov.txt")
b <- GWEIS_quantitative(plink_path, "DummyData", "Qphe_discovery.txt", "Qcov_discovery.txt")

# Extract and compute PRS
trd <- a[c("ID", "A1", "BETA")]
add <- b[c("ID", "A1", "ADD_BETA")]
gxe <- b[c("ID", "A1", "INTERACTION_BETA")]

q <- PRS_quantitative(plink_path, "DummyData", summary_input = add)
r <- PRS_quantitative(plink_path, "DummyData", summary_input = gxe) 
p <- PRS_quantitative(plink_path, "DummyData", summary_input = trd)

# Step 4: Run GCIM analysis with automatic saving and scaling
# This model specification corresponds to Model 4, as presented in the manuscript.
 result0 <- gcim_q0("Qphe_target.txt", "Qexp_target.txt", 
                  Add_PRS = q, Int_PRS = r, Cov_PRS = p)
#This model additionally incorporates the main effect of the PRS for the exposure variable. Details of this specification are provided in the model equation section of the Results. 
 result1 <- gcim_q1("Qphe_target.txt", "Qexp_target.txt", 
                  Add_PRS = q, Int_PRS = r, Cov_PRS = p)
 # Step 5: Access results and processed PRS objects
~~~

~~~
 print(result0$model_summary)
~~~

~~~
Coefficients:
                Estimate Std. Error t value Pr(>|t|)
(Intercept)     29.49287    3.67023   8.036 1.16e-13 ***
Add_PRS         -0.24579    0.36534  -0.673   0.5020
Int_PRS         -0.07551    0.40504  -0.186   0.8523
Covariate_Pheno -0.06776    0.36531  -0.185   0.8530
Conf_1          -0.31764    0.13158  -2.414   0.0168 *
Conf_2          -0.03142    0.04157  -0.756   0.4508
Conf_3           0.20054    0.23651   0.848   0.3976
Conf_4           0.38412    0.22684   1.693   0.0921 .
Conf_5          -0.38324    0.22064  -1.737   0.0841 .
Conf_6          -0.13732    0.15675  -0.876   0.3822
Conf_7           0.09351    0.06991   1.338   0.1827
Conf_8          -0.21469    0.22684  -0.946   0.3452
Conf_9          -0.29007    0.21864  -1.327   0.1863
Conf_10         -0.19900    0.18631  -1.068   0.2869
Conf_11         -0.01685    0.08666  -0.194   0.8460
Conf_12          0.06352    0.18782   0.338   0.7356
Conf_13          2.11258    0.67700   3.120   0.0021 **
Conf_14         -0.11187    0.06992  -1.600   0.1114
Int_PRS:Cov_PRS  0.13736    0.43681   0.314   0.7535
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
~~~

~~~
 print(result1$model_summary)
~~~

~~~
Coefficients:
                 Estimate Std. Error t value Pr(>|t|)
(Intercept)     28.183064   3.697016   7.623 1.37e-12 ***
Add_PRS         -0.309242   0.363647  -0.850  0.39624
Int_PRS         -0.128224   0.402497  -0.319  0.75042
Covariate_Pheno -0.008428   0.363450  -0.023  0.98152
Cov_PRS         -0.818939   0.406068  -2.017  0.04521 *
Conf_1          -0.348675   0.131386  -2.654  0.00867 **
Conf_2          -0.022893   0.041440  -0.552  0.58133
Conf_3           0.197041   0.234542   0.840  0.40196
Conf_4           0.451173   0.227383   1.984  0.04875 *
Conf_5          -0.386661   0.218801  -1.767  0.07889 .
Conf_6          -0.102698   0.156388  -0.657  0.51222
Conf_7           0.091554   0.069332   1.321  0.18834
Conf_8          -0.174943   0.225807  -0.775  0.43951
Conf_9          -0.245376   0.217937  -1.126  0.26170
Conf_10         -0.243769   0.186078  -1.310  0.19185
Conf_11         -0.015397   0.085936  -0.179  0.85800
Conf_12          0.039447   0.186630   0.211  0.83284
Conf_13          2.168696   0.671913   3.228  0.00148 **
Conf_14         -0.088428   0.070307  -1.258  0.21011
Int_PRS:Cov_PRS  0.048590   0.435383   0.112  0.91126
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
~~~

#### 1.2. Binary outcome 
The same data frames and analyses pipeline should be applied as used for the quantitative data presented above.
 
~~~
# Load required libraries
library(GxEprs)
library(GCIM)
# Set plink path
plink_path <- "<plink_path>/plink2"
# Step 1: Run GxEprs analysis for binary traits
#Note: If the exposure variable is binary, use it as is. If the exposure variable is quantitative, apply the `GWAS_quantitative` function to generate the object `a` specified as described below.
a <- GWAS_binary(plink_path, "DummyData", "Bcov_discovery_phen.txt", "Bcov_discovery_cov.txt")
b <- GWEIS_binary(plink_path, "DummyData", "Bphe_discovery.txt", "Bcov_discovery.txt")

# Step 2: Extract summary statistics
trd <- a[c("ID", "A1", "BETA")]
add <- b[c("ID", "A1", "ADD_BETA")]
gxe <- b[c("ID", "A1", "INTERACTION_BETA")]

# Step 3: Compute PRS for each component
#Note: If the exposure variable is binary, use it as is. If the exposure variable is quantitative, apply the `PRS_quantitative` function to generate the object `P` specified as described below.
p <- PRS_binary(plink_path, "DummyData", summary_input = trd)  # Covariate PRS
q <- PRS_binary(plink_path, "DummyData", summary_input = add)  # Additive PRS
r <- PRS_binary(plink_path, "DummyData", summary_input = gxe)  # Interaction PRS


# Step 4: Run GCIM analysis with automatic saving and scaling
 result0 <- gcim_b0("Bphe_target.txt", "Bexp_target.txt", 
                  Add_PRS = q, Int_PRS = r, Cov_PRS = p)
 result1 <- gcim_b1("Bphe_target.txt", "Bexp_target.txt", 
                  Add_PRS = q, Int_PRS = r, Cov_PRS = p) 
# # Step 5: Access results and processed PRS objects
 print(result0$model_summary)
~~~

~~~
Coefficients:
                Estimate Std. Error z value Pr(>|z|)
(Intercept)      4.81384    3.55634   1.354  0.17587
Add_PRS         -1.04471    1.05153  -0.994  0.32046
Int_PRS         -0.77657    1.02926  -0.754  0.45055
Covariate_Pheno  1.20406    0.66778   1.803  0.07137 .
Conf_1           0.04113    0.36059   0.114  0.90920
Conf_2           0.06771    0.09676   0.700  0.48409
Conf_3          -0.04767    0.03913  -1.218  0.22313
Conf_4           0.26105    0.23520   1.110  0.26704
Conf_5           0.15014    0.23389   0.642  0.52092
Conf_6           0.01001    0.19982   0.050  0.96003
Conf_7          -0.05855    0.16252  -0.360  0.71865
Conf_8           0.04285    0.08367   0.512  0.60854
Conf_9           0.08123    0.21203   0.383  0.70165
Conf_10          0.25170    0.17100   1.472  0.14105
Conf_11          0.09185    0.17155   0.535  0.59237
Conf_12         -0.01046    0.08043  -0.130  0.89654
Conf_13          0.04292    0.15733   0.273  0.78501
Conf_14         -0.21073    0.07682  -2.743  0.00609 **
Int_PRS:Cov_PRS -0.17201    0.48449  -0.355  0.72256

Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
~~~

~~~
 print(result1$model_summary)
~~~

~~~
Coefficients:
                 Estimate Std. Error z value Pr(>|z|)
(Intercept)      4.703893   3.618419   1.300  0.19361
Add_PRS         -1.024103   1.055763  -0.970  0.33204
Int_PRS         -0.751850   1.038001  -0.724  0.46887
Covariate_Pheno  1.208273   0.668928   1.806  0.07087 .
Cov_PRS          0.064796   0.417904   0.155  0.87678
Conf_1           0.035074   0.364660   0.096  0.92338
Conf_2           0.070014   0.098050   0.714  0.47519
Conf_3          -0.046839   0.039508  -1.186  0.23580
Conf_4           0.256237   0.237570   1.079  0.28078
Conf_5           0.148611   0.234302   0.634  0.52590
Conf_6           0.007641   0.201137   0.038  0.96970
Conf_7          -0.059300   0.162781  -0.364  0.71564
Conf_8           0.043153   0.083782   0.515  0.60651
Conf_9           0.078033   0.212945   0.366  0.71403
Conf_10          0.252354   0.170848   1.477  0.13966
Conf_11          0.096054   0.173808   0.553  0.58051
Conf_12         -0.010357   0.080341  -0.129  0.89743
Conf_13          0.042602   0.157742   0.270  0.78710
Conf_14         -0.209842   0.077039  -2.724  0.00645 **
Int_PRS:Cov_PRS -0.162503   0.489483  -0.332  0.73990
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
~~~

 ### 2. Reverse causal direction
 
To evaluate the **reverse causal direction**, re-analyze the same dataset by switching the roles of the exposure and outcome variables. This means treating the previously defined outcome variable as the new exposure, and the previous exposure variables as the new outcome. Rearrange the data using the same structure and formatting approach used for the proposed causal directions as mentioned above, ensuring consistency across analyses. The only difference should be the reassignment of variable roles.

#### 2.1. Quantitative outcome

~~~
# Load required libraries
library(GxEprs)
library(GCIM)
# Set plink path
plink_path <- "<plink_path>/plink2"

# For quantitative traits, use corresponding functions
a <- GWAS_quantitative(plink_path, "DummyData", "Qcov_disc_phen.txt", "Qcov_discovery_cov.txt")
b <- GWEIS_quantitative(plink_path, "DummyData", "Qcov_disc_out.txt", "Qcov_disc_exp.txt")

# Extract and compute PRS
trd <- a[c("ID", "A1", "BETA")]
add <- b[c("ID", "A1", "ADD_BETA")]
gxe <- b[c("ID", "A1", "INTERACTION_BETA")]

q <- PRS_quantitative(plink_path, "DummyData", summary_input = add)
r <- PRS_quantitative(plink_path, "DummyData", summary_input = gxe) 
p <- PRS_quantitative(plink_path, "DummyData", summary_input = trd)

# # Step 4: Run GCIM analysis with automatic saving and scaling
 result0 <- gcim_q0("Qexp_tar_out.txt", "Qexp_tar_exp.txt", 
                  Add_PRS = q, Int_PRS = r, Cov_PRS = p)
  result1 <- gcim_q1("Qexp_tar_out.txt", "Qexp_tar_exp.txt", 
                  Add_PRS = q, Int_PRS = r, Cov_PRS = p)
 # Step 5: Access results and processed PRS objects

~~~

~~~
 print(result0$model_summary)
~~~

~~~
Coefficients:
                 Estimate Std. Error t value Pr(>|t|)
(Intercept)      0.195159   0.875936   0.223   0.8239
Add_PRS          0.610907   0.658064   0.928   0.3545
Int_PRS          0.492336   0.653793   0.753   0.4524
Covariate_Pheno -0.002212   0.015268  -0.145   0.8850
Conf_1           0.031820   0.027421   1.160   0.2474
Conf_2          -0.006162   0.008459  -0.728   0.4673
Conf_3          -0.032466   0.048163  -0.674   0.5011
Conf_4          -0.049937   0.047016  -1.062   0.2896
Conf_5           0.084571   0.044748   1.890   0.0604 .
Conf_6           0.022523   0.032114   0.701   0.4840
Conf_7           0.001956   0.014289   0.137   0.8913
Conf_8           0.017399   0.046414   0.375   0.7082
Conf_9           0.047140   0.044522   1.059   0.2911
Conf_10         -0.024813   0.038478  -0.645   0.5198
Conf_11          0.012242   0.017624   0.695   0.4882
Conf_12          0.013844   0.038325   0.361   0.7184
Conf_13         -0.036297   0.141347  -0.257   0.7976
Conf_14          0.002498   0.014581   0.171   0.8641
Int_PRS:Cov_PRS  0.067951   0.082968   0.819   0.4139
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
~~~

~~~
 print(result1$model_summary)
~~~

~~~
Coefficients:
                  Estimate Std. Error t value Pr(>|t|)
(Intercept)      0.0588418  0.8763433   0.067   0.9465
Add_PRS          0.6821117  0.6567596   1.039   0.3004
Int_PRS          0.5977589  0.6543250   0.914   0.3622
Covariate_Pheno -0.0006245  0.0152351  -0.041   0.9673
Cov_PRS          0.1214049  0.0759022   1.599   0.1115
Conf_1           0.0320768  0.0273042   1.175   0.2416
Conf_2          -0.0044078  0.0084943  -0.519   0.6045
Conf_3          -0.0323717  0.0479567  -0.675   0.5005
Conf_4          -0.0463334  0.0468689  -0.989   0.3242
Conf_5           0.0856096  0.0445610   1.921   0.0563 .
Conf_6           0.0225979  0.0319764   0.707   0.4807
Conf_7           0.0020233  0.0142275   0.142   0.8871
Conf_8           0.0191299  0.0462282   0.414   0.6795
Conf_9           0.0386541  0.0446481   0.866   0.3878
Conf_10         -0.0236882  0.0383195  -0.618   0.5372
Conf_11          0.0119143  0.0175498   0.679   0.4981
Conf_12          0.0076775  0.0383555   0.200   0.8416
Conf_13         -0.0582612  0.1414111  -0.412   0.6808
Conf_14          0.0024458  0.0145183   0.168   0.8664
Int_PRS:Cov_PRS  0.0546324  0.0830318   0.658   0.5114
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
~~~

#### 2.2 Binary outcome 

~~~
# Load required libraries
library(GxEprs)
library(GCIM)
# Set plink path
plink_path <- "<plink_path>/plink2"
# Step 1: Run GxEprs analysis for binary traits
a <- GWAS_binary(plink_path, "DummyData", "Bcov_discovery_phen.txt", "Bcov_discovery_cov.txt")
b <- GWEIS_binary(plink_path, "DummyData", "Bcov_disc_out1.txt", "Bphe_disc_exp.txt")

# Step 2: Extract summary statistics
trd <- a[c("ID", "A1", "BETA")]
add <- b[c("ID", "A1", "ADD_BETA")]
gxe <- b[c("ID", "A1", "INTERACTION_BETA")]

# Step 3: Compute PRS for each component
q <- PRS_binary(plink_path, "DummyData", summary_input = add)  # Additive PRS
r <- PRS_binary(plink_path, "DummyData", summary_input = gxe)  # Interaction PRS
p <- PRS_binary(plink_path, "DummyData", summary_input = trd)  # Covariate PRS

# Step 4: Run GCIM analysis with automatic saving and scaling
 result0 <- gcim_b0("Bexp_tar_out.txt", "Bexp_tar_exp.txt", 
                  Add_PRS = q, Int_PRS = r, Cov_PRS = p)
 result1 <- gcim_b1("Bexp_tar_out.txt", "Bexp_tar_exp.txt", 
                  Add_PRS = q, Int_PRS = r, Cov_PRS = p) 
#Step 5: Access results and processed PRS objects
~~~

~~~
print(result0$model_summary)
~~~

~~~
Coefficients:
                Estimate Std. Error z value Pr(>|z|)
(Intercept)     -2.14559    1.77175  -1.211   0.2259
Add_PRS         -0.33958    2.28221  -0.149   0.8817
Int_PRS         -0.37669    2.28989  -0.165   0.8693
Covariate_Pheno  1.05839    0.59969   1.765   0.0776 .
Conf_1           0.44229    0.17231   2.567   0.0103 *
Conf_2           0.02413    0.05763   0.419   0.6755
Conf_3          -0.01017    0.01908  -0.533   0.5941
Conf_4          -0.14360    0.10762  -1.334   0.1821
Conf_5          -0.04562    0.11996  -0.380   0.7038
Conf_6           0.21204    0.10520   2.016   0.0438 *
Conf_7           0.13105    0.07645   1.714   0.0865 .
Conf_8          -0.07325    0.03648  -2.008   0.0446 *
Conf_9           0.10046    0.09985   1.006   0.3143
Conf_10         -0.06210    0.08605  -0.722   0.4705
Conf_11          0.07371    0.08837   0.834   0.4042
Conf_12          0.02785    0.03838   0.726   0.4681
Conf_13          0.04458    0.08003   0.557   0.5775
Conf_14          0.06213    0.03232   1.923   0.0545 .
Int_PRS:Cov_PRS -0.15908    0.23086  -0.689   0.4908
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
~~~

~~~
print(result0$model_summary)
~~~

~~~
Coefficients:
                Estimate Std. Error z value Pr(>|z|)
(Intercept)     -2.17581    1.77526  -1.226  0.22034
Add_PRS          0.10961    2.32367   0.047  0.96238
Int_PRS          0.09498    2.33450   0.041  0.96755
Covariate_Pheno  1.09816    0.60304   1.821  0.06860 .
Cov_PRS          0.22051    0.20304   1.086  0.27746
Conf_1           0.47750    0.17652   2.705  0.00683 **
Conf_2           0.02046    0.05789   0.353  0.72379
Conf_3          -0.01125    0.01919  -0.586  0.55778
Conf_4          -0.14761    0.10779  -1.369  0.17085
Conf_5          -0.06223    0.12148  -0.512  0.60848
Conf_6           0.20386    0.10582   1.926  0.05405 .
Conf_7           0.12744    0.07665   1.663  0.09638 .
Conf_8          -0.07427    0.03659  -2.030  0.04238 *
Conf_9           0.08466    0.10168   0.833  0.40506
Conf_10         -0.06056    0.08616  -0.703  0.48213
Conf_11          0.08213    0.08885   0.924  0.35531
Conf_12          0.02502    0.03868   0.647  0.51770
Conf_13          0.04305    0.08066   0.534  0.59352
Conf_14          0.06662    0.03273   2.035  0.04181 *
Int_PRS:Cov_PRS -0.15114    0.23407  -0.646  0.51848
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
~~~
 

