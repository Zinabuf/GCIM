---
# <h1 align="center">GCIM</h1>
---

The genetic causality inference model (GCIM) is a statistical method for detecting the direction of causation in Genotype-by-environment(GxE) interaction studies. 

- 
 Authors: Zinabu Fentaw, Dovini Jayasinghe, S.Hong Lee
-

_<div align="justify">GCIM is a novel statistical method that extends beyond traditional PRS×E approaches.
    It systematically evaluates both the proposed and reverse causal directions, rather than relying
    on the assumptions of prior causal directions. By explicitly testing both directions, GCIM offers researchers
    data-driven insight into the likely causal directions.</div>_

_<div align="justify">NB: The proposed direction of causation refers to the hypothesized (G×E) interaction in which the exposure affects the outcome, aligning with the researcher’s primary interest. In contrast, the reverse direction test evaluates the opposite relationship, switching the roles of exposure and outcome to assess whether the assumed causal direction is valid. Once the data has been properly prepared, begin by testing the causal direction based on the proposed direction of interest, that is the direction specified by the researcher. After this primary analysis, assess the reverse causal direction by switching the roles of the exposure and outcome variables in both the discovery and target datasets. This involves treating the original outcome as the exposure and the original exposure as the outcome, ensuring consistency in data structure and formatting across both analyses.</div>_
   
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
~~~

Load the library

~~~
library(GxEprs)
library(GCIM)
~~~

**Data Preparation:**

To ensure consistent and reliable estimation of G×E, the dataset should be split into two subsets: **a discovery dataset (80%)** and **a target dataset (20%)**. This split should maintain consistency across genetic data, outcome variables, exposure variables, and potential confounders. The genetic data must be in **PLINK binary format**, comprising three files: `.bed`, `.bim`, and `.fam`. The **outcome file** should include `FID`, `IID`, and the outcome variable. For binary outcomes, follow standard coding conventions: use **PLINK’s default coding (1 = Control, 2 = Case)** for the **discovery dataset**, and use **binary coding (0 = Control, 1 = Case)** for the **target dataset**. The **exposure and confounder file for the discovery dataset** should contain at least `m` columns (Minimum of 3) with the following format: `FID`, `IID`, `exposure`, `exposure_squared`, `confounder_1`, ..., `confounder_m`. Note that **exposure-squared is not required** in the **target dataset**. Additionally, for conducting a genome-wide association study (GWAS) of the exposure variable in the discovery dataset (for use in computing the Polygenic Risk Score of the exposure), the exposure data should also be formatted separately as:

* `FID`, `IID`, and **exposure values**
* A covariate file in the format: `FID`, `IID`, `confounder_1`, `constant_value`, `confounder_2`, ..., `confounder_m`

This standardized format ensures that all variables such as genetic, exposure, outcome, and confounders are properly aligned across discovery and target datasets, facilitating valid and replicable G×E interaction analysis.
All GWAS, GWEIS, and polygenic risk score (PRS) construction steps are performed using the [GxEprs](https://github.com/DoviniJ/GxEprs) R package 
, while the regression analyses for both binary and quantitative outcomes are conducted using the GCIM R package. 

A Guide for GCIM analyses

GCIM analyses use PLink2 to analyze discovery data, and the package is compatible with the Linux operating system. 
1. Download plink2 software from the [Plink](https://www.cog-genomics.org/plink/2.0/) website and then specify the executable Plink software path.
   
~~~
plink_path <- "<plink_path>/plink2"
~~~ 

**Example data**
<div align="justify">To conduct a GCIM analysis, we must use the same data format as GxEprs, especially in the discovery dataset. However, there is a slight difference in the target dataset, as the square of the exposure variables is not required. Here is an example of analysis using the R script in the example directories.</div>  

Data structure
To 
GWAS for quantitative exposure variables with the phenotype 

~~~
    V1   V2          V3
1 ID_1 ID_1 -0.64402046
2 ID_2 ID_2 -0.02786981
3 ID_3 ID_3  2.12865748
4 ID_4 ID_4  2.12865748
5 ID_5 ID_5 -0.95209579
6 ID_6 ID_6 -0.02786981
~~~
The covariate for adjustment 

~~~
    V1   V2           V3        V4 V5        V6       V7         V8        V9
1 ID_1 ID_1 0.4147623548 -3.831420 64 -14.03640 5.517420  0.0714337  5.662630
2 ID_2 ID_2 0.0007767262  0.614044 66 -10.85050 2.119980 -0.8828830 -0.441662
3 ID_3 ID_3 4.5311826719 -0.237792 55  -9.75369 3.183430 -2.0979300  6.873450
4 ID_4 ID_4 4.5311826719  6.698660 47  -9.07045 0.956878 -2.4840700  1.063590
5 ID_5 ID_5 0.9064863904 -1.614230 59 -12.93790 1.294610 -1.7997300  1.444040
6 ID_6 ID_6 0.0007767262 -4.389270 52 -11.85160 0.888978 -2.7231000  1.116810
        V10       V11         V12       V13        V14       V15 V16 V17
1  0.865562 -2.269570 -0.09658590 -2.354970  1.0588900  0.195302   0   7
2 -2.641770  2.789440  0.52458600  2.671340 -2.6372400 -0.998764   1  20
3 11.377700  2.969610 -1.11879000  0.873649  3.3552300 -4.578310   1  10
4 -3.132470  2.123200 -0.00976751  0.820582  0.0305345  1.630300   1  20
5 -6.828980 -2.967950 -2.91577000 -1.828810  7.1589200  2.109160   1  20
6 -3.646760 -0.594538 -1.75430000 -0.716014 -2.3906700  1.312950   1  10
~~~
 For GWEIS analyses, the outcome variables 

 ~~~
V1   V2      V3
1 ID_1 ID_1 31.6534
2 ID_2 ID_2 25.5035
3 ID_3 ID_3 26.7391
4 ID_4 ID_4 25.5271
5 ID_5 ID_5 26.7165
6 ID_6 ID_6 38.8272
~~~
GWEIS Exposure variables are 

~~~
    V1   V2          V3           V4        V5 V6        V7       V8         V9
1 ID_1 ID_1 -0.64402046 0.4147623548 -3.831420 64 -14.03640 5.517420  0.0714337
2 ID_2 ID_2 -0.02786981 0.0007767262  0.614044 66 -10.85050 2.119980 -0.8828830
3 ID_3 ID_3  2.12865748 4.5311826719 -0.237792 55  -9.75369 3.183430 -2.0979300
4 ID_4 ID_4  2.12865748 4.5311826719  6.698660 47  -9.07045 0.956878 -2.4840700
5 ID_5 ID_5 -0.95209579 0.9064863904 -1.614230 59 -12.93790 1.294610 -1.7997300
6 ID_6 ID_6 -0.02786981 0.0007767262 -4.389270 52 -11.85160 0.888978 -2.7231000
        V10       V11       V12         V13       V14        V15       V16 V17
1  5.662630  0.865562 -2.269570 -0.09658590 -2.354970  1.0588900  0.195302   0
2 -0.441662 -2.641770  2.789440  0.52458600  2.671340 -2.6372400 -0.998764   1
3  6.873450 11.377700  2.969610 -1.11879000  0.873649  3.3552300 -4.578310   1
4  1.063590 -3.132470  2.123200 -0.00976751  0.820582  0.0305345  1.630300   1
5  1.444040 -6.828980 -2.967950 -2.91577000 -1.828810  7.1589200  2.109160   1
6  1.116810 -3.646760 -0.594538 -1.75430000 -0.716014 -2.3906700  1.312950   1
  V18
1   7
2  20
3  10
4  20
5  20
6  10
~~~

Target data set 
Quantitative outcome should look like 

~~~
   V1     V2      V3
1 ID_801 ID_801 26.5723
2 ID_802 ID_802 20.2632
3 ID_803 ID_803 27.7365
4 ID_804 ID_804 18.7500
5 ID_805 ID_805 23.3025
6 ID_806 ID_806 24.9871
~~~
 The exposure variable for the target dataset should be 

 ~~~
      V1     V2          V3        V4 V5       V6      V7        V8        V9
1 ID_801 ID_801 -0.64402046 -3.826590 69 -13.8514 3.96080 -1.788050 0.0692473
2 ID_802 ID_802 -0.95209579  2.065150 60 -12.2438 4.04169 -0.905739 5.9656000
3 ID_803 ID_803 -0.02786981 -0.795863 62 -10.9195 6.91985 -2.920880 1.2601900
4 ID_804 ID_804 -0.64402046 -2.620880 67  -9.9271 4.10960 -2.354540 0.7190210
5 ID_805 ID_805  0.28020552 -3.331640 67 -11.8637 5.88272  1.072880 2.7448800
6 ID_806 ID_806  0.89635617  3.252030 61 -11.5364 5.79318 -2.311240 4.5023300
       V10       V11        V12       V13      V14       V15 V16 V17
1 -6.32556  2.853590  1.0851600 -1.303040  3.41659  1.415770   0   7
2  8.35545 -1.435760 -0.6181530  0.746918  5.11019 -0.207188   1  19
3 -5.56624 -0.552624 -0.0756095 -0.910047 -1.33896  1.726360   0   7
4 -1.82806 -1.821070  1.2157400 -3.566930 -7.91232  2.710110   0  10
5 -7.32776 -2.394770 -3.0798300 -1.436250  2.08822  1.429390   1  15
6  1.53227 -1.898840 -0.6726290  0.826352  4.01520  0.972757   1   7
~~~

1. Quantitative outcome
   
 ~~~
# Load required libraries
library(GxEprs)
library(GCIM)
# Set plink path
plink_path <- "<plink_path>/plink2"

# For quantitative traits, use corresponding functions
a <- GWAS_quantitative(plink_path, "DummyData", "Qcov_discovery_phen.txt", "Qcov_discovery_exp.txt")
b <- GWEIS_quantitative(plink_path, "DummyData", "Qphe_discovery.txt", "Qcov_discovery.txt")

# Extract and compute PRS
trd <- a[c("ID", "A1", "BETA")]
add <- b[c("ID", "A1", "ADD_BETA")]
gxe <- b[c("ID", "A1", "INTERACTION_BETA")]

q <- PRS_quantitative(plink_path, "DummyData", summary_input = add)
r <- PRS_quantitative(plink_path, "DummyData", summary_input = gxe) 
p <- PRS_quantitative(plink_path, "DummyData", summary_input = trd)

# # Step 4: Run GCIM analysis with automatic saving and scaling
 result <- gcim_q("Qphe_target.txt", "Qexp_target.txt", 
                  Add_PRS = q, Int_PRS = r, Cov_PRS = p)
 
 # Step 5: Access results and processed PRS objects
 print(result$model_summary)
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
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 4.623 on 181 degrees of freedom
Multiple R-squared:  0.1386,    Adjusted R-squared:  0.05298
F-statistic: 1.618 on 18 and 181 DF,  p-value: 0.05934
~~~

 2. Binary outcome 

~~~
# Load required libraries
library(GxEprs)
library(GCIM)
# Set plink path
plink_path <- "<plink_path>/plink2"
# Step 1: Run GxEprs analysis for binary traits
a <- GWAS_binary(plink_path, "DummyData", "Bcov_discovery_phen.txt", "Bcov_discovery_exp.txt")
b <- GWEIS_binary(plink_path, "DummyData", "Bphe_discovery.txt", "Bcov_discovery.txt")

# Step 2: Extract summary statistics
trd <- a[c("ID", "A1", "BETA")]
add <- b[c("ID", "A1", "ADD_BETA")]
gxe <- b[c("ID", "A1", "INTERACTION_BETA")]

# Step 3: Compute PRS for each component
q <- PRS_binary(plink_path, "DummyData", summary_input = add)  # Additive PRS
r <- PRS_binary(plink_path, "DummyData", summary_input = gxe)  # Interaction PRS
p <- PRS_binary(plink_path, "DummyData", summary_input = trd)  # Covariate PRS

# Step 4: Run GCIM analysis with automatic saving and scaling
 result <- gcim_b("Bphe_target.txt", "Bexp_target.txt", 
                  Add_PRS = q, Int_PRS = r, Cov_PRS = p) 
# # Step 5: Access results and processed PRS objects
 print(result$model_summary)
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

(Dispersion parameter for binomial family taken to be 1)

    Null deviance: 106.554  on 199  degrees of freedom
Residual deviance:  86.566  on 181  degrees of freedom
AIC: 124.57
~~~
    
<div align="justify">To evaluate the **reverse direction of causation**, re-analyze the same dataset by switching the roles of the exposure and outcome variables. This means treating the previously defined outcome variable as the new exposure, and the previous exposure variables as the new outcome. Rearrange the data using the same structure and formatting approach used for the proposed causal directions, ensuring consistency across analyses. The only difference should be the reassignment of variable roles.</div>  

For GWAS analyses, the data should be prepared like 

~~~
      V1     V2 V3
1   ID_1   ID_1  1
2  ID_10  ID_10  1
3 ID_100 ID_100  1
4 ID_101 ID_101  2
5 ID_102 ID_102  1
6 ID_103 ID_103  1
~~~

The exposure variables for GWAS should be 

~~~
      V1     V2         V3         V4       V5 V6       V7      V8        V9
1   ID_1   ID_1 0.62000476  0.7874038 -3.04026 45 -12.0480 2.17634 -0.940322
2  ID_10  ID_10 0.02716622 -0.1648218 -4.26064 67 -13.9720 4.60617 -4.132670
3 ID_100 ID_100 0.06735963  0.2595373 -3.75818 60 -13.7169 4.72924  0.696672
4 ID_101 ID_101 0.04674951 -0.2162164  1.51113 49 -11.6284 5.84781 -1.100270
5 ID_102 ID_102 2.02783487 -1.4240207 -3.23187 68 -12.1182 2.86926  0.405321
6 ID_103 ID_103 0.05924497 -0.2434029  3.81645 58 -12.0755 2.45074 -0.181400
         V10      V11       V12      V13       V14        V15       V16 V17
1 -0.4463510 -5.45685 -2.531610 -2.13435 -1.956230  -3.827920 -0.380636  10
2  1.7448900  7.77612 -0.285812  1.11394 -3.279560 -15.311500  3.590800   7
3 -1.2816200 -7.51913 -1.729600  2.93850 -2.630840   3.709820  2.452090  20
4  2.1046600 -6.70373 -2.120090 -1.55533 -5.620780   6.883290  5.733580  20
5  0.0536371 -2.81541 -0.490841  2.77044 -2.807250  -0.815034  4.658500   7
6  2.0906300  1.12865 -0.946887  1.83857 -0.916254   3.928980  4.353590  10
~~~
Outcome for GWEIS
~~~
 V1   V2 V3
1 ID_1 ID_1  1
2 ID_2 ID_2  2
3 ID_3 ID_3  2
4 ID_4 ID_4  1
5 ID_5 ID_5  2
6 ID_6 ID_6  2"
~~~

Exposure should also look like 

~~~
"      V1     V2 V3         V4         V5       V6 V7       V8      V9       V10
1   ID_1   ID_1  1 0.62000476  0.7874038 -3.04026 45 -12.0480 2.17634 -0.940322
2  ID_10  ID_10  1 0.02716622 -0.1648218 -4.26064 67 -13.9720 4.60617 -4.132670
3 ID_100 ID_100  1 0.06735963  0.2595373 -3.75818 60 -13.7169 4.72924  0.696672
4 ID_101 ID_101  2 0.04674951 -0.2162164  1.51113 49 -11.6284 5.84781 -1.100270
5 ID_102 ID_102  1 2.02783487 -1.4240207 -3.23187 68 -12.1182 2.86926  0.405321
6 ID_103 ID_103  1 0.05924497 -0.2434029  3.81645 58 -12.0755 2.45074 -0.181400
         V11      V12       V13      V14       V15        V16       V17 V18
1 -0.4463510 -5.45685 -2.531610 -2.13435 -1.956230  -3.827920 -0.380636  10
2  1.7448900  7.77612 -0.285812  1.11394 -3.279560 -15.311500  3.590800   7
3 -1.2816200 -7.51913 -1.729600  2.93850 -2.630840   3.709820  2.452090  20
4  2.1046600 -6.70373 -2.120090 -1.55533 -5.620780   6.883290  5.733580  20
5  0.0536371 -2.81541 -0.490841  2.77044 -2.807250  -0.815034  4.658500   7
6  2.0906300  1.12865 -0.946887  1.83857 -0.916254   3.928980  4.353590  10"
~~~

arget dataset

outcome in the target dataset 

~~~
     V1     V2 V3
1 ID_801 ID_801  1
2 ID_802 ID_802  1
3 ID_803 ID_803  0
4 ID_804 ID_804  0
5 ID_805 ID_805  0
6 ID_806 ID_806  1
~~~

Exposure variables in the targetdataset 

~~~
      V1      V2 V3          V4       V5 V6        V7      V8        V9
1 ID_1000 ID_1000  0  0.42094546 -2.91833 67 -11.79100 6.19780 -1.364380
2  ID_801  ID_801  0 -0.42082298 -4.12263 57 -13.51850 5.40198 -4.819940
3  ID_802  ID_802  1 -0.08055833  2.92534 56 -13.62360 3.21643 -0.856048
4  ID_803  ID_803  0 -1.32752645 -3.09118 61  -9.94475 3.60562 -0.917639
5  ID_804  ID_804  0  0.69800724  4.58829 49 -12.54710 4.09467 -2.589510
6  ID_805  ID_805  1 -0.65798161 -3.53948 56 -12.79500 2.91524 -2.727940
        V10      V11       V12       V13       V14       V15      V16 V17
1 -0.445164 -9.70429 -0.671265  2.258490 -2.006960  1.782850 -1.97524  10
2  0.664494 -4.92217 -0.451329  3.146770  0.427040  0.821306 -2.77705   7
3  0.750187 -2.01798 -0.350832  5.101410  2.180700 -6.043430  1.78928  19
4  0.905664 -5.09843 -1.163290 -1.881020 -1.241540  0.699574  2.24420  20
5  6.068980 12.98220 -0.704179  2.903570 -0.334968  5.042740  0.66175  10
6  3.615550  3.92957 -2.938990 -0.454737  2.310130  2.517830 -4.15592   7
~~~

3. Quantitative outcome

~~~
# Load required libraries
library(GxEprs)
library(GCIM)
# Set plink path
plink_path <- "<plink_path>/plink2"

# For quantitative traits, use corresponding functions
a <- GWAS_quantitative(plink_path, "DummyData", "Qcov_disc_phen.txt", "Qcov_discovery_exp.txt")
b <- GWEIS_quantitative(plink_path, "DummyData", "Qcov_disc_out.txt", "Qcov_disc_exp.txt")

# Extract and compute PRS
trd <- a[c("ID", "A1", "BETA")]
add <- b[c("ID", "A1", "ADD_BETA")]
gxe <- b[c("ID", "A1", "INTERACTION_BETA")]

q <- PRS_quantitative(plink_path, "DummyData", summary_input = add)
r <- PRS_quantitative(plink_path, "DummyData", summary_input = gxe) 
p <- PRS_quantitative(plink_path, "DummyData", summary_input = trd)

# # Step 4: Run GCIM analysis with automatic saving and scaling
 result <- gcim_q("Qexp_tar_out.txt", "Qexp_tar_exp.txt", 
                  Add_PRS = q, Int_PRS = r, Cov_PRS = p)
 
 # Step 5: Access results and processed PRS objects
 print(result$model_summary)
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

Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.9414 on 181 degrees of freedom
Multiple R-squared:  0.06849,   Adjusted R-squared:  -0.02415
F-statistic: 0.7393 on 18 and 181 DF,  p-value: 0.7673
~~~

4. Binary outcome 

~~~
# Load required libraries
library(GxEprs)
library(GCIM)
# Set plink path
plink_path <- "<plink_path>/plink2"
# Step 1: Run GxEprs analysis for binary traits
a <- GWAS_binary(plink_path, "DummyData", "Bcov_discovery_phen.txt", "Bcov_discovery_exp.txt")
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
 result <- gcim_b("Bexp_tar_out.txt", "Bexp_tar_exp.txt", 
                  Add_PRS = q, Int_PRS = r, Cov_PRS = p) 
# # Step 5: Access results and processed PRS objects
 print(result$model_summary)
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
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

(Dispersion parameter for binomial family taken to be 1)

    Null deviance: 272.74  on 199  degrees of freedom
Residual deviance: 247.61  on 181  degrees of freedom
AIC: 285.61
~~~
