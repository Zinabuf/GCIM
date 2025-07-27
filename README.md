---
# <h1 align="center">GCIM</h1>
---

The genetic causality inference model(GCIM) is a statistical method for detecting the direction of causation in GxE interaction studies. 
- 
 Authors: Zinabu Fentaw, Dovini Jayasinghe, S.Hong Lee
-

_GCIM is a novel statistical method that extends beyond traditional PRS×E approaches.
    It systematically evaluates both the proposed and reverse causal directions, rather than relying
    on the assumptions of prior causal directions. GCIM uses the polygenic risk score (PRS) of the
    exposure trait rather than its observed phenotype to eliminate spurious covariance between the
    exposure and the G×E component of the outcome. This helps prevent bias when a true causal effect
    exists in the proposed direction. By explicitly testing both directions, GCIM offers researchers
    data-driven insight into the likely causal directions._

_NB: The proposed direction of causation refers to the hypothesized G×E interaction in which the exposure affects the outcome, aligning with the researcher’s primary interest. In contrast, the reverse direction tests the opposite relationship, switching the roles of exposure and outcome to evaluate whether the assumed causal direction is valid._
   
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
The dataset is split into a **discovery set (80%)** and a **target set (20%)** to ensure consistency across genetic, outcome, exposure, and confounder variables. Genetic data must be provided in **PLINK binary format** (.bed, .bim, .fam). The outcome file should contain **FID, IID, and phenotype values**, with case-control coding as follows: in the discovery dataset, use PLINK's default (**1 = Control, 2 = Case**), and in the target dataset, use binary coding (**0 = Control, 1 = Case**). Exposure and confounder variables should be organized in at least **m columns** with the format: FID, IID, exposure, exposure_square, confounder1 to confounder_m, and partitioned consistently across datasets. For binary outcomes, the exposure_square column should contain constant values (e.g., 2, 2, 2, ..., 2) to ensure zero variance. This standardized structure enables proper alignment across all data types, supporting valid and reliable estimation of G×E interactions.


A Guide for GCIM analyses

GCIM analyses use PLink2 to analyze discovery data, and the package is compatible with the Linux operating system. 
1. Download plink2 software from the [Plink](https://www.cog-genomics.org/plink/2.0/) website and then specify the executable Plink file path.
   
~~~
plink_path <- "<plink_path>/plink2"
~~~
**Discovery dataset**: In the discovery dataset, please use the data input from the GxEprs data input format, then construct the PRS of the exposure using [Plink](https://www.cog-genomics.org/plink/2.0/) based on the GWAS summary statistics for the PRS of the target samples.
**Target dataset**: The target dataset for the model is also similar in data format, except, the use of PRS of the exposure variable(PRS of E) rather than the use of entire exposure values, and also includes the constant values for the fourth column, which is the square of the third column in the GxEprs model. 

**Example data**
To conduct a GCIM analysis, we must use the same data format as GxEprs, especially in the discovery dataset. However, there is a slight difference in the target dataset, as the square of the exposure variables is not required. Here is an example of analysis using the R script in the example directories. 
1. Quantitative outcome
1.1. quantitative exposure
   
 ~~~
# Load required libraries
library(GxEprs)
library(GCIM)
# Set plink path
plink_path <- "<plink_path>/plink2"

# For quantitative traits, use corresponding functions
a <- GWAS_quantitative(plink_path, "DummyData", "Qphe_discovery.txt", "Qcov_discovery.txt")
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

 Result displayed using this approach

~~~
Coefficients:
                Estimate Std. Error t value Pr(>|t|)
(Intercept)     29.59421    3.67899   8.044  1.1e-13 ***
Add_PRS         -0.23837    0.37043  -0.643  0.52072
Int_PRS         -0.14403    0.40516  -0.356  0.72263
Covariate_Pheno -0.07524    0.36446  -0.206  0.83668
Conf_1          -0.31858    0.13184  -2.416  0.01667 *
Conf_2          -0.03200    0.04150  -0.771  0.44166
Conf_3           0.21174    0.23556   0.899  0.36992
Conf_4           0.38908    0.22767   1.709  0.08917 .
Conf_5          -0.37469    0.22106  -1.695  0.09180 .
Conf_6          -0.13548    0.15685  -0.864  0.38887
Conf_7           0.09533    0.06990   1.364  0.17434
Conf_8          -0.22050    0.22720  -0.971  0.33309
Conf_9          -0.28279    0.21921  -1.290  0.19867
Conf_10         -0.19454    0.18694  -1.041  0.29942
Conf_11         -0.02127    0.08668  -0.245  0.80645
Conf_12          0.06014    0.18822   0.320  0.74971
Conf_13          2.10844    0.67541   3.122  0.00209 **
Conf_14         -0.10898    0.07029  -1.550  0.12280
Int_PRS:Cov_PRS -0.08119    0.39788  -0.204  0.83853
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 4.624 on 181 degrees of freedom
Multiple R-squared:  0.1385,    Adjusted R-squared:  0.05285
F-statistic: 1.617 on 18 and 181 DF,  p-value: 0.05972
~~~

1.2. Binary exposure variables

~~~
# Load required libraries
library(GxEprs)
library(GCIM)
# Set plink path
plink_path <- "<plink_path>/plink2"

# For quantitative traits, use corresponding functions
a <- GWAS_quantitative(plink_path, "DummyData", "Qphe_discovery.txt", "Qcov_discovery.txt")
b <- GWEIS_quantitative(plink_path, "DummyData", "Qphe_discovery.txt", "Qcov_discovery.txt")

# Extract and compute PRS
trd <- a[c("ID", "A1", "BETA")]
add <- b[c("ID", "A1", "ADD_BETA")]
gxe <- b[c("ID", "A1", "INTERACTION_BETA")]

q <- PRS_quantitative(plink_path, "DummyData", summary_input = add)
r <- PRS_quantitative(plink_path, "DummyData", summary_input = gxe) 
p <- PRS_quantitative(plink_path, "DummyData", summary_input = trd)

# # Step 4: Run GCIM analysis with automatic saving and scaling
 result <- gcim_q("Qphe_target.txt", "Bexp_target.txt", 
                  Add_PRS = q, Int_PRS = r, Cov_PRS = p) 
 # Step 5: Access results and processed PRS objects
 print(result$model_summary)
~~~

Results

~~~
Coefficients:
                 Estimate Std. Error t value Pr(>|t|)
(Intercept)     26.241837   3.859572   6.799 1.48e-10 ***
Add_PRS         -0.274757   0.396816  -0.692   0.4896
Int_PRS         -0.159843   0.408864  -0.391   0.6963
Covariate_Pheno  0.292060   0.717268   0.407   0.6844
Conf_1          -0.104441   0.364834  -0.286   0.7750
Conf_2          -0.067684   0.129205  -0.524   0.6010
Conf_3          -0.007134   0.041343  -0.173   0.8632
Conf_4          -0.009463   0.236466  -0.040   0.9681
Conf_5          -0.002942   0.258273  -0.011   0.9909
Conf_6           0.109983   0.229403   0.479   0.6322
Conf_7          -0.019449   0.166482  -0.117   0.9071
Conf_8          -0.123671   0.080512  -1.536   0.1263
Conf_9           0.544493   0.218697   2.490   0.0137 *
Conf_10         -0.059208   0.189093  -0.313   0.7546
Conf_11         -0.224780   0.195958  -1.147   0.2529
Conf_12         -0.075448   0.083313  -0.906   0.3664
Conf_13         -0.154566   0.177106  -0.873   0.3840
Conf_14          0.082395   0.070177   1.174   0.2419
Int_PRS:Cov_PRS  0.040304   0.401366   0.100   0.9201
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
~~~

 2. Binary outcome
    2.1. Quantitative exposure

~~~
# Load required libraries
library(GxEprs)
library(GCIM)
# Set plink path
plink_path <- "<plink_path>/plink2"
# Step 1: Run GxEprs analysis for binary traits
a <- GWAS_binary(plink_path, "DummyData", "Bphe_discovery.txt", "Bcov_discovery.txt")
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
result <- gcim_b("Bphe_target.txt", "Qexp_target.txt", 
                  Add_PRS = q, Int_PRS = r, Cov_PRS = p) 
# # Step 5: Access results and processed PRS objects
 print(result$model_summary)
~~~

Results

~~~
    Coefficients:
                 Estimate Std. Error z value Pr(>|z|)
(Intercept)     -9.775706   3.887043  -2.515   0.0119 *
Add_PRS         -0.712033   0.718365  -0.991   0.3216
Int_PRS         -0.451091   0.703152  -0.642   0.5212
Covariate_Pheno -0.176246   0.352987  -0.499   0.6176
Conf_1          -0.111470   0.128940  -0.865   0.3873
Conf_2           0.065496   0.041232   1.588   0.1122
Conf_3          -0.053923   0.230885  -0.234   0.8153
Conf_4          -0.003921   0.203141  -0.019   0.9846
Conf_5           0.003561   0.195466   0.018   0.9855
Conf_6           0.345131   0.160860   2.146   0.0319 *
Conf_7          -0.067287   0.064602  -1.042   0.2976
Conf_8           0.218933   0.203603   1.075   0.2822
Conf_9           0.025972   0.185094   0.140   0.8884
Conf_10         -0.185153   0.175694  -1.054   0.2920
Conf_11         -0.037268   0.067422  -0.553   0.5804
Conf_12         -0.105283   0.165413  -0.636   0.5245
Conf_13          0.532193   0.618281   0.861   0.3894
Conf_14          0.096496   0.065456   1.474   0.1404
Int_PRS:Cov_PRS  0.767477   0.440455   1.742   0.0814 .

Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

(Dispersion parameter for binomial family taken to be 1)

    Null deviance: 106.554  on 199  degrees of freedom
Residual deviance:  88.758  on 181  degrees of freedom
AIC: 126.76
~~~

2.2. Binary exposure variables

~~~
# Load required libraries
library(GxEprs)
library(GCIM)
# Set plink path
plink_path <- "<plink_path>/plink2"
# Step 1: Run GxEprs analysis for binary traits
a <- GWAS_binary(plink_path, "DummyData", "Bphe_discovery.txt", "Bcov_discovery.txt")
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
~~~~

Results 

~~~
Coefficients:
                 Estimate Std. Error z value Pr(>|z|)
(Intercept)      6.588878   3.852194   1.710  0.08719 .
Add_PRS         -1.962138   0.945819  -2.075  0.03803 *
Int_PRS         -1.425874   0.827295  -1.724  0.08479 .
Covariate_Pheno  1.219534   0.674446   1.808  0.07058 .
Conf_1           0.114997   0.362602   0.317  0.75113
Conf_2           0.110126   0.100644   1.094  0.27386
Conf_3          -0.057149   0.041885  -1.364  0.17243
Conf_4           0.356284   0.244408   1.458  0.14491
Conf_5           0.161170   0.250772   0.643  0.52042
Conf_6           0.063442   0.193255   0.328  0.74270
Conf_7          -0.105235   0.170714  -0.616  0.53760
Conf_8           0.051238   0.089461   0.573  0.56682
Conf_9           0.088239   0.220160   0.401  0.68857
Conf_10          0.298192   0.173146   1.722  0.08503 .
Conf_11          0.132276   0.184009   0.719  0.47223
Conf_12         -0.007419   0.083339  -0.089  0.92907
Conf_13          0.071013   0.158253   0.449  0.65363
Conf_14         -0.220482   0.080538  -2.738  0.00619 **
Int_PRS:Cov_PRS  0.706437   0.483414   1.461  0.14392
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

(Dispersion parameter for binomial family taken to be 1)

    Null deviance: 106.554  on 199  degrees of freedom
Residual deviance:  81.432  on 181  degrees of freedom
AIC: 119.43
~~~    
Test the reverse directions by changing the role of exposure and outcome using the same data as above.  
