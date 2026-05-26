---
# <h1 align="center">GCIM</h1>
---

The **genetic causality inference model (GCIM)** is a statistical method for detecting the causal direction of Genotype-by-environment(GxE) interaction studies 



 #### Authors: Zinabu Fentaw, Dovini Jayasinghe, S.Hong Lee


_<div align="justify">GCIM is a novel statistical method that extends beyond traditional polygenic risk scores-by-environment interactions (PRS×E) approaches.
    It systematically evaluates both the proposed and reverse causal directions, rather than relying
    on the assumptions of prior causal directions. By explicitly testing both directions, GCIM offers researchers
    data-driven insight into the likely causal directions of GxE interactions.</div>_

   
## 1. Package installation 
From GitHub 

~~~
library(devtools)
install_github("DoviniJ/GxEprs")
install_github("Zinabuf/GCIM")
~~~

## 2. Load the library

~~~
library(GxEprs)
library(GCIM)
~~~

## 3. Data Preparation
 
To ensure consistent and reliable estimation of G×E, the dataset should be split into two independent and non-overlapping subsets: a **discovery dataset (80%)** and a **target dataset (20%)**. This split should maintain consistency across genotype data, outcome variables, exposure variables, and potential confounders. 

#### Genotype data 
The genetic data must be in **PLINK binary format**, comprising three files: `.bed`, `.bim`, and `.fam`. 

**DummyData.fam:** This is a file associated with the PLINK binary format file, which contains the following columns in order. Please note that the file should not have column headings. This follows the PLINK .fam file format.

* family ID (FID)

* individual ID (IID)

* father's ID

* mother's ID

* sex

phenotype value

~~~
1 F1 I1  0  0  1 -9
2 F2 I2  0  0  2 -9
3 F3 I3  0  0  2 -9
4 F4 I4  0  0  1 -9
5 F5 I5  0  0  2 -9
6 F6 I6  0  0  1 -9
~~~

**DummyData.bim:**  This is a file associated with the PLINK binary format file, which contains the following columns in order. Please note that the file should not have column headings. This follows the PLINK .bim file format.

* chromosome code

* SNP ID

* position of centimorgans

* base-pair coordinate

* minor allele

* reference allele

~~~
1  1   rs1  0 1000  G  A
2  1  rs23  0 2000  G  A
3  1  rs45  0 3000  G  A
4  1  rs67  0 4000  G  A
5  1  rs89  0 5000  G  A
6  1 rs111  0 6000  G  A
~~~

**DummyData.bed:** This is the PLINK binary format file, which includes genotype information. This follows the PLINK .bed file format.
## 3.1. Proposed direction. 

 _<div align="justify">The **proposed causal direction** refers to the hypothesized G×E interaction in which the exposure affects the outcome, aligning with the researcher’s primary interest. The data should be carefully prepared and evaluated to test the causal direction in accordance with the researcher’s specified hypothesis of interest.</div>_
 
### 3.1.1. Input files format for the discovery dataset
**Note:** The **outcome file** should include `FID`, `IID`, and the `outcome`. For binary outcomes, follow standard coding conventions: use **PLINK’s default coding (1 = Control, 2 = Case)** in the **discovery dataset**.

#### 3.1.1.1. Genome-wide environment interaction study (GWEIS)

 Conducted a GWEIS to generate both additive and interaction polygenic risk scores (PRS). For analyses with a quantitative outcome, the input data must adhere to the same format required by the GxEprs framework 

**Qphen_disc.txt**: This is a .txt file containing the following columns in the specified order. Please note that the file should not have column headings. Therefore, the outcome file `Qphen_disc.txt` will have the following essential column:

* FID

* IID

* Outcome

 ~~~
1 F2 I2  1.6224168
2 F3 I3 -0.6784632
3 F4 I4  0.3788973
4 F5 I5 -0.4529113
5 F7 I7 -0.5511008
6 F8 I8  2.2036007
~~~

**Qexp_dis_cov.txt**: This is a .txt file containing the following columns in the specified order. Note that the file should have no column heading. The exposure variable and the covariate that are used to adjust the data frame, as expressed in GxEprs. 

* FID

* IID
  
*  Exposure variable 

* constant values.  Note: This is the input data format for GxEprs; if not specified, the model will omit the variable specified in the fourth column for the quantitative outcome. Keeping the fourth column constant is unnecessary when dealing with a binary outcome.

* Conf_1
  
  .
  .
  .
  
* Conf_n

(Note: These columns (Conf_1, ..., Conf_n) are optional. You can include any number of columns as confounders to adjust the GWEIS of the outcome.)

~~~
1 F2 I2 -0.38648517 0.1493707879 -0.2982989  0.06756177 -0.22387239 -0.11467347  -0.3169935
2 F3 I3  0.07919467 0.0062717957  0.1046483  0.16296653  0.23477644  0.17202480  0.3623103
3 F4 I4 -0.02190818 0.0004799682  0.1838144  0.30994543 -0.12036569 -0.30760305  0.2936090
4 F5 I5 -0.49888190 0.2488831453 -0.2969028 -0.18865130  0.26972598  0.23612115  -0.0617906
5 F7 I7  1.37810897 1.8991843340  0.3625800  0.01237302 -0.69589178 -0.07378856  -0.5883301
6 F8 I8 -1.39660300 1.9504999472 -0.8281950  0.31037620 -0.09381134  0.06069512   0.0726815
~~~

#### 3.1.1.2. Genome-wide association study (GWAS)

Perform a GWAS on the quantitative exposure phenotype to construct a PRS of exposure, adopting the same input data format required by the GxEprs framework. In this procedure, the exposure is treated as the outcome variable in the GWAS to obtain SNP effect estimates.

**Qexp_disc.txt:** This is a `.txt` file containing the following columns in the specified order. Please note that the file should not have column headings. Therefore, the exposure file `Qexp_disc.txt` will have the following essential column:

* FID

* IID

* Exposure

~~~
1 F2 I2 -0.38648517
2 F3 I3  0.07919467
3 F4 I4 -0.02190818
4 F5 I5 -0.49888190
5 F7 I7  1.37810897
6 F8 I8 -1.39660300
~~~

**Qcov_disc.txt:**  This covariate file is used to adjust the GWAS of the exposure variable. Note that it should not have column headings. Covariates for GWAS adjustment should be provided in a separate `.txt` file, which must include the following columns in the specified order for quantitative exposure. If the exposure variable is binary, adding constant values in the third and fourth columns is not required (as constant values can be removed), and any covariates may be included for adjustment.

* FID

* IID

* Constant_value

* Constant_value

  (Note: This is the input data format for GxEprs; if not specified, the model will omit the variable specified in the fourth column.)

  
* Conf_1
  
  .
  .
  .
  
* Conf_n

(Note: These columns (Conf_1, ..., Conf_n) are optional. You can include any number of columns as confounders to adjust the GWAS of the exposure phenotype.)

~~~
1 F2 I2  2  2 -0.2982989  0.06756177 -0.22387239 -0.11467347 -0.3169935
2 F3 I3  2  2  0.1046483  0.16296653  0.23477644  0.17202480  0.3623103
3 F4 I4  2  2  0.1838144  0.30994543 -0.12036569 -0.30760305  0.2936090
4 F5 I5  2  2 -0.2969028 -0.18865130  0.26972598  0.23612115 -0.0617906
5 F7 I7  2  2  0.3625800  0.01237302 -0.69589178 -0.07378856 -0.5883301
6 F8 I8  2  2 -0.8281950  0.31037620 -0.09381134  0.06069512  0.0726815
~~~

### 3.1.2. Input files format for the target dataset

The **outcome file** should include `FID`, `IID`, and the outcome.

**Qphen_tar.txt:** This is a '.txt` file which contains the following columns in order. Please note that the file should not have column headings.

* FID

* IID

* Outcome
  
~~~
1  F1  I1  0.4431358
2  F6  I6  0.6097717
3  F9  I9 -0.1947686
4 F14 I14 -0.4497752
5 F18 I18 -0.3722471
6 F34 I34 -0.4311244
~~~

**Qexp_tar_cov.txt:** This file should contain the exposure variable along with covariates for adjustment in the target dataset, provided separately from other inputs. Please note that the file should not have column headings. 

* FID

* IID

* Exposure 
  
* Conf_1
  
  .
  .
  .
  
* Conf_n

(Note: These columns (Conf_1, ..., Conf_n) are optional. You can include any number of columns as confounders to adjust the GxE in GCIM for the outcome.)
~~~
1  F1  I1 -0.89543441  0.3950992 -0.073442880  0.4896785 -0.2800233  0.23709983
2  F6  I6  0.03913846 -0.2444538 -0.225332496 -0.1588353  0.6133494 -0.49531486
3  F9  I9  1.04309085  0.4162231  0.112492137  0.0508843 -0.2630248 -0.01771196
4 F14 I14 -0.03991198  0.2093967 -0.218687107 -0.0126684  0.3663265  0.48036948
5 F18 I18 -0.48303534 -0.2133331 -0.003787907  0.1960149 -0.1254212 -0.49177427
6 F34 I34  0.15241566  0.2894822  0.564673054 -0.1132087 -0.5677780 -0.30682643
~~~

## 3.2. Reverse causal direction. 

 _<div align="justify">The **reverse causal direction** test assesses the validity of the assumed causal relationship by switching the roles of exposure and outcome. Following the primary analysis in the proposed causal direction, the reverse test is conducted to verify the possible causal direction, with the original outcome treated as the exposure and the original exposure treated as the outcome. Consistency in data structure and formatting should be maintained across analyses.</div>_
 
 
 ### 3.2.1. Input files format for the discovery dataset

 #### 3.2.1.1. GWEIS
 **Qexp_disc.txt:** This is a .txt file containing the following columns in the specified order. The exposure is an outcome for the reverse directions. Please note that the file should not have column headings. To generate both the additive and interaction polygenic risk scores (PRS), we performed a genome-wide environment interaction study (GWEIS) using the GxEprs data framework. When conducting a GWEIS with a quantitative outcome, the input data must follow the same format as required for the GxEprs framework. For reproducibility, the exposure variable should be organized in a dedicated file, for example: 
* FID

* IID

* Exposure

~~~
1 F2 I2 -0.38648517
2 F3 I3  0.07919467
3 F4 I4 -0.02190818
4 F5 I5 -0.49888190
5 F7 I7  1.37810897
6 F8 I8 -1.39660300
~~~
  
~~~

~~~

**Qphen_disc_cov.txt**: This is a `.txt` file containing the following columns in the specified order. Note that the file should not have a column heading. The exposure variable and the covariate that are used to adjust the data frame, as expressed in GxEprs.

* FID

* IID

* Outcome

* Constant values.
 Note: This is the input data format for GxEprs; if not specified, the model will omit the variable specified in the fourth column for the quantitative outcome. Keeping the fourth column constant is unnecessary when dealing with a binary outcome.

* Conf_1
  
  .
  .
  .
  
* Conf_n

(Note: These columns(Conf_1, ..., Conf_n)  are optional. You can include any number of columns as confounders to adjust the GWAS of the exposure phenotype.)

~~~
1  F10  I10  3.2536347 10.58613897  0.03008962 -0.05255265 -0.49929570  0.04939406 -0.01633697
2 F100 I100 -0.9757257  0.95204069  0.08582469  0.47371147  0.16039529  -0.21381246  0.18741603
3 F101 I101 -0.2124444  0.04513264  0.17514805 -0.12537552 -0.02494976  -0.19626511  0.11283626
4 F102 I102 -0.1474361  0.02173741 -0.73956266  0.34143490 -0.41246477  -0.58502939  0.19016286
5 F103 I103  0.2076099  0.04310188 -0.13795823  0.48665443 -0.35363273  0.52328458 -0.04537773
6 F104 I104  0.1599443  0.02558217 -0.10844172  0.21432447  0.34433911  0.60188311  0.75556868
~~~

3.2.1.2. GWAS
 
**Qphen_disc.txt:** This is a `.txt` file containing the following columns in the specified order. The exposure as an outcome for the reverse directions. Please note that the file should not have column headings. To generate both the additive and interaction polygenic risk scores (PRS), we performed a genome-wide environment interaction study (GWEIS) using the GxEprs data framework. When conducting a GWEIS with a quantitative outcome, the input data must follow the same format as required for the GxEprs framework. For reproducibility, the outcome data should be organized in a dedicated file, for example: FID IID quantitative outcome. 

* FID

* IID

* Outcome

~~~
1 F2 I2  1.6224168
2 F3 I3 -0.6784632
3 F4 I4  0.3788973
4 F5 I5 -0.4529113
5 F7 I7 -0.5511008
6 F8 I8  2.2036007
~~~
**Qcov_disc.txt**: This is a `.txt` file containing the following columns in the specified order. The discovery dataset has 800 individuals. Note that the file should not have a column heading. The exposure variable and the covariate that are used to adjust the data frame, as expressed in GxEprs. 

* FID

* IID

* Constant_values

* Constant_values

  (Note: This is the input data format for GxEprs; if not specified, the model will omit the variable specified in the fourth column.)

  
* Conf_1
  
  .
  .
  .
  
* Conf_n

(Note: These columns(Covf_1, ..., Conf_n) are optional. You can include any number of columns as confounders to adjust the GWAS of the exposure phenotype.)

~~~
1 F2 I2  2  2 -0.2982989  0.06756177 -0.22387239 -0.11467347 -0.3169935
2 F3 I3  2  2  0.1046483  0.16296653  0.23477644  0.17202480  0.3623103
3 F4 I4  2  2  0.1838144  0.30994543 -0.12036569 -0.30760305  0.2936090
4 F5 I5  2  2 -0.2969028 -0.18865130  0.26972598  0.23612115 -0.0617906
5 F7 I7  2  2  0.3625800  0.01237302 -0.69589178 -0.07378856 -0.5883301
6 F8 I8  2  2 -0.8281950  0.31037620 -0.09381134  0.06069512  0.0726815
~~~

### 3.2.2. Input files format for the target dataset
**Qexp_tar.txt:** This is a `.txt` file which contains the following columns in order. The target dataset consists of individuals who are independent of those in the discovery dataset. Please note that the file should not have column headings.

* FID
* IID
* Exposure
  
~~~    
1  F1  I1 -0.89543441
2  F6  I6  0.03913846
3  F9  I9  1.04309085
4 F14 I14 -0.03991198
5 F18 I18 -0.48303534
6 F34 I34  0.15241566
~~~

**Qphen_tar_cov.txt:** The exposure variable and other covariates used for adjustment should be from the target dataset and provided in a separate file. Please note that the file should not have column headings.

* FID

* IID

* Outcome 
  
* Conf_1
  
  .
  .
  .
  
* Conf_n

(Note: These columns (Conf_1, ..., Conf_n) are optional. You can include any number of columns as confounders to adjust the outcome.)
~~~
1    F1    I1  0.4431358  0.3950992 -0.07344288  0.4896785 -0.280023347  0.23709983
2 F1000 I1000  0.9743923 -0.1573689  0.32655712  0.1143439  0.284942855  -0.62625272
3  F111  I111 -0.8706305  0.3494551 -0.68904013 -0.4502063  0.001334553  -0.36046152
4  F122  I122 -0.6352419 -0.6021410  0.07734350  0.2992616  0.493944593  -0.11699865
5  F125  I125 -0.2725707  0.8709892  0.67394957  0.1530370 -0.148748060  0.26738235
6  F127  I127 -0.7433867 -0.3195486  0.38702170  0.1394358 -0.103527076  -0.04755472
~~~

## 4. Analysis workflow

GCIM analyses use PLink2 to analyze discovery data. 
1. Download plink2 software from the [Plink](https://www.cog-genomics.org/plink/2.0/) website and then specify the executable Plink software path.
   
~~~
plink_path <- "<plink_path>/plink2"
~~~ 
### 4.1. proposed causal direction 
#### 4.1.1. Quantitative outcome 
   
 ~~~
# Load required libraries
library(GxEprs)
library(GCIM)
# Set plink path
plink_path <- "<plink_path>/plink2"

# For quantitative traits, use corresponding functions
# We conducted a GWEIS of the outcome variable to generate both additive and interaction PRS.
GWEIS <- GWEIS_quantitative(plink_path, "DummyData", "Qphen_disc.txt", "Qexp_dis_cov.txt")
# We conducted a GWAS of the exposure variable to construct an exposure PRS.
GWAS <- GWAS_quantitative(plink_path, "DummyData", "Qexp_disc.txt", "Qcov_disc.txt")

# Extract and compute PRS
# Extract summary statistics, such as extracting the additive and interaction components from the GWEIS for the outcome variable. 
add <- GWEIS[c("ID", "A1", "ADD_BETA")]
int <- GWEIS[c("ID", "A1", "INTERACTION_BETA")]
# Extracting the additive component from the exposure GWAS
add_exp <- GWAS[c("ID", "A1", "BETA")]

add_prs <- PRS_quantitative(plink_path, "DummyData", summary_input = add)
int_prs <- PRS_quantitative(plink_path, "DummyData", summary_input = int)
#Note: If the exposure variable is quantitative, use it as is. If the exposure variable is binary, apply the `PRS_binary` function to generate the object `exp_prs` specified as described below. 
exp_prs <- PRS_quantitative(plink_path, "DummyData", summary_input = add_exp)

# Run GCIM analysis with automatic saving and scaling
# This model specification corresponds to Model 4, as presented in the manuscript.
result1 <- gcim_q0("Qphen_tar.txt", "Qexp_tar_cov.txt", 
                  Add_PRS = add_prs, Int_PRS = int_prs, Cov_PRS = exp_prs)
#This model additionally incorporates the main effect of the PRS for the exposure variable. Details of this specification are provided in the model equation section of the output. 
 result2 <- gcim_q1("Qphen_tar.txt", "Qexp_tar_cov.txt", 
                  Add_PRS = add_prs, Int_PRS = int_prs, Cov_PRS = exp_prs)
 # result

~~~

~~~
 print(result1$model_summary)
~~~

~~~
Coefficients:
                Estimate Std. Error t value Pr(>|t|)
(Intercept)      0.52436    0.09644   5.437 1.65e-07 ***
Add_PRS          0.23697    0.13890   1.706   0.0896 .
Int_PRS         -0.06638    0.14157  -0.469   0.6397
Covariate_Pheno  0.03804    0.08826   0.431   0.6669
Conf_1           0.36873    0.29359   1.256   0.2107
Conf_2          -0.09726    0.28749  -0.338   0.7355
Conf_3           0.30371    0.28008   1.084   0.2796
Conf_4          -0.21463    0.30101  -0.713   0.4767
Conf_5           0.52926    0.30496   1.735   0.0843 .
Int_PRS:Cov_PRS  0.38337    0.17589   2.180   0.0305 *
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
~~~

~~~
 print(result2$model_summary)
~~~

~~~
Coefficients:
                Estimate Std. Error t value Pr(>|t|)
(Intercept)      0.52731    0.09687   5.444 1.61e-07 ***
Add_PRS          0.23486    0.13927   1.686   0.0934 .
Int_PRS         -0.04416    0.15015  -0.294   0.7690
Covariate_Pheno  0.05614    0.09710   0.578   0.5638
Cov_PRS         -0.07143    0.15814  -0.452   0.6520
Conf_1           0.36499    0.29432   1.240   0.2165
Conf_2          -0.10475    0.28857  -0.363   0.7170
Conf_3           0.29173    0.28192   1.035   0.3021
Conf_4          -0.20874    0.30192  -0.691   0.4902
Conf_5           0.53095    0.30563   1.737   0.0840 .
Int_PRS:Cov_PRS  0.37272    0.17784   2.096   0.0374 *
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
~~~
Note:the output displayed in the analysis is


`Add_PRS`: additive PRS
`Int_PRS`: Interaction PRS
`Covariate_Pheno`: Exposure variable
`Cov_PRS` : Exposure PRS
`Conf_1, Conf_2, …, Conf_14` are based on the column order in the input file.
`Int_PRS:Cov_PRS`: GxE interaction 


#### 4.1.2. Binary outcome 

~~~
# Load required libraries
library(GxEprs)
library(GCIM)
# Set plink path
plink_path <- "<plink_path>/plink2"
# Run GxEprs analysis for binary traits
#  We conducted a GWEIS of the outcome variable to generate both additive and interaction PRS.
GWEIS <- GWEIS_binary(plink_path, "DummyData", "Bphen_disc.txt", "Bexp_disc_cov.txt")
# Note: If the exposure variable is binary, use it as is. If the exposure variable is quantitative, apply the `GWAS_quantitative` function to generate the `GWAS` specified as described below.
# We conducted a GWAS of the exposure variable to construct an exposure PRS.
GWAS <- GWAS_binary(plink_path, "DummyData", "Bexp_disc.txt", "Bcov_disc.txt")

# Extract summary statistics
# Extracting the additive and interaction components from the GWEIS for the outcome variable. 
add <- GWEIS[c("ID", "A1", "ADD_BETA")]
int <- GWEIS[c("ID", "A1", "INTERACTION_BETA")]
# Extracting the additive component from the exposure GWAS
add_exp <- GWAS[c("ID", "A1", "BETA")]

# Compute PRS for each component
add_prs <- PRS_binary(plink_path, "DummyData", summary_input = add)  # Additive PRS
int_prs <- PRS_binary(plink_path, "DummyData", summary_input = int)  # Interaction PRS
#Note: If the exposure variable is binary, use it as is. If the exposure variable is quantitative, apply the `PRS_quantitative` function to generate the object `exp_prs` specified as described below.
exp_prs <- PRS_binary(plink_path, "DummyData", summary_input = add_exp)  # Covariate PRS

# Run GCIM analysis with automatic saving and scaling
# A similar model specification is applied as described above for the quantitative analysis.
 result1 <- gcim_b0("Bphen_tar.txt", "Bexp_tar_cov.txt", 
                  Add_PRS = add_prs, Int_PRS = int_prs, Cov_PRS = exp_prs) 
 result2 <- gcim_b1("Bphen_tar.txt", "Bexp_tar_cov.txt", 
                  Add_PRS = add_prs, Int_PRS = int_prs, Cov_PRS = exp_prs)
# results
~~~

~~~
 print(result1$model_summary)
~~~

~~~
Coefficients:
                 Estimate Std. Error z value Pr(>|z|)
(Intercept)     -0.959830   0.171843  -5.586 2.33e-08 ***
Add_PRS          1.451542   0.656026   2.213  0.02692 *
Int_PRS          1.411390   0.643010   2.195  0.02817 *
Covariate_Pheno -0.002067   0.155615  -0.013  0.98940
Conf_1           0.740688   0.526348   1.407  0.15936
Conf_2           0.189791   0.510838   0.372  0.71024
Conf_3           0.099874   0.490649   0.204  0.83870
Conf_4          -1.714645   0.581453  -2.949  0.00319 **
Conf_5           0.544927   0.558394   0.976  0.32912
Int_PRS:Cov_PRS  0.185485   0.306038   0.606  0.54446
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
~~~

~~~
 print(result2$model_summary)
~~~

~~~
Coefficients:
                Estimate Std. Error z value Pr(>|z|)
(Intercept)     -0.97533    0.17415  -5.601 2.14e-08 ***
Add_PRS          1.43645    0.65450   2.195  0.02818 *
Int_PRS          1.42451    0.64366   2.213  0.02689 *
Covariate_Pheno  0.06795    0.17074   0.398  0.69066
Cov_PRS         -0.27794    0.27620  -1.006  0.31428
Conf_1           0.71207    0.52764   1.350  0.17716
Conf_2           0.16422    0.51130   0.321  0.74807
Conf_3           0.07750    0.49195   0.158  0.87482
Conf_4          -1.72929    0.58247  -2.969  0.00299 **
Conf_5           0.52743    0.56350   0.936  0.34927
Int_PRS:Cov_PRS  0.18738    0.30997   0.605  0.54550
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
~~~

 ### 4.2. Reverse causal direction
 
To evaluate the **reverse causal direction**, re-analyze the same dataset by switching the roles of the exposure and outcome variables. This means treating the previously defined outcome variable as the new exposure, and the previous exposure variables as the new outcome. Rearrange the data using the same structure and formatting approach used for the proposed causal directions as mentioned above, ensuring consistency across analysis pipeline. The only difference should be the reassignment of variable roles.

#### 4.2.1. Quantitative outcome

~~~
# Load required libraries
library(GxEprs)
library(GCIM)
# Set plink path
plink_path <- "<plink_path>/plink2"

# For quantitative traits, use corresponding functions
# We conducted a GWEIS of the outcome variable (treated as the exposure in the proposed causal direction) to generate both additive and interaction PRS.
GWEIS <- GWEIS_quantitative(plink_path, "DummyData", "Qexp_disc.txt", "Qphen_disc_cov.txt")
# We conducted a GWAS of the exposure variable (treated as the outcome in the proposed causal directions) to construct an exposure PRS
GWAS <- GWAS_quantitative(plink_path, "DummyData", "Qphen_disc.txt", "Qcov_disc.txt")

# Extract and compute PRS

add_exp <- GWEIS[c("ID", "A1", "ADD_BETA")]
int_exp <- GWEIS[c("ID", "A1", "INTERACTION_BETA")]

add_out <- GWAS[c("ID", "A1", "BETA")]

add_prs <- PRS_quantitative(plink_path, "DummyData", summary_input = add_exp)
int_prs <- PRS_quantitative(plink_path, "DummyData", summary_input = int_exp) 
out_prs <- PRS_quantitative(plink_path, "DummyData", summary_input = add_out)

# Run GCIM analysis with automatic saving and scaling
 result1 <- gcim_q0("Qexp_tar.txt", "Qphen_tar_cov.txt", 
                  Add_PRS = add_prs, Int_PRS = int_prs, Cov_PRS = out_prs)
 result2 <- gcim_q1("Qexp_tar.txt", "Qphen_tar_cov.txt", 
                  Add_PRS = add_prs, Int_PRS = int_prs, Cov_PRS = out_prs)
 # results
~~~

~~~
 print(result1$model_summary)
~~~

~~~
Coefficients:
                Estimate Std. Error t value Pr(>|t|)
(Intercept)      0.03861    0.07591   0.509 0.611609
Add_PRS          0.60882    0.09717   6.266 2.43e-09 ***
Int_PRS          0.34420    0.09745   3.532 0.000518 ***
Covariate_Pheno  0.03219    0.05330   0.604 0.546566
Conf_1           0.08080    0.21877   0.369 0.712284
Conf_2           0.01577    0.21526   0.073 0.941661
Conf_3          -0.04147    0.20821  -0.199 0.842337
Conf_4          -0.08103    0.22802  -0.355 0.722726
Conf_5          -0.07733    0.23000  -0.336 0.737059
Int_PRS:Cov_PRS  0.10616    0.12384   0.857 0.392391
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
~~~

~~~
 print(result2$model_summary)
~~~

~~~
Coefficients:
                Estimate Std. Error t value Pr(>|t|)
(Intercept)      0.03848    0.07602   0.506 0.613354
Add_PRS          0.60577    0.09742   6.218 3.15e-09 ***
Int_PRS          0.35804    0.09977   3.589 0.000423 ***
Covariate_Pheno  0.03570    0.05363   0.666 0.506410
Cov_PRS         -0.06836    0.10243  -0.667 0.505327
Conf_1           0.07492    0.21927   0.342 0.732961
Conf_2           0.00596    0.21608   0.028 0.978026
Conf_3          -0.02998    0.20922  -0.143 0.886218
Conf_4          -0.08915    0.22868  -0.390 0.697104
Conf_5          -0.08013    0.23037  -0.348 0.728367
Int_PRS:Cov_PRS  0.11754    0.12519   0.939 0.348987
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
~~~

#### 4.2.2 Binary outcome 

~~~
# Load required libraries
library(GxEprs)
library(GCIM)
# Set plink path
plink_path <- "<plink_path>/plink2"
# Run GxEprs analysis for binary traits

# We conducted a GWEIS of the outcome variable (treated as the exposure in the proposed causal direction) to generate both additive and interaction PRS.
GWEIS <- GWEIS_binary(plink_path, "DummyData", "Bexp_disc.txt", "Bphen_disc_cov.txt")
# We conducted a GWAS of the exposure variable (treated as the outcome in the proposed causal directions) to construct an exposure PRS
GWAS <- GWAS_binary(plink_path, "DummyData", "Bphen_disc.txt", "Bcov_disc.txt")

# Extract summary statistics
add_exp <- GWEIS[c("ID", "A1", "ADD_BETA")]
int_exp <- GWEIS[c("ID", "A1", "INTERACTION_BETA")]
add_out <- GWAS[c("ID", "A1", "BETA")]

# Compute PRS for each component
add_prs <- PRS_binary(plink_path, "DummyData", summary_input = add_exp)  # Additive PRS
int_prs <- PRS_binary(plink_path, "DummyData", summary_input = int_exp)  # Interaction PRS
out_prs <- PRS_binary(plink_path, "DummyData", summary_input = add_out)  # Covariate PRS

# Run GCIM analysis with automatic saving and scaling
# A similar model specification is applied as described above for the quantitative analysis.
 result1 <- gcim_b0("Bexp_tar.txt", "Bphen_tar_cov.txt", 
                  Add_PRS = add_prs, Int_PRS = int_prs, Cov_PRS = out_prs) 
 result2 <- gcim_b1("Bexp_tar.txt", "Bphen_tar_cov.txt", 
                  Add_PRS = add_prs, Int_PRS = int_prs, Cov_PRS = out_prs)
# results
~~~

~~~
print(result1$model_summary)
~~~

~~~
Coefficients:
                Estimate Std. Error z value Pr(>|z|)
(Intercept)     -1.15264    0.24701  -4.666 3.06e-06 ***
Add_PRS          2.06805    0.63443   3.260  0.00112 **
Int_PRS          1.94608    0.64951   2.996  0.00273 **
Covariate_Pheno  0.75477    0.35533   2.124  0.03366 *
Conf_1          -0.22743    0.53137  -0.428  0.66865
Conf_2           0.01202    0.51909   0.023  0.98153
Conf_3          -0.47566    0.51392  -0.926  0.35468
Conf_4          -0.59997    0.56149  -1.069  0.28529
Conf_5           0.23844    0.54946   0.434  0.66433
Int_PRS:Cov_PRS -0.03820    0.22323  -0.171  0.86414
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
~~~

~~~
print(result2$model_summary)
~~~

~~~
Coefficients: 
                Estimate Std. Error z value Pr(>|z|)
(Intercept)     -1.15264    0.24701  -4.666 3.06e-06 ***
Add_PRS          2.06805    0.63443   3.260  0.00112 **
Int_PRS          1.94608    0.64951   2.996  0.00273 **
Covariate_Pheno  0.75477    0.35533   2.124  0.03366 *
Conf_1          -0.22743    0.53137  -0.428  0.66865
Conf_2           0.01202    0.51909   0.023  0.98153
Conf_3          -0.47566    0.51392  -0.926  0.35468
Conf_4          -0.59997    0.56149  -1.069  0.28529
Conf_5           0.23844    0.54946   0.434  0.66433
Int_PRS:Cov_PRS -0.03820    0.22323  -0.171  0.86414
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
~~~
 
### Contact

For any inquiries or support, please feel free to contact me at:
**[zfentaw5@gmail.com](mailto:zfentaw5@gmail.com)**
