---
# <h1 align="center">GCIM</h1>
---

The **genetic causality inference model (GCIM)** is a statistical method for detecting the causal direction of Genotype-by-environment(GxE) interaction studies 



 #### Authors: Zinabu Fentaw, Dovini Jayasinghe, S.Hong Lee


_<div align="justify">GCIM is a novel statistical method that extends beyond traditional polygenic risk scores-by-environment interactions (PRS×E) approaches.
    It systematically evaluates both the proposed and reverse causal directions, rather than relying
    on the assumptions of prior causal directions. By explicitly testing both directions, GCIM offers researchers
    data-driven insight into the likely causal directions of GxE interactions.</div>_

   
## I. Package installation 
From GitHub 

~~~
library(devtools)
install_github("DoviniJ/GxEprs")
install_github("Zinabuf/GCIM")
~~~

## II. Load the library

~~~
library(GxEprs)
library(GCIM)
~~~

## III. Data Preparation
 
To ensure consistent and reliable estimation of G×E, the dataset should be split into two independent and non-overlapping subsets: a **discovery dataset (80%)** and a **target dataset (20%)**. This split should maintain consistency across genotype data, outcome variables, exposure variables, and potential confounders. 

#### Genotype data 
The genetic data must be in **PLINK binary format**, comprising three files: `.bed`, `.bim`, and `.fam`. 

**DummyData.fam:** This is a file associated with the PLINK binary format file, which contains the following columns in order. Please note that the file does not have column headings. This follows the PLINK .fam file format.

family ID (FID)

individual ID (IID)

father's ID

mother's ID

sex

phenotype value

~~~
1 ID_1 ID_1  0  0  1 -9
2 ID_2 ID_2  0  0  2 -9
3 ID_3 ID_3  0  0  2 -9
4 ID_4 ID_4  0  0  2 -9
5 ID_5 ID_5  0  0  1 -9
6 ID_6 ID_6  0  0  1 -9
~~~

**DummyData.bim:**  This is a file associated with the PLINK binary format file, which contains the following columns in order. The example dataset has 1,000 SNPs. Please note that the file does not have column headings. This follows the PLINK .bim file format

chromosome code

SNP ID

position of centimorgans

base-pair coordinate

minor allele

reference allele

~~~
1  1 SNP_1  0  768448  A  G
2  1 SNP_2  0  853954  C  A
3  1 SNP_4  0  940203  A  G
4  1 SNP_5  0  987670  T  G
5  1 SNP_6  0 1021695  A  G
6  1 SNP_7  0 1053452  A  G
~~~

**DummyData.bed:** This is the PLINK binary format file, which includes genotype information. This follows the PLINK .bed file format.
## 1. Proposed direction. 

 _<div align="justify">The **proposed causal direction** refers to the hypothesized G×E interaction in which the exposure affects the outcome, aligning with the researcher’s primary interest. The data should be carefully prepared and evaluated to test the causal direction in accordance with the researcher’s specified hypothesis of interest.</div>_
 
### 1.1. Discovery input files
**Note:** The **outcome file** should include `FID`, `IID`, and the `outcome variable`. For binary outcomes, follow standard coding conventions: use **PLINK’s default coding (1 = Control, 2 = Case)** in the **discovery dataset**.

#### 1.1.1. Genome-wide environment interaction study (GWEIS)

 Conducted a GWEIS to generate both additive and interaction polygenic risk scores (PRS). For analyses with a quantitative outcome, the input data must adhere to the same format required by the GxEprs framework 

**Qphen_disc.txt**: This is a .txt file containing the following columns in the specified order. Please note that the file does not have column headings. Therefore, the outcome file `Qphen_disc.txt` will have the following essential column:

FID

IID

Outcome

 ~~~
FID   IID    Outcome
1 ID_1 ID_1 31.6534
2 ID_2 ID_2 25.5035
3 ID_3 ID_3 26.7391
4 ID_4 ID_4 25.5271
5 ID_5 ID_5 26.7165
6 ID_6 ID_6 38.8272
~~~

**Qexp_dis_cov.txt**: This is a .txt file containing the following columns in the specified order. Note that the file has no column heading. The exposure variable and the covariate that are used to adjust the data frame, as expressed in GxEprs. 

FID

IID

standardized Exposure

constant values. Note: This is the input data format for GxEprs; if not specified, the mode will omit the variable specified in the fourth column for quantitative exposure. Keeping the fourth column constant is unnecessary when dealing with a binary outcome.
14 confounders of the discovery sample (Note: These columns are optional. Can use any number of columns as confounders to adjust the phenotype upon user requirement.)

~~~
   FID IID   Exposure   constant_value Conf_1 Conf_2 Conf_3 Conf_4  Conf_5 Conf_6
1 ID_1 ID_1 -0.64402046  2 -3.831420 64 -14.03640 5.517420  0.0714337  5.662630
2 ID_2 ID_2 -0.02786981  2  0.614044 66 -10.85050 2.119980 -0.8828830 -0.441662
3 ID_3 ID_3  2.12865748  2 -0.237792 55  -9.75369 3.183430 -2.0979300  6.873450
4 ID_4 ID_4  2.12865748  2  6.698660 47  -9.07045 0.956878 -2.4840700  1.063590
5 ID_5 ID_5 -0.95209579  2 -1.614230 59 -12.93790 1.294610 -1.7997300  1.444040
6 ID_6 ID_6 -0.02786981  2 -4.389270 52 -11.85160 0.888978 -2.7231000  1.116810
      Conf_7  Conf_8    Conf_9    Conf_10     Conf_11   Conf_12  Conf_13  Conf_14
1  0.865562 -2.269570 -0.09658590 -2.354970  1.0588900  0.195302   0   7
2 -2.641770  2.789440  0.52458600  2.671340 -2.6372400 -0.998764   1  20
3 11.377700  2.969610 -1.11879000  0.873649  3.3552300 -4.578310   1  10
4 -3.132470  2.123200 -0.00976751  0.820582  0.0305345  1.630300   1  20
5 -6.828980 -2.967950 -2.91577000 -1.828810  7.1589200  2.109160   1  20
6 -3.646760 -0.594538 -1.75430000 -0.716014 -2.3906700  1.312950   1  10
~~~

#### 1.1.2. Genome-wide environment study (GWAS)

Perform a GWAS on the quantitative exposure phenotype to construct a PRS of exposure, adopting the same input data format required by the GxEprs framework. In this procedure, the exposure is treated as the outcome variable in the GWAS to obtain SNP effect estimates.

**Qexp_disc.txt:** This is a `.txt` file containing the following columns in the specified order. Please note that the file does not have column headings. Therefore, the exposure file `Qexp_disc.txt` will have the following essential column:

FID

IID

Exposure

~~~
  FID  IID   Exposure
1 ID_1 ID_1 -0.64402046
2 ID_2 ID_2 -0.02786981
3 ID_3 ID_3  2.12865748
4 ID_4 ID_4  2.12865748
5 ID_5 ID_5 -0.95209579
6 ID_6 ID_6 -0.02786981
~~~

**Qcov_disc.txt:**  This covariate file is used to adjust the GWAS of the exposure variable. Note that it does not contain column headings. Covariates for GWAS adjustment should be provided in a separate `.txt` file, which must include the following columns in the specified order for quantitative exposure. If the exposure is binary, the third and fourth columns are not required (as constant values can be removed), and any covariates may be included for adjustment.

FID

IID

Constant value

Constant values (Note: This is the input data format for GxEprs; if not specified, the mode will omit the variable specified in the fourth column.)
14 confounders of the discovery sample (Note: These columns are optional. Can use any number of columns as confounders to adjust the phenotype upon user requirement.)

~~~
   FID  IID  Constant_value Constant_value Conf_1  Conf_2  Conf_3  Conf_4   Conf_5 Conf_6
1 ID_1 ID_1  2  2 -3.831420 64 -14.03640 5.517420  0.0714337  5.662630
2 ID_2 ID_2  2  2  0.614044 66 -10.85050 2.119980 -0.8828830 -0.441662
3 ID_3 ID_3  2  2 -0.237792 55  -9.75369 3.183430 -2.0979300  6.873450
4 ID_4 ID_4  2  2  6.698660 47  -9.07045 0.956878 -2.4840700  1.063590
5 ID_5 ID_5  2  2 -1.614230 59 -12.93790 1.294610 -1.7997300  1.444040
6 ID_6 ID_6  2  2 -4.389270 52 -11.85160 0.888978 -2.7231000  1.116810
   Conf_7    Conf_8     Conf_9     Conf_10    Conf_11    Conf_12 Conf_13 Conf_14
1  0.865562 -2.269570 -0.09658590 -2.354970  1.0588900  0.195302   0   7
2 -2.641770  2.789440  0.52458600  2.671340 -2.6372400 -0.998764   1  20
3 11.377700  2.969610 -1.11879000  0.873649  3.3552300 -4.578310   1  10
4 -3.132470  2.123200 -0.00976751  0.820582  0.0305345  1.630300   1  20
5 -6.828980 -2.967950 -2.91577000 -1.828810  7.1589200  2.109160   1  20
6 -3.646760 -0.594538 -1.75430000 -0.716014 -2.3906700  1.312950   1  10

~~~

### 1.2. Target input files

The **outcome file** should include `FID`, `IID`, and the outcome variable. For binary outcomes, use **binary coding (0 = Control, 1 = Case)** in the **target dataset**.

**Qphen_tar.txt:** This is a '.txt` file which contains the following columns in order. Please note that the file does not have column headings.

~~~
   FID     IID  Outcome
1 ID_801 ID_801 26.5723
2 ID_802 ID_802 20.2632
3 ID_803 ID_803 27.7365
4 ID_804 ID_804 18.7500
5 ID_805 ID_805 23.3025
6 ID_806 ID_806 24.9871
~~~

**Qexp_tar_cov.txt:** This file should contain the exposure variable along with covariates for adjustment in the target dataset, provided separately from other inputs.

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

## 2. Reverse causal direction. 

 _<div align="justify">The **reverse causal direction** test assesses the validity of the assumed causal relationship by switching the roles of exposure and outcome. Following the primary analysis in the proposed causal direction, the reverse test is conducted to verify the possible causal direction, with the original outcome treated as the exposure and the original exposure treated as the outcome. Consistency in data structure and formatting should be maintained across analyses.</div>_
 
 ### 2.1.. Discovery input files
 
**Qexp_disc.txt:** This is a `.txt` file containing the following columns in the specified order. The exposure as an outcome for the reverse directions. Please note that the file does not have column headings. To generate both the additive and interaction polygenic risk scores (PRS), we performed a genome-wide environment interaction study (GWEIS) using the GxEprs data framework. When conducting a GWEIS with a quantitative outcome, the input data must follow the same format as required for the GxEprs framework. For reproducibility, the outcome data should be organized in a dedicated file, for example: FID IID quantitative outcome. 

~~~
  FID   IID   Exposure
1 ID_1 ID_1 -0.64402046
2 ID_2 ID_2 -0.02786981
3 ID_3 ID_3  2.12865748
4 ID_4 ID_4  2.12865748
5 ID_5 ID_5 -0.95209579
6 ID_6 ID_6 -0.02786981
~~~
**Qphen_disc_cov.txt**: This is a `.txt` file containing the following columns in the specified order. The discovery dataset has 800 individuals. Note that the file does not have a column heading. The exposure variable and the covariate that are used to adjust the data frame, as expressed in GxEprs. 
FID
IID
standardized covariate
constant values (Note: This is the input data format for GxEprs; if not specified, the mode will omit the variable specified in the fourth column.)
14 confounders of the discovery sample (Note: These columns are optional. Can use any number of columns as confounders to adjust the phenotype upon user requirement.)

~~~
    FID    FID  Outcome Constant_value Conf_1 Conf_2 Conf_3  Conf_4  Conf_5 Conf_6
1   ID_1   ID_1 31.6534  2 -3.831420 64 -14.0364 5.51742  0.0714337  5.662630
2  ID_10  ID_10 34.1878  2 -3.408730 45 -13.1770 4.96769 -0.4055640 -1.274350
3 ID_100 ID_100 23.1237  2 -1.691620 60 -14.2821 3.01444 -0.0835603  1.622540
4 ID_101 ID_101 39.1574  2 -3.979340 52 -12.2061 5.93377 -0.4792830  3.473790
5 ID_102 ID_102 29.0466  2 -2.377090 69 -13.6426 4.53418 -2.2062600  0.222482
6 ID_103 ID_103 28.5813  2  0.268657 67 -14.7840 5.09824 -3.5448200  5.679140
   Conf_7    Conf_8    Conf_9   Conf_10   Conf_11    Conf_12  Conf_13 Conf_14
1  0.865562 -2.26957 -0.0965859 -2.354970  1.058890  0.195302   0   7
2 -9.632180 -1.74190  2.6967400 -1.321190  0.742495 -0.688296   1  13
3 -8.066920 -1.65889 -2.7431200 -0.944793  7.769740  0.140430   1  10
4 -3.995170  1.17464  2.9239000 -2.148990 -7.193830  2.918810   0  20
5 -3.869480  0.72570  0.3461710 -0.882131 -2.526320 -1.620050   1  20
6 10.636900 -2.92550  0.0215072 -1.328180 -0.774790  1.984740   1  20
~~~

**Qphen_disc.txt:** To construct PRS for the exposure variable, we first performed a GWAS on the quantitative exposure phenotype, adopting the same input data format required by the GxEprs framework. In this procedure, the exposure is treated as the outcome variable in the GWAS to obtain SNP effect estimates.

~~~
      FID   IID  Outcome
1   ID_1   ID_1 31.6534
2  ID_10  ID_10 34.1878
3 ID_100 ID_100 23.1237
4 ID_101 ID_101 39.1574
5 ID_102 ID_102 29.0466
6 ID_103 ID_103 28.5813
~~~
**Qexp_disc.txt:** To construct PRS for the exposure variable, we first performed a GWAS on the quantitative exposure phenotype, adopting the same input data format required by the GxEprs framework. In this procedure, the exposure is treated as the outcome variable in the GWAS to obtain SNP effect estimates.

~~~
      FID   IID  Constant_values Constant_value Conf_1  Conf_2  Conf_3  Conf_4 Conf_5 Conf_6
1   ID_1   ID_1  2  2 -3.831420 64 -14.0364 5.51742  0.0714337  5.662630
2  ID_10  ID_10  2  2 -3.408730 45 -13.1770 4.96769 -0.4055640 -1.274350
3 ID_100 ID_100  2  2 -1.691620 60 -14.2821 3.01444 -0.0835603  1.622540
4 ID_101 ID_101  2  2 -3.979340 52 -12.2061 5.93377 -0.4792830  3.473790
5 ID_102 ID_102  2  2 -2.377090 69 -13.6426 4.53418 -2.2062600  0.222482
6 ID_103 ID_103  2  2  0.268657 67 -14.7840 5.09824 -3.5448200  5.679140
    Conf_7   Conf_8    Conf_9   Conf_10   Conf_11  Conf_12 Conf_13 Conf_13
1  0.865562 -2.26957 -0.0965859 -2.354970  1.058890  0.195302   0   7
2 -9.632180 -1.74190  2.6967400 -1.321190  0.742495 -0.688296   1  13
3 -8.066920 -1.65889 -2.7431200 -0.944793  7.769740  0.140430   1  10
4 -3.995170  1.17464  2.9239000 -2.148990 -7.193830  2.918810   0  20
5 -3.869480  0.72570  0.3461710 -0.882131 -2.526320 -1.620050   1  20
6 10.636900 -2.92550  0.0215072 -1.328180 -0.774790  1.984740   1  20

~~~

### 2.2. Target input files
**Qexp_tar.txt:** This is a `.txt` file which contains the following columns in order. The target dataset has 200 individuals who are independent from the discovery dataset. Please note that the file does not have column headings.

~~~
      FID   IID   Exposure
1 ID_801 ID_801 -0.64402046
2 ID_802 ID_802 -0.95209579
3 ID_803 ID_803 -0.02786981
4 ID_804 ID_804 -0.64402046
5 ID_805 ID_805  0.28020552
6 ID_806 ID_806  0.89635617
~~~

**Qphen_tar_cov.txt:** The exposure variable and other covariates for the adjustments are for the target dataset and should be provided in a separate file. 

~~~
   FID      IID   Outcome   Conf_1  Conf_2 Conf_3  Conf_4  Conf_5  Conf_6
1 ID_1000 ID_1000 27.0969 -3.247510 44 -11.6394 5.47484 -3.958770 0.9431600
2  ID_801  ID_801 26.5723 -3.826590 69 -13.8514 3.96080 -1.788050 0.0692473
3  ID_802  ID_802 20.2632  2.065150 60 -12.2438 4.04169 -0.905739 5.9656000
4  ID_803  ID_803 27.7365 -0.795863 62 -10.9195 6.91985 -2.920880 1.2601900
5  ID_804  ID_804 18.7500 -2.620880 67  -9.9271 4.10960 -2.354540 0.7190210
6  ID_805  ID_805 23.3025 -3.331640 67 -11.8637 5.88272  1.072880 2.7448800
    Conf_7   Conf_8    Conf_9    Conf_10   Conf_11  Conf_12 Conf_13 Conf_14
1  9.25721 -0.542392 -1.4608900  1.836640  1.80581  0.582524   0  20
2 -6.32556  2.853590  1.0851600 -1.303040  3.41659  1.415770   0   7
3  8.35545 -1.435760 -0.6181530  0.746918  5.11019 -0.207188   1  19
4 -5.56624 -0.552624 -0.0756095 -0.910047 -1.33896  1.726360   0   7
5 -1.82806 -1.821070  1.2157400 -3.566930 -7.91232  2.710110   0  10
6 -7.32776 -2.394770 -3.0798300 -1.436250  2.08822  1.429390   1  15
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
 # Access results

~~~

~~~
 print(result1$model_summary)
~~~

~~~
Coefficients:
                Estimate Std. Error t value Pr(>|t|)
(Intercept)     29.47457    3.66878   8.034 1.17e-13 ***
Add_PRS         -0.24741    0.36531  -0.677   0.4991
Int_PRS         -0.07306    0.40494  -0.180   0.8570
Covariate_Pheno -0.06637    0.36530  -0.182   0.8560
Conf_1          -0.31775    0.13155  -2.415   0.0167 *
Conf_2          -0.03132    0.04156  -0.754   0.4521
Conf_3           0.19922    0.23642   0.843   0.4005
Conf_4           0.38424    0.22680   1.694   0.0919 .
Conf_5          -0.38378    0.22059  -1.740   0.0836 .
Conf_6          -0.13712    0.15672  -0.875   0.3828
Conf_7           0.09332    0.06990   1.335   0.1835
Conf_8          -0.21371    0.22681  -0.942   0.3473
Conf_9          -0.29075    0.21862  -1.330   0.1852
Conf_10         -0.19994    0.18630  -1.073   0.2846
Conf_11         -0.01668    0.08661  -0.193   0.8475
Conf_12          0.06303    0.18779   0.336   0.7375
Conf_13          2.11203    0.67668   3.121   0.0021 **
Conf_14         -0.11214    0.06992  -1.604   0.1105
Int_PRS:Cov_PRS  0.17233    0.43474   0.396   0.6923
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
~~~

~~~
 print(result2$model_summary)
~~~

~~~
Coefficients:
                Estimate Std. Error t value Pr(>|t|)
(Intercept)     28.10497    3.70235   7.591 1.65e-12 ***
Add_PRS         -0.31083    0.36369  -0.855  0.39387
Int_PRS         -0.12959    0.40260  -0.322  0.74792
Covariate_Pheno -0.01331    0.36327  -0.037  0.97082
Cov_PRS         -0.81484    0.40687  -2.003  0.04671 *
Conf_1          -0.34898    0.13140  -2.656  0.00862 **
Conf_2          -0.02295    0.04143  -0.554  0.58033
Conf_3           0.19151    0.23451   0.817  0.41520
Conf_4           0.45236    0.22749   1.988  0.04827 *
Conf_5          -0.38868    0.21879  -1.776  0.07734 .
Conf_6          -0.10264    0.15638  -0.656  0.51245
Conf_7           0.09179    0.06933   1.324  0.18720
Conf_8          -0.17361    0.22584  -0.769  0.44306
Conf_9          -0.24493    0.21803  -1.123  0.26277
Conf_10         -0.24231    0.18597  -1.303  0.19427
Conf_11         -0.01347    0.08591  -0.157  0.87557
Conf_12          0.03782    0.18667   0.203  0.83967
Conf_13          2.17381    0.67183   3.236  0.00144 **
Conf_14         -0.08832    0.07036  -1.255  0.21101
Int_PRS:Cov_PRS  0.06421    0.43454   0.148  0.88269
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
# Access results
~~~

~~~
 print(result1$model_summary)
~~~

~~~
Coefficients:
                Estimate Std. Error z value Pr(>|z|)
(Intercept)      4.79087    3.55061   1.349  0.17724
Add_PRS         -1.06934    1.05734  -1.011  0.31185
Int_PRS         -0.79235    1.03506  -0.766  0.44397
Covariate_Pheno  1.19664    0.66575   1.797  0.07227 .
Conf_1           0.04700    0.36126   0.130  0.89649
Conf_2           0.06718    0.09686   0.694  0.48793
Conf_3          -0.04744    0.03915  -1.212  0.22569
Conf_4           0.25827    0.23471   1.100  0.27116
Conf_5           0.14980    0.23356   0.641  0.52127
Conf_6           0.01486    0.20004   0.074  0.94079
Conf_7          -0.05879    0.16260  -0.362  0.71765
Conf_8           0.04264    0.08369   0.510  0.61040
Conf_9           0.07903    0.21139   0.374  0.70853
Conf_10          0.25093    0.17080   1.469  0.14178
Conf_11          0.09226    0.17158   0.538  0.59076
Conf_12         -0.01252    0.08042  -0.156  0.87633
Conf_13          0.03893    0.15734   0.247  0.80458
Conf_14         -0.21132    0.07689  -2.748  0.00599 **
Int_PRS:Cov_PRS -0.11467    0.49826  -0.230  0.81799
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
~~~

~~~
 print(result2$model_summary)
~~~

~~~
Coefficients:
                Estimate Std. Error z value Pr(>|z|)
(Intercept)      4.65429    3.60788   1.290  0.19704
Add_PRS         -1.04372    1.06025  -0.984  0.32491
Int_PRS         -0.76288    1.04146  -0.733  0.46386
Covariate_Pheno  1.20269    0.66735   1.802  0.07152 .
Cov_PRS          0.08265    0.41812   0.198  0.84331
Conf_1           0.03928    0.36604   0.107  0.91454
Conf_2           0.07012    0.09820   0.714  0.47520
Conf_3          -0.04637    0.03954  -1.173  0.24087
Conf_4           0.25241    0.23705   1.065  0.28695
Conf_5           0.14800    0.23394   0.633  0.52696
Conf_6           0.01198    0.20148   0.059  0.95259
Conf_7          -0.05983    0.16295  -0.367  0.71350
Conf_8           0.04304    0.08384   0.513  0.60766
Conf_9           0.07503    0.21228   0.353  0.72375
Conf_10          0.25175    0.17060   1.476  0.14003
Conf_11          0.09751    0.17379   0.561  0.57475
Conf_12         -0.01241    0.08028  -0.155  0.87711
Conf_13          0.03857    0.15783   0.244  0.80694
Conf_14         -0.21029    0.07708  -2.728  0.00637 **
Int_PRS:Cov_PRS -0.10157    0.50393  -0.202  0.84027
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
~~~

 ### 2. Reverse causal direction
 
To evaluate the **reverse causal direction**, re-analyze the same dataset by switching the roles of the exposure and outcome variables. This means treating the previously defined outcome variable as the new exposure, and the previous exposure variables as the new outcome. Rearrange the data using the same structure and formatting approach used for the proposed causal directions as mentioned above, ensuring consistency across analysis pipeline. The only difference should be the reassignment of variable roles.

#### 2.1. Quantitative outcome

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
 # Access results
~~~

~~~
 print(result1$model_summary)
~~~

~~~
Coefficients:
                 Estimate Std. Error t value Pr(>|t|)
(Intercept)      0.195039   0.876070   0.223   0.8241
Add_PRS          0.613131   0.657990   0.932   0.3527
Int_PRS          0.494480   0.653759   0.756   0.4504
Covariate_Pheno -0.002215   0.015269  -0.145   0.8848
Conf_1           0.031787   0.027421   1.159   0.2479
Conf_2          -0.006169   0.008460  -0.729   0.4668
Conf_3          -0.032482   0.048167  -0.674   0.5009
Conf_4          -0.049913   0.047018  -1.062   0.2898
Conf_5           0.084615   0.044749   1.891   0.0602 .
Conf_6           0.022504   0.032115   0.701   0.4844
Conf_7           0.001957   0.014289   0.137   0.8912
Conf_8           0.017355   0.046416   0.374   0.7089
Conf_9           0.047086   0.044526   1.057   0.2917
Conf_10         -0.024777   0.038482  -0.644   0.5205
Conf_11          0.012215   0.017624   0.693   0.4891
Conf_12          0.013804   0.038327   0.360   0.7191
Conf_13         -0.036153   0.141357  -0.256   0.7984
Conf_14          0.002512   0.014583   0.172   0.8634
Int_PRS:Cov_PRS  0.067010   0.082983   0.808   0.4204
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
~~~

~~~
 print(result2$model_summary)
~~~

~~~
Coefficients:
                  Estimate Std. Error t value Pr(>|t|)
(Intercept)      0.0636244  0.8764312   0.073    0.942
Add_PRS          0.6812089  0.6567046   1.037    0.301
Int_PRS          0.5963822  0.6542637   0.912    0.363
Covariate_Pheno -0.0006459  0.0152386  -0.042    0.966
Cov_PRS          0.1198998  0.0759055   1.580    0.116
Conf_1           0.0320716  0.0273094   1.174    0.242
Conf_2          -0.0044344  0.0084962  -0.522    0.602
Conf_3          -0.0321723  0.0479693  -0.671    0.503
Conf_4          -0.0464272  0.0468771  -0.990    0.323
Conf_5           0.0857182  0.0445706   1.923    0.056 .
Conf_6           0.0225863  0.0319837   0.706    0.481
Conf_7           0.0020104  0.0142307   0.141    0.888
Conf_8           0.0190301  0.0462377   0.412    0.681
Conf_9           0.0387226  0.0446581   0.867    0.387
Conf_10         -0.0237400  0.0383301  -0.619    0.536
Conf_11          0.0118164  0.0175533   0.673    0.502
Conf_12          0.0077673  0.0383607   0.202    0.840
Conf_13         -0.0581224  0.1414619  -0.411    0.682
Conf_14          0.0024484  0.0145233   0.169    0.866
Int_PRS:Cov_PRS  0.0538590  0.0830611   0.648    0.518
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
# Access results
~~~

~~~
print(result1$model_summary)
~~~

~~~
Coefficients:
                Estimate Std. Error z value Pr(>|z|)
(Intercept)     -2.15188    1.77293  -1.214   0.2248
Add_PRS         -0.82112    2.15921  -0.380   0.7037
Int_PRS         -0.85888    2.16276  -0.397   0.6913
Covariate_Pheno  1.07339    0.60255   1.781   0.0748 .
Conf_1           0.44001    0.17217   2.556   0.0106 *
Conf_2           0.02299    0.05770   0.398   0.6903
Conf_3          -0.01009    0.01907  -0.529   0.5966
Conf_4          -0.14219    0.10749  -1.323   0.1859
Conf_5          -0.04528    0.12053  -0.376   0.7071
Conf_6           0.20845    0.10475   1.990   0.0466 *
Conf_7           0.13222    0.07653   1.728   0.0840 .
Conf_8          -0.07370    0.03654  -2.017   0.0437 *
Conf_9           0.10079    0.09950   1.013   0.3111
Conf_10         -0.06572    0.08634  -0.761   0.4465
Conf_11          0.07486    0.08792   0.851   0.3946
Conf_12          0.02629    0.03851   0.683   0.4948
Conf_13          0.04566    0.08021   0.569   0.5692
Conf_14          0.06266    0.03230   1.940   0.0524 .
Int_PRS:Cov_PRS -0.13343    0.21118  -0.632   0.5275
~~~

~~~
print(result2$model_summary)
~~~

~~~
Coefficients:
                Estimate Std. Error z value Pr(>|z|)
(Intercept)     -2.21475    1.77687  -1.246  0.21261
Add_PRS         -0.37233    2.19574  -0.170  0.86535
Int_PRS         -0.38568    2.20229  -0.175  0.86098
Covariate_Pheno  1.11820    0.60481   1.849  0.06448 .
Cov_PRS          0.23848    0.20303   1.175  0.24015
Conf_1           0.47938    0.17698   2.709  0.00675 **
Conf_2           0.02029    0.05797   0.350  0.72634
Conf_3          -0.01138    0.01921  -0.592  0.55357
Conf_4          -0.14879    0.10773  -1.381  0.16723
Conf_5          -0.06125    0.12197  -0.502  0.61553
Conf_6           0.19845    0.10547   1.882  0.05988 .
Conf_7           0.12962    0.07675   1.689  0.09124 .
Conf_8          -0.07542    0.03670  -2.055  0.03988 *
Conf_9           0.08542    0.10129   0.843  0.39905
Conf_10         -0.06209    0.08653  -0.717  0.47307
Conf_11          0.08333    0.08844   0.942  0.34606
Conf_12          0.02379    0.03884   0.612  0.54021
Conf_13          0.04283    0.08098   0.529  0.59689
Conf_14          0.06741    0.03271   2.061  0.03932 *
Int_PRS:Cov_PRS -0.13733    0.21345  -0.643  0.51997
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
~~~
 

