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
~~~
Or the CRAN version via
~~~
install.packages("GxEprs") 
~~~

Load the library

~~~
library(GxEprs)
~~~

**Data preparations**
The dataset is divided into **discovery (80%)** and **target (20%)** subsets, ensuring consistency across genetic, outcome, exposure, and confounder data. Genetic data, stored in **PLINK binary format** (`.bed`, `.bim`, `.fam`). The outcome variable should include **FID, IID, and phenotype values**, where case-control labels follow PLINK conventions: **1 = Control, 2 = Case** in the discovery dataset, and **0 = Control, 1 = Case** in the target dataset. Exposure and confounder variables are formatted into at least **19 columns** (**FID, IID, exposure, confounder1-n **) and partitioned in the same proportions. This structured approach ensures compatibility across all data types, thereby maintaining accurate alignment to estimate GxE interactions.


A Guide for GCIM analyses

GCIM analyses use PLink2 to analyze discovery data, and the package is compatible with the Linux operating system. 
1. Download plink2 software from the [Plink](https://www.cog-genomics.org/plink/2.0/) website and then specify the executable Plink file path.
   
~~~
plink_path <- "<plink_path>/plink2"
~~~
**Discovery dataset**: In the discovery dataset, please use the data input from the GxEprs data input format, except the value of e(exposure variable)-squared, which should be prepared as constant values to remove its effect from the model. Compute the PRS of the exposure using [Plink](https://www.cog-genomics.org/plink/2.0/) based on the GWAS from the discovery dataset for the PRS of the target samples. Use input and data preparation from GxEprs.
**Target dataset**: The target dataset for the model is also similar in data format, except, the use of PRS of the exposure variable(PRS of E) rather than the use of entire exposure values, and also includes the constant values for the fourth column, which is the square of the third column in the GxEprs model. 
GWAS inputs for the covariate from the discovery dataset as 

~~~
ID_1 ID_1 -0.644020461494502 0.414762354823591 -3.83142 64 -14.0364 5.51742 0.0714337 5.66263 0.865562 -2.26957 -0.0965859 -2.35497 1.05889 0.195302 0 7
ID_2 ID_2 -0.0278698075809014 0.00077672617459647 0.614044 66 -10.8505 2.11998 -0.882883 -0.441662 -2.64177 2.78944 0.524586 2.67134 -2.63724 -0.998764 1 20
ID_3 ID_3 2.1286574811167 4.53118267191409 -0.237792 55 -9.75369 3.18343 -2.09793 6.87345 11.3777 2.96961 -1.11879 0.873649 3.35523 -4.57831 1 10
ID_4 ID_4 2.1286574811167 4.53118267191409 6.69866 47 -9.07045 0.956878 -2.48407 1.06359 -3.13247 2.1232 -0.00976751 0.820582 0.0305345 1.6303 1 20
ID_5 ID_5 -0.952095788451302 0.906486390386706 -1.61423 59 -12.9379 1.29461 -1.79973 1.44404 -6.82898 -2.96795 -2.91577 -1.82881 7.15892 2.10916 1 20
ID_6 ID_6 -0.0278698075809014 0.00077672617459647 -4.38927 52 -11.8516 0.888978 -2.7231 1.11681 -3.64676 -0.594538 -1.7543 -0.716014 -2.39067 1.31295 1 10
ID_7 ID_7 -1.2601711154081 1.5880312401089 -0.0325543 64 -8.91617 5.21113 -3.03348 -6.11197 -0.632259 5.46586 2.28142 2.31582 1.1452 -4.93502 1 10
ID_8 ID_8 -1.2601711154081 1.5880312401089 2.66342 41 -11.1959 1.72515 -1.07217 -0.149743 -0.446028 1.84385 0.969557 0.233461 -2.76736 2.46577 1 10
~~~~

Compute PRS values and merge other covariates for the adjustment in the target sample. 
The prepared data should be presented without column names. 

~~~
  FID     IID   SCORE1_AVG CONST_COL        V3 V4       V5      V6
1 ID_1000 ID_1000 -0.000513025         2 -3.247510 44 -11.6394 5.47484
2  ID_801  ID_801  0.000173346         2 -3.826590 69 -13.8514 3.96080
3  ID_802  ID_802  0.000136668         2  2.065150 60 -12.2438 4.04169
4  ID_803  ID_803 -0.000743381         2 -0.795863 62 -10.9195 6.91985
5  ID_804  ID_804 -0.001379800         2 -2.620880 67  -9.9271 4.10960
6  ID_805  ID_805 -0.000547792         2 -3.331640 67 -11.8637 5.88272
         V7        V8       V9       V10        V11       V12      V13
1 -3.958770 0.9431600  9.25721 -0.542392 -1.4608900  1.836640  1.80581
2 -1.788050 0.0692473 -6.32556  2.853590  1.0851600 -1.303040  3.41659
3 -0.905739 5.9656000  8.35545 -1.435760 -0.6181530  0.746918  5.11019
4 -2.920880 1.2601900 -5.56624 -0.552624 -0.0756095 -0.910047 -1.33896
5 -2.354540 0.7190210 -1.82806 -1.821070  1.2157400 -3.566930 -7.91232
6  1.072880 2.7448800 -7.32776 -2.394770 -3.0798300 -1.436250  2.08822
        V14 V15
1  0.582524   0
2  1.415770   0
3 -0.207188   1
4  1.726360   0
5  2.710110   0
6  1.429390   1
~~~

The same approach will be applied for the reverse direction test.
The regression outputs look like those for a quantitative outcome with a quantitative exposure variable. 

~~~
                 Coefficient   Std.Error Test.Statistic    pvalue
(Intercept)     -0.301958164 0.767536441    -0.39341215 0.6944756
E                0.023544513 0.082223120     0.28634905 0.7749365
PRS_add          0.007905589 0.076567153     0.10325040 0.9178779
PRS_gxe         -0.097510922 0.076129896    -1.28084928 0.2018764
`PRS_gxe x E`   -0.038938148 0.078699730    -0.49476850 0.6213603
`Confounder 1`  -0.002021759 0.029510412    -0.06851002 0.9454549
`Confounder 2`   0.003624143 0.009210592     0.39347557 0.6944288
`Confounder 3`  -0.014329725 0.051864775    -0.27629013 0.7826387
`Confounder 4`  -0.008979054 0.050744150    -0.17694757 0.8597463
`Confounder 5`   0.018699159 0.048476503     0.38573655 0.7001423
`Confounder 6`   0.004797342 0.035119207     0.13660165 0.8914967
`Confounder 7`   0.005113243 0.015571517     0.32837156 0.7430082
`Confounder 8`  -0.044123338 0.051005389    -0.86507209 0.3881380
`Confounder 9`  -0.019639324 0.047748447    -0.41130812 0.6813305
`Confounder 10`  0.020840701 0.041509631     0.50206905 0.6162258
`Confounder 11` -0.030173807 0.019127111    -1.57754123 0.1164076
`Confounder 12` -0.005893285 0.041673315    -0.14141628 0.8876975
`Confounder 13` -0.026170287 0.150009885    -0.17445708 0.8617002
~~~

Example data:
Quantitative outcome with Quantitative exposure:
discovery dataset: 
1. Conduct GWEIS for the outcome
2. Conduct GWAS
3. Estimate the PRS of the exposure
4. Conduct linear regression

# Quantitative outcome 
~~~
# Compute PRS of exposure variables
Conduct a GWAS for the exposure adjusted by covariates, and then compute the PRS for the exposure. 
~~~

~~~
# Conduct GWEIS based on the assigned outcome variables, exposure with confounders.
b <- GWEIS_quantitative(plink_path, "DummyData", "Qphe_discovery.txt", "Qcov_discovery.txt")
add <- b[c("ID", "A1", "ADD_BETA")]
gxe <- b[c("ID", "A1", "INTERACTION_BETA")]
~~~

~~~
# Compute the PRS of each summary 
q <- PRS_quantitative(plink_path, "DummyData", summary_input = add)
r <- PRS_quantitative(plink_path, "DummyData", summary_input = gxe)
~~~

~~~
# Conduct regression analyses for GCIM
#PRS of the exposure(PRS_exp) variables should be the third column, and include the constant values for the fourth column, and then add the  number of covariates for the adjustment. 
y <- summary_regular_quantitative("Qpt.txt", "PRS_exp_cov.txt", add_score = q, gxe_score = r, Model = 4)
~~~

~~~
y$summary
~~~
Perform the reverse analyses by changing the roles of exposure and outcome variables. 

# Binary outcome

**Compute PRS of exposure variables**
Conduct a GWAS for the exposure adjusted by covariates, and then compute the PRS for the exposure.
Conduct a log transformation of the OR if the exposure variable is binary.  Then, compute the PRS of the exposure variables using PRS score values for the target sample. 

~~~
b <- GWEIS_binary(plink_path, "DummyData", "Bpd.txt", "Bcd.txt")
add <- b[c("ID", "A1", "ADD_BETA")]
gxe <- b[c("ID", "A1", "INTERACTION_BETA")]
~~~

~~~
q <- PRS_binary(plink_path, "DummyData", summary_input = add)
r <- PRS_binary(plink_path, "DummyData", summary_input = gxe)
~~~

~~~
z <- summary_regular_binary("Bpt.txt", "PRS_exp_cov.txt", add_score = q, gxe_score = r, Model = 5)
~~~

~~~
z$summary
~~~
