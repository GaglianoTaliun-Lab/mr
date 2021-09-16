# mr
Mendelian randomization protocol and plots



#### Power calculation example
via the "Online sample size and power calculator for Mendelian randomization with a binary outcome" https://sb452.shinyapps.io/power/ based on https://pubmed.ncbi.nlm.nih.gov/24608958/. 

To illustrate, let's use Williams et al. 2020 publication https://doi.org/10.1002/ana.25880 looking at LDL-C levels (exposure) and Parkinson's disease (outcome).
In Supp Table 4, they show that an MR analysis has **91.1% power** with the following criteria: **OR per SD difference=0.95**, using an instrument that explains 1.4% variance (**R2=0.014**) in LDL-C (e.g. variance explained by a single, sentinel SNP in the PCSK9 gene).    
The power calculation results are based off of **1,019,060 Parkinson's disease samples**: 37,688 cases + 981,372 controls (totals listed in Supp. Table 1), with a **1:26 case:control ratio** (or the tool also takes the proportion of cases is 0.03846153846). They set **α (significance level) as 0.05**.

We can plug in these values into the online calculator (for binary outcomes) and see the same result (91.1% power): 

![screenshot of Online sample size and power calculator for Mendelian randomization with a binary outcome.pdf](https://github.com/GaglianoTaliun-Lab/mr/files/7180905/Online.sample.size.and.power.calculator.for.Mendelian.randomization.with.a.binary.outcome.pdf)
