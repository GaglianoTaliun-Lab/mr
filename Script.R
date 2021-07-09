#remove everything from workspace
rm(list=ls(all=TRUE))
#3 tell R to use the packages you installed
library(MendelianRandomization)
library(metafor)
par(mfrow=c(2,2)) 
par(mar=c(4,2,0,2)+.1)

#1b, read subsetted drinks per week file and find the column name for the SNP name
#It's on b37
#subsetted instrument SNPs from drinks per week sumstats
softc<-read.table("../input/outcome_snps.csv",sep=",",header=TRUE)
softc[1,]
nrow(softc)

#2 and #4.i read snp to exposure file 
#remember to use the location and file name matching the file you created
#neuroticism
#It's on b37
vtas<-read.csv("../input/my_exp_neuro_snps.csv",header=TRUE)
#check it is the correct file with the correct column names
vtas[1,]
nrow(vtas)

#Both exposure and outcome need BETAS --> fine here
#but if OR, use the following to convert, e.g.:
#beta=log(OR)
#sebeta=SE/exp(OR)


#4.ii select the genetic instrument SNPs from the outcome file 
sl<-softc[softc$RSID %in% vtas$SNP, ] # get SNP to outcome info
nrow(sl)
sl
#4.iii rename the columns in the outcome file to be the same for every file
sl$chrpos<-paste(sl$CHROM,sl$POS, sep=":")
snpto<-data.frame(SNP = sl$RSID,
                  beta.outcome = sl$BETA,
                  se.outcome = sl$SE,
                  effect_allele.outcome = sl$ALT,
                  other_allele.outcome =sl$REF,
                  eaf.outcome = sl$AF,
                  pval.outcome = sl$PVALUE,
                  chrpos.outcome = sl$chrpos             
)

#4.iv	Merge SNP to exposure and SNP to outcome files
pl<-merge(snpto,vtas,by="SNP",all.x=TRUE) #merge SNP to exposure and SNP to outcome file
nrow(pl)
pl[3,]

#4 v Align the SNPs on the same effect allele for exposure and outcome
# by changing sign of beta.exposure if effect alleles do not match
pl$effect_allele.outcome<-as.factor(pl$effect_allele.outcome)
pl$effect_allele.exposure<-as.factor(pl$effect_allele.exposure)
lev2 <- unique( c( levels(pl$effect_allele.outcome), levels(pl$effect_allele.exposure) ) )
pl$effect_allele.outcome <- factor(pl$effect_allele.outcome, levels=lev2)
pl$effect_allele.exposure <- factor(pl$effect_allele.exposure, levels=lev2)
pl$effect_allele.exposure<-gsub(" ", "",pl$effect_allele.exposure, fixed = TRUE)
pl$beta.exposure[pl$effect_allele.exposure!=pl$effect_allele.outcome]<-pl$beta.exposure[pl$effect_allele.exposure!=pl$effect_allele.outcome] * -1
pl[3,]

#4 vi.get forest plot with fixed effects, same as Mendelianrandomization IVW with fixed effects 
x<-pl$beta.exposure # beta for SNP to exposure
sigmax<-pl$se.exposure # its standard errors
y<-pl$beta.outcome # beta for SNP to outcome
sigmay<-pl$se.outcome # its standard errors
pl$Wald<-y/x #Wald estimate this is tau/alpha
pl$Waldvar<-(sigmay^2/x^2) # using Burgess's method
pl$lab<-paste(pl$SNP, pl$gene, sep=" ")
dmres<-rma.uni(yi=pl$Wald, vi=pl$Waldvar, slab=pl$lab, method="FE")
dmres
pdf("../output/Neuro_on_DrinksperW.pdf")
#forest(dmres, atransf=exp,xlab=" ", mlab="Drinks per week (OR)",at=log(c(0.07,0.5,1,3.2)),xlim=c(-0.1,1.7),cex=.7)
forest(dmres, atransf=exp,xlab=" ", mlab="Drinks per week (OR)", at=log(c(.5, 1,2)),xlim=c(-1.7,1.3),cex=.8)
dev.off()
#The black points represent the causal estimate (log-odds ratio of the Wald statistic) of each SNP on the risk of developing PD. Horizontal lines denote 95% confidence intervals. Point on the bottom of the figure represent sthe causal estimate when combining all SNPs together, using the inverse variance weighted.

#if there is an outlier... unhash the following to remove the outlying SNP(s) and replot the Forest plot
#plnew<-subset(pl, pl$SNP!="rs9384679" & pl$SNP!='rs2726491')
#pl<-plnew 
#dmres<-rma.uni(yi=pl$Wald, vi=pl$Waldvar, slab=pl$lab, method="FE")
#pdf("Forestplot-nooutlier.pdf")
#forest(dmres, atransf=exp,xlab=" ", mlab="Outcome Trait (OR)", at=log(c(.5, 1,2)),xlim=c(-1.7,1.3),cex=.8)
#dev.off()

#4 vii get estimates for outcome using Mendelian randomization package withfixed effects
MRInputObject <- mr_input(pl$beta.exposure, pl$se.exposure, pl$beta.outcome, pl$se.outcome)
mr_ivw(MRInputObject,model="fixed")
#nb exponentiate to get OR and 95% CIs matching the forest plot

#Inverse-variance weighted method
#(variants uncorrelated, fixed-effect model)

#Number of Variants : 3 

#------------------------------------------------------------------
# Method Estimate Std Error 95% CI       p-value
#    IVW    0.129     0.070 -0.009, 0.267   0.067
#------------------------------------------------------------------
#Residual standard error =  0.911 
#Residual standard error is set to 1 in calculation of confidence interval by fixed-effect assumption.
#Residual standard error is set to 1 in calculation of confidence interval when its estimate is less than 1.
#Heterogeneity test statistic = 1.6613 on 3 degrees of freedom, (p-value = 0.4358)

#4 viii get MR estimates for outcome using random effects
mr_ivw(MRInputObject)

#Inverse-variance weighted method
#(variants uncorrelated, random-effect model)

#Number of Variants : 4 

#------------------------------------------------------------------
# Method Estimate Std Error  95% CI       p-value
#    IVW    0.076     0.064 -0.048, 0.201   0.229
#------------------------------------------------------------------
#Residual standard error =  1.247 
#Residual standard error is set to 1 in calculation of confidence interval by fixed-effect assumption.
#Heterogeneity test statistic = 4.6647 on 3 degrees of freedom, (p-value = 0.1981)

#lab 4 2b get weighted median and MR-Egger estimates
mr_median(MRInputObject)

# Weighted median method 

#Number of Variants : 3 
#------------------------------------------------------------------
#                 Method Estimate Std Error  95% CI       p-value
# Weighted median method    0.136     0.084 -0.028, 0.300   0.105
#------------------------------------------------------------------

#------------------------------------------------------------------
#                 Method Estimate Std Error  95% CI       p-value
# Weighted median method    0.093     0.082 -0.067, 0.254   0.255
#------------------------------------------------------------------

mr_egger(MRInputObject)

#MR-Egger method
#(variants uncorrelated, random-effect model)

#Number of Variants =  3 

#------------------------------------------------------------------
#      Method Estimate Std Error  95% CI       p-value
#    MR-Egger    0.133     0.637 -1.115, 1.382   0.834
# (intercept)    0.000     0.011 -0.021, 0.021   0.995
#------------------------------------------------------------------
#Residual Standard Error :  1.289
#Residual standard error is set to 1 in calculation of confidence interval when its estimate is less than 1.
#Heterogeneity test statistic = 1.6612 on 1 degrees of freedom, (p-value = 0.1974)
#I^2_GX statistic: 56.8%

#Number of Variants =  4 

#------------------------------------------------------------------
#      Method Estimate Std Error  95% CI       p-value
#    MR-Egger    0.529     0.471 -0.395, 1.452   0.262
# (intercept)   -0.007     0.008 -0.022, 0.007   0.330
#------------------------------------------------------------------
#Residual Standard Error :  1.258 
#Heterogeneity test statistic = 3.1637 on 2 degrees of freedom, (p-value = 0.2056)
#I^2_GX statistic: 38.3%

#install MRPRESSO
if (!require("devtools")) { install.packages("devtools") } else {}
devtools::install_github("rondolab/MR-PRESSO")
library(MRPRESSO)

#run MR-Presso 
mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure",
          OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = pl, NbDistribution = 100000,  SignifThreshold = 0.05)
#Error in mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",  : 
#  Not enough intrumental variables

#$`Main MR results`
#       Exposure       MR Analysis Causal Estimate         Sd    T-stat   P-value
#1 beta.exposure               Raw      0.07643633 0.07923834 0.9646381 0.4058842
#2 beta.exposure Outlier-corrected              NA         NA        NA        NA

#$`MR-PRESSO results`
#$`MR-PRESSO results`$`Global Test`
#$`MR-PRESSO results`$`Global Test`$RSSobs
#[1] 7.631871

#$`MR-PRESSO results`$`Global Test`$Pvalue
#[1] 0.24512

#5 save the data for future use
write.csv(pl,"../output/Neuro_on_DrinksperW.csv")

#plot beta.exposure vs. beta.outcome for instruments
pdf('../output/exp_out-Neuro_on_DrinksperW.pdf')
plot(pl$beta.exposure, pl$beta.outcome, xlab= "Beta Exposure", ylab= "Beta Outcome", main="")
text(pl$beta.outcome ~pl$beta.exposure, labels=pl$SNP,data=pl, cex=0.5, font=2, pos=3)
abline(a=0, b=0.129, col='blue') #IVW estimate
abline(a=0, b=0.133, col='red') #MR=Egger estimate
abline(a=0, b=0.136, col="purple") #Weighted Median
legend("topleft", c("Inverse Variance Weighted", "Weighted Median", "MR-Egger"), text.col = c("blue", "purple", "red"), bty='n', cex=0.75)
dev.off()

pdf('../output/new-exp_out-Neuro_on_DrinksperW.pdf')
plot(pl$beta.exposure, pl$beta.outcome, xlab= "Beta Exposure", ylab= "Beta Outcome", main="", xlim=c(0,-0.12), ylim=c(-0.008,0.004))
text(pl$beta.outcome ~pl$beta.exposure, labels=pl$SNP,data=pl, cex=0.5, font=2, pos=3)
abline(a=0, b=0.076, col='blue') #IVW estimate
abline(a=-0.007, b=0.529, col='red') #MR=Egger estimate
abline(a=0, b=0.093, col="purple") #Weighted Median
legend("topright", c("Inverse Variance Weighted", "Weighted Median", "MR-Egger"), text.col = c("blue", "purple", "red"), bty='n', cex=0.75)
dev.off()

#plot Wald vs. Waldvar for instruments
#Plot of an ideal MR analysis will produce a symmetric scatterplot where data points with higher instrumental strength (approximated by inverse of standard error) will coalesce around the beta of the MR result.
pdf('../output/Funnel-Neuro_on_DrinksperW.pdf')
plot(pl$Wald, 1/pl$Waldvar, xlab= "Wald", ylab= "1/Waldvar", main="")
text(1/pl$Waldvar ~pl$Wald, labels=pl$SNP,data=pl, cex=0.5, font=2, pos=3)
abline(v=0.129, col='blue') #IWV estimate
abline(v=0.133, col="red") #MR-Egger estimate
abline(v=0.136, col="purple") #Weighted Median
legend("center", c("Inverse Variance Weighted", "Weighted Median", "MR-Egger"), text.col = c("blue", "purple", "red"), bty='n', cex=0.75)
dev.off()

pdf('../output/new-Funnel-Neuro_on_DrinksperW.pdf')
plot(pl$Wald, 1/pl$Waldvar, xlab= "Wald", ylab= "1/Waldvar", main="", xlim=c(-0.2,0.6))
text(1/pl$Waldvar ~pl$Wald, labels=pl$SNP,data=pl, cex=0.5, font=2, pos=3)
abline(v=0.076, col='blue') #IWV estimate
abline(v=0.529, col="red") #MR-Egger estimate
abline(v=0.093, col="purple") #Weighted Median
legend("topleft", c("Inverse Variance Weighted", "Weighted Median", "MR-Egger"), text.col = c("blue", "purple", "red"), bty='n', cex=0.75)
dev.off()
