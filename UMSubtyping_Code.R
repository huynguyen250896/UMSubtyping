1. PRE-PROCESSING


#### Library
library(CancerSubtypes)
library(tidyr)
library(tidyverse)
library(dplyr)

#### Import data
met_UM<-read.table('data_methylation_hm450.txt',sep='\t',row.names = 1, check.names = FALSE, header = TRUE)
exp_UM<-read.table('data_RNA_Seq_v2_mRNA_median_Zscores.txt',sep="\t",row.names = 1,  check.names = FALSE, header=TRUE)
cna_UM<-read.table('data_CNA.txt', sep = '\t', row.names = 1, check.names = FALSE, header = TRUE)

####Remove Entrez_Gene_Id
exp_UM = exp_UM[,-1]
cna_UM = cna_UM[,-1]
met_UM = met_UM[,-1]

#### Pre-process Clinical data
cli = read.table('UM_clinical_patient.txt', sep='\t', row.names = 1,  check.names = FALSE, header = TRUE, fill=TRUE)
rownames(cli) = paste(rownames(cli), "", sep = "-01") #rename rows
cli_patient_UM = na_if(cli, "[Not Available]") #transform [Not Available] into NA value
cli_patient_UM$event = ifelse(cli_patient_UM$`Vital Status` == 'Dead',1,0) #create binary events
cli_patient_UM$overall = cli_patient_UM$`3. Time to UM Death or Last Follow-up`

#### Retain patients sharing among three datasets
scna_UM<-colnames(cna_UM)
sexp_UM<-colnames(exp_UM)
smet_UM<-colnames(met_UM)
s_cna_exp_UM<-intersect(scna_UM,sexp_UM)
s_cna_exp_met_UM<-intersect(smet_UM,s_cna_exp_UM)
s_cna_exp_met_cli_UM<-intersect(rownames(cli_patient_UM),s_cna_exp_met_UM)

# Only retain shared samples among 3 datasets
c_cna_UM = cna_UM[,s_cna_exp_met_cli_UM]
c_met_UM = met_UM[,s_cna_exp_met_cli_UM]
c_exp_UM = exp_UM[,s_cna_exp_met_cli_UM]
cli_patient_UM = cli_patient_UM[s_cna_exp_met_cli_UM,]

#### Transform list to matrix
c_exp_UM = as.matrix(c_exp_UM)
c_cna_UM = as.matrix(c_cna_UM)
c_met_UM = as.matrix(c_met_UM)

####check missing value existing?
table(is.finite(c_exp_UM))
table(is.finite(c_cna_UM)) 
table(is.finite(c_met_UM))

#### Impute missing value
# Remove probes >50% missing
c_exp_UM <- c_exp_UM[rowSums(is.na(c_exp_UM)) < (ncol(c_exp_UM) * .5), ]
c_met_UM <- c_met_UM[rowSums(is.na(c_met_UM)) < (ncol(c_met_UM) * .5), ]
# Check missing value again
table(is.finite(c_exp_UM))
table(is.finite(c_met_UM))


#### Impute using knn method provided in the function data.imputation
c_exp_UM = data.imputation(c_exp_UM, fun = "microarray")
c_met_UM = data.imputation(c_met_UM, fun = "microarray")

#### Shared probes are reserved between EXP + CNA and EXP + MET
GENE_cna_exp_UM<-intersect(rownames(c_cna_UM),rownames(c_exp_UM))
length(GENE_cna_exp_UM) 
c_cna_UM = c_cna_UM[GENE_cna_exp_UM,]
c_exp1_UM = c_exp_UM[GENE_cna_exp_UM,]
GENE_met_exp_UM<-intersect(rownames(c_met_UM),rownames(c_exp_UM))
length(GENE_met_exp_UM) 
c_met_UM = c_met_UM[GENE_met_exp_UM,]
c_exp2_UM = c_exp_UM[GENE_met_exp_UM,]

2. SKEWNESS


#### Library
if(!require(devtools)) install.packages("devtools")
devtools::install_github("huynguyen250896/geneCor", force=TRUE)
library(geneCor)

#### Correlation between corresponding EXP + CNA
x=c_cna_UM; x=t(x)
y=c_exp1_UM; y=t(y)

#### Correlation between corresponding EXP + MET
x1=c_met_UM; x1=t(x1)
y1=c_exp2_UM; y1=t(y1)

#### Tool for Identification of CNAexp and METexp, Visualization of the distribution of expression of CNAexp genes and expression of METexp genes on a page, and Examination of the significance of each of those skewed distributions.
geneCor(cna = x, exp1 = y, alternative1="greater", met = x1, exp2 = y1, alternative2="greater", method = "spearman") #compute Spearman's Rank correlation coefficients.

#### Import cna_cor.txt and met_cor.txt 
CNAexp_UM = read.table("cna_cor.txt", sep="\t", row.names = 1, check.names = FALSE, header = TRUE)
METexp_UM = read.table("met_cor.txt", sep="\t", row.names = 1, check.names = FALSE, header = TRUE)

#### CNA/MET genes significantly correlated with gene expressison
length(rownames(CNAexp_UM))
length(rownames(METexp_UM))

#### CNAexp and METexp genes
c_cna_UM = c_cna_UM[rownames(CNAexp_UM),]
c_met_UM = c_met_UM[rownames(METexp_UM),]

#### Retain genes significantly associated with prognostic value (OS)
c_cna_UM=FSbyCox(c_cna_UM, cli_patient_UM$overall, cli_patient_UM$event, cutoff = 0.0005)
c_met_UM=FSbyCox(c_met_UM, cli_patient_UM$overall, cli_patient_UM$event, cutoff = 0.0005)
# Check dimension of the two variables 'c_cna_UM' and 'c_met_UM'
dim(c_cna_UM) 
dim(c_met_UM) 


3. CLUSTERING


cna_cor_UM = t(c_cna_UM)
met_cor_UM = t(c_met_UM)

#### Library
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("preprocessCore")
install.packages("flexclust")
install.packages(c("entropy", "pbmcapply", "doParallel", "foreach","matrixStats","Rcpp","RcppParallel",
                   "FNN", "cluster", "irlba", "mclust")) #install dependencies 
install.packages("PINSPlus")

library(PINSPlus)
library(parallel)
library(flexclust)

#### single clustering using CNAexp and METexp genes
ncore = 2

set.seed(1)
#### CNAexp
result = PerturbationClustering(data=cna_cor_UM, kMax = 10, verbose = F)
result$Discrepancy$AUC 
sub_grp = result$cluster
table(sub_grp)

 
result.5 = PerturbationClustering(data=cna_cor_UM, k = 5, verbose = F)
sub_grp.5 = result.5$cluster
table(sub_grp.5)


#### METexp
result1 = PerturbationClustering(data=met_cor_UM, kMax = 10, verbose = F)
result1$Discrepancy$AUC 
sub_grp1 = result1$cluster
table(sub_grp1)

#### Plot the result for the single clustering using CNAexp genes
plot(1:10,result$Discrepancy$AUC,type="b",xlab="Number of Subgroups", ylab="AUC")
# add axes
axis(1, at =seq(1,10,1));axis(1, at =seq(1,10,2))
#abline
abline(v=c(4,5),lty ='longdash')
title("Subgroup number\n identification of CNAexp")


#### Plot the result for the single clustering using METexp genes
plot(1:10,result1$Discrepancy$AUC,type="b",xlab="Number of Subgroups", ylab="AUC")
# add axes
axis(1, at =seq(1,10,1));axis(1, at =seq(1,10,2))
#abline
abline(v=c(2),lty ='longdash')
title("Subgroup number\n identification of METexp")

#### Integrative clustering using CNAexp and METexp genes
exp_cor_UM=rbind(c_exp1_UM, c_exp2_UM); exp_cor_UM=t(exp_cor_UM) #combining two corresponding mRNA datasets to one
dataList <- list (cna_cor_UM, met_cor_UM, exp_cor_UM)
names(dataList) = c("CNA", "MET", "EXP")

set.seed(2)
result2 <- SubtypingOmicsData(dataList = dataList, ncore = ncore, kMax = 10)
result2[["dataTypeResult"]][["CNA"]][["Discrepancy"]][["AUC"]]
result2[["dataTypeResult"]][["MET"]][["Discrepancy"]][["AUC"]]
result2[["dataTypeResult"]][["EXP"]][["Discrepancy"]][["AUC"]]
sub_grp2 = result2$cluster1 
table(sub_grp2)
 
#### Plot the result for the result2 variable
plot(1:10,result2[["dataTypeResult"]][["EXP"]][["Discrepancy"]][["AUC"]],type="b",xlab="Number of Subgroups", ylab="AUC")
# add axes
axis(1, at =seq(1,10,1));axis(1, at =seq(1,10,2))
#abline
abline(v=c(2),lty ='longdash')
title("Integrative subgroup number\n identification of mRNA + CNAexp + METexp")

#### Overlaps
set.seed(3)
# Single
chi=chisq.test(sub_grp,sub_grp1)
p.overlap=chi$p.value;p.overlap

chi.5=chisq.test(sub_grp.5,sub_grp1)
p.overlap.5=chi.5$p.value;p.overlap.5

chi.cna = chisq.test(sub_grp,sub_grp.5)
p.overlap.cna=chi.cna$p.value;p.overlap.cna

#### Plot the result
plot(chi$observed, xlab="METexp Subgroups ", ylab="CNAexp Subgroups", main="")
mtext(paste('chisq-P-value =', p.overlap), cex=0.85)
title("Overlap test between subgroups of\n METexp clustering and 4 subgroups of CNAexp clustering \n", cex.main = 0.9)

plot(chi.5$observed, xlab="METexp Subgroups ", ylab="CNAexp Subgroups", main="")
mtext(paste('chisq-P-value =', p.overlap.5), cex=0.85)
title("Overlap test between subgroups of\n METexp clustering and 5 subgroups of CNAexp clustering \n", cex.main = 0.9)

plot(chi.cna$observed, xlab="4 CNAexp Subgroups ", ylab="5 CNAexp Subgroups", main="")
mtext(paste('chisq-P-value =', p.overlap.cna), cex=0.85)
title("Overlap test between 5 subgroups of\n METexp clustering and 4 subgroups of CNAexp clustering \n", cex.main = 0.9)


#### Integrative
# CNAexp and IntSubtyping
chi1=chisq.test(sub_grp,sub_grp2)
p.overlap1=chi1$p.value;p.overlap1


chi1.5=chisq.test(sub_grp.5,sub_grp2)
p.overlap1.5=chi1.5$p.value;p.overlap1.5


#METexp and IntSubtyping
chi2=chisq.test(sub_grp1,sub_grp2)
p.overlap2=chi2$p.value;p.overlap2 


plot(chi1$observed, xlab="CNAexp Subgroups", ylab="Integrative Subgroups", main="")
mtext(paste('chisq-P-value =', p.overlap1), cex=0.85)
title("Overlap test between 4 subgroups of CNAexp clustering\n and integrative subgroups\n", cex.main = 0.9)

plot(chi1.5$observed, xlab="CNAexp Subgroups", ylab="Integrative Subgroups", main="")
mtext(paste('chisq-P-value =', p.overlap1.5), cex=0.85)
title("Overlap test between 5 subgroups of CNAexp clustering\n and integrative subgroups\n", cex.main = 0.9)

plot(chi2$observed, xlab="METexp Subgroups", ylab="Integrative Subgroups", main="")
mtext(paste('chisq-P-value =', p.overlap2), cex=0.85)
title("Overlap test between subgroups of METexp clustering\n and integrative subgroups\n", cex.main = 0.9)

4. SURVIVAL RATE


#### Library
if(!require(survival)) install.packages('survival')
library(survival)
library(survminer) 

#### Single clustering
#CNA
set.seed(4)

coxFit <- coxph(
  Surv(overall, event) ~ as.factor(sub_grp),
  data = cli_patient_UM,
  ties = "exact"
)
pcox=summary(coxFit)$logtest[3];pcox


coxFit.5 <- coxph(
  Surv(overall, event) ~ as.factor(sub_grp.5),
  data = cli_patient_UM,
  ties = "exact"
)
pcox.5=summary(coxFit.5)$logtest[3];pcox.5


mfit <- survfit(Surv(overall, event) ~ as.factor(sub_grp), data = cli_patient_UM)
mfit.5 <- survfit(Surv(overall, event) ~ as.factor(sub_grp.5), data = cli_patient_UM)

#### Plot the result

ggsurvplot(mfit, size=1, title = "4 CNAexp Clustering",
           linetype = "strata",
           risk.table = TRUE, fun ="pct", risk.table.col = "strata", break.x.by = 50,
           xlab = "Time in days",
           legend = "bottom",
           legend.title = "Group", legend.labs = c("Group 1","Group 2", "Group 3", "Group 4"),
           conf.int = FALSE, pval = paste("P-value",pcox.5,sep=" = " ),pval.size = 3.3,
           xlim = c(0,200),
           font.main = c(12, "bold"),
           font.x = c(11, "bold.italic"),
           font.y = c(11, "bold.italic"),
           palette = c("turquoise", "red", "green", "violet")) 


ggsurvplot(mfit.5, size=1, title = "5 CNAexp Clustering",
           linetype = "strata",
           risk.table = TRUE, fun ="pct", risk.table.col = "strata", break.x.by = 50,
           xlab = "Time in days",
           legend = "bottom",
           legend.title = "Group", legend.labs = c("Group 1","Group 2", "Group 3", "Group 4","Group 5"),
           conf.int = FALSE, pval = paste("P-value",pcox.5,sep=" = " ),pval.size = 3.3,
           xlim = c(0,200),
           font.main = c(12, "bold"),
           font.x = c(11, "bold.italic"),
           font.y = c(11, "bold.italic"),
           palette = c("turquoise", "red", "green", "violet", "blue")) 

#MET
coxFit1 <- coxph(
  Surv(overall, event) ~ as.factor(sub_grp1),
  data = cli_patient_UM,
  ties = "exact"
)
pcox1=summary(coxFit1)$logtest[3];pcox1 

mfit1 <- survfit(Surv(overall, event == 1) ~ as.factor(sub_grp1), data = cli_patient_UM)

#### Plot the result
ggsurvplot(mfit1, size=1, title = "METexp Clustering",
           linetype = "strata",
           risk.table = TRUE, fun ="pct", risk.table.col = "strata", break.x.by = 50,
           xlab = "Time in days",
           legend = "bottom",
           legend.title = "Group", legend.labs = c("Group 1","Group 2"),
           conf.int = FALSE, pval = paste("P-value",pcox1,sep=" = " ),pval.size = 3.5,
           xlim = c(0,200),
           font.main = c(12, "bold"),
           font.x = c(11, "bold.italic"),
           font.y = c(11, "bold.italic"),
           palette = c("turquoise", "red")) 


#Integrative
coxFit2 <- coxph(
  Surv(overall, event) ~ as.factor(sub_grp2),
  data = cli_patient_UM,
  ties = "exact"
)
pcox2=summary(coxFit2)$logtest[3];pcox2 

mfit2 <- survfit(Surv(overall, event == 1) ~ as.factor(sub_grp2), data = cli_patient_UM)

ggsurvplot(mfit2, size=1, title = "Integrated clustering",
           linetype = "strata",
           risk.table = TRUE, fun ="pct", risk.table.col = "strata", break.x.by = 50,
           xlab = "Time in days",
           legend = "bottom",
           legend.title = "Group", legend.labs = c("IntSub1","IntSub2"),
           conf.int = FALSE, pval = paste("P-value",pcox2,sep=" = " ),xlim = c(0,200),
           font.main = c(12, "bold"),
           font.x = c(11, "bold.italic"),
           font.y = c(11, "bold.italic"),
           palette = c("turquoise", "red"))


5. IDENTIFY SUBTYPE-SCECIFIC GENES 

library(devtools)
install_github("huynguyen250896/GeneCluster")
library(GeneCluster)
install.packages("maps")
library(maps)

#### CNAexp Subgroups Specific Genes
SubtypeSpecificGene(omics = cna_cor_UM, cluster = sub_grp2)

#### METexp Subgroups Specific Genes
SubtypeSpecificGene(omics = met_cor_UM, cluster = sub_grp2)

#### Integrative Subgroups Specific Genes
SubtypeSpecificGene(omics = exp_cor_UM , cluster = sub_grp2)

6. PROGNOSTIC FACTOR IDENTIFICATION

#### Age factor

cli_patient_UM$integrate_subtype = sub_grp2
cli_patient_UM_old1 = cli_patient_UM[cli_patient_UM$integrate_subtype == 1,]
cli_patient_UM_old2 = cli_patient_UM[cli_patient_UM$integrate_subtype == 2,]

#### Set 2 group of old and non-old patients in each Integrative Subgroups
cli_patient_UM_old1$Old = ifelse(cli_patient_UM_old1$Age >65, 1, 0)
cli_patient_UM_old2$Old = ifelse(cli_patient_UM_old2$Age >65, 1, 0)


# Average age of each Integrative Subgroups
mean(cli_patient_UM_old1$Age)
mean(cli_patient_UM_old2$Age)


# Average survival time of each Integrative Subgroups
mean(cli_patient_UM_old1$overall)
mean(cli_patient_UM_old2$overall)

#### Survival rate
age_factor1 = cli_patient_UM_old1$Old
age_factor2 = cli_patient_UM_old2$Old

coxFit.age1 <- coxph(
  Surv(overall, event) ~ as.factor(age_factor1),
  data = cli_patient_UM_old1,
  ties = "exact"
)
pcox.age1=summary(coxFit.age1)$logtest[3];pcox.age1


coxFit.age2 <- coxph(
  Surv(overall, event) ~ as.factor(age_factor2),
  data = cli_patient_UM_old2,
  ties = "exact"
)
pcox.age2=summary(coxFit.age2)$logtest[3];pcox.age2
