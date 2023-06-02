##############################################################################################################################################################
##############################################################################################################################################################
####################    Optimal Two-Stage Sampling for Mean Estimation in Multilevel Populations when Cluster Size is Informative.   #########################
####################                                           Supplementary Material 2                                              #########################
##############################################################################################################################################################
##############################################################################################################################################################
####################                Francesco Innocenti, Math J.J.M. Candel, Frans E.S. Tan, Gerard J.P. van Breukelen               #########################  
# Dept. Methodology & Statistics, Care and Public Health Research Institute (CAPHRI), Maastricht University, PO Box 616, 6200MD, Maastricht, the Netherlands #
####################                Correspondence to Francesco Innocenti, E-mail: francesco.innocenti@maastrichtuniversity.nl       #########################                
##############################################################################################################################################################
##############################################################################################################################################################


# file name: sample_size_calculation.R
##############################################################################################################################################################
##############################    4. 	R code for computing the maximin TSS1 design for comparing two populations       ##################################        
##############################################################################################################################################################

##### Step 1: specify input parameters
# Below are 17 input parameters to be set by the user, either as a text file to read by R,
# or as specified in this script file.
# c2_pop1 = study cost per cluster in population 1 (e.g. France).
# c1_pop1 = study cost per subject in population 1 (e.g. France).
# c2_pop2 = study cost per cluster in population 2 (e.g. Italy).
# c1_pop2 = study cost per subject in population 2 (e.g. Italy).
# iccmax = upper bound on intraclass correlation (applies to both populations).
# psimax = upper bound on the degree of informativeness of cluster size (applies to both #populations).
# tau_pop1 = coefficient of variation of cluster size distribution in population 1 
# zeta_pop1 = skewness of cluster size distribution in population 1
# tau_pop2 = coefficient of variation of cluster size distribution in population 2 
# zeta_pop2 = skewness of cluster size distribution in population 2
# Vmax = upper bound on sum of (total unexplained outcome variance population 1 + total unexplained outcome variance population 2).
# SDratiomax = upper bound on ratio (total SD population 1 / total SD population 2), must be >=1.
# SDratiomin = lower bound on ratio (total SD population 1 / total SD population 2), must be >=1.
# and SDratiomin = 1/SDratiomax must hold if both are unequal to 1.
# mus_diff =  difference between the two population means (e.g. average alcohol consumption in #France - average alcohol consumption in Italy).
# power = desired power as a percentage.
# alpha = alpha level as a percentage.
# tail = one or two tailed testing: 1 for one-tailed, 2 for two-tailed.
# The 17 compute statements below can be bypassed if input parameters are read from a text file.
rm(list=ls(all=TRUE))
c2_pop1 <- 200 
c1_pop1 <- 10
c2_pop2 <- 200
c1_pop2 <- 10
iccmax <- 0.1
psimax <- 0.35
tau_pop1 <- 0.621  
zeta_pop1 <- 0.886 
tau_pop2 <- 0.912 
zeta_pop2 <- 1.256 
Vmax <- 2
SDratiomax <- 1
SDratiomin <- 1/SDratiomax
mus_diff <- 0.5
power <- 90 
alpha <- 5  
tail <- 2
# The text file should contain the variable names in the first row, in the order as listed above,  
# and the values of these variables in each subsequent row.
# Variable names and values should be separated by blank spaces or tabs
mydata <- read.table("d:/mydata.txt",header=TRUE)
View(mydata)
c2_pop1 <- mydata$c2_pop1  
c1_pop1 <- mydata$c1_pop1
c2_pop2 <- mydata$c2_pop2
c1_pop2 <- mydata$c1_pop2
iccmax <- mydata$iccmax
psimax <- mydata$psimax
tau_pop1 <- mydata$tau_pop1  
zeta_pop1 <- mydata$zeta_pop1 
tau_pop2 <- mydata$tau_pop2 
zeta_pop2 <- mydata$zeta_pop2 
Vmax <- mydata$Vmax
SDratiomax <- mydata$SDratiomax
SDratiomin <- 1/SDratiomax
mus_diff <- mydata$mus_diff 
power <- mydata$power
alpha <- mydata$alpha  
tail <- mydata$tail
# From these 17 input parameters the program computes the following results.
# First some intermediate results needed to compute the nr of clusters and the budget needed.
# ESmin = minimum standardized difference of interest = mus_diff/square root of (Vmax/2).
# zbeta = (1-beta)-th percentile of the standard normal distribution, where (1-beta) = power.
# zalpha = (1-alpha/2)-th percentile of the standard normal for two-tailed testing,
# and zalpha = (1-alpha)-th percentile for one-tailed testing.
# n_pop1 = maximin sample size per cluster for population 1.
# n_pop2 = maximin sample size per cluster for population 2.
# costpercluster_pop1 = total cost per cluster including costs of subjects in it, for population 1.
# costpercluster_pop2 = total cost per cluster including costs of subjects in it, for population 2.
# g_pop1 = g-function of costs, iccmax, psimax, tau_pop1, and zeta_pop1, given maximin n_pop1.
# g_pop2 = g-function of costs, iccmax, psimax, tau_pop2, and zeta_pop2, given maximin n_pop2.
# hratio = square root of (g_pop1/g_pop2).
# budgetratio = maximin budget split ratio f/(1-f) between the two populations.
# fraction = fraction of budget allocated to population 1.
# mmvnumer = numerator of maximin sampling variance of muhat_pop1-muhat_pop2,
# dividing that by the budget gives the maximin sampling variance itself.
# maxvarmusdiffhat = maximum allowable sampling variance of effect (muhat_pop1-muhat_pop2),
# given true effect mu_pop1-mu_pop2 (i.e. mus_diff), desired power and alpha.
# and then final results: total budget needed, and the number of clusters per population.
# budget = total budget needed for all study costs that are proportional to sample size.
# k_pop1 = maximin nr of clusters for population 1.
# k_pop2 = maximin nr of clusters for population 2.

# compute effect size d, and quantities depending on power and alpha
ESmin <- mus_diff/sqrt(Vmax/2)
zbeta <- qnorm((power/100),mean = 0, sd =1)
zalpha <- qnorm(1-(alpha/(100*tail)),mean = 0, sd =1)

##### Step 2:compute the max allowable sampling variance of muhat_pop1-muhat_pop2 to guarantee power:
maxvarmusdiffhat <- (mus_diff / (zbeta + zalpha))**2

##### Step 3: compute the maximin sample size per cluster per population
# and compute quantities for later steps (cost per cluster including subjects)
n_pop1 <- sqrt((c2_pop1/c1_pop1)*((1-iccmax)/iccmax)*(1/(1+psimax*(tau_pop1*(zeta_pop1-tau_pop1)+1))))
n_pop2 <- sqrt((c2_pop2/c1_pop2)*((1-iccmax)/iccmax)*(1/(1+psimax*(tau_pop2*(zeta_pop2-tau_pop2)+1))))
costpercluster_pop1 <- c2_pop1 + (c1_pop1*n_pop1)
costpercluster_pop2 <- c2_pop2 + (c1_pop2*n_pop2)

##### Step 4: compute h
g_pop1 <- (sqrt(c2_pop1*iccmax*(1+psimax*(tau_pop1*(zeta_pop1-tau_pop1)+1))) + sqrt(c1_pop1*(1-iccmax)))**2
g_pop2 <- (sqrt(c2_pop2*iccmax*(1+psimax*(tau_pop2*(zeta_pop2-tau_pop2)+1))) + sqrt(c1_pop2*(1-iccmax)))**2
hratio <- sqrt(g_pop1/g_pop2)

##### Step 5: compute maximin budget ratio  
# and maximin sampling variance of difference(muhat_pop1-muhat_pop2), apart from dividing that variance by the budget C which is still unknown.
# Note: this step is valid only if (SDratiomax = 1 or > 1) AND (SDratiomin = 1 or < 1) AND (SDratiomin = 1/SDratiomax if both are unequal to 1).
# So it covers symmetric two-sided ranges, and one-sided ranges, for the SDratio as well as homogeneous variances.
budgetratio <- 99999
budgetratio <- ifelse( ((hratio <= SDratiomax) & (hratio >= SDratiomin)), hratio**2, ifelse((hratio > SDratiomax), hratio*SDratiomax, hratio*SDratiomin))
fraction <- budgetratio/(1 + budgetratio)
mmvnumer <- 99999
mmvnumer <- ifelse( ((hratio <= SDratiomax) & (hratio >= SDratiomin)),
                    g_pop2*Vmax*(1+(hratio**2)), 
                    ifelse((hratio > SDratiomax),
                           g_pop2*Vmax*(((hratio*SDratiomax)+1)**2/((SDratiomax**2)+1)),
                           g_pop2*Vmax*(((hratio*SDratiomin)+1)**2/((SDratiomin**2) + 1))   ) )

##### Step 6: compute budget such that maximin sampling variance = max allowable sampling var from step 2
budget <- mmvnumer/maxvarmusdiffhat

##### Step 7: compute budget per population 
budget_pop1 <- fraction*budget
budget_pop2 <- (1-fraction)*budget

##### Step 8: compute the maximin number of clusters per population
k_pop1 <- budget_pop1/costpercluster_pop1
k_pop2 <- budget_pop2/costpercluster_pop2

# Make a table of all scenarios and the resulting sample sizes and budget.
myresults <- data.frame(psimax,SDratiomin,SDratiomax,iccmax,c1_pop1,c2_pop1/c1_pop1,c1_pop2,c2_pop2/c1_pop2,hratio,budgetratio,n_pop1, n_pop2,k_pop1,k_pop2,budget)
colnames(myresults) <-c("psi(max)","1/u","u","icc(max)","c1_pop1","cr_pop1","c1_pop2","cr_pop2","h","f/(1-f)","n_pop1",
                        "n_pop2","k_pop1","k_pop2","C")
library(htmlTable)
htmlTable(round(myresults,2))











