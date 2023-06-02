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

# file name: simulation_study.R
##############################################################################################################################################################
#############################################       1. R code of the simulation study       ##################################################################        
##############################################################################################################################################################

rm(list=ls(all=T))
set.seed(05092016)
################################       Functions        ###################################
# Population Mean Estimator for TSS3
mu_hat_TSS3<-function(cluster_size,cluster_mean){ 
  
  Num<-cluster_size%*%cluster_mean
  Den<-sum(cluster_size)
  
  return(Num/Den)
}
# Population Mean Estimator for TSS2
mu_hat_TSS2<-function(p,cluster_size,cluster_mean){ 
  
  Cluster_Size_p<-p*cluster_size
  Num_2<-Cluster_Size_p%*%cluster_mean
  Den_2<-sum(Cluster_Size_p)
  
  return(Num_2/Den_2)
}
# Marginal Expectation for the Population Mean Estimator
marg.exp_mu_hat<-function(b0,Gamma,muN,cvN,k){ 
  b0+Gamma*((k-1)/k)*muN*cvN^2
}

# Expectation Conditional Variance for TSS3
exp.cond.var.TSS3<-function(rho,n,k,cvN,s2y){   
  return(s2y*((rho*(n-1)+1)/(n*k))*((cvN^2+1)/(1+(cvN^2/k))))
}
# Expectation Conditional Variance for TSS2
exp.cond.var.TSS2<-function(rho,n,k,cvN,s2y){ 
  return(s2y*((rho*(n*((cvN^2+1)/(1+(cvN^2/k)))-1)+1)/(n*k)))
}
# Variance Conditional Expectation for TSS2/TSS3
var.cond.exp.TSS23<-function(cvN, skewN, kurtN, rho, rhouN, k,s2y){   
  return(s2y*(rhouN^2/(1-rhouN^2))*(rho/k)*(((k-1)/k)^2*cvN^2*(kurtN-((k-3)/(k-1))+cvN*(cvN-2*skewN))+2*((k-1)/k)*cvN*(skewN-cvN)+1))
}
# Marginal variance for TSS3
marg.var_TSS3<-function(cvN, skewN, kurtN, rho, rhouN, k,s2y,n){  
  marg.tss3<-s2y*((((rho*(n-1)+1)/(n*k))*((cvN^2+1)/(1+(cvN^2/k))))+((rhouN^2/(1-rhouN^2))*(rho/k)*(((k-1)/k)^2*cvN^2*(kurtN-((k-3)/(k-1))+cvN*(cvN-2*skewN))+2*((k-1)/k)*cvN*(skewN-cvN)+1)))
  return(marg.tss3)
}
# Marginal variance for TSS2
marg.var_TSS2<-function(cvN, skewN, kurtN, rho, rhouN, k,s2y,n){ 
  marg.tss2<-s2y*(((rho*(n*((cvN^2+1)/(1+(cvN^2/k)))-1)+1)/(n*k))+((rhouN^2/(1-rhouN^2))*(rho/k)*(((k-1)/k)^2*cvN^2*(kurtN-((k-3)/(k-1))+cvN*(cvN-2*skewN))+2*((k-1)/k)*cvN*(skewN-cvN)+1)))
  return(marg.tss2)
}
# To sample from a uniform distribution
rdu<-function(q,u) {sample(1:u,q,replace=T)}
#############################     Setting of the parameters  ##########################
beta0<-10  
k<-20  # 30, 40, 100
n<-20 # 20, 50, 100
# Intraclass Correlation Coefficient
rho<-c(0.05,0.30) 
# Correlation between cluster size and cluster effect                  
rho_uN<-c(0.25,0.5)            # 0.75     
# Standard deviation of cluster effect
sigma_v<-0.5                    
# Standard deviation of the individual effect    
sigma_e<-sigma_v*sqrt((1-rho)/rho)  
#########################      Cluster size distributions   ################################
######### Negative Binomial   
# Average cluster size
theta_N<-400
r<-4           # 1 for cv=1
# Standard Deviation of cluster size
sigma_N<-sqrt(theta_N*(1+(theta_N/r)))
# Coefficient of variation of cluster size
# Note the denominator theta_N+2, because the distribution is N=X+2, where X is a Neg.Bin
tau_N<-round(sigma_N/(theta_N+2),3)
# Skewness of cluster size
zeta_N<-round((1/sqrt(r+theta_N))*(sqrt(r/theta_N)+2*sqrt(theta_N/r)),3)
# Kurtosis of cluster size
eta_N<-round(3+(6/r)+(1/sigma_N^2),3)

######### Discrete Uniform   
u<-800  
# Average cluster size 
theta_N<-(u+2)/2
# Standard Deviation of cluster size
sigma_N<-sqrt(((u-1)^2-1)/12)
# Coefficient of variation of cluster size
tau_N<-round(sqrt(((u-1)^2-1)/(u+2)^2)*(1/sqrt(3)),3)
# Skewness of cluster size
zeta_N<-0
# Kurtosis of cluster size
eta_N<-round(3-((6*((u-1)^2+1))/(5*((u-1)^2-1))),3)

######### Four-paramters Beta     
## 3 scenarios: a=1.1 b=1; a=4, b=2; a=b=2
min<-2 
max<-800 
a<-1.1   #4, 2 
b<-1     #2, 2
# Average cluster size 
theta_N<-min+(max-min)*(a/(a+b))
# Standard Deviation of cluster size
sigma_N<-sqrt((a*b)/(a+b+1))*((max-min)/(a+b))
# Coefficient of variation of cluster size
tau_N<-round(sigma_N/theta_N,3)
# Skewness of cluster size
zeta_N<-round(sqrt((a+b+1)/(a*b))*2*((b-a)/(a+b+2)),3)
# Kurtosis of cluster size
eta_N<-round(3*((a+b+1)*(2*(a+b)^2+a*b*(a+b-6)))/(a*b*(a+b+2)*(a+b+3)),3)

######### Bimodal     
#  Y= w X1 + (1-w) X2, 
# where X1=Beta 4-param.(min=2, max=1000, a=2, b=4), X2=Normal(mu=1000, sd=50)
# X1= beta 4-parameters
min<-2 ;max<-1000; a<-2; b<-4     
m1<-min+(max-min)*(a/(a+b));s1<-((a*b)/(a+b+1))*((max-min)/(a+b))^2
# X2= Normal
m2=1000;s2<-50^2
#  scenario 1: w=0.9,  scenario 2: w=0.45 
w1=0.90 #.45 
w2=1-w1
# Average cluster size 
theta_N<-round(w1*m1+w2*m2)
# Standard Deviation of cluster size
sigma_N<-sqrt((w1*s1+w2*s2)+(w1*m1^2+w2*m2^2)-theta_N^2)
# Coefficient of variation of cluster size
tau_N<-round(sigma_N/theta_N,3)
# Skewness of cluster size
zeta_N<-round((w1*((m1-theta_N)^3+3*(m1-theta_N)*s1)+w2*((m2-theta_N)^3+3*(m2-theta_N)*s2))/(sigma_N)^(3),3)
# Kurtosis of cluster size
eta_N<-round((w1*((m1-theta_N)^4+6*(m1-theta_N)^2*s1+3*s1^2)+w2*((m2-theta_N)^4+6*(m2-theta_N)^2*s2+3*s2^2))/(sigma_N)^4,3)

#########    GP List Size distribution in England     
# Available from: https://digital.nhs.uk/data-and-information/publications/statistical/patients-registered-at-a-gp-practice/patients-registered-at-a-gp-practice-october-2017-special-topic-practice-list-size-comparison-october-2013-to-october-2017
# Patients Registered at a GP Practice - October 2017: Totals (GP practice-all persons) [csv, size: 541.2 kB]
library(moments)
library(readr)
gp_reg_pat_prac_all_oct_17 <- read_csv("d:\\gp-reg-pat-prac-all-oct-17.csv")

dim(gp_reg_pat_prac_all_oct_17)
head(gp_reg_pat_prac_all_oct_17)

GP.size<-as.data.frame(gp_reg_pat_prac_all_oct_17[,10]) # GP list size distribution
sum(GP.size)     # 58719921 patients
dim(GP.size)[1]  # 7353 GPs

summary(GP.size)
# Average cluster size 
theta_N<-round(mean(as.numeric(unlist(GP.size))),0)        #  7986
# Standard Deviation of cluster size
sigma_N<-round(sd(as.numeric(unlist(GP.size))),3)
# Coefficient of variation of cluster size
tau_N<-round(sd(as.numeric(unlist(GP.size)))/mean(as.numeric(unlist(GP.size))),3)#  0.6331486
# Skewness of cluster size
zeta_N<-round(skewness(as.numeric(unlist(GP.size))) ,3) #  2.120307
# Kurtosis of cluster size
eta_N<-round(kurtosis(as.numeric(unlist(GP.size))),3)  #  14.54937


#########    High School Size distribution in Italy 
# Available from: http://dati.istruzione.it/opendata/opendata/catalogo/elements1/?area=Studenti
# "Studenti per anno di corso e fascia di eta'. Scuola statale."> "CSV" > "Distribuzione per ANNOSCOLASTICO 201617"
library(moments)
library(readr)
italian_schools <- read.csv("d:\\ALUCORSOETASTA20161720170831.csv")
# Data set for primary school, lower secondary school, high school(="SCUOLA SECONDARIA II GRADO")
dim(italian_schools)
head(italian_schools)

schools<-data.frame(italian_schools[-1],stringsAsFactors = FALSE)
schools[order(schools$CODICESCUOLA),]

###### Retriving High Schools
high_school<-data.frame(schools[schools$ORDINESCUOLA=="SCUOLA SECONDARIA II GRADO",],stringsAsFactors = FALSE)
dim(high_school)
high_school<-high_school[order(high_school$CODICESCUOLA),]

high_school_id<-matrix(0,dim(high_school)[1],1)
high_school_id[1,]<-1

for(i in 2:dim(high_school)[1]){ if(high_school[i,1]==high_school[i-1,1])
{    high_school_id[i,]<-high_school_id[i-1,]  }
  else
  {    high_school_id[i,]<-high_school_id[i-1,]+1  }}

tail(high_school_id)  # Number of High Schools = 6235
cbind(high_school_id,high_school)

high_school_size<-tapply(high_school$ALUNNI,high_school_id,sum)
sum(high_school_size) # Number of students= 2515060
summary(high_school_size)

# Average cluster size 
theta_N<-round(mean(high_school_size),0)       #  403
sigma_N<-sd(high_school_size)
# Coefficient of variation of cluster size
tau_N<-round(sigma_N/theta_N,3)                #  0.912
# Skewness of cluster size
zeta_N<-round(skewness(high_school_size) ,3)   #  1.256
# Kurtosis of cluster size
eta_N<-round(kurtosis(high_school_size),3)     #  4.315

##############################     To store the results     ###############################
# ratio of average number of individuals sampled per cluster to population average cluster size 
p<-n/theta_N       
gamma<-mu<-matrix(0,1,4)
# number of samples from cluster size distribution
sim1<-5000 
# number of samples conditional on cluster size distribution
sim2<-500 
# each column is a different draw from the cluster size distribution
cluster_sizes<-matrix(0,k,sim1)       
# each column is a different draw from the cluster size distribution,
# each array is a different combination between rho and rho_uN
MU_TSS3<-MU_TSS2<-array(0, dim=c(sim2,sim1,4))       
# conditional expectation and variance
cond.exp_TSS3<-cond.exp_TSS2<-cond.var_TSS3<-cond.var_TSS2<-array(0, dim=c(1,sim1,4))  

###########################       Simulation    #########################################
########## Step 1: sample k cluster sizes from the cluster size distribution 
# choose one distribution:
### Negative Binomial
# 5000 draws from cluster size distribution 
cluster_sizes<-replicate(sim1,rnbinom(n=k,mu = theta_N, size=r),simplify =T) +2  
# E(N) = E(Negative Binomial + 2)=E(Neg.Bin)+2
theta_N<-theta_N+2 

### Discrete Uniform 
# 5000 draws from cluster size distribution 
cluster_sizes<-replicate(sim1,rdu(q=k,u=u),simplify =T)  
### Four-parameter Beta 
# 5000 draws from cluster size distribution in the population
for(i in 1:sim1){    
  # generate k values from a Beta(a,b) in the range (0,1)  
  X<-rbeta(n=k,shape1=a,shape2=b)    
  # linear transformation which yields the Beta(a,b,min,max)
  cluster_sizes[,i]<-(max-min)*X+min 
}

### Bimodal 
w1<-0.9 # 0.45
w2<-1-w1
indicators<-matrix(0,k,sim1)
for(t in 1:sim1){
  weights<-rmultinom(n=k,size=1,prob=c(w1,w2))
  for(i in 1:k){
    for(j in 1:2){
      if(weights[j,i]==1)
      {indicators[i,t]<-j}  }  }
  # generate from beta distribution
  X1<-rbeta(n=length(which(indicators[,t]==1)),shape1=a,shape2=b) 
  Z1<-(max-min)*X1+min           
  # generate from Normal distribution
  X2<-rnorm(n=length(which(indicators[,t]==2)),mean=m2,sd=sqrt(s2))    
  # 5000 draws from cluster size distribution 
  cluster_sizes[,t]<-c(Z1,X2)  
}

### GP List Size distribution in England 
# 5000 draws from cluster size distribution 
cluster_sizes<-replicate(sim1,sample(GP.size[,1], size=k,replace = T),simplify =T)  

### High School Size distribution in Italy 
# 5000 draws from cluster size distribution 
cluster_sizes<-replicate(sim1,sample(high_school_size, size=k,replace = T),simplify =T)  


########## Step 2: sample n individuals per cluster, generate their Ys, compute mu  

#Discretize cluster size distribution 
cluster_sizes<-round(cluster_sizes) 

for(l in 1:length(rho)){
  for(t in 1:length(rho_uN)){
    # Gamma
    gamma[,ifelse(l==1,l*t,l+t)]<-sqrt(rho_uN[t]^2/(1-rho_uN[t]^2))*(sigma_v/sigma_N)      
    # Population Mean (equation (2))
    mu[,ifelse(l==1,l*t,l+t)]<-beta0+gamma[,ifelse(l==1,l*t,l+t)]*theta_N*tau_N^2            
    
    # 'for loop' to cycle across the 5000 draws from cluster size distribution 
    for(i in 1:sim1){ 
      cat("l=", l, "t=", t, "i=", i,"\n")
      # 500 samples of n individuals from the k clusters conditional on the sampled cluster sizes
      for(j in 1:sim2){ 
        
        ##2.a: generate the k cluster effects 
        cluster_effect_v<-rnorm(n=k, mean=0, sd=sigma_v)  
        
        ## 2.b: sample n individuals per sampled cluster (i.e. generate the individual effects epsilons)  
        Eps_TSS3<-Eps_TSS2<-list()  
        # number of individuals sampled per cluster
        Num_Indpc_TSS2<-NULL  
        for(s in 1:k){          
          # TSS3: sample n individuals in any cases  
          Eps_TSS3[[s]]<-rnorm(n=n, mean=0, sd=sigma_e[l])     
          
          # TSS2: 
          ifelse(round(p*cluster_sizes[s,i])==0,Num_Indpc_TSS2<-1,Num_Indpc_TSS2<-round(p*cluster_sizes[s,i]))
          Eps_TSS2[[s]]<-rnorm(n=Num_Indpc_TSS2, mean=0, sd=sigma_e[l]) # TSS2
          
        }
        
        TSS3<-data.frame(rep(1:k,each=n),rep(cluster_sizes[,i],each=n), rep(cluster_effect_v,each=n),unlist(Eps_TSS3))   
        
        TSS2<-data.frame(rep(1:k,times=lapply(Eps_TSS2,length)),rep(cluster_sizes[,i],times=lapply(Eps_TSS2,length)), rep(cluster_effect_v,times=lapply(Eps_TSS2,length)),unlist(Eps_TSS2))                          
        
        colnames(TSS3)<-colnames(TSS2)<-c("ID_Clsuter","N","v","epsilon")
        
        ## 2.c: generate the Ys  
        # Equation (1) with Assumption 4
        Y_TSS3=beta0+gamma[,ifelse(l==1,l*t,l+t)]*(TSS3$N-theta_N)+TSS3$v+TSS3$epsilon     
        Y_TSS2=beta0+gamma[,ifelse(l==1,l*t,l+t)]*(TSS2$N-theta_N)+TSS2$v+TSS2$epsilon    
        
        ## 2.d: compute mu 
        MU_TSS3[j,i,ifelse(l==1,l*t,l+t)]<-mu_hat_TSS3(cluster_size=cluster_sizes[,i],cluster_mean=as.matrix(tapply(Y_TSS3,TSS3$ID_Clsuter,mean)))     
        
        MU_TSS2[j,i,ifelse(l==1,l*t,l+t)]<-mu_hat_TSS2(p=p,cluster_size=cluster_sizes[,i],cluster_mean=as.matrix(tapply(Y_TSS2,TSS2$ID_Clsuter,mean))) 
      }
      
      rm(cluster_effect_v,Eps_TSS3,Eps_TSS2,Num_Indpc_TSS2,TSS3,TSS2,Y_TSS3,Y_TSS2)
      gc()
      
      ########## Step 3: compute the Conditional Expectation and Conditional Variance of mu_hat  
      
      ## Compute the Conditional Expectation   
      cond.exp_TSS3[,i,ifelse(l==1,l*t,l+t)]<-mean(MU_TSS3[,i,ifelse(l==1,l*t,l+t)])      
      cond.exp_TSS2[,i,ifelse(l==1,l*t,l+t)]<-mean(MU_TSS2[,i,ifelse(l==1,l*t,l+t)])      
      
      ##  Compute the Conditional Variance  
      cond.var_TSS3[,i,ifelse(l==1,l*t,l+t)]<-sum((MU_TSS3[,i,ifelse(l==1,l*t,l+t)]-mean(MU_TSS3[,i,ifelse(l==1,l*t,l+t)]))^2)/(sim2-1)  # TSS3
      cond.var_TSS2[,i,ifelse(l==1,l*t,l+t)]<-sum((MU_TSS2[,i,ifelse(l==1,l*t,l+t)]-mean(MU_TSS2[,i,ifelse(l==1,l*t,l+t)]))^2)/(sim2-1)  # TSS2
    }  }  }

########## Step 4: compute the Marginal Expectation and Variance of mu_hat 
marg.exp_TSS3<-marg.exp_TSS2<-expect.cond.var_TSS3<-expect.cond.var_TSS2<-matrix(0,1,4)
variance.cond.exp_TSS3<-variance.cond.exp_TSS2<-marginal.var_TSS3<-marginal.var_TSS2<-matrix(0,1,4)
RB_marg.exp<-RB_marg.exp_TSS3<-RB_marg.exp_TSS2<-RB_exp.cond.var_TSS3<-RB_exp.cond.var_TSS2<-matrix(0,1,4)
RB_var.cond.exp_TSS3<-RB_var.cond.exp_TSS2<-RB_marg.var_TSS3<-RB_marg.var_TSS2<-matrix(0,1,4)

## Compute the Marginal Expectation 
marg.exp_TSS3<-apply(cond.exp_TSS3,3,mean)       
marg.exp_TSS2<-apply(cond.exp_TSS2,3,mean)       

##Compute the Expectation of the Conditional Variance  
expect.cond.var_TSS3<-apply(cond.var_TSS3,3,mean)       
expect.cond.var_TSS2<-apply(cond.var_TSS2,3,mean)      

##Compute the Variance of the Conditional Expectation 
for(i in 1:4){
  variance.cond.exp_TSS3[,i]<-sum((cond.exp_TSS3[,,i]-marg.exp_TSS3[i])^2)/(sim1-1)   
  variance.cond.exp_TSS2[,i]<-sum((cond.exp_TSS2[,,i]-marg.exp_TSS2[i])^2)/(sim1-1)     
}



## Compute the Marginal Variance 
marginal.var_TSS3<-expect.cond.var_TSS3+variance.cond.exp_TSS3       
marginal.var_TSS2<-expect.cond.var_TSS2+variance.cond.exp_TSS2       



###########################     Relative Bias (RB)     ####################################
###    RB Marginal Expectation (i.e. unbiasedness)  
RB_marg.exp<-(marg.exp_mu_hat(b0=beta0,Gamma=gamma,muN=theta_N,cvN=tau_N,k=k)-mu)/mu
RB_marg.exp<-t(RB_marg.exp)
RB_marg.exp_TSS3<-(marg.exp_mu_hat(b0=beta0,Gamma=gamma,muN=theta_N,cvN=tau_N,k=k)-marg.exp_TSS3)/marg.exp_TSS3
RB_marg.exp_TSS3<-t(RB_marg.exp_TSS3)
RB_marg.exp_TSS2<-(marg.exp_mu_hat(b0=beta0,Gamma=gamma,muN=theta_N,cvN=tau_N,k=k)-marg.exp_TSS2)/marg.exp_TSS2
RB_marg.exp_TSS2<-t(RB_marg.exp_TSS2)

###    RB Expectation of the Conditional Variance 
RB_exp.cond.var_TSS3<-(rep(exp.cond.var.TSS3(rho=rho,n=n,k=k,cvN=tau_N,s2y=sigma_v^2+sigma_e^2),each=2)-expect.cond.var_TSS3)/expect.cond.var_TSS3
RB_exp.cond.var_TSS3<-as.matrix(RB_exp.cond.var_TSS3)
RB_exp.cond.var_TSS2<-(rep(exp.cond.var.TSS2(rho=rho,n=n,k=k,cvN=tau_N,s2y=sigma_v^2+sigma_e^2),each=2)-expect.cond.var_TSS2)/expect.cond.var_TSS2
RB_exp.cond.var_TSS2<-as.matrix(RB_exp.cond.var_TSS2)

###    RB Variance of the Conditional Expectation 
RB_var.cond.exp_TSS3<-(rep(var.cond.exp.TSS23(cvN=tau_N, skewN=zeta_N, kurtN=eta_N, rho=rho, rhouN=rho_uN, k=k,s2y=sigma_v^2+sigma_e^2),times=2)-variance.cond.exp_TSS3)/variance.cond.exp_TSS3
RB_var.cond.exp_TSS3<-t(RB_var.cond.exp_TSS3)
RB_var.cond.exp_TSS2<-(rep(var.cond.exp.TSS23(cvN=tau_N, skewN=zeta_N, kurtN=eta_N, rho=rho, rhouN=rho_uN, k=k,s2y=sigma_v^2+sigma_e^2),times=2)-variance.cond.exp_TSS2)/variance.cond.exp_TSS2
RB_var.cond.exp_TSS2<-t(RB_var.cond.exp_TSS2)

###    RB Marginal Variance   
MV.TSS3<-MV.TSS2<-matrix(0,1,length(rho)*length(rho_uN))

for(i in 1:length(rho)){
  for(s in 1:length(rho_uN)){
    MV.TSS3[,ifelse(i==1,i*s,i+s)]<-marg.var_TSS3(cvN=tau_N, skewN=zeta_N, kurtN=eta_N, rho=rho[i], rhouN=rho_uN[s], k=k,s2y=sigma_v^2+sigma_e[i]^2,n=n)  
    MV.TSS2[,ifelse(i==1,i*s,i+s)]<-marg.var_TSS2(cvN=tau_N, skewN=zeta_N, kurtN=eta_N, rho=rho[i], rhouN=rho_uN[s], k=k,s2y=sigma_v^2+sigma_e[i]^2,n=n)       
  }}

RB_marg.var_TSS3<-(MV.TSS3-marginal.var_TSS3)/marginal.var_TSS3
RB_marg.var_TSS3<-t(RB_marg.var_TSS3)
RB_marg.var_TSS2<-(MV.TSS2-marginal.var_TSS2)/marginal.var_TSS2
RB_marg.var_TSS2<-t(RB_marg.var_TSS2)

############################ Reporting the results ####################################
library(htmlTable) 
TSS3_results<-round(data.frame(RB_marg.exp,RB_marg.exp_TSS3,cbind(marg.exp_TSS3),
                               RB_exp.cond.var_TSS3,cbind(expect.cond.var_TSS3),
                               RB_var.cond.exp_TSS3,t(variance.cond.exp_TSS3),
                               RB_marg.var_TSS3,t(marginal.var_TSS3)),3)
rownames(TSS3_results)<-c("ICC=0.05 & corr(u,N)=0.25","ICC=0.05 & corr(u,N)=0.5","ICC=0.30 & corr(u,N)=0.25","ICC=0.30 & corr(u,N)=0.5")
colnames(TSS3_results)<-c("Rel.Bias(1) of mu_hat","Rel.Bias(2) of mu_hat","mu_hat",
                          "Rel.Bias of E(V(mu_hat|N))","E(V(mu_hat|N))",
                          "Rel.Bias of V(E(mu_hat|N))","V(E(mu_hat|N))",
                          "Rel.Bias of V(mu_hat)","V(mu_hat)")
htmlTable(TSS3_results)


TSS2_results<-round(data.frame(RB_marg.exp,RB_marg.exp_TSS2,cbind(marg.exp_TSS2),
                               RB_exp.cond.var_TSS2,cbind(expect.cond.var_TSS2),
                               RB_var.cond.exp_TSS2,t(variance.cond.exp_TSS2),
                               RB_marg.var_TSS2,t(marginal.var_TSS2)),3)
rownames(TSS2_results)<-c("ICC=0.05 & corr(u,N)=0.25","ICC=0.05 & corr(u,N)=0.5","ICC=0.30 & corr(u,N)=0.25","ICC=0.30 & corr(u,N)=0.5")
colnames(TSS2_results)<-c("Rel.Bias(1) of mu_hat","Rel.Bias(2) of mu_hat","mu_hat",
                          "Rel.Bias of E(V(mu_hat|N))","E(V(mu_hat|N))",
                          "Rel.Bias of V(E(mu_hat|N))","V(E(mu_hat|N))",
                          "Rel.Bias of V(mu_hat)","V(mu_hat)")
htmlTable(TSS2_results)
