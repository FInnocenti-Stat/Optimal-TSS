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

# file name: numerical_evaluation.R
##############################################################################################################################################################
###############    3. 	R code for numerically evaluate the ratio between V(mu_hat)_Delta and V(mu_hat)_L (i.e. section 3.3, Suppl. Material 1)  #############        
##############################################################################################################################################################


n<-5:100
k<-20:100   # set of scenarios with moderate skewness and kurtosis 
k<-100:150 # set of scenarios with extreme skewness and kurtosis
nxk<-as.matrix(expand.grid(n,k))
colnames(nxk) <- c()
nxk
dim(nxk)

# worst-case values
rho<-c(0.05,0.10,0.25)
psi<-c(0.067,0.35)
tau<-c(0.5,0.75,1)
zeta<-c(0.5,1,1.5)    # moderate skewness 
#zeta<-c(2,2.5,3)     # extreme skewness
eta<-c(1.5,3,4.5)      # moderate kurtosis
#eta<-c(6,9,12,15)   # extreme kurtosis

scenarios<-do.call(expand.grid, list(rho,psi,tau,zeta,eta))
colnames(scenarios)<-c("rho","psi","tau","zeta","eta")
dim(scenarios)
# Note that kurt >= skew^2+1
dim(scenarios[-c(which(scenarios[,5]<scenarios[,4]^2+1, arr.ind=TRUE)),])
scenarios<-scenarios[-c(which(scenarios[,5]<scenarios[,4]^2+1, arr.ind=TRUE)),]
dim(scenarios)
# 108 scenrios, small zeta and eta
# 162 scenrios, large zeta and eta

VarTSS3_kL<-VarTSS3_kdelta<-RE_TSS3<-matrix(0,dim(nxk)[1],dim(scenarios)[1])
VarTSS2_kL<-VarTSS2_kdelta<-RE_TSS2<-matrix(0,dim(nxk)[1],dim(scenarios)[1])

for(j in 1:(dim(scenarios)[1])){
  for(i in 1:(dim(nxk)[1])){
    # large k variance for TSS3 (equation (4))   
    VarTSS3_kL[i,j]<-scenarios$tau[j]^2+1+scenarios$rho[j]*((scenarios$tau[j]^2+1)*(nxk[i,1]-1)+nxk[i,1]*scenarios$psi[j]*(scenarios$tau[j]^4+scenarios$tau[j]^2*(scenarios$eta[j]-3)+2*scenarios$zeta[j]*scenarios$tau[j]*(1-scenarios$tau[j]^2)+1))
    # variance for TSS3 in Table 2 (main text)   
    VarTSS3_kdelta[i,j]<-((scenarios$tau[j]^2+1)/((scenarios$tau[j]^2/nxk[i,2])+1))+scenarios$rho[j]*(((scenarios$tau[j]^2+1)/((scenarios$tau[j]^2/nxk[i,2])+1))*(nxk[i,1]-1)+nxk[i,1]*scenarios$psi[j]*(((nxk[i,2]-1)/nxk[i,2])^2*scenarios$tau[j]^2*(scenarios$eta[j]-((nxk[i,2]-3)/(nxk[i,2]-1))+scenarios$tau[j]*(scenarios$tau[j]-2*scenarios$zeta[j]))+2*((nxk[i,2]-1)/nxk[i,2])*scenarios$tau[j]*(scenarios$zeta[j]-scenarios$tau[j])+1))
    # ratio Var_delta/Var_L
    RE_TSS3[i,j]<-VarTSS3_kdelta[i,j]/VarTSS3_kL[i,j]
    
    # large k variance for TSS2 (equation (3))   
    VarTSS2_kL[i,j]<-1+scenarios$rho[j]*(nxk[i,1]*(scenarios$tau[j]^2+1+scenarios$psi[j]*(scenarios$tau[j]^4+scenarios$tau[j]^2*(scenarios$eta[j]-3)+2*scenarios$zeta[j]*scenarios$tau[j]*(1-scenarios$tau[j]^2)+1))-1)
    # variance for TSS in Table 2 (main text)   
    VarTSS2_kdelta[i,j]<-1+scenarios$rho[j]*(nxk[i,1]*(((scenarios$tau[j]^2+1)/((scenarios$tau[j]^2/nxk[i,2])+1))+scenarios$psi[j]*(((nxk[i,2]-1)/nxk[i,2])^2*scenarios$tau[j]^2*(scenarios$eta[j]-((nxk[i,2]-3)/(nxk[i,2]-1))+scenarios$tau[j]*(scenarios$tau[j]-2*scenarios$zeta[j]))+2*((nxk[i,2]-1)/nxk[i,2])*scenarios$tau[j]*(scenarios$zeta[j]-scenarios$tau[j])+1))-1)
    # ratio Var_delta/Var_L
    RE_TSS2[i,j]<-VarTSS2_kdelta[i,j]/VarTSS2_kL[i,j]
  }}

summary(apply(RE_TSS3, MARGIN = 1, min)) # for each design
summary(apply(RE_TSS2, MARGIN = 1, min))

MIN_RE_TSS3<-apply(RE_TSS3, MARGIN = 1, min)
nxk[which(round(MIN_RE_TSS3,2)<0.95, arr.ind=TRUE),]

MIN_RE_TSS2<-apply(RE_TSS2, MARGIN = 1, min)
nxk[which(round(MIN_RE_TSS2,2)<0.95, arr.ind=TRUE),]

