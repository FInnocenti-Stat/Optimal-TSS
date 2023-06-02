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

# file name: maximin_design.R
##############################################################################################################################################################
##############################       2. A general program in R to find the maximin parameter values       ###################################################        
##############################################################################################################################################################


######################### Function for a TSS3 design: #############################

worst_case_tau_TSS3<-function(n, rho, psi, zeta_LB, zeta_UB, eta, tau_LB, tau_UB, lambda=0.001){ 
  
  #to prevent that tau_LB>tau_UB: if the following condition does not hold, the function 
  #stops and gives a warning 
  
  if(tau_LB<=tau_UB){   
    
    #generate a fine-grained grid of values for tau
    
    tau<-seq(from=tau_LB, to=tau_UB, by=lambda ) 
    
    # the range for tau is within [0,1], then Var(mu_hat) is an increasing of zeta and only 
    #zeta_UB is used
    
    if(tau_UB<=1){                         
      
      # compute Var(mu_hat)
      
      varTSS3<-as.matrix(tau^2+1+rho*((tau^2+1)*(n-1)+n*psi*(tau^4+tau^2*(eta-3)+2*zeta_UB*tau*(1-tau^2)+1)))
      
      # find the maximum value for Var(mu_hat)
      
      MAX_varTSS3<-apply(varTSS3,MARGIN=2, FUN=max)                   
      
      #find the value of tau that maximizes Var(mu_hat)
      
      tau_wc_TSS3<-tau[which(varTSS3==max(varTSS3), arr.ind=TRUE)[1]] 
      
      # return the maximum Var(mu_hat) value and the worst-case value of tau
      
      return(cbind(MAX_varTSS3,tau_wc_TSS3))
      
      # the range for tau is within (1,+inf), then Var(mu_hat) is a decreasing function of zeta 
      #and only zeta_LB is used
      
    }else if(tau_LB>1 & tau_UB>1){      
      
      varTSS3<-as.matrix(tau^2+1+rho*((tau^2+1)*(n-1)+n*psi*(tau^4+tau^2*(eta-3)+2*zeta_LB*tau*(1-tau^2)+1)))
      
      MAX_varTSS3<-apply(varTSS3,MARGIN=2, FUN=max)                   
      
      tau_wc_TSS3<-tau[which(varTSS3==max(varTSS3), arr.ind=TRUE)[1]] 
      
      return(cbind(MAX_varTSS3,tau_wc_TSS3))
      
      # the range for tau is within [0,+inf), then Var(mu_hat) is not a monotonic function of 
      #zeta, then the function searches numerically both for the worst-case tau and zeta 
      
    }else {                               
      
      # to prevent that zeta_LB>zeta_UB
      
      if(zeta_LB<=zeta_UB){ 
        
        #generate a fine-grained grid of values for zeta
        
        zeta<-seq(from=zeta_LB, to=zeta_UB, by=lambda )
        
        zetaxtau<-as.matrix(expand.grid(zeta,tau))
        colnames(zetaxtau) <- c()
        
        varTSS3<-matrix(0,(dim(zetaxtau)[1]),1)
        
        for(i in 1:(dim(zetaxtau)[1])){
          varTSS3[i]<-(zetaxtau[i,2]^2+1+rho*((zetaxtau[i,2]^2+1)*(n-1)+n*psi*(zetaxtau[i,2]^4+zetaxtau[i,2]^2*(eta-3)+2*zetaxtau[i,1]*zetaxtau[i,2]*(1-zetaxtau[i,2]^2)+1)))
        }
        MAX_varTSS3<-apply(varTSS3,MARGIN=2, FUN=max)                  
        #find that value of tau and zeta that maximize Var(mu_hat)_TSS3
        
        zeta_wc_TSS3<-zetaxtau[which(varTSS3==max(varTSS3), arr.ind=TRUE)[1],1]
        tau_wc_TSS3<-zetaxtau[which(varTSS3==max(varTSS3), arr.ind=TRUE)[1],2]
        
        return(cbind(MAX_varTSS3,tau_wc_TSS3,zeta_wc_TSS3))
        
      } else{warning("zeta_LB must be smaller than zeta_UB")}   
    }           
    
  } else{warning("tau_LB must be smaller than tau_UB")}
  
}

######################### Function for a TSS2 design: #############################
## the procedure is the same as for the previous function, so the comments will be omitted.

worst_case_tau_TSS2<-function(n, rho, psi, zeta_LB, zeta_UB, eta, tau_LB, tau_UB, lambda=0.001){ 
  
  if(tau_LB<=tau_UB){   
    
    tau<-seq(from=tau_LB, to=tau_UB, by=lambda ) 
    
    if(tau_UB<=1){                         
      
      varTSS2<-as.matrix(1+rho*(n*(tau^2+1+psi*(tau^4+tau^2*(eta-3)+2*zeta_UB*tau*(1-tau^2)+1))-1))
      MAX_varTSS2<-apply(varTSS2,MARGIN=2, FUN=max)                   
      tau_wc_TSS2<-tau[which(varTSS2==max(varTSS2), arr.ind=TRUE)[1]] 
      
      return(cbind(MAX_varTSS2,tau_wc_TSS2))
      
    }else if(tau_LB>1 & tau_UB>1){       
      varTSS2<-as.matrix(1+rho*(n*(tau^2+1+psi*(tau^4+tau^2*(eta-3)+2*zeta_LB*tau*(1-tau^2)+1))-1))
      MAX_varTSS2<-apply(varTSS2,MARGIN=2, FUN=max)                   
      tau_wc_TSS2<-tau[which(varTSS2==max(varTSS2), arr.ind=TRUE)[1]] 
      
      return(cbind(MAX_varTSS2,tau_wc_TSS2))
      
    }else {                               
      
      if(zeta_LB<=zeta_UB){  
        zeta<-seq(from=zeta_LB, to=zeta_UB, by=lambda )
        
        zetaxtau<-as.matrix(expand.grid(zeta,tau))
        colnames(zetaxtau) <- c()
        
        varTSS2<-matrix(0,(dim(zetaxtau)[1]),1)
        
        for(i in 1:(dim(zetaxtau)[1])){
          varTSS2[i]<-(1+rho*(n*(zetaxtau[i,2]^2+1+psi*(zetaxtau[i,2]^4+zetaxtau[i,2]^2*(eta-3)+2*zetaxtau[i,1]*zetaxtau[i,2]*(1-zetaxtau[i,2]^2)+1))-1))
        }
        
        MAX_varTSS2<-apply(varTSS2,MARGIN=2, FUN=max)                   
        zeta_wc_TSS2<-zetaxtau[which(varTSS2==max(varTSS2), arr.ind=TRUE)[1],1]
        tau_wc_TSS2<-zetaxtau[which(varTSS2==max(varTSS2), arr.ind=TRUE)[1],2]
        
        return(cbind(MAX_varTSS2,tau_wc_TSS2,zeta_wc_TSS2))
        
      } else{warning("zeta_LB must be smaller than zeta_UB")}   
    }           
    
  } else{warning("tau_LB must be smaller than tau_UB")}
  
}
