library(survival)
library(hier.part)
library(MASS)
library(parallel)
library(invgamma)
library(truncnorm)
library(BayesLogit)
library(mvtnorm)
library(stats)
library(rmutil) # generate laplace
library(statmod) # generate inverse gaussian
true_para = as.numeric(as.matrix(read.csv("true parameters v3.csv",header=T,sep=",",fill=T))[,2])
## for propensity score
gamma_t = true_para[1:3]

## for principal stratification
beta00_t = true_para[4:6]
beta10_t = true_para[7:9]
beta11_t = true_para[10:12]

## for outcomes 
alpha00_t = true_para[13:16]
alpha10_t = true_para[17:20]
alpha01_t = true_para[21:24]
alpha11_t = true_para[25:28]
sigsq00_t = true_para[29]
sigsq01_t = true_para[30]
sigsq10_t = true_para[31]
sigsq11_t = true_para[32]


iter=11000
burn=1000

###########################################################################
######################### MCMC gibb sample ################################
###########################################################################

lassofunc = function(dataid,simdata,iter=100,burn=0){
  set.seed(dataid)
  data_obs = simdata[[dataid]]$data_obs
  N = dim(data_obs)[1]
  m = simdata[[dataid]]$m
  n = simdata[[dataid]]$n
  ui = simdata[[dataid]]$ui
  
  ########### MCMC posterior samples
  # propensity score
  gamma_pos = matrix(-99,nrow=iter,ncol=2)
  var_gamma = rep(-99,iter)
  
  # principal strata
  beta00_pos = matrix(-99,nrow=iter,ncol=2)
  beta10_pos = matrix(-99,nrow=iter,ncol=2)
  beta11_pos = matrix(-99,nrow=iter,ncol=2)
  var_beta00 = rep(-99,iter)
  var_beta10 = rep(-99,iter)
  var_beta11 = rep(-99,iter)
  
  # outcomes
  alpha00_pos = matrix(-99,nrow=iter,ncol=3)
  alpha10_pos = matrix(-99,nrow=iter,ncol=3)
  alpha01_pos = matrix(-99,nrow=iter,ncol=3)
  alpha11_pos = matrix(-99,nrow=iter,ncol=3)
  var_alpha11 = rep(-99,iter)
  var_alpha10 = rep(-99,iter)
  var_alpha00 = rep(-99,iter)
  var_alpha01 = rep(-99,iter)
  
  sigsq00_pos = rep(-99,iter)
  sigsq01_pos = rep(-99,iter)
  sigsq10_pos = rep(-99,iter)
  sigsq11_pos = rep(-99,iter)
  
  ## hyper lasso parameters
  lambdasq = matrix(-99,nrow=iter,ncol=m)
  lasso_sigsq = matrix(-99,nrow=iter,ncol=m)
  zeta = matrix(-99,nrow=iter,ncol=m)
  zeta_11 = matrix(-99,nrow=iter,ncol=m)
  zeta_10 = matrix(-99,nrow=iter,ncol=m)
  zeta_00 = matrix(-99,nrow=iter,ncol=m)
  zeta_y_11 = matrix(-99,nrow=iter,ncol=m)
  zeta_y_10 = matrix(-99,nrow=iter,ncol=m)
  zeta_y_01 = matrix(-99,nrow=iter,ncol=m)
  zeta_y_00 = matrix(-99,nrow=iter,ncol=m)
  
  ## starting value
  # propensity score
  gamma_pos[1,] = gamma_t[1:2]
  var_gamma[1] = 100
  
  # principal strata
  beta00_pos[1,] = beta00_t[1:2]
  beta10_pos[1,] = beta10_t[1:2]
  beta11_pos[1,] = beta11_t[1:2]
  var_beta00[1] = 100
  var_beta10[1] = 100
  var_beta11[1] = 100
  
  # outcomes
  alpha00_pos[1,] = alpha00_t[1:3]
  alpha10_pos[1,] = alpha10_t[1:3]
  alpha01_pos[1,] = alpha01_t[1:3]
  alpha11_pos[1,] = alpha11_t[1:3]
  var_alpha11[1] = 100
  var_alpha10[1] = 100
  var_alpha00[1] = 100
  var_alpha01[1] = 100
  
  sigsq00_pos[1] = sigsq00_t
  sigsq01_pos[1] = sigsq01_t
  sigsq10_pos[1] = sigsq10_t
  sigsq11_pos[1] = sigsq11_t
  
  ## hyper lasso parameters
  lambdasq[1,] = rgamma(m,shape=1,rate=0.1)
  lasso_sigsq[1,] = rep(100,m)
  zeta[1,] = ui*gamma_t[3]
  zeta_11[1,] = ui*beta11_t[3]
  zeta_10[1,] = ui*beta10_t[3]
  zeta_00[1,] = ui*beta00_t[3]
  zeta_y_11[1,] = ui*alpha11_t[4]
  zeta_y_10[1,] = ui*alpha10_t[4]
  zeta_y_01[1,] = ui*alpha01_t[4]
  zeta_y_00[1,] = ui*alpha00_t[4]

  ## dummy variable/indicator matrix for cluster
  cluind = data_obs$clu
  Istar = matrix(0, nrow=length(cluind), ncol=length(unique(cluind)))
  Istar[cbind(seq_along(cluind), cluind)] = 1
  
  for (k in 2:iter){
    # Step 1: update propensity score parameters
    # 1-1 : gamma
    data_gamma = cbind(data_obs$x1,data_obs$x2)
    w_gamma = apply(data_gamma%*%gamma_pos[k-1,]+Istar%*%zeta[k-1,],1,rpg,n=1,h=1)
    k_gamma = data_obs$z-1/2
    L_gamma = k_gamma/w_gamma
    omega_gamma = diag(w_gamma)
    gamma_Sigma = solve(t(data_gamma)%*%omega_gamma%*%data_gamma + 1/var_gamma[k-1]*diag(2))
    gamma_mu = gamma_Sigma%*%(t(data_gamma)%*%omega_gamma%*%(L_gamma-Istar%*%zeta[k-1,]))
    gamma_pos[k,] = mvrnorm(1,gamma_mu,gamma_Sigma)
    
    # 1-2 : var_gamma
    var_gamma[k] = rinvgamma(1,shape=2/2+0.001,rate=0.001+t(gamma_pos[k,])%*%gamma_pos[k,]/2)
    
    # 1-3 : zeta, 1,2,...m
    tau_ps = rep(-99,m)
    for (cl in 1:m){
      # generate tau_i firstly
      mean_tau_ps = sqrt(lambdasq[k-1,cl]*lasso_sigsq[k-1,cl]/zeta[k-1,cl]^2)
      invtau_ps = rinvgauss(1,mean=mean_tau_ps,dispersion = 1/lambdasq[k-1,cl])
      tau_ps[cl] = 1/invtau_ps
      # update
      data_cl = as.matrix(data_obs[data_obs$clu==cl,3:4]) #x1,x2
      n_cl = dim(data_cl)[1]
      w_cl = apply(data_cl%*%gamma_pos[k,]+zeta[k-1,cl],1,rpg,n=1,h=1)
      k_cl = data_obs[data_obs$clu==cl,]$z-1/2
      L_cl = k_cl/w_cl
      omega_cl = diag(w_cl)
      cl_Sigma = 1/(1/(lasso_sigsq[k-1,cl]*tau_ps[cl])+sum(w_cl))
      cl_mu = cl_Sigma*(rep(1,n_cl)%*%omega_cl%*%(L_cl-data_cl%*%gamma_pos[k,]))
      zeta[k,cl] = rnorm(1,mean=cl_mu,sd=sqrt(cl_Sigma))
    }
    
    # Step 2: update principal strata parameters
    # 2-0 : sample G
    data_G = cbind(data_obs$x1,data_obs$x2,Istar%*%zeta_11[k-1,],Istar%*%zeta_10[k-1,],Istar%*%zeta_00[k-1,])
    mean_11 = data_G[,1:2]%*%beta11_pos[k-1,]+ data_G[,3]
    mean_10 = data_G[,1:2]%*%beta10_pos[k-1,]+ data_G[,4]
    mean_00 = data_G[,1:2]%*%beta00_pos[k-1,]+ data_G[,5]
    p_11 = apply(mean_11,1,pnorm,mean=0,sd=1)
    p_10 = apply(mean_10,1,pnorm,mean=0,sd=1)
    p_00 = apply(mean_00,1,pnorm,mean=0,sd=1)
    
    pi_11 = 1-p_11
    pi_10 = p_11*(1-p_10)
    pi_00 = p_11*p_10*(1-p_00)
    pi_01 = p_11*p_10*p_00
    
    fi_func = function(data_yg,alpha_pos,sigsq){
      # data_yg includes y_obs,x1,x2,zeta_yg
      yobs = data_yg[1]
      data_stra = data_yg[2:3]
      fi0 = dnorm(yobs,mean=data_stra%*%alpha_pos[1:2]+data_yg[4],sd=sqrt(sigsq))
      fi1 = dnorm(yobs,mean=data_stra%*%alpha_pos[1:2]+alpha_pos[3]+data_yg[4],sd=sqrt(sigsq))
      return(c(fi0,fi1))
    }
    fi_0_00 = apply(as.matrix(cbind(data_obs$y_obs,data_obs$x1,data_obs$x2,Istar%*%zeta_y_00[k-1,])),1,fi_func,alpha00_pos[k-1,],sigsq00_pos[k-1])[1,] 
    fi_0_01 = apply(as.matrix(cbind(data_obs$y_obs,data_obs$x1,data_obs$x2,Istar%*%zeta_y_01[k-1,])),1,fi_func,alpha01_pos[k-1,],sigsq01_pos[k-1])[1,] 
    fi_0_10 = apply(as.matrix(cbind(data_obs$y_obs,data_obs$x1,data_obs$x2,Istar%*%zeta_y_10[k-1,])),1,fi_func,alpha10_pos[k-1,],sigsq10_pos[k-1])[1,] 
    fi_0_11 = apply(as.matrix(cbind(data_obs$y_obs,data_obs$x1,data_obs$x2,Istar%*%zeta_y_11[k-1,])),1,fi_func,alpha11_pos[k-1,],sigsq11_pos[k-1])[1,] 
    
    fi_1_00 = apply(as.matrix(cbind(data_obs$y_obs,data_obs$x1,data_obs$x2,Istar%*%zeta_y_00[k-1,])),1,fi_func,alpha00_pos[k-1,],sigsq00_pos[k-1])[2,] 
    fi_1_01 = apply(as.matrix(cbind(data_obs$y_obs,data_obs$x1,data_obs$x2,Istar%*%zeta_y_01[k-1,])),1,fi_func,alpha01_pos[k-1,],sigsq01_pos[k-1])[2,] 
    fi_1_10 = apply(as.matrix(cbind(data_obs$y_obs,data_obs$x1,data_obs$x2,Istar%*%zeta_y_10[k-1,])),1,fi_func,alpha10_pos[k-1,],sigsq10_pos[k-1])[2,] 
    fi_1_11 = apply(as.matrix(cbind(data_obs$y_obs,data_obs$x1,data_obs$x2,Istar%*%zeta_y_11[k-1,])),1,fi_func,alpha11_pos[k-1,],sigsq11_pos[k-1])[2,] 
    
    p_me = matrix(0,nrow=N,ncol=4) # intermediate calculation
    p_s = matrix(0,nrow=N,ncol=4)      # 4 strata probabilities for each unit
    for (i in 1:N){
      if (data_obs$z[i]==0 &data_obs$s_obs[i]==0) {
        p_me[i,4] = pi_01[i]*fi_0_01[i]
        p_me[i,3] = pi_00[i]*fi_0_00[i]
      }else if (data_obs$z[i]==0 &data_obs$s_obs[i]==1){
        p_me[i,2] = pi_10[i]*fi_0_10[i]
        p_me[i,1] = pi_11[i]*fi_0_11[i]
      }else if (data_obs$z[i]==1 &data_obs$s_obs[i]==0){
        p_me[i,3] = pi_00[i]*fi_1_00[i]
        p_me[i,2] = pi_10[i]*fi_1_10[i]
      }else if (data_obs$z[i]==1 &data_obs$s_obs[i]==1){
        p_me[i,4] = pi_01[i]*fi_1_01[i]
        p_me[i,1] = pi_11[i]*fi_1_11[i]
      }
      
      for (j in 1:4){
        p_s[i,j] = p_me[i,j]/sum(p_me[i,])
      }
    }
    
    strata_num = rep(-99,N)
    for (i in 1:N){
      if (data_obs$z[i]==1&data_obs$s_obs[i]==1){
        strata_num[i] = sample(c(1,4),size=1,replace=TRUE,prob=p_s[i,c(1,4)])
      }else if (data_obs$z[i]==1&data_obs$s_obs[i]==0){
        strata_num[i] = sample(2:3,size=1,replace=TRUE,prob=p_s[i,2:3])
      }else if (data_obs$z[i]==0&data_obs$s_obs[i]==1){
        strata_num[i] = sample(1:2,size=1,replace=TRUE,prob=p_s[i,1:2]) 
      }else if (data_obs$z[i]==0&data_obs$s_obs[i]==0){
        strata_num[i] = sample(3:4,size=1,replace=TRUE,prob=p_s[i,3:4])
      }
    }  
    cl_g_11 = Istar%*%zeta_11[k-1,]
    cl_g_10 = Istar%*%zeta_10[k-1,]
    cl_g_00 = Istar%*%zeta_00[k-1,]
    data_obs2 = as.data.frame(cbind(data_obs,strata_num,cl_g_11,cl_g_10,cl_g_00))
    
    # 2-1: sample latent variable G*
    for (i in 1:N){
      if (data_obs2$strata_num[i]==1){
        data_obs2$gstar_11[i] = rtruncnorm(1,a = -Inf,b=0,mean=mean_11[i],sd=1)
      }else if (data_obs2$strata_num[i]!=1){data_obs2$gstar_11[i] = rtruncnorm(1,a =0,b=Inf,mean=mean_11[i],sd=1)}
      if (data_obs2$strata_num[i]==2){
        data_obs2$gstar_10[i] = rtruncnorm(1,a = -Inf,b=0,mean=mean_10[i],sd=1)
      }else if (data_obs2$strata_num[i]!=2){data_obs2$gstar_10[i] = rtruncnorm(1,a =0,b=Inf,mean=mean_10[i],sd=1)}
      if (data_obs2$strata_num[i]==3){
        data_obs2$gstar_00[i] = rtruncnorm(1,a = -Inf,b=0,mean=mean_00[i],sd=1)
      }else if (data_obs2$strata_num[i]!=3){data_obs2$gstar_00[i] = rtruncnorm(1,a =0,b=Inf,mean=mean_00[i],sd=1)}
    }
    # 2-2: update beta_11,beta_10,beta_00
    #(a) beta_11
    data_pristr = cbind(data_obs2$x1,data_obs2$x2)
    beta11_Sigma = solve(1/var_beta11[k-1]*diag(2)+t(data_pristr)%*%data_pristr)
    beta11_mu = beta11_Sigma%*%(t(data_pristr)%*%(data_obs2$gstar_11-data_obs2$cl_g_11))
    beta11_pos[k,] = mvrnorm(1,beta11_mu,beta11_Sigma)
    
    
    #(b) beta_10
    datano11 = data_obs2[data_obs2$strata_num!=1,]
    data_pristr = cbind(datano11$x1,datano11$x2)
    beta10_Sigma = solve(1/var_beta10[k-1]*diag(2)+t(data_pristr)%*%data_pristr)
    beta10_mu = beta10_Sigma%*%(t(data_pristr)%*%(datano11$gstar_10-datano11$cl_g_10))
    beta10_pos[k,] = mvrnorm(1,beta10_mu,beta10_Sigma)
    
    #(c) beta_00
    datano1110 = datano11[datano11$strata_num!=2,]
    data_pristr = cbind(datano1110$x1,datano1110$x2)
    beta00_Sigma = solve(1/var_beta00[k-1]*diag(2)+t(data_pristr)%*%data_pristr)
    beta00_mu = beta00_Sigma%*%(t(data_pristr)%*%(datano1110$gstar_00-datano1110$cl_g_00))
    beta00_pos[k,] = mvrnorm(1,beta00_mu,beta00_Sigma)
    
    # 2-3 var_beta
    var_beta00[k] = rinvgamma(1,shape=2/2+0.001,rate=0.001+t(beta00_pos[k,])%*%beta00_pos[k,]/2)
    var_beta10[k] = rinvgamma(1,shape=2/2+0.001,rate=0.001+t(beta10_pos[k,])%*%beta10_pos[k,]/2)
    var_beta11[k] = rinvgamma(1,shape=2/2+0.001,rate=0.001+t(beta11_pos[k,])%*%beta11_pos[k,]/2)
    
    # 2-3 zeta_g, 1,2,...m
    #(a) zeta_11
    tau_11 = rep(-99,m)
    for (cl in 1:m){
      # generate tau_11
      mean_tau11 = sqrt(lambdasq[k-1,cl]*lasso_sigsq[k-1,cl]/(zeta_11[k-1,cl]^2))
      invtau11 = rinvgauss(1,mean=mean_tau11,dispersion = 1/lambdasq[k-1,cl])
      tau_11[cl] = 1/invtau11
      
      # update 
      data_cl_all = data_obs2[data_obs2$clu==cl,]
      n_cl = dim(data_cl_all)[1]
      data_cl = cbind(data_cl_all$x1,data_cl_all$x2) #x1,x2
      cl_Sigma11 = 1/(1/(lasso_sigsq[k-1,cl]*tau_11[cl])+n_cl)
      cl_mu11 = cl_Sigma11*(sum(data_cl_all$gstar_11-data_cl%*%beta11_pos[k,]))
      zeta_11[k,cl] = rnorm(1,mean=cl_mu11,sd=sqrt(cl_Sigma11))
    }
    
    #(b) zeta_10
    tau_10 = rep(-99,m)
    for (cl in 1:m){
      # generate tau_10
      mean_tau10 = sqrt(lambdasq[k-1,cl]*lasso_sigsq[k-1,cl]/(zeta_10[k-1,cl]^2))
      invtau10 = rinvgauss(1,mean=mean_tau10,dispersion = 1/lambdasq[k-1,cl])
      tau_10[cl] = 1/invtau10
      # update
      data_cl_all = data_obs2[data_obs2$clu==cl&data_obs2$strata_num!=1,]
      n_cl = dim(data_cl_all)[1]
      data_cl = cbind(data_cl_all$x1,data_cl_all$x2)  #x1,x2
      cl_Sigma10 = 1/(1/(lasso_sigsq[k-1,cl]*tau_10[cl])+n_cl)
      cl_mu10 = cl_Sigma10*(sum(data_cl_all$gstar_10-data_cl%*%beta10_pos[k,]))
      zeta_10[k,cl] = rnorm(1,mean=cl_mu10,sd=sqrt(cl_Sigma10))
    }
    
    #(c) zeta_00
    tau_00 = rep(-99,m)
    for (cl in 1:m){
      # generate tau_00
      mean_tau00 = sqrt(lambdasq[k-1,cl]*lasso_sigsq[k-1,cl]/(zeta_00[k-1,cl]^2))
      invtau00 = rinvgauss(1,mean=mean_tau00,dispersion = 1/lambdasq[k-1,cl])
      tau_00[cl] = 1/invtau00
      # update
      data_cl_all = data_obs2[data_obs2$clu==cl&data_obs2$strata_num==3|data_obs2$clu==cl&data_obs2$strata_num==4,]
      n_cl = dim(data_cl_all)[1]
      data_cl = cbind(data_cl_all$x1,data_cl_all$x2)  #x1,x2
      cl_Sigma00 = 1/(1/(lasso_sigsq[k-1,cl]*tau_00[cl])+n_cl)
      cl_mu00 = cl_Sigma00*(sum(data_cl_all$gstar_00-data_cl%*%beta00_pos[k,]))
      zeta_00[k,cl] = rnorm(1,mean=cl_mu00,sd=sqrt(cl_Sigma00))
    }
    
    
    ## step 3 : update outcome parameters
    # 3-0 : add data information 
    clu_y_11 = Istar%*%zeta_y_11[k-1,]
    clu_y_10 = Istar%*%zeta_y_10[k-1,]
    clu_y_00 = Istar%*%zeta_y_00[k-1,]
    clu_y_01 = Istar%*%zeta_y_01[k-1,]
    data_obs3 =as.data.frame(cbind(data_obs2,clu_y_11,clu_y_10,clu_y_00,clu_y_01))
    
    # 3-1 : alpha_g
    update_alpha_func = function(data,sigsq_alpha,sig){
      if (dim(data)[1]==1){
        var_alpha = solve(1/sigsq_alpha*diag(3)+matrix(data[,1:3])%*%data[,1:3]/sig)
        mu_alpha = var_alpha%*%(matrix(data[,1:3])%*%(data[,5]-data[,4])/sig)
        return (mvrnorm(1,mu=mu_alpha,Sigma=var_alpha))
      } else if (dim(data)[1]>1){
        var_alpha = solve(1/sigsq_alpha*diag(3)+t(data[,1:3])%*%data[,1:3]/sig)
        mu_alpha = var_alpha%*%(t(data[,1:3])%*%(data[,5]-data[,4])/sig)
        return (mvrnorm(1,mu=mu_alpha,Sigma=var_alpha))
      }
    }
    
    # (a) : alpha_11
    data_11 = as.matrix(data_obs3[data_obs3$strata_num==1,c(3,4,5,15,7)]) #x1,x2,z,zeta_y11,y_obs
    if (dim(data_11)[1]==0){
      alpha11_pos[k,] = alpha11_pos[k-1,]
    }else{
      alpha11_pos[k,] = update_alpha_func(data_11,var_alpha11[k-1],sigsq11_pos[k-1])
    }
    # (b) : alpha_10
    data_10 = as.matrix(data_obs3[data_obs3$strata_num==2,c(3,4,5,16,7)]) #x1,x2,z,zeta_y10,y_obs
    if (dim(data_10)[1]==0){
      alpha10_pos[k,] = alpha10_pos[k-1,]
    }else{
      alpha10_pos[k,] = update_alpha_func(data_10,var_alpha10[k-1],sigsq10_pos[k-1])
    }
    # (c) : alpha_00
    data_00 = as.matrix(data_obs3[data_obs3$strata_num==3,c(3,4,5,17,7)]) #x1,x2,z,zeta_y00,y_obs
    if (dim(data_00)[1]==0){
      alpha00_pos[k,] = alpha00_pos[k-1,]
    }else{
      alpha00_pos[k,] = update_alpha_func(data_00,var_alpha00[k-1],sigsq00_pos[k-1])
    }
    # (d) : alpha_01
    data_01 = as.matrix(data_obs3[data_obs3$strata_num==4,c(3,4,5,18,7)]) #x1,x2,z,zeta_y01,y_obs
    if (dim(data_01)[1]==0){
      alpha01_pos[k,] = alpha01_pos[k-1,]
    }else{
      alpha01_pos[k,] = update_alpha_func(data_01,var_alpha01[k-1],sigsq01_pos[k-1])
    }
    
    # 3-2: var_alpha_g
    var_alpha11[k] = rinvgamma(1,shape=0.001+3/2,rate=0.001+t(alpha11_pos[k,])%*%alpha11_pos[k,]/2)
    var_alpha10[k] = rinvgamma(1,shape=0.001+3/2,rate=0.001+t(alpha10_pos[k,])%*%alpha10_pos[k,]/2)
    var_alpha01[k] = rinvgamma(1,shape=0.001+3/2,rate=0.001+t(alpha01_pos[k,])%*%alpha01_pos[k,]/2)
    var_alpha00[k] = rinvgamma(1,shape=0.001+3/2,rate=0.001+t(alpha00_pos[k,])%*%alpha00_pos[k,]/2)
    
    # 3-3: sigsq_g
    ry11 = data_11[,5]-data_11[,4]-data_11[,1:3]%*%alpha11_pos[k,]
    shape_sigsq11 = (0.002+nrow(ry11))/2
    scale_sigsq11 = (0.002*1 + t(ry11)%*%ry11)/2
    sigsq11_pos[k] = rinvgamma(1,shape=shape_sigsq11,rate=scale_sigsq11)
    
    ry10 = data_10[,5]-data_10[,4]-data_10[,1:3]%*%alpha10_pos[k,]
    shape_sigsq10 = (0.002+nrow(ry10))/2
    scale_sigsq10 = (0.002*1 + t(ry10)%*%ry10)/2
    sigsq10_pos[k] = rinvgamma(1,shape=shape_sigsq10,rate=scale_sigsq10)
    
    ry00 = data_00[,5]-data_00[,4]-data_00[,1:3]%*%alpha00_pos[k,]
    shape_sigsq00 = (0.002+nrow(ry00))/2
    scale_sigsq00 = (0.002*1 + t(ry00)%*%ry00)/2
    sigsq00_pos[k] = rinvgamma(1,shape=shape_sigsq00,rate=scale_sigsq00)
    
    ry01 = data_01[,5]-data_01[,4]-data_01[,1:3]%*%alpha01_pos[k,]
    shape_sigsq01 = (0.002+nrow(ry01))/2
    scale_sigsq01 = (0.002*1 + t(ry01)%*%ry01)/2
    sigsq01_pos[k] = rinvgamma(1,shape=shape_sigsq01,rate=scale_sigsq01)
    
    
    # 3-4: zeta_yg 1,2,...m
    # (a) zeta_y11
    tau_y11 = rep(-99,m)
    for(cl in 1:m){
      # generate tau_y11
      mean_tauy11 = sqrt(lambdasq[k-1,cl]*lasso_sigsq[k-1,cl]/(zeta_y_11[k-1,cl]^2))
      invtauy11 = rinvgauss(1,mean=mean_tauy11,dispersion = 1/lambdasq[k-1,cl])
      tau_y11[cl] = 1/invtauy11
      # update
      data_cl_11 = as.matrix(data_obs3[data_obs3$clu==cl&data_obs3$strata_num==1,c(3,4,5,7)])
      n_cl_11 = dim(data_cl_11)[1]
      var_zetay11 = 1/(1/(lasso_sigsq[k-1,cl]*tau_y11[cl])+n_cl_11/sigsq11_pos[k])
      mu_zetay11 = var_zetay11*(rep(1,n_cl_11)%*%(data_cl_11[,4]-data_cl_11[,1:3]%*%alpha11_pos[k,])/sigsq11_pos[k])
      zeta_y_11[k,cl] = rnorm(1,mu_zetay11,sqrt(var_zetay11))
    }
    # (b) zeta_y10
    tau_y10 = rep(-99,m)
    for(cl in 1:m){
      # generate tau_y10
      mean_tauy10 = sqrt(lambdasq[k-1,cl]*lasso_sigsq[k-1,cl]/(zeta_y_10[k-1,cl]^2))
      invtauy10 = rinvgauss(1,mean=mean_tauy10,dispersion = 1/lambdasq[k-1,cl])
      tau_y10[cl] = 1/invtauy10
      # update
      data_cl_10 = as.matrix(data_obs3[data_obs3$clu==cl&data_obs3$strata_num==2,c(3,4,5,7)])
      n_cl_10 = dim(data_cl_10)[1]
      var_zetay10 = 1/(1/(lasso_sigsq[k-1,cl]*tau_y10[cl])+n_cl_10/sigsq10_pos[k])
      mu_zetay10 = var_zetay10*(rep(1,n_cl_10)%*%(data_cl_10[,4]-data_cl_10[,1:3]%*%alpha10_pos[k,])/sigsq10_pos[k])
      zeta_y_10[k,cl] = rnorm(1,mu_zetay10,sqrt(var_zetay10))
    }
    # (c) zeta_y00
    tau_y00 = rep(-99,m)
    for(cl in 1:m){
      # generate tau_y00
      mean_tauy00 = sqrt(lambdasq[k-1,cl]*lasso_sigsq[k-1,cl]/(zeta_y_00[k-1,cl]^2))
      invtauy00 = rinvgauss(1,mean=mean_tauy00,dispersion = 1/lambdasq[k-1,cl])
      tau_y00[cl] = 1/invtauy00
      # update
      data_cl_00 = as.matrix(data_obs3[data_obs3$clu==cl&data_obs3$strata_num==3,c(3,4,5,7)])
      n_cl_00 = dim(data_cl_00)[1]
      var_zetay00 = 1/(1/(lasso_sigsq[k-1,cl]*tau_y00[cl])+n_cl_00/sigsq00_pos[k])
      mu_zetay00 = var_zetay00*(rep(1,n_cl_00)%*%(data_cl_00[,4]-data_cl_00[,1:3]%*%alpha00_pos[k,])/sigsq00_pos[k])
      zeta_y_00[k,cl] = rnorm(1,mu_zetay00,sqrt(var_zetay00))
    }
    # (d) zeta_y01
    tau_y01 = rep(-99,m)
    for(cl in 1:m){
      # generate tau_y01
      mean_tauy01 = sqrt(lambdasq[k-1,cl]*lasso_sigsq[k-1,cl]/(zeta_y_01[k-1,cl]^2))
      invtauy01 = rinvgauss(1,mean=mean_tauy01,dispersion = 1/lambdasq[k-1,cl])
      tau_y01[cl] = 1/invtauy01
      # update
      data_cl_01 = as.matrix(data_obs3[data_obs3$clu==cl&data_obs3$strata_num==4,c(3,4,5,7)])
      n_cl_01 = dim(data_cl_01)[1]
      var_zetay01 = 1/(1/(lasso_sigsq[k-1,cl]*tau_y01[cl])+n_cl_01/sigsq01_pos[k])
      mu_zetay01 = var_zetay01*(rep(1,n_cl_01)%*%(data_cl_01[,4]-data_cl_01[,1:3]%*%alpha01_pos[k,])/sigsq01_pos[k])
      zeta_y_01[k,cl] = rnorm(1,mu_zetay01,sqrt(var_zetay01))
    }
    
    # step 4 : update hyperparameters
    # 4-1 lambdasq
    for (i in 1:m){
      lambda_rate = tau_ps[i]+tau_00[i]+tau_10[i]+tau_11[i]+tau_y00[i]+
        tau_y10[i]+tau_y11[i]+tau_y01[i]
      lambdasq[k,i] = rgamma(1,shape=8+1,rate=0.1+lambda_rate/2)
    }
    
    # 4-2 lasso_sigsq
    for(i in 1:m){
      lasso_rate = zeta[k,i]^2/tau_ps[i]+zeta_00[k,i]^2/tau_00[i]+zeta_10[k,i]^2/tau_10[i]+zeta_11[k,i]^2/tau_11[i]+
        zeta_y_00[k,i]^2/tau_y00[i]+zeta_y_01[k,i]^2/tau_y01[i]+zeta_y_10[k,i]^2/tau_y10[i]+zeta_y_11[k,i]^2/tau_y11[i]
      lasso_sigsq[k,i] = rinvgamma(1,shape=8/2,rate=lasso_rate/2)
    }
    
    print(k)
  }
  print(dataid)
  return(list(gamma_pos=gamma_pos,var_gamma=var_gamma,beta00_pos=beta00_pos,beta10_pos=beta10_pos,beta11_pos=beta11_pos,var_beta00=var_beta00,
              var_beta10=var_beta10,var_beta11=var_beta11,alpha00_pos=alpha00_pos,alpha10_pos=alpha10_pos,alpha01_pos=alpha01_pos,
              alpha11_pos=alpha11_pos,var_alpha11=var_alpha11,var_alpha10=var_alpha10,var_alpha00=var_alpha00,var_alpha01=var_alpha01,
              sigsq00_pos=sigsq00_pos,sigsq01_pos=sigsq01_pos,sigsq10_pos=sigsq10_pos,sigsq11_pos=sigsq11_pos,lambdasq=lambdasq,lasso_sigsq=lasso_sigsq,
              zeta=zeta,zeta_11=zeta_11,zeta_10=zeta_10,zeta_00=zeta_00,zeta_y_11=zeta_y_11,zeta_y_10=zeta_y_10,zeta_y_01=zeta_y_01,zeta_y_00=zeta_y_00))
}
poster_res = mclapply(1:2,lassofunc,simdata,iter=iter,burn=burn,mc.cores=1)

###########################################################################
######################### effect estimand #################################
###########################################################################

effect_func = function(dataid,simdata,poster_res,iter=iter,burn=burn){
  set.seed(dataid)
  data_obs = simdata[[dataid]]$data_obs
  N = dim(data_obs)[1]
  effect_true = simdata[[dataid]]$effect_simu
  ## recall posterior results
  gamma_pos = poster_res[[dataid]]$gamma_pos
  beta00_pos = poster_res[[dataid]]$beta00_pos
  beta10_pos = poster_res[[dataid]]$beta10_pos
  beta11_pos = poster_res[[dataid]]$beta11_pos
  alpha00_pos = poster_res[[dataid]]$alpha00_pos
  alpha10_pos = poster_res[[dataid]]$alpha10_pos
  alpha01_pos = poster_res[[dataid]]$alpha01_pos
  alpha11_pos = poster_res[[dataid]]$alpha11_pos
  zeta = poster_res[[dataid]]$zeta
  zeta_00 = poster_res[[dataid]]$zeta_00
  zeta_10 = poster_res[[dataid]]$zeta_10
  zeta_11 = poster_res[[dataid]]$zeta_11
  zeta_y_00 = poster_res[[dataid]]$zeta_y_00
  zeta_y_10 = poster_res[[dataid]]$zeta_y_10
  zeta_y_01 = poster_res[[dataid]]$zeta_y_01
  zeta_y_11 = poster_res[[dataid]]$zeta_y_11
  
  ## dummy variable/indicator matrix for cluster
  cluind = data_obs$clu
  Istar = matrix(0, nrow=length(cluind), ncol=length(unique(cluind)))
  Istar[cbind(seq_along(cluind), cluind)] = 1
  
  gamma_ave = apply(gamma_pos[(burn+1):iter,],2,mean)
  beta00_ave = apply(beta00_pos[(burn+1):iter,],2,mean)
  beta10_ave = apply(beta10_pos[(burn+1):iter,],2,mean)
  beta11_ave = apply(beta11_pos[(burn+1):iter,],2,mean)
  alpha11_ave = apply(alpha11_pos[(burn+1):iter,],2,mean)
  alpha10_ave = apply(alpha10_pos[(burn+1):iter,],2,mean)
  alpha01_ave = apply(alpha01_pos[(burn+1):iter,],2,mean)
  alpha00_ave = apply(alpha00_pos[(burn+1):iter,],2,mean)
  zeta_ave = apply(zeta[(burn+1):iter,],2,mean)
  zeta_00_ave =  apply(zeta_00[(burn+1):iter,],2,mean)
  zeta_10_ave =  apply(zeta_10[(burn+1):iter,],2,mean)
  zeta_11_ave =  apply(zeta_11[(burn+1):iter,],2,mean)
  zeta_y_00_ave = apply(zeta_y_00[(burn+1):iter,],2,mean)
  zeta_y_10_ave = apply(zeta_y_10[(burn+1):iter,],2,mean)
  zeta_y_01_ave = apply(zeta_y_01[(burn+1):iter,],2,mean)
  zeta_y_11_ave = apply(zeta_y_11[(burn+1):iter,],2,mean)
  
  data_G = cbind(data_obs$x1,data_obs$x2,Istar%*%zeta_11_ave,Istar%*%zeta_10_ave,Istar%*%zeta_00_ave)
  mean_11 = data_G[,1:2]%*%beta11_ave+ data_G[,3]
  mean_10 = data_G[,1:2]%*%beta10_ave+ data_G[,4]
  mean_00 = data_G[,1:2]%*%beta00_ave+ data_G[,5]
  p_11 = apply(mean_11,1,pnorm,mean=0,sd=1)
  p_10 = apply(mean_10,1,pnorm,mean=0,sd=1)
  p_00 = apply(mean_00,1,pnorm,mean=0,sd=1)
  
  pi_11_est = 1-p_11
  pi_10_est = p_11*(1-p_10)
  pi_00_est = p_11*p_10*(1-p_00)
  pi_01_est = p_11*p_10*p_00
  
  data_0_y00 = as.matrix(cbind(data_obs$x1,data_obs$x2,rep(0,N),Istar%*%zeta_y_00_ave))
  data_1_y00 = as.matrix(cbind(data_obs$x1,data_obs$x2,rep(1,N),Istar%*%zeta_y_00_ave))
  Ey_0_00 = data_0_y00[,1:3]%*%alpha00_ave + data_0_y00[,4]
  Ey_1_00 = data_1_y00[,1:3]%*%alpha00_ave + data_1_y00[,4]
  ace00 = (Ey_1_00-Ey_0_00)*pi_00_est
  
  data_0_y10 = as.matrix(cbind(data_obs$x1,data_obs$x2,rep(0,N),Istar%*%zeta_y_10_ave))
  data_1_y10 = as.matrix(cbind(data_obs$x1,data_obs$x2,rep(1,N),Istar%*%zeta_y_10_ave))
  Ey_0_10 = data_0_y10[,1:3]%*%alpha10_ave + data_0_y10[,4]
  Ey_1_10 = data_1_y10[,1:3]%*%alpha10_ave + data_1_y10[,4]
  ace10 = (Ey_1_10-Ey_0_10)*pi_10_est
  
  data_0_y01 = as.matrix(cbind(data_obs$x1,data_obs$x2,rep(0,N),Istar%*%zeta_y_01_ave))
  data_1_y01 = as.matrix(cbind(data_obs$x1,data_obs$x2,rep(1,N),Istar%*%zeta_y_01_ave))
  Ey_0_01 = data_0_y01[,1:3]%*%alpha01_ave + data_0_y01[,4]
  Ey_1_01 = data_1_y01[,1:3]%*%alpha01_ave + data_1_y01[,4]
  ace01 = (Ey_1_01-Ey_0_01)*pi_01_est
  
  data_0_y11 = as.matrix(cbind(data_obs$x1,data_obs$x2,rep(0,N),Istar%*%zeta_y_11_ave))
  data_1_y11 = as.matrix(cbind(data_obs$x1,data_obs$x2,rep(1,N),Istar%*%zeta_y_11_ave))
  Ey_0_11 = data_0_y11[,1:3]%*%alpha11_ave + data_0_y11[,4]
  Ey_1_11 = data_1_y11[,1:3]%*%alpha11_ave + data_1_y11[,4]
  ace11 = (Ey_1_11-Ey_0_11)*pi_11_est
  
  ave_ate = mean(ace00+ace01+ace10+ace11)
  
  ## estimated 
  ate_est = rep(-99,(iter-burn))
  for (i in (burn+1):iter){
    j=i-burn
    
    data_G = cbind(data_obs$x1,data_obs$x2,Istar%*%zeta_11[i,],Istar%*%zeta_10[i,],Istar%*%zeta_00[i,])
    mean_11 = data_G[,1:2]%*%beta11_pos[i,]+ data_G[,3]
    mean_10 = data_G[,1:2]%*%beta10_pos[i,]+ data_G[,4]
    mean_00 = data_G[,1:2]%*%beta00_pos[i,]+ data_G[,5]
    p_11 = apply(mean_11,1,pnorm,mean=0,sd=1)
    p_10 = apply(mean_10,1,pnorm,mean=0,sd=1)
    p_00 = apply(mean_00,1,pnorm,mean=0,sd=1)
    
    pi_11_est = 1-p_11
    pi_10_est = p_11*(1-p_10)
    pi_00_est = p_11*p_10*(1-p_00)
    pi_01_est = p_11*p_10*p_00
    
    alpha00_est = alpha00_pos[i,]
    alpha01_est = alpha01_pos[i,]
    alpha10_est = alpha10_pos[i,]
    alpha11_est = alpha11_pos[i,]
    
    data_0_y00 = as.matrix(cbind(data_obs$x1,data_obs$x2,rep(0,N),Istar%*%zeta_y_00[i,]))
    data_1_y00 = as.matrix(cbind(data_obs$x1,data_obs$x2,rep(1,N),Istar%*%zeta_y_00[i,]))
    Ey_0_00 = data_0_y00[,1:3]%*%alpha00_est + data_0_y00[,4]
    Ey_1_00 = data_1_y00[,1:3]%*%alpha00_est + data_1_y00[,4]
    ace00 = (Ey_1_00-Ey_0_00)*pi_00_est
    
    data_0_y10 = as.matrix(cbind(data_obs$x1,data_obs$x2,rep(0,N),Istar%*%zeta_y_10[i,]))
    data_1_y10 = as.matrix(cbind(data_obs$x1,data_obs$x2,rep(1,N),Istar%*%zeta_y_10[i,]))
    Ey_0_10 = data_0_y10[,1:3]%*%alpha10_est + data_0_y10[,4]
    Ey_1_10 = data_1_y10[,1:3]%*%alpha10_est + data_1_y10[,4]
    ace10 = (Ey_1_10-Ey_0_10)*pi_10_est
    
    data_0_y01 = as.matrix(cbind(data_obs$x1,data_obs$x2,rep(0,N),Istar%*%zeta_y_01[i,]))
    data_1_y01 = as.matrix(cbind(data_obs$x1,data_obs$x2,rep(1,N),Istar%*%zeta_y_01[i,]))
    Ey_0_01 = data_0_y01[,1:3]%*%alpha01_est + data_0_y01[,4]
    Ey_1_01 = data_1_y01[,1:3]%*%alpha01_est + data_1_y01[,4]
    ace01 = (Ey_1_01-Ey_0_01)*pi_01_est
    
    data_0_y11 = as.matrix(cbind(data_obs$x1,data_obs$x2,rep(0,N),Istar%*%zeta_y_11[i,]))
    data_1_y11 = as.matrix(cbind(data_obs$x1,data_obs$x2,rep(1,N),Istar%*%zeta_y_11[i,]))
    Ey_0_11 = data_0_y11[,1:3]%*%alpha11_est + data_0_y11[,4]
    Ey_1_11 = data_1_y11[,1:3]%*%alpha11_est + data_1_y11[,4]
    ace11 = (Ey_1_11-Ey_0_11)*pi_11_est
    
    ate_est[j] = mean(ace00+ace01+ace10+ace11)
  }
  library(coda)
  ate_mcmc = as.mcmc(ate_est)
  hpd_ate_mcmc = HPDinterval(ate_mcmc, prob = 0.95)
  res_ate = as.data.frame(cbind(effect_true,mean(ate_est),sd(ate_est),hpd_ate_mcmc[1],hpd_ate_mcmc[2]))
  colnames(res_ate) = c("True","Mean","SD","2.5%","97.5%")
  print(dataid)
  return (list(ate_est=ate_est,ate_mcmc=ate_mcmc,hpd_ate_mcmc=hpd_ate_mcmc,res_ate=res_ate,ave_ate=ave_ate,effect_true=effect_true))
}
effect_res = mclapply(1:2,effect_func,simdata,poster_res,iter=iter,burn=burn,mc.cores=1)

