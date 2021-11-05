#setwd("/Users/LJY/Dropbox/PHD/horeshoe/Rfiles")
source("packages.R")
source("pos_horseshoe_nowcast.R")


#setwd("/Users/LJY/Dropbox/PHD/horeshoe/condo_insample_new")
simu<-readRDS("simu_48_beta_5.rds")
factors_initial_raw<-readRDS("simu_48_beta_5_ft.rds")

x_final<-simu$X
y_final<-simu$y

G=1000
M=8000
alpha_s<-2
beta_s<-6
mu_theta<-0



###########
n=dim(x_final)[1]

n1=15
n2=n-2*n1

I3<-diag(n)

I11<-diag(n1)
I12<-matrix(0,nrow=n1,ncol=(n-n1))
I1<-cbind(I11,I12)

I21<-diag(2*n1)
I22<-matrix(0,nrow=2*n1,ncol=(n-2*n1))
I2<-cbind(I21,I22)


set_i<-function(i){ #i=1,2
  if (i==1){
    I=I1
    set0=seq(n1+1,n)
  }else if(i==2){
    I=I2
    set0=seq(2*n1+1,n)
  }
  return(list("i_matrix"=I,"set0"=set0))
}

 
##for rel=3

simu_draw_all<-function(q,m,rel){
  
      index1<-set_i(1)
      set1<-unlist(index1$set0)
      index2<-set_i(2)
      set2<-unlist(index2$set0)
     
      
      x_old<-x_final[,1:((39+q)*3+m)]
      y_old<-c(y_final[1:(40+q-1)])
      x_old[set1,dim(x_old)[2]]<-10^10
      x_old[set2,dim(x_old)[2]-1]<-10^10
    
  t<-dim(x_old)[2]
  n<-dim(x_old)[1]
  k<-length(y_old)
  r<-6
  x_start<-x_old[,-c(t-1,t)]
  y_start<-y_old
  
  X<- as.matrix(x_start)
  
  Z<- matrix(0,nrow=n,ncol=t-2)
  for (i in 1:n){
    Z[i,]=(X[i,]-mean(X[i,]))/sd(X[i,])
  }
  
  
  
  S<- cor(t(Z))
  ###
  # dim(S)
  # dim(cor(Z))
  
  # s1<- cor(t(X))
  # table(round(s1,4)==round(S,4))
  
  # mean(Z[1,])
  # sd(Z[1,])
  
  
  R<- rARPACK::eigs(S,r)
  V<- R$vectors[,1:r]
  D<- diag(R$values[1:r])
  t(V)%*%V
  
  tilde.f<- matrix(0,nrow=r,ncol=t-2)
  for (i in 1:(t-2)){
    tilde.f[,i]<- t(V)%*%Z[,i]
  }
  
  theta_initial<- V
  omega_initial<- diag(diag(S-V%*%D%*%t(V)))
  
  
  
  tmp1<- matrix(0,nrow=r,ncol=r)
  tmp2<- matrix(0,nrow=r,ncol=r)
  tmp3<- matrix(0,nrow=r,ncol=r)
  for (i in 2:(t-2)){
    tmp1<- tmp1+tilde.f[,i-1]%*%t(tilde.f[,i-1])
    tmp2<- tmp2+tilde.f[,i]%*%t(tilde.f[,i-1])
    tmp3<- tmp3+tilde.f[,i]%*%t(tilde.f[,i])
  }
  
  A_initial<-tmp2%*%solve(tmp1)
  
  sigma_initial<- tmp3/(t-2)-A_initial%*%(tmp1/(t-2))%*%t(A_initial)
  
  # initials are: A_initial, theta_initial, sigma_initial, omega_initial
  
  #-------------------------------------
  
  # Kalman Filters to get initials of factors
  
  m0<- matrix(rep(0,r),ncol=1)
  C0_fun<- function(Th,O){
    solve(t(Th)%*%Th)%*%t(Th)%*%O%*%Th%*%solve(t(Th)%*%Th)
  }
  C0<- C0_fun(theta_initial,omega_initial)
  
  # dim(C0)
  # C1<- diag(c(10000,10000,10000))
  
  kal_fil_fun<- function(m,c,Th,O,A,S,D){
    t=dim(D)[2]
    r=length(m)
    factor<- matrix(0,nrow=r,ncol=t)
    for (i in 1:t){
      a<- A%*%m
      R<- A%*%c%*%t(A)+S
      e<- D[,i]-Th%*%a
      Q<- Th%*%R%*%t(Th)+O
      K<- R%*%t(Th)%*%solve(Q)
      m<- a+K%*%e
      c<- R-R%*%t(Th)%*%solve(Q)%*%Th%*%R
      factor[,i]<- m
    }
    return(factor)
  }
  
  #factors_initial<- kal_fil_fun(m=m0,c=C0,Th=theta_initial,O=omega_initial,A=A_initial,S=sigma_initial,D=Z)
  
  factor_old<-factors_initial_raw[,1:t]
  factors_initial=factor_old
  
  theta_old=theta_initial
  a_old=diag(rep(0.9,r))
  sigma_old=diag(seq(6,1))
  omega_old=omega_initial
  beta_old<-rep(2,(3*r+2))
  mu_old<-rep(2,n)
  
  theta.res<- matrix(0,nrow=G,ncol=n*r)
  mu.res<- matrix(0,nrow=G,ncol=n)
  omega.res<- matrix(0,nrow=G,ncol=n*n)
  a.res<- matrix(0,nrow=G,ncol=r)
  sigma.res<- matrix(0,nrow=G,ncol=r)
  
  beta0_old<-beta_old[1]
  beta1_old<-beta_old[2:(r+1)]
  beta2_old<-beta_old[(r+2):(2*r+1)]
  beta3_old<-beta_old[(2*r+2):(3*r+1)]
  alpha_old<-beta_old[3*r+2]
  
  beta.res<- matrix(0,nrow=2+3*r,ncol=G)
  
  eta2_old<-0.1
  eta2.res<-rep(0,G)
  
  tau_old=1
  lambda_old=rep(1,r)
  lambda.res<-matrix(0,nrow=G,ncol=r)
  
  factor1.res<- matrix(0,nrow=G,ncol=t)
  factor2.res<- matrix(0,nrow=G,ncol=t)
  factor3.res<- matrix(0,nrow=G,ncol=t)
  factor4.res<- matrix(0,nrow=G,ncol=t)
  factor5.res<- matrix(0,nrow=G,ncol=t)
  factor6.res<- matrix(0,nrow=G,ncol=t)
  
  y_pred<-rep(0,G)
  
  
  for (j in 1:(G+M)){
    
    
    theta_new<- theta_update_lag(data=x_old,factor=factor_old,Omega=omega_old,
                                 mu_x=mu_old,rel=rel)
    
    theta_old<- theta_new
    if (j>M) theta.res[j-M,]<- c(theta_new)
    
  
    
    beta_new<-beta_update(factor=factor_old,GDP=y_old,e2=eta2_old,
                          lambda=lambda_old,tau=tau_old)
    
    
    beta_old<- beta_new
    beta0_old<-beta_old[1]
    beta1_old<-beta_old[2:(r+1)]
    beta2_old<-beta_old[(r+2):(2*r+1)]
    beta3_old<-beta_old[(2*r+2):(3*r+1)]
    alpha_old<-beta_old[3*r+2]
    if (j>M) beta.res[,j-M]<-beta_new
    
    
    eta2_new<-eta2_update(factor=factor_old,GDP=y_old,beta=beta_old,
                          alpha_h=2,beta_h=0.0001)
    eta2_old<-eta2_new
    if (j>M) eta2.res[j-M]<-eta2_new
    
    Sig_A_new<-Sig_A_update(factor=factor_old,Sigma=sigma_old,A=a_old,
                            alpha_s=alpha_s, beta_s=beta_s)
    A_new<-Sig_A_new$A
    A_old<- A_new
    if (j>M) a.res[j-M,]<-diag(A_new)
    
    sigma_new<-Sig_A_new$Sigma
    sigma_old<-sigma_new
    if (j>M) sigma.res[j-M,]<-diag(sigma_new)
    
    
    mu_new<-mu_update_lag(omega=omega_old,rel=rel,data=x_old,
                          factor=factor_old,theta=theta_old)
    mu_old<- mu_new
    if (j>M) mu.res[j-M,]<-mu_new
    
    omega_new<- omega_update_lag(data=x_old,factor=factor_old,Omega=omega_old,
                                 Theta=theta_old,rel=rel,mu_x=mu_old)
    omega_old<- omega_new
    if (j>M) omega.res[j-M,]<-c(omega_new)
    
    # tau_new<-tau_update(theta=theta_old,lambda=lambda_old,tau=tau_old)
    # tau_old<- tau_new
    # if (j>M) tau.res[j-M]<-tau_new
    
    
    lambda_new<-lambda_update(lambda=lambda_old,tau=tau_old,beta=beta_old,
                              nu_lambda = 0.7)
    lambda_old<- lambda_new$lambda
    if (j>M) lambda.res[j-M,]<-as.vector(lambda_new$lambda)
    #if (j>M) jump_lam[j-M]<-as.vector(lambda_new$jump_lam)
    
    
    
    ##seperate S and Theta
    factor_new<-ft_update_lag(data=x_old,GDP=y_old,mu_x=mu_old,factor=factor_old,Theta=theta_old,
                              A=a_old,Omega=omega_old,Sigma=sigma_old,eta2=eta2_old,
                              beta0=beta0_old,beta1=beta1_old,beta2=beta2_old,
                              beta3=beta3_old,rel=rel,alpha=alpha_old)
    
    
    factor_old<- factor_new
    if (j>M) factor1.res[j-M,]<- factor_new[1,]
    if (j>M) factor2.res[j-M,]<- factor_new[2,]
    if (j>M) factor3.res[j-M,]<- factor_new[3,]
    if (j>M) factor4.res[j-M,]<- factor_new[4,]
    if (j>M) factor5.res[j-M,]<- factor_new[5,]
    if (j>M) factor6.res[j-M,]<- factor_new[6,]
    
    
    ###get nowcasts fit##
    
    if (j>M){
      t<-dim(factor_new)[2]
      
      if ((m==1)&&(rel<3)){
        y_new<-beta0_old+t(beta1_old)%*%A_new^3%*%factor_new[,t]+t(beta2_old)%*%A_new^2%*%factor_new[,t]+
          t(beta3_old)%*%A_new%*%factor_new[,t]+alpha_old*y_old[length(y_old)]
      }else if ((m==1)&&(rel==3)){
        y_new<-beta0_old+t(beta1_old)%*%A_new^2%*%factor_new[,t]+t(beta2_old)%*%A_new%*%factor_new[,t]+
          t(beta3_old)%*%factor_new[,t]+alpha_old*y_old[length(y_old)]
      }else if ((m==2)&&(rel<3)){
        y_new<-beta0_old+t(beta1_old)%*%A_new^2%*%factor_new[,t]+t(beta2_old)%*%A_new%*%factor_new[,t]+
          t(beta3_old)%*%factor_new[,t]+alpha_old*y_old[length(y_old)]
      }else if ((m==2)&&(rel==3)){
        y_new<-beta0_old+t(beta1_old)%*%A_new%*%factor_new[,t]+t(beta2_old)%*%factor_new[,t]+
          t(beta3_old)%*%factor_new[,t-1]+alpha_old*y_old[length(y_old)]
      }else if ((m==3)&&(rel<3)){
        y_new<-beta0_old+t(beta1_old)%*%A_new%*%factor_new[,t]+t(beta2_old)%*%factor_new[,t]+
          t(beta3_old)%*%factor_new[,t-1]+alpha_old*y_old[length(y_old)]
      }else if ((m==3)&&(rel==3)) {
        y_new<-beta0_old+t(beta1_old)%*%factor_new[,t]+t(beta2_old)%*%factor_new[,t-1]+
          t(beta3_old)%*%factor_new[,t-2]+alpha_old*y_old[length(y_old)]
      }
      y_pred[j-M]<-y_new
    }
    
    #---------------
    if (j%%100==0) cat(j, "\n")
  }
  
  
  
  result_all<-list("beta_all"=beta.res,"theta_all"=theta.res,
                   "lambda_all"=lambda.res,"factor1_all"=factor1.res,
                   "factor2_all"=factor2.res,"factor3_all"=factor3.res,
                   "factor4_all"=factor4.res,"factor5_all"=factor5.res,
                   "factor6_all"=factor6.res,"a_all"=a.res,
                   "sigma_all"=sigma.res,"ft_grs"=factors_initial,
                   "mu_all"=mu.res,"eta2_all"=eta2.res,
                   "y_pred"=y_pred)
  
  names<-paste("simu1_m",m,"_q",q,"_rel",rel,".rds",sep="")
  
  saveRDS(result_all, names)
  
  return(result_all)
}


df<-expand.grid(q=seq(1,20),m=c(1,2,3),rel=c(1,2))
df_sep<-df[1:16,]


result.i<-mcmapply(simu_draw_all,m=df_sep$m,rel=df_sep$rel,q=df_sep$q,mc.cores=16)

saveRDS(result.i,file="simu1_9.rds")

