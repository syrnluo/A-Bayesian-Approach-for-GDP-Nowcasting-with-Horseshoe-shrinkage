
alpha_h=2
beta_h=0.0001


###########################Theta#########
####no cons2
theta_update_lag<- function(data,factor,Omega,mu_x,rel,I,I1,I2){
  #tau scalar
  #lambda n*1 vector
  
  x=as.matrix(data)
  n=dim(x)[1]
  t_null=dim(x)[2]
  
  if (rel==1){t=t_null+1}else{t=t_null}
  
  
  mu<-as.vector(mu_x)
  omg<-as.matrix(Omega)
  
  FT=as.matrix(factor)
  r=dim(FT)[1]
  
  ifnot1<-as.numeric(rel!=1)
  
  tmp1<-matrix(0,nrow=n*r,ncol=n*r)
  tmp2<-rep(0,n*r)
  for (i in 1:(t-3+ifnot1)){
    C_i<-kronecker(diag(n),t(FT[,i]))
    tmp1=tmp1+t(C_i)%*%solve(omg)%*%C_i
    tmp2=tmp2+t(C_i)%*%solve(omg)%*%(x[,i]-mu)
  }
  
  if (rel==1){
    omg_t1<-I1%*%omg%*%t(I1)
    omg_t2<-I2%*%omg%*%t(I2)
    Ct1<-kronecker(diag(n),t(FT[,t-1]))
    Ct2<-kronecker(diag(n),t(FT[,t-2]))
    
    W=tmp1+t(Ct1)%*%t(I1)%*%solve(omg_t1)%*%I1%*%Ct1+
      t(Ct2)%*%t(I2)%*%solve(omg_t2)%*%I2%*%Ct2+diag(n*r)
    
    inv_W<-solve(W,tol=1e-20)
    
    U=inv_W%*%(tmp2+t(Ct1)%*%t(I1)%*%solve(omg_t1)%*%I1%*%(x[,t-1]-mu)+
                 t(Ct2)%*%t(I2)%*%solve(omg_t2)%*%I2%*%(x[,t-2]-mu))
  
    }else { #rel=2/3
    omg_t1<-I1%*%omg%*%t(I1)
    omg_t<-I%*%omg%*%t(I)
    Ct1<-kronecker(diag(n),t(FT[,t-1]))
    Ct<-kronecker(diag(n),t(FT[,t]))
    
    
    W=tmp1+t(Ct)%*%t(I)%*%solve(omg_t)%*%I%*%Ct+
      t(Ct1)%*%t(I1)%*%solve(omg_t1)%*%I1%*%Ct1+diag(n*r)
    
    inv_W<-solve(W,tol=1e-20)
    
    U=inv_W%*%(tmp2+t(Ct)%*%t(I)%*%solve(omg_t)%*%I%*%(x[,t]-mu)+
                 t(Ct1)%*%t(I1)%*%solve(omg_t1)%*%I1%*%(x[,t-1]-mu))
    
  }
  
  the=MASS::mvrnorm(1,U,inv_W,tol = 1e-20)
  
  res<-matrix(the,nrow=n,ncol=r,byrow = T)
  
  return(res)
}





#--------------------------------------------------------

Sig_A_update<-function(factor,Sigma,A,alpha_s,beta_s){
  FT=as.matrix(factor)
  r=dim(FT)[1]
  t=dim(FT)[2]
  A<-as.matrix(A)
  Sig<-as.matrix(Sigma)
  alpha<-alpha_s
  beta<-beta_s
  #######j=1###
  tmp1_1<-0
  tmp2_1<-0
  tmp3_1<-0
  for (i in 2:t){
    tmp1_1<-tmp1_1+FT[1,i-1]^2
    tmp2_1<-tmp2_1+FT[1,i-1]*FT[1,i]
    tmp3_1<-tmp3_1+(FT[1,i]-A[1,1]*FT[1,i-1])^2
  }
  inv_sig_1<-1+1/Sig[1,1]*tmp1_1
  mu_1<-tmp2_1/(Sig[1,1]+tmp1_1)
  alpha_sig_1<-alpha+t/2-0.5 #(t-1)/2
  beta_sig_1<-beta+0.5*tmp3_1
  l1<-0
  repeat{
    A11<-rnorm(1,mu_1, sqrt(1/inv_sig_1))
    Sig11<-MCMCpack::rinvgamma(n=1,shape=alpha_sig_1, scale=beta_sig_1)
    l1<-l1+1
    if (((abs(A11)<1)&((Sig11/(1-A11^2))>(Sig[2,2]/(1-A[2,2]^2))))|(l1>3000)){
      break
    }
  }
  if (l1>3000){
    A[1,1]<-A[1,1]
    Sig[1,1]<-Sig[1,1]
  }else{
    A[1,1]<-A11
    Sig[1,1]<-Sig11
  }
  
  ######j=2 to（r-1）####
  for (j in 2:(r-1)){
    tmp1<-0
    tmp2<-0
    tmp3<-0
    for (i in 2:t){
      tmp1<-tmp1+FT[j,i-1]^2
      tmp2<-tmp2+FT[j,i-1]*FT[j,i]
      tmp3<-tmp3+(FT[j,i]-A[j,j]*FT[j,i-1])^2
    }
    inv_sig_j<-1+1/Sig[j,j]*tmp1
    mu_j<-tmp2/(Sig[j,j]+tmp1)
    alpha_sig_j<-alpha+t/2-0.5
    beta_sig_j<-beta+0.5*tmp3
    l2<-0
    repeat{
      Ajj<-rnorm(1,mu_j, sqrt(1/inv_sig_j))
      Sigjj<-MCMCpack::rinvgamma(n=1,shape=alpha_sig_j, scale=beta_sig_j)
      l2<-l2+1
      if (((abs(Ajj)<1)&
           (Sigjj/(1-Ajj^2)<Sig[j-1,j-1]/(1-A[j-1,j-1]^2))&
           (Sigjj/(1-Ajj^2)>Sig[j+1,j+1]/(1-A[j+1,j+1]^2)))|(l2>3000)){
        break
      }
    }
    if (l2>3000){
      A[j,j]<-A[j,j]
      Sig[j,j]<-Sig[j,j]
    }else{
      A[j,j]<-Ajj
      Sig[j,j]<-Sigjj
    }
  }
  #######j=r###
  tmp1_r<-0
  tmp2_r<-0
  tmp3_r<-0
  for (i in 2:t){
    tmp1_r<-tmp1_r+FT[r,i-1]^2
    tmp2_r<-tmp2_r+FT[r,i-1]*FT[r,i]
    tmp3_r<-tmp3_r+(FT[r,i]-A[r,r]*FT[r,i-1])^2
  }
  inv_sig_r<-1+1/Sig[r,r]*tmp1_r
  mu_r<-tmp2_r/(Sig[r,r]+tmp1_r)
  alpha_sig_r<-alpha+t/2-0.5
  beta_sig_r<-beta+0.5*tmp3_r
  l3<-0
  repeat{
    Arr<-rnorm(1,mu_r, sqrt(1/inv_sig_r))
    Sigrr<-MCMCpack::rinvgamma(n=1,shape=alpha_sig_r, scale=beta_sig_r)
    l3<-l3+1
    if (((abs(Arr)<1)&((Sigrr/(1-Arr^2))<(Sig[r-1,r-1]/(1-A[r-1,r-1]^2))))|(l3>3000)){
      break
    }
  }
  if (l3>3000){
    A[r,r]<-A[r,r]
    Sig[r,r]<-Sig[r,r]
  }else{
    A[r,r]<-Arr
    Sig[r,r]<-Sigrr
  }
  
  
  Sig_A <- list("Sigma"=Sig, "A" =A)
  return(Sig_A)
}


#----------------------------------

omega_update_lag<- function(data,factor,Theta,rel,mu_x,Omega,I,I1,I2){
  #mu_x be n*1 vector
  x=as.matrix(data)
  n=dim(x)[1]
  t_null=dim(x)[2]
  if (rel==1){t=t_null+1}else{t=t_null}
  
  Omg<-as.matrix(Omega)
  
  FT=as.matrix(factor)
  r=dim(FT)[1]
  The=as.matrix(Theta)
  
  theta<-c(t(The))
  
  ifnot1<-as.numeric(rel!=1)
  
  tmp=matrix(0,nrow=n,ncol=n)
  
  for (i in 1:(t-3+ifnot1)){
    C=kronecker(diag(rep(1,n)),t(FT[,i]))
    tmp=tmp+(x[,i]-mu_x-C%*%theta)%*%t(x[,i]-mu_x-C%*%theta)
  }
  
  P=tmp+1/n*diag(n)
  o_old=Omg #nu_theta=n+2##
  o_new=MCMCpack::riwish(v=t+n,S=P) ##t-2+nu_theta
  
  if (rel==1){
    CT1=kronecker(diag(rep(1,n)),t(FT[,t-1]))
    CT2=kronecker(diag(rep(1,n)),t(FT[,t-2]))
    
    H1<-I1%*%(x[,t-1]-mu_x-CT1%*%theta)%*%t((x[,t-1]-mu_x-CT1%*%theta))%*%t(I1)
    H2<-I2%*%(x[,t-2]-mu_x-CT2%*%theta)%*%t((x[,t-2]-mu_x-CT2%*%theta))%*%t(I2)
    
    p1<-((det(I1%*%o_new%*%t(I1))*det(I2%*%o_new%*%t(I2)))/
           (det(I1%*%o_old%*%t(I1))*det(I2%*%o_old%*%t(I2))))^-0.5
    
    p2<-exp(-0.5*sum(diag(H2%*%(solve(I2%*%o_new%*%t(I2))-solve(I2%*%o_old%*%t(I2)))))-
              0.5*sum(diag(H1%*%(solve(I1%*%o_new%*%t(I1))-solve(I1%*%o_old%*%t(I1))))))
  }else{ #rel==2/3
    CT1=kronecker(diag(rep(1,n)),t(FT[,t-1]))
    CT=kronecker(diag(rep(1,n)),t(FT[,t]))
    
    H1<-I1%*%(x[,t-1]-mu_x-CT1%*%theta)%*%t((x[,t-1]-mu_x-CT1%*%theta))%*%t(I1)
    HT<-I%*%(x[,t]-mu_x-CT%*%theta)%*%t((x[,t]-mu_x-CT%*%theta))%*%t(I)
    
    p1<-((det(I1%*%o_new%*%t(I1))*det(I%*%o_new%*%t(I)))/
           (det(I1%*%o_old%*%t(I1))*det(I%*%o_old%*%t(I))))^-0.5
    
    p2<-exp(-0.5*sum(diag(H1%*%(solve(I1%*%o_new%*%t(I1))-solve(I1%*%o_old%*%t(I1)))))-
              0.5*sum(diag(HT%*%(solve(I%*%o_new%*%t(I))-solve(I%*%o_old%*%t(I))))))
  }
  
  p_all<-p1*p2
  if (runif(1)<=p_all){
    o_old<-o_new
  }else{
    o_old=o_old
  }
  return(o_old)
}

#--------------------------------------------

beta_update <- function(factor,GDP,e2,lambda,tau){
  FT=as.matrix(factor)
  t<-dim(FT)[2]
  r=dim(FT)[1]
  
  y=as.matrix(GDP)
  k=length(y)
  
  S<-tau^-2*diag(lambda^-2)
  
  
  tmp1=matrix(0,nrow=3*r+2,ncol=3*r+2)
  tmp2=matrix(0,nrow=3*r+2,ncol=1)
  for (i in 2:k){
    
    F.n<-matrix(c(1,t(FT[,3*i]),t(FT[,3*i-1]),t(FT[,3*i-2]),y[i-1]),ncol=1)
    y.n<-y[i]
    
    tmp1=tmp1+F.n%*%t(F.n)
    tmp2=tmp2+F.n*y.n
  }
  
  vcov1<-dbind(S,S)
  vcov2<-dbind(1,vcov1)
  vcov3<-dbind(S,1)
  vcov<-dbind(vcov2,vcov3)
  W=tmp1/e2+vcov
  inv_W<-solve(W,tol=1e-20)
  U=inv_W%*%(tmp2/e2)
  beta=MASS::mvrnorm(1,U,inv_W)
  return(matrix(beta,ncol=1))
}

#mm<-beta_update(factor=simu$FT,GDP=simu$y,e2=simu$eta2,lambda=simu$lambda,tau=simu$tau)

#-------------------------------------------

##without S matrix

eta2_update<- function(factor,GDP,beta,alpha_h,beta_h){
  FT=as.matrix(factor)
  y=as.matrix(GDP)
  k=length(y)
  beta<-as.vector(beta)
  tmp=0
  for (i in 2:k){
    F.n<-c(1,t(FT[,3*i]),t(FT[,3*i-1]),t(FT[,3*i-2]),y[i-1])
    tmp=tmp+(y[i]-t(beta)%*%F.n)^2
  }
  
  beta_e2=beta_h+tmp/2
  
  eta2=pscl::rigamma(n=1,alpha=alpha_h+(k-1)/2,beta=beta_e2)
  
  return(eta2)
}


######factors########


ft_update_lag<- function(data,GDP,mu_x,factor,Theta,A,Omega,Sigma,eta2,rel,
                        beta0,beta1,beta2,beta3,alpha,I,I1,I2){
  #mu_x n*1 vector
  ##eta2 is 1*1, beta is (3r+2)*1##
  x=as.matrix(data)
  n=dim(x)[1]
  t_null=dim(x)[2]
  
  if (rel==1){t=t_null+1}else{t=t_null}
  
  y=as.matrix(GDP)
  k<-length(y)
  
  FT=as.matrix(factor)
  r=dim(FT)[1]
  
  
  Theta=as.matrix(Theta)
  A=as.matrix(A)
  Sigma<-as.matrix(Sigma)
  Omega<-as.matrix(Omega)
  
  
  tilde_A=solve(A)
  s_tmp=t(A)%*%solve(Sigma)%*%A
  p<-dim(s_tmp)[1]
  ##when t=1##
  startF1=matrix(0,nrow=r,ncol=1)
  tilde_Y_1=c(x[,1]-mu_x,FT[,2])
  tilde_X_1=rbind(Theta,A)
  s_tilde_Sigma_1=cbind(rbind(solve(Omega),matrix(0,nrow=r,ncol=n)),
                        rbind(matrix(0,nrow=n,ncol=r),solve(Sigma)))
  W1=t(tilde_X_1)%*%s_tilde_Sigma_1%*%tilde_X_1
  M1=t(tilde_X_1)%*%s_tilde_Sigma_1%*%tilde_Y_1
  inv_W1=solve(W1,tol=1e-90)
  F1=MASS::mvrnorm(1,inv_W1%*%M1,inv_W1)
  startF1[,1]=F1
  tmpF=F1
  ##when t=2,3##
  startF=matrix(0,nrow=r,ncol=2)
  for (i in 2:3){
    tilde_Y_2=c(x[,i]-mu_x,FT[,(i+1)],tmpF)
    tilde_X_2=rbind(Theta,A,tilde_A)
    s_tilde_Sigma_2=cbind(rbind(solve(Omega),matrix(0,nrow=(r+p),ncol=n)),
                          rbind(matrix(0,nrow=n,ncol=r),solve(Sigma),matrix(0,nrow=p,ncol=r)),
                          rbind(matrix(0,nrow=(n+r),ncol=p),s_tmp))
    W.i=t(tilde_X_2)%*%s_tilde_Sigma_2%*%tilde_X_2
    M.i=t(tilde_X_2)%*%s_tilde_Sigma_2%*%tilde_Y_2
    inv_Wi=solve(W.i,tol=1e-90)
    F.i=MASS::mvrnorm(1, inv_Wi%*%M.i,inv_Wi)
    startF[,(i-1)]<-F.i
    tmpF=F.i
  }
  
  midF=matrix(0,nrow=r,ncol=3*k-3-3)
  for (i in 4:(3*k-3)){
    if (i%%3==1){
      y0<-y[(i+2)/3]-alpha*y[(i+2)/3-1]-beta0-t(beta1)%*%FT[,i+2]-t(beta2)%*%FT[,i+1]
      tilde_Y=c(x[,i]-mu_x,FT[,i+1],tmpF,y0)
      tilde_X=rbind(Theta,A,tilde_A,t(beta3))
      s_tilde_Sigma=cbind(rbind(solve(Omega),matrix(0,nrow=(r+p+1),ncol=n)),
                          rbind(matrix(0,nrow=n,ncol=r),solve(Sigma),matrix(0,nrow=(p+1),ncol=r)),
                          rbind(matrix(0,nrow=(n+r),ncol=p),s_tmp,matrix(0,nrow=1,ncol=p)),
                          rbind(matrix(0,nrow=(n+r+p),ncol=1),eta2^-1))
      W=t(tilde_X)%*%s_tilde_Sigma%*%tilde_X
      M=t(tilde_X)%*%s_tilde_Sigma%*%tilde_Y
      inv_W=solve(W,tol=1e-90)
      F_i=MASS::mvrnorm(1,inv_W%*%M,inv_W)
      midF[,i-3]<-F_i
      tmpF=F_i
      
    }else if(i%%3==2){
      y0<-y[(i+1)/3]-alpha*y[(i+1)/3-1]-beta0-t(beta1)%*%FT[,i+1]-t(beta3)%*%FT[,i-1]
      tilde_Y=c(x[,i]-mu_x,FT[,i+1],tmpF,y0)
      tilde_X=rbind(Theta,A,tilde_A,t(beta2))
      s_tilde_Sigma=cbind(rbind(solve(Omega),matrix(0,nrow=(r+p+1),ncol=n)),
                          rbind(matrix(0,nrow=n,ncol=r),solve(Sigma),matrix(0,nrow=(p+1),ncol=r)),
                          rbind(matrix(0,nrow=(n+r),ncol=p),s_tmp,matrix(0,nrow=1,ncol=p)),
                          rbind(matrix(0,nrow=(n+r+p),ncol=1),eta2^-1))
      W=t(tilde_X)%*%s_tilde_Sigma%*%tilde_X
      M=t(tilde_X)%*%s_tilde_Sigma%*%tilde_Y
      inv_W=solve(W,tol=1e-90)
      F_i=MASS::mvrnorm(1,inv_W%*%M,inv_W)
      midF[,i-3]<-F_i
      tmpF=F_i
      
    }else{
      y0<-y[i/3]-alpha*y[i/3-1]-beta0-t(beta2)%*%FT[,i-1]-t(beta3)%*%FT[,i-2]
      tilde_Y=c(x[,i]-mu_x,FT[,i+1],tmpF,y0)
      tilde_X=rbind(Theta,A,tilde_A,t(beta1))
      s_tilde_Sigma=cbind(rbind(solve(Omega),matrix(0,nrow=(r+p+1),ncol=n)),
                          rbind(matrix(0,nrow=n,ncol=r),solve(Sigma),matrix(0,nrow=(p+1),ncol=r)),
                          rbind(matrix(0,nrow=(n+r),ncol=p),s_tmp,matrix(0,nrow=1,ncol=p)),
                          rbind(matrix(0,nrow=(n+r+p),ncol=1),eta2^-1))
      W=t(tilde_X)%*%s_tilde_Sigma%*%tilde_X
      M=t(tilde_X)%*%s_tilde_Sigma%*%tilde_Y
      inv_W=solve(W,tol=1e-90)
      F_i=MASS::mvrnorm(1,inv_W%*%M,inv_W)
      midF[,i-3]<-F_i
      tmpF=F_i
    }
  }
  
  ifnot1<-as.numeric(rel!=1)
  
  midF2=matrix(0,nrow=r,ncol=(t-3*k+ifnot1))
  for (i in (3*k-2):(t-3+ifnot1)){
    tilde_Y_ex=c(x[,i]-mu_x,FT[,i+1],tmpF)
    tilde_X_ex=rbind(Theta,A,tilde_A)
    
    s_tilde_Sigma_ex=cbind(rbind(solve(Omega),matrix(0,nrow=(r+p),ncol=n)),
                           rbind(matrix(0,nrow=n,ncol=r),solve(Sigma),matrix(0,nrow=p,ncol=r)),
                           rbind(matrix(0,nrow=(n+r),ncol=p),s_tmp))
    
    W_ex=t(tilde_X_ex)%*%s_tilde_Sigma_ex%*%tilde_X_ex
    M_ex=t(tilde_X_ex)%*%s_tilde_Sigma_ex%*%tilde_Y_ex
    inv_Wex=solve(W_ex,tol=1e-90)
    F.ex=MASS::mvrnorm(1,inv_Wex%*%M_ex,inv_Wex)
    midF2[,(i-3*k+3)]<-F.ex
    tmpF=F.ex
  }
  
  ##When i=t-2,t-1, t##
  
  if (rel==1){
    
    nj2<-dim(I2)[1]
    tilde_Y_t2=c(I2%*%(x[,t-2]-mu_x),tmpF)
    tilde_X_t2=rbind(I2%*%Theta,tilde_A)
    Omega_T2<-I2%*%Omega%*%t(I2)
    s_tilde_Sigma_t2=cbind(rbind(solve(Omega_T2),matrix(0,nrow=p,ncol=nj2)),
                           rbind(matrix(0,nrow=nj2,ncol=p),s_tmp))
    Wt2=t(tilde_X_t2)%*%s_tilde_Sigma_t2%*%tilde_X_t2
    Mt2=t(tilde_X_t2)%*%s_tilde_Sigma_t2%*%tilde_Y_t2
    
    inv_Wt2=solve(Wt2,tol=1e-90)
    Ft2=MASS::mvrnorm(1,inv_Wt2%*%Mt2,inv_Wt2)
    tmpF=Ft2
    
    nj1<-dim(I1)[1]
    tilde_Y_t1=c(I1%*%(x[,t-1]-mu_x),tmpF)
    tilde_X_t1=rbind(I1%*%Theta,tilde_A)
    Omega_T1<-I1%*%Omega%*%t(I1)
    s_tilde_Sigma_t1=cbind(rbind(solve(Omega_T1),matrix(0,nrow=p,ncol=nj1)),
                           rbind(matrix(0,nrow=nj1,ncol=p),s_tmp))
    Wt1=t(tilde_X_t1)%*%s_tilde_Sigma_t1%*%tilde_X_t1
    Mt1=t(tilde_X_t1)%*%s_tilde_Sigma_t1%*%tilde_Y_t1
    
    inv_Wt1=solve(Wt1,tol=1e-90)
    Ft1=MASS::mvrnorm(1,inv_Wt1%*%Mt1,inv_Wt1)
    
    res=cbind(startF1,startF,midF,midF2,Ft2,Ft1)
  }else{ #rel==2/3
    
    nj1<-dim(I1)[1]
    tilde_Y_t1=c(I1%*%(x[,t-1]-mu_x),tmpF)
    tilde_X_t1=rbind(I1%*%Theta,tilde_A)
    Omega_T1<-I1%*%Omega%*%t(I1)
    s_tilde_Sigma_t1=cbind(rbind(solve(Omega_T1),matrix(0,nrow=p,ncol=nj1)),
                           rbind(matrix(0,nrow=nj1,ncol=p),s_tmp))
    Wt1=t(tilde_X_t1)%*%s_tilde_Sigma_t1%*%tilde_X_t1
    Mt1=t(tilde_X_t1)%*%s_tilde_Sigma_t1%*%tilde_Y_t1
    
    inv_Wt1=solve(Wt1,tol=1e-90)
    Ft1=MASS::mvrnorm(1,inv_Wt1%*%Mt1,inv_Wt1)
    tmpF=Ft1
    
    
    nj<-dim(I)[1]
    tilde_Y_t=c(I%*%(x[,t]-mu_x),tmpF)
    tilde_X_t=rbind(I%*%Theta,tilde_A)
    Omega_T<-I%*%Omega%*%t(I)
    s_tilde_Sigma_t=cbind(rbind(solve(Omega_T),matrix(0,nrow=p,ncol=nj)),
                           rbind(matrix(0,nrow=nj,ncol=p),s_tmp))
    Wt=t(tilde_X_t)%*%s_tilde_Sigma_t%*%tilde_X_t
    Mt=t(tilde_X_t)%*%s_tilde_Sigma_t%*%tilde_Y_t
    
    inv_Wt=solve(Wt,tol=1e-90)
    Ft=MASS::mvrnorm(1,inv_Wt%*%Mt,inv_Wt)
    res=cbind(startF1,startF,midF,midF2,Ft1,Ft)
  }
  
  return(res)
}


######mu########
mu_update_lag<-function(omega,rel,data,factor,theta,I,I1,I2){
  The<-as.matrix(theta)
  FT<-as.matrix(factor)
  Omg<-as.matrix(omega)
  X<-as.matrix(data)
  n<-dim(Omg)[1]
  r<-dim(FT)[1]
  
  t_null=dim(X)[2]
  
  if (rel==1){t=t_null+1}else{t=t_null}
 
  ifnot1<-as.numeric(rel!=1)
  
  tmp=matrix(0,nrow=n,ncol=1)
  for (i in 1:(t-3+ifnot1)){
    K=The%*%FT[,i]
    tmp<-tmp+solve(Omg)%*%(X[,i]-K)
  }
  
  if (rel==1){
    inv_sig<-(t-3)*solve(Omg)+diag(n)+t(I2)%*%solve(I2%*%Omg%*%t(I2))%*%I2+
      t(I1)%*%solve(I1%*%Omg%*%t(I1))%*%I1
    
    mean_mu<-solve(inv_sig)%*%(tmp+t(I2)%*%solve(I2%*%Omg%*%t(I2))%*%I2%*%
                                 (X[,t-2]-The%*%FT[,t-2])+
                                 t(I1)%*%solve(I1%*%Omg%*%t(I1))%*%I1%*%
                                 (X[,t-1]-The%*%FT[,t-1]))
    
  }else{#rel==2/3
    inv_sig<-(t-3+ifnot1)*solve(Omg)+diag(n)+t(I)%*%solve(I%*%Omg%*%t(I))%*%I+
      t(I1)%*%solve(I1%*%Omg%*%t(I1))%*%I1
    
    mean_mu<-solve(inv_sig)%*%(tmp+t(I)%*%solve(I%*%Omg%*%t(I))%*%I%*%
                                 (X[,t]-The%*%FT[,t])+
                                 t(I1)%*%solve(I1%*%Omg%*%t(I1))%*%I1%*%
                                 (X[,t-1]-The%*%FT[,t-1]))
  }
  new_mu<-MASS::mvrnorm(1,mu=mean_mu,Sigma=solve(inv_sig))
  
  return(new_mu)
}




######metroplis-Hastings#####

###tau###
# alpha_tau=2.01 
# beta_tau=0.1

##final for with constraint2/ no S matrix, no tau on y
# tau_update<-function(theta,lambda,tau){
#   
#   The<-as.matrix(theta)
#   r<-dim(The)[2]
#   n<-dim(The)[1]
#   
#   tmp<-0
#   for (j in 1:r){
#     tmp=tmp+t(The[j:n,j])%*%The[j:n,j]/(lambda[j]^2)
#   }
#   
#   shape=0.5*(n*r)+alpha_tau
#   scale=0.5*tmp+beta_tau
#   
#   
#   #tau2<-MCMCpack::rinvgamma(1,shape=shape,scale=scale)
#   
#   time=0
#   repeat {
#     tau2<-MCMCpack::rinvgamma(1,shape=shape,scale=scale)
#     time=time+1
#     if ((tau2<=1)|(time>3000)){break}
#   }
#   if (time<=3000){
#     res=sqrt(tau2)
#   }else{
#     res=tau
#   }
#   return(res)
# }

#final for no constraint2/ S matrix, no tau on y-edit
# tau_update<-function(theta,lambda){
#   
#   The<-as.matrix(theta)
#   r<-dim(The)[2]
#   n<-dim(The)[1]
#   
#   tmp<-0
#   for (j in 1:r){
#     tmp=tmp+t(The[,j])%*%The[,j]/(lambda[j]^2)
#   }
#   
#   shape=0.5*(n*r)+alpha_tau
#   scale=0.5*tmp+beta_tau
#   
#   
#   repeat {
#     tau2<-MCMCpack::rinvgamma(1,shape=shape,scale=scale)
#     if (tau2<=1){break}
#   }
#   
#   return(sqrt(tau2))
# }


##final for no constraint2/ no S matrix, no tau on y
# tau_update<-function(theta,lambda){
#   
#   The<-as.matrix(theta)
#   r<-dim(The)[2]
#   n<-dim(The)[1]
#   
#   tmp<-0
#   for (j in 1:r){
#     tmp=tmp+t(The[,j])%*%The[,j]/(lambda[j]^2)
#   }
#   
#   shape=0.5*(n*r)+alpha_tau
#   scale=0.5*tmp+beta_tau
#   
#   
#   repeat {
#   tau2<-MCMCpack::rinvgamma(1,shape=shape,scale=scale)
#   if (tau2<=1){break}
#   }
#   
#   return(sqrt(tau2))
# }


# tau_update<-function(data,GDP,mu_x,factor,Theta,lambda,Omega,beta,I,e2){
#   x=as.matrix(data)
#   n=dim(x)[1]
#   t=dim(x)[2]
#   
#   y=as.matrix(GDP)
#   k<-length(y)
#   
#   FT=as.matrix(factor)
#   r=dim(FT)[1]
#   
#   The<-as.matrix(Theta)
#   mu<-as.vector(mu_x)
#   Lam<-diag(lambda)
#   omg<-as.matrix(Omega)
#   beta<-as.vector(beta)
#   
#   beta0=beta[1]
#   beta1=beta[2:(r+1)]
#   beta2=beta[(r+2):(2*r+1)]
#   beta3=beta[(2*r+2):(3*r+1)]
#   beta4=beta[(3*r+2)]
#   
#   
#   tmp1=0
#   tmp2=0
#   for (i in 1:(t-1)){
#     tmp1=tmp1+t(x[,i]-mu)%*%solve(omg)%*%The%*%Lam%*%FT[,i]
#     tmp2=tmp2+t(The%*%Lam%*%FT[,i])%*%solve(omg)%*%The%*%Lam%*%FT[,i]
#   }
#   
#   
#   tmp3=0
#   tmp4=0
#   for (j in 2:k){
#     tmp3=tmp3+(y[j]-beta0-beta4*y[j-1])*(t(beta1)%*%Lam%*%FT[,3*j]+
#                                            t(beta2)%*%Lam%*%FT[,3*j-1]+
#                                            t(beta3)%*%Lam%*%FT[,3*j-2])
#     tmp4=tmp4+(t(beta1)%*%Lam%*%FT[,3*j]+t(beta2)%*%Lam%*%FT[,3*j-1]+
#                  t(beta3)%*%Lam%*%FT[,3*j-2])^2
#   }
#   
#   inv_sig_tau=tmp2+t(I%*%The%*%Lam%*%FT[,t])%*%solve(I%*%omg%*%t(I))%*%
#     I%*%The%*%Lam%*%FT[,t]+(1/e2)*tmp4+1
#   
#   mu_tau=(1/inv_sig_tau)*(tmp1+t(I%*%(x[,t]-mu))%*%
#                             solve(I%*%omg%*%t(I))%*%I%*%The%*%Lam%*%FT[,t]+
#                             (1/e2)*tmp3)
#   
#   repeat {
#     tau_new<-rnorm(1,mean=mu_tau,sd=sqrt(1/inv_sig_tau))
#     if (tau_new>=0){break}
#   }
#   
#   
#   return(tau_new)
# }


# #alpha_tau=2.01, beta_tau=0.1
# tau_update<-function(theta,lambda){
# 
#   The<-as.matrix(theta)
#   r<-dim(The)[2]
#   n<-dim(The)[1]
# 
#   tmp<-0
#    for (j in 1:r){
#     tmp=tmp+t(The[,j])%*%The[,j]/(lambda[j]^2)
#    }
#   
#   shape=0.5
#   scale=0.5*tmp
# 
#   tau2_old=tau^2
#   tau2_new<-MCMCpack::rinvgamma(1,shape=shape,scale=scale)
# 
#   p_old=tau2_old/(1+tau2_old)
#   p_new=tau2_new/(1+tau2_new)
# 
#   p_all<-p_new/p_old
#   if (runif(1)<=p_all){
#     tau2_new<-tau2_new
#   }else{
#     tau2_new=tau2_old
#   }
#   return(sqrt(tau2_new))
# }


####lambda####
###no sum=r constraint##
lambda_update<-function(lambda,tau,beta,nu_lambda){
  #lambda/sd_lam R*1 vector,
  
  r=length(lambda)
  
  lambda_raw<-lambda
  
  # beta0=beta[1]
  # beta1=beta[2:(r+1)]
  # beta2=beta[(r+2):(2*r+1)]
  # beta3=beta[(2*r+2):(3*r+1)]
  # beta4=beta[(3*r+2)]
  
  
  for (j in 1:r){
    
    lambdaj_old<-lambda[j]
    lambda_old<-lambda
    
    tmp=0
    for (i in 1:3){
      beta_i<-beta[((i-1)*r+2):(i*r+1)]
      tmp=tmp+(beta_i[j])^2
    }
    
    scale_beta=2*tau^2/tmp
    
    lambdaj_new_raw=rgamma(1, shape=2.5, scale = scale_beta)
    lambdaj_new=sqrt(1/lambdaj_new_raw)
    
    lambda_new_all<-lambda
    lambda_new_all[j]<-lambdaj_new
    
    cr=as.numeric(min(abs(lambda_new_all))>=0.00000001)
    
    
    p_old=(1/(1+(lambdaj_old/nu_lambda^(j))^2))
    p_new=(1/(1+(lambdaj_new/nu_lambda^(j))^2))
    
    
    p_all<-p_new/p_old
    if ((runif(1)<=p_all)&(cr==1)){
      lambda<-lambda_new_all
    }else{
      lambda<-lambda_old
    }
  }
  
  
  if (rmse(lambda,lambda_raw)<=0.001){
    jump_lam=0
  }else{
    jump_lam=1
  }
  
  return(list("lambda"=lambda,"jump_lam"=jump_lam))
}


