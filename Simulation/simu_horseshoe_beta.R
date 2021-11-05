
set.seed(floor(rnorm(1,200,300)))
setwd("/Users/LJY/Dropbox/PHD/horeshoe/simulations")

n=40
t=302
k=100

#alpha_s<-3
#beta_s=2

alpha_s<-2
beta_s=6
r=6
mu_x<-10


simulate<-function(n,t,k,r,alpha_s,beta_s,mu_x,omg_type){
  
  # generate Thera & Omega
  n1<-ceiling(n/3)
  n3<-n-2*n1
  Omg1<- MCMCpack::riwish(v=n1+0.001,S=1/n*diag(n1))
  Omg2<- MCMCpack::riwish(v=n1+0.001,S=1/n*diag(n1))
  Omg3<- MCMCpack::riwish(v=n3+0.001,S=1/n*diag(n3))
  cor1<-cor(Omg1)
  cor2<-cor(Omg2)
  cor3<-cor(Omg3)
  if (omg_type==1){
    Omg<- MCMCpack::riwish(v=n,S=1/n*diag(n))
  }else if(omg_type==2){
    Omg<-adiag(Omg1,Omg2,Omg3)
  }else if(omg_type==3){
    cor_all<-adiag(cor1,cor2,cor3)
    cor_all[cor_all==0]<-0.1
    cor_true<-as.matrix(nearPD(cor_all,corr=T,keepDiag = T)$mat)
    sd<-c(diag(Omg1),diag(Omg2),diag(Omg3))
    Omg<-cor2cov(cor.mat=cor_true, sd=sqrt(sd))
  }else{
    cor_all<-adiag(cor1,cor2,cor3)
    cor_all[cor_all==0]<-0.9
    cor_true<-as.matrix(nearPD(cor_all,corr=T,keepDiag = T)$mat)
    sd<-c(diag(Omg1),diag(Omg2),diag(Omg3))
    Omg<-cor2cov(cor.mat=cor_true, sd=sqrt(sd))
  }
  
  Omg_cor<-cor(Omg)
  
  tau<-1
  lambda_raw<-c(4,4,0.1,0.1,0.1,0.1)
  #lambda=lambda_raw*sqrt(r/sum(lambda_raw^2))
  lambda<-lambda_raw
  

  
  S=tau*diag(lambda) 
  
  
 the<-MASS::mvrnorm(1,mu=rep(0,(n*r)),Sigma = diag(n*r))
 The<-matrix(the, ncol=r,nrow=n,byrow = T)
  
  # The<-matrix(0,ncol=r,nrow=n)
  # 
  # for (j in 1:r){
  #   repeat{
  #     The[j:n,j]=MASS::mvrnorm(1,mu=rep(0,(n+1-j)),Sigma = diag(n+1-j))
  #     if(The[j,j]>0){
  #       break
  #     }
  #   }
  # }
  
  
  #------------------------
  
  # generate A & Sigma
  
  Sig<-diag(c(5.5,3,1,0.5,0.25,0.1))
  #Sig<-diag(c(7.5,5.8,5.4,3.2))
 
  
  # Sigma_n<- MCMCpack::rinvgamma(r,shape=alpha_s,scale=beta_s)
  # Sig<-diag(sort(Sigma_n,decreasing = T))
  # 
  # a<-rep(0,r)
  # repeat{
  #   a[1]<-rnorm(1)
  #   if(abs(a[1])<1){
  #     break
  #   }
  # }
  # for (j in 2:r){
  #   repeat{
  #     a[j]<-rnorm(1)
  #     if((abs(a[j])<1)&((Sig[j,j]/(1-a[j]^2))<(Sig[j-1,j-1]/(1-a[j-1]^2)))){
  #       break
  #     }
  #   }
  # }
  
  
  a<-c(0.9,-0.8,0.75,0.7,-0.65,0.6)
  #a<-c(0.1,-0.3,0.3,-0.06)
  A<-diag(a)
  
  #------------------------
  
  eta2<-1
  
  beta_0<- 0.5
  beta_null<-c(5,-5,mvrnorm(1,mu=rep(0,r-2),Sigma = diag(rep(0.01,r-2))))
  #beta_null<-c(1,-1,0,0,0,0)
  beta_1<-beta_null
  beta_2<-beta_null
  beta_3<-beta_null
  
  # beta_null<-matrix(mvrnorm(1,mu=rep(0,r),Sigma = diag(r)),ncol=1)
  # beta_null<-c(2,1.4,1,0.5,0.35,0.175)
  # #beta_null<-c(2,1.4,1,0.5)
  # beta_1<-S%*%beta_null
  # beta_2<-S%*%beta_null
  # beta_3<-S%*%beta_null
  alpha<-0.15
  beta<-c(beta_0,beta_1,beta_2,beta_3,alpha)
  #------------------------------
  
  # Simulate F
  F0<- MASS::mvrnorm(1,mu=rep(0,r),
                     Sigma=solve(t(The)%*%The)%*%t(The)%*%Omg%*%The%*%solve(t(The)%*%The))
  
  F0
  
  FT<- cbind(F0,matrix(0,nrow=r,ncol=t))
  for (i in 1:t){
    FT[,i+1]<- A%*%FT[,i]+MASS::mvrnorm(1,mu=rep(0,r),Sigma=Sig) ##?#
  }
  FT<- FT[,-1]
  
  
  #---------------------------
  
  #simulate x_t 
  X<- matrix(0,nrow=n,ncol=t)
  for (i in 1:t){
    X[,i]=mu_x+The%*%FT[,i]+MASS::mvrnorm(1,mu=rep(0,n),Sigma=Omg)   #?##
  }
  
  
  
  #simulate y_k
  y<- NULL
  y[1]<-beta_0+t(beta_1)%*%FT[,3]+t(beta_2)%*%FT[,2]+
    t(beta_3)%*%FT[,1]+rnorm(1,0,eta2)
  for (i in 2:k){
    y[i]=beta_0+t(beta_1)%*%FT[,3*i]+t(beta_2)%*%FT[,(3*i-1)]+
      t(beta_3)%*%FT[,(3*i-2)]+alpha*y[i-1]+rnorm(1,0,eta2)
  }
  
  simuR <- list("Sig"=Sig,"A" = A,"Omg"=Omg,"The"=The,"tau"=tau,"lambda"=lambda,
                "X"=X,"y"=y,"FT"=FT,"Omg_cor"=Omg_cor,"beta"=beta,"eta2"=eta2)
  return(simuR)
}
#---------------------------
simu_ran<-simulate(n,t,k,r,alpha_s,beta_s,mu_x,1)
#simu_0cor<-simulate(n,t,k,r,alpha_s,beta_s,mu_x,2)

plot(simu_ran$X[1,1:100],type="l",xlab="Time",ylab="X1") #check how x_t looks
plot(ts(simu_ran$y),ylab="GDP") #check how y_k looks

plot(ts(simu_ran$FT[1,1:150]),ylab="First Factor")
#plot(ts(simu_ran$FT[6,1:150]),ylab="6th Factor")

saveRDS(simu_ran,file="simu_ran40_cons_beta2_6_s.rds")
#saveRDS(simu_0cor,file="simu_0cor_60.rds")
#saveRDS(simu3,file="simu3.rds")
#saveRDS(simu4,file="simu4.rds")

#saveRDS(simu_ran,file="simu_ran.rds")

####200 realizations##
X_all<-matrix(0,nrow=60*200,ncol=302)
Y<-matrix(0,nrow=100,ncol=200)
A_all<-matrix(0,nrow=3*200,ncol=3)
FT_all<-matrix(0,nrow=3*200,ncol=302)
Theta_all<-matrix(0,nrow=60,ncol=3*200)
for (i in 1:200){
  simu<-simulate(n=60,t=302,k=100,r=3,alpha_s,beta_s,mu_theta)
  X_all[(1+(i-1)*60):(60+(i-1)*60),]<-simu$X
  Y[,i]<-simu$y
  A_all[(3*i-2):(3*i),]<-simu$A
  FT_all[(3*i-2):(3*i),]<-simu$FT
  Theta_all[,(3*i-2):(3*i)]<-simu$The
}

MASS::write.matrix(X_all,file="X_all.txt",sep=" ")
MASS::write.matrix(Y,file="Y_all.txt",sep=" ")
MASS::write.matrix(A_all,file="A_all.txt",sep=" ")
MASS::write.matrix(FT_all,file="FT_all.txt",sep=" ")
MASS::write.matrix(Theta_all,file="Theta_all.txt",sep=" ")





MASS::write.matrix(simu$Omg,file="Omega.txt",sep=" ")
MASS::write.matrix(simu$The,file="Theta.txt",sep=" ")
MASS::write.matrix(simu$Omg_cor,file="Omg_cor.txt",sep=" ")


MASS::write.matrix(simu$Sig,file="Sigma.txt",sep=" ")
MASS::write.matrix(simu$A,file="A.txt",sep=" ")





MASS::write.matrix(simu$FT, file="factors.txt",sep=" ")
MASS::write.matrix(simu$X, file="var_s.txt",sep=" ")
MASS::write.matrix(simu$y, file="sim_gdp.txt",sep=" ")


The<-read.table("Theta.txt")
Omg<-read.table("Omega.txt")
Omg_cor<-read.table("Omg_cor.txt")

A<-read.table("A.txt")
Sig<-read.table("Sigma.txt")

FT<-read.table("factors.txt")
