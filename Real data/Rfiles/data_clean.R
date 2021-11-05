setwd("/Users/LJY/Dropbox/PHD/horeshoe/Real_data/Nowcasting-master/data/US")

require(tidyverse)
require(magittr)
require(readxl)
require(lubridate)
library(zoo)

path <- "/Users/LJY/Dropbox/PHD/horeshoe/Real_data/Nowcasting-master/data/US"

url_xlsx <- list.files(path, pattern = "*.xls", recursive = TRUE)

read_xlsx_files <- function(x){
  df <- read_xls(path = paste(path, x, sep = "/"),sheet=1)
  return(df)
}

df <- lapply(url_xlsx, read_xlsx_files ) %>%
  bind_rows()


#####
df1<- read_xlsx_files(url_xlsx[length(url_xlsx)])
df2<- read_xlsx_files(url_xlsx[length(url_xlsx)-5])

df1 %>% filter(Date>="1992-06-01") %>% data.frame()->df1_cut
#df2 %>% filter(Date>="1992-06-01") %>% data.frame()->df2_cut


#####focus on df1#

max(unique(df1_cut$Date)) ##last month 2017-01-01

df1_cut %>% 
  select_if(function(x) any(is.na(x))) %>% 
  summarise_each(funs(sum(is.na(.))))

df1_cut %>% filter(is.na(GACDISA066MSFRBNY)) %>% view()

df1_cut %>% slice(1:100) %>% view()


####get response####

df1_cut %>% dplyr::select(Date,GDPC1) %>% drop_na() %>% 
  mutate(quarter=quarter(Date, with_year = T),
         index=seq(1,length(quarter)))->response_raw
  

response_raw %>% mutate(index_lag=c(seq(2,length(index)+1))) %>% 
  dplyr::select(GDPC1,index_lag) %>% 
  dplyr::rename(index=index_lag,GDP_lag=GDPC1)->response_lag


response_raw %>% left_join(.,response_lag,by="index") %>% 
  drop_na() %>% 
  mutate(gdp=((1+(GDPC1-GDP_lag)/GDP_lag)^4-1)*100) %>% 
  dplyr::select(quarter,gdp) %>% filter(quarter>=1992.4) %>% 
  slice(-1) %>% 
  mutate(quarter_index=seq(1,length(quarter))) ->response_final


####quarter 1992.4-gdp=4.06783540

##save response##
setwd("/Users/LJY/Dropbox/PHD/horeshoe/Real_data/clean_data")
saveRDS(response_final,"GDP_old.rds")

#response_final<-readRDS("GDP_old.rds")

######transform monthly series#########
colSums(apply(df1_cut,2,is.na))

df1_cut %>% filter(Date>="1992-11-01") %>% select(TTLCONS)->TTLCONS_sep

df1_cut %>% filter(Date>="1992-11-01")%>% 
  select(-JTSJOL,-PPIFIS,-GACDISA066MSFRBNY,-TTLCONS,
         -PCEC96, -GDPC1,-ULCNFB,-A261RX1Q020SBEA)%>% 
  mutate(TTLCONS=c(458080,TTLCONS_sep$TTLCONS[-1])) %>% 
  do(na.locf(.))->df_cut




#################transform groups#################
##1: no transform
df_cut %>% select(Date, GACDFSA066MSFRBPHI)->trans1

##2: level change
df_cut %>% select(Date, PAYEMS, UNRATE, TCU,PERMIT)->trans2  

##3: MOM change
df_cut %>% select(Date,TTLCONS,BOPTEXP,BOPTIMP,DGORDER,BUSINV,IR,IQ,RSAFS,
                  INDPRO,WHLSLRIMSA,CPIAUCSL,CPILFESL,HOUST,HSN1F,DSPIC96,
                  PCEPI,PCEPILFE)->trans3

##functions for transformation#

trans2_fun<-function(vector){
  vector %>% data.frame() %>% 
    mutate(value_lag=lag(.)) %>% 
    mutate(transformed_value=.-value_lag)->res_raw
  return(as.vector(res_raw$transformed_value))
}

trans3_fun<-function(vector){
  vector %>% data.frame() %>%
    mutate(value_lag=lag(.)) %>% 
    mutate(transformed_value=100*(.-value_lag)/value_lag)->res_raw
  return(as.vector(res_raw$transformed_value))
}

##
trans2 %>% mutate_at(vars(-Date),trans2_fun) %>% 
  filter(Date>="1992-12-01")->trans2_final


trans3 %>% mutate_at(vars(-Date),trans3_fun) %>% 
  filter(Date>="1992-12-01")->trans3_final

##combine transformed series##

library(plyr)
 
trans_final<-join_all(list(trans1,trans2_final,trans3_final),
                      by='Date', type='left')


trans_final %>% filter(Date>="1992-12-01")->X_final

X_final %>% select(GACDFSA066MSFRBPHI,PAYEMS,UNRATE,IR,IQ,RSAFS,INDPRO,
                   TCU,CPIAUCSL,CPILFESL,HOUST,PERMIT,DGORDER,WHLSLRIMSA,
                   HSN1F,DSPIC96,PCEPI,PCEPILFE,TTLCONS,BOPTEXP,BOPTIMP,
                   BUSINV) %>% as.matrix()->X_trans




###Save transformed series
setwd("/Users/LJY/Dropbox/PHD/horeshoe/Real_data/clean_data")
saveRDS(t(X_trans),"X_old.rds")

X_old<-readRDS("X_old.rds")

###############get factors####################
x<-X_old
y<-response_final$gdp



G=1000
M=8000
alpha_s<-2
beta_s<-6
mu_theta<-0

n<-dim(x)[1]
x_old<-x[,1:288]
#x_old[31:40,97]<-1000000000000

y_old<-c(y[1:96])





X<- as.matrix(x_old[,-dim(x_old)[2]])
n<-dim(X)[1]
t<-dim(X)[2]
k<-length(y_old)
Z<- matrix(0,nrow=n,ncol=t)
for (i in 1:n){
  Z[i,]=(X[i,]-mean(X[i,]))/sd(X[i,])
}
S<- cor(t(Z))

Z_all<- matrix(0,nrow=n,ncol=t+1)
for (i in 1:n){
  Z_all[i,]=(x_old[i,]-mean(x_old[i,]))/sd(x_old[i,])
}


r=6
R<- rARPACK::eigs(S,r)
V<- R$vectors[,1:r]
D<- diag(R$values[1:r])
t(V)%*%V

tilde.f<- matrix(0,nrow=r,ncol=t)
for (i in 1:t){
  tilde.f[,i]<- t(V)%*%X[,i]
}

# tmp<- matrix(0,nrow=3,ncol=3)
# tmp_s<- matrix(0,nrow=n,ncol=3)
# for (i in 1:t){
#   tmp<- tmp+tilde.f[,i]%*%t(tilde.f[,i])
#   tmp_s<- tmp_s+X[,i]%*%t(tilde.f[,i])
# }

theta_initial<-V
#theta_initial<-tmp_s%*%solve(tmp)


omega_initial<- diag(diag(S-V%*%D%*%t(V)))
omega_kf<-omega_initial
# if (rel==2){
#   diag(omega_kf)[31:40]<-10^10
# }else{
#   diag(omega_kf)[16:40]<-10^10
# }


tmp1<- matrix(0,nrow=r,ncol=r)
tmp2<- matrix(0,nrow=r,ncol=r)
tmp3<- matrix(0,nrow=r,ncol=r)
for (i in 2:t){
  tmp1<- tmp1+tilde.f[,i-1]%*%t(tilde.f[,i-1])
  tmp2<- tmp2+tilde.f[,i]%*%t(tilde.f[,i-1])
  tmp3<- tmp3+tilde.f[,i]%*%t(tilde.f[,i])
}

A_initial<-tmp2%*%solve(tmp1)

sigma_initial<- tmp3/(t-1)-A_initial%*%(tmp1/(t-1))%*%t(A_initial)

#m0<-tilde.f[,t]
m0<- matrix(rep(0,r),ncol=1)
C0_fun<- function(Th,O){
  solve(t(Th)%*%Th)%*%t(Th)%*%O%*%Th%*%solve(t(Th)%*%Th)
}
C0<- C0_fun(theta_initial,omega_kf)
kal_fil_fun<- function(m,c,Th,O,A,S,D){
  t=dim(D)[2]
  r=length(m)
  factor<- matrix(0,nrow=r,ncol=t)
  for (i in 1:t){
    a<- A%*%m
    R<- A%*%c%*%t(A)+S
    e<- D[,i]-Th%*%a
    Q<- Th%*%R%*%t(Th)+O
    #Q<-nearPD(Th%*%R%*%t(Th)+O)$mat
    K<- R%*%t(Th)%*%solve(Q)
    m<- a+K%*%e
    c<- R-R%*%t(Th)%*%solve(Q)%*%Th%*%R
    factor[,i]<- as.vector(m)
  }
  return(factor)
}

factors_initial<- kal_fil_fun(m=m0,c=C0,Th=theta_initial,O=omega_kf,A=A_initial,S=sigma_initial,D=Z_all)
t<-t+1

factor_old<-factors_initial
#factor_old<-simu$FT[,1:t]



saveRDS(factors_initial, "real_old_ft_initial.rds")


