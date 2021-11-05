
library(ggforce)
library(stringr)
library(readr)
library(tidyverse)
library(ggpubr)

setwd("/Users/LJY/Dropbox/PHD/horeshoe/simulations")
simu<-read_rds("simu_ran40_cons_beta2_6_s.rds")



setwd("/Users/LJY/Dropbox/PHD/horeshoe/results")
horse_40_r2<-readRDS("para_all40_r2_s.rds") 
#horse_40_r1_basic=horse_40_r1[[1]]
horse_40_r2_basic=horse_40_r2
FT_grs<-horse_40_r2$ft_grs

####plot beta##
####plot beta####
beta.res<-horse_40_r2_basic$beta_all
beta_res<-data.frame(t(abs(beta.res)))

beta_res %>% gather(.,type0,value,X1:X20,factor_key=TRUE) %>% 
  mutate(num=as.numeric(str_extract(type0, "[0-9]+"))) %>% 
  mutate(type=ifelse(num==1,"beta0",
                     ifelse(num==20,"alpha",
                            ifelse((num>=2)&(num<=7),"beta1",
                                   ifelse((num>=8)&(num<=13),"beta2","beta3")))),
         index=ifelse(!(num%in%c(1,20)),(num-1)%%6,7))->beta_res_all

#beta_res_all$index<-as.factor(beta_res_all$index)



vline_beta <- summarise(group_by(beta_res_all,index,type),
                   mean = mean(value),
                   q1=quantile(value,0.025),
                   q3=quantile(value,0.975))

vline_beta %>%filter(index%in%c(1,2))->vline1_beta
vline_beta %>% filter(index%in%c(3,4,5,0))->vline2_beta
vline_beta %>% filter(type%in%c("beta0"))->vline3_beta
vline_beta %>% filter(type%in%c("alpha"))->vline4_beta


# beta_res_all %>% filter(index%in%c(1,0)) %>% 
#   mutate(Index=as.factor(index+1)) %>% 
#   ggplot(., aes(value,..scaled..)) + 
#   geom_density(aes(group=Index,color=Index,fill=Index),alpha=0.4)+
#   geom_vline(data=vline1_beta, aes(xintercept = mean),
#              col='red',size=0.4,alpha=0.4)+
#   facet_grid(type~.)+ggtitle("Posterior Density of Important Beta/Simu2")
# 
# beta_res_all %>% filter(index%in%c(2,3,4,5)) %>%
#   mutate(Index=as.factor(index+1)) %>% 
#   ggplot(., aes(value,..scaled..)) + 
#   geom_density(aes(group=Index,color=Index,fill=Index),alpha=0.4)+
#   geom_vline(data=filter(vline2_beta,index==3), aes(xintercept = mean),
#              col='red',size=0.4,alpha=0.4)+
#   facet_grid(type~.)+ggtitle("Posterior Density of unimportant Beta/Simu2")

##density plot seperate
beta_true<-data.frame(index=c(seq(1,5),0),beta_true=simu$beta[2:7])

beta_res_all %>%left_join(.,beta_true,by="index") %>%
  filter(index%in%c(1,2)) %>% 
  ggplot(., aes(value,..scaled..)) + 
  geom_density(color="blue",fill="blue",alpha=0.4)+
  geom_vline(aes(xintercept = abs(beta_true)),
             col='red',size=0.4,alpha=0.4)+
  facet_grid(type~as.factor(index),scales = 'free_x')+
  ggtitle("Posterior Density of important Beta/Simu1")

beta_res_all %>%left_join(.,beta_true,by="index") %>%
  filter(index%in%c(3,4,5,0)) %>% 
  ggplot(., aes(value,..scaled..)) + 
  geom_density(color="blue",fill="blue",alpha=0.4)+
  facet_grid(type~as.factor(index),scales = 'free')+
  xlim(0,10^-7)+
  ggtitle("Posterior Density of unimportant Beta/Simu1")




b0_p<-beta_res_all %>% filter(type%in%c("beta0"))%>%
  ggplot(., aes(value,..scaled..)) + 
  geom_density(aes(color=type,fill=type),alpha=0.4)+
  geom_vline(aes(xintercept = 0.5),
             col='red',size=0.4,alpha=0.4)

alpha_p<-beta_res_all %>% filter(type%in%c("alpha"))%>%
  ggplot(., aes(value,..scaled..)) + 
  geom_density(aes(color=type,fill=type),alpha=0.4)+
  geom_vline(aes(xintercept =0.15),
             col='red',size=0.4,alpha=0.4)



b0a_den<-ggarrange(b0_p, alpha_p,ncol = 1, nrow = 2)
annotate_figure(b0a_den,
                top = text_grob("Posterior Density for Beta0 and Alpha/Simu2", 
                                face = "bold", size = 10))


####plot lambda###
lambda.res<-horse_40_r2_basic$lambda_all
lambda.res %>% data.frame() %>%
  rename(lam1="X1",lam2="X2",lam3="X3",lam4="X4",lam5="X5",lam6="X6") %>%
  gather(.,lambda,value,lam1:lam6,factor_key = T)->lambda_res

##posterior density plot##

vline <- summarise(group_by(lambda_res,lambda),
                   mean = mean(value),
                   q1=quantile(value,0.025),
                   q3=quantile(value,0.975))
vline %>%
  filter(lambda %in% c("lam1", "lam2"))->vline1

vline %>%
  filter(lambda %in% c("lam3", "lam4","lam5", "lam6"))->vline2


lambda_res %>% filter(lambda%in%c("lam1","lam2")) %>% 
  ggplot(., aes(value,..scaled..)) + 
  geom_density(aes(group=lambda,color=lambda,fill=lambda),alpha=0.4)+
  geom_vline(data=vline1, aes(xintercept = mean),
             col='red',size=0.4,alpha=0.4)->lam_den1

lambda_res %>% filter(!(lambda%in%c("lam1","lam2"))) %>% 
  ggplot(., aes(value,..scaled..)) + 
  geom_density(aes(group=lambda,color=lambda,fill=lambda),alpha=0.4)+
  geom_vline(data=filter(vline2,lambda=="lam3"), aes(xintercept = mean),
             col='red',size=0.4,alpha=0.4)->lam_den2

lam_den<-ggarrange(lam_den1, lam_den2,ncol = 1, nrow = 2)
annotate_figure(lam_den,
                top = text_grob("Posterior Density of Lambda for Simulation 2", 
                                face = "bold", size = 10))


# annotate_figure(lam_den,
#                 top = text_grob("Posterior Density for Simulation 1", face = "bold", size = 14),
#                 bottom = text_grob("Data source: \n ToothGrowth data set", color = "blue",
#                                    hjust = 1, x = 1, face = "italic", size = 10),
#                 left = text_grob("Figure arranged using ggpubr", color = "green", rot = 90),
#                 right = text_grob(bquote("Superscript: ("*kg~NH[3]~ha^-1~yr^-1*")"), rot = 90),
#                 fig.lab = "Figure 1", fig.lab.face = "bold"
# )


##histogram##
# lambda_res%>%
#   filter(lambda %in% c("lam1", "lam2"))%>%
#   ggplot(data=., aes(value)) +
#   geom_histogram(aes(y =..density..),col="red",
#                  fill="orange",
#                  binwidth=0.5,
#                  alpha=1) +
#   geom_density(col=4)+
#   geom_vline(data=vline1, aes(xintercept = mean),col='blue',size=0.4)+
#   geom_vline(data=vline1, aes(xintercept = q1),col='green',size=0.4)+
#   geom_vline(data=vline1, aes(xintercept = q3),col='green',size=0.4)+
#   facet_wrap(lambda~., ncol=2)
# 
# lambda_res%>%
#   filter(lambda %in% c("lam3", "lam4","lam5", "lam6"))%>%
#   ggplot(data=., aes(value)) +
#   geom_histogram(aes(y =..density..),col="red",
#                  fill="orange",
#                  binwidth=0.5*10^-8,
#                  alpha=1) +
#   geom_density(col=4)+facet_wrap(lambda~., ncol=2)+
#   geom_vline(data=vline2, aes(xintercept = mean),col='blue',size=0.4)+
#   geom_vline(data=vline2, aes(xintercept = q1),col='green',size=0.4)+
#   geom_vline(data=vline2, aes(xintercept = q3),col='green',size=0.4)


####################plot kappa=1/(1+lambda^2)##############
kappa_res<-lambda_res %>% 
  mutate(kappa=1/(1+value^2))

lambda_true<-cbind(c("lam1","lam2","lam3","lam4","lam5","lam6"),
                   c("beta1","beta2","beta3","beta4","beta5","beta6"),
                   simu$beta[2:7])

names(lambda_true)<-c("lambda","beta","beta_value")

kappa_res_all<-lambda_true %>% data.frame() %>% 
  rename(lambda="X1",beta="X2",beta_true="X3")%>% 
  mutate_at("beta_true",as.numeric) %>% 
  right_join(.,kappa_res,by="lambda") %>% 
  rename(lam_value="value") 




#boxplot
kappa_res_all %>% filter(lambda%in%c("lam1","lam2")) %>%
  ggplot(aes(x=beta, y=kappa,fill="red"),alpha=0.4) + 
  geom_boxplot()->kappa_p1



kappa_res_all %>%ggplot(aes(x=beta, y=kappa),alpha=0.4)+ 
  geom_boxplot()+ggtitle("Boxplots for Kappa/Simu2")

###scatter plot on posterior mean ??

find_mode<-function(var){ #var is a vector
  var<-as.numeric(var)
  density_estimate<-density(var)
  mode_value <- density_estimate$x[which.max(density_estimate$y)]
  return(mode_value)
}

kappa_res_all %>% group_by(beta,beta_true) %>% 
  summarise(kappa_mean=mean(kappa)) %>% 
  ggplot(aes(x=beta_true, y=kappa_mean))+
  geom_point()+ggtitle("Posterior Mean for Kappa/Simu2")

kappa_res_all %>% group_by(beta,beta_true) %>% 
  summarise(kappa_mode=find_mode(kappa)) %>% 
  ggplot(aes(x=beta_true, y=kappa_mode))+
  geom_point()+ggtitle("MAP for Kappa/Simu2")


##posterior density plot##


vline <- summarise(group_by(kappa_res_all,beta),
                   mean = mean(kappa),
                   q1=quantile(kappa,0.025),
                   q3=quantile(kappa,0.975),
                   mode=find_mode(kappa))


kappa_res_all %>% 
  ggplot(., aes(kappa,..scaled..)) + 
  geom_density(aes(group=beta,color=beta,fill=beta),alpha=0.4)+
  geom_vline(data=vline, aes(xintercept = mode),
             col='red',size=0.4,alpha=0.4)+
  ggtitle("Posterior Density of Kappa for Simulation 4")

####heatmap##
kappa_res_all %>% group_by(beta) %>%
  summarise(mode_kappa=find_mode(kappa),
            index=1) %>% 
ggplot(aes(x=index,y = beta)) +
  geom_tile(aes(fill = mode_kappa ),alpha=0.6) +
  geom_text(aes(label = round(mode_kappa, 2)))+
  scale_fill_continuous(low = "violetred", high = "aquamarine")+
  ggtitle("Mode for Kappa/Simulation1")

##########plot factors-first two####
factor_bay1<-horse_40_r2_basic$factor1_all
factor_bay2<-horse_40_r2_basic$factor2_all

factor1_bay<-apply(abs(factor_bay1),2,mean)
factor2_bay<-apply(abs(factor_bay2),2,mean)

t<-length(factor1_bay)

factor_simu<-data.frame(cbind(factor1_bay,factor2_bay,FT_grs[1,],FT_grs[2,],
                               simu$FT[1,1:t],simu$FT[2,1:t]))

# factor_simu %>% mutate(f2_bay_f=abs(V6)+0.2*abs(factor2_bay-V6),
#                        f1_bay_f=abs(V5)+0.2*abs(factor1_bay-V5)) %>%
#   select(-factor2_bay,-factor1_bay) %>%
#   rename(factor2_bay=f2_bay_f,
#          factor1_bay=f1_bay_f)%>%
#   rename(f1_bay=factor1_bay,f2_bay=factor2_bay,
#          f1_grs=V3,f2_grs=V4,f1_true=V5,f2_true=V6) %>%
#   mutate(month=seq(1,length(f1_bay)))%>%
#   gather(.,type,value,f1_grs:f1_bay,factor_key = T) %>%
#   mutate(factor=str_extract(type, ".+?(?=_)"),
#          method=str_sub(type, 4))->factor_p


factor_simu %>% rename(f1_bay=factor1_bay,f2_bay=factor2_bay,
                        f1_grs=V3,f2_grs=V4,f1_true=V5,f2_true=V6) %>% 
  mutate(month=seq(1,length(f1_bay)))%>% 
  gather(.,type,value,f1_bay:f2_true,factor_key = T) %>% 
  mutate(factor=str_extract(type, ".+?(?=_)"),
         method=str_sub(type, 4))->factor_p



factor_p %>% filter(method%in%c("bay","grs"))%>% 
  ggplot(aes(x=month,y=abs(value),color=method))+
  geom_line(size = 0.4)+facet_wrap(.~factor,scales = "free_y",ncol=1)+
  ggtitle("Factor comparison for simulation 2")


####check in-sample fit###
r=6

factor_bay1<-horse_40_r2_basic$factor1_all
factor_bay2<-horse_40_r2_basic$factor2_all
factor_bay3<-horse_40_r2_basic$factor3_all
factor_bay4<-horse_40_r2_basic$factor4_all
factor_bay5<-horse_40_r2_basic$factor5_all
factor_bay6<-horse_40_r2_basic$factor6_all
beta.res<-horse_40_r2_basic$beta_all

k<-floor(dim(factor_bay1)[2]/3)
  
y_insample_all<-matrix(0,nrow=1000,ncol=k)
for (i in 1:1000){
  factor_all_i<-rbind(factor_bay1[i,],factor_bay2[i,],factor_bay3[i,],
                    factor_bay4[i,],factor_bay5[i,],factor_bay6[i,])
  beta_i<-beta.res[,i]
  
  beta0<-beta_i[1]
  beta1<-beta_i[2:(r+1)]
  beta2<-beta_i[(r+2):(2*r+1)]
  beta3<-beta_i[(2*r+2):(3*r+1)]
  alpha<-beta_i[3*r+2]
  
  y_insample_all[i,1]<-beta0+t(beta1)%*%factor_all_i[,3]+
    t(beta2)%*%factor_all_i[,2]+t(beta3)%*%factor_all_i[,1]
  
  for (j in 2:k){
    y_insample_all[i,j]<-beta0+t(beta1)%*%factor_all_i[,3*j]+
      t(beta2)%*%factor_all_i[,3*j-1]+t(beta3)%*%factor_all_i[,3*j-2]+
      alpha*simu$y[j-1]
  }
}

y_insample_mean<-apply(y_insample_all,2,mean)
plot(y_insample_mean,simu$y[1:32])



######################visual 100 sample path###########


clean_100p<-function(index,mu_beta,sparse){ #sparse=1 means sparsity
  
  setwd("/Users/LJY/Dropbox/PHD/horeshoe/condo_insample_new")
  
  if (sparse==0){
    data_name<-paste("simu_",index,"_beta_",mu_beta,".rds",sep="")
    res_name<-paste("result_",index,"_beta_",mu_beta,".rds",sep="")
  
  }else{
    data_name<-paste("simu_",index,"_beta_",mu_beta,"_sparse.rds",sep="")
    res_name<-paste("result_",index,"_beta_",mu_beta,"_sparse.rds",sep="")
  }
  
  simu<-readRDS(data_name)
  result<-readRDS(res_name)
  
  r<-dim(simu$A)[1]
  n<-dim(simu$X)[1]
  t<-dim(result$ft_grs)[2]
  k<-dim(result$y_insample)[2]
  
  FT<-simu$FT[,1:t]
  
  ##beta###
  
  beta_all<-data.frame("beta_bay_abs"=apply(abs(result$beta_all),1,mean),
                       "beta_true"=simu$beta)
  
  beta_all<-beta_all %>% mutate(beta_num=c(0,rep(1,6),rep(2,6),rep(3,6),4),
                                beta_num_in=c(0,rep(seq(1,6),3),7),
                                index=index,mu_beta=mu_beta,sparse=sparse)
  
  ##
  
  kappa_test<-apply(1/(1+(result$lambda_all)^2),2,mean)
  kappa_mod<-apply(1/(1+(result$lambda_all)^2),2,find_mode)
  kappa_med<-apply(1/(1+(result$lambda_all)^2),2,median)
  
  kappa_all<-data.frame("lambda_bay"=lambda_test,
                        "lambda_mod"=lambda_mod,
                        "kappa_bay"=kappa_test,
                        "kappa_mod"=kappa_mod,
                        "kappa_med"=kappa_med,
                        "beta_bay"=beta_test[2:(r+1)],
                        "beta_true"=simu$beta[2:(r+1)])
  kappa_all<-kappa_all %>% mutate(index=index,mu_beta=mu_beta,
                                  lambda_num=seq(1,6),sparse=sparse)
  
  ###RISFEE for insample fit##
  y_insample<-apply(result$y_insample,2,mean)
  factors1_bay<-apply(abs(result$factor1_all),2,mean)
  factors2_bay<-apply(abs(result$factor2_all),2,mean)
  factors3_bay<-apply(abs(result$factor3_all),2,mean)
  factors4_bay<-apply(abs(result$factor4_all),2,mean)
  factors5_bay<-apply(abs(result$factor5_all),2,mean)
  factors6_bay<-apply(abs(result$factor6_all),2,mean)
  
  RISFEE<-data.frame("y_ris"=sum(abs(simu$y[1:k]-y_insample))/sum(abs(simu$y[1:k])),
                     "f1_bay_ris"=sum(abs(abs(simu$FT[1,1:t])-factors1_bay))/sum(abs(simu$FT[1,1:t])),
                     "f2_bay_ris"=sum(abs(abs(simu$FT[2,1:t])-factors1_bay))/sum(abs(simu$FT[2,1:t])),
                     "f3_bay_ris"=sum(abs(abs(simu$FT[3,1:t])-factors1_bay))/sum(abs(simu$FT[3,1:t])),
                     "f4_bay_ris"=sum(abs(abs(simu$FT[4,1:t])-factors1_bay))/sum(abs(simu$FT[4,1:t])),
                     "f5_bay_ris"=sum(abs(abs(simu$FT[5,1:t])-factors1_bay))/sum(abs(simu$FT[5,1:t])),
                     "f6_bay_ris"=sum(abs(abs(simu$FT[6,1:t])-factors1_bay))/sum(abs(simu$FT[6,1:t])))

  RISFEE<-RISFEE %>% mutate(index=index,mu_beta=mu_beta,sparse=sparse)
  
  
  ###factors/only first 2##
  
  factors_compare<-data.frame("Bay_1"=factors1_bay,"Bay2"=factors2_bay,
                              "grs_1"=result$ft_grs[1,],"grs_2"=result$ft_grs[2,],
                              "true_1"=simu$FT[1,1:t],"true_2"=simu$FT[2,1:t])
  
  factors_compare<-factors_compare %>% 
    mutate(month=seq(1,t),index=index,mu_beta=mu_beta,sparse=sparse)
  
  RISFEE_grs<-data.frame("f1_grs_ris"=sum(abs(abs(simu$FT[1,1:t])-abs(result$ft_grs[1,])))/sum(abs(simu$FT[1,1:t])),
                         "f2_grs_ris"=sum(abs(abs(simu$FT[2,1:t])-abs(result$ft_grs[2,])))/sum(abs(simu$FT[2,1:t])),
                         "f3_grs_ris"=sum(abs(abs(simu$FT[3,1:t])-abs(result$ft_grs[3,])))/sum(abs(simu$FT[3,1:t])),
                         "f4_grs_ris"=sum(abs(abs(simu$FT[4,1:t])-abs(result$ft_grs[4,])))/sum(abs(simu$FT[4,1:t])),
                         "f5_grs_ris"=sum(abs(abs(simu$FT[5,1:t])-abs(result$ft_grs[5,])))/sum(abs(simu$FT[5,1:t])),
                         "f6_grs_ris"=sum(abs(abs(simu$FT[6,1:t])-abs(result$ft_grs[6,])))/sum(abs(simu$FT[6,1:t])))
  
  RISFEE_grs<-RISFEE_grs %>% mutate(index=index,mu_beta=mu_beta,sparse=sparse)
  
  return(list("beta_res"=beta_all,
              "kappa_res"=kappa_all,
              "RISF"=RISFEE,
              "RISF_grs"=RISFEE_grs,
              "factor_compare"=factors_compare))
  }


testFunction <- function (index,mu_beta,sparse) {
  return(tryCatch(clean_100p(index,mu_beta,sparse), error=function(e) NULL))
}

df_p<-expand_grid("mu_beta"=c(5,1),"index"=seq(1,100),sparse=c(0,1))
p100_all<-mapply(testFunction,index=df_p$index,mu_beta=df_p$mu_beta,
                 sparse=df_p$sparse)

df_add<-expand_grid("mu_beta"=0.1,"index"=seq(1,96),sparse=c(0,1))
p100_add<-mapply(testFunction,index=df_add$index,mu_beta=df_add$mu_beta,
                 sparse=df_add$sparse)



setwd("/Users/LJY/Dropbox/PHD/horeshoe/sample100_res")
saveRDS(p100_all,"sample100_all.rds") ##5/1
saveRDS(p100_add,"sample100_add.rds") ##0.1

p100_all<-read_rds("sample100_all.rds")

#p100_all_f<-mapply(testFunction,index=df_p$index,mu_beta=df_p$mu_beta)


beta_p100_all<-NULL
kappa_p100_all<-NULL
RISFEE_p100_all<-NULL
RISF_grs_p100_all<-NULL
ft_compare_p100_all<-NULL

for (i in 1:length(p100_all)){
  p_sep<-p100_all[i][[1]]
  
  if (as.numeric(is.null(p_sep))==0){
  beta_p100_all<-rbind(beta_p100_all,p_sep$beta_res)
  kappa_p100_all<-rbind(kappa_p100_all,p_sep$kappa_res)
  RISFEE_p100_all<-rbind(RISFEE_p100_all,p_sep$RISF)
  RISF_grs_p100_all<-rbind(RISF_grs_p100_all,p_sep$RISF_grs)
  ft_compare_p100_all<-rbind(ft_compare_p100_all,p_sep$factor_compare)
  }
}


#################
####for beta####
beta_p100_add %>% filter(beta_num_in%in%c(1,2)) %>% 
  group_by(mu_beta,sparse,beta_num_in) %>% 
  summarise(mean_beta=mean(beta_bay_abs),
            prob_less=sum(beta_bay_abs<=0.001)/length(beta_bay_abs))


beta_p100_add<-beta_p100_add %>% mutate(n=seq(1,length(beta_bay_abs)))



beta_p100_add %>% filter(sparse==0) %>% 
  dplyr::rename(Bay="beta_bay_abs",True="beta_true") %>% 
  gather(.,type,value,Bay:True,factor_key=TRUE) %>% 
  filter(beta_num_in%in%c(1,2)) %>% 
  ggplot(., aes(value,..scaled..)) + 
  geom_density(color="blue",fill="blue",alpha=0.4)+
  facet_grid(type~as.factor(beta_num_in),scales = 'free_x')+
  ggtitle("Posterior Density of important Beta/No sparsity")

beta_p100_add %>% filter(sparse==1) %>% 
  dplyr::rename(Bay="beta_bay_abs",True="beta_true") %>% 
  gather(.,type,value,Bay:True,factor_key=TRUE) %>% 
  filter(beta_num_in%in%c(1,2)) %>% 
  ggplot(., aes(value,..scaled..)) + 
  geom_density(color="blue",fill="blue",alpha=0.4)+
  facet_grid(type~as.factor(beta_num_in),scales = 'free_x')+
  ggtitle("Posterior Density of important Beta/Sparsity")

beta_p100_all_cl %>% filter(sparse==0)%>% 
  dplyr::rename(Bay="beta_bay_abs",Bay_ori="beta_bay_ori",True="beta_true") %>% 
  gather(.,type,value,Bay:True,factor_key=TRUE) %>% 
  filter(!(beta_num_in%in%c(1,2,0,7)),type!="Bay_ori") %>% 
  ggplot(., aes(value,..scaled..)) + 
  geom_density(color="blue",fill="blue",alpha=0.4)+
  facet_grid(type~as.factor(beta_num_in),scales = 'free_x')+
  ggtitle("Posterior Density of unimportant Beta/No sparsity")

beta_p100_all_cl %>% filter(sparse==1)%>% 
  dplyr::rename(Bay="beta_bay_abs",Bay_ori="beta_bay_ori",True="beta_true") %>% 
  gather(.,type,value,Bay:True,factor_key=TRUE) %>% 
  filter(!(beta_num_in%in%c(1,2,0,7)), type!="Bay_ori") %>% 
  ggplot(., aes(value,..scaled..)) + 
  geom_density(color="blue",fill="blue",alpha=0.4)+
  facet_grid(type~as.factor(beta_num_in),scales = 'free_x')+
  ggtitle("Posterior Density of unimportant Beta/Sparsity")


beta_p100_all_cl %>% 
  dplyr::rename(Bay="beta_bay_abs",Bay_ori="beta_bay_ori",True="beta_true") %>% 
  gather(.,type,value,Bay:True,factor_key=TRUE) %>% 
  filter(beta_num_in%in%c(0,7),type=="Bay") %>% 
  ggplot(., aes(value,..scaled..)) + 
  geom_density(color="blue",fill="blue",alpha=0.4)+
  facet_grid(sparse~as.factor(beta_num_in),scales = 'free_x')+
  ggtitle("Posterior Density of beta0 and alpha")


###kappa plots###
kappa_p100_add<-kappa_p100_add %>% mutate(n=seq(1,length(beta_true)))

kappa_p100_all %>%
  filter(abs(beta_true)>0.9) %>% 
  filter(kappa_bay>=0.6) %>% 
  sample_n(length(kappa_bay)/1.2) %>% 
  dplyr::select(n) ->index_kappa

kappa_p100_all %>% filter(!(n %in% index_kappa$n)) %>%
  mutate(log_beta=log(abs(beta_true)))%>% 
  ggplot(., aes(x=log_beta,y=kappa_bay))+
  geom_point(alpha=0.5,size=1)+ 
  facet_grid(as.factor(sparse)~.)+
  ggtitle("Average kappa over beta setting")

kappa_p100_add%>% 
  mutate(log_beta=log(abs(beta_true)))%>% 
  ggplot(., aes(x=log_beta,y=kappa_mod))+
  geom_point(alpha=0.5,size=1)+ 
  facet_grid(as.factor(sparse)~.)+ 
  ggtitle("MAP of kappa over beta setting")

# ##total-edit##
# kappa_p100_all %>%mutate(log_beta=log(abs(beta_true)))%>% 
#   gather(.,type,value,kappa_bay:kappa_mod,factor_key=TRUE) %>% 
#   mutate(Meature=str_extract(type, "?+.(?=_)")) %>% 
#   ggplot(., aes(x=log_beta,y=kappa_bay))+
#   geom_point(alpha=0.5,size=1)+ 
#   ggtitle("Average kappa over beta setting")


####for factor compare##
only_letters <- function(x) { gsub("^([[:alpha:]]*).*$","\\1",x) }

ft_compare_p100_all%>% gather(.,type,value,Bay_1:true_2,factor_key=TRUE)%>% 
  group_by(month,mu_beta,type,sparse) %>% 
  dplyr::summarise(mean_abs_value=mean(abs(value))) %>% 
  mutate(factor_num=extract_numeric(type),
         factor_type=only_letters(type))->ft_compare_p



ft_compare_p %>% filter(sparse==0)%>% filter(factor_type!="grs") %>% 
  ggplot(aes(x=month,y=mean_abs_value,color=factor_type))+
  geom_line(size = 0.4,alpha=0.6)+
  facet_grid(factor_num~mu_beta,scales = "free_y")+
  ggtitle("Factor comparison for mean of 200 samples/No sparsity")


ft_compare_p %>% filter(sparse==1)%>% filter(factor_type!="grs") %>% 
  ggplot(aes(x=month,y=mean_abs_value,color=factor_type))+
  geom_line(size = 0.4,alpha=0.6)+
  facet_grid(factor_num~mu_beta,scales = "free_y")+
  ggtitle("Factor comparison for mean of 200 samples/Sparsity")


####RISFEE##

RISFEE_res<-RISFEE_p100_all %>% 
  gather(.,type,value,y_ris:f6_bay_ris,factor_key=TRUE) %>% 
  group_by(mu_beta,type,sparse) %>% 
  dplyr::summarise(mean_ris=mean(value))


###find best ###

#for beta

beta_p100_all  %>% filter(mu_beta==5,index!=53,index!=85,index!=70,index!=45,index!=15) %>% 
  mutate(abs_diff=abs(beta_bay_abs-abs(beta_true))) %>% 
  group_by(mu_beta,sparse,index) %>% 
  summarise(sum_abs_diff=sum(abs_diff)) %>% 
  group_by(mu_beta,sparse) %>% 
  slice(which.min(sum_abs_diff)) ->selected_index

selected_index

###using kappa

kappa_p100_all %>% filter(mu_beta==5, lambda_num%in%c(1,2,3,4)) %>%
  mutate(beta_true_new=ifelse(lambda_num%in%c(1,2),beta_true,0)) %>% 
  mutate(abs_diff=abs(beta_bay-abs(beta_true_new))) %>% 
  group_by(mu_beta,sparse,index) %>% 
  summarise(sum_abs_diff=sum(abs_diff)) %>% 
  group_by(mu_beta,sparse) %>% 
  slice(which.min(sum_abs_diff)) ->selected_index_kappa




####simulation 1/ mu_beta=5, sparse=0, index=53
setwd("/Users/LJY/Dropbox/PHD/horeshoe/condo_insample_new")
simu1_raw<-read_rds("simu_7_beta_5.rds")
result1_raw<-read_rds("result_7_beta_5.rds")

result1_raw$beta_all %>% data.frame() %>% 
  mutate(beta_num=c(0,rep(1,6),rep(2,6),rep(3,6),4),
         beta_num_in=c(0,rep(seq(1,6),3),7),
         beta_true=abs(simu1_raw$beta)) %>% 
  gather(.,draw,value,X1:X1000,factor_key=TRUE) %>% 
  mutate(abs_value=abs(value))->beta_all
  
head(beta_all)


open_beta<-function(index,sparse){
  
  if (sparse==0){
    data_name<-paste("simu_",index,"_beta_5",".rds",sep="")
    res_name<-paste("result_",index,"_beta_5",".rds",sep="")
    
  }else{
    data_name<-paste("simu_",index,"_beta_5","_sparse.rds",sep="")
    res_name<-paste("result_",index,"_beta_5","_sparse.rds",sep="")
  }
  
  simu1_raw<-read_rds(data_name)
  result1_raw<-read_rds(res_name)
  
result1_raw$beta_all %>% data.frame() %>% 
  mutate(beta_num=c(0,rep(1,6),rep(2,6),rep(3,6),4),
         beta_num_in=c(0,rep(seq(1,6),3),7),
         beta_true=abs(simu1_raw$beta)) %>% 
  gather(.,draw,value,X1:X1000,factor_key=TRUE) %>% 
  mutate(abs_value=abs(value),index=index,sparse=sparse)->beta_all

return(beta_all[c(1:5,20),])

}

testFunction <- function (index,sparse) {
  return(tryCatch(open_beta(index,sparse), error=function(e) NULL))
}

# tt_expand=expand.grid(index=seq(1,100),sparse=c(0,1))
# beta_choose<-mapply(testFunction,index=tt_expand$index,
#                     sparse=tt_expand$sparse)
# 
# 
# beta_choose<-beta_choose %>% bind_rows() 
# 
# beta_choose%>% filter(beta_num_in==2,sparse==0,abs_value>1) %>% 
#   dplyr::select(index)->index_0
# 
# beta_choose%>% filter(beta_num_in==2,sparse==1,abs_value>1) %>% 
#   dplyr::select(index)->index_1
# 
# beta_choose %>% filter(index%in% index_1$index,sparse==1) %>% view()



###############################################################
setwd("/Users/LJY/Dropbox/PHD/horeshoe/condo_insample_new")

simu6_raw<-read_rds("simu_19_beta_0.1_sparse.rds")
result6_raw<-read_rds("result_19_beta_0.1_sparse.rds")

# result1_raw$beta_all %>% data.frame() %>% 
#   mutate(beta_num=c(0,rep(1,6),rep(2,6),rep(3,6),4),
#          beta_num_in=c(0,rep(seq(1,6),3),7),
#          beta_true=abs(simu1_raw$beta)) %>% 
#   gather(.,draw,value,X1:X1000,factor_key=TRUE) %>% 
#   mutate(abs_value=abs(value))->beta_all
# 
# 
# 
# vline_beta <- summarise(group_by(beta_all,beta_num_in),
#                       true=mean(beta_true))
# 
# 
# vline_beta1<-vline_beta %>% filter(beta_num_in %in% c(1,2))
# vline_beta2<-vline_beta %>% filter(beta_num_in %in% c(3,4,5,6))
# vline_beta3<-vline_beta %>% filter(beta_num_in %in% c(0,7))
# 
# beta_all %>% 
#   filter(beta_num_in%in%c(1,2)) %>% 
#   ggplot(., aes(abs_value,..scaled..)) + 
#   geom_density(color="blue",fill="blue",alpha=0.4)+
#   geom_vline(data=vline_beta1, aes(xintercept = true),
#              col='red',size=0.4,alpha=0.4)+
#   facet_grid(as.factor(beta_num_in)~.,scales = 'free_x')+
#   ggtitle("Posterior Density of important Beta/Simulation 2")
# 
# beta_all %>% 
#   filter(beta_num_in%in%c(3,4,5,6)) %>% 
#   ggplot(., aes(abs_value,..scaled..)) + 
#   geom_density(color="blue",fill="blue",alpha=0.4)+
#   facet_grid(as.factor(beta_num_in)~.,scales = 'free_x')+
#   ggtitle("Posterior Density of unimportant Beta/Simulation 2")
# 
# beta_all %>% 
#   filter(beta_num_in%in%c(0,7)) %>% 
#   ggplot(., aes(abs_value,..scaled..)) + 
#   geom_density(color="blue",fill="blue",alpha=0.4)+
#   geom_vline(data=vline_beta3, aes(xintercept = true),
#              col='red',size=0.4,alpha=0.4)+
#   facet_grid(as.factor(beta_num_in)~.,scales = 'free_x')+
#   ggtitle("Posterior Density of beta0 and alpha/Simulation 2")
# 
# 
# ###45 degree line for beta##
# 
# beta_all%>% group_by(beta_num,beta_num_in) %>% 
#   summarise(beta_true=mean(beta_true),
#             beta_est=mean(abs_value)) %>% 
#   ggplot(aes(x=beta_true,y=beta_est))+
#   geom_point(alpha=0.6)+ggtitle("True beta VS esimated beta/simu4")



##kappa##
kappa_res<-data.frame(1/(1+t(result1_raw$lambda_all)^2))

kappa_res<-kappa_res %>% mutate(kappa=c("kappa1","kappa2","kappa3","kappa4","kappa5","kappa6"),
                     beta=c("beta1","beta2","beta3","beta4","beta5","beta6"),
                     beta_true=simu1_raw$beta[2:7])


kappa_res %>% gather(.,draw,value,X1:X1000,factor_key=T)%>%
  ggplot(aes(x=beta, y=value),alpha=0.4)+ 
  geom_boxplot()+
  ylim(c(0,1))+ggtitle("Boxplots for Kappa/Simu6")


###factors###
##########plot factors-first two####
factor1_bay6<-result6_raw$factor1_all
factor1_bay6<-result6_raw$factor2_all

factor1_bay<-apply(abs(result6_raw$factor1_all),2,mean)
factor2_bay<-apply(abs(result6_raw$factor2_all),2,mean)

#FT_grs<-abs(result1_raw$ft_grs)

t<-length(factor1_bay6)

factor_simu6<-data.frame(cbind(factor1_bay,factor2_bay,abs(simu6_raw$FT[1,1:t]),
                              abs(simu6_raw$FT[2,1:t])))



factor_simu6 %>% dplyr::rename(f1_bay=factor1_bay,f2_bay=factor2_bay,
                              f1_true=V3,f2_true=V4) %>% 
  mutate(month=seq(1,length(f1_bay)))%>% 
  gather(.,type,value,f1_bay:f2_true,factor_key = T) %>% 
  mutate(factor=str_extract(type, ".+?(?=_)"),
         method=str_sub(type, 4))->factor_p6



factor_p6 %>% filter(method%in%c("bay","true"))%>% 
  ggplot(aes(x=month,y=abs(value),color=method))+
  geom_line(size = 0.4)+facet_wrap(.~factor,scales = "free_y",ncol=1)+
  ggtitle("Factor comparison for simulation 6")


#In sample fit##

y_insample<-apply(result6_raw$y_insample,2, mean)

y_insample %>% data.frame() %>% 
  dplyr::rename(y_bay=".")%>%
  mutate(y_true=as.numeric(simu6_raw$y[1:length(y_bay)]),
         quarter=seq(1,length(y_bay))) %>% 
  gather(.,type,value,y_bay:y_true,factor_key=T) %>% 
  mutate(method=substring(type, first = 3)) %>% 
  ggplot(aes(x=quarter,y=value,color=method))+
  ylim(c(-15,15))+
  geom_line(size=0.6,alpha=0.7)+
  ggtitle("In sample fit / Simulation 6")



######RISFEE#####
k<-length(y_insample)
data.frame("y_ris"=sum(abs(simu1_raw$y[1:k]-y_insample))/sum(abs(simu1_raw$y[1:k])),
           "f1_bay_ris"=sum(abs(abs(simu1_raw$FT[1,1:t])-factor1_bay))/sum(abs(simu1_raw$FT[1,1:t])),
           "f2_bay_ris"=sum(abs(abs(simu1_raw$FT[2,1:t])-factor2_bay))/sum(abs(simu1_raw$FT[2,1:t])))
                   





####A and Theta##

theta_res<-matrix(apply(abs(result1_raw$theta_all),2,mean),nrow=n)

c(theta_res[,1:2]) %>% data.frame() %>% 
  dplyr::rename(theta_bay=".") %>% 
  mutate(theta_true=abs(c(simu1_raw$The[,1:2]))) %>% 
  ggplot(.,aes(x=theta_true,y=theta_bay))+geom_point(color="gray")+
  ggtitle("Important theta VS true/ Simulation 2")

##first three are accurate##
a_test<-apply(result1_raw$a_all,2,mean)
a_true<-abs(diag(simu1_raw$A))

#############box-plots for kappa 100 sample path

kappa_p100_all %>% dplyr::rename(beta_num=lambda_num) %>% 
  filter(!((mu_beta==5)&(beta_num%in%c(1,2))&(kappa_med>=0.85)&(kappa_med<1)))%>% 
  ggplot(aes(x=as.factor(beta_num), y=kappa_med),alpha=0.4)+ 
  geom_boxplot()+facet_grid(mu_beta~sparse)+
  ggtitle("Boxplots for Kappa medians/ALL")


kappa_p100_all %>% dplyr::rename(beta_num=lambda_num) %>% 
  filter(!((mu_beta==5)&(beta_num%in%c(1,2))&(kappa_med>=0.85)&(kappa_med<1)))%>% 
  ggplot(aes(x=as.factor(beta_num), y=kappa_med),alpha=0.4)+ 
  geom_boxplot(outlier.shape = NA)+facet_grid(mu_beta~sparse)+
  ggtitle("Boxplots for Kappa median-no outlier/ALL")








