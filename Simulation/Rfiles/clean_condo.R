
setwd("/Users/LJY/Dropbox/PHD/horeshoe/condo_insample_new")

simu1<-readRDS("simu_37_beta_5.rds")
factors_initial_raw1<-readRDS("simu_37_beta_5_ft.rds")

simu2<-readRDS("simu_76_beta_5_sparse.rds")
factors_initial_raw2<-readRDS("simu_76_beta_5_sparse_ft.rds")


simu3<-readRDS("simu_25_beta_1.rds")
factors_initial_raw3<-readRDS("simu_25_beta_1_ft.rds")

simu4<-readRDS("simu_56_beta_1_sparse.rds")
factors_initial_raw4<-readRDS("simu_56_beta_1_sparse_ft.rds")

simu5<-readRDS("simu_75_beta_0.1.rds")
factors_initial_raw5<-readRDS("simu_75_beta_0.1_ft.rds")

simu6<-readRDS("simu_19_beta_0.1_sparse.rds")
factors_initial_raw6<-readRDS("simu_19_beta_0.1_sparse_ft.rds")



setwd("/Users/LJY/Dropbox/PHD/horeshoe/nowcasts") 
##/nowcasts_new for simu1 and simu2

nowcast_clean<-function(m,q,rel,type){
  
  res_name=paste("simu",type,"_m",m,"_q",q,"_rel",rel,".rds",sep="")
  res_raw<-readRDS(res_name)
  
  ###kappa###
  kappa<-apply(1/(1+res_raw$lambda_all^2),2,mean)
  y_pred=mean(res_raw$y_pred)
  #beta<-apply(res_raw$beta_all,1,mean)
  
  
  kappa %>% data.frame() %>% dplyr::rename(kappa=".") %>% 
    mutate(lam_num=seq(1,6),month=m,quarter=q,rel=rel)->res_kappa
  
  t(res_raw$beta_all) %>% data.frame()  %>% 
    mutate(iteration=seq(1,1000),month=m,quarter=q,rel=rel)->res_beta
  #nowcast
  y_pred %>% data.frame()%>% dplyr::rename(y_pred=".") %>% 
    mutate(month=m,quarter=q,rel=rel)->res_y
  
  return(res=list("kappa"=res_kappa,"pred_y"=res_y,"beta"=res_beta))
}


testFunction <- function (m,q,rel,type) {
  return(tryCatch(nowcast_clean(m,q,rel,type), error=function(e) NULL))
}

df<-expand.grid(q=seq(1,20),m=c(1,2,3),rel=c(1,2,3))

res_simu1_raw<-mapply(testFunction,m=df$m,q=df$q,
                      rel=df$rel,type=1)

res_simu2_raw<-mapply(testFunction,m=df$m,q=df$q,
                      rel=df$rel,type=2)

##
res_simu3_raw<-mapply(testFunction,m=df$m,q=df$q,
                 rel=df$rel,type=3)

res_simu4_raw<-mapply(testFunction,m=df$m,q=df$q,
                      rel=df$rel,type=4)

res_simu5_raw<-mapply(testFunction,m=df$m,q=df$q,
                      rel=df$rel,type=5)

res_simu6_raw<-mapply(testFunction,m=df$m,q=df$q,
                      rel=df$rel,type=6)

###simu1
kappa_simu1_all<-NULL
y_simu1_all<-NULL
beta_simu1_all<-NULL
for (i in 1:180){
  res_sep<-res_simu1_raw[,i]
  kappa_simu1_all<-rbind(kappa_simu1_all,res_sep$kappa)
  y_simu1_all<-rbind(y_simu1_all,res_sep$pred_y)
  beta_simu1_all <-rbind(beta_simu1_all,res_sep$beta)
}

##simu2
kappa_simu2_all<-NULL
y_simu2_all<-NULL
beta_simu2_all<-NULL
for (i in 1:180){
  res_sep<-res_simu2_raw[,i]
  kappa_simu2_all<-rbind(kappa_simu2_all,res_sep$kappa)
  y_simu2_all<-rbind(y_simu2_all,res_sep$pred_y)
  beta_simu2_all <- rbind(beta_simu2_all,res_sep$beta)
}

#simu3
kappa_simu3_all<-NULL
y_simu3_all<-NULL
beta_simu3_all<-NULL
for (i in 1:180){
  res_sep<-res_simu3_raw[,i]
  kappa_simu3_all<-rbind(kappa_simu3_all,res_sep$kappa)
  y_simu3_all<-rbind(y_simu3_all,res_sep$pred_y)
  beta_simu3_all <-rbind(beta_simu3_all,res_sep$beta)
}

##simu4
kappa_simu4_all<-NULL
y_simu4_all<-NULL
beta_simu4_all<-NULL
for (i in 1:180){
  res_sep<-res_simu4_raw[,i]
  kappa_simu4_all<-rbind(kappa_simu4_all,res_sep$kappa)
  y_simu4_all<-rbind(y_simu4_all,res_sep$pred_y)
  beta_simu4_all <-rbind(beta_simu4_all,res_sep$beta)
}

#simu5
kappa_simu5_all<-NULL
y_simu5_all<-NULL
beta_simu5_all<-NULL
for (i in 1:180){
  res_sep<-res_simu5_raw[,i]
  kappa_simu5_all<-rbind(kappa_simu5_all,res_sep$kappa)
  y_simu5_all<-rbind(y_simu5_all,res_sep$pred_y)
  beta_simu5_all <-rbind(beta_simu5_all,res_sep$beta)
}



kappa_simu6_all<-NULL
y_simu6_all<-NULL
beta_simu6_all<-NULL
for (i in 1:180){
  res_sep<-res_simu6_raw[,i]
  kappa_simu6_all<-rbind(kappa_simu6_all,res_sep$kappa)
  y_simu6_all<-rbind(y_simu6_all,res_sep$pred_y)
  beta_simu6_all <-rbind(beta_simu6_all,res_sep$beta)
}

###for true####
#simu1
y_true1<-simu1$y

y_true1 %>% data.frame() %>% 
  dplyr::rename(y_true=".") %>% 
  mutate(quarter=seq(1,length(y_true)))%>% 
  filter(quarter<=60)->y_true1_f

#simu2
y_true2<-simu2$y

y_true2 %>% data.frame() %>% 
  dplyr::rename(y_true=".") %>% 
  mutate(quarter=seq(1,length(y_true)))%>% 
  filter(quarter<=60)->y_true2_f


#simu3
y_true3<-simu3$y

y_true3 %>% data.frame() %>% 
  dplyr::rename(y_true=".") %>% 
  mutate(quarter=seq(1,length(y_true)))%>% 
  filter(quarter<=60)->y_true3_f

#simu4
y_true4<-simu4$y

y_true4 %>% data.frame() %>% 
  dplyr::rename(y_true=".") %>% 
  mutate(quarter=seq(1,length(y_true)))%>% 
  filter(quarter<=60)->y_true4_f

#simu5
y_true5<-simu5$y

y_true5 %>% data.frame() %>% 
  dplyr::rename(y_true=".") %>% 
  mutate(quarter=seq(1,length(y_true)))%>% 
  filter(quarter<=60)->y_true5_f

##simu6
y_true6<-simu6$y

y_true6 %>% data.frame() %>% 
  dplyr::rename(y_true=".") %>% 
  mutate(quarter=seq(1,length(y_true)))%>% 
  filter(quarter<=60)->y_true6_f


###add in sample quarters##
cbp2 <- c("#999999", "#E69F00", "#009E73",
           "#CC79A7","#000000")

##simu1##
y_simu1_all%>%  mutate(quarter=quarter+40)%>%
  mutate_at(c("month","rel"),as.factor) %>% 
  left_join(.,y_true1_f,by="quarter")%>% 
  spread(.,rel,y_pred)%>% 
  dplyr::rename(rel1="1",rel2="2",rel3="3") %>% 
  gather(.,type,value,y_true:rel3,factor_key=TRUE)->y_simu1_mod


y_true1_f %>% filter((quarter<41)&(quarter>35)) %>% 
  gather(.,type,value,y_true,factor_key=TRUE) %>%
  slice(rep(1:n(), each=3)) %>% 
  mutate(month=rep(c(1,2,3),length(unique(value)))) %>% 
  dplyr::select(month,quarter,type,value) %>% arrange(month) %>% 
  rbind(.,y_simu1_mod)->y_simu1_p


y_simu1_p%>% 
  ggplot(.,aes(x=quarter,y=value,color=type,linetype=type,shape=type))+
  geom_line(size=0.5)+
  geom_point(size=1)+
  scale_linetype_manual(labels = c("True", "RL1",'RL2','RL3'),values = c(1,2,1,2))+
  scale_shape_manual(labels = c("True", "RL1",'RL2','RL3'),values=c(17,17,16,16))+
  scale_colour_manual(labels = c("True", "RL1",'RL2','RL3'),values=cbp2)+
  geom_vline(xintercept=41,color="red",linetype='dotted',alpha=0.6)+
  facet_grid(month~.)+ggtitle("Simulation 1")+
  ylab("GDP")+
  theme(text=element_text(size=10), axis.text=element_text(size=7), 
        axis.title=element_text(size=10), plot.title=element_text(size=10))->p_nowcast_s1

##simu2
y_simu2_all%>%  mutate(quarter=quarter+40)%>%
  mutate_at(c("month","rel"),as.factor) %>% 
  left_join(.,y_true2_f,by="quarter")%>% 
  spread(.,rel,y_pred)%>% 
  dplyr::rename(rel1="1",rel2="2",rel3="3") %>% 
  gather(.,type,value,y_true:rel3,factor_key=TRUE)->y_simu2_mod


y_true2_f %>% filter((quarter<41)&(quarter>35)) %>% 
  gather(.,type,value,y_true,factor_key=TRUE) %>%
  slice(rep(1:n(), each=3)) %>% 
  mutate(month=rep(c(1,2,3),length(unique(value)))) %>% 
  dplyr::select(month,quarter,type,value) %>% arrange(month) %>% 
  rbind(.,y_simu2_mod)->y_simu2_p


y_simu2_p%>% 
  ggplot(.,aes(x=quarter,y=value,color=type,linetype=type,shape=type))+
  geom_line(size=0.5)+
  geom_point(size=1)+
  scale_linetype_manual(labels = c("True", "RL1",'RL2','RL3'),values = c(1,2,1,2))+
  scale_shape_manual(labels = c("True", "RL1",'RL2','RL3'),values=c(17,17,16,16))+
  scale_colour_manual(labels = c("True", "RL1",'RL2','RL3'),values=cbp2)+
  geom_vline(xintercept=41,color="red",linetype='dotted',alpha=0.6)+
  facet_grid(month~.)+ggtitle("Simulation 2")+
  ylab("GDP")+
  theme(text=element_text(size=10), axis.text=element_text(size=7), 
        axis.title=element_text(size=10), plot.title=element_text(size=10))->p_nowcast_s2


##simu3
y_simu3_all%>%  mutate(quarter=quarter+40)%>%
  mutate_at(c("month","rel"),as.factor) %>% 
  left_join(.,y_true3_f,by="quarter")%>% 
  spread(.,rel,y_pred)%>% 
  dplyr::rename(rel1="1",rel2="2",rel3="3") %>% 
  gather(.,type,value,y_true:rel3,factor_key=TRUE)->y_simu3_mod


y_true3_f %>% filter((quarter<41)&(quarter>35)) %>% 
  gather(.,type,value,y_true,factor_key=TRUE) %>%
  slice(rep(1:n(), each=3)) %>% 
  mutate(month=rep(c(1,2,3),length(unique(value)))) %>% 
  dplyr::select(month,quarter,type,value) %>% arrange(month) %>% 
  rbind(.,y_simu3_mod)->y_simu3_p


y_simu3_p%>% 
  ggplot(.,aes(x=quarter,y=value,color=type,linetype=type,shape=type))+
  geom_line(size=0.5)+
  geom_point(size=1)+
  scale_linetype_manual(labels = c("True", "RL1",'RL2','RL3'),values = c(1,2,1,2))+
  scale_shape_manual(labels = c("True", "RL1",'RL2','RL3'),values=c(17,17,16,16))+
  scale_colour_manual(labels = c("True", "RL1",'RL2','RL3'),values=cbp2)+
  geom_vline(xintercept=41,color="red",linetype='dotted',alpha=0.6)+
  facet_grid(month~.)+ggtitle("Simulation 3")+
  ylab("GDP")+
  theme(text=element_text(size=10), axis.text=element_text(size=7), 
        axis.title=element_text(size=10), plot.title=element_text(size=10))->p_nowcast_s3

##simu4##
y_simu4_all%>%  mutate(quarter=quarter+40)%>%
  mutate_at(c("month","rel"),as.factor) %>% 
  left_join(.,y_true4_f,by="quarter")%>% 
  spread(.,rel,y_pred)%>% 
  dplyr::rename(rel1="1",rel2="2",rel3="3") %>% 
  gather(.,type,value,y_true:rel3,factor_key=TRUE)->y_simu4_mod


y_true4_f %>% filter((quarter<41)&(quarter>35)) %>% 
  gather(.,type,value,y_true,factor_key=TRUE) %>%
  slice(rep(1:n(), each=3)) %>% 
  mutate(month=rep(c(1,2,3),length(unique(value)))) %>% 
  dplyr::select(month,quarter,type,value) %>% arrange(month) %>% 
  rbind(.,y_simu4_mod)->y_simu4_p


y_simu4_p%>% 
  ggplot(.,aes(x=quarter,y=value,color=type,linetype=type,shape=type))+
  geom_line(size=0.5)+
  geom_point(size=1)+
  scale_linetype_manual(labels = c("True", "RL1",'RL2','RL3'),values = c(1,2,1,2))+
  scale_shape_manual(labels = c("True", "RL1",'RL2','RL3'),values=c(17,17,16,16))+
  scale_colour_manual(labels = c("True", "RL1",'RL2','RL3'),values=cbp2)+
  geom_vline(xintercept=41,color="red",linetype='dotted',alpha=0.6)+
  facet_grid(month~.)+ggtitle("Simulation 4")+
  ylab("GDP")+
  theme(text=element_text(size=10), axis.text=element_text(size=7), 
        axis.title=element_text(size=10), plot.title=element_text(size=10))->p_nowcast_s4


##simu5##

y_simu5_all%>%  mutate(quarter=quarter+40)%>%
  mutate_at(c("month","rel"),as.factor) %>% 
  left_join(.,y_true5_f,by="quarter")%>% 
  spread(.,rel,y_pred)%>% 
  dplyr::rename(rel1="1",rel2="2",rel3="3") %>% 
  gather(.,type,value,y_true:rel3,factor_key=TRUE)->y_simu5_mod


y_true5_f %>% filter((quarter<41)&(quarter>35)) %>% 
  gather(.,type,value,y_true,factor_key=TRUE) %>%
  slice(rep(1:n(), each=3)) %>% 
  mutate(month=rep(c(1,2,3),length(unique(value)))) %>% 
  dplyr::select(month,quarter,type,value) %>% arrange(month) %>% 
  rbind(.,y_simu5_mod)->y_simu5_p


y_simu5_p%>% 
  ggplot(.,aes(x=quarter,y=value,color=type,linetype=type,shape=type))+
  geom_line(size=0.5)+
  geom_point(size=1)+
  scale_linetype_manual(labels = c("True", "RL1",'RL2','RL3'),values = c(1,2,1,2))+
  scale_shape_manual(labels = c("True", "RL1",'RL2','RL3'),values=c(17,17,16,16))+
  scale_colour_manual(labels = c("True", "RL1",'RL2','RL3'),values=cbp2)+
  geom_vline(xintercept=41,color="red",linetype='dotted',alpha=0.6)+
  facet_grid(month~.)+ggtitle("Simulation 5")+
  ylab("GDP")+
  theme(text=element_text(size=10), axis.text=element_text(size=7), 
        axis.title=element_text(size=10), plot.title=element_text(size=10))->p_nowcast_s5

####simu6
y_simu6_all%>%  mutate(quarter=quarter+40)%>%
  mutate_at(c("month","rel"),as.factor) %>% 
  left_join(.,y_true6_f,by="quarter")%>% 
  spread(.,rel,y_pred)%>% 
  dplyr::rename(rel1="1",rel2="2",rel3="3") %>% 
  gather(.,type,value,y_true:rel3,factor_key=TRUE)->y_simu6_mod


y_true6_f %>% filter((quarter<41)&(quarter>35)) %>% 
  gather(.,type,value,y_true,factor_key=TRUE) %>%
  slice(rep(1:n(), each=3)) %>% 
  mutate(month=rep(c(1,2,3),length(unique(value)))) %>% 
  dplyr::select(month,quarter,type,value) %>% arrange(month) %>% 
  rbind(.,y_simu6_mod)->y_simu6_p


y_simu6_p%>% 
  ggplot(.,aes(x=quarter,y=value,color=type,linetype=type,shape=type))+
  geom_line(size=0.5)+
  geom_point(size=1)+
  scale_linetype_manual(labels = c("True", "RL1",'RL2','RL3'),values = c(1,2,1,2))+
  scale_shape_manual(labels = c("True", "RL1",'RL2','RL3'),values=c(17,17,16,16))+
  scale_colour_manual(labels = c("True", "RL1",'RL2','RL3'),values=cbp2)+
  geom_vline(xintercept=41,color="red",linetype='dotted',alpha=0.6)+
  facet_grid(month~.)+ggtitle("Simulation 6")+
  ylab('GDP')+
  theme(text=element_text(size=10), axis.text=element_text(size=7), 
        axis.title=element_text(size=10), plot.title=element_text(size=10))->p_nowcast_s6


b0a_den<-ggarrange(p_nowcast_s1,p_nowcast_s2,p_nowcast_s3,
                   p_nowcast_s4,p_nowcast_s5,p_nowcast_s6, 
                   ncol = 2, nrow = 3,common.legend = TRUE,legend="bottom")
annotate_figure(b0a_den,
                top = text_grob("One-step ahead nowcasts for 6 simulations", 
                                face = "bold", size = 10))




####MANE##

y_true6_f %>% mutate(y_rw=c(0,y_true[-length(y_true)]))->y_true6_all


y_simu6_all %>% mutate(quarter=quarter+40) %>% 
  left_join(.,y_true6_all,by="quarter")->y_simu6_result_combine
  

y_simu6_result_combine  %>% group_by(month,rel) %>% 
  summarise(MANE=mean(abs(y_pred-y_true)),MANE_rw=mean(abs(y_rw-y_true))) %>% 
  mutate(MANE_ratio=MANE/MANE_rw)->mane_bay_simu6

mane_bay_simu6 %>% mutate(mane_reduced=1-MANE_ratio) %>%
  dplyr::select(month, rel, mane_reduced) %>% group_by(month,rel) %>% 
  summarise(mean_mane=mean(mane_reduced))


p<-ggplot(data=mane_bay_simu6, aes(x=rel, y=MANE_ratio,fill=as.factor(rel))) +
  geom_bar(stat="identity")+facet_grid(.~month)

p+coord_cartesian(ylim = c(0,1.5))+scale_fill_grey()+
  ggtitle("MANE Ratio for Simu6")+ ylab("MANE Ratio")+labs(fill = "release")+
  geom_hline(yintercept=1, linetype="dashed", color = "red")


#####beta and kappa##
library(ggbreak)
##simu1
beta_simu1_all %>% group_by(month,quarter,rel) %>% 
  summarise(across(everything(), mean)) %>% 
  select(-iteration,-X1,-X20) %>% 
  pivot_longer(.,4:21,values_to="beta",names_to="index") %>% 
  mutate(index=parse_number(index)-1) %>% 
  mutate(lam_num=index%%6,beta_num=floor(index/6)+1) %>% 
  mutate(lam_num=ifelse(lam_num==0,6,lam_num))->beta_simu1_f

beta_simu1_f %>% left_join(kappa_simu1_all,by=c("month","quarter","rel","lam_num")) %>% 
  arrange(lam_num,beta_num) %>% 
  mutate(index_beta=10*lam_num+beta_num) %>% 
  ggplot(.,aes(x=abs(beta),y=kappa))+
  geom_point(alpha=0.4)+
  ggtitle("Simulation 1")+
  xlim(0, 8)+
  xlab("beta")+
  ylim(0,1)+
  theme(text=element_text(size=10), axis.text=element_text(size=7), 
        axis.title=element_text(size=7), plot.title=element_text(size=7))->k_b_simu1



simu1$beta[-c(1,20)] %>% data.frame() %>% 
  rename(beta=".") %>% mutate(lam_num=rep(seq(1,6),3)) %>% 
  right_join(kappa_simu1_all,by="lam_num") %>% 
  ggplot(.,aes(x=abs(beta),y=kappa))+
  geom_point(alpha=0.4)+
  ggtitle("Simulation 1")+
  scale_x_break(c(0.5,4.5))+
  theme(text=element_text(size=10), axis.text=element_text(size=7), 
        axis.title=element_text(size=7), plot.title=element_text(size=7))->k_b_simu1



###simu2##

beta_simu2_all %>% group_by(month,quarter,rel) %>% 
  summarise(across(everything(), mean)) %>% 
  select(-iteration,-X1,-X20) %>% 
  pivot_longer(.,4:21,values_to="beta",names_to="index") %>% 
  mutate(index=parse_number(index)-1) %>% 
  mutate(lam_num=index%%6,beta_num=floor(index/6)+1) %>% 
  mutate(lam_num=ifelse(lam_num==0,6,lam_num))->beta_simu2_f

beta_simu2_f %>% left_join(kappa_simu2_all,by=c("month","quarter","rel","lam_num")) %>% 
  arrange(lam_num,beta_num) %>% 
  mutate(index_beta=10*lam_num+beta_num) %>% 
  ggplot(.,aes(x=abs(beta),y=kappa))+
  geom_point(alpha=0.4)+
  xlab("beta")+
  ggtitle("Simulation 2")+
  xlim(0, 8)+
  ylim(0,1)+
  theme(text=element_text(size=10), axis.text=element_text(size=7), 
        axis.title=element_text(size=7), plot.title=element_text(size=7))->k_b_simu2

###simu3##

beta_simu3_all %>% group_by(month,quarter,rel) %>% 
  summarise(across(everything(), mean)) %>% 
  select(-iteration,-X1,-X20) %>% 
  pivot_longer(.,4:21,values_to="beta",names_to="index") %>% 
  mutate(index=parse_number(index)-1) %>% 
  mutate(lam_num=index%%6,beta_num=floor(index/6)+1) %>% 
  mutate(lam_num=ifelse(lam_num==0,6,lam_num))->beta_simu3_f

beta_simu3_f %>% left_join(kappa_simu3_all,by=c("month","quarter","rel","lam_num")) %>% 
  arrange(lam_num,beta_num) %>% 
  mutate(index_beta=10*lam_num+beta_num) %>% 
  ggplot(.,aes(x=abs(beta),y=kappa))+
  geom_point(alpha=0.4)+
  xlab("beta")+
  ggtitle("Simulation 3")+
  xlim(0, 1.5)+
  ylim(0,1)+
  theme(text=element_text(size=10), axis.text=element_text(size=7), 
        axis.title=element_text(size=7), plot.title=element_text(size=7))->k_b_simu3



###simu4##

beta_simu4_all %>% group_by(month,quarter,rel) %>% 
  summarise(across(everything(), mean)) %>% 
  select(-iteration,-X1,-X20) %>% 
  pivot_longer(.,4:21,values_to="beta",names_to="index") %>% 
  mutate(index=parse_number(index)-1) %>% 
  mutate(lam_num=index%%6,beta_num=floor(index/6)+1) %>% 
  mutate(lam_num=ifelse(lam_num==0,6,lam_num))->beta_simu4_f

beta_simu4_f %>% left_join(kappa_simu4_all,by=c("month","quarter","rel","lam_num")) %>% 
  arrange(lam_num,beta_num) %>% 
  mutate(index_beta=10*lam_num+beta_num) %>% 
  ggplot(.,aes(x=abs(beta),y=kappa))+
  geom_point(alpha=0.4)+
  xlab("beta")+
  ggtitle("Simulation 4")+
  xlim(0, 1.5)+
  ylim(0,1)+
  theme(text=element_text(size=10), axis.text=element_text(size=7), 
        axis.title=element_text(size=7), plot.title=element_text(size=7))->k_b_simu4

##simu5##

beta_simu5_all %>% group_by(month,quarter,rel) %>% 
  summarise(across(everything(), mean)) %>% 
  select(-iteration,-X1,-X20) %>% 
  pivot_longer(.,4:21,values_to="beta",names_to="index") %>% 
  mutate(index=parse_number(index)-1) %>% 
  mutate(lam_num=index%%6,beta_num=floor(index/6)+1) %>% 
  mutate(lam_num=ifelse(lam_num==0,6,lam_num))->beta_simu5_f

beta_simu5_f %>% left_join(kappa_simu5_all,by=c("month","quarter","rel","lam_num")) %>% 
  arrange(lam_num,beta_num) %>% 
  mutate(index_beta=10*lam_num+beta_num) %>% 
  ggplot(.,aes(x=abs(beta),y=kappa))+
  geom_point(alpha=0.4)+
  xlab("beta")+
  ggtitle("Simulation 5")+
  xlim(0, 0.6)+
  ylim(0,1)+
  theme(text=element_text(size=10), axis.text=element_text(size=7), 
        axis.title=element_text(size=7), plot.title=element_text(size=7))->k_b_simu5

###simu6##

beta_simu6_all %>% group_by(month,quarter,rel) %>% 
  summarise(across(everything(), mean)) %>% 
  select(-iteration,-X1,-X20) %>% 
  pivot_longer(.,4:21,values_to="beta",names_to="index") %>% 
  mutate(index=parse_number(index)-1) %>% 
  mutate(lam_num=index%%6,beta_num=floor(index/6)+1) %>% 
  mutate(lam_num=ifelse(lam_num==0,6,lam_num))->beta_simu6_f

beta_simu6_f %>% left_join(kappa_simu6_all,by=c("month","quarter","rel","lam_num")) %>% 
  arrange(lam_num,beta_num) %>% 
  mutate(index_beta=10*lam_num+beta_num) %>% 
  ggplot(.,aes(x=abs(beta),y=kappa))+
  geom_point(alpha=0.4)+
  xlab("beta")+
  ggtitle("Simulation 6")+
  xlim(0, 0.6)+
  ylim(0,1)+
  theme(text=element_text(size=10), axis.text=element_text(size=7), 
        axis.title=element_text(size=7), plot.title=element_text(size=7))->k_b_simu6


beta_kappa<-ggarrange(k_b_simu1,k_b_simu2,k_b_simu3,
                      k_b_simu4,k_b_simu5,k_b_simu6,
                     ncol = 2, nrow = 3,common.legend = TRUE,legend="bottom")

annotate_figure(beta_kappa)



###kappa result#
kappa_simu1_all %>% 
  ggplot(aes(x=as.factor(lam_num),y=kappa))+
  geom_boxplot()+ylim(c(0,1))+
  xlab("Number")+
  ggtitle("Simulation 1")+
  theme(text=element_text(size=10), axis.text=element_text(size=7), 
        axis.title=element_text(size=7), plot.title=element_text(size=7))->p_kappa_s1


kappa_simu2_all%>% 
  ggplot(aes(x=as.factor(lam_num),y=kappa))+
  geom_boxplot()+ylim(c(0,1))+
  xlab("Number")+
  ggtitle("Simulation 2")+
  theme(text=element_text(size=10), axis.text=element_text(size=7), 
        axis.title=element_text(size=7), plot.title=element_text(size=7))->p_kappa_s2

kappa_simu3_all %>% 
  ggplot(aes(x=as.factor(lam_num),y=kappa))+
  geom_boxplot()+ylim(c(0,1))+
  xlab("Number")+
  ggtitle("Simulation 3")+
  theme(text=element_text(size=10), axis.text=element_text(size=7), 
        axis.title=element_text(size=7), plot.title=element_text(size=7))->p_kappa_s3

kappa_simu4_all %>% 
  ggplot(aes(x=as.factor(lam_num),y=kappa))+
  geom_boxplot()+ylim(c(0,1))+
  xlab("Number")+
  ggtitle("Simulation 4")+
  theme(text=element_text(size=10), axis.text=element_text(size=7), 
        axis.title=element_text(size=7), plot.title=element_text(size=7))->p_kappa_s4


kappa_simu5_all %>% 
  ggplot(aes(x=as.factor(lam_num),y=kappa))+
  geom_boxplot()+ylim(c(0,1))+
  xlab("Number")+
  ggtitle("Simulation 5")+
  theme(text=element_text(size=10), axis.text=element_text(size=7), 
        axis.title=element_text(size=7), plot.title=element_text(size=7))->p_kappa_s5

kappa_simu6_all %>% 
  ggplot(aes(x=as.factor(lam_num),y=kappa))+
  geom_boxplot(size=0.1)+ylim(c(0,1))+
  xlab("Number")+
  ggtitle("Simulation 6")+
  theme(text=element_text(size=10), axis.text=element_text(size=7), 
        axis.title=element_text(size=7), plot.title=element_text(size=7))->p_kappa_s6


b0a_kappa<-ggarrange(p_kappa_s1,p_kappa_s2,p_kappa_s3,
                   p_kappa_s4,p_kappa_s5,p_kappa_s6, 
                   ncol = 2, nrow = 3,common.legend = TRUE,legend="bottom")
annotate_figure(b0a_kappa,
                top = text_grob("Kappa estimations for 6 simulations", 
                                face = "bold", size = 7))


####plot factors##
setwd("/Users/LJY/Dropbox/PHD/horeshoe/condo_insample_new")

simu1_raw<-read_rds("simu_37_beta_5.rds")
result1_raw<-read_rds("result_37_beta_5.rds")

simu2_raw<-read_rds("simu_76_beta_5_sparse.rds")
result2_raw<-read_rds("result_76_beta_5_sparse.rds")

simu3_raw<-read_rds("simu_25_beta_1.rds")
result3_raw<-read_rds("result_25_beta_1.rds")

simu4_raw<-read_rds("simu_56_beta_1_sparse.rds")
result4_raw<-read_rds("result_56_beta_1_sparse.rds")

simu5_raw<-read_rds("simu_75_beta_0.1.rds")
result5_raw<-read_rds("result_75_beta_0.1.rds")

simu6_raw<-read_rds("simu_19_beta_0.1_sparse.rds")
result6_raw<-read_rds("result_19_beta_0.1_sparse.rds")


###simu1####
factor1_bay1<-apply(abs(result1_raw$factor1_all),2,mean)
factor2_bay1<-apply(abs(result1_raw$factor2_all),2,mean)

#FT_grs<-abs(result1_raw$ft_grs)

t<-length(factor1_bay1)

factor_simu1<-data.frame(cbind(factor1_bay1,factor2_bay1,abs(simu1_raw$FT[1,1:t]),
                               abs(simu1_raw$FT[2,1:t])))



factor_simu1 %>% dplyr::rename(f1_bay=factor1_bay1,f2_bay=factor2_bay1,
                               f1_true=V3,f2_true=V4) %>% 
  mutate(month=seq(1,length(f1_bay)))%>% 
  gather(.,type,value,f1_bay:f2_true,factor_key = T) %>% 
  mutate(factor=str_extract(type, ".+?(?=_)"),
         method=str_sub(type, 4))->factor_p1



factor_p1 %>% filter(method%in%c("bay","true"))%>% 
  ggplot(aes(x=month,y=abs(value),color=method))+
  geom_line(size = 0.4)+facet_wrap(.~factor,scales = "free_y",ncol=1)+
  scale_colour_manual(labels = c("BAY", "True"),values=cbp2[c(2,1)])+
  ylab("|Ft|")+xlab('Month')+
  ggtitle("Simulation 1")+
  theme(text=element_text(size=10), axis.text=element_text(size=7), 
        axis.title=element_text(size=7), 
        plot.title=element_text(size=7))->p_ft_s1


###simu2####
factor1_bay2<-apply(abs(result2_raw$factor1_all),2,mean)
factor2_bay2<-apply(abs(result2_raw$factor2_all),2,mean)

#FT_grs<-abs(result1_raw$ft_grs)

t<-length(factor1_bay2)

factor_simu2<-data.frame(cbind(factor1_bay2,factor2_bay2,abs(simu2_raw$FT[1,1:t]),
                               abs(simu2_raw$FT[2,1:t])))



factor_simu2 %>% dplyr::rename(f1_bay=factor1_bay2,f2_bay=factor2_bay2,
                               f1_true=V3,f2_true=V4) %>% 
  mutate(month=seq(1,length(f1_bay)))%>% 
  gather(.,type,value,f1_bay:f2_true,factor_key = T) %>% 
  mutate(factor=str_extract(type, ".+?(?=_)"),
         method=str_sub(type, 4))->factor_p2



factor_p2 %>% filter(method%in%c("bay","true"))%>% 
  ggplot(aes(x=month,y=abs(value),color=method))+
  geom_line(size = 0.4)+facet_wrap(.~factor,scales = "free_y",ncol=1)+
  scale_colour_manual(labels = c("BAY", "True"),values=cbp2[c(2,1)])+
  ylab("|Ft|")+xlab('Month')+
  ggtitle("Simulation 2")+
  theme(text=element_text(size=10), axis.text=element_text(size=7), 
        axis.title=element_text(size=7), 
        plot.title=element_text(size=7))->p_ft_s2


###simu3####
factor1_bay3<-apply(abs(result3_raw$factor1_all),2,mean)
factor2_bay3<-apply(abs(result3_raw$factor2_all),2,mean)

#FT_grs<-abs(result1_raw$ft_grs)

t<-length(factor1_bay3)

factor_simu3<-data.frame(cbind(factor1_bay3,factor2_bay3,abs(simu3_raw$FT[1,1:t]),
                               abs(simu3_raw$FT[2,1:t])))



factor_simu3 %>% dplyr::rename(f1_bay=factor1_bay3,f2_bay=factor2_bay3,
                               f1_true=V3,f2_true=V4) %>% 
  mutate(month=seq(1,length(f1_bay)))%>% 
  gather(.,type,value,f1_bay:f2_true,factor_key = T) %>% 
  mutate(factor=str_extract(type, ".+?(?=_)"),
         method=str_sub(type, 4))->factor_p3



factor_p3 %>% filter(method%in%c("bay","true"))%>% 
  ggplot(aes(x=month,y=abs(value),color=method))+
  geom_line(size = 0.4)+facet_wrap(.~factor,scales = "free_y",ncol=1)+
  scale_colour_manual(labels = c("BAY", "True"),values=cbp2[c(2,1)])+
  ylab("|Ft|")+xlab('Month')+
  ggtitle("Simulation 3")+
  theme(text=element_text(size=10), axis.text=element_text(size=7), 
        axis.title=element_text(size=7), 
        plot.title=element_text(size=7))->p_ft_s3

###simu4####
factor1_bay4<-apply(abs(result4_raw$factor1_all),2,mean)
factor2_bay4<-apply(abs(result4_raw$factor2_all),2,mean)

#FT_grs<-abs(result1_raw$ft_grs)

t<-length(factor1_bay4)

factor_simu4<-data.frame(cbind(factor1_bay4,factor2_bay4,abs(simu4_raw$FT[1,1:t]),
                               abs(simu4_raw$FT[2,1:t])))



factor_simu4 %>% dplyr::rename(f1_bay=factor1_bay4,f2_bay=factor2_bay4,
                               f1_true=V3,f2_true=V4) %>% 
  mutate(month=seq(1,length(f1_bay)))%>% 
  gather(.,type,value,f1_bay:f2_true,factor_key = T) %>% 
  mutate(factor=str_extract(type, ".+?(?=_)"),
         method=str_sub(type, 4))->factor_p4



factor_p4 %>% filter(method%in%c("bay","true"))%>% 
  ggplot(aes(x=month,y=abs(value),color=method))+
  geom_line(size = 0.4)+facet_wrap(.~factor,scales = "free_y",ncol=1)+
  scale_colour_manual(labels = c("BAY", "True"),values=cbp2[c(2,1)])+
  ylab("|Ft|")+xlab('Month')+
  ggtitle("Simulation 4")+
  theme(text=element_text(size=10), axis.text=element_text(size=7), 
        axis.title=element_text(size=7), 
        plot.title=element_text(size=7))->p_ft_s4

###simu5####
factor1_bay5<-apply(abs(result5_raw$factor1_all),2,mean)
factor2_bay5<-apply(abs(result5_raw$factor2_all),2,mean)

#FT_grs<-abs(result1_raw$ft_grs)

t<-length(factor1_bay5)

factor_simu5<-data.frame(cbind(factor1_bay5,factor2_bay5,abs(simu5_raw$FT[1,1:t]),
                               abs(simu5_raw$FT[2,1:t])))



factor_simu5 %>% dplyr::rename(f1_bay=factor1_bay5,f2_bay=factor2_bay5,
                               f1_true=V3,f2_true=V4) %>% 
  mutate(month=seq(1,length(f1_bay)))%>% 
  gather(.,type,value,f1_bay:f2_true,factor_key = T) %>% 
  mutate(factor=str_extract(type, ".+?(?=_)"),
         method=str_sub(type, 4))->factor_p5



factor_p5 %>% filter(method%in%c("bay","true"))%>% 
  ggplot(aes(x=month,y=abs(value),color=method))+
  geom_line(size = 0.4)+facet_wrap(.~factor,scales = "free_y",ncol=1)+
  scale_colour_manual(labels = c("BAY", "True"),values=cbp2[c(2,1)])+
  ylab("|Ft|")+xlab('Month')+
  ggtitle("Simulation 5")+
  theme(text=element_text(size=10), axis.text=element_text(size=7), 
        axis.title=element_text(size=7), 
        plot.title=element_text(size=7))->p_ft_s5


###simu6
factor1_bay6<-apply(abs(result6_raw$factor1_all),2,mean)
factor2_bay6<-apply(abs(result6_raw$factor2_all),2,mean)

#FT_grs<-abs(result1_raw$ft_grs)

t<-length(factor1_bay6)

factor_simu6<-data.frame(cbind(factor1_bay6,factor2_bay6,abs(simu6_raw$FT[1,1:t]),
                               abs(simu6_raw$FT[2,1:t])))



factor_simu6 %>% dplyr::rename(f1_bay=factor1_bay6,f2_bay=factor2_bay6,
                               f1_true=V3,f2_true=V4) %>% 
  mutate(month=seq(1,length(f1_bay)))%>% 
  gather(.,type,value,f1_bay:f2_true,factor_key = T) %>% 
  mutate(factor=str_extract(type, ".+?(?=_)"),
         method=str_sub(type, 4))->factor_p6



factor_p6 %>% filter(method%in%c("bay","true"))%>% 
  ggplot(aes(x=month,y=abs(value),color=method))+
  geom_line(size = 0.4)+facet_wrap(.~factor,scales = "free_y",ncol=1)+
  scale_colour_manual(labels = c("BAY", "True"),values=cbp2[c(2,1)])+
  ylab("|Ft|")+xlab('Month')+
  ggtitle("Simulation 6")+
  theme(text=element_text(size=10), axis.text=element_text(size=7), 
        axis.title=element_text(size=7), 
        plot.title=element_text(size=7))->p_ft_s6

b0a_ft<-ggarrange(p_ft_s1,p_ft_s2,p_ft_s3,
                     p_ft_s4,p_ft_s5,p_ft_s6, 
                     ncol = 2, nrow = 3,common.legend = TRUE,legend="bottom")
annotate_figure(b0a_ft,
                top = text_grob("First two factors for 6 simulations", 
                                face = "bold", size = 7))

####in-sample fit- result1_raw, beta_simu1_all,y_true1_f###

beta_simu1_all %>% filter(month==1,quarter==1,rel==1) %>% 
  select(-iteration,-month,-quarter,-rel)->beta_simu1

y_true1_f$y_true

y_insampe_1<-NULL
  for (i in 1:1000){
    ft_all_i<-rbind(result1_raw$factor1_all[i,],result1_raw$factor2_all[i,],
                  result1_raw$factor3_all[i,],result1_raw$factor4_all[i,],
                  result1_raw$factor5_all[i,],result1_raw$factor6_all[i,])
    
    beta_i<-beta_simu1[i,]
    y_i<-rep(0,31)
    
    for (k in 2:32){
      y_i[k-1]<-unlist(beta_i[1]+as.matrix(beta_i[2:7])%*%ft_all_i[,3*k]+
        as.matrix(beta_i[8:13])%*%ft_all_i[,(3*k-1)]+
        as.matrix(beta_i[14:19])%*%ft_all_i[,3*k-2]+
          beta_i[20]*y_true1_f$y_true[k-1])
    }
    
    y_i %>% unlist() %>% data.frame->y_i_f
    
    y_insampe_1<-cbind(y_insampe_1,y_i_f)
  }
   



