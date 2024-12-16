library(dynpred)
cindexFun <- function(dat,x,add,f){
  x <- union(x,add)
  x <- paste(x,collapse = "+")
  f <- paste0(f,x)
  print(f)
  f <- as.formula(f)
  ix = cindex(f, data=dat)$cindex
  print(ix)
  ix
}
cindexFun_rfs <- function(dat,x,add){
  x <- union(x,add)
  x <- paste(x,collapse = "+")
  f <- paste0("Surv(PFS_Time, PFS_Status)~",x)
  print(f)
  f <- as.formula(f)
  dynpred::cindex(f, data=dat)$cindex
}
cindexFun2 <- function(dat,x,f){
  print(x)
  a <- cindexFun(dat,x=x,add=NULL,f)
  b <- cindexFun(dat,x=x,add="Index",f)
  c(a,b)
}


mrs <- readRDS("MRS.UCAN_ICAM.WGS.rds")
mrs$UU.Tumor.WGS <- mrs$UU.Tumor.WGS[rownames(mrs$UU.NAT.WGS),,drop=F]
mrs$UM.Tumor.WGS <- mrs$UM.Tumor.WGS[rownames(mrs$UM.NAT.WGS),,drop=F]


# data
phe.ucan <- read.csv("../../../2.SourceData/01.phe/20230727-CRC-SW.20230704-v5.clinical.csv",row.names = 1)
phe.icam <- readRDS("../../../2.SourceData/00.SourceData/phenotype.ICAM.rds")
phe.ucan$Tumor_Site <- ifelse(phe.ucan$Tumor_Site=="Right_Colon","Right","Left")
phe.ucan$Age <- phe.ucan$Age_At_Diagnosis
phe.icam$AJCC_path_stage <- as.character(phe.icam$AJCC_path_stage)

x <- c("Age","Tumor_Site","Tumor_Stage","Gender","Hyper_Mutated_Status","MSI_Status","CMS_Subtype")


# confounders
adj.uu <- phe.ucan[rownames(mrs$UU.Tumor.WGS),x]
adj.um <- phe.ucan[rownames(mrs$UM.Tumor.WGS),x]
adj.icam <- phe.icam[rownames(mrs$ICAM.WGS),c("Age","Location","AJCC_path_stage","Gender","HM","MSI","CMS")]
colnames(adj.icam) <- c("Age","Tumor_Site","Tumor_Stage","Gender","Hyper_Mutated_Status","MSI_Status","CMS_Subtype")
adj.icam$Tumor_Site <- factor(adj.icam$Tumor_Site)
adj.icam$Tumor_Stage <- factor(adj.icam$Tumor_Stage)
adj.icam$CMS_Subtype[adj.icam$CMS_Subtype=="mixed"] <- NA


dat.uu <- phe.ucan[rownames(adj.uu),c("OS_5Years_with_censored_update","OS_Status_5Years_with_censored_update","RFS_Status_5Years_with_censored","RFS_5Years_with_censored")]
dat.um <- phe.ucan[rownames(adj.um),c("OS_5Years_with_censored_update","OS_Status_5Years_with_censored_update","RFS_Status_5Years_with_censored","RFS_5Years_with_censored")]
colnames(dat.uu)[1:4] <- colnames(dat.um)[1:4] <- c("OS_time","OS_status","RFS_status","RFS_time")
dat.icam <- phe.icam[rownames(adj.icam),c("OS_status","OS_time","PFS_status","PFS_time")]


d1 <- cbind(dat.uu,adj.uu[rownames(dat.uu),],mrs$UU.Tumor.WGS[rownames(dat.uu),,drop=F])
d2 <- cbind(dat.um,adj.um[rownames(dat.um),],mrs$UM.Tumor.WGS[rownames(dat.um),,drop=F])
d3 <- cbind(dat.icam,adj.icam[rownames(dat.icam),],mrs$ICAM.WGS[rownames(dat.icam),,drop=F])
d4 <- cbind(dat.uu,adj.uu[rownames(dat.uu),],mrs$UU.NAT.WGS[rownames(dat.uu),,drop=F])
d5 <- cbind(dat.um,adj.um[rownames(dat.um),],mrs$UM.NAT.WGS[rownames(dat.um),,drop=F])


Cindex.addon <- function(d1,d2,d3,d4,d5,x){
  f1 <- "Surv(OS_time,OS_status)~"
  f2 <- "Surv(RFS_time,RFS_status)~"
  f3 <- "Surv(PFS_time,PFS_status)~"

  res1 <- list(
    UU.T = cindexFun2(d1,x=x,f=f1),
    UM.T = cindexFun2(d2,x=x,f=f1),
    ICAM.T = cindexFun2(d3,x=x,f=f1),
    UU.N = cindexFun2(d4,x=x,f=f1),
    UM.N = cindexFun2(d5,x=x,f=f1)
  )
  res2 <- list(
    UU.T = cindexFun2(d1,x=x,f=f2),
    UM.T = cindexFun2(d2,x=x,f=f2),
    ICAM.T = cindexFun2(d3,x=x,f=f3),
    UU.N = cindexFun2(d4,x=x,f=f2),
    UM.N = cindexFun2(d5,x=x,f=f2)
  )
  
  res1 <- do.call(rbind,res1)
  res2 <- do.call(rbind,res2)
  res1 <- as.data.frame(res1) ; res2 <- as.data.frame(res2)
  res1$Group <- "OS"; res2$Group <- "RFS"
  res <- rbind(res1,res2)
  res$Cohort <- rep(c("UU","UM","AC-ICAM","UU","UM"),2)
  res$MRS <- rep(c("MRS-T","MRS-T","MRS-T","MRS-N","MRS-N"),2)
  colnames(res)[1:2] <- c("Ref","Add.MRS")
  
  a <- res %>%
    gather(x,Cindex,-Cohort,-MRS,-Group)
  b <- a[a$x=="Ref",]
  a$x[a$x=="Add.MRS"] <- as.character(a$MRS[a$x=="Add.MRS"])
  b$x <- b$MRS
  a <- rbind(a,b)
  a$x <- factor(a$x,levels = c("MRS-N","MRS-T","Ref"))
  a$Cohort <- factor(a$Cohort,levels = c("UU","UM","AC-ICAM"))
  a$g <- paste0(a$Group," (",a$Cohort,")")
  a$g <- factor(a$g,levels = c("OS (UU)","OS (UM)","OS (AC-ICAM)","RFS (UU)","RFS (UM)","RFS (AC-ICAM)"))
  a$g2 <- ifelse(a$Cindex %in% a$Cindex[a$x=="Ref"],"g1","g2")
  print(a)
  
  ggplot(a,aes(x,Cindex))+
    geom_point(size=3,aes(color=g2))+
    coord_flip()+
    facet_wrap(g~.,scales = "free",nrow = 6,strip.position="left")+
    geom_line(aes(group=x),arrow=arrow(ends = "first",length = unit(0.1,"inches")))+
    scale_color_manual(values = c("#16A085","#E74C3C"))+
    #scale_y_continuous(limits = c(0.63,0.88))+
    theme_bw()+
    theme(
      axis.text = element_text(color = 1),
      axis.title = element_text(color = 1),
      panel.grid = element_blank(),
      legend.position = ""
    )
}


res.all <- Cindex.addon(
  d1,
  d2,
  d3,
  d4,
  d5,
  x = x
)+labs(title = "All")

res.left <- Cindex.addon(
  d1[d1$Tumor_Site=="Left",],
  d2[d2$Tumor_Site=="Left",],
  d3[d3$Tumor_Site=="Left",],
  d4[d4$Tumor_Site=="Left",],
  d5[d5$Tumor_Site=="Left",],
  x = x[x!="Tumor_Site"]
)+labs(title = "Left")


res.right <- Cindex.addon(
  d1[d1$Tumor_Site=="Right",],
  d2[d2$Tumor_Site=="Right",],
  d3[d3$Tumor_Site=="Right",],
  d4[d4$Tumor_Site=="Right",],
  d5[d5$Tumor_Site=="Right",],
  x = x[x!="Tumor_Site"]
)+labs(title = "Right")


ggarrange(res.all,res.left,res.right,nrow = 1,ncol = 3,align = "hv")

ggsave("point.Cindex.pdf",width = 6,height = 7)






