library(pheatmap)
library(tidyverse)
library(ggpubr)
library(survival)
library(survminer)
library(openxlsx)
library(forestploter)


# load data
# profiles (UU-WGS)
# load data
load("Data.rd")


prognosisFun <- function(datx,f,adj){
  x1 <- "Enterobacteriaceae"
  r1 <- HR(datx[which(datx$CMS=="CMS1"),],f,x=x1,adj=adj)
  r2 <- HR(datx[which(datx$CMS=="CMS2"),],f,x=x1,adj=adj)
  r3 <- HR(datx[which(datx$CMS=="CMS3"),],f,x=x1,adj=adj)
  r4 <- HR(datx[which(datx$CMS=="CMS4"),],f,x=x1,adj=adj)
  
  hr.cms <- rbind(r1,r2,r3,r4)
  hr.cms$x <- (c("CMS1","CMS2","CMS3","CMS4"))
  hr.cms
}
#-------------------------------------------------------------------------------
HR <- function(data,f,x,adj=NULL,outputAll=FALSE){
  data <- data[!is.na(data[,x]),]
  formula <- paste0(f,x)
  if(!is.null(adj)){
    lab.adj <- paste(adj,collapse = "+")
    formula <- paste0(formula,"+",lab.adj)
  }
  formula <- as.formula(formula)
  print(formula)
  #x <- unlist(strsplit(as.character(formula)[3]," "))[1]
  fit <- coxph(formula,data=data)
  a <- summary(fit)
  # output
  out <- cbind(
    a$coefficients[,c(1,4,5),drop=F],
    a$conf.int[,c(1,3,4),drop=F]
  )
  colnames(out) <- c("coef","z","P","HR","lower.95","upper.95")
  out <- out %>%
    as.data.frame() %>%
    rownames_to_column(var = "x") %>%
    mutate(
      HR2 = paste0(round(HR,2)," (",round(lower.95,2)," to ",round(upper.95,2),")")
    )
  if(!outputAll){
    tmp <- lapply(x,function(i,j) grep(i,j),j=out$x)
    tmp <- unlist(tmp)
    out <- out[tmp,]
    if(is.factor(data[,x])){
      n <- table(data[,x])
      out$N1 <- n[1]
      out$N2 <- n[2]
      out$N <- paste0(out$N2," vs. ",out$N1)
    }else{
      out$N1 <- NA; out$N2 <- NA
      out$N <- nrow(data)
    }
  }
  out
}
HR_Forest <- function(dat,title=NULL){
  dat$` ` <- paste(rep(" ", 20), collapse = " ")
  colnames(dat)[colnames(dat)=="HR2"] <- "HR (95% CI)"
  dat$P <- round(dat$P,4)
  forest(
    data = dat[,c(1,ncol(dat)-1,ncol(dat),8,4),drop=F],
    est=dat$HR,
    lower = dat$lower.95,
    upper = dat$upper.95,
    ci_column = 3,
    ref_line = 1,
    size=0.6,
    x_trans = "log",
    xlab = c("Hazard ratio"),
    xlim = c(0.5,5),
    title = title
  )
}


f1 <- "Surv(OS_time,OS_status)~"
f2 <- "Surv(RFS_time,RFS_status)~"
f3 <- "Surv(PFS_time,PFS_status)~"

a <- dat.um.t
b <- dat.icam.t2
colnames(b)[3:4] <- colnames(a)[3:4]
a <- a[,colnames(b)]
a$OS_status <- a$OS_status -1 
a$RFS_status <- a$RFS_status -1
a$Group <- "UM"
b$Group <- "ICAM"
dat2 <- rbind(a,b)

dat1 <- dat.uu.t
dat3 <- dat.icam.t1

dat1$Treatment_After_Surgery <- phe.ucan[rownames(dat1),"Treatment_After_Surgery"]


dat1$CMS <- factor(dat1$CMS,levels = paste0("CMS",1:4),labels = c("CMS-other","CMS2","CMS-other","CMS-other"))
dat2$CMS <- factor(dat2$CMS,levels = paste0("CMS",1:4),labels = c("CMS-other","CMS2","CMS-other","CMS-other"))
dat3$CMS <- factor(dat3$CMS,levels = paste0("CMS",1:4),labels = c("CMS-other","CMS2","CMS-other","CMS-other"))


prognosisFun <- function(datx,f,adj){
  x1 <- "Enterobacteriaceae"
  r1 <- HR(datx[which(datx$CMS=="CMS-other"),],f,x=x1,adj=adj)
  r2 <- HR(datx[which(datx$CMS=="CMS2"),],f,x=x1,adj=adj)
  
  hr.cms <- rbind(r1,r2)
  hr.cms$x <- (c("CMS-other","CMS2"))
  hr.cms
}


f <- f1
r1 <- prognosisFun(dat1,f,c(adj.ucan))
r2 <- prognosisFun(dat2,f,c(adj.icam2,"Group"))
r3 <- prognosisFun(dat3,f,c(adj.icam2))

p1 = HR_Forest(r1)
p2 = HR_Forest(r2)
p3 = HR_Forest(r3)


pdf("Fig4f.forest.combined.V2.pdf",width = 10,height = 6)
gridExtra::grid.arrange(p1,p2,p3,nrow = 3)
dev.off()
