# load pathway data
library(data.table)
library(openxlsx)
library(tidyverse)
library(rmeta)

# prognosis taxa
load("HR.taxa.rd")
a <- list(
  Tumor = r.uu.t,
  NAT = r.uu.n
)
lst.taxa <- union(
  a$Tumor$x[a$Tumor$fdr<0.05],
  a$NAT$x[a$NAT$P<0.01]
)
a1 <- a$Tumor[a$Tumor$x %in% lst.taxa,]
a2 <- a$NAT[a$NAT$x %in% lst.taxa,]
a1$prognostic <- ifelse(a1$HR>1,"Shorter","Longer")
a1$prognostic[a1$fdr>=0.05] <- "ns"
a2$prognostic <- ifelse(a2$HR>1,"Shorter","Longer")
a2$prognostic[a2$P>=0.01] <- "ns"
all(a1$x==a2$x)
table(a1$prognostic,a2$prognostic)

out <- cbind(a1,a2)

## kk2 annotation
ann <- read.csv("../../1.data.process/kk2.v2/kraken2.db.annotation.csv")
ann1 <- ann[!duplicated(ann$genus),]
ann2 <- ann[!duplicated(ann$species),]
rownames(ann1) <- ann1$genus
rownames(ann2) <- ann2$species
ann <- rbind(ann1,ann2)
ann$p <- gsub("_\\w$","",ann$p)
ann$p <- gsub("p__","",ann$p)
rownames(ann) <- make.names(rownames(ann))
all(out[,1] %in% rownames(ann))
ann <- ann[out[,1],]
out <- cbind(ann[,c(1,2,3,4,5,6)],out)
#write.xlsx(out,"HR.taxa.xlsx")

source("prognostic_fun.R")
l.pos.t <- a1$x[a1$prognostic=="Shorter"]
l.neg.t <- a1$x[a1$prognostic=="Longer"]

a <- readRDS("../../1.data.process/kk2.v2/Data.UCAN.rmHOST.rds")
a <- a$RA$WGS
dat.reads <- read.csv("../../1.data.process/kk2.v2/reads.summary.csv",row.names = 1)
dat.reads.t <- dat.reads[dat.reads$sampleType=="Tumor",]
dat.reads.n <- dat.reads[dat.reads$sampleType=="NAT",]
pr.uu.t <- a[rownames(dat.reads.t)[dat.reads.t$cohort=="UU"],]
pr.uu.n <- a[rownames(dat.reads.n)[dat.reads.n$cohort=="UU"],]
pr.um.t <- a[rownames(dat.reads.t)[dat.reads.t$cohort=="UM"],]
pr.um.n <- a[rownames(dat.reads.n)[dat.reads.n$cohort=="UM"],]

a <- readRDS("../../1.data.process/kk2.v2/Data.ICAM.rmHOST.rds")
pr.icam <- a$RA$WGS
colnames(pr.uu.t) <- make.names(colnames(pr.uu.t))
colnames(pr.uu.n) <- make.names(colnames(pr.uu.n))
colnames(pr.um.t) <- make.names(colnames(pr.um.t))
colnames(pr.um.n) <- make.names(colnames(pr.um.n))
colnames(pr.icam) <- make.names(colnames(pr.icam))

mrs.uu.t <- Index(pr.uu.t, pos = l.pos.t, neg = l.neg.t)
mrs.um.t <- Index(pr.um.t, pos = l.pos.t, neg = l.neg.t)
mrs.icam.t <- Index(pr.icam, pos = l.pos.t, neg = l.neg.t)

l.pos.n <- a2$x[a2$prognostic=="Shorter"]
l.neg.n <- a2$x[a2$prognostic=="Longer"]
mrs.uu.n <- Index(pr.uu.n, pos = l.pos.n, neg = l.neg.n)
mrs.um.n <- Index(pr.um.n, pos = l.pos.n, neg = l.neg.n)


# summary
mrs <- list(
  UU.Tumor.WGS = mrs.uu.t,
  UU.NAT.WGS = mrs.uu.n,
  UM.Tumor.WGS = mrs.um.t,
  UM.NAT.WGS = mrs.um.n,
  ICAM.WGS = mrs.icam.t
)
saveRDS(mrs,"MRS.UCAN_ICAM.WGS.rds")

# rm stage IV patients
a <- dat.reads[dat.reads$cohort%in%c("UU","UM"),]
a <- a[a$sampleType!="Blood",]
a$Stage <- phe.ucan[a$patientID,"Tumor_Stage"]
a <- a[a$Stage !="IV",]

b <- phe.icam

mrs2 <- list(
  UU.Tumor.WGS = mrs.uu.t[phe.uu.t$Tumor_ID,,drop=F],
  UU.NAT.WGS = mrs.uu.n[phe.uu.n$Normal_ID,,drop=F],
  UM.Tumor.WGS = mrs.um.t[phe.um.t$Tumor_ID,,drop=F],
  UM.NAT.WGS = mrs.um.n[phe.um.t$Normal_ID,,drop=F],
  ICAM.WGS = mrs.icam.t[rownames(phe.icam),,drop=F]
)

sapply(mrs2,dim)
saveRDS(mrs2,"MRS.UCAN_ICAM.WGS.withoutIV.rds")







