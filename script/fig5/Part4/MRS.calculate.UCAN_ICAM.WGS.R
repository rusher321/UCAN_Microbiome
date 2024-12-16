# load pathway data
library(data.table)
library(openxlsx)
library(tidyverse)
library(rmeta)

# profile
pr.taxa <- readRDS("../../../2.SourceData/00.SourceData/Kraken2.UCAN.WGS.Splited.rds")
pr.uu.t <- pr.taxa$UU.Tumor
pr.um.t <- pr.taxa$UM.Tumor
pr.uu.n <- pr.taxa$UU.NAT
pr.um.n <- pr.taxa$UM.NAT
pr.icam.t <- readRDS("../../../2.SourceData/00.SourceData/Kraken2.ICAM.WGS.rds")

phe.ucan <- read.csv("../../../2.SourceData/01.phe/20230727-CRC-SW.20230704-v5.clinical.csv",row.names = 1)
a <- phe.ucan[rownames(pr.uu.n),]
pr.uu.n <- pr.uu.n[a$Normal_DNA_Control_Type=="Tumor_Adjacent_Tissues",,drop=F]
SID.rmStage4 <- rownames(phe.ucan)[phe.ucan$Tumor_Stage!="IV"]


## preprocess
# UCAN 
rownames(pr.uu.t) <- paste0(rownames(pr.uu.t),"-T")
rownames(pr.uu.n) <- paste0(rownames(pr.uu.n),"-N")
rownames(pr.um.t) <- paste0(rownames(pr.um.t),"-T")
rownames(pr.um.n) <- paste0(rownames(pr.um.n),"-N")

# combine tumor and NAT to re-normalization
pr.uu <- rbind(pr.uu.t,pr.uu.n)
pr.um <- rbind(pr.um.t,pr.um.n)
filterPer = function (x, row, percent, include = T) 
{
  if (include) {
    index <- apply(x, row, function(x) {
      (sum(x != 0)/length(x)) > percent
    })
  }
  else {
    index <- apply(x, row, function(x) {
      (sum(x != 0)/length(x)) < percent
    })
  }
  if (row == 1) {
    out <- x[index, ]
  }
  else {
    out <- x[, index]
  }
  return(out)
}

renorm <- function(x,occ=0.2){x <- filterPer(x, 2, occ);x/rowSums(x)}
a <- pr.uu[,grep("g__",colnames(pr.uu))]
a <- renorm(a,occ = 0.5)
b <- pr.uu[,grep("s__",colnames(pr.uu))]
b <- renorm(b)
pr.uu <- cbind(a,b)
a <- pr.um[,grep("g__",colnames(pr.um))]
a <- renorm(a,occ = 0.5)
b <- pr.um[,grep("s__",colnames(pr.um))]
b <- renorm(b)
pr.um <- cbind(a,b)
pr.uu.t <- pr.uu[grep("-T",rownames(pr.uu)),]
pr.um.t <- pr.um[grep("-T",rownames(pr.um)),]
rownames(pr.uu.t) <- gsub("-T","",rownames(pr.uu.t))
rownames(pr.um.t) <- gsub("-T","",rownames(pr.um.t))
pr.uu.n <- pr.uu[grep("-N",rownames(pr.uu)),]
pr.um.n <- pr.um[grep("-N",rownames(pr.um)),]
rownames(pr.uu.n) <- gsub("-N","",rownames(pr.uu.n))
rownames(pr.um.n) <- gsub("-N","",rownames(pr.um.n))

# remove stage IV
pr.uu.t <- pr.uu.t[rownames(pr.uu.t) %in% SID.rmStage4,]
pr.uu.n <- pr.uu.n[rownames(pr.uu.n) %in% SID.rmStage4,]
pr.um.t <- pr.um.t[rownames(pr.um.t) %in% SID.rmStage4,]
pr.um.n <- pr.um.n[rownames(pr.um.n) %in% SID.rmStage4,]


# ICAM-WGS
# remove stage IV
phe.icam <- readRDS("../../../2.SourceData/00.SourceData/phenotype.ICAM.rds")
phe.icam <- phe.icam[rownames(pr.icam.t),]
pr.icam.t <- pr.icam.t[phe.icam$AJCC_path_stage!=4,]

# prognosis taxa
lst.taxa <- read.xlsx("HZ-List.Prognosis.Taxa.V2.xlsx",sheet = 3,rowNames = T)


source("prognostic_fun.R")
l.pos.t <- rownames(lst.taxa)[lst.taxa$tumor==1 & lst.taxa$Group=="shorter OS"]
l.neg.t <- rownames(lst.taxa)[lst.taxa$tumor==1 & lst.taxa$Group=="longer OS"]
mrs.uu.t <- Index(pr.uu.t, pos = l.pos.t, neg = l.neg.t)
mrs.um.t <- Index(pr.um.t, pos = l.pos.t, neg = l.neg.t)
mrs.icam.t <- Index(pr.icam.t, pos = l.pos.t, neg = l.neg.t)

l.pos.n <- rownames(lst.taxa)[lst.taxa$normal==1 & lst.taxa$Group=="shorter OS"]
l.neg.n <- rownames(lst.taxa)[lst.taxa$normal==1 & lst.taxa$Group=="longer OS"]
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


# debug
a <- readRDS("BACKUP/Index_tumor.rds")
b <- readRDS("BACKUP/Index_normal.rds")
c <- readRDS("BACKUP/icam_phe_rm4.rd")
ix.uu.t <- a$uu
ix.um.t <- a$um
ix.uu.n <- b$uu
ix.um.n <- b$um
ix.icam.t <- c

rownames(ix.uu.t) <- gsub(".*-(U.*)-.*","\\1",rownames(ix.uu.t))
rownames(ix.um.t) <- gsub(".*-(U.*)-.*","\\1",rownames(ix.um.t))
rownames(ix.uu.n) <- gsub(".*-(U.*)-.*","\\1",rownames(ix.uu.n))
rownames(ix.um.n) <- gsub(".*-(U.*)-.*","\\1",rownames(ix.um.n))
rownames(ix.icam.t) <- gsub("SER-SILU-CC-P0","",ix.icam.t$Patient_ID)

ix.uu.t <- ix.uu.t[rownames(mrs.uu.t),]
ix.uu.n <- ix.uu.n[rownames(mrs.uu.n),]
ix.um.t <- ix.um.t[rownames(mrs.um.t),]
ix.um.n <- ix.um.n[rownames(mrs.um.n),]
ix.icam.t <- ix.icam.t[rownames(mrs.icam.t),]

all(ix.uu.t$Index1==mrs.uu.t$Index)
all(ix.uu.n$Index1==mrs.uu.n$Index)
all(ix.um.t$Index1==mrs.um.t$Index)
all(ix.um.n$Index1==mrs.um.n$Index)


ix <- data.frame(old=ix.icam.t[,"MRS_index"],new=mrs.icam.t[rownames(ix.icam.t),])
ix <- data.frame(old=ix.uu.n[,"Index1"],new=mrs.uu.n$Index)

