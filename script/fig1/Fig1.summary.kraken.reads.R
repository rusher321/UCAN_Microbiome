library(openxlsx)
library(tidyverse)
library(ggvenn)
library(ggpubr)
library(data.table)



# load raw data
load("RawData.rd")

#-------------------------------------------------------------------------------
# FigS1.a boxplot

Boxplot_Reads <- function(d,x,y){
  ggboxplot(d,x=x,y=y,fill=x,ylab = "Reads Number",xlab = "")+
    scale_y_log10(limit=c(10,5e9))+
    scale_fill_manual(values = c("#9EBCE3","#F1D1BB","#D5B226","#862320"))+
    theme(legend.position = "")
}

# WGS
p1 <- Boxplot_Reads(dat.wgs,x="Group1",y="nonHuman")
p2 <- Boxplot_Reads(dat.wgs,x="Group1",y="Microbiota")
p3 <- Boxplot_Reads(dat.wgs,x="Group1",y="nonContam")
p4 <- Boxplot_Reads(dat.wgs,x="Group1",y="TissueResident")

# WTS
p5 <- Boxplot_Reads(dat.wts,x="Group1",y="nonHuman")
p6 <- Boxplot_Reads(dat.wts,x="Group1",y="Microbiota")
p7 <- Boxplot_Reads(dat.wts,x="Group1",y="nonContam")
p8 <- Boxplot_Reads(dat.wts,x="Group1",y="TissueResident")

ggarrange(p1,p5,p2,p6,p3,p7,p4,p8,nrow = 4,ncol = 2,align = "hv",widths = c(2.8,2))
ggsave("FigS1a.boxplot.pipeline.reads.pdf",width = 4,height = 10)



#-------------------------------------------------------------------------------
# FigS1.a barplot

# Ratio of discard and remaining
fun3 <- function(d){
  d$Ratio.prokary <- d$Microbiota/d$nonHuman
  d$Ratio.nonContam <- d$nonContam/d$AllGenus
  d$Ratio.tissueResident <- d$TissueResident/d$AllGenus
  d
}

dat.wgs <- fun3(dat.wgs)
dat.wts <- fun3(dat.wts)

# function for barplot
barplot_reads <- function(dat,y){
  dat <- dat[,c("Group1",y)] %>%
    group_by(Group1) %>%
    summarise_all(mean)
  colnames(dat) <- c("g","remain")
  a <- dat
  a$lab <- paste0(round(a$remain*100),"%")
  dat <- dat %>%
    mutate(discard=1-remain) %>%
    gather(Group,Ratio,-g) %>%
    mutate(
      Group = factor(Group,levels=c("discard","remain"),labels=c("Propotion reads discard","Propotion reads remaining"))
    )
  print(dat)
  ggbarplot(dat,x="g",y="Ratio",fill="Group")+
    geom_text(data=a,aes(x=g,y=0.2,label=lab),size=5)+
    scale_y_continuous(expand = expansion(mult = c(0,0)))+
    scale_fill_manual(values = c("gray50","#ed7d31"))+
    xlab("")+ylab("Propotion of reads")+
    guides(fill=guide_legend(nrow=2))+
    theme(
      axis.text = element_text(size = 13),
      axis.title = element_text(size = 14),
      legend.text = element_text(size = 13),
      legend.title = element_blank()
    )
}

p9  <- barplot_reads(dat.wgs,y="Ratio.prokary")
p10 <- barplot_reads(dat.wgs,y="Ratio.nonContam")
p11 <- barplot_reads(dat.wgs,y="Ratio.tissueResident")
p12 <- barplot_reads(dat.wts,y="Ratio.prokary")
p13 <- barplot_reads(dat.wts,y="Ratio.nonContam")
p14 <- barplot_reads(dat.wts,y="Ratio.tissueResident")

ggarrange(p1,NULL,p2,p9,p3,p10,p4,p11,nrow = 4,ncol = 2,align = "hv",common.legend = T,legend = "top")
ggsave("summary.reads.WGS.pdf",width = 6,height = 13)
ggarrange(p5,NULL,p6,p12,p7,p13,p8,p14,nrow = 4,ncol = 2,align = "hv",common.legend = T,legend = "top")
ggsave("summary.reads.RNA.pdf",width = 4,height = 13)




###---------------------------------------------------------------------------------
# FigS1b
# WGS vs. RNA in 775 UU patients
n.wgs <- dat.wgs[dat.wgs$cohort=="UU" & dat.wgs$sampleType=="Tumor",]
all(rownames(n.wgs) %in% rownames(dat.wts))
n.wts <- dat.wts[rownames(n.wgs),]

n.wgs$IDs <- rownames(n.wgs)
n.wts$IDs <- rownames(n.wts)
n.wgs$Group <- "WGS"
n.wts$Group <- "WTS"
dat <- rbind(n.wgs[,colnames(n.wts)],n.wts)
dat <- dat[order(dat$IDs),]

plt <- function(d,x,y,ylab){
  ggboxplot(d,x=x,y=y,fill=x)+
    ylab(ylab)+xlab("")+
    stat_compare_means(comparisons = list(c(1,2)))+
    scale_y_log10(limits=c(1,500000000))+
    scale_fill_manual(values = c("#FF7F50","#87CEEB"))+
    theme(
      #axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
      axis.text = element_text(size = 13),
      axis.title = element_text(size = 14),
      legend.position = ""
    )
}

p1 <- plt(dat,x="Group",y="Microbiota",ylab="# Kraken2-mapped prokaryotic reads")
p2 <- plt(dat,x="Group",y="TissueResident",ylab="# Reads annotated to\ntissue-resident species")


# species number
a <- readRDS("../../../1.data.process/kk2.v2/Data.UCAN.rmHOST.rds")
pr1 <- a$counts$WGS[rownames(n.wgs),]
pr2 <- a$counts$WTS[rownames(n.wgs),]
pr1 <- pr1[,grep("s__",colnames(pr1))]
pr2 <- pr2[,grep("s__",colnames(pr2))]

nSpecies <- data.frame(
  WGS = apply(pr1,1,function(x)  sum(x>=100)),
  WTS = apply(pr2,1,function(x)  sum(x>=100))
)

nSpecies <- nSpecies %>%
  gather(Group,n)
nSpecies$Group <- factor(nSpecies$Group,levels = c("WGS","WTS"))

p3 <-
  ggboxplot(nSpecies,x="Group",y="n",fill="Group",ylab = "Identified tissue-resident species")+
  stat_compare_means(comparisons = list(c(1,2)))+
  scale_fill_manual(values = c("#FF7F50","#87CEEB"))+
  theme(legend.position = "")

ggarrange(p1,p2,p3,nrow = 1,ncol = 3,align = "hv")
ggsave("FigS1cd.boxplot.WGS_WTS.reads.pdf",width = 7,height = 4)

wilcox.test(nSpecies$n~nSpecies$Group,paired=T)$p.value
#1.665334e-89

n <- n[order(n$SampleID),]
wilcox.test(n$Reads.KrakenMapped~n$Group,paired=T)$p.value
wilcox.test(n$Reads.Common~n$Group,paired=T)$p.value





#-------------------------------------------------------------------------------
#                               Fig 1
#-------------------------------------------------------------------------------
# Fig1b paired tumor and blood samples

dat <- dat.wgs[dat.wgs$Group1 %in% c("UU","Blood"),]
dat$PatientID <- gsub(".*-(U.*)-.*","\\1",rownames(dat))
dat <- dat[dat$PatientID %in% dat$PatientID[dat$Group1=="Blood"],]
dat <- dat[order(dat$PatientID),]
dat$sampleType <- factor(dat$sampleType,levels = c("Tumor","Blood"))
dat <- dat[order(dat$sampleType),]

ggplot(dat,aes(sampleType,y=Microbiota,fill=sampleType))+
  geom_line(aes(group=PatientID),alpha=0.1)+
  geom_point(size=3,alpha=0.7,aes(fill=sampleType),shape=21,color="white")+
  geom_boxplot(outlier.shape=NA,width=0.75)+
  ylab("Kraken2-mapped prokaryotic reads")+
  scale_y_log10(limits=c(10,5e9))+stat_compare_means(comparisons = list(c(1,2)))+
  scale_fill_manual(values = c("#9EBCE3","#862320"))+
  scale_color_manual(values = c("#9EBCE3","#862320"))+
  theme_pubr()+
  theme(legend.position = "")
ggsave("Fig1b.boxplot.reads.pdf",width = 3,height = 4)

wilcox.test(dat$Microbiota~dat$sampleType,paired=T)$p.value
# 9.146918e-33

#------------------------------------------------------------------------------
# Fig1c pie plot
ann <- read.csv("../../../1.data.process/kk2.v2/kraken2.db.annotation.csv",header = T)
ann <- ann[!duplicated(ann$genus),]
ann <- ann[ann$d != "d__Eukaryota",]
ann$p <- gsub("_\\w$","",ann$p)
ann$p <- gsub("p__","",ann$p)

dat$Ratio.contam <- 1-dat$Ratio.nonContam

dat.ratio <- dat[,c("Ratio.contam","Ratio.nonContam")]
colnames(dat.ratio) <- c("Contaminants","NonContaminants")

dat.tumor <-  data.frame(mean=apply(dat.ratio[dat$sampleType=="Tumor",],2,mean))
dat.tumor$x <- rownames(dat.tumor)
dat.tumor$ypos <- cumsum(dat.tumor$mean)- 0.5*dat.tumor$mean
dat.tumor$lab <- paste0(dat.tumor$x,"\n",round(dat.tumor$mean*100,2),"%")
dat.tumor$x <- factor(dat.tumor$x,levels = c("Contaminants","NonContaminants"))


dat.blood <-  data.frame(mean=apply(dat.ratio[dat$sampleType=="Blood",],2,mean))
dat.blood$x <- rownames(dat.blood)
dat.blood$ypos <- cumsum(dat.blood$mean)- 0.5*dat.blood$mean
dat.blood$lab <- paste0(dat.blood$x,"\n",round(dat.blood$mean*100,2),"%")
dat.blood$x <- factor(dat.blood$x,levels = c("Contaminants","NonContaminants"))


p1 <- 
  ggplot(dat.blood, aes(x="", y=mean, fill=x)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=3.7) +
  theme_void() + 
  theme(legend.position="none") +
  geom_text(aes(y = ypos, label = lab), color = "white", size=3) +
  scale_fill_manual(values = c(Contaminants="#AB8944",NonContaminants="#3778AD"))+
  labs(title = "Blood")
p1
p2 <- 
  ggplot(dat.tumor, aes(x="", y=mean, fill=x)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=1.7) +
  theme_void() + 
  theme(legend.position="none") +
  geom_text(aes(y = ypos, label = lab), color = "white", size=3) +
  scale_fill_manual(values = c(Contaminants="#AB8944",NonContaminants="#3778AD"))+
  labs(title = "Tissue")
p2

ggarrange(p1,p2,align = "hv",nrow = 2,ncol = 1)
ggsave(filename = "Fig1c.1.pie.pdf",width = 3,height = 6)



# phylum
fun4 <- function(p,l,ann){
  p <- p[,l]
  p <- as.data.frame(t(p))
  all(rownames(p) %in% ann$genus)
  rownames(ann) <- ann$genus
  ann <- ann[rownames(p),]
  p$phylum <- ann$p
  p <- p %>%
    group_by(phylum) %>%
    summarise_all(sum) %>%
    mutate_if(is.numeric,function(x) x/sum(x,na.rm = T)) %>%
    column_to_rownames(var = "phylum")
  a <- data.frame(mean=apply(p,1,function(x) mean(x,na.rm = T)))
}

dat.blood <- fun4(pr.wgs[rownames(dat)[dat$sampleType=="Blood"],],rownames(lst.contam),ann)
dat.tumor <- fun4(pr.wgs[rownames(dat)[dat$sampleType=="Tumor"],],rownames(lst.contam),ann)
dat.tumor <- dat.tumor[rownames(dat.blood),,drop=F]
dat <- cbind(dat.tumor,dat.blood)
colnames(dat) <- c("UU-Tissue","UU-Blood")
dat <- dat[order(dat$`UU-Blood`,decreasing = T),]
dat$Phylum <- factor(rownames(dat),levels = rownames(dat))
dat <- dat %>% gather(Group,abun,-Phylum)
dat$Group <- factor(dat$Group,levels = c("UU-Blood","UU-Tissue"))

p1 <- 
  ggplot(dat[dat$Group=="UU-Blood",],aes(x="",y=abun,fill=Phylum))+
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  xlab("")+ylab("Mean relative abundance")+
  scale_fill_brewer(palette="Set1")+
  labs(title = "Blood")
p2 <- 
  ggplot(dat[dat$Group=="UU-Tissue",],aes(x="",y=abun,fill=Phylum))+
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  xlab("")+ylab("Mean relative abundance")+
  scale_fill_brewer(palette="Set1")+
  labs(title = "Tissue")
ggarrange(p1,p2,nrow = 1,ncol = 2,align = "hv")
ggsave("Fig1c.2.pie.pdf",width = 10,height = 4)




#------------------------------------------------------------------------------
# Fig1d

# NAT vs. Tumor vs. Blood (UU WGS)
n <- dat.wgs[dat.wgs$cohort=="UU",]
n$Group <- factor(n$sampleType,levels = c("NAT","Tumor","Blood"))
ggplot(n,aes(Group,y=TissueResident,fill=Group))+
  geom_jitter(size=3,alpha=0.6,aes(fill=Group),shape=21,color="white",width = 0.2)+
  stat_compare_means(method = "wilcox.test",comparisons = list(c(1,2),c(2,3),c(1,3)))+
  geom_boxplot(outlier.shape=NA,width=0.5)+
  scale_y_log10(limits=c(10,5e9))+stat_compare_means(comparisons = list(c(1,2)))+
  ylab("Reads annotated to tissue resident taxa")+
  scale_fill_manual(values = c("#399938","#9EBCE3","#862320"))+
  scale_color_manual(values = c("#399938","#9EBCE3","#862320"))+
  theme_pubr()+
  theme(legend.position = "")
ggsave("Fig1d.boxplot.reads.pdf",width = 2.5,height = 4)

a <- n[n$Group %in% c("Tumor","Blood"),]
a$Group <- droplevels(a$Group)
wilcox.test(a$TissueResident~a$Group)$p.value
# 2.087532e-168
a <- n[n$Group %in% c("NAT","Blood"),]
a$Group <- droplevels(a$Group)
wilcox.test(a$TissueResident~a$Group)$p.value
# 5.248637e-103


# coverage
dat.cov <- read.csv("../../../1.data.process/bedtools.cov.summary.UCAN.allReads_db2.csv",row.names = 1,header = T,check.names = F)
genome2species <- read.table("../../../1.data.process/Selected.genomes.db2.txt",row.names = 1,header = T,sep = "\t")
genome2species <- genome2species[,"species",drop=F]
colnames(genome2species) <- "Species"
genome2species <- genome2species[rownames(dat.cov),,drop=F]
rownames(genome2species) <- rownames(dat.cov)
genome2species$Species[1] <- "s__Fusobacterium nucleatum subspecies animalis Clade1"
all(rownames(genome2species)==rownames(dat.cov))
rownames(dat.cov) <- genome2species$Species
dat.cov <- as.data.frame(t(dat.cov))
rownames(dat.cov) <- gsub("B$","",rownames(dat.cov))
dat.cov <- dat.cov[,colnames(pr1)]
dat.cov <- dat.cov[rownames(dat.wgs)[dat.wgs$cohort=="UU"],]
dat.cov$Group <- dat.wgs[rownames(dat.cov),"sampleType"]

# average coverage
a <- dat.cov %>%
  group_by(Group) %>%
  summarise_all(mean) %>%
  gather(species,coverage,-Group)
a$Group <- factor(a$Group,levels = c("NAT","Tumor","Blood"))

ggplot(a,aes(Group,y=coverage,fill=Group))+
  geom_jitter(size=3,alpha=0.6,aes(fill=Group),shape=21,color="white",width = 0.2)+
  stat_compare_means(method = "wilcox.test",comparisons = list(c(1,2),c(2,3),c(1,3)))+
  geom_boxplot(outlier.shape=NA,width=0.5)+
  #scale_y_log10(limits=c(10,5e9))+
  stat_compare_means(comparisons = list(c(1,2)))+
  ylab("Reads annotated to tissue resident taxa")+
  scale_fill_manual(values = c("#399938","#9EBCE3","#862320"))+
  scale_color_manual(values = c("#399938","#9EBCE3","#862320"))+
  theme_pubr()+
  theme(legend.position = "")
ggsave("Fig1d.boxplot.coverage.pdf",width = 2.5,height = 4)

b <- a[a$Group %in% c("Tumor","Blood"),]
b$Group <- droplevels(b$Group)
wilcox.test(b$coverage~b$Group)$p.value
# 1.62675e-119
b <- a[a$Group %in% c("NAT","Blood"),]
b$Group <- droplevels(b$Group)
wilcox.test(b$coverage~b$Group)$p.value
# 1.389106e-119
b <- a[a$Group %in% c("NAT","Tumor"),]
b$Group <- droplevels(b$Group)
wilcox.test(b$coverage~b$Group)$p.value
# 3.815928e-09

#------------------------------------------------------------------------------
# Fig1e

# average coverage
a <- dat.cov
a$Group <- factor(a$Group,levels = c("NAT","Tumor","Blood"),labels = c("Tissue","Tissue","Blood"))
meanCoverage <- a %>%
  group_by(Group) %>%
  summarise_all(mean) %>%
  gather(species,coverage,-Group)
# average read number
b <- readRDS("../../../1.data.process/kk2.v2/Data.UCAN.rmHOST.rds")
b <- b$counts$WGS[rownames(dat.cov),]
b <- b[,grep("s__",colnames(b))]
b$Group <- a$Group
meanReads <- b %>%
  group_by(Group) %>%
  summarise_all(mean) %>%
  gather(species,reads,-Group)

dat <- full_join(meanCoverage,meanReads,by=c("species","Group"))
# top10
a <- dat[order(dat$reads,decreasing = T),]
a <- a$species[!duplicated(a$species)]
dat$lab <- ifelse(dat$species %in% a[1:10] & dat$Group=="Tissue",as.character(dat$species),"")

ggscatter(dat,x="reads",y="coverage",fill="Group",color="white",size = 4,shape = 21,alpha=0.8)+
  scale_x_log10()+
  ggrepel::geom_text_repel(aes(label=lab),max.overlaps = 1000,size=5)+
  scale_fill_manual(values = c("#9BB7DB","#862320"))+
  ylab("Coverage of genomes")+
  xlab("Average of Kraken-mapped reads")
ggsave("Fig1e.scatterplot.nReads.coverage.species.pdf",width = 5,height = 5)


#------------------------------------------------------------------------------
# Fig1f
n <- dat.wgs[dat.wgs$sampleType=="Tumor",]
n$cohort <- factor(n$cohort,levels = c("UU","UM","AC-ICAM"))

Boxplot_Reads2 <- function(d,x,y){
  ggboxplot(d,x=x,y=y,fill=x,ylab = "Reads Number",xlab = "")+
    scale_y_log10(limit=c(10,5e9))+
    stat_compare_means(method = "wilcox.test",comparisons = list(c(1,2),c(2,3),c(1,3)))+
    scale_fill_manual(values = c("#9EBCE3","#F1D1BB","#D5B226"))+
    theme(legend.position = "")
}

# WGS
p1 <- Boxplot_Reads2(n,x="cohort",y="nonHuman")
p2 <- Boxplot_Reads2(n,x="cohort",y="Microbiota")
p3 <- Boxplot_Reads2(n,x="cohort",y="TissueResident")

ggarrange(p1,p2,p3,nrow = 1,ncol = 3,align = "hv")
ggsave("Fig1f.boxplot.reads.pdf",width = 9,height = 4)
#-------------------------------------------------------------------------------



