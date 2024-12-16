library(tidyverse)
library(ggpubr)
library(aplot)
library(openxlsx)

# source data
lst <- read.xlsx("HZ-List.Prognosis.Taxa.V2.xlsx",sheet = "Zhun-modified",rowNames = T)

# groups
grp.TN <- readRDS("../../../2.SourceData/00.SourceData/MaAsLin2.NAT_Tumor.UCAN.rds")
grp.TN <- grp.TN$UU.All
grp.RL <- readRDS("../../../2.SourceData/00.SourceData/MaAsLin2.Location.UCAN.rds")
grp.RL.T <- grp.RL$UU.Tumor
grp.RL.N <- grp.RL$UU.Normal
grp.HM <- readRDS("../../../2.SourceData/00.SourceData/MaAsLin2.HM_nHM.UCAN.rds")
grp.HM <- grp.HM$UU.Tumor

# add groups
grp.TN <- grp.TN[rownames(lst),]
grp.HM <- grp.HM[rownames(lst),]
grp.RL.N <- grp.RL.N[rownames(lst),]
grp.RL.T <- grp.RL.T[rownames(lst),]

colnames(grp.HM) <- colnames(grp.TN)
colnames(grp.RL.N) <- colnames(grp.TN)
colnames(grp.RL.T) <- colnames(grp.TN)

# HR
load("../Validation.PrognosticTaxa/Backup/uu_hr.Rd")
load("../Validation.PrognosticTaxa/Backup/uu_normal_hr.rd")
hrOS.T <- r1_uu_wgs 
hrOS.N <- r1_uu_normal_wgs
rownames(hrOS.T) <- hrOS.T$x
rownames(hrOS.N) <- hrOS.N$x

# rename
l <- read.xlsx("../../../2.SourceData/02.profile/Kraken2/Mannual.K2.taxa.reName.xlsx",rowNames = T)
rownames(l)[grep("s__",rownames(l),invert = T)] <- paste0("g__",rownames(l)[grep("s__",rownames(l),invert = T)])
rownames(l) <- make.names(rownames(l))
l <- l[rownames(hrOS.T),]
hrOS.N <- hrOS.N[rownames(hrOS.T),]
rownames(hrOS.T) <- rownames(hrOS.N) <- l$StandardName

hrOS.T <- hrOS.T[rownames(lst),]
hrOS.N <- hrOS.N[rownames(lst),]



# heatmap
# 1
a <- data.frame(
  Tumor = hrOS.T$HR,
  NAT = hrOS.N$HR,
  Taxa = rownames(hrOS.T)
) %>%
  gather(y,HR,-Taxa) %>%
  mutate(
    y = factor(y,levels=c("NAT","Tumor"),labels=c("HR (NAT)","HR (Tumor)")),
    Taxa = factor(Taxa,levels=rownames(lst))
  )
b <- lst[,6:7] %>%
  rownames_to_column(var = "Taxa") %>%
  gather(y,lab,-Taxa)
b <- b[b$lab!=0,]
b$y <- factor(b$y,levels=c("normal","tumor"),labels=c("HR (NAT)","HR (Tumor)"))

p1 <- 
ggplot()+
  geom_tile(a,mapping=aes(x=Taxa,y=y,fill=HR))+
  geom_point(data=b,mapping=aes(x=Taxa,y=y),shape=5,size=3)+
  scale_x_discrete(expand = expansion(mult = c(0,0)))+
  scale_y_discrete(expand = expansion(mult = c(0,0)))+
  scale_fill_gradient2(low = "#55aaaa",high = "#aa9074",mid = "white",midpoint = 1)+
  theme_bw()+
  theme(
    #axis.text.x = element_text(angle = 270,hjust = 0, vjust = 0.5),
    axis.text.x = element_blank(),
    axis.text = element_text(color=1),
    panel.grid = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank()
  )


# T vs. N
a <- grp.TN %>%
  rownames_to_column(var = "Taxa") %>%
  mutate(
    Taxa = factor(Taxa,levels=rownames(lst)),
    y = "Tumor vs. NAT"
  )
p2 <- 
ggplot()+
  geom_tile(data=a,mapping=aes(x=Taxa,y=y,fill=Effectsize))+
  geom_point(data=a[which(a$FDR<0.05),],mapping=aes(x=Taxa,y=y),shape=1,size=3)+
  scale_x_discrete(expand = expansion(mult = c(0,0)))+
  scale_y_discrete(expand = expansion(mult = c(0,0)))+
  scale_fill_gradient2(low = "#848CCF",high = "#BE5683",mid = "white",midpoint = 0,breaks=c(-0.5,0,0.5,1),limits=c(-0.5,1))+
  theme_bw()+
  theme(
    #axis.text.x = element_text(angle = 270,hjust = 0, vjust = 0.5),
    axis.text.x = element_blank(),
    axis.text = element_text(color=1),
    panel.grid = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank()
  )

# Location
a <- data.frame(
  Tumor = grp.RL.T$Effectsize,
  NAT = grp.RL.N$Effectsize,
  Taxa = rownames(grp.RL.T)
) %>%
  gather(y,Effectsize,-Taxa) %>%
  mutate(
    y = factor(y,levels=c("NAT","Tumor"),labels=c("Right vs. Left (NAT)","Right vs. Left (Tumor)")),
    Taxa = factor(Taxa,levels=rownames(lst))
  )
b <- data.frame(
  Tumor = grp.RL.T$FDR,
  NAT = grp.RL.N$FDR,
  Taxa = rownames(grp.RL.T)
) %>%
  gather(y,FDR,-Taxa) %>%
  mutate(
    y = factor(y,levels=c("NAT","Tumor"),labels=c("Right vs. Left (NAT)","Right vs. Left (Tumor)")),
    Taxa = factor(Taxa,levels=rownames(lst))
  )

p3 <- 
  ggplot()+
  geom_tile(a,mapping=aes(x=Taxa,y=y,fill=Effectsize))+
  geom_point(data=b[b$FDR<0.05,],mapping=aes(x=Taxa,y=y),shape=1,size=3)+
  scale_x_discrete(expand = expansion(mult = c(0,0)))+
  scale_y_discrete(expand = expansion(mult = c(0,0)))+
  scale_fill_gradient2(low = "#248888",high = "#F0D879",mid = "white",midpoint = 1)+
  theme_bw()+
  theme(
    #axis.text.x = element_text(angle = 270,hjust = 0, vjust = 0.5),
    axis.text.x = element_blank(),
    axis.text = element_text(color=1),
    panel.grid = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank()
  )
p3

# HM
a <- grp.HM %>%
  rownames_to_column(var = "Taxa") %>%
  mutate(
    Taxa = factor(Taxa,levels=rownames(lst)),
    y = "nHM vs. HM (Tumor)"
  )
p4 <- 
  ggplot()+
  geom_tile(data=a,mapping=aes(x=Taxa,y=y,fill=Effectsize))+
  geom_point(data=a[which(a$FDR<0.05),],mapping=aes(x=Taxa,y=y),shape=1,size=3)+
  scale_x_discrete(expand = expansion(mult = c(0,0)))+
  scale_y_discrete(expand = expansion(mult = c(0,0)))+
  scale_fill_gradient2(low = "#4aa3ba",high = "#da5c53",mid = "white",midpoint = 0,breaks=c(-0.5,0,0.5,1),limits=c(-0.9,0.9))+
  theme_bw()+
  theme(
    #axis.text.x = element_text(angle = 270,hjust = 0, vjust = 0.5),
    axis.text.x = element_blank(),
    axis.text = element_text(color=1),
    panel.grid = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank()
  )


# annotation
a <- lst[,c(1,4)]
l1 <- c("p__Firmicutes","p__Actinobacteriota","p__Proteobacteria","p__Bacteroidota","p__Verrucomicrobiota")
l2 <- c("f__Lachnospiraceae","f__Ruminococcaceae","f__Clostridiaceae","Others")
a$Family[!a$Family %in% l2] <- "Others"
a <- a %>%
  rownames_to_column(var = "Taxa") %>%
  gather(y,value,-Taxa) %>%
  mutate(
    Taxa = factor(Taxa,levels=rownames(lst)),
    value = factor(value,levels=c(l1,l2)),
    y = factor(y,levels=c("Family","Phylum"))
  )
a <- a[a$y=="Phylum",]
p5 <- 
ggplot()+
  geom_tile(a,mapping=aes(x=Taxa,y=y,fill=value))+
  scale_x_discrete(expand = expansion(mult = c(0,0)))+
  scale_y_discrete(expand = expansion(mult = c(0,0)))+
  scale_fill_brewer(palette = "Set3")+
  theme_bw()+
  theme(
    axis.text.x = element_text(angle = 270,hjust = 0, vjust = 0.5),
    axis.text = element_text(color=1),
    panel.grid = element_blank()
  )

p5 %>% 
  insert_top(p4,height = 1) %>%
  insert_top(p3,height = 2) %>%
  insert_top(p2,height = 1) %>%
  insert_top(p1,height = 2)
ggsave("Fig5.heatmap.HR.V2.pdf",width = 24,height = 7)






