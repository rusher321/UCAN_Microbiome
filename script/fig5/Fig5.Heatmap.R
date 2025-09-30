library(tidyverse)
library(ggpubr)
library(aplot)
library(openxlsx)

# Prognostic data
loadData1 <- function(){
  load("HR.taxa.rd")
  list(
    Tumor = r.uu.t,
    NAT = r.uu.n
  )
}
a <- loadData1()
hrOS.T <- a$Tumor
hrOS.N <- a$NAT
rownames(hrOS.T) <- hrOS.T$x
rownames(hrOS.N) <- hrOS.N$x


# 筛选
a <- hrOS.T$x[hrOS.T$fdr<0.05]
b <- hrOS.N$x[hrOS.N$P<0.01]
l <- union(a,b)
hrOS.T <- hrOS.T[l,]
hrOS.N <- hrOS.N[l,]

# 排序
ann <- read.csv("../../1.data.process/kk2.v2/kraken2.db.annotation.csv",header = T)
ann1 <- ann[!duplicated(ann$genus),]
ann2 <- ann[!duplicated(ann$species),]
rownames(ann1) <- ann1$genus
rownames(ann2) <- ann2$species
ann <- rbind(ann1,ann2)
ann$p <- gsub("_\\w$","",ann$p)
ann <- ann[,c(1,2,3,4,5,9)]
ann$f <- gsub("_\\w$","",ann$f)
ann$taxa <- rownames(ann)
rownames(ann) <- make.names(rownames(ann))
ann <- ann[l,]
ann$hrOS.T <- ifelse(hrOS.T$fdr<0.05,1,0)
ann$hrOS.N <- ifelse(hrOS.N$P<0.01,1,0)
ann$Association <- ifelse(hrOS.T$HR>1,"Shorter OS","Longer OS")
ann$Association[ann$hrOS.T==0 & hrOS.N$HR>1] <- "Shorter OS"

ann <- ann[order(ann$genus),]
ann <- ann[order(ann$f),]
ann <- ann[order(ann$p),]
ann <- ann[order(ann$hrOS.T,decreasing = F),]
ann <- ann[order(ann$hrOS.N,decreasing = T),]
ann <- ann[order(ann$Association,decreasing = T),]

#write.xlsx(ann,"TMP.order.HR.taxa.xlsx",rowNames=T)
ann <- read.xlsx("TMP.order.HR.taxa.xlsx",rowNames = T)

l <- rownames(ann)
l2 <- ann$taxa

hrOS.T <- hrOS.T[l,]
hrOS.N <- hrOS.N[l,]

# Tumor vs. NAT
a <- readRDS("../Part1/MaAsLin.Tumor_NAT.rds")
grp.TN <- a$UU

# Right vs. left
a <- readRDS("../Part1/MaAsLin.Location.rds")
grp.RL.T <- a$UU.Tumor
grp.RL.N <- a$UU.NAT

# HM vs. nHM
a <- readRDS("../Part2/glm.HM.rds")
grp.HM <- a$UU

# add groups
grp.TN <- grp.TN[l2,]
grp.HM <- grp.HM[l2,]
grp.RL.N <- grp.RL.N[l2,]
grp.RL.T <- grp.RL.T[l2,]

colnames(grp.HM) <- colnames(grp.TN)
colnames(grp.RL.N) <- colnames(grp.TN)
colnames(grp.RL.T) <- colnames(grp.TN)


# heatmap
# 1
a <- data.frame(
  Tumor = hrOS.T$HR,
  NAT = hrOS.N$HR,
  Taxa = rownames(hrOS.T),
  IF.Tumor = ann$hrOS.T,
  IF.NAT = ann$hrOS.N
) 
a$Taxa <- factor(a$Taxa,levels = l,labels = l2)
b <- a[,c("Tumor","NAT","Taxa")] %>%
  gather(y,HR,-Taxa) %>%
  mutate(
    y = factor(y,levels=c("NAT","Tumor"),labels=c("HR (NAT)","HR (Tumor)"))
  )
c <- a[,c("Taxa","IF.Tumor","IF.NAT")] %>%
  gather(y,lab,-Taxa)
c <- c[c$lab!=0,]
c$y <- factor(c$y,levels=c("IF.NAT","IF.Tumor"),labels = c("HR (NAT)","HR (Tumor)"))

p1 <- 
ggplot()+
  geom_tile(b,mapping=aes(x=Taxa,y=y,fill=HR))+
  geom_point(data=c,mapping=aes(x=Taxa,y=y),shape=5,size=3)+
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
p1

# T vs. N
a <- grp.TN %>%
  rownames_to_column(var = "Taxa") %>%
  mutate(
    Taxa = factor(Taxa,levels=l2),
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
p2

# Location
a <- data.frame(
  Tumor = grp.RL.T$Effectsize,
  NAT = grp.RL.N$Effectsize,
  Taxa = rownames(grp.RL.T)
) %>%
  gather(y,Effectsize,-Taxa) %>%
  mutate(
    y = factor(y,levels=c("NAT","Tumor"),labels=c("Right vs. Left (NAT)","Right vs. Left (Tumor)")),
    Taxa = factor(Taxa,levels=l2)
  )
b <- data.frame(
  Tumor = grp.RL.T$FDR,
  NAT = grp.RL.N$FDR,
  Taxa = rownames(grp.RL.T)
) %>%
  gather(y,FDR,-Taxa) %>%
  mutate(
    y = factor(y,levels=c("NAT","Tumor"),labels=c("Right vs. Left (NAT)","Right vs. Left (Tumor)")),
    Taxa = factor(Taxa,levels=l2)
  )

p3 <- 
  ggplot()+
  geom_tile(a,mapping=aes(x=Taxa,y=y,fill=Effectsize))+
  geom_point(data=b[b$FDR<0.05,],mapping=aes(x=Taxa,y=y),shape=1,size=3)+
  scale_x_discrete(expand = expansion(mult = c(0,0)))+
  scale_y_discrete(expand = expansion(mult = c(0,0)))+
  #scale_fill_gradient2(high = "#248888",low = "#F0D879",mid = "white",midpoint = 0,limits=c(-3,3))+
  scale_fill_gradientn(
    colours = c("#F0D879", "white", "#248888"),
    values = scales::rescale(c(-3, 0, 1.3)),  # 替换min_value和max_value为你的实际数据范围
    limits = c(-3, 1.3)  # 设置相同的限制确保颜色映射正确
  )+
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
    Taxa = factor(Taxa,levels=l2),
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
p4

# annotation
a <- ann[,c("taxa","p"),drop=F]
l1 <- c("p__Firmicutes","p__Actinobacteriota","p__Proteobacteria","p__Bacteroidota","p__Verrucomicrobiota")
a <- a %>%
  mutate(
    Taxa = factor(taxa,levels=l2),
    y = "Phylum"
  )
p5 <- 
ggplot()+
  geom_tile(a,mapping=aes(x=Taxa,y=y,fill=p))+
  scale_x_discrete(expand = expansion(mult = c(0,0)))+
  scale_y_discrete(expand = expansion(mult = c(0,0)))+
  scale_fill_brewer(palette = "Set3")+
  theme_bw()+
  theme(
    axis.text.x = element_text(angle = 270,hjust = 0, vjust = 0.5),
    axis.text = element_text(color=1),
    panel.grid = element_blank()
  )
p5

p5 %>% 
  insert_top(p4,height = 1) %>%
  insert_top(p3,height = 2) %>%
  insert_top(p2,height = 1) %>%
  insert_top(p1,height = 2)
ggsave("Fig5d.heatmap.HR.V2.pdf",width = 24,height = 7)


#----------------------------------------------------------------------------------
a <- data.frame(
  x = c("1","2","3","4","5","6"),
  y = c(sum(ann$hrOS.T),sum(ann$hrOS.N),sum(grp.TN$FDR<0.05),sum(grp.RL.T$FDR<0.05),sum(grp.RL.N$FDR<0.05),sum(grp.HM$FDR<0.05)),
  g = c("1","1","2","3","3","4")
)
a$x <- factor(a$x,levels = rev(a$x))

ggbarplot(a,x="x",y="y",fill="gray")+
  coord_flip()+
  facet_grid(g~.,scale="free",space="free")+
  geom_text(aes(label=y),hjust=-0.2)+
  scale_y_continuous(expand = expansion(mult = c(0,0.2)))+
  xlab("")+
  ylab("Number of significantly different tax")+
  theme(
    strip.text = element_blank(),
    axis.text.y = element_blank()
  )
ggsave("Fig5d.barplot.pdf",width = 3,height = 3)


