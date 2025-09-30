library(openxlsx)
library(tidyverse)
library(ggpubr)
library(metacoder)

# read HRs for prognostic taxa generated from UU tumors and NATs
dat <- read.xlsx("TMP.order.HR.taxa.xlsx",rowNames = T)

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

l <- rownames(dat)
l2 <- dat$taxa

hr.tumor <- hrOS.T[l,]
hr.nat <- hrOS.N[l,]


#-----------------------------------------------------------------------------
# heattree

# preprocess
a <- dat
l1 <- a[grep("g__",rownames(a)),]
l2 <- a[grep("s__",rownames(a)),]
l1$Lineage <- apply(l1[,1:6],1,paste,collapse="|")
l2$Lineage <- apply(l2[,1:7],1,paste,collapse="|")
lstTaxa <- rbind(l1[,"Lineage",drop=F],l2[,"Lineage",drop=F])

dat$Lineage <- lstTaxa[rownames(dat),1]

# add significance
dat$If_T <- dat$hrOS.T
dat$If_N <- dat$hrOS.N
dat$Prognositic <- dat$Association


# add other infor.
dat$If_Both <- ifelse(dat$If_T==1 & dat$If_N==1, "Intersect","Specific")
dat$nodeColor <- ifelse(dat$Prognositic=="Shorter OS","#AA9276","#7CC9CA")
#write.csv(dat,"data.for.plot.csv")

# prognostic in Tumors
# configure parameters



dat <- dat[,11:16]


# ploting

# Tumor
obj <- parse_tax_data(dat,
                      class_cols = "Lineage", 
                      class_sep = "|",
                      class_regex = "^(.+)__(.+)$", 
                      class_key = c(tax_rank = "info", 
                                    tax_name = "taxon_name"))
obj$data$tax_abund <- calc_taxon_abund(obj, "tax_data")
a <- as.data.frame(obj$data$tax_data)
rownames(a) <- a$taxon_id
obj$data$tax_abund$Node_color <- a[obj$data$tax_abund$taxon_id,7]
obj$data$tax_abund$Node_color[is.na(obj$data$tax_abund$Node_color)] <- "gray"
obj$data$tax_abund$If_Both <- a[obj$data$tax_abund$taxon_id,6]
obj$data$tax_abund$If_Both[is.na(obj$data$tax_abund$If_Both)] <- "no"
obj$data$tax_abund$Node_color <- a[obj$data$tax_abund$taxon_id,7]
obj$data$tax_abund$Node_color[is.na(obj$data$tax_abund$Node_color)] <- "gray"
obj$data$tax_abund$Node_color[obj$data$tax_abund$If_T==0] <- "white"
obj$data$tax_abund$Label_color <- ifelse(obj$data$tax_abund$If_Both == "Intersect","#E74C3C","black")
obj$data$tax_abund$Label_color[obj$data$tax_abund$If_T==0] <- "white"
b <- obj$data$class_data
b <- b[!duplicated(b$taxon_id),]
b <- b[order(b$taxon_id),]
b$node_size <- 1
b$node_size[b$tax_rank=="g"] <- 3
b$node_size[b$tax_rank=="f"] <- 6
b$node_size[b$tax_rank=="o"] <- 8
b$node_size[b$tax_rank=="c"] <- 8
b$node_size[b$tax_rank=="p"] <- 15
b$node_size[b$tax_rank=="d"] <- 20
c <- b$node_size
names(c) <- b$taxon_id
d <- c
d[names(d)=="s"] <- 0.01

set.seed(10)
heat_tree(obj,  
          node_label = obj$taxon_names(),
          node_size = c,  
          node_label_size = d,
          layout = "davidson-harel", 
          initial_layout = "fruchterman-reingold",
          node_color = Node_color, 
          #node_color_interval = c(0.6, 1.4), 
          node_size_range = c(0.01, 0.02),
          node_color_range = c("#7CC9CA", "#AA9276"),
          node_label_color = Label_color,
          node_label_size_range = c(0.008,0.04),
          edge_size_range = c(0.001, 0.008),
          output_file = "Heat_tree.prognostic.Tumor.V2.pdf"
)


set.seed(10)
heat_tree(obj,  
          node_label = obj$taxon_names(),
          node_size = c,  
          node_label_size = d,
          layout = "davidson-harel", 
          initial_layout = "fruchterman-reingold",
          node_color = Node_color, 
          #node_color_interval = c(0.6, 1.4), 
          node_size_range = c(0.01, 0.02),
          node_color_range = c("#7CC9CA", "#AA9276"),
          node_label_color = Label_color,
          node_label_size_range = c(0.001,0.006),
          edge_size_range = c(0.001, 0.008),
          output_file = "Heat_tree.prognostic.Tumor.TMP.pdf"
)



# NAT
obj <- parse_tax_data(dat,
                      class_cols = "Lineage", 
                      class_sep = "|",
                      class_regex = "^(.+)__(.+)$", 
                      class_key = c(tax_rank = "info", 
                                    tax_name = "taxon_name"))
obj$data$tax_abund <- calc_taxon_abund(obj, "tax_data")
a <- as.data.frame(obj$data$tax_data)
rownames(a) <- a$taxon_id
obj$data$tax_abund$Node_color <- a[obj$data$tax_abund$taxon_id,7]
obj$data$tax_abund$Node_color[is.na(obj$data$tax_abund$Node_color)] <- "gray"
obj$data$tax_abund$If_Both <- a[obj$data$tax_abund$taxon_id,6]
obj$data$tax_abund$If_Both[is.na(obj$data$tax_abund$If_Both)] <- "no"
obj$data$tax_abund$Node_color <- a[obj$data$tax_abund$taxon_id,7]
obj$data$tax_abund$Node_color[is.na(obj$data$tax_abund$Node_color)] <- "gray"
obj$data$tax_abund$Node_color[obj$data$tax_abund$If_N==0] <- "white"
obj$data$tax_abund$Label_color <- ifelse(obj$data$tax_abund$If_Both == "Intersect","#E74C3C","black")
obj$data$tax_abund$Label_color[obj$data$tax_abund$If_N==0] <- "white"
b <- obj$data$class_data
b <- b[!duplicated(b$taxon_id),]
b <- b[order(b$taxon_id),]
b$node_size <- 1
b$node_size[b$tax_rank=="g"] <- 3
b$node_size[b$tax_rank=="f"] <- 6
b$node_size[b$tax_rank=="o"] <- 8
b$node_size[b$tax_rank=="c"] <- 8
b$node_size[b$tax_rank=="p"] <- 15
b$node_size[b$tax_rank=="d"] <- 20
c <- b$node_size
names(c) <- b$taxon_id
d <- c
d[names(d)=="s"] <- 0.01

set.seed(10)
heat_tree(obj,  
          node_label = obj$taxon_names(),
          node_size = c,  
          node_label_size = d,
          layout = "davidson-harel", 
          initial_layout = "fruchterman-reingold",
          node_color = Node_color, 
          #node_color_interval = c(0.6, 1.4), 
          node_size_range = c(0.01, 0.02),
          node_color_range = c("#7CC9CA", "#AA9276"),
          node_label_color = Label_color,
          node_label_size_range = c(0.008,0.04),
          edge_size_range = c(0.001, 0.008),
          output_file = "Heat_tree.prognostic.NAT.V2.pdf"
)
