library(openxlsx)
library(tidyverse)
library(ggpubr)
library(metacoder)

# read HRs for prognostic taxa generated from UU tumors and NATs
load("uu_hr.Rd") 
load("uu_normal_hr.rd")
hr.tumor <- r1_uu_wgs
hr.nat <- r1_uu_normal_wgs

# re-name the taxa
ann.taxa <- read.xlsx("../../../2.SourceData/02.profile/Kraken2/Mannual.K2.taxa.reName.xlsx",rowNames = T)
rownames(ann.taxa)[grep("s__",rownames(ann.taxa),invert = T)] <- paste0("g__",rownames(ann.taxa)[grep("s__",rownames(ann.taxa),invert = T)])
rownames(ann.taxa) <- make.names(rownames(ann.taxa))
all(hr.tumor$x %in% rownames(ann.taxa)) # TRUE
hr.tumor$x <- ann.taxa[hr.tumor$x,"StandardName"]
hr.nat$x <- ann.taxa[hr.nat$x,"StandardName"]

# make sure the same order of taxa
rownames(hr.tumor) <- hr.tumor$x
rownames(hr.nat) <- hr.nat$x
hr.nat <- hr.nat[rownames(hr.tumor),]


# debug
sum(hr.tumor$fdr<0.05) # 136
sum(hr.nat$P<0.01)     # 35

#_-----------------------------------------------------------------------------
# only focus on the prognostic taxa
# FDR<0.05 for tumors
# P  <0.01 for NATs
lst.tumor <- rownames(hr.tumor)[hr.tumor$fdr<0.05]
lst.nat   <- rownames(hr.nat)[hr.nat$P<0.01]
a <- union(lst.tumor,lst.nat)

hr.tumor <- hr.tumor[a,]
hr.nat   <- hr.nat[a,]

#-----------------------------------------------------------------------------
# consistency analysis
a <- hr.tumor[order(hr.tumor$HR),]
b <- hr.nat[rownames(a),]
dat <- data.frame(
  row.names = rownames(a),
  Tumor = a$HR,
  NAT = b$HR
)

ggplot(dat,aes(Tumor,NAT))+
  geom_point()+
  geom_vline(xintercept = 1,linetype="dashed")+
  geom_hline(yintercept = 1,linetype="dashed")+
  theme_pubr()+
  xlab("HR for Tumors")+
  ylab("HR for NATs")+
  labs(title = "148 prognostic taxa")

# the result showed that the associations between tissue-resident taxa and prognosis
# are consistent in tumor and normal



#-----------------------------------------------------------------------------
# heattree

# preprocess
# standarize the lineage
a <- ann.taxa[,"Lineage",drop=F] %>% separate(Lineage,c("d","p","c","o","f","g"),sep = "\\|")
rownames(a) <- ann.taxa$StandardName
l1 <- a[grep("g__",rownames(a)),]
l2 <- a[grep("s__",rownames(a)),]
l <- l1
l$Genus <- rownames(l)
rownames(l) <- l$g
l2$g <- l[l2$g,"Genus"]
l2$s <- rownames(l2)
l1$g <- rownames(l1)
# combine
l1$Lineage <- apply(l1,1,paste,collapse="|")
l2$Lineage <- apply(l2,1,paste,collapse="|")
lstTaxa <- rbind(l1[,"Lineage",drop=F],l2[,"Lineage",drop=F])

dat$Lineage <- lstTaxa[rownames(dat),1]

# add significance
hr.tumor <- hr.tumor[rownames(dat),]
hr.nat <- hr.nat[rownames(dat),]
dat$If_T <- ifelse(hr.tumor$fdr<0.05,1,0)
dat$If_N <- ifelse(hr.nat$P<0.01,1,0)
dat$Prognositic <- NA
dat$Prognositic[dat$NAT>1 & dat$If_N==1] <- "Shorter OS"
dat$Prognositic[dat$NAT<1 & dat$If_N==1] <- "Longer OS"
dat$Prognositic[dat$Tumor>1 & dat$If_T==1] <- "Shorter OS"
dat$Prognositic[dat$Tumor<1 & dat$If_T==1] <- "Longer OS"


# add other infor.
dat$If_Both <- ifelse(dat$If_T==1 & dat$If_N==1, "Intersect","Specific")
dat$nodeColor <- ifelse(dat$Prognositic=="Shorter OS","#AA9276","#7CC9CA")
#write.csv(dat,"data.for.plot.csv")

# prognostic in Tumors
# configure parameters



# ploting

# Tumor
set.seed(1)
obj <- parse_tax_data(dat,
                      class_cols = "Lineage", 
                      class_sep = "|",
                      class_regex = "^(.+)__(.+)$", 
                      class_key = c(tax_rank = "info", 
                                    tax_name = "taxon_name"))
obj$data$tax_abund <- calc_taxon_abund(obj, "tax_data")
a <- as.data.frame(obj$data$tax_data)
rownames(a) <- a$taxon_id
obj$data$tax_abund$Node_color <- a[obj$data$tax_abund$taxon_id,9]
obj$data$tax_abund$Node_color[is.na(obj$data$tax_abund$Node_color)] <- "gray"
obj$data$tax_abund$If_Both <- a[obj$data$tax_abund$taxon_id,8]
obj$data$tax_abund$If_Both[is.na(obj$data$tax_abund$If_Both)] <- "no"
obj$data$tax_abund$Node_color <- a[obj$data$tax_abund$taxon_id,9]
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
heat_tree(obj,  
          node_label = obj$taxon_names(),
          node_size = c,  
          node_label_size = d,
          layout = "davidson-harel", 
          initial_layout = "fruchterman-reingold",
          node_color = Node_color, 
          #node_color_interval = c(0.6, 1.4), 
          node_size_range = c(0.01, 0.03),
          node_color_range = c("#7CC9CA", "#AA9276"),
          node_label_color = Label_color,
          #node_label_size_range = c(0.001,0.04),
          edge_size_range = c(0.004, 0.005)#,
          #output_file = "Heat_tree.prognostic.Tumor.V2.pdf"
)

# NAT
set.seed(1)
obj <- parse_tax_data(dat,
                      class_cols = "Lineage", 
                      class_sep = "|",
                      class_regex = "^(.+)__(.+)$", 
                      class_key = c(tax_rank = "info", 
                                    tax_name = "taxon_name"))
obj$data$tax_abund <- calc_taxon_abund(obj, "tax_data")
a <- as.data.frame(obj$data$tax_data)
rownames(a) <- a$taxon_id
obj$data$tax_abund$Node_color <- a[obj$data$tax_abund$taxon_id,9]
obj$data$tax_abund$Node_color[is.na(obj$data$tax_abund$Node_color)] <- "gray"
obj$data$tax_abund$If_Both <- a[obj$data$tax_abund$taxon_id,8]
obj$data$tax_abund$If_Both[is.na(obj$data$tax_abund$If_Both)] <- "no"
obj$data$tax_abund$Node_color <- a[obj$data$tax_abund$taxon_id,9]
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
heat_tree(obj,  
          node_label = obj$taxon_names(),
          node_size = c,  
          node_label_size = d,
          layout = "davidson-harel", 
          initial_layout = "reingold-tilford",
          node_color = Node_color, 
          #node_color_interval = c(0.6, 1.4), 
          node_size_range = c(0.01, 0.03),
          node_color_range = c("#7CC9CA", "#AA9276"),
          node_label_color = Label_color,
          #node_label_size_range = c(0.001,0.04),
          edge_size_range = c(0.004, 0.005),
          output_file = "Heat_tree.prognostic.NAT.V2.pdf"
)



#-----------------------------------------------------------------------------








heat_tree(obj,  
          node_label = taxon_names,
          node_size = n_leaves,
          node_color = Node_color,
          #node_color_interval = c(0.6, 1.4), 
          node_color_range = c("#7CC9CA", "#AA9276"),
          node_size_range = c(0.01, 0.03),
          node_label_size_range = c(0.013,0.04),
          node_label_color = Label_color,
          edge_size_range = c(0.001, 0.001),
          initial_layout = "re",
          layout = "da",
          output_file = "Heat_tree.prognostic.Tumor.pdf"
)
heat_tree(obj,  node_label = obj$taxon_names(),
          node_size = obj$n_obs(),  
          layout = "davidson-harel", 
          initial_layout = "reingold-tilford",
          node_color = Node_color, 
          node_color_interval = c(0.6, 1.4), 
          node_color_range = c("#7CC9CA", "#AA9276"),
          node_label_color = Label_color,
          node_label_size_range = c(0.005,0.005),
          node_size_range = c(0.01,0.02),
          edge_size_range = c(0.001, 0.001),
          output_file = "Heat_tree.prognostic.Tumor.V2.pdf"
)

