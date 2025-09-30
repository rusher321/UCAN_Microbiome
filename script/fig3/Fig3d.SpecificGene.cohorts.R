########  Load data  ###########################################################

library(openxlsx)
library(tidyverse)

# -----------------------------------------------------------------------------
# Read associaion results for driver genes and DDR genes
# -----------------------------------------------------------------------------
a        <- readRDS("drivers.rds")   # driver-gene associations
b        <- readRDS("DDRs.rds")      # DDR-gene associations
dat.uu1  <- b$UU                          # DDR results UU cohort only

# Combine driver + DDR estimates / P-values for each cohort
dat.uu   <- cbind(a$UU, b$UU[rownames(a$UU), ])
dat.um   <- cbind(a$UM, b$UM[rownames(a$UM), ])
dat.icam <- cbind(a$ICAM, b$ICAM[rownames(a$ICAM), ])

# Align rows / columns across cohorts
dat.uu   <- dat.uu[rownames(dat.um), colnames(dat.um)]

# Ensure ICAM contains exactly the same columns as UU
for (x in 1:ncol(dat.uu)) {
  if (!colnames(dat.uu)[x] %in% colnames(dat.icam)) {
    dat.icam[, colnames(dat.uu)[x]] <- NA
  }
}
dat.icam <- dat.icam[, colnames(dat.uu)]
dat.icam <- dat.icam[rownames(dat.uu), ]

# Impute missing values 
for (x in 1:ncol(dat.icam)) {
  n <- colnames(dat.icam)[x]
  if (grepl("Pval", n))     dat.icam[is.na(dat.icam[, x]), x] <- 1
  if (grepl("FDR", n))      dat.icam[is.na(dat.icam[, x]), x] <- 1
  if (grepl("estimate", n)) dat.icam[is.na(dat.icam[, x]), x] <- 0
}


# -----------------------------------------------------------------------------
# Build consensus significance matrix
# -----------------------------------------------------------------------------
# Extract P-value matrices
p.uu   <- dat.uu[, grep("Pval", colnames(dat.uu))]
p.um   <- dat.um[, grep("Pval", colnames(dat.um))]
p.icam <- dat.icam[, grep("Pval", colnames(dat.icam))]
colnames(p.uu) <- colnames(p.um) <- colnames(p.icam) <- gsub("\\..*", "", colnames(p.uu))

# Extract effect-size matrices
r.uu   <- dat.uu[, grep("estimate", colnames(dat.uu))]
r.um   <- dat.um[, grep("estimate", colnames(dat.um))]
r.icam <- dat.icam[, grep("estimate", colnames(dat.icam))]
colnames(r.uu) <- colnames(r.um) <- colnames(r.icam) <- gsub("\\..*", "", colnames(r.uu))

# Convert to signed significance (-1 / 0 / 1)
a <- apply(p.uu, 2, function(x) ifelse(x < 0.05, 1, 0)) * sign(r.uu)
b <- apply(p.um, 2, function(x) ifelse(x < 0.05, 1, 0)) * sign(r.um)
c <- apply(p.icam, 2, function(x) ifelse(x < 0.05, 1, 0)) * sign(r.icam)

# Require UU signal and concordant direction in ≥1 replication cohort
b[a == 0] <- 0
c[a == 0] <- 0

d.uu <- a
for (i in 1:nrow(a)) {
  for (j in 1:ncol(a)) {
    if (a[i, j] == 1  && (b[i, j] == 1  || c[i, j] == 1))  d.uu[i, j] <- 1
    else if (a[i, j] == -1 && (b[i, j] == -1 || c[i, j] == -1)) d.uu[i, j] <- -1
    else d.uu[i, j] <- 0
  }
}

# Retain taxa with ≥1 consensus hit
d.uu <- d.uu[apply(d.uu, 1, function(x) any(x != 0)), ]
d.uu <- d.uu[, apply(d.uu, 2, function(x) any(x != 0))]

# Filter companion matrices identically
a <- a[rownames(d.uu), colnames(d.uu)]
b <- b[rownames(d.uu), colnames(d.uu)]
c <- c[rownames(d.uu), colnames(d.uu)]
a[d.uu == 0] <- 0
b[d.uu == 0] <- 0
c[d.uu == 0] <- 0

# Label matrix (sum of signs)
lab <- a + b + c
lab[lab == 0] <- ""

# Order genes / taxa by number of hits
tmp  <- data.frame(n = apply(d.uu, 2, function(x) sum(x != 0))); tmp  <- tmp[order(tmp$n, decreasing = TRUE), , drop = FALSE]
tmp2 <- data.frame(n = apply(d.uu, 1, function(x) sum(x != 0))); tmp2 <- tmp2[order(tmp2$n, decreasing = TRUE), , drop = FALSE]

d.uu <- d.uu[rownames(tmp2), rownames(tmp)]
lab  <- lab[rownames(tmp2), rownames(tmp)]

# -----------------------------------------------------------------------------
# Long-format for ggplot
# -----------------------------------------------------------------------------
d1 <- d.uu %>% rownames_to_column(var = "taxa") %>% gather(gene, symbols, -taxa)
d2 <- lab %>% rownames_to_column(var = "taxa") %>% gather(gene, lab, -taxa)

dat <- full_join(d1, d2, by = c("gene", "taxa"))

# Validation status
dat$Validation <- NA
for (i in 1:nrow(dat)) {
  x1 <- b[dat$taxa[i], dat$gene[i]]
  x2 <- c[dat$taxa[i], dat$gene[i]]
  if (x1 != 0) dat$Validation[i] <- "UM"
  if (x2 != 0) dat$Validation[i] <- "ICAM"
  if (dat$lab[i] %in% c("3", "-3")) dat$Validation[i] <- "Both"
}
dat$Validation[is.na(dat$Validation)] <- "Not reproducible"

# Factor levels for plotting
dat$taxa      <- factor(dat$taxa, levels = rev(rownames(tmp2)))
dat$gene      <- factor(dat$gene, levels = rownames(tmp))
dat$Validation <- factor(dat$Validation, levels = c("Both", "UM", "ICAM", "Not reproducible"))

# -----------------------------------------------------------------------------
# Add UU FDR status
# -----------------------------------------------------------------------------
f.uu <- dat.uu[, grep("FDR", colnames(dat.uu))]
colnames(f.uu) <- gsub("\\..*", "", colnames(f.uu))
f <- f.uu[rownames(d.uu), colnames(d.uu)]
p <- p.uu[rownames(d.uu), colnames(d.uu)]

dat$P.uu   <- NA
dat$FDR.uu <- NA
for (x in 1:nrow(dat)) {
  dat$P.uu[x]   <- p[as.character(dat$taxa[x]), as.character(dat$gene[x])]
  dat$FDR.uu[x] <- f[as.character(dat$taxa[x]), as.character(dat$gene[x])]
}

dat$Sig.uu <- ifelse(dat$P.uu < 0.05, "P<0.05", "ns")
dat$Sig.uu[dat$FDR.uu < 0.05] <- "FDR<0.05"
dat$Sig.uu <- factor(dat$Sig.uu, levels = c("ns", "P<0.05", "FDR<0.05"))
dat <- dat[order(dat$Sig.uu), ]

# Move APC gene to last column
l <- levels(dat$gene)
dat$gene <- factor(dat$gene, levels = c(l[l != "APC"], "APC"))

# Association direction
dat$Association <- ifelse(dat$Sig.uu != "ns", "positive", "ns")
dat$Association[dat$gene == "APC" & dat$Association == "positive"] <- "negative"

# Gene group (DDR vs Driver)
lst.ddr <- gsub("\\..*", "", colnames(dat.uu1))
dat$Group <- ifelse(dat$gene %in% lst.ddr, "DDR gene", "Driver gene")
dat$Group <- factor(dat$Group, levels = c("Driver gene", "DDR gene"))

# Fill colour combines significance + FDR
dat$Fill <- as.character(dat$Association)
dat$Fill[dat$Association == "positive" & dat$FDR.uu < 0.05] <- "positive2"
dat$Fill[dat$Association == "negative" & dat$FDR.uu < 0.05] <- "negative2"
dat$Fill <- factor(dat$Fill, levels = c("negative2", "negative", "ns", "positive", "positive2"))

# Remove Escherichia (poor annotation)
dat <- dat[grep("Escherichia", dat$taxa, invert = TRUE), ]
dat$taxa <- droplevels(dat$taxa)

# -----------------------------------------------------------------------------
# Main figure 3d: heat-map
# -----------------------------------------------------------------------------
p1 <- ggplot(dat, aes(x = gene, y = taxa, fill = Validation)) +
  geom_tile(aes(fill = Fill), color = "black") +
  geom_point(data = dat[dat$Validation %in% c("UM", "Both"), ], shape = 17) +
  geom_point(data = dat[dat$Validation %in% c("ICAM", "Both"), ], shape = 1, size = 3) +
  facet_grid(~Group, scales = "free", space = "free") +
  scale_fill_manual(values = c(positive2 = "#ea8c8c", positive = "#FFE5D9",
                               negative2 = "#3769c1", negative = "#99CCFF", ns = "white")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust = 0))
p1

# -----------------------------------------------------------------------------
# Sidebar: # species per gene
# -----------------------------------------------------------------------------
a <- tmp
a$gene  <- factor(rownames(a), levels = levels(dat$gene))
a$Group <- ifelse(a$gene %in% lst.ddr, "DDR gene", "Driver gene")
a$Group <- factor(a$Group, levels = c("Driver gene", "DDR gene"))
b <- data.frame(n = apply(f.uu, 2, function(x) sum(x < 0.05)))
a$N <- b[rownames(a), ]

p2 <- ggplot(a, aes(gene, n)) +
  geom_bar(stat = "identity", width = 0.5) +
  facet_grid(~Group, scales = "free", space = "free") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.grid = element_blank(), legend.position = "")
p2

# -----------------------------------------------------------------------------
# Top bar: # genes per species
# -----------------------------------------------------------------------------
a <- tmp2
a <- a[rownames(a) %in% levels(dat$taxa), , drop = FALSE]
a$taxa  <- factor(rownames(a), levels = rev(rownames(a)))
a$Group <- "x"
b <- data.frame(n = apply(f.uu, 1, function(x) sum(x < 0.05)))
a$N <- b[rownames(a), ]

p3 <- ggplot(a, aes(taxa, n)) +
  geom_bar(stat = "identity", width = 0.5) +
  facet_grid(~Group) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  coord_flip() +
  xlab("") +
  theme_bw() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid = element_blank())
p3

# -----------------------------------------------------------------------------
# Gene prevalence in HM vs nHM (UU cohort)
# -----------------------------------------------------------------------------
dat.gene1 <- read.table("../../../../2.SourceData/03.Mutation/CRC-SW.Ensemble.1063_DNBSEQ.20210706.lite.coding_splice.nonsyn.mut.n.tsv",
                        header = TRUE, sep = "\t")
geneList1 <- read.xlsx("../../../../2.SourceData/03.Mutation/DR pathway gene.xlsx", rowNames = TRUE)
rownames(geneList1) <- geneList1$Genes
dat.gene2 <- read.xlsx("../../../../2.SourceData/03.Mutation/nonsyn.mut.n_dNdScv-1.xlsx", sheet = 1)
geneList2 <- read.xlsx("../../../../2.SourceData/03.Mutation/96 driver from cong.xlsx")[, c("gene_name", "pathway.f2")]
colnames(geneList2) <- c("Genes", "Pathway")
rownames(geneList2) <- geneList2$Genes

process1 <- function(dat.gene, geneList) {
  dat.gene <- dat.gene %>% spread(Hugo_Symbol, Num_Nonsyn_Mut) %>% column_to_rownames(var = "Tumor_Sample_Barcode")
  dat.gene[is.na(dat.gene)] <- 0
  rownames(dat.gene) <- gsub("B$", "", rownames(dat.gene))
  dat.gene <- dat.gene[, rownames(geneList)]
  dat.gene[dat.gene > 0] <- 1
  dat.gene
}

dat.gene1 <- process1(dat.gene1, geneList1)
dat.gene2 <- process1(dat.gene2, geneList2)
dat.gene  <- cbind(dat.gene1, dat.gene2)

# Restrict to UU samples
phe <- read.csv("../../../../2.SourceData/01.phe/20230727-CRC-SW.20230704-v5.clinical.csv", row.names = 1)
phe <- phe[phe$Sample_Center == "UU", ]
dat.gene <- dat.gene[phe$Tumor_ID, ]
dat.gene$HM <- factor(phe$HM_Status, levels = c("HM", "nHM"))

pre <- dat.gene %>%
  group_by(HM) %>%
  summarise_all(function(x) mean(x, na.rm = TRUE)) %>%
  gather(Gene, Prevalence, -HM)

pre <- pre[pre$Gene %in% levels(dat$gene), ]
pre$Gene  <- factor(pre$Gene, levels = levels(dat$gene))
pre$Group <- ifelse(pre$Gene %in% lst.ddr, "DDR gene", "Driver gene")
pre$Group <- factor(pre$Group, levels = c("Driver gene", "DDR gene"))

p4 <- ggplot(pre, aes(Gene, HM, fill = Prevalence)) +
  geom_tile() +
  facet_grid(~Group, scales = "free", space = "free") +
  scale_fill_gradient(low = "white", high = "black") +
  scale_y_discrete(expand = expansion(mult = c(0, 0))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust = 0),
        axis.text = element_text(color = 1),
        axis.title = element_text(color = 1),
        panel.grid = element_blank())
p4

