######## load data #################################################

library(openxlsx)
library(ggtree)
library(ggplot2)
library(ggtreeExtra)
library(patchwork)
library(ggridges)
library(phyloseq)
library(reshape)
library(ggstar)
library(ggnewscale)
library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)
set.seed(1024)

# ------------------------------------------------------------------
# Kraken2 taxonomy annotation
# ------------------------------------------------------------------
ann <- read.csv("../../1.data.process/kk2.v2/kraken2.db.annotation.csv", header = TRUE)
ann1 <- ann[!duplicated(ann$genus), ]
ann2 <- ann[!duplicated(ann$species), ]
rownames(ann1) <- ann1$genus
rownames(ann2) <- ann2$species
ann  <- rbind(ann1, ann2)

# Remove suffixes from phylum names
a <- gsub("p__", "", ann$p)
a <- gsub("_[A-Z]{1,2}$", "", a)
ann$p <- paste0("p__", a)

# ------------------------------------------------------------------
# MaAsLin2 differential-results objects
# ------------------------------------------------------------------
res.location <- readRDS("MaAsLin.Location.rds")   # Left vs Right
res.tumor    <- readRDS("MaAsLin.Tumor_NAT.rds")  # Tumor vs NAT

# ==========================================================================
#  Fig 2c – Left vs Right (WGS, species level, circular heat-map)
# ==========================================================================

# ------------------------------------------------------------------
# Significant species common to UU Tumor & UU NAT (FDR < 0.05)
# ------------------------------------------------------------------
f1 <- res.location$UU.Tumor[res.location$UU.Tumor$FDR < 0.05, ]
f2 <- res.location$UU.NAT[res.location$UU.NAT$FDR   < 0.05, ]
f.location <- intersect(rownames(f1), rownames(f2))
f.location <- f.location[grep("s__", f.location)]

# Keep species with concordant effect directions
f1 <- res.location$UU.Tumor[f.location, ]
f2 <- res.location$UU.NAT[f.location, ]
f.location <- rownames(f1)[sign(f1$estimate) == sign(f2$estimate)]

# Build effect-size matrix (all cohorts)
qdat.location <- data.frame(
  row.names = f.location,
  UU_T = res.location$UU.Tumor[f.location, 3],
  UU_N = res.location$UU.NAT[f.location, 3],
  UM_T = res.location$UM.Tumor[f.location, 3],
  UM_N = res.location$UM.NAT[f.location, 3]
)

# Merge taxonomy and enrichment direction
ann1 <- ann[f.location, ]
ann1$Enrich  <- ifelse(qdat.location$UU_T > 0, "Right", "Left")
ann1$Enrich  <- factor(ann1$Enrich, levels = c("Right", "Left"))
ann1$estimate <- qdat.location$UU_T
ann1 <- ann1[order(ann1$species), ]

# Phylum-level
l <- c("p__Firmicutes", "p__Campylobacterota", "p__Bacteroidota",
       "p__Actinobacteriota", "p__Verrucomicrobiota",
       "p__Proteobacteria", "p__Fusobacteriota")
ann1$p[!ann1$p %in% l] <- "Others"
ann1$p <- factor(ann1$p, levels = c(l, "Others"))
ann1 <- ann1[order(ann1$p), ]
ann1 <- ann1[order(ann1$Enrich, decreasing = TRUE), ]

# Remove poorly annotated / CAG / UBA / numeric genera
ann1 <- ann1[grep("g__CAG-\\d+$",      ann1$genus, invert = TRUE), ]
ann1 <- ann1[grep("g__UBA\\d+$",       ann1$genus, invert = TRUE), ]
ann1 <- ann1[grep("g__MGYG-HGUT-\\d+$", ann1$genus, invert = TRUE), ]
ann1 <- ann1[grep("\\d",              ann1$genus, invert = TRUE), ]

qdat.location <- qdat.location[rownames(ann1), ]

# Pretty row-names for MGYG-HGUT species
a <- grep("MGYG-HGUT", rownames(ann1))
if (length(a)) {
  b <- gsub(".*-", "", rownames(ann1)[a])
  rownames(ann1)[a] <- paste0(gsub("g__", "s__", ann1$genus[a]), " (", b, ")")
  rownames(qdat.location) <- rownames(ann1)
}

# Final ordering: phylum → enrichment
ann1 <- ann1[order(ann1$p, decreasing = TRUE), ]
ann1 <- ann1[order(ann1$Enrich, decreasing = TRUE), ]
qdat.location <- qdat.location[rownames(ann1), ]

# ------------------------------------------------------------------
# Circular heat-map
# ------------------------------------------------------------------
label.col <- ifelse(ann1$Enrich == "Right", "#248888", "#d57b2c")
col_fun1  <- colorRamp2(c(-2, 0, 2.8), c("#F0D879", "white", "#248888"))
col_phylum <- c(
  p__Actinobacteriota  = "#7FC97F",
  p__Bacteroidota      = "#BEAED4",
  p__Campylobacterota  = "#666666",
  p__Others            = "#FFFF99",
  p__Firmicutes        = "#386CB0",
  p__Fusobacteriota    = "#FDC086",
  p__Proteobacteria    = "#BF5B17",
  p__Verrucomicrobiota = "#F0027F"
)

pdf("Fig2c.pdf", width = 8, height = 8)
circos.clear()
circos.par(gap.after = 10)
circos.heatmap(ann1[, "p", drop = FALSE], col = col_phylum, track.height = 0.04)
circos.heatmap(mat = qdat.location, col = col_fun1, cluster = FALSE,
               rownames.side = "inside", rownames.col = label.col,
               bg.border = "black", rownames.cex = 0.7, track.height = 0.2)
lgd <- Legend(title = "WGS", col_fun = col_fun1)
grid.draw(lgd)
dev.off()

# ==========================================================================
#  Fig 2e – Tumour vs NAT (WGS, species level, circular heat-map)
# ==========================================================================

# ------------------------------------------------------------------
# Significant species common to UU & UM tumours (FDR < 0.05 & p < 0.05)
# ------------------------------------------------------------------
f1 <- res.tumor$UU[res.tumor$UU$FDR < 0.05, ]
f2 <- res.tumor$UM[res.tumor$UM$pval < 0.05, ]
f.tumor <- intersect(rownames(f1), rownames(f2))
f.tumor <- f.tumor[grep("s__", f.tumor)]

# Concordant effect direction
f1 <- res.tumor$UU[f.tumor, ]
f2 <- res.tumor$UM[f.tumor, ]
f.tumor <- rownames(f1)[sign(f1$Effectsize) == sign(f2$Effectsize)]

# Effect-size matrix (all sub-cohorts)
qdat <- data.frame(
  row.names = f.tumor,
  UU   = res.tumor$UU[f.tumor, 3],
  UU_L = res.tumor$UU_L[f.tumor, 3],
  UU_R = res.tumor$UU_R[f.tumor, 3],
  UM   = res.tumor$UM[f.tumor, 3],
  UM_L = res.tumor$UM_L[f.tumor, 3],
  UM_R = res.tumor$UM_R[f.tumor, 3]
)

# Taxonomy & enrichment direction
ann2 <- ann[f.tumor, ]
ann2$Enrich <- ifelse(qdat$UU > 0, "Tumor", "NAT")
ann2$Enrich <- factor(ann2$Enrich, levels = c("Tumor", "NAT"))
ann2 <- ann2[order(ann2$species), ]

l <- c("p__Proteobacteria", "p__Fusobacteriota", "p__Firmicutes",
       "p__Campylobacterota", "p__Bacteroidota", "p__Actinobacteriota",
       "p__Verrucomicrobiota")
ann2$p[!ann2$p %in% l] <- "Others"
ann2$p <- factor(ann2$p, levels = c(l, "Others"))
ann2 <- ann2[order(ann2$p), ]
ann2 <- ann2[order(ann2$Enrich), ]

# Remove poorly annotated taxa (same filters as above)
ann2 <- ann2[grep("g__CAG-\\d+$",      ann2$genus, invert = TRUE), ]
ann2 <- ann2[grep("g__UBA\\d+$",       ann2$genus, invert = TRUE), ]
ann2 <- ann2[grep("g__MGYG-HGUT-\\d+$", ann2$genus, invert = TRUE), ]
ann2 <- ann2[grep("\\d",              ann2$genus, invert = TRUE), ]

qdat <- qdat[rownames(ann2), ]

# Pretty row-names for MGYG-HGUT species
a <- grep("MGYG-HGUT", rownames(ann2))
if (length(a)) {
  b <- gsub(".*-", "", rownames(ann2)[a])
  rownames(ann2)[a] <- paste0(gsub("g__", "s__", ann2$genus[a]), " (", b, ")")
  rownames(qdat) <- rownames(ann2)
}

# Final ordering
ann2 <- ann2[order(ann2$p, decreasing = TRUE), ]
ann2 <- ann2[order(ann2$Enrich, decreasing = TRUE), ]
qdat <- qdat[rownames(ann2), ]

# Add Left/Right information from Fig2c
ann2$Location <- ann1[rownames(ann2), "Enrich"]
ann2$Location <- factor(ann2$Location, levels = c("Right", "Left", "None"))
ann2$Location[is.na(ann2$Location)] <- "None"

# ------------------------------------------------------------------
# Circular heat-map
# ------------------------------------------------------------------
label.col <- ifelse(ann2$Enrich == "NAT", "#848CCF", "#BE5683")
col_fun1  <- colorRamp2(c(-2, 0, 2.8), c("#848CCF", "white", "#BE5683"))
col_dir   <- structure(c("#248888", "white", "#F0D879"),
                       names = c("Right", "None", "Left"))

pdf("Fig2e.pdf", width = 8, height = 8)
circos.clear()
circos.par(gap.after = 10)
circos.heatmap(ann2[, "Location", drop = FALSE], col = col_dir, track.height = 0.04)
circos.heatmap(ann2[, "p", drop = FALSE], col = col_phylum, track.height = 0.04)
circos.heatmap(mat = qdat, col = col_fun1, cluster = FALSE,
               rownames.side = "inside", rownames.col = label.col,
               bg.border = "black", track.height = 0.2, rownames.cex = 0.7)
lgd <- Legend(title = "WGS", col_fun = col_fun1)
grid.draw(lgd)
dev.off()