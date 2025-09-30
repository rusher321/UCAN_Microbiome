library(vegan)
library(openxlsx)
library(ggpubr)
library(ggbeeswarm)

# ------------------------------------------------------------------
# Figure 2 – microbiota diversity and tissue-resident reads
# ------------------------------------------------------------------
load("Fig2.data.rd")          # pre-compiled list: dat.wgs, etc.

# ------------------------------------------------------------------
# Generic plotting function: boxplot + beeswarm + facet by cohort
# ------------------------------------------------------------------
plt <- function(dat, x, y, color){
  ggboxplot(dat, x = x, y = y, color = x, notch = TRUE, outlier.shape = NA) +
    geom_beeswarm(aes_string(color = x), alpha = 0.1, size = 1, dodge.width = 0.8) +
    facet_grid(~ cohort) +
    stat_compare_means(comparisons = list(c(1, 2))) +
    scale_y_continuous(expand = expansion(mult = c(.1, .1)), limits = c(0, 4.5)) +
    scale_color_manual(values = color) +
    theme_bw() +
    theme(
      axis.title  = element_text(color = 1),
      axis.text   = element_text(color = 1),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
      panel.grid  = element_blank(),
      strip.background = element_rect(color = "white", linewidth = 0.1),
      legend.position = ""
    )
}

# ------------------------------------------------------------------
# Fig 2a-b  Right- vs Left-sided tumours  (WGS)
# ------------------------------------------------------------------
a <- dat.wgs[dat.wgs$sampleType == "Tumor", ]
p1 <- plt(a, x = "Location", y = "Shannon", color = c("#248888", "#f0d879"))
p2 <- plt(a, x = "Location", y = "TissueResident", color = c("#248888", "#f0d879")) +
  scale_y_log10(limits = c(100, 5e8))

# ------------------------------------------------------------------
# Fig 2c-d  Right- vs Left-sided NAT  (WGS)
# ------------------------------------------------------------------
a <- dat.wgs[dat.wgs$sampleType == "NAT", ]
p3 <- plt(a, x = "Location", y = "Shannon",   color = c("#248888", "#f0d879"))
p4 <- plt(a, x = "Location", y = "TissueResident", color = c("#248888", "#f0d879")) +
  scale_y_log10(limits = c(100, 5e8))

# ------------------------------------------------------------------
# Fig 2e-g  Tumour vs NAT  (paired patients, all / right / left)
# ------------------------------------------------------------------
a <- dat.wgs[dat.wgs$PatientID %in% dat.wgs$PatientID[dat.wgs$sampleType == "NAT"], ]
p5 <- plt(a, x = "sampleType", y = "Shannon", color = c("#848ccf", "#be5683"))
p6 <- plt(a[a$Location == "Right", ], x = "sampleType", y = "Shannon", color = c("#848ccf", "#be5683"))
p7 <- plt(a[a$Location == "Left", ],  x = "sampleType", y = "Shannon", color = c("#848ccf", "#be5683"))

# ------------------------------------------------------------------
# save
# ------------------------------------------------------------------
ggarrange(p2, p4, p1, p3, p5, p6, p7,
          nrow = 1, ncol = 7, align = "hv")
ggsave("Fig2abd.boxplot.Shannon.pdf", width = 14, height = 3.3)