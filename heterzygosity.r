install.packages("data.table")
library(data.table)
install.packages('R.utils')
library(R.utils)
install.packages("ggplot2")
library(ggplot2)
install.packages("scales")
library("scales")




##Sites - O(HOM) = O(HET)
# this is the observed and expected homozygosity i want heterozygosity.  
hom <- fread("/Users/stephaniecook/Desktop/HudsonAlpha/HetBar/peanut.line8.d8.poly.snp.vcf.het")


# based on https://www.biostars.org/p/266502/
hom[,"E(HET)" := hom$N_SITES - hom$`E(HOM)`]
hom[,"O(HET)" := hom$N_SITES - hom$`O(HOM)`]

hom[, "H0" := hom$`O(HOM)`/hom$N_SITES]


# JIGT_PCRfree_NA_NA_ATGTAAGT_Arachis_hypogaea_Line8
# JPQS_PCRfree_Leaf_126_CCGCGGTT_Arachis_P-1_PI502120
# JPQT_PCRfree_Leaf_127_TTATAACC_Arachis_P-3_AU17
# JPQU_PCRfree_Leaf_128_GGACTTGG_Arachis_P-5_Georgia_Green
# JQAA_PCRfree_Leaf_129_AAGTCCAA_Arachis_P-7_AP3
# JQAB_PCRfree_Leaf_130_ATCCACTG_Arachis_P-9_587
# JQAE_PCRfree_Leaf_131_GCTTGTCA_Arachis_P-11_C76_16
# JQAF_PCRfree_Leaf_132_CAAGCTAG_Arachis_P-13_AT3085RO
# JQIH_PCRfree_Leaf_134_AGTTCAGG_Arachis_P-17_Tifrunner


pdf("/Users/biotrain2022/Desktop/thingsForUAB/HetBar/comparedHet.pdf", height = 6, width = 12)
ggplot(hom, aes(x = INDV, y = `O(HET)`, fill = INDV))+geom_bar(stat = "identity") +
  theme(panel.background = element_rect(fill="black"),
        axis.line = element_line(color = "white"),
        panel.grid = element_blank(),
        axis.ticks = element_line(color = "white"),
        axis.text = element_text(family = "Helvetica", size = 18, color = "black"),
        axis.title = element_text(family = "Helvetica", size = 20, color = "black"), 
        legend.position = "none") +
  scale_fill_manual(values = c("purple", rep("lightgrey", 2), "red", 
                                rep("lightgrey", 2),"skyblue", rep("lightgrey", 2))) +
  labs(x = "Library", y = "Observed Heterozygosity") +
  scale_x_discrete(labels = c("Line8", "PI502120", "AU17", 
                                "GAGreen", "AP3", "587", "C76", "AT3085RO", "Tiffrunner"))
dev.off()



