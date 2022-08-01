library(data.table)
install.packages("ggplot2")
library(ggplot2)
install.packages("dplyr")
library(dplyr)
# sliding window analysis of heterozygosity in peanut

# make windows
make_windows <- function(step = 2e5, windowSize = 4e6, minWindSize = 2e6, fai){
  winds <- with(fai[1:20,], data.table(chr = V1, length = V2))
  winds <- winds[,list(start = seq(from = 1, to = length, by = step)), by = c("chr", "length")]
  winds[,end := start + windowSize - 1]
  winds$end[winds$end > winds$length] <- winds$length[winds$end > winds$length]
  winds[,`:=`(length = NULL, width = end - start + 1)]
  winds <- subset(winds, width > minWindSize)

  setkey(winds, chr, start, end)
  return(winds)
}


# path to vcf =  "/Users/stephaniecook/Desktop/HudsonAlpha/hetLines/peanut.line8.d8.poly.snp.vcf.gz"

metad <- as.data.table(tstrsplit(strsplit("JIGT_PCRfree_NA_NA_ATGTAAGT_Arachis_hypogaea_Line8
JPQS_PCRfree_Leaf_126_CCGCGGTT_Arachis_P-1_PI502120
JPQT_PCRfree_Leaf_127_TTATAACC_Arachis_P-3_AU17
JPQU_PCRfree_Leaf_128_GGACTTGG_Arachis_P-5_GAGreen
JQAA_PCRfree_Leaf_129_AAGTCCAA_Arachis_P-7_AP3
JQAB_PCRfree_Leaf_130_ATCCACTG_Arachis_P-9_587
JQAE_PCRfree_Leaf_131_GCTTGTCA_Arachis_P-11_C76
JQAF_PCRfree_Leaf_132_CAAGCTAG_Arachis_P-13_AT3085RO
JQIH_PCRfree_Leaf_134_AGTTCAGG_Arachis_P-17_Tifrunner", "\n")[[1]],"_"))[, c(1,8)]

setnames(metad, c("lib", "ID"))



f <- "/Users/stephaniecook/Desktop/HudsonAlpha/hetLines/peanut.line8.d8.poly.snp.vcf.gz"
raw <- fread(f, select = c("#CHROM", "POS", metad$lib), 
               col.names = c("chr", "pos", metad$ID))

fai <- fread("/Users/stephaniecook/Desktop/HudsonAlpha/hetLines/line8genome_bothChr9s.fa.fai")

long<- melt(raw, id.vars = c("chr", "pos"), variable.name = "ID", value.name = "geno")
long <- subset(long, geno == "0/1")

hets <- with(long, data.table(chr = chr, start = pos, end = pos, ID = ID))
setkey(hets,  chr, start, end)

winds <- make_windows(step = 1e5, windowSize = 2e6, minWindSize = 1e6, fai)

fo <- foverlaps(winds, hets)
foc <- fo[,list(nHets = sum(!is.na(start))), by = c("ID", "chr", "i.start", "i.end", "width")]
foc[,H0 := nHets / width]


foc[,chrn := gsub("Arahy.", "", chr)]
foc[, subg := ifelse(as.numeric(chrn)> 10, "B", "A")]
#foc[, chrn2:= ifelse(subg == "B", as.character(as.numeric(chrn)-10), chrn)]
foc[, chrnn := ifelse(as.numeric(chrn)> 10, as.numeric(chrn) -10, as.numeric(chrn))]

#foc[, chrnn := ifelse(subg == "B", as.numeric(chrn)10, chrn)]


foc[, medH0:= median(H0, na.rm = TRUE), by = c("chr", "i.start", "i.end")]
foc[, difH0 := H0 - medH0]

foc[, tr := ifelse(difH0 <= 0, 0, difH0)]

foc[, bin:= ifelse(tr < 0.0001, 0, 
                   ifelse(tr < 0.0002, 1, ifelse(tr < 0.0005, 2, 
                                                 ifelse(tr < 0.002, 3, 4))))]


focLine8 <- foc[foc$ID %in% "Line8"]
foc587 <- foc[foc$ID %in% "587"]
focAP3 <- foc[foc$ID %in% "AP3"]
focPI502120 <- foc[foc$ID %in% "PI502120"]
focC76 <- foc[foc$ID %in% "C76"]
focAT3085RO <- foc[foc$ID %in% "AT3085RO"]
focAU17 <- foc[foc$ID %in% "AU17"]
focGAGreen <- foc[foc$ID %in% "GAGreen"]
focTifrunner <- foc[foc$ID %in% "Tifrunner"]


# happy with this one, so far. the only thing i dont like is the x axis
pdf("/Users/biotrain2022/Desktop/thingsForUAB/hetLines/allHetdiffAltered2.pdf", height = 8, width = 16)


ggplot(focLine8, aes(x = i.start, y = tr, col = ID))+geom_line()+
  geom_line(mapping = aes(x = i.start, y = tr, col = ID), focTifrunner) +
  geom_line(mapping = aes(x = i.start, y = tr, col = ID), focAP3) +
  geom_line(mapping = aes(x = i.start, y = tr, col = ID), foc587) +
  geom_line(mapping = aes(x = i.start, y = tr, col = ID), focAT3085RO) +
  geom_line(mapping = aes(x = i.start, y = tr, col = ID), focAU17) +
  geom_line(mapping = aes(x = i.start, y = tr, col = ID), focC76) +
  geom_line(mapping = aes(x = i.start, y = tr, col = ID), focGAGreen) +
  geom_line(mapping = aes(x = i.start, y = tr, col = ID), focPI502120) +
  scale_color_manual(values = c(rep("grey", 6), "purple", "grey","darkgrey"))+
  theme(panel.background = element_rect(fill = "black"), panel.grid.minor = element_blank(), 
        panel.grid.major = element_line(size = .01, color = rgb(1,1,1,.1)), 
        axis.title = element_text(size = 20, family = "Helvetica"), 
        axis.text.x.bottom = element_text(size = 18, family = "Helvetica", angle = 90), 
        axis.text.y = element_text(size = 18, family = "Helvetica"),
        legend.background = element_rect(fill = "darkgrey"),
        legend.title = element_text(size = 12, family = "Helvetica")) + 
  labs(x = "Positions", y = "Heterozygosity")+
  facet_wrap(subg~chrnn, scale = "fixed", nrow = 2, ncol = 10)

dev.off()

ggsave(altered, filename = '/Users/biotrain2022/Desktop/thingsForUAB/hetLines/allHetdiffAltered2.png', width = 16, height = 8)
  



# pdf("~/Desktop/allHetdiffAltered.pdf", height = 12, width = 24)
# ggplot(focLine8, aes(x = i.start, y = tr, col = id))+geom_line()+
#   geom_line(mapping = aes(x = i.start, y = tr, col = id), focTifrunner) +
#   geom_line(mapping = aes(x = i.start, y = tr, col = id), focAP3) +
#   geom_line(mapping = aes(x = i.start, y = tr, col = id), foc587) +
#   geom_line(mapping = aes(x = i.start, y = tr, col = id), focAT3085RO) +
#   geom_line(mapping = aes(x = i.start, y = tr, col = id), focAU17) +
#   geom_line(mapping = aes(x = i.start, y = tr, col = id), focC76) +
#   geom_line(mapping = aes(x = i.start, y = tr, col = id), focGAGreen) +
#   geom_line(mapping = aes(x = i.start, y = tr, col = id), focPI502120) +
#   scale_color_manual(values = c(rep("grey", 6), "pink", "grey","darkgrey"))+
#   theme(panel.background = element_rect(fill = "black"), panel.grid.minor = element_blank(), 
#         panel.grid.major = element_line(size = .01, color = rgb(1,1,1,.1)),
#         axis.title = element_text(size = 20, family = "Helvetica")) + 
#   labs(x = "Positions", y = "Heterozygosity") + 
#   facet_grid(subg~chrnn, scale = "free")
# dev.off()

# # playing with this one
# pdf("~/Desktop/allHetdiffAltered.pdf", height = 12, width = 24)
# ggplot(foc, aes(x = i.start, y = tr, col = H0))+geom_point(pch = ".")+
#   stat_smooth()+
#   facet_grid(id~chr, scale = "free", space = "free")
# dev.off()



# pdf("~/Desktop/allHetBin.pdf", height = 12, width = 24)
# ggplot(foc, aes(x = i.start, y = bin, col = H0))+geom_jitter(pch = ".", height = 0.2)+
#   facet_grid(id~chr, scale = "free", space = "free")
# dev.off()
# 
# pdf("~/Desktop/allHetBin.pdf", height = 4, width = 12)
# ggplot(foc, aes(x = i.start, y = id, col = as.factor(bin)))+geom_jitter(pch = ".", height = 0.2)+
#    facet_grid(subg~chrnn, scale = "free", space = "free")
#  dev.off()
# 
# pdf("~/Desktop/allHet3.pdf", height = 12, width = 24)
# ggplot(foc, aes(x = i.start, y = log10(H0), col = H0))+geom_point(pch = ".")+
#   facet_grid(id~chr, scale = "free", space = "free")
# dev.off()
# 
# 
# pdf("~/Desktop/allHetdiff.pdf", height = 12, width = 24)
# ggplot(foc, aes(x = i.start, y = difH0, col = H0))+geom_point(pch = ".")+
#   stat_smooth()+
#   facet_grid(id~chr, scale = "free", space = "free")
# dev.off()
# 
# 
# pdf("~/Desktop/allHetdiff.pdf", height = 12, width = 24)
# ggplot(foc, aes(x = i.start, y = log10(tr + 1), col = H0))+geom_point(pch = ".")+
#   facet_grid(id~chr, scale = "free", space = "free")
# dev.off()
