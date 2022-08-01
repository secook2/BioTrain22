library(data.table)

install.packages("ggplot2")
library(ggplot2)

wind <- fread("/Users/biotrain2022/Desktop/thingsForUAB/Coverage_Map_for_CHR9/line8_10kb_wind.bed", 
              col.names = c("chr", "pos", "depth"))

gag <- fread("/Users/biotrain2022/Desktop/thingsForUAB/Coverage_Map_for_CHR9/GAgreen_cov.txt",
             col.names = c("chr", "pos", "depth"))

C76 <- fread("/Users/biotrain2022/Desktop/thingsForUAB/Coverage_Map_for_CHR9/C76_cov.txt", 
             col.names = c("chr", "pos", "depth"))

card <- fread("/Users/biotrain2022/Desktop/thingsForUAB/Coverage_Map_for_CHR9/cardenasii_cov.txt",
              col.names = c("chr", "pos", "depth"))

line8 <- fread("/Users/biotrain2022/Desktop/thingsForUAB/Coverage_Map_for_CHR9/line8_cov.txt",
               col.names = c("chr", "pos", "depth"))

cardall <- merge(card, wind, by = c("chr", "pos"), all.y = T)
cardall[is.na(cardall)] <- 0


C76[, genome := "C76"]
card[, genome := "cardenasii"]
gag[, genome := "GAgreen"]
line8[, genome := "line8"]



chr9sgag <- subset(gag, grepl("09", chr))
chr9sC76 <- subset(C76, grepl("09", chr))
chr9sline8 <- subset(line8, grepl("09", chr))
chr9scard <- subset(card, grepl("09", chr))

chr9sgag$depth[chr9sgag$depth > 60] <- 60
chr9sC76$depth[chr9sC76$depth > 60] <- 60
chr9sline8$depth[chr9sline8$depth > 60] <- 60
chr9scard$depth[chr9scard$depth > 60] <- 60

with(chr9sC76, plot(x = pos, y = depth, pch = "."))
with(chr9sgag, plot(x = pos, y = depth, pch = "."))
with(chr9sline8, plot(x = pos, y = depth, pch = "."))
with(chr9scard, plot(x = pos, y = depth, pch = "."))

comb <- rbindlist(list(chr9sC76, chr9sgag, chr9scard, chr9sline8))

ggplot(comb, aes(x = pos, y =depth)) + 
  geom_point(pch = ".") +
  stat_smooth()+
  facet_grid(.~chr, scale = "free", space = "free")

ggplot(chr9sC76, aes(x = pos, y =depth)) + 
  geom_point(pch = ".") +
  stat_smooth()+
  facet_grid(.~chr, scale = "free", space = "free")


comb$depth[comb$depth > 150] <- 150
comb[,corDepth := depth/mean(depth, na.rm = T), by = "chr"]

pdf("/Users/biotrain2022/Desktop/thingsForUAB/Coverage_Map_for_CHR9/coverage.pdf", height = 6, width = 12)
coverage <- ggplot(comb[comb$genome %in% "C76"], aes(x = pos/1e6, col = corDepth, y = depth))+
  geom_point(pch = ".", color = "skyblue") +
  geom_point(comb[comb$genome %in% "line8"], mapping = aes(x = pos/1e6, col = corDepth, y = depth), 
             pch = ".", color = "purple") +
  geom_point(comb[comb$genome %in% "GAgreen"], mapping = aes(x = pos/1e6, col = corDepth, y = depth), 
             pch = ".", color = "red") +
  geom_point(comb[comb$genome %in% "cardenasii"], mapping = aes(x = pos/1e6, col = corDepth, y = depth), 
             pch = ".", color = "lightgreen") + 
  facet_grid(genome ~ chr, scale = "free", space = "free")+
  geom_smooth(data=comb[comb$genome %in% "C76"], aes(x = pos/1e6, col = corDepth, y = depth), se=F, color = "white") +
  geom_smooth(data=comb[comb$genome %in% "line8"], aes(x = pos/1e6, col = corDepth, y = depth), se=F, color = "white")+
  geom_smooth(data=comb[comb$genome %in% "GAgreen"], aes(x = pos/1e6, col = corDepth, y = depth), se=F, color = "white") +
  geom_smooth(data=comb[comb$genome %in% "cardenasii"], aes(x = pos/1e6, col = corDepth, y = depth), se=F, color = "white")+
  #stat_smooth(span = .2, color = "white") + 
  labs(x = "Physical position (Mbp)", y = "Coverage (X)",
       main = "CCS coverage") +
  theme(panel.background = element_rect(fill="black"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "white", size = .1, linetype = 2),
        panel.grid.minor.y = element_blank(),
        axis.title = element_text(size = 20, family = "Helvetica", color = "black"), 
        axis.text.x.bottom = element_text(size = 14, family = "Helvetica", color = "black"), 
        axis.text.y = element_text(size = 18, family = "Helvetica", color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "none")
dev.off()


ggsave(coverage, filename = '/Users/biotrain2022/Desktop/thingsForUAB/Coverage_Map_for_CHR9/coverage.png', width = 6, height = 12)



# needs modifications?
pdf("/Users/biotrain2022/Desktop/6.27/new.pdf", height = 6, width = 12)
ggplot(comb, aes(x = pos/1e6, col = corDepth, y = depth))+
  geom_point(pch = ".", color = "skyblue")+
  facet_grid(genome ~ chr, scale = "free", space = "free")+
  stat_smooth(span = .2, color = "white") + 
  #scale_y_log10()+
  labs(x = "Physical position (Mbp)", y = "Coverage (X)",
       main = "CCS coverage") +
  theme(panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "white", size = .1, linetype = 2),
        panel.grid.minor.y = element_blank(),
        axis.title = element_text(size = 20, family = "Helvetica", color = "white"), 
        axis.text.x.bottom = element_text(size = 14, family = "Helvetica", color = "white"), 
        axis.text.y = element_text(size = 18, family = "Helvetica", color = "white"),
        legend.position = "none")
dev.off()
