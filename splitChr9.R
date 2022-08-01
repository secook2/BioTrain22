install.packages("BiocManager")
BiocManager::install("Biostrings")
library(Biostrings)
install.packages("ggplot2")
library(ggplot2)
install.packages("dbscan")
library(dbscan)
library(data.table)

# reading in the chr9 split file
d <- fread("/Users/biotrain2022/Desktop/things for powerpoint and paper/Best_Map/chr9_split_asm5p1.paf",
           fill = T)

# creating a new column for tmp it is splitting a character vector in this
# case it is splitting a binary file
d[,tmp := as.numeric(d[,tstrsplit(V1, "_")][[2]])]

# creating new columns for line 8 start cr start and end
d[,`:=`(l8start = V8/1e6, l8end = V9/1e6, 
        crstart = (tmp + V3)/1e6,
        crend = (tmp + V4)/1e6)]

# creating a new data set to be able to actually see what each column shows
dproc <- with(d, data.table(line8chr = V6, line8pos = (l8end + l8start)/2,
  cardpos = (crend + crstart)/2,
  score = V12, percId = V10/V11, alignLength = V11))


# subsets the dproc on the perc score alignLength
tp <- subset(dproc, percId > .9 & score >= 55 & alignLength > 
               8e3 & alignLength < 12e3)

# creates new columns for the ords and Returns 
# the sample ranks of the values in a vector
tp[,`:=`(ord1 = frank(line8pos, ties.method = "dense"), 
         ord2 = frank(cardpos, ties.method = "dense")),
   by = "line8chr"]


# creates new column for cluster and estimates the density around each 
# data point
tp[,clus := dbscan(frNN(cbind(ord1, ord2), eps = sqrt(2)*10), 
                   minPts = 10)$cluster, by = "line8chr"]

# subsets tp by the clusters and sets the x and y coloring based on line 8 chr
pdf("/Users/biotrain2022/Desktop/things for powerpoint and paper/Best_Map/bestMap.pdf", height = 6, width = 8)
# i still need help with this one unless this is good
ggplot(subset(tp, clus > 0), 
       aes(x = cardpos, y = line8pos, col = line8chr))+
  # creates scatter plot with points of .5
  geom_point(cex = .5) +
  # labels the x and y axis and key by the best map
  labs(x = "Cardenasii Genome Position (Mbp)", 
        y = "Line8 Genome Position (Mbp)", col = "Best Map") +
  #axis(1, )
  coord_fixed(ratio = 1,  xlim = c(0,130), ylim = c(0,130), expand = TRUE) +
  # themes the plot
  theme(panel.background = element_rect(fill='white'), 
        axis.title = element_text(size = 20, family = "Helvetica", color = "black"), 
        axis.text.x.bottom = element_text(size = 18, family = "Helvetica", color = "black"), 
        axis.text.y = element_text(size = 18, family = "Helvetica", color = "black"),
        panel.grid.major.x = element_line(color = "black", size = .2,
                                          linetype = 2),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "black",
                                          size = .2, linetype = 2),
        panel.grid.minor.y = element_blank(),
        legend.position = "none") +
  scale_x_continuous(breaks = seq(0,125, by = 25)) + 
  scale_y_continuous(breaks = seq(0,125, by = 25)) +
  scale_color_manual(values = c("Arahy.09" = "red", "Arahy.09_alt" = "skyblue"))
dev.off()



