library(ggplot2)
library(reshape2)
library(ggsci)

#setwd("/hpc/cog_bioinf/cuppen/project_data/Jose_COLO829/shared_COLO829/subsampling_purity/")
#df <- read.table("subsampling_table.txt", header = T, sep = '\t')

setwd("~/work/colo/")
df <- read.table("table_recall.txt", header = T, sep = '\t')

###ADD CALCULATIONS###
df$precision <- (df$tp) / (df$tp + df$fp)
df$recall <- (df$tp) / (df$tp + df$fn)
#df$f1 <- 2 * ((df$precision*df$recall)/(df$precision+df$recall))

df$purity <- gsub("purity", "", df$purity)
df$purity <- factor(df$purity, 
                    levels = c('0','10', '20', '25', '50', '75', '100'))
#df$group <- paste(df$tech, df$type)
df[is.na(df)] <- 0

dfm <- melt(df, id.vars = c("tech", "n", "purity", "type", "tp", "fp", "fn"))
dfm$group <- paste(dfm$tech, dfm$type)
medians <- aggregate(list(dfm$value), by=list(dfm$tech, dfm$purity, dfm$type, dfm$group, dfm$variable), median)
colnames(medians) <- c("tech", "purity", "type", "group", "variable", "value")
#dfm$variable <- factor(dfm$variable, levels=c("recall", "precision", "f1"))
dfm$variable <- factor(dfm$variable, levels=c("recall", "precision", "f1"))
# 
# p <- ggplot(dfm, aes(x=purity, y=value, group = group)) + 
#   facet_grid(. ~ variable) +
#   geom_jitter(aes(shape=type, col=tech), fill = NA, size = 3, height=0, width=0.05, alpha=.6) +
#   geom_line(data=medians, aes(x=purity, y=value, col=tech, linetype = type)) +
#   scale_color_manual(values = c("darkolivegreen3", "deepskyblue4", "firebrick3"), 
#                      limits =c("illumina", "nanopore", "pacbio")) +
#   scale_shape_manual(values = c(21, 19), limits = c("raw", "somatic")) +
#   scale_linetype_manual(values = c("dashed", "solid"), limits = c("raw", "somatic")) +
#   scale_y_continuous(breaks = seq(0,1, 0.1)) +
#   #scale_x_continuous(breaks = seq(0,1, 0.1), limits = c(-0.000001,1.000001)) +  
#   xlab("Purity") + ylab("") +
#   theme_bw() +
#   theme(
#     text = element_text(face = 'bold', size = 15)
#   )
# p

som <- dfm[dfm$type == "somatic",]
medsom <- medians[medians$type == "somatic",]
psom <- ggplot(som, aes(x=purity, y=value, group = group)) + 
  facet_grid(. ~ variable) +
  geom_jitter(aes(col=tech), fill = NA, size = 3, height=0, width=0.05, alpha=.9) +
  geom_line(data=medsom, aes(x=purity, y=value, col=tech)) +
  # scale_color_npg() + 
  scale_color_manual(values = c("springgreen3", "deepskyblue4", "firebrick3"), 
                    limits =c("illumina", "nanopore", "pacbio")) +
  #scale_shape_manual(values = c(21, 19), limits = c("raw", "somatic")) +
  scale_y_continuous(breaks = seq(0,1, 0.1), limits = c(0,1), labels = seq(0,100,10)) +
  #scale_x_continuous(breaks = seq(0,1, 0.1), limits = c(-0.000001,1.000001)) +  
  xlab("Purity") + ylab("%") +
  theme_bw() +
  theme(
    text = element_text(face = 'bold', size = 15)
  )
psom

precall <-  ggplot(som[som$variable == "recall", ], aes(x=purity, y=value, group = group)) + 
  geom_jitter(aes(col=tech), fill = NA, size = 3, height=0, width=0.05, alpha=.9) +
  geom_line(data=medsom[medsom$variable == "recall", ], aes(x=purity, y=value, col=tech)) +
  scale_color_manual(values = c("springgreen3", "deepskyblue4", "firebrick3"), 
                      limits =c("illumina", "nanopore", "pacbio")) +
  # scale_color_npg() + 
  #scale_shape_manual(values = c(21, 19), limits = c("raw", "somatic")) +
  scale_y_continuous(breaks = seq(0,1, 0.1), limits = c(0,1), labels = seq(0,100,10)) +
  #scale_x_continuous(breaks = seq(0,1, 0.1), limits = c(-0.000001,1.000001)) +  
  xlab("Purity") + ylab("% recall") +
  theme_bw() +
  theme(
    text = element_text(face = 'bold', size = 15)
  )
precall
all <- ggplot(dfm, aes(x=purity, y=value, group = group)) + 
  facet_grid(. ~ variable) +
  geom_jitter(aes(col=tech, shape = type), fill = NA, size = 3, height=0, width=0.05, alpha=.9) +
  geom_line(data=medians, aes(x=purity, y=value, col=tech, linetype = type)) +
  # scale_color_npg() + 
  scale_color_manual(values = c("springgreen3", "deepskyblue4", "firebrick3"), 
                      limits =c("illumina", "nanopore", "pacbio")) +
  scale_shape_manual(values = c(21, 19), limits = c("raw", "somatic")) +
  scale_y_continuous(breaks = seq(0,1, 0.1), limits = c(0,1), labels = seq(0,100,10)) +
  scale_linetype_manual(values = c("dashed", "solid"), limits = c("raw", "somatic")) +
  #scale_x_continuous(breaks = seq(0,1, 0.1), limits = c(-0.000001,1.000001)) +  
  xlab("Purity") + ylab("%") +
  theme_bw() +
  theme(
    text = element_text(face = 'bold', size = 15)
  )
all

pdf("colo_purity.pdf", width = 8, height = 6)
#p
psom
precall
dev.off()
write.csv(df, file = "purity_calculations.csv")
write.csv(medians, file = "purity_medians.csv")

  

