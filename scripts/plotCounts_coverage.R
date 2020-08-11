library(ggplot2)
library(reshape2)
library(ggsci)

df <- read.table("coverage_table.tsv", header = T, sep = '\t')

###ADD CALCULATIONS###
df$precision <- (df$tp) / (df$tp + df$fp)
df$recall <- (df$tp) / (df$tp + df$fn)

df$purity <- gsub("purity", "", df$purity)
df$purity <- factor(df$purity, 
                    levels = c('0','10', '20', '25', '50', '75', '100'))
df$coverage <- gsub("cov", "", df$coverage)
df$coverage <- factor(df$coverage, 
                    levels = c('1', '5', '10', '30', '50', '98'))
df[is.na(df)] <- 0

dfm <- melt(df, id.vars = c("tech", "coverage", "purity", "type", "tp", "fp", "fn"))
dfm$group <- paste(dfm$coverage, dfm$type)
medians <- aggregate(list(dfm$value), 
                     by=list(dfm$tech, dfm$coverage, dfm$purity, dfm$type, dfm$group, dfm$variable), median)
colnames(medians) <- c("tech", "coverage", "purity", "type", "group", "variable", "value")
dfm$variable <- factor(dfm$variable, levels=c("recall", "precision"))

som <- dfm[dfm$type == "somatic",]
medsom <- medians[medians$type == "somatic",]
psom <- ggplot(som, aes(x=purity, y=value, group = coverage)) + 
  facet_grid(. ~ variable) +
  geom_jitter(aes(col=coverage), fill = NA, size = 3, height=0, width=0.05, alpha=.9) +
  geom_line(data=medsom, aes(x=purity, y=value, col=coverage)) +
  scale_color_npg() + 
  scale_y_continuous(breaks = seq(0,1, 0.1), limits = c(0,1), labels = seq(0,100,10)) +
  xlab("Purity") + ylab("%") +
  theme_bw() +
  theme(
    text = element_text(face = 'bold', size = 15)
  )
psom

precall <-  ggplot(som[som$variable == "recall", ], 
                   aes(x=purity, y=value, group = coverage)) + 
  geom_jitter(aes(col=coverage), fill = NA, size = 3, height=0, width=0.05, alpha=.9) +
  
  geom_line(data=medsom[medsom$variable == "recall", ], 
            aes(x=purity, y=value, col=coverage), show.legend = F) +
  scale_color_npg() + 
  scale_y_continuous(breaks = seq(0,1, 0.1), limits = c(0,1), labels = seq(0,100,10)) +
  xlab("Purity") + ylab("% recall") +
  theme_bw() +
  theme(
    text = element_text(face = 'bold', size = 15)
  )
precall

all <- ggplot(dfm, aes(x=purity, y=value, group = group)) + 
  facet_grid(. ~ variable) +
  geom_jitter(aes(col=coverage, shape = type), fill = NA, size = 3, height=0, width=0.05, alpha=.9) +
  geom_line(data=medians, aes(x=purity, y=value, col=coverage, linetype = type)) +
  scale_color_npg() + 
  scale_shape_manual(values = c(21, 19), limits = c("raw", "somatic")) +
  scale_y_continuous(breaks = seq(0,1, 0.1), limits = c(0,1), labels = seq(0,100,10)) +
  scale_linetype_manual(values = c("dashed", "solid"), limits = c("raw", "somatic")) +
  xlab("Purity") + ylab("%") +
  theme_bw() +
  theme(
    text = element_text(face = 'bold', size = 15)
  )
all
pdf("colo_coverage.pdf", width = 8, height = 6)
psom
precall
all
dev.off()
write.csv(df, file = "coverage_calculations.csv")
write.csv(medians, file = "coverage_medians.csv")
