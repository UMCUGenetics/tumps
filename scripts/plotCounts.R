library(ggplot2)
library(reshape2)

setwd("/hpc/cog_bioinf/cuppen/project_data/Jose_COLO829/shared_COLO829/subsampling_purity/")
df <- read.table("subsampling_table.txt", header = T, sep = '\t')
df$purity <- gsub("purity", "", df$purity)
df$purity <- factor(df$purity, 
                    levels = c('0','10', '20', '25', '50', '75', '100'))
df$group <- paste(df$tech, df$type)
medians <- aggregate(list(df$value), by=list(df$tech, df$purity, df$type, df$group), median)
colnames(medians) <- c("tech", "purity", "type", "group", "value")

#Original
p <- ggplot(df, aes(x=purity, y=value, group = group)) + 
    geom_jitter(aes(shape=type, col=tech), fill = NA, size = 3, height=0, width=0.05) +
    geom_line(data=medians, aes(x=purity, y=value, col=tech, linetype = type)) +
    scale_color_manual(values = c("darkolivegreen3", "deepskyblue4", "firebrick3"), 
                       limits =c("illumina", "nanopore", "pacbio")) +
    scale_shape_manual(values = c(21, 19), limits = c("raw", "somatic")) +
    scale_linetype_manual(values = c("dashed", "solid"), limits = c("raw", "somatic")) +
    scale_y_continuous(breaks = seq(0,100, 10)) +
    xlab("% tumor purity") + ylab("% recall") +
    theme_bw() +
    theme(
        text = element_text(face = 'bold', size = 15)
    )

#p2 only somatic
p2 <- ggplot(df[df$type == "somatic",], aes(x=purity, y=value, group = group)) +
    geom_jitter(aes(shape=type, col=tech), fill = NA, size = 3, height=0, width=0.05) +
    geom_line(data=medians[medians$type == "somatic",], aes(x=purity, y=value, col=tech, linetype = type)) +
    scale_color_manual(values = c("darkolivegreen3", "deepskyblue4", "firebrick3"),
                       limits =c("illumina", "nanopore", "pacbio")) +
    scale_shape_manual(values = c(21, 19), limits = c("raw", "somatic"), guide=F) +
    scale_linetype_manual(values = c("dashed", "solid"), limits = c("raw", "somatic"), guide=F) +
    scale_y_continuous(breaks = seq(0,100, 10)) +
    xlab("% tumor purity") + ylab("% recall") +
    theme_bw() +
    theme(
        text = element_text(face = 'bold', size = 15)
    )

#p3 facet_grid
p3 <- ggplot(df, aes(x=purity, y=value, group = group)) +
    facet_grid(. ~ tech) +
    geom_jitter(aes(shape=type, col=tech), fill = NA, size = 3, height=0, width=0.05) +
    geom_line(data=medians, aes(x=purity, y=value, col=tech, linetype = type)) +
    scale_color_manual(values = c("darkolivegreen3", "deepskyblue4", "firebrick3"),
                       limits =c("illumina", "nanopore", "pacbio")) +
    scale_shape_manual(values = c(21, 19), limits = c("raw", "somatic")) +
    scale_linetype_manual(values = c("dashed", "solid"), limits = c("raw", "somatic")) +
    scale_y_continuous(breaks = seq(0,100, 10)) +
    xlab("% tumor purity") + ylab("% recall") +
    theme_bw() +
    theme(
        text = element_text(face = 'bold', size = 15)
    )

#p4 facet_grid: only somatic

p4 <- ggplot(df[df$type == "somatic",], aes(x=purity, y=value, group = group)) +
    facet_grid(. ~ tech) +
    geom_jitter(aes(shape=type, col=tech), fill = NA, size = 3, height=0, width=0.05) +
    geom_line(data=medians[medians$type == "somatic",], aes(x=purity, y=value, col=tech, linetype = type)) +
    scale_color_manual(values = c("darkolivegreen3", "deepskyblue4", "firebrick3"),
                       limits =c("illumina", "nanopore", "pacbio")) +
    scale_shape_manual(values = c(21, 19), limits = c("raw", "somatic"), guide=F) +
    scale_linetype_manual(values = c("dashed", "solid"), limits = c("raw", "somatic"), guide=F) +
    scale_y_continuous(breaks = seq(0,100, 10)) +
    xlab("% tumor purity") + ylab("% recall") +
    theme_bw() +
    theme(
        text = element_text(face = 'bold', size = 15)
    )
 


pdf("subsampling_plot.pdf", width = 16, height = 8)
p
p2
p3
p4
dev.off()
