suppressPackageStartupMessages({
    library(tidyverse)
    library(scales)
})

# read in variables
ivals      <- read.csv(snakemake@input[[1]])
len_barplt <- snakemake@output[['len_barplt']]

# rank interval length per chr
ivals <- ivals %>%
  group_by(chr) %>%
  mutate(position = rank(-interval_length),
         chr_no = str_remove(chr, "chr"))

# reorder chrs
ivals$chr <- factor(ivals$chr, levels = c(paste0("chr",seq(1, 38, by = 1)),"chrM","chrX"))

tiff(len_barplt,width = 1000, height = 1000)
ggplot(ivals, 
       aes(chr,
           interval_length,
           group = position,
           fill = chr)) +
  geom_bar(stat = "identity", position=position_dodge()) +
  scale_fill_manual(
    values=rep(c("#7A0019","#FFCC33"), ceiling(length(ivals$chr)/2))[1:length(ivals$chr)]) +
  scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6), 
                     expand = expansion(mult = c(0,0.1))) +
  theme_classic() +
  ylab("Interval length") +
  xlab("Chromosome") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    aspect.ratio=1
  )
dev.off()
