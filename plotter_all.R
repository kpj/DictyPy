library(ggplot2)
library(reshape)
library(grid)
library(plyr)


df <- read.csv("out.csv", header=TRUE)

ggplot(df, aes(x=amino_acid, y=codon_usage, fill=codon)) +
  geom_bar(stat="identity") +
  coord_flip() +
  facet_grid(group ~ .) +
  theme(panel.margin = unit(2, "lines"))
