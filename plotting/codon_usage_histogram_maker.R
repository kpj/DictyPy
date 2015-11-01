library(ggplot2)
library(rjson)
library(reshape2)


# read data
fname <- "results/gene_codon_usages.json"
data <- fromJSON(file=fname, method='C')

# preprocess data
df <- data.frame(marker=numeric(0), counts=numeric(0), edges=numeric(0))
for(cur in data) {
  df <- rbind(df, data.frame(marker=cur$marker, counts=cur$counts, edges=cur$edges))
}

# plot data
ggplot(df, aes(x=edges, y=counts, fill=counts)) +
  geom_bar(stat="identity") +
  geom_text(aes(label=ifelse(counts == 0, "", counts), angle=90, hjust=-0.25), size=3) +
  xlab("codon usage") +
  facet_grid(marker ~ .) +
  ggsave(filename=paste0("codon_usage_hist.png"), width=13)
