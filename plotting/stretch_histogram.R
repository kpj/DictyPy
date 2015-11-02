library(ggplot2)
library(rjson)
library(reshape2)


# read data
fname <- "results/longest_stretches.json"
data <- fromJSON(file=fname, method='C')

# preprocess data
df <- data.frame(codon=numeric(0), counts=numeric(0), edges=numeric(0))
for(cur in data) {
  df <- rbind(df, data.frame(codon=cur$codon, counts=cur$counts, edges=cur$edges))
}

# plot data
ggplot(df, aes(x=edges, y=counts, fill=counts)) +
  geom_bar(stat="identity") +
  geom_text(aes(label=ifelse(counts == 0, "", counts), angle=90, hjust=-0.25), size=3) +
  xlab("stretch length") +
  facet_grid(codon ~ .) +
  ggsave(filename=paste0("stretch_hist.png"), width=13)
