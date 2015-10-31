library(ggplot2)
library(rjson)
library(reshape2)


# read data
fname <- "results/longest_stretches.json"
data <- fromJSON(file=fname, method='C')

# preprocess data
df <- data.frame(codon=numeric(0), length=numeric(0))
for(cur in data$aaa_len) {
  df <- rbind(df, data.frame(codon="AAA", len=cur))
}
for(cur in data$caa_len) {
  df <- rbind(df, data.frame(codon="CAA", len=cur))
}

# plot data
ggplot(df, aes(x=len, fill=..count..)) +
  geom_histogram(binwidth=2) +
  xlab("stretch length") +
  scale_y_log10() +
  facet_grid(codon ~ .) +
  ggsave(filename=paste0("stretch_hist.png"))
