library(ggplot2)
library(rjson)
library(reshape2)


# read data
fname <- "results/longest_stretches.json"
data <- fromJSON(file=fname, method='C')

# preprocess data
df <- data.frame()
for(cur in data) {
  for(entry in cur$data) {
    df <- rbind(df, data.frame(codon=cur$codon, x_pos=entry$x, y_pos=entry$y, count=entry$z, a=ifelse(entry$z == 0, 0, 1)))
  }
}

# plot data
ggplot(df, aes(x=x_pos, y=y_pos)) +
  geom_raster(aes(fill=ifelse(count == 0, 0, log(count)))) + # alpha=a
  xlab("codon number in stretch") +
  ylab("relative stretch position in gene") +
  facet_grid(codon ~ .) +
  ggsave(filename=paste0("stretch_pos_hist.png"), width=13)
