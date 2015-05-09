library(ggplot2)


e.vals <- as.numeric(scan("e_values.dat", what="", sep="\n"))
e.vals <- e.vals[e.vals<=1]
df <- data.frame(e=e.vals)

ggplot(df, aes(x=e, fill=..count..)) +
  geom_histogram(binwidth=0.02) +
  xlab("e value") +
  ggsave(filename="evalue_hist.png")
