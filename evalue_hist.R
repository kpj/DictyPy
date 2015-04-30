library(ggplot2)


e.vals <- as.numeric(scan("e_values.dat", what="", sep="\n"))
e.vals <- e.vals[e.vals<=0.15]
df <- data.frame(e=e.vals)

ggplot(df, aes(x=e, fill=..count..)) +
  geom_histogram(binwidth=0.001) +
  xlab('e value')
